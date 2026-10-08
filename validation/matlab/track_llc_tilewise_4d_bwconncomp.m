function [tracks, info] = track_llc_tilewise_4d_bwconncomp(mask, seamEdges, varargin)
%TRACK_LLC_TILEWISE_4D_BWCONNCOMP Track MHWs by tile-wise 4-D connected components.
%
%   [tracks, info] = track_llc_tilewise_4d_bwconncomp(mask, seamEdges)
%
%   mask is ni x nj x ntile x nz x ntime.
%
%   The algorithm is:
%   1. For each tile, run bwconncomp on mask(:,:,tile,:,:) as a 4-D array
%      [i j depth time] using full 3x3x3x3 connectivity.
%   2. Use exact LLC seamEdges to union components across tile boundaries.
%      Across a seam, a cell is connected to the seam-neighbor cell at
%      depth/time offsets dz,dt in -1:1, matching the full 4-D connectivity
%      used inside each tile.
%   3. Assemble final components as tracks with global idx5.
%
%   Output tracks keep:
%       tracks(k).idx5      cell array of global 5-D indices by day
%       tracks(k).day       global days
%       tracks(k).nvoxels   voxel counts by day

validateattributes(mask, {'numeric','logical'}, {'nonempty'}, mfilename, 'mask', 1);
if ndims(mask) ~= 5
    error('%s:MaskMustBe5D', mfilename, ...
        'mask must be ni x nj x ntile x nz x ntime.');
end
validateattributes(seamEdges, {'numeric'}, {'ncols', 2}, mfilename, 'seamEdges', 2);

dims5 = double(size(mask));
dims4 = dims5(1:4);
ntimeLocal = dims5(5);
dimsTile = [dims4(1), dims4(2), dims4(4), ntimeLocal];
n2 = prod(dims4(1:3));
n4 = prod(dims4);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Time', [], @(x) isempty(x) || numel(x) == ntimeLocal);
addParameter(p, 'MaskSize', dims5, @(x) isnumeric(x) && numel(x) == 5);
addParameter(p, 'MinDuration', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'SaveFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
opt = p.Results;
opt.SaveFile = char(opt.SaveFile);
opt.Verbose = logical(opt.Verbose);

maskSize = double(opt.MaskSize(:)');
if ~isequal(maskSize(1:4), dims4)
    error('%s:MaskSizeMismatch', mfilename, ...
        'MaskSize(1:4) must match size(mask,1:4).');
end
ntimeTotal = maskSize(5);

if isempty(opt.Time)
    dayLabels = (1:ntimeLocal)';
else
    dayLabels = opt.Time(:);
end
if any(dayLabels < 1 | dayLabels > ntimeTotal)
    error('%s:TimeOutOfRange', mfilename, ...
        'Time values must be inside [1, MaskSize(5)].');
end

conn4 = true(3, 3, 3, 3);

seamEdges = uint32(seamEdges);
isSeamCell = false(n2, 1);
isSeamCell(double(seamEdges(:))) = true;

parent = uint32([]);
rank = zeros(0, 1, 'uint8');
compTile = zeros(0, 1, 'uint16');
compPixels = {};

seamKey = zeros(0, 1, 'uint64');
seamComp = zeros(0, 1, 'uint32');

if opt.Verbose
    fprintf('Running tile-wise 4-D bwconncomp...\n');
end

for tile = 1:dims4(3)
    maskTile = reshape(mask(:,:,tile,:,:), dimsTile);
    CC = bwconncomp(maskTile ~= 0, conn4);
    nHere = CC.NumObjects;
    firstID = numel(parent) + 1;
    ids = uint32(firstID:(firstID + nHere - 1))';

    parent = [parent; ids]; %#ok<AGROW>
    rank = [rank; zeros(nHere, 1, 'uint8')]; %#ok<AGROW>
    compTile = [compTile; repmat(uint16(tile), nHere, 1)]; %#ok<AGROW>
    compPixels = [compPixels; CC.PixelIdxList(:)]; %#ok<AGROW>

    for c = 1:nHere
        pix = uint64(CC.PixelIdxList{c}(:));
        [ii, jj, zz, tt] = ind2sub(dimsTile, double(pix));
        hidx = uint64(sub2ind(dims4(1:3), ii(:), jj(:), repmat(tile, numel(ii), 1)));
        keep = isSeamCell(double(hidx));
        if any(keep)
            key = makeKey(hidx(keep), uint64(zz(keep)), uint64(tt(keep)), uint64(n2), uint64(dims4(4)));
            seamKey = [seamKey; key(:)]; %#ok<AGROW>
            seamComp = [seamComp; repmat(ids(c), numel(key), 1)]; %#ok<AGROW>
        end
    end

    if opt.Verbose
        fprintf('  tile %d/%d: %d components\n', tile, dims4(3), nHere);
    end
end

if opt.Verbose
    fprintf('Unioning components across exact LLC seams...\n');
end

if ~isempty(seamKey)
    [seamKey, order] = sort(seamKey);
    seamComp = seamComp(order);
    [seamKey, ia] = unique(seamKey, 'stable');
    seamComp = seamComp(ia);
    [parent, rank, nSeamPairs] = unionExactSeamsVectorized(parent, rank, seamEdges, ...
        seamKey, seamComp, n2, dims4(4), ntimeLocal, opt.Verbose);
else
    nSeamPairs = 0;
end

root = zeros(numel(parent), 1, 'uint32');
for c = 1:numel(parent)
    [parent, root(c)] = findRoot(parent, uint32(c));
end
[rootList, ~, compID] = unique(double(root));
nComp = numel(rootList);

if opt.Verbose
    fprintf('Assembling %d raw tile-wise 4-D components into tracks...\n', nComp);
end

rawTracks = repmat(emptyTrack(), nComp, 1);
for c = 1:numel(compPixels)
    outID = compID(c);
    tile = double(compTile(c));
    pix = compPixels{c};
    if isempty(pix)
        continue
    end

    [ii, jj, zz, tt] = ind2sub(dimsTile, double(pix(:)));
    idx4 = uint64(sub2ind(dims4, ii(:), jj(:), repmat(tile, numel(ii), 1), zz(:)));
    globalDay = uint64(dayLabels(tt(:)));
    idx5 = idx4 + (globalDay - uint64(1)) .* uint64(n4);

    [daysHere, ~, dayID] = unique(tt(:));
    for d = 1:numel(daysHere)
        rows = dayID == d;
        dayLabel = int32(dayLabels(daysHere(d)));
        rawTracks(outID).day(end+1,1) = dayLabel;
        rawTracks(outID).idx5{end+1,1} = idx5(rows);
        rawTracks(outID).nvoxels(end+1,1) = uint32(nnz(rows));
    end
end

for k = 1:numel(rawTracks)
    rawTracks(k) = mergeTrackDays(rawTracks(k));
end

duration = arrayfun(@(s) numel(s.day), rawTracks);
tracks = rawTracks(duration >= opt.MinDuration);

info = struct();
info.method = 'tilewise_4d_bwconncomp_exact_llc_seam_union';
info.mask_size = maskSize;
info.grid_size = dims4;
info.ntime_local = ntimeLocal;
info.ntime_total = ntimeTotal;
info.min_duration = opt.MinDuration;
info.connectivity = 'full 3x3x3x3 within tile; exact seam union with dz,dt in -1:1';
info.n_tile_components = numel(compPixels);
info.n_raw_components = nComp;
info.n_tracks = numel(tracks);
info.n_seam_edges = size(seamEdges, 1);
info.n_active_seam_cells = numel(seamKey);
info.n_seam_component_pairs = nSeamPairs;

if ~isempty(opt.SaveFile)
    save(opt.SaveFile, 'tracks', 'info', '-v7.3');
end
end

function [parent, rank, nPairTotal] = unionExactSeamsVectorized(parent, rank, seamEdges, seamKey, seamComp, n2, nz, ntimeLocal, verbose)
[hAll, zAll, tAll] = parseKeyVector(seamKey, uint64(n2), uint64(nz));
hAllDouble = double(hAll);
[hSorted, hOrder] = sort(hAllDouble);
[hUnique, firstPos] = unique(hSorted, 'first');
[~, lastPos] = unique(hSorted, 'last');

hFirst = zeros(n2, 1, 'uint32');
hLast = zeros(n2, 1, 'uint32');
hFirst(hUnique) = uint32(firstPos);
hLast(hUnique) = uint32(lastPos);

nPairTotal = 0;
nComp = uint64(numel(parent));

for e = 1:size(seamEdges, 1)
    a = double(seamEdges(e,1));
    b = double(seamEdges(e,2));

    pa = hFirst(a);
    if pa == 0
        continue
    end
    rowsA = hOrder(double(pa):double(hLast(a)));
    zA = double(zAll(rowsA));
    tA = double(tAll(rowsA));

    pairA = zeros(0, 1, 'uint32');
    pairB = zeros(0, 1, 'uint32');

    for dz = -1:1
        z2 = zA + dz;
        validZ = z2 >= 1 & z2 <= nz;
        if ~any(validZ)
            continue
        end

        for dt = -1:1
            t2 = tA + dt;
            valid = validZ & t2 >= 1 & t2 <= ntimeLocal;
            if ~any(valid)
                continue
            end

            srcRows = rowsA(valid);
            targetKey = makeKey(repmat(uint64(b), nnz(valid), 1), ...
                uint64(z2(valid)), uint64(t2(valid)), uint64(n2), uint64(nz));
            [tf, loc] = ismember(targetKey, seamKey);
            if any(tf)
                pairA = [pairA; seamComp(srcRows(tf))]; %#ok<AGROW>
                pairB = [pairB; seamComp(loc(tf))]; %#ok<AGROW>
            end
        end
    end

    if ~isempty(pairA)
        keep = pairA ~= pairB;
        pairA = pairA(keep);
        pairB = pairB(keep);
        if ~isempty(pairA)
            pairKey = uint64(pairA) + nComp .* (uint64(pairB) - uint64(1));
            [~, ia] = unique(pairKey);
            pairA = pairA(ia);
            pairB = pairB(ia);
            nPairTotal = nPairTotal + numel(pairA);
            for q = 1:numel(pairA)
                [parent, rank] = unionOne(parent, rank, pairA(q), pairB(q));
            end
        end
    end

    if verbose && (e == 1 || e == size(seamEdges, 1) || mod(e, 1000) == 0)
        fprintf('  seam edge %d/%d; component pairs so far %d\n', ...
            e, size(seamEdges, 1), nPairTotal);
    end
end
end

function key = makeKey(hidx, z, t, n2, nz)
key = uint64(hidx(:)) + (uint64(z(:)) - uint64(1)) .* n2 + ...
    (uint64(t(:)) - uint64(1)) .* n2 .* nz;
end

function [hidx, z, t] = parseKeyVector(key, n2, nz)
x = uint64(key(:)) - uint64(1);
hidx = mod(x, n2) + uint64(1);
x = idivide(x, n2, 'floor');
z = mod(x, nz) + uint64(1);
t = idivide(x, nz, 'floor') + uint64(1);
end

function [parent, rank] = unionOne(parent, rank, a, b)
[parent, ra] = findRoot(parent, a);
[parent, rb] = findRoot(parent, b);
if ra == rb
    return
end

ia = double(ra);
ib = double(rb);
if rank(ia) < rank(ib)
    parent(ia) = rb;
elseif rank(ia) > rank(ib)
    parent(ib) = ra;
else
    parent(ib) = ra;
    rank(ia) = rank(ia) + 1;
end
end

function [parent, root] = findRoot(parent, x)
root = x;
while parent(double(root)) ~= root
    root = parent(double(root));
end
while parent(double(x)) ~= x
    next = parent(double(x));
    parent(double(x)) = root;
    x = next;
end
end

function tr = mergeTrackDays(tr)
if isempty(tr.day)
    return
end

[days, ~, dayID] = unique(tr.day);
idx5New = cell(numel(days), 1);
nvoxNew = zeros(numel(days), 1, 'uint32');
for d = 1:numel(days)
    parts = tr.idx5(dayID == d);
    idx = vertcat(parts{:});
    idx5New{d} = idx(:);
    nvoxNew(d) = uint32(numel(idx));
end

tr.day = days(:);
tr.idx5 = idx5New;
tr.nvoxels = nvoxNew;
end

function s = emptyTrack()
s = struct('day', zeros(0,1,'int32'), ...
           'idx5', {cell(0,1)}, ...
           'nvoxels', zeros(0,1,'uint32'));
end
