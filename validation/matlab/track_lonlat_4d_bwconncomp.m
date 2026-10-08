function [tracks, info] = track_lonlat_4d_bwconncomp(mask, varargin)
%TRACK_LONLAT_4D_BWCONNCOMP Track MHWs on a regular lon-lat-depth-time grid.
%
%   [tracks, info] = track_lonlat_4d_bwconncomp(mask)
%
%   mask is nlon x nlat x nz x ntime.
%
%   This is the regular-grid counterpart of
%   track_llc_tilewise_4d_bwconncomp. The algorithm is equivalent, except
%   that the only nonlocal spatial seam is the periodic longitude boundary:
%
%   1. Run bwconncomp on mask(:,:,:,:) as a 4-D array
%      [lon lat depth time] using full 3x3x3x3 connectivity.
%   2. Union components across the longitude seam between i = 1 and i =
%      nlon. Across that seam, cells are connected with dlat,dz,dt in -1:1,
%      matching the full 4-D connectivity used inside the domain.
%   3. Assemble final components as tracks with global idx4.
%
%   Output tracks keep:
%       tracks(k).idx4      cell array of global 4-D indices by day
%       tracks(k).day       global days
%       tracks(k).nvoxels   voxel counts by day
%
%   Name-value options:
%       'Time'              local-to-global time labels, default 1:ntime
%       'MaskSize'          full [nlon nlat nz ntime_total], default size(mask)
%       'MinDuration'       minimum number of occupied days, default 5
%       'SaveFile'          optional .mat output file
%       'AddIdx5Alias'      also add tracks(k).idx5 = idx4, default false
%       'Verbose'           print progress, default true

validateattributes(mask, {'numeric','logical'}, {'nonempty'}, mfilename, 'mask', 1);

dims4 = double(size(mask));
if numel(dims4) < 4
    dims4(end+1:4) = 1;
elseif numel(dims4) > 4
    error([mfilename ':MaskMustBe4D'], ...
        'mask must be nlon x nlat x nz x ntime.');
end
ntimeLocal = dims4(4);
if dims4(1) < 2
    error([mfilename ':LongitudeTooShort'], ...
        'Longitude dimension must contain at least two cells for periodic boundary tracking.');
end

n2 = prod(dims4(1:2));
n3 = prod(dims4(1:3));

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Time', [], @(x) isempty(x) || numel(x) == ntimeLocal);
addParameter(p, 'MaskSize', dims4, @(x) isnumeric(x) && numel(x) == 4);
addParameter(p, 'MinDuration', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'SaveFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'AddIdx5Alias', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
opt = p.Results;
opt.SaveFile = char(opt.SaveFile);
opt.AddIdx5Alias = logical(opt.AddIdx5Alias);
opt.Verbose = logical(opt.Verbose);

maskSize = double(opt.MaskSize(:)');
if ~isequal(maskSize(1:3), dims4(1:3))
    error([mfilename ':MaskSizeMismatch'], ...
        'MaskSize(1:3) must match size(mask,1:3).');
end
ntimeTotal = maskSize(4);

if isempty(opt.Time)
    dayLabels = (1:ntimeLocal)';
else
    dayLabels = opt.Time(:);
end
if any(dayLabels < 1 | dayLabels > ntimeTotal)
    error([mfilename ':TimeOutOfRange'], ...
        'Time values must be inside [1, MaskSize(4)].');
end

conn4 = true(3, 3, 3, 3);

if opt.Verbose
    fprintf('Running regular-grid 4-D bwconncomp...\n');
end

CC = bwconncomp(reshape(mask ~= 0, dims4), conn4);
nHere = CC.NumObjects;
parent = uint32(1:nHere)';
rank = zeros(nHere, 1, 'uint8');
compPixels = CC.PixelIdxList(:);

compMap = zeros(prod(dims4), 1, 'uint32');
nActiveBoundaryCells = 0;
for c = 1:nHere
    pix = double(compPixels{c}(:));
    compMap(pix) = uint32(c);
    [ii, ~, ~, ~] = ind2sub(dims4, pix);
    nActiveBoundaryCells = nActiveBoundaryCells + nnz(ii == 1 | ii == dims4(1));
end

lonSeamEdges = buildLongitudeSeamEdges(dims4(1:2));

if opt.Verbose
    fprintf('Unioning components across periodic longitude seam...\n');
end

[parent, rank, nSeamPairs] = unionPeriodicLongitudeComponents(parent, rank, ...
    compPixels, compMap, dims4, opt.Verbose);

root = zeros(numel(parent), 1, 'uint32');
for c = 1:numel(parent)
    [parent, root(c)] = findRoot(parent, uint32(c));
end
[rootList, ~, compID] = unique(double(root));
nComp = numel(rootList);

if opt.Verbose
    fprintf('Assembling %d raw regular-grid 4-D components into tracks...\n', nComp);
end

rawTracks = repmat(emptyTrack(), nComp, 1);
for c = 1:numel(compPixels)
    outID = compID(c);
    pix = compPixels{c};
    if isempty(pix)
        continue
    end

    [ii, jj, zz, tt] = ind2sub(dims4, double(pix(:)));
    idx3 = uint64(sub2ind(dims4(1:3), ii(:), jj(:), zz(:)));
    globalDay = uint64(dayLabels(tt(:)));
    idx4 = idx3 + (globalDay - uint64(1)) .* uint64(n3);

    [daysHere, ~, dayID] = unique(tt(:));
    for d = 1:numel(daysHere)
        rows = dayID == d;
        dayLabel = int32(dayLabels(daysHere(d)));
        rawTracks(outID).day(end+1,1) = dayLabel;
        rawTracks(outID).idx4{end+1,1} = idx4(rows);
        rawTracks(outID).nvoxels(end+1,1) = uint32(nnz(rows));
    end
end

for k = 1:numel(rawTracks)
    rawTracks(k) = mergeTrackDays(rawTracks(k));
end

duration = arrayfun(@(s) numel(s.day), rawTracks);
tracks = rawTracks(duration >= opt.MinDuration);

if opt.AddIdx5Alias
    tracks = addIdx5Alias(tracks);
end

info = struct();
info.method = 'regular_lonlat_4d_bwconncomp_periodic_lon_union';
info.mask_size = maskSize;
info.grid_size = dims4(1:3);
info.ntime_local = ntimeLocal;
info.ntime_total = ntimeTotal;
info.min_duration = opt.MinDuration;
info.connectivity = 'full 3x3x3x3 within regular lon-lat grid; periodic lon seam union with dlat,dz,dt in -1:1';
info.index_field = 'idx4';
info.has_idx5_alias = opt.AddIdx5Alias;
info.n_components_before_periodic_union = nHere;
info.n_raw_components = nComp;
info.n_tracks = numel(tracks);
info.n_lon_seam_edges = size(lonSeamEdges, 1);
info.n_active_lon_boundary_cells = nActiveBoundaryCells;
info.n_lon_seam_component_pairs = nSeamPairs;

if ~isempty(opt.SaveFile)
    save(opt.SaveFile, 'tracks', 'info', '-v7.3');
end
end

function [parent, rank, nPairTotal] = unionPeriodicLongitudeComponents(parent, rank, compPixels, compMap, dims4, verbose)
nlon = dims4(1);
nlat = dims4(2);
nz = dims4(3);
ntime = dims4(4);
nComp = uint64(numel(parent));

pairA = zeros(0, 1, 'uint32');
pairB = zeros(0, 1, 'uint32');

for c = 1:numel(compPixels)
    pix = double(compPixels{c}(:));
    if isempty(pix)
        continue
    end

    [ii, jj, zz, tt] = ind2sub(dims4, pix);
    boundary = ii == 1 | ii == nlon;
    if ~any(boundary)
        continue
    end

    ii = ii(boundary);
    jj = jj(boundary);
    zz = zz(boundary);
    tt = tt(boundary);

    iiWrap = ones(size(ii));
    iiWrap(ii == 1) = nlon;

    hits = zeros(0, 1, 'uint32');
    for dj = -1:1
        jj2 = jj + dj;
        validJ = jj2 >= 1 & jj2 <= nlat;
        if ~any(validJ)
            continue
        end

        for dz = -1:1
            zz2 = zz + dz;
            validZ = zz2 >= 1 & zz2 <= nz;
            if ~any(validZ)
                continue
            end

            for dt = -1:1
                tt2 = tt + dt;
                valid = validJ & validZ & tt2 >= 1 & tt2 <= ntime;
                if ~any(valid)
                    continue
                end

                targetPix = sub2ind(dims4, iiWrap(valid), jj2(valid), zz2(valid), tt2(valid));
                hit = compMap(targetPix);
                hits = [hits; hit(hit > 0 & hit ~= uint32(c))]; %#ok<AGROW>
            end
        end
    end

    hits = unique(hits);
    if ~isempty(hits)
        pairA = [pairA; repmat(uint32(c), numel(hits), 1)]; %#ok<AGROW>
        pairB = [pairB; hits(:)]; %#ok<AGROW>
    end

    if verbose && (c == 1 || c == numel(compPixels) || mod(c, 1000) == 0)
        fprintf('  scanned component %d/%d for periodic seam pairs\n', c, numel(compPixels));
    end
end

if isempty(pairA)
    nPairTotal = 0;
    return
end

pairKey = uint64(pairA) + nComp .* (uint64(pairB) - uint64(1));
[~, ia] = unique(pairKey);
pairA = pairA(ia);
pairB = pairB(ia);
nPairTotal = numel(pairA);

for q = 1:nPairTotal
    [parent, rank] = unionOne(parent, rank, pairA(q), pairB(q));
end
end

function seamEdges = buildLongitudeSeamEdges(dims2)
nlon = dims2(1);
nlat = dims2(2);
parts = cell(nlat, 1);
for jj = 1:nlat
    jj2 = max(1, jj - 1):min(nlat, jj + 1);
    left = repmat(sub2ind(dims2, 1, jj), numel(jj2), 1);
    right = sub2ind(dims2, repmat(nlon, numel(jj2), 1), jj2(:));
    parts{jj} = [left, right; right, left];
end
seamEdges = uint32(unique(vertcat(parts{:}), 'rows'));
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
        fprintf('  lon seam edge %d/%d; component pairs so far %d\n', ...
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
x = floor(x ./ n2);
z = mod(x, nz) + uint64(1);
t = floor(x ./ nz) + uint64(1);
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
idx4New = cell(numel(days), 1);
nvoxNew = zeros(numel(days), 1, 'uint32');
for d = 1:numel(days)
    parts = tr.idx4(dayID == d);
    idx = vertcat(parts{:});
    idx4New{d} = idx(:);
    nvoxNew(d) = uint32(numel(idx));
end

tr.day = days(:);
tr.idx4 = idx4New;
tr.nvoxels = nvoxNew;
end

function tracks = addIdx5Alias(tracks)
if isempty(tracks)
    tracks = struct('day', {}, 'idx4', {}, 'nvoxels', {}, 'idx5', {});
    return
end
for k = 1:numel(tracks)
    tracks(k).idx5 = tracks(k).idx4;
end
end

function s = emptyTrack()
s = struct('day', zeros(0,1,'int32'), ...
           'idx4', {cell(0,1)}, ...
           'nvoxels', zeros(0,1,'uint32'));
end
