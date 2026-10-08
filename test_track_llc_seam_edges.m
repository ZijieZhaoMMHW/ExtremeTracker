function report = test_track_llc_seam_edges(varargin)
%TEST_TRACK_LLC_SEAM_EDGES Exhaustive synthetic tests for LLC seam unions.
%
% This test does not alter seamEdges, the tracker, or any MHW data.  It
% creates isolated two-voxel masks for every selected seam pair and checks
% the public tracker against the prescribed full 3x3x3x3 connectivity.
%
% Examples:
%   report = test_track_llc_seam_edges('MaxEdges', 20);
%   report = test_track_llc_seam_edges;  % all seam pairs

p = inputParser;
addParameter(p, 'SeamFile', 'llc90_seamEdges_from_bounds.mat', @(x) ischar(x) || isstring(x));
addParameter(p, 'BatchSize', 16, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'MaxEdges', inf, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
opt = p.Results;
opt.Verbose = logical(opt.Verbose);

loaded = load(opt.SeamFile, 'seamEdges');
seamEdges = double(loaded.seamEdges);
gridSize = [90 90 13];
nHorizontal = prod(gridSize);

validateSeamTable(seamEdges, gridSize);
nEdges = min(size(seamEdges, 1), floor(opt.MaxEdges));
seamEdges = seamEdges(1:nEdges, :);

if opt.Verbose
    fprintf('Testing %d LLC seam pairs in batches of %d.\n', nEdges, opt.BatchSize);
end

report = struct();
report.seam_file = char(opt.SeamFile);
report.n_edges_tested = nEdges;
report.batch_size = opt.BatchSize;
report.static_validation = 'passed';
report.same_depth_same_time = runCase(seamEdges, gridSize, opt.BatchSize, 'same', opt.Verbose);
report.forward_depth_time_diagonal = runCase(seamEdges, gridSize, opt.BatchSize, 'forwardDiagonal', opt.Verbose);
report.reverse_depth_time_diagonal = runCase(seamEdges, gridSize, opt.BatchSize, 'reverseDiagonal', opt.Verbose);
report.depth_gap_two_control = runCase(seamEdges, gridSize, opt.BatchSize, 'depthGapTwo', opt.Verbose);
report.time_gap_two_control = runCase(seamEdges, gridSize, opt.BatchSize, 'timeGapTwo', opt.Verbose);

if opt.Verbose
    fprintf('All selected LLC seam tests passed.\n');
end
end

function validateSeamTable(seamEdges, gridSize)
nHorizontal = prod(gridSize);
assert(isnumeric(seamEdges) && size(seamEdges, 2) == 2, ...
    'test_track_llc_seam_edges:InvalidSeamTable', 'seamEdges must be N-by-2.');
assert(all(isfinite(seamEdges(:))) && all(seamEdges(:) == round(seamEdges(:))), ...
    'test_track_llc_seam_edges:InvalidIndices', 'seamEdges must contain finite integer indices.');
assert(all(seamEdges(:) >= 1 & seamEdges(:) <= nHorizontal), ...
    'test_track_llc_seam_edges:OutOfRange', 'A seam index lies outside the LLC horizontal grid.');
assert(~any(seamEdges(:, 1) == seamEdges(:, 2)), ...
    'test_track_llc_seam_edges:SelfEdge', 'A seam edge cannot connect a cell to itself.');

canonical = sort(seamEdges, 2);
assert(size(unique(canonical, 'rows'), 1) == size(canonical, 1), ...
    'test_track_llc_seam_edges:DuplicateEdge', 'Duplicate undirected seam pairs were found.');

[~, ~, tileA] = ind2sub(gridSize, seamEdges(:, 1));
[~, ~, tileB] = ind2sub(gridSize, seamEdges(:, 2));
assert(~any(tileA == tileB), ...
    'test_track_llc_seam_edges:SameTileEdge', 'A seam pair must span two LLC tiles.');
end

function result = runCase(seamEdges, gridSize, batchSize, caseName, verbose)
nEdges = size(seamEdges, 1);
nBatches = ceil(nEdges / batchSize);
result = struct('case_name', caseName, 'n_edges', nEdges, 'n_batches', nBatches, ...
    'expected_components_per_edge', expectedComponents(caseName));

for batchIndex = 1:nBatches
    first = (batchIndex - 1) * batchSize + 1;
    last = min(batchIndex * batchSize, nEdges);
    edgeBlock = seamEdges(first:last, :);
    nHere = size(edgeBlock, 1);

    % Four time steps per seam pair guarantee empty temporal buffers between
    % independent synthetic cases, even for full dt = -1:1 connectivity.
    mask = false(gridSize(1), gridSize(2), gridSize(3), 3, 4 * nHere);
    for edgeIndex = 1:nHere
        time0 = 1 + 4 * (edgeIndex - 1);
        mask = placeCase(mask, edgeBlock(edgeIndex, :), gridSize, time0, caseName);
    end

    [tracks, info] = track_llc_tilewise_4d_bwconncomp(mask, seamEdges, ...
        'MinDuration', 1, 'Verbose', false);
    expected = expectedComponents(caseName) * nHere;
    actual = info.n_raw_components;
    assert(actual == expected && numel(tracks) == expected, ...
        'test_track_llc_seam_edges:ConnectivityFailure', ...
        ['%s failed for seam rows %d:%d: expected %d disconnected/connected ', ...
         'components, obtained %d.'], caseName, first, last, expected, actual);

    if verbose && (batchIndex == 1 || batchIndex == nBatches || mod(batchIndex, 20) == 0)
        fprintf('  %s: batch %d/%d passed.\n', caseName, batchIndex, nBatches);
    end
end
end

function mask = placeCase(mask, pair, gridSize, time0, caseName)
[ia, ja, tileA] = ind2sub(gridSize, pair(1));
[ib, jb, tileB] = ind2sub(gridSize, pair(2));

switch caseName
    case 'same'
        mask(ia, ja, tileA, 1, time0) = true;
        mask(ib, jb, tileB, 1, time0) = true;
    case 'forwardDiagonal'
        mask(ia, ja, tileA, 1, time0) = true;
        mask(ib, jb, tileB, 2, time0 + 1) = true;
    case 'reverseDiagonal'
        mask(ib, jb, tileB, 1, time0) = true;
        mask(ia, ja, tileA, 2, time0 + 1) = true;
    case 'depthGapTwo'
        mask(ia, ja, tileA, 1, time0) = true;
        mask(ib, jb, tileB, 3, time0) = true;
    case 'timeGapTwo'
        mask(ia, ja, tileA, 1, time0) = true;
        mask(ib, jb, tileB, 1, time0 + 2) = true;
    otherwise
        error('test_track_llc_seam_edges:UnknownCase', 'Unknown test case: %s', caseName);
end
end

function count = expectedComponents(caseName)
switch caseName
    case {'same', 'forwardDiagonal', 'reverseDiagonal'}
        count = 1;
    case {'depthGapTwo', 'timeGapTwo'}
        count = 2;
    otherwise
        error('test_track_llc_seam_edges:UnknownCase', 'Unknown test case: %s', caseName);
end
end
