function meanDensity = radialdensity2(x, y, isPositive, inMask, radius_um)
% Radial density with batched KD-tree queries
% meanDensity = radialdensity(x, y, isPositive, inMask, radius_um)
% 1) Build one KD-tree on ALL cells in the region (inMask).
% 2) For EVERY cell i (all N), count positive neighbors within r (exclude self if positive).
% 3) Per-cell density = (#pos neighbors)/area; return the MEAN over all cells.
%
% Memory-safe: queries are processed in batches so neighbor lists aren't held for all N cells at once.

    if nargin < 5, radius_um = 50; end
    if isempty(inMask), inMask = true(size(x)); end

    % adjust scale factor depending on microns/pixels:
    px2um = 0.65;  % <- scale factor 
    x = single(x(:) * px2um);
    y = single(y(:) * px2um);

    inMask     = logical(inMask(:));
    isPositive = logical(isPositive(:));

    % Restrict to the region (e.g., tumor or per-cluster)
    x = x(inMask); y = y(inMask);
    isPositive = isPositive(inMask);

    N = numel(x);
    if N == 0
        meanDensity = NaN; return;
    end

    coords = [x y];
    area_mm2 = pi * (radius_um^2) / 1e6;

    % Build one KD-tree (on ALL cells in region)
    mdl = createns(coords, 'NSMethod','kdtree');

    % ---- Batch over ALL cells as queries  ----
    BATCH = 2e4;

    totalPositiveNeighbors = 0;

    for sIdx = 1:BATCH:N
        eIdx = min(sIdx + BATCH - 1, N);
        qIdx = sIdx:eIdx;

        % Neighbor lists for this batch (cell array of indices into coords)
        nbrs = rangesearch(mdl, coords(qIdx, :), radius_um);

        % Count positives per query, excluding the query itself if it is positive
        % (The query index in global masked space is qIdx(k))
        qIsPos = isPositive(qIdx);

        % Sum positives among neighbors for each query
        posCounts = cellfun(@(v) sum(isPositive(v)), nbrs);

        % Exclude self-count for positive queries (self at distance 0)
        posCounts = posCounts - double(qIsPos);

        totalPositiveNeighbors = totalPositiveNeighbors + sum(posCounts);

        % Free batch cell array promptly
        nbrs = []; %#ok<NASGU>
    end

    % Mean per-cell radial density (over ALL cells), in positives per mm^2
    meanDensity = (totalPositiveNeighbors / double(N)) / area_mm2;
end
