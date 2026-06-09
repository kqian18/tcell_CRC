function clusterStats = countMarkersInClusters(dataStruct, slideNames, slideInfo, opts)
% countMarkersInClusters  Count immune markers within lymphoid aggregates.
%
%   clusterStats = countMarkersInClusters(dataStruct, slideNames, slideInfo)
%   clusterStats = countMarkersInClusters(dataStruct, slideNames, slideInfo, opts)
%
%   INPUTS
%     dataStruct  - struct where each field is a slide table (rows = cells,
%                   columns include IC, Region, tX/tY_centroid, and marker
%                   positivity columns ending in 'p', e.g. CD4p, FOXP3p)
%     slideNames  - cell array of slide name strings (must match field names
%                   in dataStruct)
%     slideInfo   - table with a 'slideName' column for joining metadata
%     opts        - (optional) struct with fields:
%                     .radialDensity  : true/false, compute radial density
%                                       (default: false)
%                     .radialRadius   : radius for radialdensity2 (default: 50)
%
%   OUTPUT
%     clusterStats - table with one row per cluster per slide, containing
%                    raw counts, fractions, region labels, and (optionally)
%                    radial density for every marker combination in markerDefs

% -------------------------------------------------------------------------
% OPTIONS
% -------------------------------------------------------------------------
if nargin < 4 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'radialDensity'), opts.radialDensity = false; end
if ~isfield(opts, 'radialRadius'),  opts.radialRadius  = 50;    end

% -------------------------------------------------------------------------
% MARKER DEFINITIONS — edit only this section to add/remove markers
%
% Each entry defines one marker combination to quantify.
%   .label    : output column name suffix
%               → <label>_Pos, Fraction_<label>, RadialDensity_<label>
%   .channels : cell array of column names that must ALL be positive (AND-ed)
% -------------------------------------------------------------------------
markerDefs = struct( ...
    'label',    {'CD4',           'CD4FOXP3',                    'CD4PD1',                     'CD4ICOS',                    'CD4CD45RO',                       'CD4GZMB'}, ...
    'channels', {{'CD4p','CD3dp'}, {'CD4p','FOXP3p','CD3dp'}, {'CD4p','PD_1p','CD3dp'}, {'CD4p','ICOSp','CD3dp'}, {'CD4p','CD45ROp','CD3dp'}, {'CD4p','GranzymeBp','CD3dp'}} ...
);
% -------------------------------------------------------------------------

nMarkers  = numel(markerDefs);
clusterStats = table();

for i = 1:numel(slideNames)
    currSlide = slideNames{i};
    fprintf('Processing: %s\n', currSlide);

    data = dataStruct.(currSlide);

    clusterIDs = unique(data.IC);
    clusterIDs(clusterIDs == 0) = [];

    if isempty(clusterIDs)
        warning('No clusters found for slide: %s', currSlide);
        continue;
    end

    nClusters = numel(clusterIDs);

    % Preallocate
    counts       = zeros(nClusters, nMarkers);
    radialDens   = NaN(nClusters, nMarkers);   % only populated if requested
    cluster_size = zeros(nClusters, 1);
    regionLabels = strings(nClusters, 1);

    for c = 1:nClusters
        mask          = data.IC == clusterIDs(c);
        cluster_cells = data(mask, :);
        cluster_size(c) = sum(mask);

        % Region label: majority vote
        if sum(cluster_cells.Region == 1) > sum(cluster_cells.Region ~= 1)
            regionLabels(c) = "Tumor";
        else
            regionLabels(c) = "Non-Tumor";
        end

        x = cluster_cells.tX_centroid;
        y = cluster_cells.tY_centroid;

        for m = 1:nMarkers
            positiveMask = buildPositiveMask(cluster_cells, markerDefs(m).channels);
            counts(c, m) = sum(positiveMask);

            if opts.radialDensity && counts(c, m) > 0
                radialDens(c, m) = radialdensity2( ...
                    x, y, positiveMask, true(size(positiveMask)), opts.radialRadius);
            end
        end
    end

    % Build per-slide table
    slideStats            = table();
    slideStats.slideName  = repmat(string(currSlide), nClusters, 1);
    slideStats.ClusterID  = clusterIDs;
    slideStats.TotalCells = cluster_size;

    for m = 1:nMarkers
        lbl = markerDefs(m).label;
        slideStats.([lbl, '_Pos'])      = counts(:, m);
        slideStats.(['Fraction_', lbl]) = counts(:, m) ./ cluster_size;
        if opts.radialDensity
            slideStats.(['RadialDensity_', lbl]) = radialDens(:, m);
        end
    end

    slideStats.RegionLabel = regionLabels;
    clusterStats = [clusterStats; slideStats]; %#ok<AGROW>
end

% Join slide metadata
if ~isempty(clusterStats)
    clusterStats = join(clusterStats, slideInfo, 'Keys', 'slideName');
end

end


% -------------------------------------------------------------------------
% Helper: build logical mask where ALL listed channels are positive (> 0)
% -------------------------------------------------------------------------
function mask = buildPositiveMask(cellTable, channels)
    mask = true(height(cellTable), 1);
    for k = 1:numel(channels)
        mask = mask & (cellTable.(channels{k}) > 0);
    end
end