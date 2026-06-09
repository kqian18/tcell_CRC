function outsideStats = countMarkersOutsideClusters(dataStruct, slideNames, slideInfo, opts)
% countMarkersOutsideClusters  Count immune markers outside lymphoid aggregates.
%
%   outsideStats = countMarkersOutsideClusters(dataStruct, slideNames, slideInfo)
%   outsideStats = countMarkersOutsideClusters(dataStruct, slideNames, slideInfo, opts)
%
%   Counts cells with IC == 0 (not part of any cluster), split by region
%   (Tumor vs Non-Tumor). Produces one row per slide per region.
%
%   INPUTS
%     dataStruct  - struct where each field is a slide table (rows = cells,
%                   columns include IC, Region, tX/tY_centroid, and marker
%                   positivity columns ending in 'p', e.g. CD4p, FOXP3p)
%     slideNames  - cell array of slide name strings (must match field names
%                   in dataStruct)
%     slideInfo   - table with a 'slideName' column for joining metadata
%     opts        - (optional) struct with fields:
%                     .radialDensity : true/false, compute radial density
%                                      (default: false)
%                     .radialRadius  : radius for radialdensity2
%                                      (default: 50)
%
%   OUTPUT
%     outsideStats - table with one row per slide per region. Columns:
%                    slideName, TotalCells, <label>_Pos, Fraction_<label>,
%                    [RadialDensity_<label>], RegionLabel, plus any
%                    slideInfo metadata columns.
%
% 
%  -------------------------------------------------------------------------
markerDefs = struct( ...
    'label',    {'CD4',           'CD4FOXP3',                    'CD4PD1',                     'CD4ICOS',                    'CD4CD45RO',                       'CD4GZMB'}, ...
    'channels', {{'CD4p','CD3dp'}, {'CD4p','FOXP3p','CD3dp'}, {'CD4p','PD_1p','CD3dp'}, {'CD4p','ICOSp','CD3dp'}, {'CD4p','CD45ROp','CD3dp'}, {'CD4p','GranzymeBp','CD3dp'}} ...
);
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% OPTIONS
% -------------------------------------------------------------------------
if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'radialDensity'), opts.radialDensity = false; end
if ~isfield(opts, 'radialRadius'),  opts.radialRadius  = 50;    end

nMarkers     = numel(markerDefs);
outsideStats = table();

regionDefs = struct( ...
    'label', {"Tumor", "Non-Tumor"}, ...
    'maskFn', { @(data) data.IC == 0 & data.Region == 1, ...
                @(data) data.IC == 0 & data.Region ~= 1  } ...
);

for i = 1:numel(slideNames)
    currSlide = slideNames{i};
    fprintf('Processing: %s\n', currSlide);

    data = dataStruct.(currSlide);

    for r = 1:numel(regionDefs)
        regionLabel  = regionDefs(r).label;
        current_mask = regionDefs(r).maskFn(data);

        if sum(current_mask) == 0
            warning('No outside-cluster %s cells for slide: %s', regionLabel, currSlide);
            continue;
        end

        outside_cells = data(current_mask, :);
        total_cells   = sum(current_mask);

        x = outside_cells.tX_centroid;
        y = outside_cells.tY_centroid;

        counts     = zeros(1, nMarkers);
        radialDens = NaN(1, nMarkers);

        for m = 1:nMarkers
            posMask      = buildPositiveMask(outside_cells, markerDefs(m).channels);
            counts(m)    = sum(posMask);
            if opts.radialDensity && counts(m) > 0
                radialDens(m) = radialdensity2( ...
                    x, y, posMask, true(size(posMask)), opts.radialRadius);
            end
        end

        % Build one-row table for this slide + region
        rowStats            = table();
        rowStats.slideName  = string(currSlide);
        rowStats.TotalCells = total_cells;

        for m = 1:nMarkers
            lbl = markerDefs(m).label;
            rowStats.([lbl, '_Pos'])      = counts(m);
            rowStats.(['Fraction_', lbl]) = counts(m) / total_cells;
            if opts.radialDensity
                rowStats.(['RadialDensity_', lbl]) = radialDens(m);
            end
        end

        rowStats.RegionLabel = regionLabel;
        outsideStats = [outsideStats; rowStats]; %#ok<AGROW>
    end
end

if ~isempty(outsideStats)
    outsideStats = join(outsideStats, slideInfo, 'Keys', 'slideName');
end

end


% -------------------------------------------------------------------------
% Helper: logical mask where ALL listed channels are positive (> 0)
% -------------------------------------------------------------------------
function mask = buildPositiveMask(cellTable, channels)
    mask = true(height(cellTable), 1);
    for k = 1:numel(channels)
        mask = mask & (cellTable.(channels{k}) > 0);
    end
end