function plot_ROI(T, regionID, marker)

% T: table with columns x and y and region
% regionID: numeric or string indicating region of interest
% marker: cell array of phenotyped marker column names (i.e. {'PD1p'} or
% {'PD1p', 'CD8p'}

    %plot all cells
    scatter(T.tX_centroid, T.tY_centroid, 10, [0.8, 0.8, 0.8], 'filled');
    hold on;
    
    if ~isempty(regionID)
        mask = T.Region == regionID;
        scatter(T.tX_centroid(mask), T.tY_centroid(mask), 10, 'r', 'filled');
    end
    
    if nargin >=3 && ~isempty(marker)
        % combine all marker masks with logical AND
        phenotypeMask = true(height(T), 1);
        for i = 1:numel(marker)
            phenotypeMask = phenotypeMask & T.(marker{i}) == 1;
        end
        scatter(T.tX_centroid(phenotypeMask), T.tY_centroid(phenotypeMask), 12, 'b', 'filled');
    end 
    
    hold off;
    
    axis equal;
    xlabel('X');
    ylabel('Y');
    
    if ~isempty(regionID)
        title(['Region: ', num2str(regionID)]);
    elseif nargin >=3 && ~isempty(markerCols)
        title(['Phenotype: ', strjoin(marker, ' & ')]);
    else
        title('Sample: ', inputname(1));
        
end