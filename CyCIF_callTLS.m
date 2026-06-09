% ==== User Input ====
output1 = input('Please input the output name: ', 's');
flag1 = input('Output maps (true/false): ');
epsilon = 100;
minPts = 300;

% Define immune marker list
immune_markers = {'CD3dp','CD4p','CD8ap','CD20p','CD45p','CD45ROp','CD68p'};

%% Start parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool('local');
end

%% Process each slide
for i = 1:length(slideName)
    disp(['Processing: ', slideName{i}]);
    data0 = eval(['dataStruct1.', slideName{i}]);
    %data0.cellID = (1:height(data0))';

    coords = [data0.tX_centroid, data0.tY_centroid];
    bCells = data0.CD20p > 0;
    idx_b = find(bCells);

    % Step 1: CD20+ cells with >25% CD20+ in 25-NN
    cd20_mask_local = false(height(data0), 1);
    parfor j = 1:length(idx_b)
        pt = coords(idx_b(j), :);
        dists = sqrt(sum((coords - pt).^2, 2));
        [~, nn_idx] = sort(dists);
        nn_idx = nn_idx(2:26);  % exclude self
        cd20_frac = sum(data0.CD20p(nn_idx) > 0) / 25;
        if cd20_frac > 0.25
            cd20_mask_local(j) = true;
        end
    end
    
    cd20_mask = false(height(data0),1);
    cd20_mask(idx_b(cd20_mask_local)) = true;

    % Step 2: Expand to include immune marker neighbors
    
    idx_seed = find(cd20_mask);
    expanded_mask_local = false(length(idx_seed), 1);
    neighbor_indices_cell = cell(length(idx_seed),1);

    parfor j = 1:length(idx_seed)
        pt = coords(idx_seed(j), :);
        dists = sqrt(sum((coords - pt).^2, 2));
        [~, nn_idx] = sort(dists);
        nn_idx = nn_idx(2:101);  % 100 nearest neighbors

        % Check for immune marker positivity
        immune_pos = false(length(nn_idx), 1);
        for k = 1:length(immune_markers)
            immune_pos = immune_pos | (data0{nn_idx, immune_markers{k}} > 0);
        end

        if sum(immune_pos)/100 >= 0.20
            expanded_mask_local(j)= true;
            neighbor_indices_cell{j} = nn_idx(immune_pos); % store contributing neighbors
        else
            neighbor_indices_cell{j} = [];
        end
    end
    
    expanded_mask = false(height(data0),1);
    expanded_mask(idx_seed(expanded_mask_local)) = true;
    for j = 1:length(neighbor_indices_cell)
        expanded_mask(neighbor_indices_cell{j}) = true;
    end
    
    coords = [data0.tX_centroid, data0.tY_centroid]; 
    assert(length(expanded_mask) == size(coords,1), 'Mismatch between expanded mask and coords');
    
    fprintf('CD20+ seeds considered: %d\n', length(idx_seed)); 
    fprintf('CD20+ seeds expanded: %d\n', sum(expanded_mask_local));
    fprintf('Total cells in expanded_mask/DBSCAN input rows: %d\n', sum(expanded_mask));
    fprintf('Length of expanded_mask: %d\n', length(expanded_mask)); 
    fprintf('coords rows: %d\n', size(coords,1));
    
   

    % Step 3: DBSCAN clustering
%     cluster_labels = zeros(height(data0),1);
%     if sum(expanded_mask) > 0
%         X = coords(expanded_mask, :) * 0.65; % scaling
%         [labels, ~] = dbscan(X, epsilon, minPts);
%         
%         
%         disp('unique labels (pre-conversion):');
%         disp(unique(labels));
%         
%         % convert noise to 0
%         labels(labels < 0) = 0;
%         
%         
%         cluster_labels(expanded_mask) = labels;
%     end


    % pcsegdist instead of DBSCAN
    cluster_labels = zeros(height(data0), 1);
    
    if sum(expanded_mask) > 0
        X = coords(expanded_mask,:) * 0.65;
        pc = pointCloud([X, ones(size(X,1),1)]);
        minDist = 50;
        labels = pcsegdist(pc, minDist);
        label_counts = tabulate(labels);
        valid_labels = label_counts(label_counts(:,2) >= 300,1);
        labels(~ismember(labels,valid_labels)) = 0;
        
        cluster_labels(expanded_mask) = labels;
        
        fprintf('clusters: %d\n', length(valid_labels) - any(valid_labels ==0));
    end

    % Store result
    data0.(output1) = cluster_labels;
    eval(['dataStruct1.', slideName{i}, ' = data0;']);
    disp(tabulate(data0.(output1)));
end

%% Output segmentation maps
if flag1
    for i = 1:length(slideName)
        disp(['Plotting: ', slideName{i}]);
        data1 = eval(['dataStruct1.', slideName{i}]);
        figure('units','normalized','outerposition',[0 0 1 1]);
        CycIF_tumorview(data1, output1, 9, 100000);
        daspect([1 1 1]);
        set(gcf, 'color', 'w');
        title(slideName{i});
        filename = strcat(slideName{i}, '_', output1, '.png');
        saveas(gcf, filename);
        close;
    end
end
