%% Count markers in lymphoid aggregates/clusters or in non-IC regions
% Both functions accept the same arguments and opts
opts.radialDensity = false;   % or true to enable

newClusterStats  = countMarkersInClusters(dataStruct1, slideNames, slideInfo, opts);
notClusterStats  = countMarkersOutsideClusters(dataStruct1, slideNames, slideInfo, opts);
%% visualize
groupOrder = ["dMMR", "tdpMMR", "tipMMR"];
%allClusterStats.MMR = categorical(allClusterStats.MMR_status, groupOrder, 'Ordinal', true);
tumorClusterStats = newClusterStats_sub(newClusterStats_sub.RegionLabel=="Tumor", :);
tumorNot = notClusterStats_CD4_sub(notClusterStats_CD4_sub.RegionLabel=="Tumor",:);
notTumorCluster = newClusterStats_sub(newClusterStats_sub.RegionLabel=="Non-Tumor", :);
notTumorNot = notClusterStats_CD4_sub(notClusterStats_CD4_sub.RegionLabel=="Non-Tumor",:);
figure;
myboxplot3(tumorNot.Fraction_CD4GZMB, tumorNot.MMR_status);
%violinplot(tumorClusterStats.Fraction_TCF1, tumorClusterStats.MMR);
%violinplot(allClusterStats.RadialDensity_CD8PD1TCF1, allClusterStats.MMR)
%boxplot(allClusterStats.Fraction_Pos, allClusterStats.MMR_status);
ylabel('CD4+ GZMB+ cells');
title('Fraction CD4+ GZMB+ in tumor (exclusive of clusters)');

%saveas(gcf, 'CD8PD1TCF1_tumor_nonTLS.pdf');


%% calculate metrics
clusterSizeVar = "TotalCells";
weightVar = 'TotalCells';   % cluster size
fracVars  = {'Fraction_CD4', ...
    'Fraction_CD4PD1','Fraction_CD4ICOS','Fraction_CD4CD45RO', ...
    'Fraction_CD4FOXP3','Fraction_CD4GZMB'};   % fractions
% radVars   = {'RadialDensity_CD8','RadialDensity_TCF1'...
%     'RadialDensity_CD8PD1TCF1','RadialDensity_CD8TCF1', 'RadialDensity_CD8GZMB', ...
%     'RadialDensity_CD4GZMB', 'RadialDensity_CD4'};        % radial density
countVars = {'CD4_Pos','CD4FOXP3_Pos', ...
    'CD4ICOS_Pos','CD4CD45RO_Pos','CD4PD1_Pos', 'CD4GZMB_Pos'};                  % counts

%% aggregate data per slide
% varsCluster = newClusterStats.Properties.VariableNames;
% varsNon     = notClusterStats_CD4.Properties.VariableNames;
% 
% commonVars = intersect(varsCluster, varsNon, 'stable');  % columns that exist in BOTH
% 
% newClusterStats_sub = newClusterStats(:, commonVars);
% notClusterStats_CD4_sub = notClusterStats_CD4(:, commonVars);


% allStats2 = [newClusterStats_sub; notClusterStats_CD4_sub]; 
varsToKeep = [{'slideName','MMR_status', 'RegionComp', 'Compartment'}, ...
              fracVars, countVars, {weightVar}];

T = allStats2(:, varsToKeep);

% Split cluster vs noncluster
isCluster = T.Compartment == "Cluster";

C = T(isCluster,:);
N = T(~isCluster,:);   % already slide-level

% ---- aggregate CLUSTERS ----
[G, keyTbl] = findgroups(C(:,{'slideName','MMR_status', 'RegionComp', 'Compartment'}));
w = double(C.(weightVar));

aggC = keyTbl;

% Fractions + radial density → cell-weighted mean
for v = fracVars % [fracVars, radVars]
    x = double(C.(v{1}));
    num = splitapply(@(x,w) sum(x.*w,'omitnan'), x, w, G);
    den = splitapply(@(w) sum(w,'omitnan'), w, G);
    aggC.(v{1}) = num ./ den;
end

% Counts → sum
for v = countVars
    x = double(C.(v{1}));
    aggC.(v{1}) = splitapply(@(x) sum(x,'omitnan'), x, G);
end

%aggC.Fraction_TCF1inCD8 = aggC.CD8TCF1_Pos ./ aggC.CD8_Pos;
%aggC.RatioCD4_CD8 = aggC.CD4_Pos ./ aggC.CD8_Pos; 

%N.Fraction_TCF1inCD8 = N.CD8TCF1_Pos ./ N.CD8_Pos;
%N.RatioCD4_CD8 = N.CD4_Pos ./ N.CD8_Pos;


% Drop weight column (it was only needed for aggregation)
if ismember(weightVar, N.Properties.VariableNames)
    N(:, weightVar) = [];
end

% ---- stack back together ----
agg = [aggC; N];

% sanity check
assert(all(countcats(categorical(agg.slideName)) <= 4), ...
    'Aggregation failed: more than one row per slide × RegionComp');

%% plots 

% per MMR
%plotMetric_4Comp_boxplots_perMMR(agg, 'Fraction_CD8');
plotMetric_4Comp_boxplots_perMMR(agg_sub, 'RatioCD4_CD8');
%plotMetric_4Comp_boxplots_perMMR(agg, 'CD8_Pos');

% plot pie charts

plotCD4CD8_pooledPies_perMMR(agg_sub);

% Cluster vs NonCluster
fc1 = log2FC_clusterVsNoncluster_byRegion(agg, metric);
plotLog2FC_byMMR(fc1, 'log2FC_Tumor', ...
    sprintf('%s: Cluster vs NonCluster (Tumor)', metric));
plotLog2FC_byMMR(fc1, 'log2FC_NonTumor', ...
    sprintf('%s: Cluster vs NonCluster (Non-Tumor)', metric));

% Tumor vs NonTumor
fc2 = log2FC_tumorVsNontumor_byCompartment(agg, metric);
plotLog2FC_byMMR(fc2, 'log2FC_Cluster', ...
    sprintf('%s: Tumor vs Non-Tumor (Cluster)', metric));
plotLog2FC_byMMR(fc2, 'log2FC_NonCluster', ...
    sprintf('%s: Tumor vs Non-Tumor (NonCluster)', metric));

%% plot log2FC
metric = 'Fraction_TCF1';
groupLabels = {'tdpMMR', 'dMMR', 'tipMMR'};
groupOrder = categorical(fc.MMR_status, groupLabels, 'Ordinal', true);

fc = log2FC_clusterVsNoncluster_byRegion(agg_sub, metric);

figure;
%boxchart(categorical(fc.MMR_status), fc.log2FC_NonTumor); %or fc.log2FC_NonTumor
myboxplot3(fc.log2FC_Tumor, groupOrder);
yline(0,'k--');
xlabel('MMR Status');
ylabel(sprintf('log_2(Cluster / NonCluster): %s', metric), 'Interpreter','none');
title(sprintf('%s — Tumor', metric), 'Interpreter','none');
set(gca,'Box','off');

%% tumor nontumor fold change
metric = 'Fraction_TCF1';
fc2 = log2FC_tumorVsNontumor_byCompartment(agg_sub, metric);

% Cluster compartment only
figure;
%boxchart(categorical(fc2.MMR_status), fc2.log2FC_Cluster);
myboxplot3(fc2.log2FC_NonCluster, categorical(fc2.MMR_status));
yline(0,'k--');
xlabel('MMR Status');
ylabel(sprintf('log_2(Tumor / Non-Tumor): %s', metric), 'Interpreter','none');
title([metric ' — Not Cluster'], 'Interpreter','none');
set(gca,'Box','off');

%% visualize tumor clusters
groupOrder = ["tdpMMR", "dMMR", "tipMMR"];
allClusterStats.MMR = categorical(allClusterStats.MMR_status, groupOrder, 'Ordinal', true);
tumorClusterStats = allClusterStats(allClusterStats.RegionLabel=="Tumor", :);
figure;
myboxplot3(tumorClusterStats.CD8PD1TCF1_Pos, tumorClusterStats.MMR);
%violinplot(tumorClusterStats.Fraction_TCF1, tumorClusterStats.MMR);
%violinplot(allClusterStats.RadialDensity_CD8PD1TCF1, allClusterStats.MMR)
%boxplot(allClusterStats.Fraction_Pos, allClusterStats.MMR_status);
ylabel('TCF1+ Cells');
title('TCF1+ in lymphoid aggregates in tumor compartment)');

%saveas(gcf, 'CD8PD1TCF1_tumor_nonTLS.pdf');

%% are there more TCF1+ cells in larger clusters?
T = tumorClusterStats;
figure;
scatter(T.TotalCells, T.TCF1_pos, 10, 'filled');
xlabel('Total cells in cluster');
ylabel('TCF1+ cells in cluster');
title('TCF1+ count vs. cluster size');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on; 

%% plot across MMR groups
T = allClusterStats(allClusterStats.RegionLabel=='Non-Tumor', :); 
T.MMR_status = categorical(T.MMR_status);

groups = categories(T.MMR_status);
colors = lines(numel(groups));

figure; hold on;
for i = 1:numel(groups)
    mask = T.MMR_status == groups{i};
    scatter(T.TotalCells(mask), T.CD8PD1TCF1_Pos(mask), ...
        30, colors(i,:), 'filled', 'DisplayName', char(groups{i}));
end
set(gca,'XScale','log','YScale','log')
xlabel('Total cells in cluster')
ylabel('TCF1^+ cells')
title('TCF1^+ vs Cluster Size by MMR group')
legend('Location','best')
grid on

%% subplots
groups = categories(T.MMR_status);

figure;
for i = 1:numel(groups)
    subplot(1, numel(groups), i)
    mask = T.MMR_status == groups{i};

    scatter(T.TotalCells(mask), T.CD8PD1TCF1_Pos(mask), 15, 'filled')
    set(gca,'XScale','log','YScale','log')
    title(groups{i})
    xlabel('Total cells')
    ylabel('TCF1^+ cells')
    grid on
end
sgtitle('TCF1^+ vs Cluster Size by MMR Group')

%% scatter + regression
figure; hold on;

groups = categories(T.MMR_status);
colors = lines(numel(groups));

for i = 1:numel(groups)

    % ---- SELECT CURRENT GROUP ----
    mask = T.MMR_status == groups{i};

    x = T.TotalCells(mask);
    y = T.CD8PD1TCF1_Pos(mask);

    % Remove zeros or negatives for log-space operations
    valid = x > 0 & y > 0;
    x = x(valid);
    y = y(valid);

    % ---- SCATTER PLOT ----
    scatter(x, y, 18, colors(i,:), 'filled', ...
        'DisplayName', char(groups{i}));

    % ---- FIT LOG-LOG REGRESSION ----
    x_log = log10(x);
    y_log = log10(y);

    mdl = fitlm(x_log, y_log);

    % Create smooth fitted line
    xx = linspace(min(x), max(x), 200);
    yy = 10.^(predict(mdl, log10(xx)'));

    % Plot fit line
    plot(xx, yy, 'Color', colors(i,:), 'LineWidth', 1.8);

    % ---- COMPUTE STATISTICS ----
    r2  = mdl.Rsquared.Ordinary;
    [rho, pval] = corr(x, y, 'Type','Spearman');

    % ---- ADD TEXT LABEL NEAR UPPER LEFT OF EACH GROUP CLOUD ----
    tx = min(x) * 1.2;
    ty = max(y) / 1.3;

    %txt = sprintf('%s: R^2=%.2f, rho=%.2f (p=%.1g)', char(groups{i}), r2, rho, pval);

    txt = sprintf('%s:  R^2=%.2f', char(groups{i}), r2);
    
    text(tx, ty, txt, 'Color', colors(i,:), 'FontSize', 10, 'FontWeight','bold');
end

set(gca,'XScale','log','YScale','log')
xlabel('Total cells in cluster')
ylabel('TCF1^+ cells')
title('TCF1^+ vs Cluster Size by MMR Group (log–log)')
legend('Location','best')
grid on


%% spearman
for i = 1:numel(groups)
    mask = T.MMR_status == groups{i};
    [rho, pval] = corr(T.TotalCells(mask), T.Marker2_Pos(mask), 'Type','Spearman');
    fprintf('%s: Spearman rho = %.3f (p = %.3g)\n', groups{i}, rho, pval);
end


%% statistical tests - KW, Levene and Brown-Forsythe

x = tumorClusterStats.Marker2_Pos;
g = categorical(tumorClusterStats.MMR_status, {'tdpMMR', 'dMMR', 'tipMMR'});

%Kruskal-Wallis

[p_kw, tbl_kw, stats_kw] = kruskalwallis(x,g,'off');
fprintf('Krusal-Wallis p = %.3g\n', p_kw);


% Levene's test (mean-based)
[p_levene, stats_levene] = vartestn(x,g,...
    'TestType', 'LeveneAbsolute', ...
    'Display', 'off');

% Brown-Forsythe test (median-based)
[p_bf, stats_bf] = vartestn(x,g,...
    'TestType', 'BrownForsythe', ...
    'Display', 'off');

%% GLM
T = tumorClusterStats;
T.MMR_status = categorical(T.MMR_status, {'dMMR', 'tdpMMR', 'tipMMR'}); % dMMR reference group
validMask = T.TotalCells > 0 & T.CD8PD1TCF1_Pos >=0;
T = T(validMask,:);

offsetTerm = log(T.TotalCells);
mdl = fitglm(T, ...
    'ResponseVar', 'CD8PD1TCF1_Pos', ...
    'PredictorVars', 'MMR_status', ...
    'Distribution', 'poisson', ...
    'Offset', offsetTerm);

disp(mdl);

%anova_mdl = anova(mdl, 'summary');
%disp(anova_mdl);

%% aggregate (so the measurement is per-tumor/group)
T = allClusterStats(allClusterStats.RegionLabel=='Non-Tumor', :);

T.MMR_status = categorical(T.MMR_status);

% keep only rows with sensible values
valid = T.TotalCells > 0 & T.Marker2_Pos >= 0 & ~isundefined(T.MMR_status);
T = T(valid, :);

% ----- PER-TUMOR AGGREGATION -----
% group by SampleID and MMR_status
tumorSummary = groupsummary(T, {'slideName','MMR_status'}, ...
    'sum', {'Marker2_Pos','TotalCells'});

% nicer variable names
tumorSummary.Properties.VariableNames( ...
    end-1:end) = {'Sum_TCF1', 'Sum_TotalCells'};

% per-tumor rate of TCF1+ per cell
tumorSummary.TCF1_rate = tumorSummary.Sum_TCF1 ./ tumorSummary.Sum_TotalCells;

%% compare tumors across MMR groups
x = tumorSummary.CD8PD1TCF1_rate;
g = tumorSummary.MMR_status;

[p_kw, tbl_kw, stats_kw] = kruskalwallis(x, g, 'off');
fprintf('Kruskal–Wallis on per-tumor TCF1_rate: p = %.3g\n', p_kw);

%% pairwise
if p_kw < 0.05
    post = multcompare(stats_kw, 'CType', 'dunn-sidak', 'Display', 'off');
    postTbl = array2table(post, ...
        'VariableNames', {'Group1','Group2','LowerCI','Diff','UpperCI','pValue'});
    disp(postTbl);
end

%% visualize per-tumor distributions
figure;
groupLabels = {'tdpMMR', 'dMMR', 'tipMMR'};
%groupOrder = categorical(tumorSummary.MMR_status, groupLabels, 'Ordinal', true);

myboxplot3_mod(tumorSummary.TCF1_rate, tumorSummary.MMR_status, groupLabels);
ylabel('Per-tumor TCF1^+ rate');
title('Per-tumor TCF1^+ rate by MMR group');

% can optionally log the rate

%% summary stats test

groups = categories(g);
nG = numel(groups);

N = zeros(nG,1);
Mean = zeros(nG,1);
Median = zeros(nG,1);
Std = zeros(nG,1);
VarVal = zeros(nG,1);

for i = 1:nG
    mask = (g == groups{i});
    xi = x(mask);
    
    N(i) = numel(xi);
    Mean(i) = mean(xi);
    Median(i) = median(xi);
    VarVal(i) = var(xi,1);
    Std(i) = std(xi,1);
end

summarybyMMR = table(groups, N, Mean, Median, Std, VarVal, ...
    'VariableNames', {'MMR_status', 'N', 'Mean', 'Median', 'Std', 'Var'});


%% visualize barchart
% Make sure the group label is categorical and ordered
subsetCluster.Group = strcat(subsetCluster.MMR_status, " - ", subsetCluster.RegionLabel);
groupCats = categories(categorical(subsetCluster.Group));

figure; hold on;

colors = lines(length(groupCats));  % Or use your own colormap

for i = 1:length(groupCats)
    thisGroup = groupCats{i};
    mask = subsetCluster.Group == thisGroup;
    boxchart(...
        repmat(i, sum(mask), 1), ...
        subsetCluster.Fraction_CD8PD1TCF1(mask), ...
        'BoxFaceColor', colors(i,:));
end

set(gca, 'XTick', 1:length(groupCats), 'XTickLabel', groupCats);
xtickangle(45);
ylabel('Fraction CD8+PD1+TCF1+ per Cluster');
title('Cluster Fraction by Region and MMR Status');

%%
subsetData = allClusterStats(allClusterStats.RegionLabel == "Tumor", :);
figure;
%myboxplot2_log(subsetData.RadialDensity_CD8PD1TCF1, subsetData.RegionLabel)
myboxplot3(subsetData.Fraction_CD8PD1TCF1, subsetData.MMR_status)
ylabel('CD8+ PD1+ TCF1+ Density');
title('CD8+ PD1+ TCF1+ in lymphoid aggregates in tumor');

%% plot cluster size v radial density

figure;
gscatter(allClusterStats.Fraction_CD8PD1TCF1, ...
    allClusterStats.RadialDensity_CD8PD1TCF1, ...
    allClusterStats.RegionLabel);
xlabel('Fraction CD8+PD1+TCF1+ per cluster');
ylabel('Density of CD8+PD1+TCF1+ (cells/mm2) in clusters');

%% viz MMR TCF density vs fraction
T = allClusterStats_sub;
T.MMR_status = categorical(T.MMR_status);

mmrGroups = categories(T.MMR_status);

figure;
for i = 1:numel(mmrGroups)
    mask = T.MMR_status == mmrGroups{i};

    subplot(1, numel(mmrGroups), i);
    h = gscatter(T.Fraction_CD8GZMB(mask), ...
             T.RadialDensity_CD8GZMB(mask), ...
             T.RegionLabel(mask));
    for j = 1:numel(h)
        h(j).MarkerFaceColor = 'none';
    end
    
    %ylim([0 500]);
         
    xlabel('Fraction CD8+GZMB+ per cluster');
    ylabel('Density of CD8+GZMB+ (cells/mm^2) in clusters');
    %ylabel('Total cells in clusters');
    title(char(mmrGroups{i}));
    legend('Location', 'best');
    grid on;
end


%% visualize
groupOrder = ["tdpMMR", "dMMR", "tipMMR"];
allClusterStats.MMR = categorical(allClusterStats.MMR_status, groupOrder, 'Ordinal', true);

logVals = log10(allClusterStats.RadialDensity_CD8PD1TCF1) + 1e-4;

figure;
myboxplot3(allClusterStats.RadialDensity_CD8PD1TCF1, allClusterStats.MMR);
%violinplot(allClusterStats.RadialDensity_CD8PD1TCF1, allClusterStats.MMR)
%boxplot(allClusterStats.Fraction_Pos, allClusterStats.MMR_status);

% groupedScatter = gscatter(...
%     repmat(double(allClusterStats.MMR), 1, 1) + 0.1*randn(size(allClusterStats.MMR)), ...
%     allClusterStats.RadialDensity_CD8PD1TCF1, ...
%     allClusterStats.MMR);
ylabel('CD8+ PD1+ TCF1+ (cells/mm2)');
title('% CD8+ PD1+ TCF1+ in Lymphoid Aggregate (CD20+ cluster)');


%% boxplot of cluster sizes by MMR

figure; 
boxplot(allClusterStats.TotalCells, allClusterStats.MMR_status);
ylabel('Cluster size (# CD20+ cells)');


%% Ranked barplots for specific markers 
%allClusterStats = join(allClusterStats,slideInfo,'Keys','slideName');
summaryStats = groupsummary(allClusterStats, 'slideName', 'mean', 'Fraction_CD8PD1TCF1');
% Step 2: Sort data
sortedData = sortrows(summaryStats, 'mean_Fraction_CD8PD1TCF1', 'descend');
sortedData = join(sortedData,slideInfo,'Keys','slideName');

% Step 3: Generate categorical labels combining IO and mutation status
%sortedData.Labels = strcat(sortedData.MMR_status, '-', sortedData.slideName);
sortedData.Labels = strcat(sortedData.MMR_status);

% Step 4: Create the ranked bar plot
figure;
bar(sortedData.mean_Fraction_CD8PD1TCF1, 'FaceColor', 'b'); % Blue bars

% Step 5: Formatting
xticks(1:height(sortedData)); % Set x-ticks for each slide
xticklabels(sortedData.Labels); % Use IO - mutstatus as labels
xtickangle(45); % Rotate labels for readability
ylabel('Mean CD8+PD1+TCF1+ Fraction per Cluster');
xlabel('Samples');
title('Ranked CD8+PD1+TCF1+ Cluster Fractions by Slide');
grid on;

% Step 6: Adjust figure for better readability
set(gca, 'FontSize', 10); % Adjust font size for readability

%% Ranked by MMR status
% Step 1: Summary and join with metadata
summaryStats = groupsummary(allClusterStats, 'slideName', 'mean', 'Fraction_CD8PD1TCF1');
summaryStats = join(summaryStats, slideInfo, 'Keys', 'slideName');

% Step 2: Define MMR group order
groupOrder = ["tdpMMR", "dMMR", "tipMMR"];  % Your desired order
sortedGrouped = table();

% Step 3: Loop over groups, sort each group internally
for i = 1:numel(groupOrder)
    g = groupOrder(i);
    gData = summaryStats(summaryStats.MMR_status == g, :);
    gData = sortrows(gData, 'mean_Fraction_CD8PD1TCF1', 'descend');
    sortedGrouped = [sortedGrouped; gData];  % Append in order
    
end

% Step 4: Categorical labels (can include slideName too if needed)
sortedGrouped.Labels = strcat(sortedGrouped.MMR_status);

% Step 5: Bar plot
figure;
bar(sortedGrouped.mean_Fraction_CD8PD1TCF1, 'FaceColor', 'b');

% Step 6: Formatting
xticks(1:height(sortedGrouped));
xticklabels(sortedGrouped.Labels);
xtickangle(45);
ylabel('Mean CD8+PD1+TCF1+ Fraction per Cluster');
xlabel('Samples (grouped by MMR)');
title('CD8+PD1+TCF1+ Fraction Ranked Within Each MMR Group');
grid on;
set(gca, 'FontSize', 10);




