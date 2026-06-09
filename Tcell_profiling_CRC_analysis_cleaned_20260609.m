%% Analysis CRC T cell profiling CyCIF
%  last modified: Kristin Qian 2026/06/09

%% run CyCIF_importgroup_KQ.m
% imported single cell tables into data structure 

%% Generate summary tables

%run CyCIF_alldata_kq to concatenate all tables into alldata;

tic;

sumAll = varfun(@mean,alldata,'GroupingVariables','slideName');
sumAll = join(sumAll,slideInfo,'Keys','slideName');
toc;


sumTumor = varfun(@mean,alldata(alldata.Region==1,:),'GroupingVariables','slideName');
sumTumor = join(sumTumor,slideInfo,'Keys','slideName');
toc;

% allNormal (or allMucosa/not Tumor)
tic; 
sumMucosa = varfun(@mean,alldata(~alldata.Region==1,:),'GroupingVariables','slideName');
sumMucosa = join(sumMucosa,slideInfo,'Keys','slideName');
toc;

% epithelial: panCK/E_cadherin pos, CD45 neg, and SMA neg.  
%immune = cd45+, 
%stroma = SMA pos, panCK/E_cadherin neg, CD45 neg 

% in tumor
sumTumorStroma = varfun(@mean,alldata(alldata.CD45p==0 & ...
    alldata.Region==1 & alldata.SMAp == 1 & ...
    (alldata.Pan_CKp == 0 | alldata.E_cadherinp ==0),:),'GroupingVariables','slideName');
sumTumorStroma = join(sumTumorStroma,slideInfo,'Keys','slideName');
toc;

sumTumorImmune = varfun(@mean, alldata(alldata.CD45p == 1 & alldata.Region==1,:), 'GroupingVariables', 'slideName');
sumTumorImmune = join(sumTumorImmune, slideInfo, 'Keys', 'slideName');


sumTumorEpi = varfun(@mean,alldata(alldata.CD45p == 0 & alldata.SMAp == 0 & alldata.Region==1 & (alldata.Pan_CKp == 1 | alldata.E_cadherinp == 1), :),'GroupingVariables','slideName');
sumTumorEpi = join(sumTumorEpi, slideInfo, 'Keys', 'slideName');
toc;

% non tumor
% epithelial: panCK/E_cadherin pos, CD45 neg, and SMA neg.  
%immune = cd45+, 
%stroma = SMA pos, panCK/E_cadherin neg, CD45 neg

sumStroma = varfun(@mean,alldata(alldata.CD45p==0 & ~alldata.Region== 1 & ...
    (alldata.Pan_CKp == 0 | alldata.E_cadherinp ==0) & ...
    alldata.SMAp == 1,:),'GroupingVariables','slideName');
sumStroma = join(sumStroma,slideInfo,'Keys','slideName');
toc;

sumImmune = varfun(@mean, alldata(alldata.CD45p == 1 & ~alldata.Region==1,:), 'GroupingVariables', 'slideName');
sumImmune = join(sumImmune, slideInfo, 'Keys', 'slideName');
toc;

sumEpi = varfun(@mean,alldata(~alldata.Region==1 & alldata.CD45p == 0 & alldata.SMAp == 0 & (alldata.Pan_CKp == 1 | alldata.E_cadherinp == 1), :),'GroupingVariables','slideName');
sumEpi = join(sumEpi, slideInfo, 'Keys', 'slideName');
toc;
% lymphoid aggregates
sumLA = varfun(@mean,alldata(alldata.Region==2,:),'GroupingVariables','slideName');
sumLA = join(sumLA,slideInfo,'Keys','slideName');
toc;

sumTLS = varfun(@mean, alldata(alldata.Region==3,:), 'GroupingVariables', 'slideName');
sumTLS = join(sumTLS, slideInfo, 'Keys', 'slideName');

sumImmune = varfun(@mean, alldata(alldata.onlyImmune == 1,:), 'GroupingVariables', 'slideName');
sumImmune = join(sumImmune, slideInfo, 'Keys', 'slideName');
toc;

sumtumorImmune = varfun(@mean, alldata(alldata.onlyImmune == 1 & alldata.Region==1,:), 'GroupingVariables', 'slideName');
sumtumorImmune = join(sumtumorImmune, slideInfo, 'Keys', 'slideName');

sumstromaImmune = varfun(@mean, alldata(alldata.onlyImmune == 1 & alldata.Region==0,:), 'GroupingVariables', 'slideName');
sumstromaImmune = join(sumstromaImmune, slideInfo, 'Keys', 'slideName');

sumPanCK = varfun(@mean,alldata(alldata.Region==1 & alldata.panCKp == 1, :),'GroupingVariables','slideName');
sumPanCK = join(sumPanCK, slideInfo, 'Keys', 'slideName');
toc;

sumstromaPanCK = varfun(@mean,alldata(alldata.Region==0 & alldata.panCKp == 1, :),'GroupingVariables','slideName');
sumstromaPanCK = join(sumstromaPanCK, slideInfo, 'Keys', 'slideName');

sumCD45posTumor = varfun(@mean,alldata(alldata.Region==1 & alldata.CD45p == 1, :),'GroupingVariables','slideName');
sumCD45posTumor = join(sumCD45posTumor, slideInfo, 'Keys', 'slideName');

sumCD45pos = varfun(@mean,alldata(alldata.CD45p == 1, :),'GroupingVariables','slideName');
sumCD45pos = join(sumCD45pos, slideInfo, 'Keys', 'slideName');

%% Boxplot (all intensity)

sum2 = sumTumor;
title1 = '';

% Define the two types you want to include
selectedTypes = {'mutant', 'WT'};

% Filter the rows in summarytable where Type matches one of the selectedTypes
filteredRows = ismember(sum2.RNF43_659_status, selectedTypes);

% Apply the filtering
filteredTable = sum2(filteredRows, :);


figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(label1)
    marker1 = strcat('mean_',label1{i});
    subplot(3,7,i);
    %myboxplot3(log(sum2{:,marker1}),sum2.RNF43_659_status);
    myboxplot2d(log(filteredTable{:,marker1}),filteredTable.RNF43_659_status);
    % Use filteredTable in your function
    %myboxplot2d(log(filteredTable{:, marker1}), filteredTable.RNF43_659_status);
    ylabel('Intensity(log)');
    title(label1{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle(title1);


%% Boxplot (all gated markers)

%sum2 = sumTumor;
title1 = '';

% Define the two types you want to include
%selectedTypes = {'yes', 'tdpMMR'};

% Filter the rows in summarytable where Type matches one of the selectedTypes
%filteredRows = ismember(sum2.RNF43_659_status, selectedTypes);

% Apply the filtering
%filteredTable = sum2(filteredRows, :);

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(sublabels)
    %marker1 = strcat('mean_',markerCols{i});
    marker1 = sublabels{i};
    subplot(2,3,i);
    %myboxplot3(sumAll_sub{:,marker1}*100,sumAll_sub.MMR_status);
    myboxplot3_nosig(sumTumor_subset{:,marker1}*100,sumTumor_subset.MMR_status);
    ytickformat('percentage');
    title(sublabels{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle(title1);

%% Boxplot double gates 

% Filter the rows in summarytable where Type matches one of the selectedTypes

% Apply filtering if needed
%filteredTable = sumPanCK(filteredRows, :);
%sum2 = sumTumor;

groupLabels = {'tdpMMR', 'dMMR', 'tipMMR'};
groupOrder = categorical(sumTumor4.MMR_status, groupLabels, 'Ordinal', true);

figure('units','normalized','outerposition',[0.2 0.2 0.7 0.6]);
for i = 1:length(labels)
    %marker1 = strcat('mean_',label4{i});
    marker1 = strcat(labels{i});
    subplot(3,6,i);
    myboxplot3(sumTumor4{:,marker1}*100,sumTumor4.MMR_status);
    ytickformat('percentage');
    title(labels{i},'Interpreter','none');
end

set(gcf,'color','w');
sgtitle(title1);


%% what percentage of cytotoxic T cells are CD103+?

base = alldata.CD4p & alldata.GranzymeBp & ~alldata.FOXP3p & alldata.Region==1;     % CD8+ GZMB+ cells
sub  = base & alldata.CD103p;          % ...that are also CD103+

[G, slides] = findgroups(alldata.slideName);

den = splitapply(@sum, base, G);
num = splitapply(@sum, sub,  G);

pct = 100 * num ./ den;

out = table(slides, num, den, pct, ...
    'VariableNames', {'slideName','Num_CD103_in_CD4GZMB','Den_CD4GZMB','Pct'});

%% what percentage of GZMB+ cells are also GNLY+?
% Assumes T has a column identifying the slide, e.g. T.Slide or T.SlideName
% Adjust 'Slide' below to match your actual column name.

slides = unique(alldata.slideName);   % or T.SlideName, T.slide_id, etc.

for i = 1:numel(slideName)

    % Filter to this slide
    s = slideName(i);
    if iscell(s)
        mask = strcmp(alldata.slideName, s{1});
        slide_label = s{1};
    else
        mask = alldata.slideName == s;
        slide_label = num2str(s);
    end

    % ---- Diagram 1: CD3+CD8+ subset ----
    subset = mask & alldata.CD3dp == 1 & alldata.CD8ap == 1;
    total  = sum(subset);

    if total == 0
        fprintf('Slide %s: no CD3+CD8+ cells, skipping.\n', slide_label);
    else
        gzmb_only = sum( subset & alldata.GranzymeBp == 1 & alldata.GNLYp == 0 );
        gnly_only = sum( subset & alldata.GranzymeBp == 0 & alldata.GNLYp == 1 );
        both      = sum( subset & alldata.GranzymeBp == 1 & alldata.GNLYp == 1 );

        fprintf('Slide %s | CD3+CD8+ | GZMB only: %d  GNLY only: %d  DP: %d  total: %d\n', ...
            slide_label, gzmb_only, gnly_only, both, total);

        figure('Name', sprintf('Venn CD3+CD8+ — %s', slide_label), ...
               'Color','white','Position',[100 100 520 420]);
        drawVenn(gzmb_only, gnly_only, both, total, ...
            'GZMB+','GNLY+', ...
            sprintf('GZMB vs GNLY | CD3^+CD8^+ | %s', slide_label), ...
            [0.27 0.58 0.80],[0.95 0.60 0.30]);
        fig = gcf;
        exportgraphics(fig, fullfile(sprintf('%s_CD3CD8_venn.png', slide_label)), 'Resolution', 300);
    end
    

    % ---- Diagram 2: All cells on this slide ----
    total2    = sum(mask);
    gzmb_only2 = sum( mask & alldata.GranzymeBp == 1 & alldata.GNLYp == 0 );
    gnly_only2 = sum( mask & alldata.GranzymeBp == 0 & alldata.GNLYp == 1 );
    both2      = sum( mask & alldata.GranzymeBp == 1 & alldata.GNLYp == 1 );

    fprintf('Slide %s | All cells | GZMB only: %d  GNLY only: %d  DP: %d  total: %d\n', ...
        slide_label, gzmb_only2, gnly_only2, both2, total2);

%     figure('Name', sprintf('Venn all cells — %s', slide_label), ...
%            'Color','white','Position',[640 100 520 420]);
%     drawVenn(gzmb_only2, gnly_only2, both2, total2, ...
%         'GZMB+','GNLY+', ...
%         sprintf('GZMB vs GNLY | all cells | %s', slide_label), ...
%         [0.27 0.58 0.80],[0.95 0.60 0.30]);
end

%% CD3+CD8+ subset across all slides 

subset = alldata.CD3dp == 1 & alldata.CD8ap == 1 & alldata.Region==1;
n_total = sum(subset);

gzmb_only = sum(subset & alldata.GranzymeBp == 1 & alldata.GNLYp == 0);
gnly_only = sum(subset & alldata.GranzymeBp == 0 & alldata.GNLYp == 1);
both      = sum(subset & alldata.GranzymeBp == 1 & alldata.GNLYp == 1);
%neither   = n_total - gzmb_only - gnly_only - both;

props = [gzmb_only, gnly_only, both] / n_total * 100;
labels = {'GZMB+ only','GNLY+ only','GZMB^+ & GNLY^+'};
colors = [0.27 0.58 0.80; 0.95 0.60 0.30; 0.45 0.72 0.45];

figure('Color','white','Position',[100 100 500 420]);
b = bar(props, 'FaceColor','flat');
for k = 1:3
    b.CData(k,:) = colors(k,:);
end
xticklabels(labels)
xtickangle(30)
ylabel('% of CD3^+CD8^+ cells')
title(sprintf('GZMB / GNLY co-expression — all cells (n=%d)', n_total))
box off
ylim([0 25])

% Add value labels on top of each bar
for k = 1:4
    text(k, props(k) + 1.5, sprintf('%.1f%%', props(k)), ...
        'HorizontalAlignment','center','FontSize',10);
end

%exportgraphics(gcf, fullfile(outDir, 'coexpression_merged.png'), 'Resolution',300);

%% GZMB/GNLY positivity by MMR status
statuses = unique(alldata.MMR_status);   % replace 'Status' with your actual column name

figure('Color','white','Position',[100 100 900 420]);
colors = [0.27 0.58 0.80; 0.95 0.60 0.30; 0.45 0.72 0.45; 0.88 0.88 0.88];

for i = 1:numel(statuses)
    st = statuses{i};
    mask   = strcmp(alldata.MMR_status, st);
    subset = mask & alldata.CD3dp == 1 & alldata.CD8ap == 1 & alldata.Region == 1;
    n_total = sum(subset);

    gzmb_only = sum(subset & alldata.GranzymeBp == 1 & alldata.GNLYp == 0);
    gnly_only = sum(subset & alldata.GranzymeBp == 0 & alldata.GNLYp == 1);
    both      = sum(subset & alldata.GranzymeBp == 1 & alldata.GNLYp == 1);
    %neither   = n_total - gzmb_only - gnly_only - both;

    props = [gzmb_only, gnly_only, both] / n_total * 100;

    subplot(1, numel(statuses), i)
    b = bar(props, 'FaceColor','flat');
    for k = 1:3
        b.CData(k,:) = colors(k,:);
    end
    xticklabels({'GZMB^+ only','GNLY^+ only','GZMB^+ & GNLY^+',})
    xtickangle(30)
    ylabel('% of CD3^+CD8^+ cells')
    title(sprintf('%s (n=%d)', st, n_total))
    ylim([0 30])
    box off

    for k = 1:3
        text(k, props(k) + 1.5, sprintf('%.1f%%', props(k)), ...
            'HorizontalAlignment','center','FontSize',9);
    end
end

%% GZMB/GNLY co-expression
statuses = unique(alldata.MMR_status);
n_groups = numel(statuses);

gzmb_coex = zeros(1, n_groups);   % % of GZMB+ also GNLY+
gnly_coex = zeros(1, n_groups);   % % of GNLY+ also GZMB+
ns        = zeros(1, n_groups);

for i = 1:n_groups
    st     = statuses{i};
    mask   = strcmp(alldata.MMR_status, st);
    subset = mask & alldata.CD3dp == 1 & alldata.CD8ap == 1 & alldata.Region==1;

    gzmb_only = sum(subset & alldata.GranzymeBp == 1 & alldata.GNLYp == 0);
    gnly_only = sum(subset & alldata.GranzymeBp == 0 & alldata.GNLYp == 1);
    both      = sum(subset & alldata.GranzymeBp == 1 & alldata.GNLYp == 1);

    gzmb_coex(i) = 100 * both / max(both + gzmb_only, 1);
    gnly_coex(i) = 100 * both / max(both + gnly_only, 1);
    ns(i)        = sum(subset);
end

status_labels = cellfun(@(s,n) sprintf('%s\n(n=%d)', s, n), ...
    statuses, num2cell(ns(:)), 'UniformOutput', false);

figure('Color','white','Position',[100 100 600 420]);
hold on

x = 1:n_groups;
b1 = bar(x - 0.2, gzmb_coex, 0.35, 'FaceColor',[0.27 0.58 0.80]);
b2 = bar(x + 0.2, gnly_coex, 0.35, 'FaceColor',[0.95 0.60 0.30]);

% Value labels
for i = 1:n_groups
    text(i - 0.2, gzmb_coex(i) + 1.5, sprintf('%.1f%%', gzmb_coex(i)), ...
        'HorizontalAlignment','center','FontSize',9,'Color',[0.27 0.58 0.80]*0.7);
    text(i + 0.2, gnly_coex(i) + 1.5, sprintf('%.1f%%', gnly_coex(i)), ...
        'HorizontalAlignment','center','FontSize',9,'Color',[0.95 0.60 0.30]*0.7);
end

xticks(x)
xticklabels(statuses)
ylabel('Co-expression rate (%)')
title('GNLY / GZMB co-expression in CD3^+CD8^+ cells')
legend({'% of GZMB^+ also GNLY^+','% of GNLY^+ also GZMB^+'}, 'Location','best')
ylim([0 100])
box off
hold off

%exportgraphics(gcf, fullfile(outDir, 'coexpression_rates_by_status.png'), 'Resolution',300);


%% Hierarchical proportion breakdown of T cells and checkpoint markers

PLOT_FIGS = true;
RUN_BY_STATUS = true;
outDir    = 'proportion_output';
if ~exist(outDir,'dir'), mkdir(outDir); end

% Color palette
col_blue   = [0.27 0.58 0.80];
col_orange = [0.95 0.60 0.30];
col_green  = [0.45 0.72 0.45];
col_pink   = [0.85 0.40 0.60];
col_purple = [0.55 0.40 0.80];
col_gray   = [0.70 0.70 0.70];

% Setup: either loop over status groups or run once on all data
if RUN_BY_STATUS
    status_groups = unique(alldata.MMR_status);
else
    status_groups = {'ALL'};   % single dummy group
end

for g = 1:numel(status_groups)

    if RUN_BY_STATUS
        grp_label = status_groups{g};
        grp_mask  = strcmp(alldata.MMR_status, grp_label);
    else
        grp_label = 'all';
        grp_mask  = true(height(alldata), 1);
    end

    % Subset the data for this group
    grp = alldata(grp_mask, :);   % work with grp.CD3, grp.CD8 etc. inside the loop
    grp = grp(grp.Region == 1, :);
    % ARM 1 -------------------------------------------------------
    gzmb_mask  = grp.GranzymeBp == 1;
    n_gzmb     = sum(gzmb_mask);
    
    % CD3+ vs CD3- among GZMB+
    gzmb_cd3pos = sum(gzmb_mask & grp.CD3dp == 1);
    gzmb_cd3neg = sum(gzmb_mask & grp.CD3dp == 0);

    % Among GZMB+CD3+: CD4 vs CD8
    gzmb_cd3_mask = gzmb_mask & grp.CD3dp == 1;
    gzmb_cd4  = sum(gzmb_cd3_mask & grp.CD4p == 1);
    gzmb_cd8  = sum(gzmb_cd3_mask & grp.CD8ap == 1);
    gzmb_dn   = sum(gzmb_cd3_mask & grp.CD4p == 0 & grp.CD8ap == 0);  % double negative

    % Among GZMB+CD8+: PD1 fraction
    gzmb_cd8_mask = gzmb_mask & grp.CD8ap == 1;
    n_gzmb_cd8    = sum(gzmb_cd8_mask);
    gzmb_cd8_pd1pos = sum(gzmb_cd8_mask & grp.PD_1p == 1);
    gzmb_cd8_pd1neg = sum(gzmb_cd8_mask & grp.PD_1p == 0);
    
    % Among GZMB+CD8+: NKG2A fraction
    gzmb_cd8_nkg2apos = sum(gzmb_cd8_mask & grp.NKG2Ap == 1);
    gzmb_cd8_nkg2aneg = sum(gzmb_cd8_mask & grp.NKG2Ap == 0);

    % Among GZMB+CD4+: NKG2A fraction
    gzmb_cd4_mask = gzmb_mask & grp.CD4p == 1;
    n_gzmb_cd4    = sum(gzmb_cd4_mask);
    gzmb_cd4_nkg2apos = sum(gzmb_cd4_mask & grp.NKG2Ap == 1);
    gzmb_cd4_nkg2aneg = sum(gzmb_cd4_mask & grp.NKG2Ap == 0);    
   
    fprintf('=== ARM 1: GZMB+ breakdown ===\n')
    fprintf('GZMB+ total:            %d\n', n_gzmb)
    fprintf('  CD3+:                 %d (%.1f%%)\n', gzmb_cd3pos, 100*gzmb_cd3pos/n_gzmb)
    fprintf('  CD3-:                 %d (%.1f%%)\n', gzmb_cd3neg, 100*gzmb_cd3neg/n_gzmb)
    fprintf('  Among CD3+:\n')
    fprintf('    CD4+:               %d (%.1f%%)\n', gzmb_cd4, 100*gzmb_cd4/gzmb_cd3pos)
    fprintf('    CD8+:               %d (%.1f%%)\n', gzmb_cd8, 100*gzmb_cd8/gzmb_cd3pos)
    fprintf('    DN:                 %d (%.1f%%)\n', gzmb_dn,  100*gzmb_dn/gzmb_cd3pos)
    fprintf('  Among CD8+:\n')
    fprintf('    PD1+:               %d (%.1f%%)\n', gzmb_cd8_pd1pos, 100*gzmb_cd8_pd1pos/max(n_gzmb_cd8,1))
    fprintf('    PD1-:               %d (%.1f%%)\n', gzmb_cd8_pd1neg, 100*gzmb_cd8_pd1neg/max(n_gzmb_cd8,1))
    fprintf('    NKG2A+:               %d (%.1f%%)\n', gzmb_cd8_nkg2apos, 100*gzmb_cd8_nkg2apos/max(n_gzmb_cd8,1))
    fprintf('    NKG2A-:               %d (%.1f%%)\n', gzmb_cd8_nkg2aneg, 100*gzmb_cd8_nkg2aneg/max(n_gzmb_cd8,1))
    fprintf('  Among CD4+:\n')
    fprintf('    NKG2A+:             %d (%.1f%%)\n', gzmb_cd4_nkg2apos, 100*gzmb_cd4_nkg2apos/max(n_gzmb_cd4,1))
    fprintf('    NKG2A-:             %d (%.1f%%)\n', gzmb_cd4_nkg2aneg, 100*gzmb_cd4_nkg2aneg/max(n_gzmb_cd4,1))

    if PLOT_FIGS
    % --- Figure 1a: GZMB+ top-level CD3 split ---
    figure('Color','white','Position',[100 600 500 380]);
    drawPropBar( ...
        [gzmb_cd3pos, gzmb_cd3neg], ...
        {'CD3+','CD3-'}, ...
        [col_blue; col_gray], ...
        sprintf('GZMB^+ cells (n=%d)', n_gzmb), ...
        '% of GZMB^+ cells', outDir, sprintf('%s_arm1_gzmb_cd3split', grp_label));

    % --- Figure 1b: Among GZMB+CD3+: CD4 vs CD8 ---
    figure('Color','white','Position',[620 600 500 380]);
    drawPropBar( ...
        [gzmb_cd4, gzmb_cd8, gzmb_dn], ...
        {'CD4+','CD8+','DN'}, ...
        [col_orange; col_blue; col_gray], ...
        sprintf('GZMB^+CD3^+ cells (n=%d)', gzmb_cd3pos), ...
        '% of GZMB^+CD3^+ cells', outDir, sprintf('%s_arm1_gzmb_cd4cd8split', grp_label));

    % --- Figure 1c: GZMB+CD8+ PD1 split ---
    figure('Color','white','Position',[100 150 500 380]);
    drawPropBar( ...
        [gzmb_cd8_pd1pos, gzmb_cd8_pd1neg], ...
        {'PD1+','PD1-'}, ...
        [col_purple; col_gray], ...
        sprintf('GZMB^+CD8^+ cells (n=%d)', n_gzmb_cd8), ...
        '% of GZMB^+CD8^+ cells', outDir, sprintf('%s_arm1_gzmb_cd8_pd1', grp_label));
    

    % --- Figure 1d: GZMB+CD8+ NKG2A split ---
    figure('Color','white','Position',[620 150 500 380]);
    drawPropBar( ...
        [gzmb_cd8_nkg2apos, gzmb_cd8_nkg2aneg], ...
        {'NKG2A+','NKG2A-'}, ...
        [col_green; col_gray], ...
        sprintf('GZMB^+CD8^+ cells (n=%d)', n_gzmb_cd8), ...
        '% of GZMB^+CD8^+ cells', outDir, sprintf('%s_arm1_gzmb_cd8_nkg2a', grp_label));
    
    % --- Figure 1e: GZMB+CD4+ NKG2A split ---
    figure('Color','white','Position',[620 150 500 380]);
    drawPropBar( ...
        [gzmb_cd4_nkg2apos, gzmb_cd4_nkg2aneg], ...
        {'NKG2A+','NKG2A-'}, ...
        [col_green; col_gray], ...
        sprintf('GZMB^+CD4^+ cells (n=%d)', n_gzmb_cd4), ...
        '% of GZMB^+CD4^+ cells', outDir, sprintf('%s_arm1_gzmb_cd4_nkg2a', grp_label));
    end
    

    % ARM 2:  -------------------------------------------------------
    % CD4 Treg vs Thelper (FOXP3)
    cd4_mask = grp.CD3dp == 1 & grp.CD4p == 1;
    n_cd4      = sum(cd4_mask);

    treg       = sum(cd4_mask & grp.FOXP3p == 1);
    thelper    = sum(cd4_mask & grp.FOXP3p == 0);

    % Among Thelper: GZMB fraction
    thelper_mask  = cd4_mask & grp.FOXP3p == 0;
    th_gzmb       = sum(thelper_mask & grp.GranzymeBp == 1);
    th_nogzmb     = sum(thelper_mask & grp.GranzymeBp == 0);

    % Among Treg: GZMB fraction
    treg_mask     = cd4_mask & grp.FOXP3p == 1;
    treg_gzmb     = sum(treg_mask & grp.GranzymeBp == 1);
    treg_nogzmb   = sum(treg_mask & grp.GranzymeBp == 0);

    fprintf('\n=== ARM 2: CD4+ Treg vs Thelper ===\n')
    fprintf('CD4+ total:             %d\n', n_cd4)
    fprintf('  FOXP3+ (Treg):        %d (%.1f%%)\n', treg,    100*treg/n_cd4)
    fprintf('  FOXP3- (Thelper):     %d (%.1f%%)\n', thelper, 100*thelper/n_cd4)
    fprintf('  Among Thelper:\n')
    fprintf('    GZMB+:              %d (%.1f%%)\n', th_gzmb,   100*th_gzmb/max(thelper,1))
    fprintf('  Among Treg:\n')
    fprintf('    GZMB+:              %d (%.1f%%)\n', treg_gzmb, 100*treg_gzmb/max(treg,1))

    if PLOT_FIGS
        % --- Figure 2a: CD4+ Treg vs Thelper ---
        figure('Color','white','Position',[100 600 500 380]);
        drawPropBar( ...
            [thelper, treg], ...
            {'Thelper (FOXP3-)','Treg (FOXP3+)'}, ...
            [col_orange; col_pink], ...
            sprintf('CD3^+CD4^+ cells (n=%d)', n_cd4), ...
            '% of CD4^+ cells', outDir, sprintf('%s_arm2_cd4_treg_thelper', grp_label));

        % --- Figure 2b: Thelper GZMB ---
        figure('Color','white','Position',[620 600 500 380]);
        drawPropBar( ...
            [th_gzmb, th_nogzmb], ...
            {'GZMB+','GZMB-'}, ...
            [col_blue; col_gray], ...
            sprintf('Thelper cells (n=%d)', thelper), ...
            '% of Thelper cells', outDir, sprintf('%s_arm2_thelper_gzmb', grp_label));
    end

    % ARM 3 -------------------------------------------------------
    % Checkpoint markers — PD1, PDL1, CTLA4
    % Among CD4+: CTLA4 fraction
    cd4_ctla4     = sum(cd4_mask & grp.CTLA4p == 1);
    cd4_noctla4   = sum(cd4_mask & grp.CTLA4p == 0);

    % Among CD8+: PD1 fraction
    cd8_mask      = grp.CD3dp == 1 & grp.CD8ap == 1;
    n_cd8         = sum(cd8_mask);
    cd8_pd1       = sum(cd8_mask & grp.PD_1p == 1);
    cd8_nopd1     = sum(cd8_mask & grp.PD_1p == 0);

    % PD1 / CTLA4 co-expression among CD4+
    pd1_ctla4_both  = sum(cd4_mask & grp.PD_1p == 1 & grp.CTLA4p == 1);
    pd1_only        = sum(cd4_mask & grp.PD_1p == 1 & grp.CTLA4p == 0);
    ctla4_only      = sum(cd4_mask & grp.PD_1p == 0 & grp.CTLA4p == 1);
    neither_ck      = sum(cd4_mask & grp.PD_1p == 0 & grp.CTLA4p == 0);

    fprintf('\n=== ARM 3: Checkpoint markers ===\n')
    fprintf('CD4+ total:             %d\n', n_cd4)
    fprintf('  CTLA4+:               %d (%.1f%%)\n', cd4_ctla4,   100*cd4_ctla4/n_cd4)
    fprintf('  CTLA4-:               %d (%.1f%%)\n', cd4_noctla4, 100*cd4_noctla4/n_cd4)
    fprintf('CD8+ total:             %d\n', n_cd8)
    fprintf('  PD1+:                 %d (%.1f%%)\n', cd8_pd1,     100*cd8_pd1/n_cd8)
    fprintf('  PD1-:                 %d (%.1f%%)\n', cd8_nopd1,   100*cd8_nopd1/n_cd8)
    fprintf('CD4+ PD1/CTLA4 overlap:\n')
    fprintf('  PD1+ only:            %d (%.1f%%)\n', pd1_only,       100*pd1_only/n_cd4)
    fprintf('  CTLA4+ only:          %d (%.1f%%)\n', ctla4_only,     100*ctla4_only/n_cd4)
    fprintf('  PD1+ & CTLA4+:        %d (%.1f%%)\n', pd1_ctla4_both, 100*pd1_ctla4_both/n_cd4)
    fprintf('  Neither:              %d (%.1f%%)\n', neither_ck,     100*neither_ck/n_cd4)

    if PLOT_FIGS
        % --- Figure 3a: CD4+ CTLA4 ---
        figure('Color','white','Position',[100 600 500 380]);
        drawPropBar( ...
            [cd4_ctla4, cd4_noctla4], ...
            {'CTLA4+','CTLA4-'}, ...
            [col_purple; col_gray], ...
            sprintf('CD4^+ cells (n=%d)', n_cd4), ...
            '% of CD4^+ cells', outDir, sprintf('%s_arm3_cd4_ctla4', grp_label));

        % --- Figure 3b: CD8+ PD1 ---
        figure('Color','white','Position',[620 600 500 380]);
        drawPropBar( ...
            [cd8_pd1, cd8_nopd1], ...
            {'PD1+','PD1-'}, ...
            [col_purple; col_gray], ...
            sprintf('CD8^+ cells (n=%d)', n_cd8), ...
            '% of CD8^+ cells', outDir, sprintf('%s_arm3_cd8_pd1', grp_label));

        % --- Figure 3c: PD1 vs CTLA4 Venn in CD4+ ---
        figure('Color','white','Position',[100 150 520 420]);
        drawVenn(pd1_only, ctla4_only, pd1_ctla4_both, n_cd4, ...
            'PD1+','CTLA4+', ...
            'PD1 vs CTLA4 in CD4^+ cells', ...
            col_purple, col_pink);
        exportgraphics(gcf, fullfile(outDir, sprintf('%s_arm3_cd4_pd1_ctla4_venn.png', grp_label)),'Resolution',300);
    end

end

%% Checkpoint markers

PLOT_FIGS = true;
RUN_BY_STATUS = true;
outDir    = 'proportion_output';
if ~exist(outDir,'dir'), mkdir(outDir); end

% Color palette
col_blue   = [0.27 0.58 0.80];
col_orange = [0.95 0.60 0.30];
col_green  = [0.45 0.72 0.45];
col_pink   = [0.85 0.40 0.60];
col_purple = [0.55 0.40 0.80];
col_gray   = [0.70 0.70 0.70];

% Setup: either loop over status groups or run once on all data
if RUN_BY_STATUS
    status_groups = unique(alldata.MMR_status);
else
    status_groups = {'ALL'};   % single dummy group
end

for g = 1:numel(status_groups)

    if RUN_BY_STATUS
        grp_label = status_groups{g};
        grp_mask  = strcmp(alldata.MMR_status, grp_label);
    else
        grp_label = 'all';
        grp_mask  = true(height(alldata), 1);
    end

    % Subset the data for this group
    grp = alldata(grp_mask, :);   % work with grp.CD3, grp.CD8 etc. inside the loop
    grp = grp(grp.Region == 1, :);
    
    % Checkpoint markers — PD1, PDL1, CTLA4
    % Among CD8+: PD1 fraction
    cd8_mask      = grp.CD3dp == 1 & grp.CD8ap == 1;
    n_cd8         = sum(cd8_mask);
    cd8_pd1       = sum(cd8_mask & grp.PD_1p == 1);
    cd8_nopd1     = sum(cd8_mask & grp.PD_1p == 0);
    
    % Among CD8+: CTLA4 fraction
    cd8_ctla4     = sum(cd8_mask & grp.CTLA4p == 1);
    cd8_noctla4   = sum(cd8_mask & grp.CTLA4p == 0);

    % PD1 / CTLA4 co-expression among CD8+
    pd1_ctla4_both  = sum(cd8_mask & grp.PD_1p == 1 & grp.CTLA4p == 1);
    pd1_only        = sum(cd8_mask & grp.PD_1p == 1 & grp.CTLA4p == 0);
    ctla4_only      = sum(cd8_mask & grp.PD_1p == 0 & grp.CTLA4p == 1);
    neither_ck      = sum(cd8_mask & grp.PD_1p == 0 & grp.CTLA4p == 0);
    
    % LAG3 on CD8
    cd8_lag3 = sum(cd8_mask & grp.LAG3p == 1);
    cd8_nolag3 = sum(cd8_mask & grp.LAG3p == 0);
    
    % PD1 / LAG3 co-expression among CD8+
    pd1_lag3_both  = sum(cd8_mask & grp.PD_1p == 1 & grp.LAG3p == 1);
    pd1_only        = sum(cd8_mask & grp.PD_1p == 1 & grp.LAG3p == 0);
    lag3_only      = sum(cd8_mask & grp.PD_1p == 0 & grp.LAG3p == 1);
    neither_ck      = sum(cd8_mask & grp.PD_1p == 0 & grp.LAG3p == 0);

    fprintf('\n=== ARM 3: Checkpoint markers ===\n')
    fprintf('CD8+ total:             %d\n', n_cd8)
    fprintf('  CTLA4+:               %d (%.1f%%)\n', cd8_ctla4,   100*cd8_ctla4/n_cd8)
    fprintf('  CTLA4-:               %d (%.1f%%)\n', cd8_noctla4, 100*cd8_noctla4/n_cd8)
    fprintf('CD8+ total:             %d\n', n_cd8)
    fprintf('  PD1+:                 %d (%.1f%%)\n', cd8_pd1,     100*cd8_pd1/n_cd8)
    fprintf('  PD1-:                 %d (%.1f%%)\n', cd8_nopd1,   100*cd8_nopd1/n_cd8)
    fprintf('CD8+ total:             %d\n', n_cd8)
    fprintf('  LAG3+:                 %d (%.1f%%)\n', cd8_lag3,     100*cd8_lag3/n_cd8)
    fprintf('  LAG3-:                 %d (%.1f%%)\n', cd8_nolag3,   100*cd8_nolag3/n_cd8)
    fprintf('CD8+ PD1/CTLA4 overlap:\n')
    fprintf('  PD1+ only:            %d (%.1f%%)\n', pd1_only,       100*pd1_only/n_cd8)
    fprintf('  CTLA4+ only:          %d (%.1f%%)\n', ctla4_only,     100*ctla4_only/n_cd8)
    fprintf('  PD1+ & CTLA4+:        %d (%.1f%%)\n', pd1_ctla4_both, 100*pd1_ctla4_both/n_cd8)
    fprintf('  Neither:              %d (%.1f%%)\n', neither_ck,     100*neither_ck/n_cd8)
    fprintf('CD8+ PD1/LAG3 overlap:\n')
    fprintf('  PD1+ only:            %d (%.1f%%)\n', pd1_only,       100*pd1_only/n_cd8)
    fprintf('  LAG3+ only:          %d (%.1f%%)\n', lag3_only,     100*lag3_only/n_cd8)
    fprintf('  PD1+ & LAG3+:        %d (%.1f%%)\n', pd1_lag3_both, 100*pd1_lag3_both/n_cd8)
    fprintf('  Neither:              %d (%.1f%%)\n', neither_ck,     100*neither_ck/n_cd8)

    if PLOT_FIGS
        % --- Figure 3a: CD8+ CTLA4 ---
        figure('Color','white','Position',[100 600 500 380]);
        drawPropBar( ...
            [cd8_ctla4, cd8_noctla4], ...
            {'CTLA4+','CTLA4-'}, ...
            [col_purple; col_gray], ...
            sprintf('CD8^+ cells (n=%d)', n_cd8), ...
            '% of CD8^+ cells', outDir, sprintf('%s_arm3_cd8_ctla4', grp_label));

        % --- Figure 3b: CD8+ PD1 ---
        figure('Color','white','Position',[620 600 500 380]);
        drawPropBar( ...
            [cd8_pd1, cd8_nopd1], ...
            {'PD1+','PD1-'}, ...
            [col_purple; col_gray], ...
            sprintf('CD8^+ cells (n=%d)', n_cd8), ...
            '% of CD8^+ cells', outDir, sprintf('%s_arm3_cd8_pd1', grp_label));

        % --- Figure 3c: PD1 vs CTLA4 Venn in CD8+ ---
        figure('Color','white','Position',[100 150 520 420]);
        drawVenn(pd1_only, ctla4_only, pd1_ctla4_both, n_cd8, ...
            'PD1+','CTLA4+', ...
            'PD1 vs CTLA4 in CD8^+ cells', ...
            col_purple, col_pink);
        exportgraphics(gcf, fullfile(outDir, sprintf('%s_arm3_cd8_pd1_ctla4_venn.png', grp_label)),'Resolution',300);
        
        % --- Figure 3d: PD1 vs LAG3 Venn in CD8+ ---
        figure('Color','white','Position',[100 150 520 420]);
        drawVenn(pd1_only, lag3_only, pd1_lag3_both, n_cd8, ...
            'PD1+','LAG3+', ...
            'PD1 vs LAG3 in CD8^+ cells', ...
            col_purple, col_pink);
        exportgraphics(gcf, fullfile(outDir, sprintf('%s_arm3_cd8_pd1_lag3_venn.png', grp_label)),'Resolution',300);
        
        % --- Figure 3e: CD8+ LAG3 ---
        figure('Color','white','Position',[620 600 500 380]);
        drawPropBar( ...
            [cd8_lag3, cd8_nolag3], ...
            {'LAG3+','LAG3-'}, ...
            [col_purple; col_gray], ...
            sprintf('CD8^+ cells (n=%d)', n_cd8), ...
            '% of CD8^+ cells', outDir, sprintf('%s_arm3_cd8_lag3', grp_label));
    end

end

%% Grouped by MMR status
% Preallocate — rows = status groups, cols = categories
props_gzmb_cd3 = zeros(numel(status_groups), 2);   % CD3+ | CD3-
props_cd4_split = zeros(numel(status_groups), 3);  % CD4+ | CD8+ | DN
props_cd8_pd1   = zeros(numel(status_groups), 2);  % PD1+ | PD1-
% etc. for each arm
props_cd4_th_treg = zeros(numel(status_groups),2); % Th | Treg
props_th_gzmb = zeros(numel(status_groups),2); % GZMB+ | GZMB - 
props_treg_gzmb = zeros(numel(status_groups),2); % GZMB + | GZMB -

props_ctla4 = zeros(numel(status_groups), 2); % CD4+ CTLA4
props_cd8pd1 = zeros(numel(status_groups), 2); %CD8 PD1+ PD1-
props_cd4 = zeros(numel(status_groups), 4); %CD4 PD1+ CTLA4+


for g = 1:numel(status_groups)
    grp_label = status_groups{g};
    grp = alldata(strcmp(alldata.MMR_status, grp_label), :);

    % ARM 1
    gzmb_mask   = grp.GranzymeBp == 1;
    n_gzmb      = sum(gzmb_mask);
    gzmb_cd3pos = sum(gzmb_mask & grp.CD3dp == 1);
    gzmb_cd3neg = sum(gzmb_mask & grp.CD3dp == 0);
    props_gzmb_cd3(g,:) = [gzmb_cd3pos, gzmb_cd3neg] / max(n_gzmb,1) * 100;

    gzmb_cd3_mask = gzmb_mask & grp.CD3dp == 1;
    gzmb_cd4 = sum(gzmb_cd3_mask & grp.CD4p == 1);
    gzmb_cd8 = sum(gzmb_cd3_mask & grp.CD8ap == 1);
    gzmb_dn  = sum(gzmb_cd3_mask & grp.CD4p == 0 & grp.CD8ap == 0);
    props_cd4_split(g,:) = [gzmb_cd4, gzmb_cd8, gzmb_dn] / max(sum(gzmb_cd3_mask),1) * 100;

    gzmb_cd8_mask   = gzmb_mask & grp.CD8ap == 1;
    n_gzmb_cd8      = sum(gzmb_cd8_mask);
    gzmb_cd8_pd1pos = sum(gzmb_cd8_mask & grp.PD_1p == 1);
    gzmb_cd8_pd1neg = sum(gzmb_cd8_mask & grp.PD_1p == 0);
    props_cd8_pd1(g,:) = [gzmb_cd8_pd1pos, gzmb_cd8_pd1neg] / max(n_gzmb_cd8,1) * 100;
    
    % arm 2
    % CD4 Treg vs Thelper (FOXP3)
    cd4_mask = grp.CD3dp == 1 & grp.CD4p == 1;
    n_cd4      = sum(cd4_mask);

    treg       = sum(cd4_mask & grp.FOXP3p == 1);
    thelper    = sum(cd4_mask & grp.FOXP3p == 0);
    
    props_cd4_th_treg(g,:) = [treg, thelper]/max(n_cd4,1) *100;

    % Among Thelper: GZMB fraction
    thelper_mask  = cd4_mask & grp.FOXP3p == 0;
    n_thelper = sum(thelper_mask);
    th_gzmb       = sum(thelper_mask & grp.GranzymeBp == 1);
    th_nogzmb     = sum(thelper_mask & grp.GranzymeBp == 0);
    
    props_th_gzmb(g,:) = [th_gzmb, th_nogzmb]/max(n_thelper,1) *100;

    % Among Treg: GZMB fraction
    treg_mask     = cd4_mask & grp.FOXP3p == 1;
    n_treg = sum(treg_mask);
    treg_gzmb     = sum(treg_mask & grp.GranzymeBp == 1);
    treg_nogzmb   = sum(treg_mask & grp.GranzymeBp == 0);
    props_treg_gzmb(g,:) = [treg_gzmb, treg_nogzmb]/max(n_treg,1) *100;
    
    %arm 3
    % Among CD4+: CTLA4 fraction
    cd4_ctla4     = sum(cd4_mask & grp.CTLA4p == 1);
    cd4_noctla4   = sum(cd4_mask & grp.CTLA4p == 0);
    props_ctla4(g,:) = [cd4_ctla4, cd4_noctla4]/max(n_cd4,1) *100;

    % Among CD8+: PD1 fraction
    cd8_mask      = grp.CD3dp == 1 & grp.CD8ap == 1;
    n_cd8         = sum(cd8_mask);
    cd8_pd1       = sum(cd8_mask & grp.PD_1p == 1);
    cd8_nopd1     = sum(cd8_mask & grp.PD_1p == 0);
    props_cd8pd1(g,:) = [cd8_pd1, cd8_nopd1]/max(n_cd8,1) *100;

    % PD1 / CTLA4 co-expression among CD4+
    pd1_ctla4_both  = sum(cd4_mask & grp.PD_1p == 1 & grp.CTLA4p == 1);
    pd1_only        = sum(cd4_mask & grp.PD_1p == 1 & grp.CTLA4p == 0);
    ctla4_only      = sum(cd4_mask & grp.PD_1p == 0 & grp.CTLA4p == 1);
    neither_ck      = sum(cd4_mask & grp.PD_1p == 0 & grp.CTLA4p == 0);
    props_cd4(g,:) = [pd1_ctla4_both, pd1_only, ctla4_only, neither_ck]/max(n_cd4,1) *100;
end

% Now plot all groups together
figure('Color','white','Position',[100 100 1000 350]);
tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');
grp_colors = [col_blue; col_orange; col_green];   % one per status group

% Panel 1: GZMB+ CD3 split
nexttile; hold on;
drawGroupedBar(props_gzmb_cd3, status_groups, {'CD3+','CD3-'}, ...
    grp_colors, 'GZMB^+ — CD3 split');

% Panel 2: GZMB+CD3+ lineage split
nexttile;
drawGroupedBar(props_cd4_split, status_groups, {'CD4+','CD8+','DN'}, ...
    grp_colors, 'GZMB^+CD3^+ — lineage');

% Panel 3: GZMB+CD8+ PD1
nexttile;
drawGroupedBar(props_cd8_pd1, status_groups, {'PD1+','PD1-'}, ...
    grp_colors, 'GZMB^+CD8^+ — PD1');

exportgraphics(gcf, fullfile(outDir, 'arm1_grouped_by_status.eps'), 'Resolution',300);

% arm2
% Now plot all groups together
figure('Color','white','Position',[100 100 1000 350]);
tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');
grp_colors = [col_blue; col_orange; col_green];   % one per status group

% Panel 1: CD4 split
nexttile; hold on;
drawGroupedBar(props_cd4_th_treg, status_groups, {'Th','Treg'}, ...
    grp_colors, 'CD4^+ — lineage');

% Panel 2: Thelper split
nexttile;
drawGroupedBar(props_th_gzmb, status_groups, {'GZMB+','GZMB-'}, ...
    grp_colors, 'Thelper');

% Panel 3: Treg split
nexttile;
drawGroupedBar(props_treg_gzmb, status_groups, {'GZMB+','GZMB-'}, ...
    grp_colors, 'Treg');

exportgraphics(gcf, fullfile(outDir, 'arm2_grouped_by_status.eps'), 'Resolution',300);

% arm3
% Now plot all groups together
figure('Color','white','Position',[100 100 1000 350]);
tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');
grp_colors = [col_blue; col_orange; col_green];   % one per status group

% Panel 1: CD4 CTLA4 fraction
nexttile; hold on;
drawGroupedBar(props_ctla4, status_groups, {'CTLA4+','CTLA4-'}, ...
    grp_colors, 'CD4^+ — CTLA4');

% Panel 2: PD1 in CD8
nexttile;
drawGroupedBar(props_cd8pd1, status_groups, {'PD1+','PD1-'}, ...
    grp_colors, 'CD8^+');

% Panel 3: CD4 split PD1 or CTLA4
nexttile;
drawGroupedBar(props_cd4, status_groups, {'PD1^+CTLA4^+','PD1^+CTLA4-', 'PD1^-CTLA4^+', 'Neither'}, ...
    grp_colors, 'CD4');

exportgraphics(gcf, fullfile(outDir, 'arm3_grouped_by_status.eps'), 'Resolution',300);

%% Mac/DC PDL1 graph
% PDL1+ fraction on DCs and Macrophages per slide, colored by MMR status

% Define populations — adjust marker names to match your columns
panck_mask  = alldata.panCKp == 1;
myeloid_mask = alldata.CD68p  == 1 | alldata.CD11cp == 1 | alldata.HLA_DRB1p == 1;

n_slides = numel(slideName); 
pdl1_panck      = zeros(n_slides, 1);
pdl1_myeloid     = zeros(n_slides, 1);
slide_status = cell(n_slides, 1);

for i = 1:n_slides
    s     = slides{i};
    smask = strcmp(alldata.slideName, s);
    n_panckpdl1  = sum(smask & panck_mask);
    n_myeloid = sum(smask & myeloid_mask);
    pdl1_panck(i)      = sum(smask & panck_mask  & alldata.PD_L1p == 1) / max(n_panckpdl1,  1);
    pdl1_myeloid(i)     = sum(smask & myeloid_mask & alldata.PD_L1p == 1) / max(n_myeloid, 1);
    slide_status{i} = alldata.MMR_status{find(smask, 1)};
end

% Status colors
status_list   = unique(slide_status);
status_colors = [col_blue; col_orange; col_green];

%% PDL1+ fraction — bidirectional DC vs Mac per slide

dc_mask  = alldata.CD11cp == 1 & alldata.HLA_DRB1p == 1;
mac_mask = alldata.CD68p  == 1;

n_slides     = numel(slideName);
pdl1_dc      = zeros(n_slides, 1);
pdl1_mac     = zeros(n_slides, 1);
slide_status = cell(n_slides, 1);

for i = 1:n_slides
    s     = slides{i};
    smask = strcmp(alldata.slideName, s);
    n_dc  = sum(smask & dc_mask);
    n_mac = sum(smask & mac_mask);
    pdl1_dc(i)      = sum(smask & dc_mask  & alldata.PD_L1p == 1) / max(n_dc,  1);
    pdl1_mac(i)     = sum(smask & mac_mask & alldata.PD_L1p == 1) / max(n_mac, 1);
    slide_status{i} = alldata.MMR_status{find(smask, 1)};
end

% Status colors
status_list   = unique(slide_status);
status_colors = [col_blue; col_orange; col_green];

% Sort by DC fraction ascending
[~, sort_idx]  = sort(pdl1_mac, 'ascend');
sorted_dc      = pdl1_dc(sort_idx);
sorted_mac     = pdl1_mac(sort_idx);
sorted_slides  = slides(sort_idx);
sorted_status  = slide_status(sort_idx);

% sort by combined DC + mac fraction
% [~, sort_idx] = sort(pdl1_dc + pdl1_mac, 'ascend');
% sorted_dc      = pdl1_dc(sort_idx);
% sorted_mac     = pdl1_mac(sort_idx);
% sorted_slides  = slides(sort_idx);
% sorted_status  = slide_status(sort_idx);

figure('Color','white','Position',[100 100 700 800]);
hold on;

for i = 1:n_slides
    col = status_colors(strcmp(status_list, sorted_status{i}), :);
    
    % DC — draws rectangle from -dc_val to 0
    dc_val  = sorted_dc(i);
    mac_val = sorted_mac(i);
    
    % left bar (DC)
    patch([-dc_val, 0, 0, -dc_val], ...
          [i-bar_h, i-bar_h, i+bar_h, i+bar_h], ...
          col, 'EdgeColor',[0.2 0.2 0.2], 'LineWidth',0.5);
    
    % right bar (Mac)
    patch([0, mac_val, mac_val, 0], ...
          [i-bar_h, i-bar_h, i+bar_h, i+bar_h], ...
          col, 'EdgeColor',[0.2 0.2 0.2], 'LineWidth',0.5);
end

xline(0, 'k-', 'LineWidth', 1.2);

max_val = max([pdl1_dc; pdl1_mac]) * 1.1;
xlim([-max_val, max_val]);

% Show absolute values on x axis
tick_vals = -0.4:0.1:0.4;
xticks(tick_vals);
xticklabels(arrayfun(@(x) sprintf('%.1f', abs(x)), tick_vals, 'UniformOutput',false));

yticks(1:n_slides);
yticklabels(sorted_slides);

% Direction labels
text(-max_val*0.5, n_slides+1.5, '← DCs (CD11c^+HLADR^+)', ...
    'FontSize',9,'Color',[0.3 0.3 0.3],'HorizontalAlignment','center');
text( max_val*0.5, n_slides+1.5, 'Macrophages (CD68^+) →', ...
    'FontSize',9,'Color',[0.3 0.3 0.3],'HorizontalAlignment','center');

xlabel('PDL1^+ cell fraction');
title('PDL1^+ fraction — DCs vs Macrophages', 'FontSize',11,'FontWeight','bold');

% Legend — same patch approach as original
% for g = 1:numel(status_list)
%     patch(nan, nan, status_colors(g,:), 'EdgeColor','k', ...
%         'LineWidth',0.5, 'DisplayName', status_list{g});
% end
% legend('Location','southeast','FontSize',8,'Box','off');

box off;
hold off;

exportgraphics(gcf, fullfile(outDir,'pdl1_bidirectional_mac.eps'),'Resolution',300);
%% PanCK myeloid dominance bidirectional
[~, sort_idx]  = sort(pdl1_panck, 'ascend');
sorted_myeloid      = pdl1_myeloid(sort_idx);
sorted_panck     = pdl1_panck(sort_idx);
sorted_slides  = slides(sort_idx);
sorted_status  = slide_status(sort_idx);

figure('Color','white','Position',[100 100 700 800]);
hold on;

for i = 1:n_slides
    col = status_colors(strcmp(status_list, sorted_status{i}), :);
    
    % DC — draws rectangle from -dc_val to 0
    myeloid_val  = sorted_myeloid(i);
    panck_val = sorted_panck(i);
    
    % left bar (DC)
    patch([-myeloid_val, 0, 0, -myeloid_val], ...
          [i-bar_h, i-bar_h, i+bar_h, i+bar_h], ...
          col, 'EdgeColor',[0.2 0.2 0.2], 'LineWidth',0.5);
    
    % right bar (Mac)
    patch([0, panck_val, panck_val, 0], ...
          [i-bar_h, i-bar_h, i+bar_h, i+bar_h], ...
          col, 'EdgeColor',[0.2 0.2 0.2], 'LineWidth',0.5);
end

xline(0, 'k-', 'LineWidth', 1.2);

max_val = max([pdl1_myeloid; pdl1_panck]) * 1.1;
xlim([-max_val, max_val]);

% Show absolute values on x axis
tick_vals = -0.4:0.1:0.4;
xticks(tick_vals);
xticklabels(arrayfun(@(x) sprintf('%.1f', abs(x)), tick_vals, 'UniformOutput',false));

yticks(1:n_slides);
yticklabels(sorted_slides);

% Direction labels
text(-max_val*0.5, n_slides+1.5, '← Myeloid', ...
    'FontSize',9,'Color',[0.3 0.3 0.3],'HorizontalAlignment','center');
text( max_val*0.5, n_slides+1.5, 'PanCK^+ →', ...
    'FontSize',9,'Color',[0.3 0.3 0.3],'HorizontalAlignment','center');

xlabel('PDL1^+ cell fraction');
title('PDL1^+ fraction — PanCK vs Myeloid', 'FontSize',11,'FontWeight','bold');

% Legend — same patch approach as original
% for g = 1:numel(status_list)
%     patch(nan, nan, status_colors(g,:), 'EdgeColor','k', ...
%         'LineWidth',0.5, 'DisplayName', status_list{g});
% end
% legend('Location','southeast','FontSize',8,'Box','off');

box off;
hold off;

exportgraphics(gcf, fullfile(outDir,'pdl1_bidirectional_panck_myeloid.eps'),'Resolution',300);

%% subset part of alldata

cols = {'slideName', 'CTLA4p', 'CD8CTLA4p', 'CD4CTLA4p', 'CD8Chkpt', 'CD8CTLA4PD1p', 'CD4CTLA4PD1p', 'Region'};

subsetData = alldata(:,cols);

sumTumor_subset = varfun(@mean,subsetData(subsetData.Region==1,:),'GroupingVariables','slideName');
sumTumor_subset = join(sumTumor_subset,slideInfo,'Keys','slideName');

%% join with slideInfo

sumCD4 = varfun(@mean,sumCD4(sumCD4.Region==1, :),'GroupingVariables','slideName');
sumCD4.slideName = string(sumCD4.slideName);
slideInfo.slideName = string(slideInfo.slideName);
sumCD4 = join(sumCD4,slideInfo,'Keys','slideName');

%% opts
opts = struct;
opts.GroupLookup = mmrMap;          % REQUIRED
opts.Cols = struct('CD8','CD8ap','TIM3','TIM3p','LAG3','LAG3p','PD1','PD_1p','Region','Region'); % adapt names
opts.Thresh = struct('CD8',0,'TIM3',0,'LAG3',0,'PD1',0);  % thresholds if using expression columns
opts.RegionUse = 1;                 % tumor only; set [] for no region filter
opts.DoVenn = true;                 % try venn(); fallback to UpSet automatically
opts.DoGZMB = true;                % set true to include GZMB in UpSet + summaries
opts.OutputDir = 'venn_outputs';
opts.GroupOrder = {'tdpMMR','dMMR','tipMMR'};  % optional

%% Plot marker/immune clusters (ICs)
outDir = 'IC_tumorviews';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

slideNames = fieldnames(dataStruct1);

for i = 1:numel(slideNames)
    slideName = slideNames{i};
    tmp = dataStruct1.(slideName);
   
    tmp.bg = zeros(height(tmp),1);

    sw1 = 2;
    sw2 = 50000;
    c1 = [0.7 0.7 0.7];
    s1 = 3;

    fig = figure('Visible', 'on');
    fig.Units = 'inches';
    fig.Position = [1 1 6 6];
    
    CycIF_tumorview(tmp, 'bg', sw1, sw2, c1, s1);

    hold on;
    inCluster = tmp.IC > 0;
    hCluster = scatter(tmp.Xt(inCluster), tmp.Yt(inCluster), 3, ...
        'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.05);

    hCluster.DisplayName = 'IC';
    % hold on;
    pos = tmp.CD4p ==1;
    pos2 = tmp.CD8ap == 1;
%     hPos = scatter(tmp.Xt(pos), tmp.Yt(pos), 0.2, 'blue', 'filled');
%     hPos.DisplayName = 'CD4+';
    
    % density contour 
    if sum(pos) > 10
        hCD4 = dscatter(tmp.Xt(pos), tmp.Yt(pos), ...
            'PLOTTYPE', 'CONTOUR');
        hCD4.DisplayName = 'CD4+ density';
        hCD4.LineColor = [0.00 0.45 0.74];   % MATLAB default blue
        hCD4.LineWidth = 0.75;
        %hCD4.LineStyle = '--';
    end
    
    if sum(pos2) > 10
        hCD8 = dscatter(tmp.Xt(pos2), tmp.Yt(pos2), ...
            'PLOTTYPE', 'CONTOUR');
        hCD8.DisplayName = 'CD8+ density';
        hCD8.LineColor = [0.85 0.33 0.10];
        hCD8.LineWidth = 0.75;
    end
    
    % region outline
    % Outline the specified region (e.g. Region == 1 for tumor)
    
    regionMask = tmp.Region == 1;
    x = tmp.Xt(regionMask);
    y = tmp.Yt(regionMask);

    if numel(x) > 2
        K = convhull(x, y);
        hTumor = plot(x(K), y(K), 'k--', 'LineWidth', 1);  % black outline
        
        hTumor.DisplayName = 'Tumor boundary';
    end
    
    % --- build legend handles safely ---
    hList = hCluster;
    lList = {'IC'};

    if exist('hCD4','var') && isgraphics(hCD4)
        hList(end+1) = hCD4;
        lList{end+1} = 'CD4+ density';
    end

    if exist('hCD8','var') && isgraphics(hCD8)
        hList(end+1) = hCD8;
        lList{end+1} = 'CD8+ density';
    end
    if ~isempty(hTumor) && isgraphics(hTumor)
        hList(end+1) = hTumor; lList{end+1} = 'Tumor boundary';
    end
    

    lgd = legend(hList, lList, 'Location','southwest', ...
    'Orientation','horizontal');
    lgd.FontSize = 10;
    lgd.ItemTokenSize = [12 8];
    lgd.EdgeColor = [0.3 0.3 0.3];  % subtle gray border
    lgd.Color = 'white';            % opaque background
    lgd.AutoUpdate = 'off';   % lock it
    lgd.Box = 'on';

    
    title(slideName, 'Interpreter', 'none');
    
    hold off;

    % optional
    % maskplot = (tmp.IC > 0) & (tmp.CD4p ==1);
    % dscatter(tmp.Xt(maskplot), tmp.Yt(maskplot), ...
    %     'PLOTTYPE', 'CONTOUR');
    

    outFile = fullfile(outDir, [slideName '_CD4_overlay.png']);
    exportgraphics(fig, outFile, 'Resolution', 300);
    
    close(fig);
end 
%% compare T cells on tumor (panck+) vs in surrounding area of tumor (panck-) but all within annotated tumor region

xLabels = string(sumTumor_subset.slideName);
y1 = sumTumor_subset.mean_panCKnegCD8pGZMBp;
y2 = sumTumor_subset.mean_panCKpCD8pGZMBp;

figure; 
bar([y1,y2]);
xticks(1:length(xLabels));
xticklabels(xLabels);
xtickangle(45);
legend({'panCK-CD8+GZMB+', 'panCK+CD8+GZMB+'}, 'Location', 'northwest');

%% barplot of intraepithelial vs stromal/peritumoral T cells by MMR status
sumTumor3.MMR_status = string(cellstr(sumTumor3.MMR_status));
groups = unique(sumTumor3.MMR_status);
ratioTable = table();

for i = 1:numel(groups)
    groupVal = groups(i);
    thisGroup = sumTumor3(sumTumor3.MMR_status == groupVal, :);
    
    xLabels = string(thisGroup.slideName);
    y1 = thisGroup.mean_panCKnegCD8pTCF1p;
    y2 = thisGroup.mean_CD8apTCF1p;
    
    ratios = y2 ./ y1;
    ratios(y1 == 0) = NaN;
    
    % ratio table
    T = table(thisGroup.slideName, repmat(groupVal, size(thisGroup,1),1), ratios, ...
        'VariableNames', {'slideName', 'MMR_status', 'intra_to_stromal_ratio'});
    
    ratioTable = [ratioTable; T];
    
    figure;
    bar([y1,y2]);
    xticks(1:length(xLabels));
    xticklabels(xLabels);
    xtickangle(45);
    legend({'stromal/peritumoral', 'intraepithelial'}, 'Location', 'northwest');
    
    title(['CD8+ TCF1+ T cells in the tumor compartment: ' groupVal]);
end
%% ratios boxplot
figure;
myboxplot3(ratioTable.intra_to_stromal_ratio, ratioTable.MMR_status);
ylabel('Intraepithelial-to-Stromal Ratio');
title('CD8+ TCF1+ localization in tumor compartment');

%% scatter plot comparing two marker sets
% i.e. CD4-CD8+FOXP3+ vs CD4+CD8+FOXP3+
CD4neg_counts = zeros(length(slideName), 1);
CD4pos_counts = zeros(length(slideName), 1);

for i = 1:length(slideName)
    T = dataStruct1.(slideName{i});
    
    CD4neg_counts(i) = sum(T.CD4negCD8pFOXP3);
    CD4pos_counts(i) = sum(T.CD4pCD8apFOXP3p);
    
end

figure;
scatter(CD4neg_counts, CD4pos_counts, 60, 'filled')
xlabel('CD4-CD8+FOXP3+ cells')
ylabel('CD4+CD8+FOXP3+ cells')
title('Sample-level comparison of FOXP3+ CD8 T cell subsets')
grid on

for i = 1:length(slideName)
    text(CD4neg_counts(i), CD4pos_counts(i), slideName{i}, ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8); 
end 


[r,p] = corr(CD4neg_counts, CD4pos_counts, 'type', 'Pearson');
text(max(CD4neg_counts)*0.6, max(CD4pos_counts)*0.9, sprintf('r = %.2f, p = %.3g', r, p))

%% top tipMMR samples expressing a cell type 

tumor = sumTumorsubset;
% Step 1: Multiply mean_CD4negCD8+FOXP3+ by 100 to get proportions
tumor.Proportion_CD4negCD8pFOXP3 = tumor.mean_CD4negCD8pFOXP3 * 100;

% Step 2: Sort data by cell type proportion
sortedData = sortrows(tumor, 'Proportion_CD4negCD8pFOXP3', 'descend');

% Step 3: Generate categorical labels combining IO and mutation status
sortedData.Labels = strcat(sortedData.MMR_status); 

% Step 4: Create the ranked bar plot
figure;
bar(sortedData.Proportion_CD4negCD8pFOXP3, 'FaceColor', 'b'); % Blue bars

% Step 5: Formatting
xticks(1:height(sortedData)); % Set x-ticks for each slide
xticklabels(sortedData.Labels); % Use IO - mutstatus as labels
xtickangle(45); % Rotate labels for readability
ylabel('CD4negCD8pFOXP3 Proportion (%)');
xlabel('Samples (Mutation Status)');
title('Ranked CD4negCD8pFOXP3 Proportion by Mutation Status');
grid on;

% Step 6: Adjust figure for better readability
set(gca, 'FontSize', 10); % Adjust font size for readability

%% Kaplan Meier curve

threshold = median(FOXP3_counts.CD4posCD8pFOXP3p);
group = FOXP3_counts.CD4posCD8pFOXP3p > threshold; % 1 = high, 0 = low
time = FOXP3_counts.PFSDays;
censoring = ~FOXP3_counts.PFSCensor; % 1 event occurred, 0 = censored

keep = ["tipMMR", "dMMR"];
idx = ismember(FOXP3_counts.MMR_status, keep);
time_sub = time(idx);
censor_sub = censoring(idx);
group_sub = FOXP3_counts.MMR_status(idx);

group_label = cell(height(FOXP3_counts), 1);
group_label(FOXP3_counts.CD4posCD8pFOXP3p > median(FOXP3_counts.CD4posCD8pFOXP3p)) = {'Above Median'};
group_label(FOXP3_counts.CD4posCD8pFOXP3p <= median(FOXP3_counts.CD4posCD8pFOXP3p)) = {'Below Median'};

% Call MatSurv to generate the KM plot
[p, fh, stats] = MatSurv(time_sub, censor_sub, group_sub, ...
                         'DispP', true, 'DispHR', true, 'NoRiskTable', true, ...
                         'XLabel', 'Time (PFS)', 'YLabel', 'Survival Probability', ...
                         'Title', 'Kaplan-Meier: CD4+CD8+FOXP3+');

% Save the figure
%saveas(fh, 'all_Ki67p_MatSurv.png');

%% stacked bar chart 

% --- Define Phenotype Rules (mixed AND/OR logic) ---
phenotypeRules = {
    'Cytotoxic T cell',          struct('AND', {{'CD8ap', 'CD3dp'}}, 'OR', {{}});
    'Regulatory T cell',         struct('AND', {{'CD4p', 'CD3dp', 'FOXP3p'}}, 'OR', {{}});

    'Tpex (Progenitor Exhausted)', struct('AND', {{'CD8ap', 'CD3dp', 'PD_1p', 'TCF1p'}}, 'OR', {{'SLAMF6p'}});
    'Exhausted T cell',          struct('AND', {{'CD8ap', 'CD3dp'}}, 'OR', {{'PD_1p', 'LAG3p', 'TIM3p'}});
    'Active T cell',             struct('AND', {{'CD3dp'}}, 'OR', {{'CD134p', 'ICOSp'}});
    'Checkpoint+ T cell',        struct('AND', {{'CD3dp'}}, 'OR', {{'PD_1p', 'LAG3p', 'TIM3p', 'CTLA4p'}});

    'Macrophage',                struct('AND', {{'CD68p', 'CD206p'}}, 'OR', {{}});
    'Dendritic cell',            struct('AND', {{'CD11cp', 'HLA_DRB1p'}}, 'OR', {{}});
    'B cell',                    struct('AND', {{'CD20p'}}, 'OR', {{}});
    %'Other Immune',              struct('AND', {{'CD45p'}}, 'OR', {{}});
    
   % 'Proliferating Tumor',       struct('AND', {{'panCKp', 'Ki_67_2p'}}, 'OR', {{}});
   % 'Non-Proliferating Tumor',   struct('AND', {{'panCKp'}}, 'OR', {{}});
    'Other',                     struct('AND', {{}}, 'OR', {{}});
};



phenotypeNames = phenotypeRules(:, 1);
phenotypeStructs = phenotypeRules(:, 2);
numPhenotypes = size(phenotypeRules, 1);

% --- Set colors ---
%colors = lines(numPhenotypes);
colors = [
    0.0000, 0.4470, 0.7410;  % blue
    0.8500, 0.2050, 0.6980;  % orange
    0.9290, 0.6940, 0.1250;  % yellow
    0.5940, 0.1240, 0.6560;  % purple
    0.4660, 0.6740, 0.1880;  % green
    0.3010, 0.7450, 0.9330;  % light blue
    0.6350, 0.0780, 0.1840;  % dark red
    0.6000, 0.6000, 0.6000;  % gray
    0.9800, 0.5000, 0.4500;  % coral
    0.0000, 0.7490, 0.7490;  % cyan
    0.5000, 0.5000, 0.0000;  % olive
    0.8700, 0.4900, 0.0000;  % brownish-orange
    0.3600, 0.2000, 0.5000;  % deep purple
];


% --- Define MMR groups ---
MMR_groups = {'tdpMMR', 'dMMR', 'tipMMR'};
MMR_agg_counts = zeros(length(MMR_groups), numPhenotypes);

% --- Loop through samples ---
tableNames = strcat('dataStruct1.', slideName);
for t = 1:length(tableNames)
    currentTable = eval(tableNames{t});
    currentName = slideName{t};
    
    % Lookup MMR status from slideInfo
    matchIdx = strcmp(slideInfo.slideName, currentName);
    if sum(matchIdx) ~= 1
        warning(['Could not uniquely match slideName: ', currentName]);
        continue;
    end
    MMR = slideInfo.MMR_status{matchIdx};
    
    mmrIdx = find(strcmp(MMR_groups, MMR));
    if isempty(mmrIdx)
        warning(['Unknown MMR_status: ', MMR]);
        continue;
    end
    
    % subset to immune cells IN THE TUMOR
    subsetMask = (currentTable.Region==1) & (currentTable.CD45p | currentTable.CD3dp | ... 
        currentTable.CD8ap | currentTable.CD45ROp | currentTable.CD4p | ...
        currentTable.CD20p | currentTable.CD68p | currentTable.CD11cp);
    
    currentTable = currentTable(subsetMask, :); 

    % Count cells per phenotype
    % Initialize counts for this sample
    phenotypeCounts = zeros(1, numPhenotypes);
    
    for j = 1:numPhenotypes
        rule = phenotypeStructs{j};

        % --- AND logic ---
        if isfield(rule, 'AND') && ~isempty(rule.AND)
            andMask = true(height(currentTable), 1);
            for m = 1:length(rule.AND)
                marker = rule.AND{m};
                if ismember(marker, currentTable.Properties.VariableNames)
                    andMask = andMask & currentTable{:, marker} == 1;
                else
                    warning('Missing marker in AND: %s', marker);
                    andMask = false(height(currentTable), 1);
                end
            end
        else
            andMask = true(height(currentTable), 1);  % no AND condition
        end

        % --- OR logic ---
        if isfield(rule, 'OR') && ~isempty(rule.OR)
            orMask = false(height(currentTable), 1);
            for m = 1:length(rule.OR)
                marker = rule.OR{m};
                if ismember(marker, currentTable.Properties.VariableNames)
                    orMask = orMask | currentTable{:, marker} == 1;
                else
                    warning('Missing marker in OR: %s', marker);
                end
            end
        else
            orMask = true(height(currentTable), 1);  % no OR condition
        end

        % --- Final mask ---
        finalMask = andMask & orMask;

        % --- Count and store ---
        phenotypeCounts(j) = sum(finalMask);
    end


    % Add to aggregate
    MMR_agg_counts(mmrIdx, :) = MMR_agg_counts(mmrIdx, :) + phenotypeCounts;
end

t_and_b_indices = ismember(phenotypeNames, {
    'Cytotoxic T cell', 'Regulatory T cell', ...
    'Tpex (Progenitor Exhausted)', 'Exhausted T cell', ...
    'Active T cell', 'Checkpoint+ T cell', ...
    'B cell'});

tb_myeloid_indices = ismember(phenotypeNames, {
    'Cytotoxic T cell', 'Regulatory T cell', ...
    'Tpex (Progenitor Exhausted)', 'Exhausted T cell', ...
    'Active T cell', 'Checkpoint+ T cell', ...
    'B cell', ...
    'Macrophage', 'Dendritic cell'
});

% --- Normalize to % per group ---
MMR_agg_percent_tumor = MMR_agg_counts ./ sum(MMR_agg_counts, 2) * 100;

% Normalize only by T + B cell total
T_B_sums = sum(MMR_agg_counts(:, t_and_b_indices), 2);
MMR_agg_percent_t_b = MMR_agg_counts ./ T_B_sums * 100;

% Normalize to T+B+myeloid
tb_myeloid_sums = sum(MMR_agg_counts(:, tb_myeloid_indices), 2);
MMR_tbMyeloid_percent = MMR_agg_counts ./ tb_myeloid_sums * 100;

%% --- Plot stacked bar chart ---
figure;
%h = bar(MMR_agg_percent_tumor, 'stacked'); % percent
% absolute numbers
h = bar(MMR_agg_counts, 'stacked');
%colormap(colors);
for i = 1:numPhenotypes
    h(i).FaceColor = colors(i, :);
end
xticks(1:3);
xticklabels(MMR_groups);
ylabel('Number of cells');
legend(phenotypeNames, 'Interpreter', 'none', 'Location', 'northeastoutside');
title('Immune cell absolute counts in the tumor by MMR Status');
set(gcf, 'Color', 'w');

saveas(gcf, 'phenotype_tumor_immune_counts_stackedbar.pdf');

%% Plot 

figure;
%MMR_agg_percent_plot = MMR_agg_percent(:, t_and_b_indices);  % only T+B
%MMR_tbMyeloid_percent_plot = MMR_tbMyeloid_percent(:, tb_myeloid_indices);
MMR_agg_percent_tb_tumor = MMR_agg_percent_t_b(:, t_and_b_indices);
MMR_tbMyeloid_percent_tumor = MMR_tbMyeloid_percent(:, tb_myeloid_indices);

% absolute counts
MMR_tbMyeloid_counts = MMR_agg_counts(:, tb_myeloid_indices);
MMR_tb_counts = MMR_agg_counts(:, t_and_b_indices);

h = bar(MMR_tbMyeloid_counts, 'stacked');
color_idx = find(tb_myeloid_indices);  % only T+B colors

for i = 1:length(color_idx)
    h(i).FaceColor = colors(color_idx(i), :);
end

xticks(1:3);
xticklabels(MMR_groups);
ylabel('Number of Cells');
legend(phenotypeNames(tb_myeloid_indices), 'Interpreter', 'none', 'Location', 'northeastoutside');
title('T + B + Myeloid Phenotype Abundance in Tumor by MMR Status');
set(gcf, 'Color', 'w');

saveas(gcf, 'phenotype_tbmyeloid_tumor_counts_stackedbar.pdf');