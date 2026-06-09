function drawGroupedBar(props, group_labels, cat_labels, grp_colors, titleStr)
    % props: n_groups x n_cats matrix of percentages
    % Plots grouped bars — one cluster per category, one bar per status group
    n_groups = size(props, 1);
    n_cats   = size(props, 2);
    x        = 1:n_cats;
    w        = 0.8 / n_groups;

    hold on;
    for g = 1:n_groups
        offset = (g - (n_groups+1)/2) * w;
        for k = 1:n_cats
            bar(x(k) + offset, props(g,k), w*0.9, ...
                'FaceColor', grp_colors(g,:), 'EdgeColor','none');
        end
    end

    % Legend entries (one per group)
    for g = 1:n_groups
        bar(nan, nan, 'FaceColor', grp_colors(g,:), 'EdgeColor','none', ...
            'DisplayName', group_labels{g});
    end
    legend('Location','northeast','FontSize',8);

    xticks(x);
    xticklabels(cat_labels);
    ylabel('% of parent population');
    title(titleStr, 'FontSize',11,'FontWeight','bold');
    ylim([0 115]);
    box off;
    hold off;
end