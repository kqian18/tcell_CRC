function drawPropBar(counts, labels, colors, titleStr, ylabelStr, outDir, fname)
    n_total = sum(counts);
    props   = counts / n_total * 100;
    n_bars  = numel(counts);

    ax = gca; cla(ax); hold(ax,'on');
    for k = 1:n_bars
        b = bar(ax, k, props(k), 'FaceColor', colors(k,:), 'EdgeColor','none', 'BarWidth',0.6);
        text(ax, k, props(k) + 1.5, sprintf('%.1f%%', props(k)), ...
            'HorizontalAlignment','center','FontSize',10,'FontWeight','bold', ...
            'Color', colors(k,:) * 0.65);
    end
    xticks(ax, 1:n_bars);
    xticklabels(ax, labels);
    xtickangle(ax, 30);
    ylabel(ax, ylabelStr);
    title(ax, titleStr, 'FontSize',12,'FontWeight','bold');
    ylim(ax, [0 115]);
    box(ax,'off');
    hold(ax,'off');
    exportgraphics(gcf, fullfile(outDir, [fname '.png']), 'Resolution',300);
end