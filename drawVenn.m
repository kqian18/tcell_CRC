function drawVenn(n_left, n_right, n_both, n_total, ...
                  label_left, label_right, titleStr, col_left, col_right)

    ax = gca; cla(ax); axis(ax,'off'); hold(ax,'on');
    ax.Position = [0.05 0.08 0.90 0.82];

    % Circle geometry
    r   = 0.28;
    cx1 = 0.38;   cx2 = 0.62;   cy = 0.50;
    overlap_frac = n_both / max(n_left + n_right + n_both, 1);

    % Draw circles with transparency
    th = linspace(0, 2*pi, 300);
    fill(ax, cx1 + r*cos(th), cy + r*sin(th), col_left,  ...
        'FaceAlpha', 0.35, 'EdgeColor', col_left,  'LineWidth', 1.8);
    fill(ax, cx2 + r*cos(th), cy + r*sin(th), col_right, ...
        'FaceAlpha', 0.35, 'EdgeColor', col_right, 'LineWidth', 1.8);

    % Percentages
    pct = @(n) sprintf('%d\n(%.1f%%)', n, 100*n/n_total);

    % Left-only label
    text(ax, cx1 - 0.13, cy, pct(n_left), ...
        'HorizontalAlignment','center','FontSize',11,'FontWeight','bold', ...
        'Color', col_left * 0.6);

    % Overlap label
    text(ax, (cx1+cx2)/2, cy, pct(n_both), ...
        'HorizontalAlignment','center','FontSize',11,'FontWeight','bold', ...
        'Color', [0.25 0.25 0.25]);

    % Right-only label
    text(ax, cx2 + 0.13, cy, pct(n_right), ...
        'HorizontalAlignment','center','FontSize',11,'FontWeight','bold', ...
        'Color', col_right * 0.6);

    % Circle labels (above circles)
    text(ax, cx1 - 0.10, cy + r + 0.06, label_left, ...
        'HorizontalAlignment','center','FontSize',12,'FontWeight','bold', ...
        'Color', col_left * 0.6);
    text(ax, cx2 + 0.10, cy + r + 0.06, label_right, ...
        'HorizontalAlignment','center','FontSize',12,'FontWeight','bold', ...
        'Color', col_right * 0.6);

    % Title
    title(ax, titleStr, 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.15 0.15 0.15]);

    % Footer: total n
    text(ax, 0.50, 0.04, sprintf('n = %d total cells shown', n_total), ...
        'HorizontalAlignment','center','FontSize',9,'Color',[0.55 0.55 0.55], ...
        'Units','normalized');

    % Overlap annotation (percentage of GZMB+ that are also GNLY+ and other way around)
    if n_left + n_both > 0
    coex1 = 100 * n_both / (n_left + n_both);
    coex2 = 100 * n_both / (n_right + n_both);
    text(ax, 0.50, 0.10, ...
        sprintf('%.1f%% of GZMB^+ are GNLY^+  |  %.1f%% of GNLY^+ are GZMB^+', coex1, coex2), ...
        'HorizontalAlignment','center','FontSize',9,'Color',[0.40 0.40 0.40], ...
        'Units','normalized');
    end

    xlim(ax,[0 1]); ylim(ax,[0 1]); axis(ax,'equal','off');
    hold(ax,'off');
end