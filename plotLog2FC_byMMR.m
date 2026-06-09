function plotLog2FC_byMMR(fc, yVar, titleStr)

    mmr = categorical(fc.MMR_status);
    mmrCats = categories(mmr);

    for i = 1:numel(mmrCats)
        m = mmrCats{i};
        T = fc(mmr == m, :);
        if height(T) == 0, continue; end

        figure;
        boxchart(T.(yVar));
        %myboxplot2d(T.(yVar));
        yline(0,'k--');
        ylabel('log_2 fold change');
        title(sprintf('%s — %s', titleStr, m), 'Interpreter','none');
        set(gca,'Box','off','XTick',[]);
    end
end
