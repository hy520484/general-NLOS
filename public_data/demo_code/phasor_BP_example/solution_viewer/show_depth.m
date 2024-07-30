function show_depth(depth,save_name)
    min_depth = min(depth(:));
    max_depth = max(depth(:));
    figure('Position',[50,50,320,400]);
    axes('Position',[0.02 0.212 0.96 0.768]);
    tmp_depth_img = imagesc(depth);
    colormap parula
    caxis([min_depth max_depth])
    set(tmp_depth_img,'alphadata',~isnan(depth));
    current_figure = gca;
    set(current_figure,'color',[0 0 0]);
    current_figure.TickLength = [0 0];
    axis equal
    current_figure.XTick = [];
    current_figure.YTick = [];
    
    % add a colorbar
    c = colorbar;
    c.Location = 'south';
    c.Position = [0.02 0.11 0.96 0.08];
    c.Limits = [min_depth,max_depth];
    c.TickLength = 0;
    c.Ticks = [min_depth + 0.09 * (max_depth - min_depth ) ,max_depth -  0.09 * (max_depth - min_depth )];
    c.TickLabels = {num2str(min_depth,'%1.2f'),num2str(max_depth,'%1.2f')};
    c.FontName = 'Times New Roman';
    c.FontSize = 23;
    c.AxisLocation = 'out';
    c.Box = 'on';
    saveas(gcf,[save_name,'.svg'])
    pause
end