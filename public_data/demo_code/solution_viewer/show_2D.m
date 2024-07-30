function show_2D(img,colormap_name,save_name)
% This function saves an image with a colorbar.
min_img = min(img(:));
max_img = max(img(:));
% show the figure
figure('Position',[50,50,320,400]);
axes('Position',[0.02 0.212 0.96 0.768]);
imagesc(img), colormap(colormap_name), axis equal, axis off
% add a colorbar
c = colorbar;
c.Location = 'south';
c.Position = [0.02 0.11 0.96 0.08];
c.Limits = [min_img,max_img];
c.TickLength = 0;
c.Ticks = [min_img + 0.09 * (max_img - min_img ) ,max_img -  0.09 * (max_img - min_img )];
c.TickLabels = {num2str(min_img,'%.2f'),num2str(max_img,'%.2f')};
c.FontName = 'Times New Roman';
c.FontSize = 23;
c.AxisLocation = 'out';
c.Box = 'on';
if nargin == 3
    saveas(gcf,[save_name,'.png'])
end
end