function y = three(albedo, fig_num)
% This function shows the three views of the reconstructed albedo.

min_albedo = 0; 
max_albedo = max(albedo(:));

figure(fig_num)
% front view
    subplot(1,3,1)
    temp_img = rot90(squeeze(max(albedo,[],1)));
    imagesc(temp_img(:,end:-1:1));
    xlabel('y');
    ylabel('z');
    colormap('gray');
    caxis([min_albedo, max_albedo])
    axis equal;
    title('front view')
    axis off

% top view  
    subplot(1,3,2)
    temp_img = squeeze(max(albedo,[],3));
    imagesc(temp_img(end:-1:1,end:-1:1));
    xlabel('y');
    ylabel('x');
    colormap('gray');
    caxis([min_albedo, max_albedo])
    axis equal;
    title('top view')
    axis off
    
% side view
    subplot(1,3,3)
    temp_img = rot90(squeeze(max(albedo,[],2)));
    imagesc(temp_img);
    xlabel('x');
    ylabel('z');
    colormap('gray');
    caxis([min_albedo, max_albedo])
    axis equal;
    title('side view')
    axis off

    y = [];
end