function view_directional_albedo(u, method, threshold_rate)
%% maximum albedo value projection
alb = cpt_alb(u);
max_alb = max(alb(:));
u = u / max_alb;
alb = alb / max_alb;
max_alb_ind = get_max_alb_ind(alb);
u_proj = rot90(get_proj_u(u,max_alb_ind));
alb_proj = rot90(get_proj_alb(alb,max_alb_ind));
%% converting view
u_proj = u_proj(:,end:-1:1,:);
alb_proj = alb_proj(:,end:-1:1);
u_proj(:,:,1) = - u_proj(:,:,1);
u_proj(:,:,2) =   u_proj(:,:,2);
u_proj(:,:,3) = - u_proj(:,:,3);
%% showing projections
show_2D(u_proj(:,:,1),'hot',[method,'-x']);
show_2D(u_proj(:,:,2),'hot',[method,'-y']);
show_2D(u_proj(:,:,3),'hot',[method,'-z']);
show_2D(alb_proj, 'hot',[method,'-albedo']);
%% Three views
three(alb, 100);
saveas(gcf,[method,'-three.svg'])
alb(alb < threshold_rate) = 0;
three(alb,200);
saveas(gcf,[method,'-threshold.svg'])
end