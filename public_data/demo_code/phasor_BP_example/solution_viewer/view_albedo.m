function view_albedo(u, method, threshold_rate)
%% maximum albedo value projection
u(u<0) = 0;
u = u / max(u(:));
max_u_ind = get_max_alb_ind(u);
alb_proj = rot90(get_proj_alb(u,max_u_ind));
%% converting view
alb_proj = alb_proj(:,end:-1:1);
%% showing projections
show_2D(alb_proj, 'hot',[method,'-albedo']);
%% Three views
three(u, 100);
saveas(gcf,[method,'-three.svg'])
u(u < threshold_rate) = 0;
three(u, 200);
saveas(gcf,[method,'-threshold.svg'])
end