function [depth_gt, alb_gt, normal_gt] = get_oracle(u_gt, x_cor)
%% maximum albedo value projection
alb = cpt_alb(u_gt);
max_alb = max(alb(:));
u_gt = u_gt / max_alb;
alb = alb / max_alb;
max_alb_ind = get_max_alb_ind(alb);
u_proj = rot90(get_proj_u(u_gt,max_alb_ind));
alb_proj = rot90(get_proj_alb(alb,max_alb_ind));
%% converting view
u_proj = u_proj(:,end:-1:1,:);
alb_proj = alb_proj(:,end:-1:1);
u_proj(:,:,1) = - u_proj(:,:,1);
u_proj(:,:,2) =   u_proj(:,:,2);
u_proj(:,:,3) = - u_proj(:,:,3);
% get normal
normal_gt = u_proj ./ repmat(alb_proj,1,1,3);
% get albedo
alb_gt = alb_proj;
% get depth 
max_alb_ind = rot90(max_alb_ind);
max_alb_ind = max_alb_ind(:,end:-1:1);
depth_gt = get_depth(max_alb_ind, x_cor);
end