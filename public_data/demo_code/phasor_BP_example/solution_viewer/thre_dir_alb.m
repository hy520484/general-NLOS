function u = thre_dir_alb(u, threshold_rate)
if length(size(u)) == 4 % directional albedo
    alb = cpt_alb(u);
    max_alb = max(alb(:));
    u = u / max_alb;
    alb = alb / max_alb;
    u(repmat(alb,1,1,1,3) < threshold_rate) = 0;
else % albedo
    max_alb = max(u(:));
    u = u / max_alb;
    u(u < threshold_rate) = 0;
end
end

