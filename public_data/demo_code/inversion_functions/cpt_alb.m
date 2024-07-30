function alb = cpt_alb(u)
% This function computes the albedo of the directional albedo.
alb = sqrt(sum(u.^2,4));
end

