function proj = get_proj_alb(alb,ind)
% This function generates the front view of the albedo with the provided
% indices.
[~,n,o] = size(alb);
proj = zeros(n,o);
for j = 1:n
    for k = 1:o
        if isnan(ind(j,k)), continue, end
        proj(j,k) = alb(ind(j,k),j,k);
    end
end
end
