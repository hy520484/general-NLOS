function proj = get_proj_alb(alb,ind)
[~,n,o] = size(alb);
proj = zeros(n,o);
for j = 1:n
    for k = 1:o
        if isnan(ind(j,k)), continue, end
        proj(j,k) = alb(ind(j,k),j,k);
    end
end
end
