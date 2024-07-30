function proj = get_proj_u(u,ind)
[~,n,o,~] = size(u);
proj = zeros(n,o,3);
for j = 1:n
    for k = 1:o
        if isnan(ind(j,k)), continue, end
        proj(j,k,:) = u(ind(j,k),j,k,:);
    end
end
end
