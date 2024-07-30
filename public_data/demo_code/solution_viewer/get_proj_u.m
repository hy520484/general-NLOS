function proj = get_proj_u(u,ind)
% This function generates the front view of the directional albedo with the
% provided indices.
[~,n,o,~] = size(u);
proj = zeros(n,o,3);
for j = 1:n
    for k = 1:o
        if isnan(ind(j,k)), continue, end
        proj(j,k,:) = u(ind(j,k),j,k,:);
    end
end
end
