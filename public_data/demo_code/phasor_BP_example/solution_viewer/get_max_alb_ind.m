function max_alb_ind = get_max_alb_ind(alb)
[~,n,o] = size(alb);
max_alb_ind = NaN * zeros(n,o);
for j = 1:n
    for k = 1:o
        temp_alb = alb(:,j,k);
        max_temp_alb = max(temp_alb);
        if max_temp_alb == 0, continue, end
        temp_alb_ind = find(temp_alb == max_temp_alb);
        temp_alb_ind = temp_alb_ind(1);
        max_alb_ind(j,k) = temp_alb_ind;
    end
end
end