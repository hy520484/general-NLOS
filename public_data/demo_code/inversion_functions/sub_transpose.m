function Y = sub_transpose(Bu, num_nei)
% A block-by-block transpose of the block dataset
[pxy,Nnn2D] = size(Bu);                                             
N = Nnn2D / num_nei;                                                   
Y = NaN * zeros(num_nei,pxy * N);                                             
for j = 1:N                                                   
    temp_block = Bu(:, num_nei*(j-1) + 1 : num_nei*(j-1) + num_nei);               
    Y(:,pxy*(j-1) + 1 : pxy*(j-1) + pxy) = temp_block';                    
end
end