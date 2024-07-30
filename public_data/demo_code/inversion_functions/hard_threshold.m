function x = hard_threshold( x,lambda )
% Hard-thresholding operator
x(abs(x)<lambda) = 0;
end

