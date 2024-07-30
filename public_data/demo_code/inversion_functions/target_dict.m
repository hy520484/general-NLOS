function [uDC,Dxy,Dz] = target_dict( u, amp_alb, lambda_pu, lambda_RMSE, pxo,pyo,pzo, sxo,syo,szo,wso, nno, sparse_value_o, sparse_threshold_o)
% This function computes the dictionaries and the representation of the
% directional albedo.
alb = sqrt(sum(u.^2,4));
[denoised_alb, Dxy, Dz] = dict_learn_alb( alb, amp_alb, lambda_pu,lambda_RMSE,...
                                          pxo, pyo, pzo, sxo,syo, szo,wso, nno, sparse_value_o, sparse_threshold_o);
uDC = repmat(denoised_alb,1,1,1,3) .* u ./ repmat(alb,1,1,1,3);
uDC(isnan(uDC)) = 0;
end