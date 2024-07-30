function [denoised_alb, Dxy, Dz] = dict_learn_alb( alb, amp_alb, lambda_pu,...
                                                   lambda_RMSE, pxo, pyo, pzo, sxo,syo, szo,...
                                                   wso, nno, sparse_value_o, sparse_threshold_o)
% This function carries out the dictionary learning for the reconstructed
% target.

% basic parameters
[m,n,o] = size(alb); wlo = 2 * wso + 1; 
N_iter = 10; N_xy_iter = 5; N_z_iter = 2; 
pxyz = pxo * pyo * pzo;
mpx = m - pxo + 1; npy = n - pyo + 1; opz = o - pzo + 1;
Px = 1:sxo:mpx; if Px(length(Px))~=mpx, Px = [Px,mpx]; end, lx = length(Px);
Py = 1:syo:npy; if Py(length(Py))~=npy, Py = [Py,npy]; end, ly = length(Py);
Pz = 1:szo:opz; if Pz(length(Pz))~=opz, Pz = [Pz,opz]; end, lz = length(Pz);

% amplification
alb = alb * amp_alb;

% constructing block                                                    
Bu = NaN * zeros(pxyz,lx*ly*lz*nno);
part_block = 0;
for j = 1:lx
    for k = 1:ly
        for l = 1:lz
            jj = Px(j); kk = Py(k); ll = Pz(l);
            ref_patch = alb(jj:jj+pxo-1,kk:kk+pyo-1,ll:ll+pzo-1);
            ref_patch = ref_patch(:);
            if sum(abs(ref_patch) > sparse_value_o) < sparse_threshold_o
                continue
            end
            count_num = 0;
            temp_block = NaN * zeros(pxyz,wlo^3);
            MSE = NaN * zeros(1,wlo^3);
            for x = -wso:wso
                for y = -wso:wso
                    for z = -wso:wso
                        tj = jj + x; tk = kk + y; tl = ll + z;
                        if tj >= 1 && tj <= mpx &&...        
                           tk >= 1 && tk <= npy &&...
                           tl >= 1 && tl <= opz
                            count_num = count_num + 1;
                            nei_patch = alb(tj:tj+pxo-1,tk:tk+pyo-1,tl:tl+pzo-1);       
                            nei_patch = nei_patch(:);
                            temp_block(:,count_num) = nei_patch;
                            MSE(count_num) = sum((ref_patch - nei_patch).^2)/pxyz;
                        end
                    end
                end
            end

            % check current block
            [MSE_ascend_value,MSE_ascend_order] = sort(MSE,'ascend');
            f = find(MSE_ascend_value < lambda_RMSE.^2);

            if length(f) < nno 
                continue
            else
               temp_block = temp_block(:,MSE_ascend_order(1:nno));
               Bu(:,part_block + 1 : part_block + nno) = temp_block;
               part_block = part_block + nno;
            end

        end
    end
end
Bu(:,isnan(Bu(1,:))) = [];

% initializing dictionaries
Dxy = kron(dctmtx(pyo)',dctmtx(pxo)'); Dxy = kron(dctmtx(pzo)',Dxy);
Dz = dctmtx(nno)';

% updating dictionaries
for j = 1:N_iter
    Coef_xy = sub_transpose(Dxy' * Bu,nno);
    for j_sub = 1:N_z_iter
        C = hard_threshold(Dz' * Coef_xy, lambda_pu);
        [U,~,V] = svd(Coef_xy * C'); Dz = U * V';
    end
    Coef_z = sub_transpose(Dz' * sub_transpose(Bu,nno),pxyz);
    for j_sub = 1:N_xy_iter
        C = hard_threshold(Dxy' * Coef_z, lambda_pu);
        [U,~,V] = svd(Coef_z * C'); Dxy = U * V';
    end      
end

% aggregation
ebuff = zeros(m,n,o);
wbuff = zeros(m,n,o);
for j = 1:lx
    parfor k = 1:ly
        for l = 1:lz
            jj = Px(j); kk = Py(k); ll = Pz(l);
            ref_patch = alb(jj:jj+pxo-1,kk:kk+pyo-1,ll:ll+pzo-1);
            ref_patch = ref_patch(:);     
            if sum(abs(ref_patch) > sparse_value_o) < sparse_threshold_o
                continue
            end
            temp_est_mtx = zeros(m,n,o);
            temp_w_mtx = zeros(m,n,o);
            count_num = 0;
            temp_block = NaN * zeros(pxyz,wlo^3);
            MSE = NaN * zeros(1,wlo^3);
            nei_x = NaN * zeros(1,wlo^3);
            nei_y = NaN * zeros(1,wlo^3);
            nei_z = NaN * zeros(1,wlo^3);
            for x = -wso:wso
                for y = -wso:wso
                    for z = -wso:wso
                        tj = jj + x; tk = kk + y; tl = ll + z;
                        if tj >= 1 && tj <= mpx &&...        
                           tk >= 1 && tk <= npy &&...
                           tl >= 1 && tl <= opz
                            count_num = count_num + 1;
                            nei_patch = alb(tj:tj+pxo-1,tk:tk+pyo-1,tl:tl+pzo-1);
                            nei_patch = nei_patch(:);
                            temp_block(:,count_num) = nei_patch;
                            MSE(count_num) = sum((ref_patch - nei_patch).^2)/pxyz;
                            nei_x(count_num) = x; nei_y(count_num) = y; nei_z(count_num) = z;
                        end
                    end
                end
            end
            
            [MSE_ascend_value,MSE_ascend_order] = sort(MSE,'ascend');
            f = find(MSE_ascend_value < lambda_RMSE.^2);
            
            if length(f) < nno
               temp_Dz = dctmtx(length(f))';
            else
               temp_Dz = Dz;
               f = f(1:nno);
            end
            
            count_num = length(f);
            noisy_block = temp_block(:,MSE_ascend_order(f));
            nei_x = nei_x(MSE_ascend_order(f));
            nei_y = nei_y(MSE_ascend_order(f));
            nei_z = nei_z(MSE_ascend_order(f)); 
            temp_coef = Dxy' * noisy_block;
            temp_coef = temp_Dz' * temp_coef';
            temp_coef(abs(temp_coef) < lambda_pu) = 0;
            temp_weight = 1;
            temp_coef = temp_Dz * temp_coef;
            denoised_block = Dxy * temp_coef';
            
            for nei = 1:count_num
                        temp_patch = reshape(denoised_block(:,nei),pxo,pyo,pzo);
                        tnx = jj + nei_x(nei); tny = kk + nei_y(nei); tnz = ll + nei_z(nei);
                        temp_w_mtx(tnx:tnx+pxo-1,tny:tny+pyo-1,tnz:tnz+pzo-1) =...
                            temp_w_mtx(tnx:tnx+pxo-1,tny:tny+pyo-1,tnz:tnz+pzo-1) + temp_weight;
                        temp_est_mtx(tnx:tnx+pxo-1,tny:tny+pyo-1,tnz:tnz+pzo-1) =...
                            temp_est_mtx(tnx:tnx+pxo-1,tny:tny+pyo-1,tnz:tnz+pzo-1) +  temp_weight * temp_patch;
            end
            
            ebuff = ebuff + temp_est_mtx;
            wbuff = wbuff + temp_w_mtx; 

        end
    end
end

denoised_alb = ebuff./wbuff;

% post processing
denoised_alb(isnan(denoised_alb)) = 0;
denoised_alb(denoised_alb < 0) = 0;

% amplification
denoised_alb = denoised_alb / amp_alb;
end