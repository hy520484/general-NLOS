function [noisy_img, Dxyz] = DDTF(noisy_img,noise_level,pxd,pyd,pzd,sxd,syd,szd)
% This function carries out the dictionary learning process for blocks in a 3D
% tensor.

% basic parameters
[m,n,o] = size(noisy_img);
mpx = m - pxd + 1; npy = n - pyd + 1; opz = o - pzd + 1;
Px = 1:sxd:mpx; if Px(length(Px))~=mpx, Px = [Px,mpx]; end, lx = length(Px);
Py = 1:syd:npy; if Py(length(Py))~=npy, Py = [Py,npy]; end, ly = length(Py);
Pz = 1:szd:opz; if Pz(length(Pz))~=opz, Pz = [Pz,opz]; end, lz = length(Pz);

% DCT dictionary
Dxyz = kron(dctmtx(pzd) , kron(dctmtx(pyd) , dctmtx(pxd)))';

% constructing the patch dataset
Pu = NaN * zeros(pxd * pyd * pzd, lx * ly * lz);
row = 0;
count_num = zeros(m,n,o);
for j = 0:pxd-1
    for k = 0:pyd-1
        for l = 0:pzd-1
            row = row + 1;
            temp_noisy = noisy_img(Px+j,Py+k,Pz+l);
            count_num(Px+j,Py+k,Pz+l) = count_num(Px+j,Py+k,Pz+l) + 1;
            Pu(row,:) = temp_noisy(:)';
        end
    end
end

for num_iter = 1:30
    % updating the coefficients
    noisy_coef = Dxyz' * Pu;
    noisy_coef = hard_threshold(noisy_coef, 2.6 * noise_level);
    % updating the dictionary
    [U,~,V] = svd(Pu * noisy_coef'); Dxyz = U * V';
end

% frequency hard-thresholding
noisy_coef = Dxyz' * Pu;
noisy_coef = hard_threshold(noisy_coef, 2.6 * noise_level);
Pu = Dxyz * noisy_coef;

% aggregation
noisy_img = zeros(m,n,o);
row = 0;
for j = 0:pxd-1
    for k = 0:pyd-1
        for l = 0:pzd-1
            row = row + 1;
            noisy_img(Px+j,Py+k,Pz+l) = noisy_img(Px+j,Py+k,Pz+l) + reshape(Pu(row,:),lx,ly,lz); 
        end
    end
end

% post-processing
noisy_img = noisy_img ./ count_num;
noisy_img(noisy_img < 0) = 0;
end