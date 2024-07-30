%% Part I: speed of light & temporal resolution
% 1.1 speed of light
c = 3E8;
% 1.2 temporal resolution
delta_t = 16E-12;

%% Part II: size and location of the reconstruction domain
% 2.1 size of a basic voxel
bx = c * delta_t; by = 0.01; bz = 0.01;
% 2.2 location of the reconstruction domain
xb = 0.85; yb = 0.00; zb = 0.00;
xe = 1.00; ye = 1.27; ze = 1.27;
% 2.3 number of voxels extended for reconstruction
ext_x = 0; ext_y = 0; ext_z = 0;
% 2.4 number of basic voxels in each inversion voxel
nbx = 1; nby = 2; nbz = 2;

% compute the size of the extended reconstruction domain 
xb = xb - nbx * bx * ext_x; yb = yb - nby * by * ext_y; zb = zb - nbz * bz * ext_z;
xe = xe + nbx * bx * ext_x; ye = ye + nby * by * ext_y; ze = ze + nbz * bz * ext_z;

% compute the index of voxels in three dimensions
xbc = floor(xb / bx); ybc = floor(yb / by); zbc = floor(zb / bz);
xec = ceil (xe / bx); yec = ceil (ye / by); zec = ceil (ze / bz);

% compute the number of voxels in each dimension 
xg = xbc:nbx:xec; lxs = length(xg);
yg = ybc:nby:yec; lys = length(yg);
zg = zbc:nbz:zec; lzs = length(zg);

% number of inversion voxels in total
num_voxel = lxs * lys * lzs;

%% Part III: Location of full rectangular grids on the relay surface
% full imaginary illuminated points of rectangular grids
% endding with *_f

% 3.1 grid indices in horizonal and vertical directions
grid_hor_f = 126:-2:0;
grid_ver_f = 126:-2:0;

% number of points in each direction
np_hor_f = length(grid_hor_f);
np_ver_f = length(grid_ver_f);

% generating rectangular focal points (C_i: Coordinates of illuminated points)
num_meas_f = np_hor_f * np_ver_f;
C_i_f = zeros(num_meas_f,3);
count_meas_f = 0;
for j = 1:np_ver_f
    for k = 1:np_hor_f
        count_meas_f = count_meas_f + 1;
        C_i_f(count_meas_f,2) = grid_hor_f(k);
        C_i_f(count_meas_f,3) = grid_ver_f(j);
    end
end
clear count_meas_f
C_d_f = C_i_f; % Confocal case

%% Part IV: Location of grid points where illumination occurs
% endding without *_f
Sig_flag = sum(sub_Sig,2) > 0;
sub_Sig(Sig_flag == 0,:) = [];
C_i(Sig_flag == 0,:) = [];
C_d(Sig_flag == 0,:) = [];
num_meas = size(C_i,1);
% recording real measurements
[real_meas, real_order] = ismember([C_i_f,C_d_f], [C_i,C_d], 'rows'); 
real_order(real_order == 0) = []; [~,real_order] = sort(real_order,'ascend');
[inter_meas,inter_order] = ismember([C_i,C_d], [C_i_f,C_d_f], 'rows'); 
inter_order(inter_order == 0) = []; [~,inter_order] = sort(inter_order,'ascend');