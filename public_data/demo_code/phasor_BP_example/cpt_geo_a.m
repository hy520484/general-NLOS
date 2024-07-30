function [Arr_bin_i,Arr_bin_d, Des_value_i, Des_value_d, min_bin, max_bin,...
          xmin_i, ymin_i, zmin_i, xmin_d, ymin_d, zmin_d] =...
          cpt_geo_a(c,delta_t, xbc,ybc,zbc, bx,by,bz, xec,yec,zec, C_i, C_d)

ixbc = min(C_i(:,1)); ixec = max(C_i(:,1));
iybc = min(C_i(:,2)); iyec = max(C_i(:,2));
izbc = min(C_i(:,3)); izec = max(C_i(:,3));

xmin_i = xbc - ixec; xmax_i = xec - ixbc;
ymin_i = ybc - iyec; ymax_i = yec - iybc;
zmin_i = zbc - izec; zmax_i = zec - izbc;

dxbc = min(C_d(:,1)); dxec = max(C_d(:,1));
dybc = min(C_d(:,2)); dyec = max(C_d(:,2));
dzbc = min(C_d(:,3)); dzec = max(C_d(:,3));

xmin_d = xbc - dxec; xmax_d = xec - dxbc;
ymin_d = ybc - dyec; ymax_d = yec - dybc;
zmin_d = zbc - dzec; zmax_d = zec - dzbc;

[JJ_i,KK_i,LL_i] = meshgrid(ymin_i:ymax_i,xmin_i:xmax_i,zmin_i:zmax_i);
JJ_i = by * JJ_i; KK_i = bx * KK_i; LL_i = bz * LL_i;

[JJ_d,KK_d,LL_d] = meshgrid(ymin_d:ymax_d,xmin_d:xmax_d,zmin_d:zmax_d);
JJ_d = by * JJ_d; KK_d = bx * KK_d; LL_d = bz * LL_d;

n_vec_i = sqrt(JJ_i.^2 + KK_i.^2 + LL_i.^2);
n_vec_d = sqrt(JJ_d.^2 + KK_d.^2 + LL_d.^2);

Arr_bin_i = n_vec_i / c / delta_t;
Arr_bin_d = n_vec_d / c / delta_t;

Des_value_i = 1 ./ n_vec_i .^2;
Des_value_d = 1 ./ n_vec_d .^2;

min_bin = ceil(min(Arr_bin_i(:)) + min(Arr_bin_d(:)));
max_bin = ceil(max(Arr_bin_i(:)) + max(Arr_bin_d(:)));
end