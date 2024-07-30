function target = NLOS_adjoint_a(Sig_flag, xmin_i,ymin_i,zmin_i,xmin_d,...
                  ymin_d,zmin_d, xg,yg,zg, C_i, C_d, Arr_bin_i,Arr_bin_d, Des_value_i, Des_value_d, Sig)

lxs = length(xg); lys = length(yg); lzs = length(zg);
num_measure = size(Sig,1);

target = zeros(lxs,lys,lzs);

for s = 1:num_measure
       if Sig_flag(s) == 0, continue, end
       temp_sig = Sig(s,:); 
       id_ind_i = C_i(s,:); id_ind_d = C_d(s,:);
       
       temp_Arr_i = Arr_bin_i(xg - id_ind_i(1) - xmin_i + 1, yg - id_ind_i(2) - ymin_i + 1, zg - id_ind_i(3) - zmin_i + 1);
       temp_Arr_d = Arr_bin_d(xg - id_ind_d(1) - xmin_d + 1, yg - id_ind_d(2) - ymin_d + 1, zg - id_ind_d(3) - zmin_d + 1);
       temp_Arr = ceil(temp_Arr_i + temp_Arr_d);

       temp_Des_i = Des_value_i(xg - id_ind_i(1) - xmin_i + 1, yg - id_ind_i(2) - ymin_i + 1, zg - id_ind_i(3) - zmin_i + 1);
       temp_Des_d = Des_value_d(xg - id_ind_d(1) - xmin_d + 1, yg - id_ind_d(2) - ymin_d + 1, zg - id_ind_d(3) - zmin_d + 1);
       
       target = target + temp_sig(temp_Arr) .* temp_Des_i .* temp_Des_d;

end

end