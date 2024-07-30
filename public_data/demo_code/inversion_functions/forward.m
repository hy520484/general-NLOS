function Sig = forward(C_i,C_d,u,max_bin,xmin_i,ymin_i,zmin_i,xmin_d,ymin_d,zmin_d,xg,yg,zg,Arr_bin_i,Arr_bin_d,Des_value_i,Des_value_d)
% This function generates the simulated signal with the provided directional albedo.
num_measure = size(C_i,1);
Sig = NaN * zeros(num_measure,max_bin);

parfor s = 1:num_measure
       id_ind_i = C_i(s,:); id_ind_d = C_d(s,:);
       temp_Arr_i = Arr_bin_i(xg - id_ind_i(1) - xmin_i + 1, yg - id_ind_i(2) - ymin_i + 1, zg - id_ind_i(3) - zmin_i + 1);
       temp_Arr_d = Arr_bin_d(xg - id_ind_d(1) - xmin_d + 1, yg - id_ind_d(2) - ymin_d + 1, zg - id_ind_d(3) - zmin_d + 1);
       temp_Arr = ceil(temp_Arr_i + temp_Arr_d);

       temp_Des_i = Des_value_i(xg - id_ind_i(1) - xmin_i + 1, yg - id_ind_i(2) - ymin_i + 1, zg - id_ind_i(3) - zmin_i + 1);
       temp_Des_d = Des_value_d(xg - id_ind_d(1) - xmin_d + 1, yg - id_ind_d(2) - ymin_d + 1, zg - id_ind_d(3) - zmin_d + 1,:);
       
       temp = sum(temp_Des_d .* u, 4) .* temp_Des_i;
       
       Sig(s,:) = full(sparse(1,temp_Arr(:),temp(:),1,max_bin));
end

end