function x = NLOS_reconstruction_II(Sig_weight,Sig_weight_f, Sig_flag,Sig_flag_f, xmin_i,ymin_i,zmin_i,xmin_d,ymin_d,zmin_d,xmin_i_f,ymin_i_f,zmin_i_f,xmin_d_f,ymin_d_f,zmin_d_f, xg,yg,zg, max_bin,max_bin_f, C_i,C_d,C_i_f,C_d_f, Arr_bin_i,Arr_bin_d, Arr_bin_i_f,Arr_bin_d_f, Des_value_i,Des_value_d, Des_value_i_f,Des_value_d_f, init_guess, RHS, max_iter, max_r_tol, mu)
% Least-squares reconstruction with the virtual confocal signal

% basic parameters
lxs = length(xg); lys = length(yg); lzs = length(zg); 
num_measure = size(C_i,1); num_measure_f = size(C_i_f,1); 

% CG
init_guess = init_guess(:); RHS = RHS(:);
nb = norm(RHS); r_err = NaN * zeros(1,max_iter);
x = init_guess;
r = RHS - A(x);

r_err(1) = norm(r)/nb;
rr = r' * r; p = r;

if r_err(1) >= max_r_tol
    for k = 2:max_iter
        Ap = A(p);
        pAp = p' * Ap;
        alpha = rr / pAp;
        x_new = x + alpha * p;
        r_new = r - alpha * Ap;
        r_err(k) =  norm(r_new)/nb;
        if r_err(k) < max_r_tol
            x = x_new;
            break
        end
        rr_new = r_new' * r_new;
        beta = rr_new / rr;
        p = r_new + beta * p;
        x = x_new;
        r = r_new;
        rr = rr_new;
    end
end

x = reshape(x,lxs,lys,lzs,3);

function y = A(x)
y = zeros(lxs,lys,lzs,3);
x = reshape(x,lxs,lys,lzs,3);

parfor s = 1:num_measure
    if Sig_flag(s) == 0, continue, end
    id_ind_i = C_i(s,:); id_ind_d = C_d(s,:);

    temp_Arr_i = Arr_bin_i(xg - id_ind_i(1) - xmin_i + 1, yg - id_ind_i(2) - ymin_i + 1, zg - id_ind_i(3) - zmin_i + 1);
    temp_Arr_d = Arr_bin_d(xg - id_ind_d(1) - xmin_d + 1, yg - id_ind_d(2) - ymin_d + 1, zg - id_ind_d(3) - zmin_d + 1);
    temp_Arr = ceil(temp_Arr_i + temp_Arr_d);

    temp_Des_i = Des_value_i(xg - id_ind_i(1) - xmin_i + 1, yg - id_ind_i(2) - ymin_i + 1, zg - id_ind_i(3) - zmin_i + 1);
    temp_Des_d = Des_value_d(xg - id_ind_d(1) - xmin_d + 1, yg - id_ind_d(2) - ymin_d + 1, zg - id_ind_d(3) - zmin_d + 1,:);
    
    temp = sum(temp_Des_d .* x, 4) .* temp_Des_i;
    
    mid_sig = full(sparse( 1, temp_Arr(:), temp(:), 1, max_bin ));
    
    temp_y = repmat(mid_sig(temp_Arr),1,1,1,3) .* temp_Des_d .* repmat(temp_Des_i,1,1,1,3);
    y = y + Sig_weight(s) * temp_y;
end
parfor s = 1:num_measure_f
    if Sig_flag_f(s) == 0, continue, end
    id_ind_i = C_i_f(s,:); id_ind_d = C_d_f(s,:);

    temp_Arr_i = Arr_bin_i_f(xg - id_ind_i(1) - xmin_i_f + 1, yg - id_ind_i(2) - ymin_i_f + 1, zg - id_ind_i(3) - zmin_i_f + 1);
    temp_Arr_d = Arr_bin_d_f(xg - id_ind_d(1) - xmin_d_f + 1, yg - id_ind_d(2) - ymin_d_f + 1, zg - id_ind_d(3) - zmin_d_f + 1);
    temp_Arr = ceil(temp_Arr_i + temp_Arr_d);

    temp_Des_i = Des_value_i_f(xg - id_ind_i(1) - xmin_i_f + 1, yg - id_ind_i(2) - ymin_i_f + 1, zg - id_ind_i(3) - zmin_i_f + 1);
    temp_Des_d = Des_value_d_f(xg - id_ind_d(1) - xmin_d_f + 1, yg - id_ind_d(2) - ymin_d_f + 1, zg - id_ind_d(3) - zmin_d_f + 1,:);
    
    temp = sum(temp_Des_d .* x, 4) .* temp_Des_i;
    
    mid_sig = full(sparse( 1, temp_Arr(:), temp(:), 1, max_bin_f ));
    
    temp_y = repmat(mid_sig(temp_Arr),1,1,1,3) .* temp_Des_d .* repmat(temp_Des_i,1,1,1,3);
    y = y + Sig_weight_f(s) * temp_y;
end

y = y + mu * x;
y = y(:);
end

end