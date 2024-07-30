% This file implements the proposed "confocal-complemented signal-object
% collaborative regularization" (CC-SOCR) algorithm described in the paper
% "Non-line-of-sight imaging with arbitrary illumination and detection pattern".

%--------------------------------------------------------------------------
% Loading data
id = 21;
switch id
    case {1}
        Instance_name = 'bunny';
        Relay_setting = 'full';
        threshold_rate = 0.07;
    case {2}
        Instance_name = 'bunny';
        Relay_setting = 'horizontal';
        threshold_rate = 0.07;
    case {3}
        Instance_name = 'bunny';
        Relay_setting = 'sticks';
        threshold_rate = 0.07;
    case {4}
        Instance_name = 'bunny';
        Relay_setting = 'vertical';
        threshold_rate = 0.07;
    case {5}
        Instance_name = 'bunny';
        Relay_setting = 'window';
        threshold_rate = 0.07;
    case {6}
        Instance_name = 'figure_4';
        Relay_setting = 'full';
        threshold_rate = 0.15;
    case {7}
        Instance_name = 'figure_4';
        Relay_setting = 'horizontal';
        threshold_rate = 0.15;
    case {8}
        Instance_name = 'figure_4';
        Relay_setting = 'random';
        threshold_rate = 0.15;
    case {9}
        Instance_name = 'figure_4';
        Relay_setting = 'sparse';
        threshold_rate = 0.15;
    case {10}
        Instance_name = 'figure_4';
        Relay_setting = 'vertical';
        threshold_rate = 0.15;
    case {11}
        Instance_name = 'non_planar_letters_NT';
        Relay_setting = 'full';
        threshold_rate = 0.20;
    case {12}
        Instance_name = 'non_planar_letters_NT';
        Relay_setting = 'horizontal';
        threshold_rate = 0.20;
    case {13}
        Instance_name = 'non_planar_letters_NT';
        Relay_setting = 'oval';
        threshold_rate = 0.20;        
    case {14}
        Instance_name = 'non_planar_letters_NT';
        Relay_setting = 'random';
        threshold_rate = 0.20;
    case {15}
        Instance_name = 'non_planar_letters_NT';
        Relay_setting = 'sparse';
        threshold_rate = 0.20;
    case {16}
        Instance_name = 'non_planar_letters_NT';
        Relay_setting = 'vertical';
        threshold_rate = 0.20;
    case {17}
        Instance_name = 'pyramid';
        Relay_setting = 'edge';
        threshold_rate = 0.25;
    case {18}
        Instance_name = 'statue';
        Relay_setting = 'full';
        threshold_rate = 0.15;
    case {19}
        Instance_name = 'statue';
        Relay_setting = 'heart';
        threshold_rate = 0.15;
    case {20}
        Instance_name = 'statue';
        Relay_setting = 'horizontal';
        threshold_rate = 0.15;
    case {21}
        Instance_name = 'statue';
        Relay_setting = 'NLOS';
        threshold_rate = 0.15;
    case {22}
        Instance_name = 'statue';
        Relay_setting = 'random';
        threshold_rate = 0.15;
    case {23}
        Instance_name = 'statue';
        Relay_setting = 'sparse';
        threshold_rate = 0.15;
    case {24}
        Instance_name = 'statue';
        Relay_setting = 'sticks';
        threshold_rate = 0.15;
    case {25}
        Instance_name = 'statue';
        Relay_setting = 'vertical';
        threshold_rate = 0.15;
end
cd(['./data/',Instance_name]);
% loading data
load(['./',Relay_setting,'/C_i.mat'])
load(['./',Relay_setting,'/C_d.mat'])
load(['./',Relay_setting,'/sub_Sig.mat'])
% experimental setup
experimental_setup
% parameters
parameters
cd ..
cd ..
save_results_dir = ['results_',Instance_name,'_',Relay_setting];
mkdir(save_results_dir);

%--------------------------------------------------------------------------
%% Initializing parallel pool
parpool(8, 'IdleTimeout', inf);

%--------------------------------------------------------------------------
% The CC-SOCR algorithm starts from here
tic;
addpath('./inversion_functions')
%% Computing the forward operator
[Arr_bin_i,Arr_bin_d, Des_value_i, Des_value_d, min_bin, max_bin, xmin_i, ymin_i, zmin_i, xmin_d, ymin_d, zmin_d] = cpt_geo(c,delta_t, xbc,ybc,zbc, bx,by,bz, xec,yec,zec, C_i, C_d);
[Arr_bin_i_f,Arr_bin_d_f, Des_value_i_f, Des_value_d_f, min_bin_f, max_bin_f, xmin_i_f, ymin_i_f, zmin_i_f, xmin_d_f, ymin_d_f, zmin_d_f] = cpt_geo(c,delta_t, xbc,ybc,zbc, bx,by,bz, xec,yec,zec, C_i_f, C_d_f);
max_bin_f = max(max_bin, max_bin_f);
%% 0. Pre-processing raw signal
num_bin = size(sub_Sig, 2);
if num_bin < max_bin_f
    ext_bin = max_bin_f - num_bin;
    sub_Sig = [sub_Sig,zeros(num_meas, ext_bin)];
else
    sub_Sig = sub_Sig(:,1:max_bin_f);
end
sub_Sig(sub_Sig < 0) = 0; sub_Sig = sub_Sig / max(sub_Sig(:)) * 255; sub_Sig = double(sub_Sig);
%% Stage I: The initialization stage
%% 1.1 Initializing the estimated signal b
b = sub_Sig; b(abs(b) < s_b_imp) = 0;
%% 1.2 Initializing the reconstruction u
BP_b = NLOS_adjoint(ones(1,num_meas), xmin_i,ymin_i,zmin_i,xmin_d,ymin_d,zmin_d, xg,yg,zg, C_i, C_d, Arr_bin_i,Arr_bin_d, Des_value_i, Des_value_d, b);
u = NLOS_reconstruction_I(ones(1,num_meas), ones(1,num_meas), xmin_i,ymin_i,zmin_i,xmin_d,ymin_d,zmin_d, xg,yg,zg, max_bin, C_i,C_d, Arr_bin_i,Arr_bin_d, Des_value_i,Des_value_d, zeros(lxs,lys,lzs,3), BP_b, num_iter_CG_small, r_tol_CG_coarse, 0);
b_u = forward(C_i,C_d,u,max_bin_f,xmin_i,ymin_i,zmin_i,xmin_d,ymin_d,zmin_d,xg,yg,zg,Arr_bin_i,Arr_bin_d,Des_value_i,Des_value_d);
% choosing parameters
alb = cpt_alb(u); alb = alb(:);
s_u = s_u_imp * sum((b(:) - b_u(:)).^2) / sum(alb);
alb(alb == 0) = [];
mu = mu_imp * s_u * mean(1 ./ alb);
% split Bregman iteration
q = zeros(lxs,lys,lzs,3);
R_uv_init = NaN * zeros(1,Init_Bregman_iter);
for init_Breg_loop = 1:Init_Bregman_iter
% update v
temp_var = u - q;
alb = cpt_alb(temp_var);
v = repmat(max(0,1 - s_u / mu / 2 ./alb),1,1,1,3) .* temp_var;
% update u
temp_var = q + v;
u = NLOS_reconstruction_I(ones(1,num_meas), ones(1,num_meas), xmin_i,ymin_i,zmin_i,xmin_d,ymin_d,zmin_d, xg,yg,zg, max_bin, C_i,C_d, Arr_bin_i,Arr_bin_d, Des_value_i,Des_value_d, v, BP_b + mu * temp_var, num_iter_CG_small, r_tol_CG_coarse, mu);
R_uv_init(init_Breg_loop) = norm(u(:) - v(:)) / norm(u(:));
if R_uv_init(init_Breg_loop) < R_uv_init_tol
fprintf('At iteration %d, R(u,v) = %1.6f\n', init_Breg_loop, R_uv_init(init_Breg_loop));
    break
end
% update q
q = q + v - u;
end
u = v;
b_u = forward(C_i,C_d,u,max_bin_f,xmin_i,ymin_i,zmin_i,xmin_d,ymin_d,zmin_d,xg,yg,zg,Arr_bin_i,Arr_bin_d,Des_value_i,Des_value_d);
d_u = forward(C_i_f,C_d_f,u,max_bin_f,xmin_i_f,ymin_i_f,zmin_i_f,xmin_d_f,ymin_d_f,zmin_d_f,xg,yg,zg,Arr_bin_i_f,Arr_bin_d_f,Des_value_i_f,Des_value_d_f);
%% 1.3 Initializing the dictionaries D_n, D_s, and the coefficients C
amp_alb = 255 / max(alb(:));
[uDC,Ds,Dn] = target_dict( u, amp_alb, lambda_pu_imp, lambda_RMSE, pxo,pyo,pzo, sxo,syo,szo, wso, nno, sparse_value_o, sparse_threshold_o);
%% 1.4 Initializing the Wiener coefficients S
Sig_S = reshape(Wiener(reshape(b_u,num_meas,1,[]), reshape( (sub_Sig + lambda_sb * b)/(1+lambda_sb),num_meas,1,[]), sigma_b / sqrt(1 + lambda_sb),pxb,pyb,pzb, sxd,syd,szd), num_meas, []);
%% 1.5 Initializing the virtual confocal signal d
d = d_u;
d(real_meas == 0,:) = hard_threshold(d(real_meas == 0,:) , s_d_imp);
temp_b = b(inter_meas == 1,:); temp_b = temp_b(inter_order,:);
d(real_meas == 1,:) = hard_threshold((lambda_d * d(real_meas == 1,:) + lambda_bd * temp_b) / (lambda_d + lambda_bd) , s_d_imp);
%% 1.6 Dictionary learning of the simulated signal and the virtual confocal signal $\Psi$, Q
[PsiQ, Dxyz] = DDTF( reshape((d + lambda_sd * d_u) / (1 + lambda_sd), np_hor_f, np_ver_f, max_bin_f ), lambda_fd_imp, pxq, pyq, pzq, sxq, syq, szq);
PsiQ = reshape( PsiQ, num_meas_f,[]);
%% Stage II: The iteration stage
for num_loop = 1:num_main_loop
    %% 2.1 Updating the estimated signal b
    b(inter_meas == 0,:) = hard_threshold( (1 * b_u(inter_meas == 0,:) + lambda_b * sub_Sig(inter_meas == 0,:) + lambda_b * lambda_pb * lambda_sb * Sig_S(inter_meas == 0,:)                                  )/(1 + lambda_b + lambda_b * lambda_pb * lambda_sb),...
                           sqrt(s_b / (1 + lambda_b + lambda_b * lambda_pb * lambda_sb)));
    temp_d = d(real_meas == 1,:); temp_d = temp_d(real_order,:);
    b(inter_meas == 1,:) = hard_threshold( (1 * b_u(inter_meas == 1,:) + lambda_b * sub_Sig(inter_meas == 1,:) + lambda_b * lambda_pb * lambda_sb * Sig_S(inter_meas == 1,:) + lambda_bd * temp_d)/(1 + lambda_b + lambda_b * lambda_pb * lambda_sb + lambda_bd),...
                           sqrt(s_b / (1 + lambda_b + lambda_b * lambda_pb * lambda_sb + lambda_bd)));
    %% 2.2 Updating the Wiener coefficients S
    Sig_S = reshape(Wiener(reshape(b_u,num_meas,1,[]),  reshape( (sub_Sig + lambda_sb * b)/(1+lambda_sb),num_meas,1,[]),  sigma_b / sqrt(1 + lambda_sb),  pxb,pyb,pzb, sxd,syd,szd),   num_meas, []);
    %% 2.3 Updating the reconstructed target u
    if num_loop == 1
        err_b = sum((b_u(:) - b(:)).^2);               
        err_uDC = sum((u(:) - uDC(:)).^2);
        err_d = sum((d_u(:) - d(:)).^2);
        err_PsiQ = sum((d_u(:) - PsiQ(:)).^2);
        alb = cpt_alb(u); L1_norm_alb = sum(alb(:));    
        lambda_u = lambda_u_imp * err_b / err_uDC;
        lambda_d = lambda_d_imp * err_b / err_d;
        s_u = s_u * (1 + lambda_u_imp + lambda_d_imp * (1 + lambda_pd * lambda_sd));
        mu = mu * (1 + lambda_u_imp + lambda_d_imp * (1 + lambda_pd * lambda_sd));
    end
    q = zeros(lxs,lys,lzs,3);
    R_uv_main = NaN * zeros(1,Main_Bregman_iter);
    for main_Breg_loop = 1:Main_Bregman_iter
        % update v
        temp_var = u - q;
        alb = cpt_alb(temp_var);
        v = repmat(max(0,1 - s_u / mu / 2 ./alb),1,1,1,3) .* temp_var;  
        % update u
        alpha = lambda_d + lambda_d * lambda_pd * lambda_sd;
        y_full = (lambda_d * d + lambda_d * lambda_pd * lambda_sd * PsiQ) / (lambda_d + lambda_d * lambda_pd * lambda_sd);
        beta = lambda_u + mu;
        y_error = (lambda_u * uDC + mu * (v + q)) / (lambda_u + mu);
        % BP_b
        BP_b =  NLOS_adjoint(ones(1,num_meas), xmin_i,ymin_i,zmin_i,xmin_d,ymin_d,zmin_d, xg,yg,zg, C_i, C_d, Arr_bin_i,Arr_bin_d, Des_value_i, Des_value_d, b);        
        BP_b = BP_b + alpha * NLOS_adjoint(ones(1,num_meas_f), xmin_i_f,ymin_i_f,zmin_i_f,xmin_d_f,ymin_d_f,zmin_d_f, xg,yg,zg, C_i_f, C_d_f, Arr_bin_i_f,Arr_bin_d_f, Des_value_i_f, Des_value_d_f, y_full);
        BP_b = BP_b + beta * y_error;
        % reconstruct the target
        u = NLOS_reconstruction_II(ones(1,num_meas),alpha * ones(1,num_meas_f), ones(1,num_meas),ones(1,num_meas_f), xmin_i,ymin_i,zmin_i,xmin_d,ymin_d,zmin_d, xmin_i_f,ymin_i_f,zmin_i_f,xmin_d_f,ymin_d_f,zmin_d_f, xg,yg,zg, max_bin,max_bin_f, C_i,C_d, C_i_f,C_d_f, Arr_bin_i,Arr_bin_d,Arr_bin_i_f,Arr_bin_d_f, Des_value_i,Des_value_d, Des_value_i_f,Des_value_d_f, v, BP_b, num_iter_CG_small, r_tol_CG_coarse, beta);                           
        R_uv_main(main_Breg_loop) = norm(u(:) - v(:)) / norm(u(:));
        if R_uv_main(main_Breg_loop) < R_uv_main_tol
           fprintf('At iteration %d, R(u,v) = %1.6f\n', main_Breg_loop, R_uv_main(main_Breg_loop));
           break
        end
        % update q
        q = q + v - u;
    end
    u = v;
    cd(save_results_dir), save(['u',num2str(num_loop),'.mat'],'u'), cd ..
    
    if num_loop == num_main_loop
        break
    end
    
    b_u = forward(C_i,C_d,u,max_bin_f,xmin_i,ymin_i,zmin_i,xmin_d,ymin_d,zmin_d,xg,yg,zg,Arr_bin_i,Arr_bin_d,Des_value_i,Des_value_d);
    d_u = forward(C_i_f,C_d_f,u,max_bin_f,xmin_i_f,ymin_i_f,zmin_i_f,xmin_d_f,ymin_d_f,zmin_d_f,xg,yg,zg,Arr_bin_i_f,Arr_bin_d_f,Des_value_i_f,Des_value_d_f); 
    %% 2.4 Updating the dictionaries Dn, Ds and the coefficients C
    [uDC,Ds,Dn] = target_dict( u, amp_alb, lambda_pu_imp, lambda_RMSE, pxo,pyo,pzo, sxo,syo,szo, wso, nno, sparse_value_o, sparse_threshold_o);
    %% 2.5 Updating the virtual confocal signal d
    s_d = s_d_imp * lambda_d * (1 + lambda_pd);
    d(real_meas == 0,:) = hard_threshold((d_u(real_meas == 0,:) + lambda_pd * PsiQ(real_meas == 0,:)) / (1 + lambda_pd), sqrt(s_d / (lambda_d + lambda_d * lambda_pd)));
    temp_b = b(inter_meas == 1,:); temp_b = temp_b(inter_order,:);
    d(real_meas == 1,:) = hard_threshold((lambda_d * d_u(real_meas == 1,:) + lambda_d * lambda_pd * PsiQ(real_meas == 1,:) + lambda_bd * temp_b) / (lambda_d + lambda_d * lambda_pd + lambda_bd), sqrt(s_d / (lambda_d + lambda_d * lambda_pd + lambda_bd)));
    cd(save_results_dir), save(['d',num2str(num_loop),'.mat'],'d'), cd ..
    %% 2.6 Updating the dictionary of the simulated signal and the virtual confocal signal
    [PsiQ, Dxyz] = DDTF( reshape((d + lambda_sd * d_u) / (1 + lambda_sd), np_hor_f, np_ver_f, max_bin_f ), lambda_fd_imp, pxq, pyq, pzq, sxq, syq, szq);
    PsiQ = reshape(PsiQ, num_meas_f,[]);
end
%% Stage III: Show the reconstructed albedo and the surface normal
addpath('solution_viewer')
u(repmat(u(:,:,:,1),1,1,1,3) > 0) = 0;
alb = cpt_alb(u);
max_alb = max(alb(:));
u = u / max_alb;
alb = alb / max_alb;
u(repmat(alb,1,1,1,3) < threshold_rate) = 0;
alb(alb < threshold_rate) = 0;
max_alb_ind = get_max_alb_ind(alb);
u_proj = rot90(get_proj_u(u,max_alb_ind));
alb_proj = rot90(get_proj_alb(alb,max_alb_ind));
% viewing from the relay surface
u_proj = u_proj(:,end:-1:1,:);
alb_proj = alb_proj(:,end:-1:1);
u_proj(:,:,1) = - u_proj(:,:,1);
u_proj(:,:,2) =   u_proj(:,:,2);
u_proj(:,:,3) = - u_proj(:,:,3);
% showing projections
cd(save_results_dir)
show_2D(u_proj(:,:,1),'hot','x-component');
show_2D(u_proj(:,:,2),'hot','y-component');
show_2D(u_proj(:,:,3),'hot','z-component');
show_2D(alb_proj, 'hot','albedo');
% three views
three(alb, 5);
saveas(gcf,'three view.png')
%% Show illumination points
figure
plot3(bx * C_i(:,1), by * C_i(:,2), bz * C_i(:,3),'k.')
view(33,29)
grid on
xlabel('depth (m)')
ylabel('horizontal (m)')
zlabel('vertical (m)')
title('Coordinates of the illumination points')
saveas(gcf,'illumination_pattern.png')
cd ..
toc;