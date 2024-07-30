%% Part I
% CG solver
num_iter_CG_small = 20;
r_tol_CG_coarse = 5E-3;
% number of iterations
num_main_loop = 1;
Init_Bregman_iter = 20; R_uv_init_tol = 5E-3;
Main_Bregman_iter = 3; R_uv_main_tol = 5E-3;
%% Part II
% noise levels
sigma_b = 40; lambda_fd_imp = 40;
% target dictionary learning parameters
lambda_pu_imp = 2.6 * 60;
lambda_RMSE = 15; pxo = 3; pyo = 3; pzo = 3; sxo = 1; syo = 1; szo = 1;
wso = 3; nno = 3; sparse_value_o = 70; sparse_threshold_o = 6;
% patch size of complete signal
pxb = 1; pyb = 1; pzb = 3; sxd = 1; syd = 1; szd = 1;
pxq = 3; pyq = 3; pzq = 3; sxq = 1; syq = 1; szq = 1;
%% Part III
s_u_imp = 100; mu_imp = 0.05;
lambda_u_imp = 10;
lambda_b = 1; lambda_sb = 0.25; lambda_pb = 16;
lambda_d_imp = 2; lambda_sd = 1; lambda_pd = 4;
s_b_imp = 255 * 0.01; s_b = lambda_b * s_b_imp^2;
s_d_imp = 255 * 0.01;
lambda_bd = 4; lambda_d = 2;