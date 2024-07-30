%% PF-BP
lambda_times = 2;
cycle_times  = 4;
addpath(genpath(pwd))
load('sub_Sig.mat')
experimental_setup
[num_meas, num_bin] = size(sub_Sig);
[Arr_bin_i,Arr_bin_d, Des_value_i, Des_value_d, min_bin, max_bin,...
xmin_i, ymin_i, zmin_i, xmin_d, ymin_d, zmin_d] =...
cpt_geo_a(c,delta_t, xbc,ybc,zbc, bx,by,bz, xec,yec,zec, C_i, C_d);
sub_Sig = Sig_alignment(sub_Sig,max_bin);
wavelength = lambda_times * 2 * by;
gridsize = wavelength*1/4;
P = phasor_pulse(wavelength, cycle_times, c * delta_t, sub_Sig);

R_P = NLOS_adjoint_a(ones(num_meas,1), xmin_i,ymin_i,zmin_i, xmin_d,ymin_d,zmin_d, xg,yg,zg,...
                      C_i, C_d, Arr_bin_i,Arr_bin_d, Des_value_i, Des_value_d, real(P));
I_P = NLOS_adjoint_a(ones(num_meas,1), xmin_i,ymin_i,zmin_i, xmin_d,ymin_d,zmin_d, xg,yg,zg,...
                      C_i, C_d, Arr_bin_i,Arr_bin_d, Des_value_i, Des_value_d, imag(P));
u = abs(complex(R_P, I_P));
clear R_P I_P
thresh_value = 0.15;
threshold_rate = 0.15;
u = thre_dir_alb(u, threshold_rate);
view_albedo(u, 'phasor-back-projection', thresh_value)