%% ODT Reconstruction
clear
useGPU(1);
load('simulation_shep_model_CG_n0_1.33_dn_0.10_Nviews_41_Nx_256_Ny_256_Nz_72_dx_0.10_dy_0.10_dz_0.10_GaussianBeam_0_photon_1_illum_circle_-45_45.mat');
%% Load simulated data

main_param_LSm;% Set the parameters
load('input_Ryt.mat');% Load initial guess of Rytov reconstruction
LSm;