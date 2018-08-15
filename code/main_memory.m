addpath(genpath('tomography/myFun'));
addpath(genpath('InvPbLib'));

disp('Press enter');
pause;

load('/Users/anpham/Google Drive/ODT/code/mat files/fastLipp/fine_grid_shep_noisy_n0_1.33_dn_0.13_Nviews_31_Nx_1024_GaussianBeam_0_noise_splitGauss_1_1_dB.mat',...
    'M','SP','SP_discrete','f','n','simP','uM_true');


%main_seagle;
main_lippman;