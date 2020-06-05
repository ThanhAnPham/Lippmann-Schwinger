%% Assessment methods from [1]
%[1] Ayoub, A. B., Pham, T. A., Lim, J., Unser, M., & Psaltis, D. (2020).
%    A method for assessing the fidelity of optical diffraction tomography 
%    reconstruction methods using structured illumination.
%    Optics Communications, 454, 124486.
clear;
close all;
useGPU(0);
set(0,'DefaultTextInterpreter','latex');
%% Load param
load('simulation_shep_model_CG_n0_1.33_dn_0.10_Nviews_41_Nx_256_Ny_256_Nz_72_dx_0.10_dy_0.10_dz_0.10_GaussianBeam_0_photon_1_illum_circle_-45_45_light.mat');

f = gpuCpuConverter((simP.k*simP.dx)^2*((n/simP.n0).^2 - 1));% Scattering Potential
sz = size(f);

%% Set the sensors position (outside the volume)
SP.Nmeas = [150,150];%Number of sensors
%XY positions: origin at the center of the sample
[SP.x,SP.y] = meshgrid((simP.dx:simP.dx:simP.dx*SP.Nmeas(1)) - simP.dx*SP.Nmeas(1)/2,...
    (simP.dy:simP.dy:simP.dy*SP.Nmeas(2)) - simP.dy*SP.Nmeas(2)/2);

SP.y = SP.y(:); SP.x = SP.x(:);

%one pixel outside of the sample.
SP.z = simP.dz*(sz(3)/2+1)*ones(length(SP.y),1);

%% Build the incident field (EPFL logo)

epfl = 2.5*flip(double(imread('epfl.png') > 0),1);
epfl = imresize(epfl(1 + end/2 -72:end/2 + 72,1 + end/2 -72:end/2 + 72),sz(1:2),'nearest');
uin0 = exp(1i*epfl);
uin0 = padarray(uin0,(SP.Nmeas - size(uin0))/2,1,'both');

par.distUin = simP.Lz_vec - simP.Lz/2;par.Nz = length(par.distUin);

uin = gpuCpuConverter(setFields3D(uin0,[],par,simP,'lipp','uin'));

S_sp2vol = LinOpSelectorPatch(size(uin),1 + size(uin)/2 - sz/2, size(uin)/2 + sz/2);
uin = S_sp2vol*uin;

%% Build the forward model based on Lippmann-Schwinger

Niter = 200;% # iterations to compute the forward model
xtol = 1e-9;% tolerance for convergence criteria
par.distuM = -SP.z(end);% Focal plane at the center of the volume

% Green function in the volume
G = FreeSpaceKernel3D(simP.k*simP.dx,sz,'Vainikko');
% Green function to get the scattered field at the sensor positions
Gtil = Vainikkotilde(simP.k,sz,SP,[simP.dx,simP.dy,simP.dz],ones(1,3));

H = OpLipp3D(G,Gtil,LinOpIdentity(sz),uin,Niter,xtol,[],false);

P = LinOpTiltedPropagator(simP.lambda0,simP.n0,par.distuM,simP.dx,SP.Nmeas,[0,0],'AS');
%% Compute the scattered field on the sensors positions
tic
ysc = P*H*f;%scattered field
y = uin0 + ysc;%total field
toc
%% Ayoub method for assessment of the quality of reconstruction

%reverse the fields (conjugate)
yr = conj(y);
uinr = S_sp2vol*gpuCpuConverter(setFields3D(yr,[],par,simP,'lipp','uin'));

%% Load reconstructions from different methods and set the reverse forward model
load('recon_born.mat');
load('recon_ryt.mat');
load('recon_lipp.mat');

H.MGtild = Gtil;
H.uin = uinr;
%% Estimate the uin0 from them
tic
uin0_born = conj(yr + P*H*f_born);
uin0_ryt = conj(yr + P*H*f_ryt);
uin0_lipp = conj(yr + P*H*f_lipp);
toc

%% Plots and Metrics
figure(1);set(gcf,'Units','normalized','Position',[0.5,0,0.5,1]);
subplot(421);
imagesc(abs(uin0));axis image;colorbar;title('Ground truth (ampl)');caxis([0,2]);
subplot(422);
imagesc(angle(uin0));axis image;colorbar;title('Ground truth (phase)');caxis([-pi,pi]);
subplot(423);
imagesc(abs(uin0_born));axis image;colorbar;title('Born (ampl)');caxis([0,2]);
subplot(424);
imagesc(angle(uin0_born));axis image;colorbar;title('Born (phase)');caxis([-pi,pi]);
subplot(425);
imagesc(abs(uin0_ryt));axis image;colorbar;title('Ryt (ampl)');caxis([0,2]);
subplot(426);
imagesc(angle(uin0_ryt));axis image;colorbar;title('Ryt (phase)');caxis([-pi,pi]);
subplot(427);
imagesc(abs(uin0_lipp));axis image;colorbar;title('LSm (ampl)');caxis([0,2]);
subplot(428);
imagesc(angle(uin0_lipp));axis image;colorbar;title('LSm (phase)');caxis([-pi,pi]);
set(findall(gcf,'Type','Axes'),'YDir','normal');

fprintf('-------------------------Table of comparison---------------------------------\n');
fprintf('RMSE Born : %1.2e\n',norm(angle(uin0_born(:)) - angle(uin0(:)))/sqrt(nnz(uin0)));
fprintf('RMSE Rytov : %1.2e\n',norm(angle(uin0_ryt(:)) - angle(uin0(:)))/sqrt(nnz(uin0)));
fprintf('RMSE LSm : %1.2e\n',norm(angle(uin0_lipp(:)) - angle(uin0(:)))/sqrt(nnz(uin0)));
fprintf('-----------------------------------------------------------------------------\n');