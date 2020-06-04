%% Script for parameters initialisation
%% Load a data mat file
%   - simP : param structure with fields
%            -ObjectType : Name of the samples (string for filename saving)
%            -lambda0 : wavelength
%            -n0 : background refractive index (RI)
%            -k : wavenumber 2 pi/simP.lambda0*simP.n0
%            -thetas : angles vector (radian)
%            -Ntheta : Numer of angles (length of simP.thetas)
%            -Ltheta : angle scanning (max(simP.thetas) - min(simP.thetas))
%            -Nx, Nz : Number of pixels (experimental : Nx determined by camera, Nz arbitrary)
%            -Lx, Lz : Physical size (same as Nx,Nz)
%            -dx, dz : Physical step size (dx = Lx/Nx, etc.)
%            Optional
%            -GaussianBeam, sigmaBeam : boolean for Gaussian beam approx.
%                                       (for simulation mainly) and its associated std

%   - uM.GT : data measurement in matrix form M x Ntheta
%   - uin : incident field either of same size than uM.GT or full volume
%           (measured or simulated rsp.).
%   - SP.x / .z : physical positions of each camera pixel wrt the physical
%                 length (origin at the left bottom)
%   Optional
%   - log_imag_ur2D : for classical Rytov and cie, unwrapped phase diff.
%   - f : Scattering potential ground truth
%   - n : Refractive index ground truth

%Example
%load('/Users/anpham/Google Drive/ODT/code/mat files/Intensity exp 2/to do/fine_grid_bead_noisy_n0_1.33_dn_0.06_Nviews_61_Nx_1024_GaussianBeam_0_noise_splitGauss_30_30_dB_radius_1.683e-06.mat',...
%    'M', 'ObjType', 'SP', 'SP_discrete', 'f', 'iuin', 'n', 'simP', 'uM', 'uM_true', 'uem')
%load('2017-06-10GlassFiberLarger_mean_1_doOpp_0_conj_0.mat');
%load('2017-06-10GlassFiber_mean_1_doOpp_0_invTheta_0_extractPos_470.mat');
%load('2017-06-10GlassFiberLarger_mean_1_doOpp_0_conj_0.mat');
%load('/Volumes/GoogleDrive/My Drive/ODT/code/mat files/3D mat/simulation_bead_noisy_n0_1.33_dn_0.06_Nviews_41_Nx_256_Ny_256_Nz_256_GaussianBeam_1_photon_1_radius_3.500e+00.mat',...
%    'M', 'SP', 'f', 'uin', 'n', 'simP', 'uM', 'uem','log_imag_ur');
%load('/Volumes/2ToMac/real_data/2017-06-10GlassFiber/2017-06-10GlassFiber_mean_1_doOpp_0_invTheta_0_extractPos_512.mat')
clear par;

%% Parameters
par.param_fname = strcat(mfilename(),'.m');

%Preprocessing of data
simP.Lxroi = simP.dx*simP.Nx;%Experimental data has artifact at the border, it determines the ROI before any propagation, etc.
simP.Lyroi = simP.dy*simP.Ny;

%BPM(i) : area where field is propagated (avoids periodic bound' artifacts)
%Otherwise : area including the sensors positions

N = [simP.Nx,simP.Ny,simP.Nz];

par.Lxext = simP.dx*N(1);
par.Lyext = simP.dy*N(2);
par.Lzext = simP.dz*N(3);

%Effective area on which algo optimizes
par.Ntheta = simP.Ntheta;
par.Nx = simP.Nroi(1);
par.Ny = simP.Nroi(2);
par.Nz = simP.Nroi(3);

%Discretisation on area delimited by par.Lxext x par.Lzext
%Note for BPM : set par.Nxext = par.Nx in general.
par.Nxext = simP.Nx;%par.Nx;%par.Lxext/par.dx;
par.Nyext = simP.Ny;%par.Ny;%par.Lxext/par.dx;
par.Nzext = simP.Nz;%par.Nz;%par.Lzext/par.dz;

par.dx = par.Lxext/par.Nxext;
par.dy = par.Lyext/par.Nyext;
par.dz = par.Lzext/par.Nzext;

par.Lx = par.Nx*par.dx;
par.Ly = par.Ny*par.dy;
par.Lz = par.Nz*par.dz;

%maximal angle scanning (if scalar, assumes symmetric : 180 -> [-90,90], can be an interval [-30,45])
par.max_scan_theta = {deg2rad(360),deg2rad(360)};
%Exclude any angle (if too noisy, etc), set [] if none;
%index corresponds to the full dataset
par.excluded_angles = [];
par.Ntheta = simP.Ntheta;
par.noreflection = 0;

par.SP = SP;
%% Algo parameters (FISTA or ADMM for complex and intensity rsp.)
par.Lsub = 8; %Number of angles considered for gradient evaluation (stochastic or not)
par.stochGrad = 1;%Stochastic Gradient descent or switch of one with equally angle spaced Lsub angles
par.gam = 2e3;%Step size for FBS (FISTA)
par.redStep = false;%reduced step at each iteration of FISTA
par.accel = 'fista';
par.momRestart = false;
par.alpha = 1;
par.Niter = 200;%Number of iterations
par.ItUpOut = 10;%Save every ItUpOut iterations
par.xtol = 1e-10;%Tolerance for stopping criterion
par.doPlot = 1;%plots within algorithms
par.iterVerb = 1;%Iterations display
par.folder = './';%'/Volumes/2ToMac/3D';%Folder where data are saved (in subdirectory with algorithm name)
par.saveX = 25;

%LIPPMANN/Rytj only (for instance)
par.discreteGtild = 0;%true for FFT based Gtild, false for matrix based Gtild (more flexible)
%par.storeGtild = false;%Define if Gtild matrix is stored for matrix based Gtild
par.yscat = 1;
par.memoizeApply = 1;

%Intensity only algo
par.rho = 1;%ADMM

%Inner parameters for iterative forward (lippmann)
par.Niterin = 500;%number of iterations for iterative forward

%Intensity only INNER FISTA
par.NiterInnerFISTA = 50;
par.doSpecialConv = false;
%Nupdate > 0 activates displays every 50 iterations and iterates saved in
%F.OutOpIn every Nupdate for Lipp
par.Nupdatein = 10;
par.xtolin = eps;%Tolerance for stopping criterion
par.reset = 1;%reset gam (step)

%TV parameters
par.RegOpti = 'FGP';%ADMM'
par.RegType = 'TV';%'Hess','TV'
par.penalty = 'l21';%penalty type 'l20','l10'
par.regul = 5e-2;%Regularization parameter
par.bounds = [0,inf];%[0 inf];%set constraint (positivity [0 inf], or e.g. [0,0.1])
par.prox_iter = 10;%number of iterations for ADMM
par.LTV = 12; %Parameter for FGP
par.rhoTV = 4e2*[0.5,1];%4e1*[0.5,1];%Penalty parameters (rho) for ADMM

par.field_scale = 1;%Scale the fields (equivalent to higher regul)
%Define how Uin is built
%   -'sim' : analytical expression for 
%   -'real' : propagated from measured uin
par.uin = 'sim';%'sim';%'real';%'cylindrical';%'planewave';

%ONLY used when par.uin == 'real'
%par.distUin = -par.Lx/2 + par.dz*20;%BPM
par.shift = 0*simP.dz*2;
%par.distUin = -par.Lz/2;% - par.dx;%BPM
par.distUin = max(unique(SP.z)) - simP.Lz/2 + (0:par.Nz-1)*par.dz + par.shift - par.Lz/2;%others
%Data measurements

par.sigma = 0;
uM.curr_y = uM.GTc;%or any suitable variable
curr_uem = uem;

uM.curr_y = addNoise(uM.curr_y,par.sigma,'Gauss');
curr_uem = addNoise(curr_uem,par.sigma,'Gauss');

par.distuM = 0*((par.Lzext/2) + par.shift);%propagation distance,-par.dz
%Special requirements (cut, etc. (for real data mainly))
%e.g. uM.curr_y = uM.curr_y(max(1 + end/2 - 708,1):min(end/2+708,end),:);

par.Zsen = [];%par.Lz + par.dz;%allows us to place the sensor anywhere we want. Should be outside the reconstruction domain

par.Mroi = [1,1];%par.Lxext/simP.Lx;%Not working well
[par, M_curr,curr_SP] = setParamODT3D(par,simP,SP);

par.snr = snr(uM.GTc - uem,uM.GTc - uM.curr_y);
fprintf('Measurement SNR: %1.2f dB\n',par.snr);
%% Special param set, loop over it
% might require few adaptations in the algorithm script (e.g. main_lippmann)

%regularization params
regul_set = [1e-2];

regul_set = regul_set/prod(par.siz);

%% Create a missing uin (if Intensity exp only, no uin available)
clear new_uwr curr_par

if ~exist('uin','var')
    uin = [];
end
if ~exist('uem','var')
    uem = uin;
end

fprintf('dx : %1.2e, dy : %1.2e, dz : %1.2e, ROI Lx x Ly x Lz : %1.2e x %1.2e x %1.2e, Lx x Ly x Lz (ext): %1.2e x %1.2e x %1.2e\n',...
    par.dx,par.dy,par.dz,par.Lx,par.Ly,par.Lz,par.Lxext,par.Lyext,par.Lzext);

fprintf('Number of angles : %i\n',par.Ntheta);