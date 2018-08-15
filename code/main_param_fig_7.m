%% Script for parameters initialisation
%% Load a data mat file
%   - simP : param structure with fields
%            -ObjectType : Name of the samples (string for filename saving)
%            -lambda0 : wavelength
%            -k : wavenumber 2 pi/simP.lambda0
%            -n0 : background refractive index (RI)
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
%load('fine_grid_shep_noisy_n0_1.33_dn_0.13_Nviews_31_Nx_1024_GaussianBeam_0_noise_splitGauss_0_0_dB.mat',...
%    'M', 'ObjType', 'SP', 'SP_discrete', 'f', 'n','simP', 'uM', 'uM_true', 'uem')
clear par;
%Dataset provided by
%Jean-Michel Geffrin, Pierre Sabouroux and Christelle Eyraud, “Free space experimental scattering database
%continuation: experimental set-up and measurement precision,” Inverse Probl. 21, (2005).
load('FoamDielExtTM_3_GHz_conj_1.mat');

par.param_fname = strcat(mfilename(),'.m');
%% Parameters

%Preprocessing of data
simP.Lxroi = simP.Lx;%Experimental data has artifact at the border, it determines the ROI before any propagation, etc.

%BPM(i) : area where field is propagated (avoids periodic bound' artifacts)
%Otherwise : area including the sensors positions
par.Lxext = 1*simP.Lx;
par.Lzext = 1*simP.Lz;% for simulated data

%Effective area on which algo optimizes
par.Ntheta = simP.Ntheta;
par.Nx = 256;
par.Nz = 256;

%Discretisation on area delimited by par.Lxext x par.Lzext
par.Nxext = par.Nx;
par.Nzext = par.Nz;

par.dx = par.Lxext/par.Nxext;%par.Lx/par.Nx;
par.dz = par.Lzext/par.Nzext;%par.Lz/par.Nz;

par.Lx = par.Nx*par.dx;%par.Lzext/1.5;%
par.Lz = par.Nz*par.dz;%par.Lxext/1.5;%

%maximal angle scanning (if scalar, assumes symmetric : 180 -> [-90,90], can be an interval [-30,45])
par.max_scan_theta = deg2rad(360);
%Exclude any angle (if too noisy, etc), set [] if none;
%index corresponds to the full dataset
par.excluded_angles = [];%[80:83];%[80,82,85:87,101:103];

par.noreflection = false;

par.SP = SP;
%% REAL DATA
%Where the data measurements are placed (at the end of the physical window
%in the optical axis by default
%par.SP.z = par.Lzext*ones(length(SP.z),1);
%% Algo parameters (FISTA or ADMM for complex and intensity rsp.)
par.Lsub = 8; %Number of angles considered for gradient evaluation (stochastic or not)
par.stochGrad = 0;%Stochastic Gradient descent or switch of one with equally angle spaced Lsub angles
par.gam = 5e-3;%Step size for FBS (FISTA)
par.redStep = 0;%reduced step at each iteration of FISTA
par.Niter = 2e2;%Number of iterations
par.xtol = 1e-8;%Tolerance for stopping criterion
par.ItUpOut = 25;%Save every ItUpOut iterations
par.doPlot = 1;%plots within algorithms
par.iterVerb = 1;%Iterations display
par.folder = pwd;%Folder where data are saved (in subdirectory with algorithm name)
par.discreteGtild = false;%true for FFT based Gtild, false for matrix based Gtild
par.storeGtild = true;%Define if Gtild matrix is stored for matrix based Gtild
par.yscat = true;

%Inner parameters for iterative forward (lippmann) OR Intensity only for
%the inner FISTA
par.Niterin = 200;%number of iterations for iterative forward/inner FISTA
par.Nupdatein = 0;%Nupdate > 0 activates displays every 50 iterations
par.xtolin = 1e-10;%Tolerance for stopping criterion
%and iterates saved in F.OutOpIn every Nupdate for SEAGLE

%TV ADMM parameters
par.penalty = 'l21';%'l20','l21'
par.regul = 1.8e-2;%Regularization parameter
par.bounds = [0,inf];%[0 inf];%set constraint (positivity [0 inf], or e.g. [0,0.1])
par.prox_iter = 7.5e2;%number of iterations for ADMM
par.rhoTV = 4e1*[0.5,1];%Penalty parameters (rho)

par.field_scale = 1;%Scale the fields (equivalent to higher regul)
%Define how Uin is built
%   -'sim' : analytical expression for 
%   -'real' : propagated from measured uin
par.uin = 'cylindrical';%'sim';%'real'

%O0NLY used when par.uin == 'real'

par.shift = 0*par.dz*20;

par.distUin = (0:par.Nzext-1)*par.dz - par.Lzext/2 + par.shift;%others
%Data measurements
uM.curr_y = uM.GT;%or any suitable variable
par.distuM = 0;%propagation distance,
%Special requirements (cut, etc. (for real data mainly))
%e.g. uM.curr_y = uM.curr_y(max(1 + end/2 - 708,1):min(end/2+708,end),:);

par.Mroi = 1;%par.Lxext/simP.Lx;
[par, M_curr,SP] = setParamODT(par,simP,SP);
%% Special param set, loop over it
% might require few adaptations in the algorithm script (e.g. main_lippmann)

%regularization params
regul_set = 1.6e-2;%[1e-2:1e-3:1.5e-2];%[1.6e-2:1e-3:2e-2];%[1.8e-2,4e-3,1e-2];

%Step size
gam_set = [5e-3];
%% Initial guess
%main_backproj;
%% Create a missing uin (if Intensity only, no uin available in reality)
if ~exist('uin','var')
    uin = [];
end

fprintf('dx : %1.2e, dz : %1.2e, ROI Lx x Lz : %1.2e x %1.2e, Lx x Lz : %1.2e x %1.2e\n',...
    par.dx,par.dz,par.Lx,par.Lz,par.Lxext,par.Lzext);