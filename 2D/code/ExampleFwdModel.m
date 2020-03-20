%-------------------------------------------------------------------------
% This script contains an example of use of the iterative Lippman-Schwinger
% forward model. Given the incident field and the scattering potentiel, it
% computes the total field.
%-------------------------------------------------------------------------
clear; close;

%% Load inputs
addpath(genpath('../'));
load('uin');            % Load incident field
load('f');              % Load scattering potential f
lamb=406-9;             % wavelength
nb=1.33;                % refractive index background
siz=size(f);            % size of the region of interest (containing the support of f)
dz=16*lamb/siz(1);      % axial discretization step (o have a ROI of 16*lamb
kdz=2*pi/lamb*nb*dz;    % wavenumber

%% Instanciation and application of the forward model
% - Instanciate operators
G = LinOpFreeSpaceKernel(kdz,siz,2*siz);     % Convolution operator with the Green's function -> G*x convolves the input x with the Green's function
Op = LinOpLippmannUi(G,f);                   % Operator (I - G*diag(f))   --> Op*x computes (I - G*diag(f))*x

% - Compute the total field utot = (I - G*diag(f))^{-1} uin
%   Done via the resolution of
%         utot = arg min_u ||Op*u - uin||^2   
%    <=>  Op'*Op*utot = Op'*uin
% Some precomputations
A = Op'*Op;    
b = Op'*uin;
% Instanciate a conjugate gradient to solve A*utot=b            
CG = OptiConjGrad(A,b);                    % Conjugate gradient algorithm
CG.maxiter = 50;                           % maximal number of iterations
CG.ItUpOut = round(CG.maxiter/10);         % Display iterations every maxiter/10 iterations
CG.CvOp = TestCvgStepRelative(1e-8);       % Stopping criteria with difference between two iterates
CG.run(uin);                               % Run the CG initialized with uin 
utot = CG.xopt;                            % Get the solution

%% Display
figure;
subplot(1,2,1);imagesc(abs(uin)); axis image; axis off; title('Amplitude uin'); colorbar;
subplot(1,2,2);imagesc(angle(uin)); axis image; axis off; title('Phase uin'); colorbar;
figure; 
imagesc(f); axis image; axis off; title('Scattering potential'); colorbar;
figure;
subplot(1,2,1);imagesc(abs(utot)); axis image; axis off; title('Amplitude utot'); colorbar;
subplot(1,2,2);imagesc(angle(utot)); axis image; axis off; title('Phase utot'); colorbar;