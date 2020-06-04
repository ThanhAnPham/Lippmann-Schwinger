%% Scattered field from a bead using Lippmann-Schwinger model
clear;
close all;
useGPU(0);

%% Create two beads
simP.lambda0 = 0.406;% wavelength
simP.Nx = 96;% Number of voxels X
simP.Ny = 96;% Number of voxels Y
simP.Nz = 96;% Number of voxels Z
simP.Lx = 6*simP.lambda0;% Physical length X
simP.Ly = 6*simP.lambda0;% Physical length Y
simP.Lz = 6*simP.lambda0;% Physical length Z
simP.n0 = 1.333;% Background refractive index
radii = [simP.lambda0,simP.lambda0];% Beads radii
centers = cat(1,[simP.Nx,simP.Ny,simP.Nz/2]/2,...
    [simP.Nx,simP.Ny,simP.Nz*3/2]/2);% Beads center
simP.dn = [0.1,0.1];% refractive index (diff) of the beads

simP.dx = simP.Lx/simP.Nx;% step size X
simP.dy = simP.Ly/simP.Ny;% step size Y
simP.dz = simP.Lz/simP.Nz;% step size Z
simP.k0 = 2*pi/simP.lambda0;% Freespace wavenumber
simP.k = simP.k0*simP.n0;% Wavenumber
assert(simP.dx==simP.dy & simP.dy==simP.dz,'Step sizes must be equal');

sz = [simP.Nx,simP.Ny,simP.Nz];

n = simP.n0*ones(sz);%Refractive index distribution

for kk = 1:length(radii)
    Vectx = (0:simP.Nx-1) - floor(centers(kk,1));
    Vecty = (0:simP.Ny-1) - floor(centers(kk,2));
    Vectz = (0:simP.Nz-1) - floor(centers(kk,3));
    [X,Y,Z] = meshgrid(Vectx,Vecty,Vectz);
    
    R = sqrt((X*simP.dx).^2 + (Y*simP.dy).^2 + (Z*simP.dz).^2);
    posObj = R < radii(kk);
    cn = zeros(sz);
    cn(posObj) = simP.dn(kk);
    n = n + cn;
end

f = gpuCpuConverter((simP.k*simP.dx)^2*((n/simP.n0).^2 - 1));% Scattering Potential
%% Build the incident field (ideal plane wave)
simP.thetas = [pi/8,0];%Angle of incidence (rad.)

Vectx = (1:simP.Nx) - simP.Nx/2;
Vecty = (1:simP.Ny) - simP.Ny/2;
Vectz = (1:simP.Nz) - simP.Nz/2;
[X,Y,Z] = meshgrid(Vectx*simP.dx,Vecty*simP.dy,Vectz*simP.dz);


simP.kvec = simP.k*[sin(simP.thetas(1));sin(simP.thetas(2))];%Wave vector
simP.kvec(3) = (-1)^((abs(simP.thetas(1)) > pi/2) + (abs(simP.thetas(2)) > pi/2))*sqrt(simP.k^2 - sum(simP.kvec(1:2).^2));
uin = gpuCpuConverter(exp(1i*(simP.kvec(1)*X + simP.kvec(2)*Y + simP.kvec(3)*Z)));

%% Set the sensors position (outside the volume)
SP.Nmeas = [150,150];%Number of sensors
%XY positions: origin at the center of the sample
[SP.x,SP.y] = meshgrid((simP.dx:simP.dx:simP.dx*SP.Nmeas(1)) - simP.dx*SP.Nmeas(1)/2,...
    (simP.dy:simP.dy:simP.dy*SP.Nmeas(2)) - simP.dy*SP.Nmeas(2)/2);

SP.y = SP.y(:); SP.x = SP.x(:);

%one pixel outside of the sample.
SP.z = simP.dz*(sz(3)/2+1)*ones(length(SP.y),1);

%% Build the forward model based on Lippmann-Schwinger

Niter = 200;% # iterations to compute the forward model
xtol = 1e-9;% tolerance for convergence criteria

% Green function in the volume
G = FreeSpaceKernel3D(simP.k*simP.dx,sz,'Vainikko');
% Green function to get the scattered field at the sensor positions
Gtil = Vainikkotilde(simP.k,sz,SP,[simP.dx,simP.dy,simP.dz],ones(1,3));

H = OpLipp3D(G,Gtil,LinOpIdentity(sz),uin,Niter,xtol,[],false);

%% Compute the scattered field on the sensors positions
tic
y = H*f;
toc
%% Orthoviews of 3D total field
Orthoviews(abs(H.u));
%% (Optional) Repropagate the scattered field at the (center + shiftZ)
shiftZ = simP.Nz/4;
par.distuM = -SP.z(end) + shiftZ*simP.dz;%-(simP.Nz)/2*simP.dz;%par.distUin = par.distuM;

Vectx = (1:SP.Nmeas(1)) - SP.Nmeas(1)/2;
Vecty = (1:SP.Nmeas(2)) - SP.Nmeas(2)/2;
[X,Y] = meshgrid(Vectx*simP.dx,Vecty*simP.dy);

%incident field at the repropagated z-position.
uem = exp(1i*(simP.kvec(1)*X + simP.kvec(2)*Y + simP.kvec(3)*(1 + shiftZ)*simP.dz));

%Propagator
P = LinOpTiltedPropagator(simP.lambda0,simP.n0,par.distuM,simP.dx,SP.Nmeas,[0,0],'AS');

yprop = P*y;

figure(1);
subplot(221);
imagesc(abs(y));axis image;colorbar;title('scattered field');
subplot(222);
imagesc(abs(yprop));axis image;colorbar;title('centered scattered field');
subplot(223);
imagesc(abs(yprop+uem));axis image;colorbar;title('centered total field');
subplot(224);
imagesc(angle((yprop+uem)./uem));axis image;colorbar;title('centered relative phase');

fprintf('Approximate expected relative phase %1.2f\n',simP.k0*dot(simP.dn,2*radii));