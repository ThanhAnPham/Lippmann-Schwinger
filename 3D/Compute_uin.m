%% Figure 5
clear
close all;
useGPU(0);
simP.lambda0 = 0.406;% wavelength
simP.Nx = 256;% Number of voxels X
simP.Ny = 256;% Number of voxels Y
simP.Nz = 128;% Number of voxels Z
simP.Lx = 6*simP.lambda0;% Physical length X
simP.Ly = 6*simP.lambda0;% Physical length Y
simP.Lz = 3*simP.lambda0;% Physical length Z
simP.n0 = 1.333;% Background refractive index
simP.thetas = [pi/8,pi/4];

simP.dx = simP.Lx/simP.Nx;% step size X
simP.dy = simP.Ly/simP.Ny;% step size Y
simP.dz = simP.Lz/simP.Nz;% step size Z
simP.k0 = 2*pi/simP.lambda0;% Freespace wavenumber
simP.k = simP.k0*simP.n0;% Wavenumber
assert(simP.dx==simP.dy & simP.dy==simP.dz,'Step sizes must be equal');



simP.kvec = simP.k*[sin(simP.thetas(1));sin(simP.thetas(2))];%Wave vector
simP.kvec(3) = (-1)^((abs(simP.thetas(1)) > pi/2) + (abs(simP.thetas(2)) > pi/2))*sqrt(simP.k^2 - sum(simP.kvec(1:2).^2));

Vectx = (1:simP.Nx) - simP.Nx/2;
Vecty = (1:simP.Ny) - simP.Ny/2;
Vectz = (1:simP.Nz) - simP.Nz/2;
[X,Y,Z] = meshgrid(Vectx*simP.dx,Vecty*simP.dy,Vectz*simP.dz);

sz = [simP.Nx,simP.Ny,simP.Nz];

uin = gpuCpuConverter(exp(1i*(simP.kvec(1)*X + simP.kvec(2)*Y + simP.kvec(3)*Z)));
uin0 = uin(:,:,1);

dProp = (simP.Nz-1)*simP.dz;

uin_pGT = uin(:,:,1 + round(dProp/simP.dz));
%% Angular spectrum

P = LinOpTiltedPropagator(simP.lambda0,...
    simP.n0,dProp,simP.dx,sz(1:2),[0,0],'AS');
uin_AS =P*uin0;
%% Angular spectrum with tilt transfer

[X2D,Y2D] = meshgrid(Vectx*simP.dx,Vecty*simP.dy);
tiltFactorB = exp(1i*simP.k*(sin(simP.thetas(1))*X2D + sin(simP.thetas(2))*Y2D));
curr_uinf = fft2(uin0./tiltFactorB);
P = LinOpTiltedPropagator(simP.lambda0,...
    simP.n0,dProp,simP.dx,sz(1:2),-rad2deg(simP.thetas),'AS');

uin_AStt = (ifft2(P.ephi.*curr_uinf)).*tiltFactorB;%*exp(1i*simP.k*simP.distUin(kk)*2);
%% Display

figure(1);set(gcf,'Units','normalized','Position',[0,0,.5,1]);
subplot(321);
imagesc(abs(uin_pGT));axis image;colorbar;title('Ground truth (amplitude)');
subplot(322);
imagesc(angle(uin_pGT));axis image;colorbar;title('Ground truth (phase)');
subplot(323);
imagesc(abs(uin_AS));axis image;colorbar;title('AS propagated (amplitude)');
subplot(325);
imagesc(abs(uin_AStt));axis image;colorbar;title('tilted AS propagated (amplitude)');
subplot(324);
imagesc(abs(uin_AS - uin_pGT));axis image;colorbar;title('diff AS-GT propagated (amplitude)');
subplot(326);
imagesc(abs(uin_AStt - uin_pGT));axis image;colorbar;title('diff tilted AS-GT propagated (amplitude)');