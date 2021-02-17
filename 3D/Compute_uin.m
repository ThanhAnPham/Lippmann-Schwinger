%% Figure 5
clear
close all;
useGPU(0);
simP.lambda0 = 0.532;% wavelength
simP.Nx = 300;% Number of voxels X
simP.Ny = 300;% Number of voxels Y
simP.Nz = 96;% Number of voxels Z
simP.dx = 0.0992561403508772;% Same than real data, step size X
simP.dy = simP.dx;% step size Y
simP.dz = simP.dx;% step size Z
simP.Lx = simP.Nx*simP.dx;%6*simP.lambda0;% Physical length X
simP.Ly = simP.Ny*simP.dy;%6*simP.lambda0;% Physical length Y
simP.Lz = simP.Nz*simP.dz;% Physical length Z
simP.n0 = 1.338;% Background refractive index
simP.thetas = [0.558053163587562, 0.24846127635636]; %reproduce paper
         
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
uin0 = uin(:,:,end);

dProp = -(simP.Nz/2-1)*simP.dz;

uin_pGT = uin(:,:,end + round(dProp/simP.dz));
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
Noi = 96;
S = LinOpSelectorPatch(size(uin_AS),1 + size(uin_AS)/2 - Noi/2,size(uin_AS)/2 + Noi/2);
figure(1);set(gcf,'Units','normalized','Position',[0,0,.5,1]);
subplot(321);
imagesc(abs(S*uin_pGT));axis image;colorbar;title('Ground truth (amplitude)');
subplot(322);
imagesc(angle(S*uin_pGT));axis image;colorbar;title('Ground truth (phase)');
subplot(323);
imagesc(abs(S*uin_AS));axis image;colorbar;title('AS propagated (amplitude)');
subplot(325);
imagesc(abs(S*uin_AStt));axis image;colorbar;title('tilted AS propagated (amplitude)');
subplot(324);
imagesc(S*abs(uin_AS - uin_pGT));axis image;colorbar;title('diff AS-GT propagated (amplitude)');
subplot(326);
imagesc(S*abs(uin_AStt - uin_pGT));axis image;colorbar;title('diff tilted AS-GT propagated (amplitude)');