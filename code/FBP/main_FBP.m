%% MAIN FBP
%Does FBP on phase difference between wavefield with and without object
%=>fhat_roi ~ simP.k0*dn*simP.dz
clear roi

len_oi = simP.Nz*3/4;%avoid border effect of backpropagation and stuff
doInterp = false;

par.fbp.Lx = simP.Lx; par.fbp.Lz = simP.Lz;

curr_thetas = find(abs(par.thetas) <= deg2rad(45));
Ntheta = length(curr_thetas);

roix = SP.z==simP.Lz;

roi.x = (1:len_oi) + simP.Nx/2 - len_oi/2;
roi.z = (1:len_oi) + simP.Nz/2 - len_oi/2;

%Diffraction step
siz = nnz(roix);

dkx = 2*pi/simP.Lx;

Kxx = (dkx * [0:(siz/2-1), (-siz/2):-1]').^2;
k0 = 2*pi/simP.lambda0;
k = k0*simP.n0;

dphi = real(Kxx./(k+ sqrt((k)^2 - Kxx))); % diffraction phase factor

if ~exist('uem','var') || any(isnan(uem(:)))
    uin_fbp = setUin(par,simP,'fbp');
    if ~ismatrix(uin_fbp)
        uin_fbp = squeeze(uin_fbp(:,end,:));
    end
else
    uin_fbp = uem;
end
%Back propagation to the center for unwrapping
Zplane = 0.5;%*(1 - len_oi/simP.Nz);

uem_oi = ifft(fft(uin_fbp(roix,curr_thetas))...
    .*repmat(exp(-1i*dphi*(-Zplane)*simP.Lz),1,Ntheta));
uob_oi = ifft(fft(uM.curr_y(roix,curr_thetas))...
    .*repmat(exp(-1i*dphi*(-Zplane)*simP.Lz),1,Ntheta));

uem_oi = uem_oi(roi.x,:);
uob_oi = uob_oi(roi.x,:);


figure(8);subplot(222);
imagesc(radon((n - simP.n0)*k0,simP.thetas*180/pi-90));colorbar;title('GT sinogram');%'GT' sinogram

dnhat_roi = radonRoiEst2D(uem_oi, uob_oi, par.thetas(curr_thetas),...
    k0*0.2*simP.dx,doInterp,par.negative,0,1,[],false);
dnhat_roi = (-1)^par.negative*dnhat_roi/(k0*simP.dx);

n_hat.fbp = simP.n0*ones(simP.Nx,simP.Nz);
n_hat.fbp(roi.x,roi.z) = simP.n0 + dnhat_roi;

par.fbp.dx = simP.dx; par.fbp.dz = par.fbp.dx;

clear uob_oi uem_oi dphi Kxx roi k k0 dnhat_roi Zplane len_oi doInterp dkx roix siz Ntheta curr_thetas
%%
figure(13);
imagesc(simP.Lx_vec,simP.Lz_vec,n_hat.fbp);%n_hat);
%if max(n_hat.fbp(:)) > 1.5*(simP.n0 + simP.dn)
%caxis([simP.n0,simP.n0 + max(simP.dn)]);%might help to see
%end
if par.negative
    colormap(flipud(gray));
else
    colormap gray;
end
axis image;colorbar;
fprintf('SNR : %1.2f\n',snr(n,n - n_hat.fbp));
clear uin_fbp