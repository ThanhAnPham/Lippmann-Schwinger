%% Main Rytov FBP

if hasandis(simP,'realData')
    uM.rytovFBP = uM.GT;
else
    uM.rytovFBP = uM.curr_y(1 + end/2:end,:);%    uM_true(1 + end/2:end,:);
end
Muin = zeros(simP.Nx,simP.Ntheta);
if ~exist('M','var') || ~isa(M,'LinOpSelector') %early simu prob
    M = setM(simP.Nx,simP.Nz,simP.dx,simP.dz,SP);
end

if exist('uin','var') && ~isempty(uin)
    if ~ismatrix(uin)
        for kk = 1:simP.Ntheta
            tmp = M*uin(:,:,kk);
            Muin(:,kk) = tmp(1+end/2:end);
        end
    else
        Muin = uin;
    end
else
    Muin = uem(1 + end/2:end,:);
end

clear tmp;

%%
if ~exist('log_imag_ur2D','var')
    if hasandis(simP,'realData')
        par.rytovFBP.distUin = -simP.Lz/2;
        par.rytovFBP.distuM = par.rytovFBP.distUin;
    else
        par.rytovFBP.distUin = -simP.Lz/2;
        par.rytovFBP.distuM = par.rytovFBP.distUin;
    end
    [Muin,uM.rytovFBP] = setFields(Muin,uM.rytovFBP,par.rytovFBP,simP,'rytovFBP'); %,-par.rytovFBP.Lzext/2,-par.rytovFBP.Lzext/2);
    
    [~,loguM] = unwrapMeas(Muin,uM.rytovFBP,simP.thetas,false);%,{utot, uin});
    
else
    loguM = 1i*log_imag_ur2D;
end

%par.rytovFBP.max_scan_theta = deg2rad(100);
%par.rytovFBP = setScanningAngle(par.rytovFBP,simP);
%par.rytovFBP.thetas = simP.thetas(par.rytovFBP.curr_thetas);
uM.rytovFBP = uM.rytovFBP(:,par.curr_thetas);
Muin = Muin(:,par.curr_thetas);
loguM = loguM(:,par.curr_thetas);

if exist('k_x','var')
    [n_rytovFBP,f_rytovFBP,ff_rytovFBP,ff_mask] = FBProp2D('rytov',uM.rytovFBP,Muin,imag(loguM),simP,par.thetas,k_x(par.curr_thetas),k_y(par.curr_thetas));
else
    [n_rytovFBP,f_rytovFBP,ff_rytovFBP,ff_mask] = FBProp2D('rytov',uM.rytovFBP,Muin,imag(loguM),simP,par.thetas);
end
par.rytovFBP.dx = simP.dx; par.rytovFBP.dz = simP.dx;
par.rytovFBP.Lxext = simP.dx*size(n_rytovFBP,1); par.rytovFBP.Lzext = simP.dx*size(n_rytovFBP,2);
par.rytovFBP.Lx_vec = simP.dx:simP.dx:par.rytovFBP.Lxext; par.rytovFBP.Lz_vec = par.rytovFBP.Lx_vec;
%%
if hasandis(simP,'realData') && par.negative
    %n_hat.rytovFBP = n_rytovFBP;
    n_hat.rytovFBP = min(n_rytovFBP,simP.n0);%n_rytovFBP - min(n_rytovFBP(:)) + simP.n0;%
else
    n_hat.rytovFBP = max(n_rytovFBP,simP.n0);% circshift(n_rytovFBP,[0,simP.Nz/2]);%
end
%% n_hat.rytovFBP(n_hat.rytovFBP <= 1.538) = simP.n0;
figure(13);
subplot(121);
if isfield(simP,'Lx_vec')
    imagesc(par.rytovFBP.Lx_vec,par.rytovFBP.Lz_vec,n_hat.rytovFBP);
else
    imagesc(n_hat.rytovFBP);
end
axis image;colorbar;
%caxis([simP.n0,simP.n0+simP.dn]);
subplot(122);imagesc(real(f_rytovFBP));axis image;colorbar;

%fprintf('SNR : %1.2f dB\n',snr(n,n - n_hat.rytovFBP));