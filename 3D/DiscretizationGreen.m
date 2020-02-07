%--------------------------------------------------------------------------
% This script reproduce Figure 3 in
%
% [1] Three-Dimensional Optical Diffraction Tomographywith Lippmann-Schwinger Model,
%     IEEE Transactions on Computational Imaging (2020),
%     T-a. Pham, E. Soubies, A. Ayoub, J. Lim, D. Psaltis, and M. Unser.
% 
% Copyright 2020
%   Thanh-an Pham (thanh-an.pham@epfl.ch)
%   Emmanuel Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
clear, close
addpath('Utils/');

%% Parameters
p=4;                   % padding factor
nlist=10:10:100;       % number of points per side
L=1;                   % lenght of Omega [0,L]
kb=1.5;                % wavevector
sig=0.04;              % parameter of the Gaussian source1

%% Gaussian source distribution
v_fun   = @(r) 1/((sig)^3*(2*pi)^1.5)*exp(-r.^2./(2*(sig)^2));
sol_fun = @(r) -1./(4*pi*r).*exp(-sig^2*kb^2/2).*(real(exp(-1i*kb*r).*erf_((2*sig^2*1i*kb-2*r)/(2*sqrt(2*sig^2)),2e4))-1i*sin(kb*r));

%% Compute Green's kernel with discretization proposed in [1] (Theorem III.1 and Proposition III.2)
% Called : Truncated Green's function
it=1;
for n=nlist
    % -- Some precomputations
    h=L/n;
    fourierStep=2/(2*h*n*p);
    fr=fftshift(linspace(-1/(2*h),1/(2*h)-fourierStep,n*p));
    % -- Get the Fourier of the Green's kernel
    gc=zeros([2*n,2*n,2*n]);
    cc=[0:n-1,n*p-n:n*p-1];  % indexes in the cropped region (fftshift)
    [x1,x2,x3]=meshgrid(cc,cc,cc);
    for s1=0:p/2-1
        fr1=fr(s1+1:p/2:end);
        for s2=0:p/2-1
            fr2=fr(s2+1:p/2:end);
            for s3=0:p/2-1
                fr3=fr(s3+1:p/2:end);
                [k1,k2,k3]=meshgrid(fr1,fr2,fr3);
                normk=2*pi*sqrt(k1.^2+k2.^2+k3.^2);
                gchat=(-1 + exp(1i*sqrt(3)*L*kb)*(cos(sqrt(3)*L*normk) - 1i * kb * sqrt(3)*L *sinc(sqrt(3)*L*normk/pi)))./(kb^2 - abs(normk).^2);
                gchat(normk==kb)=1i*(sqrt(3)*L/(2*kb) -exp(1i*sqrt(3)*L*kb)*sin(sqrt(3)*L*kb)/(2*kb^2));                                           % singularity at kb
                gc=gc+fftshift(ifftn(gchat).*exp(2*1i*pi/(n*p)*(s1.*x1+s2.*x2+s3.*x3))/(p/2)^3);
            end
        end
    end
    gchat_vainiko{it}=fftshift(fftn(fftshift(gc)));
    it=it+1;
end

%% Compute Green's kernel by "cropping" the singularity 
% Called : Naive crop
it=1;
for n=nlist
    h=L/n;
    x=linspace(-L,L-h,2*n);
    [X,Y,Z]=meshgrid(x,x,x);
    R=sqrt(X.^2+Y.^2+Z.^2);
    thres=eps;
    gc = (exp(1i*kb*R)./(4*pi*R))/n^3;
    gc(isinf(gc))=max(gc(~isinf(gc)));
    gchat_crop{it}= fftn(fftshift(gc));
    it=it+1;
end

%% Compare the errors after convolution with Green's kernel
it=1;
for n=nlist
    h=L/n;
    % - Discretization points of Omega and discrete source distribution
    x=linspace(-L/2,L/2-h,n);
    [X,Y,Z]=meshgrid(x,x,x);
    R=sqrt(X.^2+Y.^2+Z.^2);
    v = v_fun(R);
    % - Analytical solution
    sol = sol_fun(R);
    sol(R==0)=sol_fun(eps);   % extends by continuity at zero
    % - Convolution with Green's function
    vp = padarray(v,size(v)/2,'both');
    % Vainiko
    Gv_Trunc=fftshift(ifftn(fftshift(gchat_vainiko{it}).*fftn(fftshift(vp))));Gv_Trunc=Gv_Trunc(n/2+1:end-n/2,n/2+1:end-n/2,n/2+1:end-n/2);
    % Crop
    Gv_crop=(ifftn(gchat_crop{it}.*fftn((vp))));Gv_crop=Gv_crop(n/2+1:end-n/2,n/2+1:end-n/2,n/2+1:end-n/2);
    % - Compute the error
    tt=~isnan(sol(:)); 
    err_Vainiko(it)=norm(Gv_Trunc(tt)-sol(tt))/norm(sol(tt));
    err_Crop(it)=norm(Gv_crop(tt)-sol(tt))/norm(sol(tt));
    disp(['-- n = ',num2str(n)]);
    disp(['Error with Truncated Green''s: ',num2str(err_Vainiko(it))]);
    disp(['Error with Naive crop: ',num2str(err_Crop(it))]);
    it=it+1;
end

%% Plots and error calculation
figure;
subplot(2,3,1);
imdisp(abs(sol(:,:,n/2+1)),'Analytical Solution (abs)',0);colormap default;
subplot(2,3,4);
imdisp(abs(sol(:,:,n/2+1)),'Analytical Solution (abs)',0);colormap default;
subplot(2,3,2);
imdisp(abs(Gv_Trunc(:,:,n/2+1)),'Truncated Green''s (abs)',0);colormap default;
subplot(2,3,3);
imdisp(abs(Gv_Trunc(:,:,n/2+1)-sol(:,:,n/2+1)),'Error',0);colormap default;colorbar;
subplot(2,3,5);
imdisp(abs(Gv_crop(:,:,n/2+1)),'Naive crop (abs)',0);colormap default;
subplot(2,3,6);
imdisp(abs(Gv_crop(:,:,n/2+1)-sol(:,:,n/2+1)),'Error',0);colormap default;colorbar;

figure; hold all;set(gca,'FontSize',14);set(gca,'yscale','log');grid;
semilogy(nlist,err_Crop,'-x','Linewidth',1.5);
semilogy(nlist,err_Vainiko,'-o','Linewidth',1.5);
legend('Naive crop','Truncated Green''s Function');
xlabel('n');ylabel('Relative error');