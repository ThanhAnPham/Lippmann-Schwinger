function [out,dphi] = FBProj(dphi,angs,varargin)
%dphi : phase difference
%thetas : in degrees


doInterp =  0;
doGlobUnwrap = 0;

if nargin > 2
    doInterp = varargin{1};
    if nargin > 3
        doGlobUnwrap = varargin{2};
        
    end
end


Ntheta = size(dphi,2);
dphi(dphi < 0) = 0;

Nx = size(dphi,1);

if doInterp
    thetas = deg2rad(angs + 90);
    x = 1:Nx;
    % Below done by Kamilov, but better without in 2D for BPM
    t = x; % coordinate orthogonal to the propagation direction
    Nt = Nx; % number of samples in the ortho direction
    
    [Theta, X] = meshgrid(thetas, x); % Meshgrid for interp1
    rdata2 = zeros(Nt, Ntheta); % initialize the cosine transformed data
    
    %%% Transform (x, theta) to (t, theta)
    for ind_theta = 1:Ntheta
        theta = thetas(ind_theta);
        xt = t/cos(theta);%t(ind_t)/cos(theta); % = t/cos(theta)
        rdata2(:,ind_theta) =...
            interp2(Theta,X,phaseData,theta,xt,'linear',0);
    end
end

%dphi = medfilt2(dphi, [3,3]);%smooth

if doGlobUnwrap
    dphi = unwrap(dphi,[],2);
    dphi = max(0,dphi);
    subplot(224);
    imagesc(dphi);
end

out = iradon(dphi,angs,'linear','han',Nx);
out(out < 0) = 0;
end