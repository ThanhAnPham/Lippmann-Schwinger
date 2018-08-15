function fhat_roi = radonRoiEst2D(gem, gob, thetas,thet_jump,doInterp,isnegative, varargin)

GB = false;
fact = 1;
doGlobUnwrap = false;

if nargin>=6
    doPlot = varargin{1};
    if nargin>=7
        fact = varargin{2};
        if nargin>=8
            if ~isempty(varargin{3})
                %attempt
                simP = varargin{3}{1};
                GB = simP.GaussianBeam;
                roi_x = varargin{3}{2};
            end
            if nargin>=9
                doGlobUnwrap = varargin{4};
            end
        end
    end
else
    doPlot = false;
end

if isempty(doInterp)
    doInterp = false;
end

[Nx, Ntheta] = size(gem);

angs = rad2deg(thetas(end:-1:1))-90; % for iradon in degrees

phaseData = zeros(Nx, Ntheta);
%figure(6);

for ind_theta = 1:Ntheta
    imgem = angle(gem(:,ind_theta));
    imgob = angle(gob(:,ind_theta));
    
    %%% Unwrap the phase
    %phiem1 = unwrap(imgem,thet_jump);
    %phiob1 = unwrap(imgob,thet_jump);
    
    if round(1e3*thetas(ind_theta))/1e3 > 0
        phiem1 = unwrap(imgem,thet_jump);
        phiob1 = unwrap(imgob,thet_jump);
    elseif round(1e3*thetas(ind_theta))/1e3~=0
        
        phiem1 = flipud(unwrap(flipud(imgem),thet_jump));
        phiob1 = flipud(unwrap(flipud(imgob),thet_jump));
        
    end
    
    phi = (-1)^isnegative*(phiob1 - phiem1);
    phi(isnan(phi)) = 0;
    
    if GB
        angx = simP.thetas(Ntheta);
        zz = simP.dz*simP.Nz/2;
        xx = simP.dx*(-simP.Nx/2+1:simP.Nx/2)';
        Bx = simP.sigmaBeam*simP.Lx; % beam scale (m)
        Bx = Bx/cos(angx); % projection of the propagation plane to XY plane
        Bz = simP.sigmaBeam*simP.Lz; Bz = -Bz/sin(angx);%"-" because from (z',x') -> (z,x)
        Bz = -Bz/sin(angx);%"-" because from (z',x') -> (z,x)
        ain = exp(-(xx(roi_x)./Bx + zz./Bz).^2); % illumination beam amplitude
    else
        ain = 1;
    end
    
    phaseData(:,ind_theta) = ain.*phi;%attempt
    
    if doPlot || ind_theta==Ntheta %disp at the end
        subplot(2,2,[1,3]);
        plot(imgob);hold on;
        plot(imgem);
        plot(phiob1);
        plot(phiem1);
        plot(phi);hold off;
        title(sprintf('Angle %1.2f',thetas(ind_theta)*180/pi));
        legend('wrapped ob','wrapped em',...
            'unwrapped ob','unwrapped em','diff');
        ylabel('phase(x)');xlabel('x');
        subplot(224);
        imagesc(phaseData);colorbar;
        pause(0.001);
    end
end

[fhat_roi,phi] = FBProj(phaseData,angs,doInterp,doGlobUnwrap);

subplot(224);
imagesc(phi);colorbar;pause(0.05)
fhat_roi = fhat_roi/fact;
end