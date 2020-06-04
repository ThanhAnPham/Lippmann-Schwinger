function [uin_out,y] = setFields3D(uin,y,param,simP,model,varargin)
%SETFIELDS [uin_out,y] = setFields3D(uin,y,param,simP,model,varargin)
%   Propagates uin /y at a given distance

mode = 'both';%determines which field(s) is/are propagated, 'both','uin','uM'
dx = simP.dx;
dy = simP.dy;
doSpecialCase = 0;

if hasandis(param,'reflect')
    param.distUin = -param.distUin;
    param.distuM = -param.distuM;
end

if nargin > 5
    if ~isempty(varargin{1})
        mode = varargin{1};
    end
    if nargin > 6
        if ~isempty(varargin{2})
            dx = varargin{2}(1);
            dy = varargin{2}(2);
        end
    end
    if nargin > 7
        if ~isempty(varargin{3})
            doSpecialCase = varargin{3};
        end
    end
end
switch mode
    case 'uin'
        param.distuM = 0;
    case 'uM'
        param.distUin = 0;
    case 'uem'
        param.distUin = param.distuM;
        param.distuM = 0;
end

if isempty(uin)
    Nx = size(y,1);
    Ny = size(y,2);
else
    Nx = size(uin,1);
    Ny = size(uin,2);
end
padfact = ones(1,2);
dkx = 2*pi/(padfact(1)*Nx*dx);
dky = 2*pi/(padfact(2)*Ny*dy);
Nx = padfact(1)*Nx;
Ny = padfact(2)*Ny;

if mod(Nx,2)==0
    Kxx = (dkx * [0:(Nx/2-1), (-Nx/2):-1]').^2;
else
    Kxx = (dkx * [0:(Nx-1)/2, (-(Nx-1)/2):-1]').^2;
end

if mod(Ny,2)==0
    Kyy = (dky * [0:(Ny/2-1), (-Ny/2):-1]').^2;
else
    Kyy = (dky * [0:(Ny-1)/2, (-(Ny-1)/2):-1]').^2;
end

[Kyy, Kxx] = meshgrid(Kyy, Kxx);

k = 2*pi/(simP.lambda0)*simP.n0;

dphi = real(k - sqrt(k^2 - Kxx - Kyy)); % diffraction phase factor


switch model
    case {'BPM','bpm','bpmi','rytovFBP','rytFBP','lfr'} %propagates uin -0.5 Lz and uM 0.5 Lz
        
        if length(param.distUin) > 1
            error('Number of propagations of Uin too high, should be 1 for BPM and cie\n');
        end
        
        if param.distuM ~=0
            
            if doSpecialCase
                
                sz_y = size(y); sz_y = sz_y(1:2);
                pfact = 1;
                Padder = LinOpSelectorPatch(pfact*sz_y,1 + (pfact-1)*sz_y/2,(pfact+1)*sz_y/2)';
                PadderExt = LinOpSelectorPatch([pfact*sz_y],[1 + (pfact-1)*sz_y/2],[(pfact+1)*sz_y/2])';
                
                xx = (dx:dx:(sz_y(1)*dx)) - (sz_y(1)*dx)/2 - dx;
                yy = (dy:dy:(sz_y(2)*dy)) - (sz_y(2)*dy)/2 - dx;
                [YYB,XXB] = meshgrid(xx,yy);
                useGPU(0);
                
                yout = zeros(size(y));
                %%
                for ll = 1:param.Ntheta
                    %tic
                    tiltFactorB = exp(1i*k*(sin(param.thetas(ll,1))*YYB + sin(param.thetas(ll,2))*XXB));
                    curr_yf = fft2(Padder*(y(:,:,param.curr_ind_thetas(ll))./tiltFactorB));
                    cthetas = -rad2deg(param.thetas(ll,:))';%asin(simP.k*sin(simP.thetas(ll,[2,1]))*simP.lambda0/(2*pi).*[simP.Lx,simP.Ly])';
                    P = LinOpTiltedPropagator(simP.lambda0,...
                        simP.n0,0,dx,pfact*sz_y,cthetas([2,1]),'AS');
                    lf = @(z) abs(z*sqrt(P.mu)./(sqrt((  P.n0/ P.lambda)^2- (P.mu + P.mv)))) <= P.sizein(1)*P.dxy...
                        & abs(z*sqrt(P.mv)./(sqrt((  P.n0/ P.lambda)^2- (P.mu + P.mv)))) <= P.sizein(2)*P.dxy;
                    ephi = arrayfun(@(x) P.mod.*exp(P.phase*x).*lf(x),-param.distuM,'UniformOutput',false);
                    ephi = cat(3,ephi{:});
                    
                    yout(:,:,ll) = (PadderExt'*ifft2(ephi.*curr_yf)).*tiltFactorB;
                    
                    %toc
                end
                y = yout;
            else
                y = ifft2(fft2(padarray(y,([padfact,1] - 1).*[Nx,Ny,0]/2,0,'both')).*repmat(exp(-1i*dphi*param.distuM),...
                    [ones(1,2 + (ndims(y)==4)),size(y,max(3,ndims(y)))]));
            end
            %y = exp(1i*k*param.distuM)*y;
        end
        if any(param.distUin ~=0)
             if doSpecialCase
                
                sz_y = size(y); sz_y = sz_y(1:2);
                pfact = 1;
                Padder = LinOpSelectorPatch(pfact*sz_y,1 + (pfact-1)*sz_y/2,(pfact+1)*sz_y/2)';
                PadderExt = LinOpSelectorPatch([pfact*sz_y],[1 + (pfact-1)*sz_y/2],[(pfact+1)*sz_y/2])';
                
                xx = (dx:dx:(sz_y(1)*dx)) - (sz_y(1)*dx)/2 - dx;
                yy = (dy:dy:(sz_y(2)*dy)) - (sz_y(2)*dy)/2 - dx;
                [YYB,XXB] = meshgrid(xx,yy);
                useGPU(0);
                
                uin_out = zeros(size(y));
                %%
                for ll = 1:param.Ntheta
                    %tic
                    tiltFactorB = exp(1i*k*(sin(param.thetas(ll,1))*YYB + sin(param.thetas(ll,2))*XXB));
                    curr_yf = fft2(Padder*(uin(:,:,param.curr_ind_thetas(ll))./tiltFactorB));
                    cthetas = -rad2deg(param.thetas(ll,:))';%asin(simP.k*sin(simP.thetas(ll,[2,1]))*simP.lambda0/(2*pi).*[simP.Lx,simP.Ly])';
                    P = LinOpTiltedPropagator(simP.lambda0,...
                        simP.n0,0,dx,pfact*sz_y,cthetas([2,1]),'AS');
                    lf = @(z) abs(z*sqrt(P.mu)./(sqrt((  P.n0/ P.lambda)^2- (P.mu + P.mv)))) <= P.sizein(1)*P.dxy...
                        & abs(z*sqrt(P.mv)./(sqrt((  P.n0/ P.lambda)^2- (P.mu + P.mv)))) <= P.sizein(2)*P.dxy;
                    ephi = arrayfun(@(x) P.mod.*exp(P.phase*x).*lf(x),-param.distUin,'UniformOutput',false);
                    ephi = cat(3,ephi{:});
                    
                    uin_out(:,:,ll) = (PadderExt'*ifft2(ephi.*curr_yf)).*tiltFactorB;
                    
                    %toc
                end
             else
            uin_out = ifft2(fft2(padarray(uin,([padfact,1] - 1).*[Nx,Ny,0]/2,0,'both')).*repmat(exp(-1i*dphi*min(param.distUin)),...
                [ones(1,2 + (ndims(uin)==4)),size(uin,max(3,ndims(uin)))]));
            %uin_out = exp(1i*k*param.distUin)*uin_out;
             end
        else
            uin_out = uin;
        end
    otherwise
        if param.distuM ~=0
            y = ifft2(fft2(y).*repmat(exp(-1i*dphi*param.distuM),[ones(1,2 + (ndims(y)==4)),size(y,max(3,ndims(y)))]));
            y = exp(1i*k*param.distuM)*y;
        end
        if strcmp(mode,'both') || strcmp(mode,'uin') || strcmp(mode,'uem')
            if any(param.distUin~=0)
                if length(param.distUin)~=param.Nz && ~strcmp(mode,'uem') && ~param.discreteGtild
                    error('param.distUin needs to be of the same length than Nz');
                end
                if ndims(uin)>2
                    Ntheta = size(uin,ndims(uin));
                else
                    Ntheta = 1;
                end
                
                Nz = length(param.distUin);
                repel = repmat({':'},1,ndims(uin) - 1 + 1*(Ntheta==1));
                expkz = permute(repmat(exp(1i*k*reshape(param.distUin,length(param.distUin),1)),[1,Nx,Ny]),[2,3,1]);
                if ndims(uin)==4 %mode uem with several planes
                    uin_out = zeros([Nx,Ny,Nz,size(uin,ndims(uin)),Ntheta]);
                else
                    uin_out = zeros([Nx,Ny,Nz,Ntheta]);
                end
                if false
                    fprintf('THE PROPAGATOR IS ALSO LOW PASS\n')
                    mask = dphi < k;
                else
                    mask = 1;
                end
                propagator = arrayfun(@(x) mask.*exp(-1i*dphi*x),param.distUin,'UniformOutput',false);
                propagator = cat(3,propagator{:});
                for kk = 1:Ntheta
                    curr_fuin = repmat(fft2((uin(repel{:},kk))),[1,1,Nz]);
                    curr_out = ifft2(curr_fuin.*propagator);
                    uin_out(repel{:},:,kk) = curr_out.*expkz;
                end
                uin_out = squeeze(uin_out);
            else
                uin_out = uin;
            end
        else
            uin_out = [];
        end
        
end

end