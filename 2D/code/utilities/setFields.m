function [uin_out,y] = setFields(uin,y,param,simP,model,varargin)
%SETFIELDS Propagtes uin /y at a given distance

mode = 'both';%determines which field(s) is/are propagated, 'both','uin','uM'
if nargin > 5
    mode = varargin{1};
end
switch mode
    case 'uin'
        param.distuM = 0;
    case 'uM'
        param.distUin = 0;
end

if isempty(uin)
    Nx = size(y,1);
else
    Nx = size(uin,1);
end
padfact = 1;
dkx = 2*pi/(padfact*Nx*simP.dx);
Nx = padfact*Nx;
if mod(Nx,2)==0
    Kxx = (dkx * [0:(Nx/2-1), (-Nx/2):-1]').^2;
else
    Kxx = (dkx * [0:(Nx-1)/2, (-(Nx-1)/2):-1]').^2;
end

k = 2*pi/(simP.lambda0)*simP.n0;

dphi = real(Kxx./(k + real(sqrt(k^2 - Kxx)))); % diffraction phase factor

switch model
    case {'bpm','bpmi','rytovFBP'} %propagates uin -0.5 Lz and uM 0.5 Lz
        
        if strcmp(mode,'uin') && ~ismatrix(uin)
            warning('uin not in the expected format\nPropagating a slice in dim 1 for each angle')
            uin = squeeze(uin(:,end/2,:));
        end
        
        if length(param.distUin) > 1
            error('Number of propagations of Uin too high, should be 1 for BPM and cie\n');
        end
        
        if param.distuM ~=0
            y = ifft(fft(padarray(y,(padfact - 1)*Nx/2,0,'both')).*repmat(exp(-1i*dphi*param.distuM),1,size(y,2)));
            %y = exp(1i*k*param.distuM)*y;
        end
        if any(param.distUin ~=0)
            uin_out = ifft(fft(padarray(uin,(padfact - 1)*Nx/2,0,'both')).*repmat(exp(-1i*dphi*min(param.distUin)),1,size(uin,2)));
            %uin_out = exp(1i*k*param.distUin)*uin_out;
        else
            uin_out = uin;
        end
    otherwise
        if param.distuM ~=0
            y = exp(1i*k*param.distuM)*ifft(fft(y).*repmat(exp(-1i*dphi*param.distuM),1,size(y,2)));
        end
        if any(param.distUin~=0) && ~isempty(uin)
            uin_out = setUinRealData(uin,simP,param.distUin);
        else
            uin_out = uin;
        end
end
end