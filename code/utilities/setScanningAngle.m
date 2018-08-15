function par = setScanningAngle(par,simP,varargin)
%par
%   fields : max_scan_theta
%simP
%   fields : Ltheta, theta
orig_thetas = simP.thetas;

if nargin > 2
    simP.thetas(varargin{1}) = [];%exclude
    simP.Ntheta = length(simP.thetas);
end

if isscalar(par.max_scan_theta)
    doReduced = par.max_scan_theta < simP.Ltheta;
elseif length(par.max_scan_theta)==2
    doReduced = diff(par.max_scan_theta) < simP.Ltheta;
else
    doReduced = false;
end
if doReduced
    if isscalar(par.max_scan_theta)
        par.curr_thetas = find(abs(simP.thetas) <= par.max_scan_theta/2);
    else
        par.curr_thetas = find(simP.thetas >= par.max_scan_theta(1)...
            & simP.thetas <= par.max_scan_theta(2));
    end
%     if isfield(par,'Ntheta') && length(par.curr_thetas) > par.Ntheta
%         par.curr_thetas = par.curr_thetas(1) - 1 ...
%             + unique(round(1:length(par.curr_thetas)/par.Ntheta:length(par.curr_thetas)));
%     end
else
    par.curr_thetas = true(simP.Ntheta,1);
end

par.thetas = simP.thetas(par.curr_thetas);
par.Ntheta = length(par.thetas);
par.curr_thetas = ismember(orig_thetas,par.thetas);
par.curr_thetas(varargin{1}) = 0;

end