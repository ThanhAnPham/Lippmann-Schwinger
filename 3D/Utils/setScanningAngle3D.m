function par = setScanningAngle3D(par,simP,varargin)
%% SETSCANNINGANGLE3D: Select a set of angles
%par
%   fields : max_scan_theta
%simP
%   fields : Ltheta, theta
orig_thetas = simP.thetas;

if iscell(simP.Ltheta)
    Ltheta = cell2mat(simP.Ltheta);
else
    Ltheta = simP.Ltheta;
end
if nargin > 2
    simP.thetas(varargin{1},:) = [];%exclude
    simP.Ntheta = size(simP.thetas,1);
end
for kk = 1:length(par.max_scan_theta)
    if isscalar(par.max_scan_theta{kk})
        doReduced = par.max_scan_theta{kk} < Ltheta(kk);
    elseif length(par.max_scan_theta{kk})==2
        doReduced = diff(par.max_scan_theta{kk}) < Ltheta(kk);
    else
        doReduced = false;
    end
    if doReduced
        if isscalar(par.max_scan_theta{kk})
            par.curr_thetas(:,kk) = abs(simP.thetas(:,kk)) <= par.max_scan_theta{kk}/2;
        else
            par.curr_thetas(:,kk) = simP.thetas(:,kk) >= par.max_scan_theta{kk}(1)...
                & simP.thetas(:,kk) <= par.max_scan_theta{kk}(2);
        end
    elseif kk > 1
        par.curr_thetas(:,kk) = true(size(par.curr_thetas,1),1);
    else
        par.curr_thetas(:,kk) = true(simP.Ntheta,1);
    end
end


if isfield(par,'Ntheta') && par.Ntheta < nnz(all(par.curr_thetas,2))
    curr_thet = find(all(par.curr_thetas,2));
    par.curr_thetas(curr_thet(round(linspace(1,length(curr_thet),length(curr_thet) - par.Ntheta))),:) = 0;
end

par.thetas = simP.thetas(all(par.curr_thetas,2),:);
par.Ntheta = size(par.thetas,1);
par.curr_thetas = ismember(orig_thetas,par.thetas,'rows');
if nnz(par.curr_thetas)> par.Ntheta %there is / are twice the same illumination angles
    curr_thet = orig_thetas;
    [u_thet,ia] = unique(curr_thet,'rows') ;
    curr_thet = zeros(size(curr_thet)) ;
    curr_thet(ia,:) = u_thet;
    par.curr_thetas = ismember(curr_thet,par.thetas,'rows');
    if nnz(par.curr_thetas) < par.Ntheta
        fprintf('Duplicate angles...Taking some to match the asked Ntheta\n');
        sel = find(~par.curr_thetas);
        par.curr_thetas(sel(round(linspace(1,length(sel),par.Ntheta - nnz(par.curr_thetas))))) = true;
    end
end
par.curr_thetas(varargin{1}) = 0;
par.curr_ind_thetas = find(par.curr_thetas);
end