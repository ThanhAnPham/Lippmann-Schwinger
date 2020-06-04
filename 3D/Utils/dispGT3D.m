function dispGT3D(par,simP,n,nhat,negative,varargin)
%DISPGT3D 

curr_siz = size(nhat);

if nargin > 5
    colval = varargin{1};
end

if ~exist('colval','var') || isempty(colval)
    if isempty(n)
        colval = [min(nhat(:)),max(nhat(:))];
    else
        colval = sort([simP.n0,simP.n0 + max(simP.dn)]);
    end
end
if nargin > 6
    figure(varargin{2});
else
    figure;
end
if nargin > 7
    pos = varargin{3};
else
    pos = round(size(nhat)/2);
end
% XY
subplot(221);
if isfield(par,'Lx_vec')
    imagesc(par.Ly_vec,par.Lx_vec,nhat(:,:,pos(3)));
    if ~isempty(n)
        dispGTover3D(simP,n,1,2,0.5,curr_siz([2,1]).*[par.dy,par.dx]*0.5);
    end
    daspect([par.dx,par.dy,par.dz]/par.dx);
else
    imagesc(nhat(:,:,pos(3)));
end
title('XY');
axis image;
colorbar;

caxis(colval);
% XZ
subplot(222);
if isfield(par,'Lx_vec')
    imagesc(par.Lz_vec,par.Lx_vec,squeeze(nhat(:,pos(2),:)));
    if ~isempty(n)
        dispGTover3D(simP,n,1,3,0.5,curr_siz([3,1]).*[par.dz,par.dx]*0.5);
    end
    daspect([par.dz,par.dx,par.dy]/par.dx);
else
    imagesc(squeeze(nhat(:,pos(2),:)));
end
title('XZ');
axis image;
colorbar;
caxis(colval);

% ZY
subplot(223);
if isfield(par,'Ly_vec')
    imagesc(par.Ly_vec,par.Lz_vec,squeeze(nhat(pos(1),:,:))');
    if ~isempty(n)
        dispGTover3D(simP,n,3,2,0.5,curr_siz([2,3]).*[par.dy,par.dz]*0.5);
    end
    daspect([par.dy,par.dz,par.dx]/par.dx);
else
    imagesc(squeeze(nhat(pos(1),:,:)));
end
title('ZY');
axis image;
colorbar;
caxis(colval);
if negative
    colormap(flipud(viridis));
end
drawnow
end

