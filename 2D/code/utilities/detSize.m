function siz = detSize(simP,par,varargin)
%DETSIZE Determine the correct size

%if simP.Lx < par.Lxext
if hasandis(simP,'realData') %&& (nargin > 2 && (contains(varargin{1},'lipp') || contains(varargin{1},'rytj')))
    siz = par.siz_ext;
else
    siz = round([simP.Lx/par.dx,simP.Lz/par.dz]);
end
%else
%    siz = round([simP.Lx/par.dx,simP.Lz/par.dz]);
%end
end

