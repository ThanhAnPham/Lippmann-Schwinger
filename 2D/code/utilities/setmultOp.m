function out = setmultOp(uin,M,par,uem,varargin)
% SETMULTOP : set the stack of (Linear) Operators

Ntheta = par.Ntheta;
out = cell(Ntheta,1);


G = LinOpFreeSpaceKernel(par.kdz,par.siz,2*par.siz);

if par.discreteGtild
    Gtil = LinOpGtild_FFT(par.kdz, par.siz, par.siz_ext,...
        2*par.siz_ext,M.sel);
else
    if nargin > 4 && ~isempty(varargin{1})
        if nargin > 6 && ~isempty(varargin{2})
            ddx = varargin{2}.dx/2;
            ddz = varargin{2}.dz/2;
        else
            ddx = 0;
            ddz = 0;
        end
        curr_SP.x = varargin{1}.x - ddx;
        curr_SP.z = varargin{1}.z - ddz;%dx,dz should be identical
        curr_SP.Nmeas = varargin{1}.Nmeas;
    else %imprecise because of rounderror
        curr_SP.x = par.SP.x*par.dx;% - par.dx/2;
        curr_SP.z = par.SP.z*par.dz;%dx,dz should be identical
        curr_SP.Nmeas = par.SP.Nmeas;
    end
    if curr_SP.Nmeas~=length(curr_SP.x)
        fprintf('Different sensor positions for each view\n');
        %Ntheta = length(curr_SP.x)/curr_SP.Nmeas;
        Gtil = cell(Ntheta,1);
        for kk = 1:Ntheta
            curr_ind = find(par.curr_thetas,kk);curr_ind = curr_ind(end);
            subSP.x = curr_SP.x(1 + (curr_ind-1)*curr_SP.Nmeas:curr_ind*curr_SP.Nmeas);
            subSP.z = curr_SP.z(1 + (curr_ind-1)*curr_SP.Nmeas:curr_ind*curr_SP.Nmeas);
            Gtil{kk} = LinOpGtild(par.kdz,par.siz,subSP,...
                par.dx,par.storeGtild,par.dx/2);%for convention of Green function  H_0^(1,2)
        end
    else %must be planar
        curr_SP = centerPos(curr_SP,max(curr_SP.x),max(curr_SP.z));
        if nargin > 4 && hasandis(varargin{2},'realData')
            curr_SP.z = par.distuM*ones(curr_SP.Nmeas,1);
        end
        if par.noreflection
            curr_SP.x = curr_SP.x(curr_SP.z > 0);
            curr_SP.z = curr_SP.z(curr_SP.z > 0);
        end
        Gtil = LinOpGtild(par.kdz,par.siz,curr_SP,...
            par.dx,par.storeGtild,par.dx/2);
    end
end

for kk = 1:Ntheta
    if iscell(Gtil)
        out{kk} = OpLipp(G,Gtil{kk},uin(:,:,kk),par.Niterin,par.xtolin,~par.yscat,uem(:,kk));
    else
        out{kk} = OpLipp(G,Gtil,uin(:,:,kk),par.Niterin,par.xtolin,~par.yscat,uem(:,kk));
    end
end
out = StackMap(out);%useful for computing the whole forward
end