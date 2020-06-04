function [M,SP_out] = setM3D(Nx,Ny,Nz,dx,dy,dz,SP,simP)
%% SETM Set the mask for measurement
x = SP.x;
y = SP.y;
z = SP.z;

if any(x < 0)
    x = x + simP.Lx/2;
end
if any(y < 0)
    y = y + simP.Ly/2;
end

M = false(Nx,Ny,Nz);

SP_out.x = round((1e14*x)/(1e14*dx),7);
SP_out.y = round((1e14*y)/(1e14*dy),7);
SP_out.z = round((1e14*z)/(1e14*dz),7);



for kk = 1:length(x)
    %ceil, floor or round ?
    if SP_out.x(kk) <= Nx  && SP_out.x(kk) > 0 ...
            && SP_out.y(kk) <= Ny && SP_out.y(kk) > 0 ...
            && SP_out.z(kk) <= Nz &&  SP_out.z(kk) > 0
        M(min(max(ceil(SP_out.x(kk)),1),Nx),...
            min(max(ceil(SP_out.y(kk)),1),Ny),...
            min(max(ceil(SP_out.z(kk)),1),Nz)) = true;
    end
end

ind_rm = SP_out.y <= 0 | SP_out.y > Ny;
SP_out.x(ind_rm) = [];
SP_out.y(ind_rm) = [];
SP_out.z(ind_rm) = [];

ind_rm = SP_out.x <= 0 | SP_out.x > Nx;
SP_out.x(ind_rm) = [];
SP_out.y(ind_rm) = [];
SP_out.z(ind_rm) = [];

SP_out.unique_pos = unique(ceil([SP_out.x(:),SP_out.y(:), SP_out.z(:)]),'rows','stable');%might contain even data out of ROI

M = LinOpSelector(M);
sizeout = [length(unique(ceil(SP_out.x))),length(unique(ceil(SP_out.y)))];
if length(unique(z)) > 1
    sizeout = [sizeout,length(unique(z))];
end
if prod(M.sizeout) > 0
    M = LinOpShape(M.sizeout,sizeout)*M;
end

if prod(SP.Nmeas)==length(SP.x)
    SP_out.Nmeas = [length(unique(SP_out.unique_pos(:,1))),length(unique(SP_out.unique_pos(:,2)))];
else %special case (Guess: Fresnel)
    SP_out.Nmeas = SP.Nmeas;
end
SP_out.nCamera = length(unique(SP_out.z));
end

