function [M,SP_out] = setM(Nx,Nz,dx,dz,SP)
%SETM SET MASK FOR MEASUREMENT SAMPLING
%Could be done with imresize actually...
x = SP.x;
z = SP.z;

M = false(Nx,Nz);

SP_out.x = round((1e14*x)/(1e14*dx),7);
SP_out.z = round((1e14*z)/(1e14*dz),7);

for kk = 1:length(x)
    %ceil, floor or round ?
    if SP_out.x(kk) <= Nx  && SP_out.x(kk) > 0 ...
            && SP_out.z(kk) <= Nz &&  SP_out.z(kk) > 0
        M(min(max(ceil(SP_out.x(kk)),1),Nx), min(max(ceil(SP_out.z(kk)),1),Nz)) = true;
    end
end
M = LinOpSelector(M);

SP_out.unique_pos = unique(ceil([SP_out.x(:), SP_out.z(:)]),'rows','stable');%might contain even data out of ROI
if SP.Nmeas==length(SP.x)
    SP_out.Nmeas = prod(M.sizeout);
else %Special case
    SP_out.Nmeas = SP.Nmeas;
end
SP_out.nCamera = length(unique(SP_out.z));
end

