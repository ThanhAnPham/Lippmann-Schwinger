function out = setInput3D(in,par_in,par_out,siz_out,n0,varargin)
%% SETINPUT3D
%   setInput3D(in,par_in,par_out,siz_out,n0,varargin)


%if isfield(par1,'Lxext')
%    if nargin > 5 && varargin{1}
%        ratiox = par2.dx*siz_out(1)/par1.Lxext*par1.Lxext/par1.Lx;
%        ratioy = par2.dy*siz_out(2)/par1.Lyext*par1.Lyext/par1.Ly;
%        ratioz = par2.dz*siz_out(3)/par1.Lzext*par1.Lzext/par1.Lz;
%    else
%        ratiox = par2.dx*siz_out(1)/par1.Lxext;
%        ratioy = par2.dy*siz_out(2)/par1.Lyext;
%        ratioz = par2.dz*siz_out(3)/par1.Lzext;
%    end
%elseif isfield(par1,'siz')
%ratiox = par2.dx*siz_out(1)/(par1.siz(1)*par1.dx);
%ratioy = par2.dy*siz_out(2)/(par1.siz(2)*par1.dy);
%ratioz = par2.dz*siz_out(3)/(par1.siz(3)*par1.dz);
%end


%dx1 = par1.dx; dz1 = par1.dz; dx2 = par2.dx; dz2 = par2.dz;
[Nx1,Ny1,Nz1] = size(in);
%Nx2 = par2.Nx; Nz2 = par2.Nz;
ratiox = siz_out(1)*par_out.dx/(par_in.dx*Nx1);
ratioy = siz_out(2)*par_out.dy/(par_in.dy*Ny1);
ratioz = siz_out(3)*par_out.dz/(par_in.dz*Nz1);

new_Nx = round(ratiox*Nx1); new_Nx = new_Nx + mod(new_Nx,2);
new_Ny = round(ratioy*Ny1); new_Ny = new_Ny + mod(new_Ny,2);
new_Nz = round(ratioz*Nz1); new_Nz = new_Nz + mod(new_Nz,2);
out = in;
%Select the "ROI" (or expand)
if ratiox < 1
    out = out(1 + end/2 - new_Nx/2:end/2 + new_Nx/2,:,:);
else
    out = padarray(out,[(new_Nx - Nx1)/2,0,0],n0,'both');
end

if ratioy < 1
    out = out(:,1 + end/2 - new_Ny/2:end/2 + new_Ny/2,:);
else
    out = padarray(out,[0,(new_Ny - Ny1)/2,0],n0,'both');
end

if ratioz < 1
    out = out(:,:,1 + end/2 - new_Nz/2:end/2 + new_Nz/2);
else
    out = padarray(out,[0,0,(new_Nz - Nz1)/2],n0,'both');
end
if exist('imresize3','file')
    out = imresize3(out,siz_out,'linear');
else
    out = myimresize3(out,siz_out,'linear');
end