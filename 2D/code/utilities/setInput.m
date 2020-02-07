function out = setInput(in,par1,par2,siz_out,n0,varargin)

% if isfield(par1,'Lxext')
%     if nargin > 5 && varargin{1}
%         ratiox = par2.dx*siz_out(1)/par1.Lxext*par1.Lxext/par1.Lx;
%         ratioz = par2.dz*siz_out(2)/par1.Lzext*par1.Lzext/par1.Lz;
%     else
%         ratiox = par2.dx*siz_out(1)/par1.Lxext;
%         ratioz = par2.dz*siz_out(2)/par1.Lzext;
%     end
% else
[Nx1,Nz1] = size(in);
ratiox = par2.dx*siz_out(1)/(par1.dx*Nx1);
ratioz = par2.dz*siz_out(2)/(par1.dz*Nz1);
%end
%dx1 = par1.dx; dz1 = par1.dz; dx2 = par2.dx; dz2 = par2.dz;

%Nx2 = par2.Nx; Nz2 = par2.Nz;

new_Nx = round(ratiox*Nx1);% new_Nx = new_Nx + mod(new_Nx,2);
new_Nz = round(ratioz*Nz1);% new_Nz = new_Nz + mod(new_Nz,2);
out = in;
%Select the "ROI" (or expand)
if ratiox < 1
    out = out(1 + end/2 - new_Nx/2:end/2 + new_Nx/2,:);
else
    out = padarray(out,[(new_Nx - Nx1)/2,0],n0,'both');
end

if ratioz < 1
    out = out(:,1 + end/2 - new_Nz/2:end/2 + new_Nz/2);
else
    out = padarray(out,[0,(new_Nz - Nz1)/2],n0,'both');
end

%out = imresize(out,siz_out,'box');
out = imresize(out,siz_out);