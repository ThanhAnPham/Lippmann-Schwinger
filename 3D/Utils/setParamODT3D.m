function [par,M_curr,SP] = setParamODT3D(par,simP,SP)
%% SETPARAMODT3D Set the parameters for 3D ODT reconstruction

par.siz = [par.Nx,par.Ny,par.Nz];par.siz_ext = [par.Nxext,par.Nyext,par.Nzext];

par.ForwardModel = 'Vainikko';%CroppedGreen,Vainikko
par.NOvs = 2*par.siz;%Oversampling for Green fct (2 * support)

par = setScanningAngle3D(par,simP,par.excluded_angles);

par.ratio_x = par.dx/simP.dx;
par.ratio_x = par.dy/simP.dy;
par.ratio_z = par.dz/simP.dz;
%For display purpose
if hasandis(simP,'realData') && ~hasandis(simP,'fresnel')
    
    %Where the data measurements are placed (at the end of the physical window
    %in the optical axis by default
    if isfield(par,'Zsen') && ~isempty(par.Zsen)
        par.SP.z = par.Zsen*ones(prod(SP.Nmeas),1);
    else
        par.SP.z = (par.Lzext/2 + par.distuM)*ones(prod(SP.Nmeas),1);
    end
    %par.SP.z = (2*par.distuM - par.shift)*ones(prod(SP.Nmeas),1);
    par.Lz_vec = par.dz:par.dz:par.Lz;
    par.Lx_vec = par.dx:par.dx:par.Lx;
    par.Ly_vec = par.dy:par.dy:par.Ly;
else
    par.Lz_vec = par.dz:par.dz:par.Lz;
    par.Lx_vec = par.dx:par.dx:par.Lx;
    par.Ly_vec = par.dy:par.dy:par.Ly;
end
%free-space and background wavenumbers
par.k0 = 2*pi/simP.lambda0;
par.k = par.k0*simP.n0;
%par.kdz = par.k*par.dz;

if simP.Lxroi < simP.Lx
    if hasandis(par.SP,'centered')
        ind_rm = par.SP.x <= -simP.Lxroi/2 | par.SP.x > simP.Lxroi/2;
    else
        ind_rm = par.SP.x - simP.Lx/2 <= -simP.Lxroi/2 | par.SP.x - simP.Lx/2 > simP.Lxroi/2;
    end
    par.SP.x(ind_rm) = [];
    par.SP.y(ind_rm) = [];
    par.SP.z(ind_rm) = [];
    
    par.SP.x = par.SP.x - (simP.Lx - simP.Lxroi)/2;
    par.SP.Nmeas = [length(unique(par.SP.x)), length(unique(par.SP.y))];
    
    if hasandis(SP,'centered')
        ind_rm = SP.x <= - simP.Lxroi/2 | SP.x > simP.Lxroi/2;
    else
        ind_rm = SP.x - simP.Lx/2 <= - simP.Lxroi/2 | SP.x - simP.Lx/2 > simP.Lxroi/2;
    end
    SP.x(ind_rm) = [];
    SP.y(ind_rm) = [];
    SP.z(ind_rm) = [];
    SP.x = SP.x - (simP.Lx - simP.Lxroi)/2;
    SP.Nmeas = [length(unique(SP.x)), length(unique(SP.y))];
end

if simP.Lyroi < simP.Ly
    if hasandis(par.SP,'centered')
        ind_rm = par.SP.y - simP.Ly/2 <= -simP.Lyroi/2 | par.SP.y - simP.Ly/2 > simP.Lyroi/2;
    else
        ind_rm = par.SP.y - simP.Ly/2 <= -simP.Lyroi/2 | par.SP.y - simP.Ly/2 > simP.Lyroi/2;
    end
    par.SP.x(ind_rm) = [];
    par.SP.y(ind_rm) = [];
    par.SP.z(ind_rm) = [];
    par.SP.y = par.SP.y - (simP.Ly - simP.Lyroi)/2;
    par.SP.Nmeas = [length(unique(par.SP.x)), length(unique(par.SP.y))];
    
    if hasandis(SP,'centered')
        ind_rm = SP.y <= - simP.Lyroi/2 | SP.y  > simP.Lyroi/2;
    else
        ind_rm = SP.y - simP.Ly/2 <= - simP.Lyroi/2 | SP.y - simP.Ly/2 > simP.Lyroi/2;
    end
    SP.x(ind_rm) = [];
    SP.y(ind_rm) = [];
    SP.z(ind_rm) = [];
    SP.y = SP.y - (simP.Ly - simP.Lyroi)/2;
    SP.Nmeas = [length(unique(SP.x)), length(unique(SP.y))];
end

if par.Mroi(1) < 1
    par.SPorig = par.SP;
    ind_rm = abs(par.SP.x - simP.Lxroi/2) > par.Mroi(1)*simP.Lxroi/2;
    par.SP.x(ind_rm) = [];
    par.SP.z(ind_rm) = [];
    par.SP.Nmeas = [length(unique(par.SP.x)), length(unique(par.SP.y))];
end
if par.Mroi(2) < 1
    ind_rm = abs(par.SP.y - simP.Lyroi/2) > par.Mroi(2)*simP.Lyroi/2;
    par.SP.x(ind_rm) = [];
    par.SP.y(ind_rm) = [];
    par.SP.z(ind_rm) = [];
    
    par.SP.Nmeas = [length(unique(par.SP.x)), length(unique(par.SP.y))];
end

if par.Lxext > simP.Lxroi %mainly for BPM
    %Correct the x positions if window OI is larger
    %Assumes extra length/2 is added on both sides
    par.SP.x = par.SP.x + (par.Lxext - simP.Lxroi)/2;
elseif par.Lxext < simP.Lx %mainly for BPM
    par.SP.x = par.SP.x - (simP.Lxroi - par.Lxext)/2;
end

if par.Lyext > simP.Lyroi %mainly for BPM
    %Correct the y positions if window OI is larger
    %Assumes extra length/2 is added on both sides
    par.SP.y = par.SP.y + (par.Lyext - simP.Lyroi)/2;
elseif par.Lyext < simP.Ly %mainly for BPM
    par.SP.y = par.SP.y - (simP.Lyroi - par.Lyext)/2;
end

if hasandis(simP,'fresnel')
    M_curr = LinOpIdentity(par.siz);
else
    [M_curr, par.SP] = setM3D(par.Nxext,par.Nyext,par.Nzext,...
        par.dx,par.dy,par.dz,par.SP,simP);
end
if isfield(par,'SPorig')%means par.Mroi < 1
    [~, par.SPorig] = setM3D(par.Nxext,par.Nyext,par.Nzext,...
        par.dx,par.dy,par.dz,par.SPorig,simP);
else %means par.Mroi == 1
    par.SPorig = par.SP;
end

if max(par.bounds) <= 0
    par.negative = true;%For colormap purpose, invert if sample is negative
else
    par.negative = false;
end

par.n0 = simP.n0;
if ~isfield(simP,'realData')
    simP.realData = false;
end
end