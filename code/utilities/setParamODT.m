function [par,M_curr,SP] = setParamODT(par,simP,SP)
par.siz = [par.Nx,par.Nz];par.siz_ext = [par.Nxext,par.Nzext];

par.NOvs = 2*par.siz;%Oversampling for Green fct (2 * support)

par = setScanningAngle(par,simP,par.excluded_angles);

par.ratio_x = par.dx/simP.dx;
par.ratio_z = par.dz/simP.dz;
%For display purpose
if hasandis(simP,'realData')
    par.Lz_vec = par.dz:par.dz:par.Lzext;
    par.Lx_vec = par.dx:par.dx:par.Lxext;
else
    par.Lz_vec = par.dz:par.dz:simP.Lz;
    par.Lx_vec = par.dx:par.dx:simP.Lx;
end
%free-space and background wavenumbers
par.k0 = 2*pi/simP.lambda0;
par.k = par.k0*simP.n0;
par.kdz = par.k*par.dz;

if isfield(par.SP,'M')%retrocompatibility with old simulated stuffs
    par.SP.Nmeas = par.SP.M;
    par.SP = rmfield(par.SP,'M');
elseif hasandis(simP,'fresnel')
    par.SP.Nmeas = length(par.SP.x)/simP.Ntheta;
else
    par.SP.Nmeas = length(par.SP.x);
end

if isfield(SP,'M') %retrocompatibility with old simulated stuffs
    SP.Nmeas = SP.M;
    SP = rmfield(SP,'M');
elseif hasandis(simP,'fresnel')
    SP.Nmeas = length(SP.x)/simP.Ntheta;
else
    SP.Nmeas = length(SP.x);
end

if simP.Lxroi < simP.Lx
    ind_rm = par.SP.x - simP.Lx/2 < -simP.Lxroi/2 | par.SP.x - simP.Lx/2 > simP.Lxroi/2;
    par.SP.x(ind_rm) = [];
    par.SP.z(ind_rm) = [];
    par.SP.x = par.SP.x - (simP.Lx - simP.Lxroi)/2;
    par.SP.Nmeas = length(par.SP.x);
    
    ind_rm = SP.x - simP.Lx/2 < - simP.Lxroi/2 | SP.x - simP.Lx/2 > simP.Lxroi/2;
    SP.x(ind_rm) = [];
    SP.z(ind_rm) = [];
    SP.x = SP.x - (simP.Lx - simP.Lxroi)/2;
    SP.Nmeas = length(SP.x);
end

if par.Mroi < 1
    par.SPorig = par.SP;
    ind_rm = abs(par.SP.x - simP.Lxroi/2) > par.Mroi*simP.Lxroi/2;
    par.SP.x(ind_rm) = [];
    par.SP.z(ind_rm) = [];
    par.SP.Nmeas = length(par.SP.x);
end

if par.Lxext > simP.Lx %mainly for BPM
    %Correct the x positions if window OI is larger
    %Assumes extra length/2 is added on both sides
    par.SP.x = par.SP.x + (par.Lxext - simP.Lx)/2;
elseif par.Lxext < simP.Lx %mainly for BPM
    par.SP.x = par.SP.x - (simP.Lx - par.Lxext)/2;
end

[M_curr, par.SP] = setM(par.Nxext,par.Nzext,...
    par.dx,par.dz,par.SP);

if isfield(par,'SPorig')%means par.Mroi < 1
    [~, par.SPorig] = setM(par.Nxext,par.Nzext,...
    par.dx,par.dz,par.SPorig);
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