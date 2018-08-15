%% MAIN 2D ManuLipp

fpath = createFolder(par.folder,'lipp',par.param_fname);%copy main_param in (newly created) folder

fprintf('%1.2f < 0.5 ? \n',par.dx/simP.lambda0*simP.n0);
%% Set Measurements
%Assume same ratio for instance in all dim (but can be converted for anisotropic case)
%It will require to know in which dimension lie the measurements

par.curr_siz = detSize(simP,par,'lipp');

[uin_lipp, noise_lipp] = setUin(par,simP,'lipp',uin,SP);

curr_SP = SP;
if par.discreteGtild
    uem_lipp = squeeze(zeros([M_curr.sizeout,par.Ntheta]));
    for kk = 1:par.Ntheta
        uem_lipp(:,kk) = M_curr*uin_lipp(:,:,kk);
    end
else
    if hasandis(simP,'realData')
        if ~hasandis(simP,'fresnel')
            curr_SP = SP;
            %such that it is at the end of a physical area of 2*par.distuM
            curr_SP.z = par.Lzext*ones(curr_SP.Nmeas,1);%2*par.distuM*ones(curr_SP.Nmeas,1);
        end
    else
        par.SP.Nmeas = length(par.SP.x);
    end
    if exist('uem','var')
        uem_lipp = uem;
    elseif ~isempty(uin)
        uem_lipp = uin;
    else
        error('Please provide Uem or Uin at the sensor positions');
    end
end

[uM.lipp,M_lipp] = setMeasurements(uM.curr_y,par,simP,par.noreflection,'lipp',...
    M_curr,curr_SP, uem);
%%

H = setmultOp(uin_lipp,M_curr,par,uem_lipp,curr_SP, simP);

%%
F = cell(par.Ntheta,1);
for kk = 1:par.Ntheta
    F{kk} = CostL2([],uM.lipp(:,kk))*H.mapsCell{kk};
end
F = CostSummation(F,1/par.Lsub);
F = F.makePartialSummation(par.Lsub);
F.partialGrad = par.stochGrad;
%%
Reg = CostMixNorm21([par.siz,2],3)*LinOpGrad(par.siz,[],'circular');

pos_ROIz = (par.Nzext - par.Nz)/2;
pos_ROIx = (par.Nxext - par.Nx)/2;

if exist('f','var')
    n_gt.lipp = setInput(n,simP,par,par.curr_siz,simP.n0);
    f_res.lipp = setInput(f,simP,par,par.siz_ext,0)*par.ratio_x^2;
    f_res.lipp = f_res.lipp(1 + pos_ROIx:end-pos_ROIx, 1 + pos_ROIz:end-pos_ROIz);
    
    dispG(f_res.lipp,n_gt.lipp,uin_lipp,uM.lipp,H.mapsCell{1}.G);
else
    f_res.lipp = [];
    n_gt.lipp = zeros(par.curr_siz);
end

Reg.setProxAlgo(par.bounds,par.prox_iter,par.xtolin,...
    myOutputOpti(1,f_res.lipp,5,par.kdz,simP.n0,'Lipp',par.negative));
Reg.optim.ItUpOut = 0;

%%

tic
u_gt = H*f_res.lipp;
toc

compareFields(u_gt, uM.lipp);
%% Run FISTA for Lippmann
figure;
for kk = 1:length(regul_set)
    par.regul = regul_set(kk);
    par.gam = gam_set(1);
    fprintf('Starting regul %1.2e, gam %1.2e, Niterin %i\n',par.regul,par.gam,par.Niterin);
    OutOp = OutputOptiLipp(1,f_res.lipp,par.iterVerb,[par.kdz,simP.k],simP.n0,'lipp',par.negative);
    OutOp.doPlot = par.doPlot;
    
    if isa(Reg,'CostMultiplication')
        Reg = par.regul*Reg.cost2;
    else
        Reg = par.regul*Reg;
    end
    FBS = OptiFBS(F,Reg);
    
    FBS.OutOp = OutOp;
    FBS.reducedStep = par.redStep;
    FBS.maxiter = par.Niter;               % set maxiter (default 50)
    FBS.CvOp = TestCvgStepRelative(par.xtol);
    FBS.fista = true;                % activate FISTA
    FBS.gam = par.gam;
    FBS.mingam = FBS.gam/12;
    FBS.ItUpOut = par.ItUpOut;
    
    %input = simP.n0*ones(par.siz);%choose this initialisation for Fresnel
    input = setInput(n_hat.fbp,simP,par,par.siz,simP.n0);%Shepp-Logan
    
    input = par.kdz^2*((input/simP.n0).^2 - 1);
    
    FBS.run(input);
    
    n_hat.lipp = padarray(sqrt(FBS.xopt/par.kdz^2 + 1)*simP.n0,...
        (par.curr_siz - par.siz)/2, simP.n0);
    
    fname = sprintf('%s',simP.ObjectType);
    if strcmp(simP.ObjectType,'mult_bead') && isfield(simP,'nBeads')
        fname = sprintf('%s_%i_low',fname,simP.nBeads);
    end
    if isfield(simP,'dn')
        fname = sprintf('%s_dn_%1.2f',fname,simP.dn(end));
    end
    fname = sprintf('%s_lipp_regul_%1.2e_Nx_%i_Nxext_%i_Nz_%i_Nzext_%i_iterin_%i_snr_%1.2f_gam_%1.2e.mat',...
        fname,regul_set(kk),par.Nx,par.Nxext,par.Nz,par.Nzext,...
        F.mapsCell{1}.H2.Niter,snr(n_gt.lipp,n_gt.lipp - n_hat.lipp),FBS.gam);
    fprintf('Saving data as %s...\nin %s\n',fname,fpath);
    OutOp = FBS.OutOp;
    save(fullfile(fpath,fname),'OutOp','input','par','n_hat','-v7.3');%,'FBS'
end