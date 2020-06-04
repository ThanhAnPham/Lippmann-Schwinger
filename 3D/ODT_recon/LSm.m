%% MAIN 3D FastLipp

fpath = createFolder(par.folder,'lipp',par.param_fname);%copy main_param in (newly created) folder

fprintf('%1.2f < 0.5 ? \n', par.dx/simP.lambda0*simP.n0);
%% Set 3D incident fields

cuin = uin;

uin_lipp = setUin3D(par,simP,'lipp',single(cuin),SP);

%% Set Measurements
if par.discreteGtild
    if par.SP.nCamera > 1
        uem_lipp = curr_uem(:,:,:,par.curr_thetas);
    else
        uem_lipp = curr_uem(:,:,par.curr_thetas);
    end
else
    if exist('curr_uem','var')
        uem_lipp = curr_uem(:,:,par.curr_thetas);
    elseif ~isempty(cuin)
        uem_lipp = cuin(:,:,par.curr_thetas);
    else
        error('Please provide Uem or Uin at the sensor positions');
    end
end

[uM.lipp,M_lipp,uem_lipp] = setMeasurements3D(uM.curr_y,par,simP,par.noreflection,'lipp',...
    M_curr,curr_SP, uem_lipp);

%% Set the forward model and regularization

H = setmultOp3Dnew(uin_lipp,M_lipp,par,'lipp',uem_lipp,curr_SP, simP);

F = createCostsArray(H,uM.lipp,true);
F = F.makePartialSummation(par.Lsub);
F.partialGrad = par.stochGrad;

switch par.RegType
    case 'TV'
        if isfield(par,'mu') && par.mu <= 1
            Reg = CostMixNorm21([par.siz,3],4)*LinOpGrad(par.siz,[],'circular',1./[par.mu,par.mu,1 - par.mu]);
        else
            Reg = CostMixNorm21([par.siz,3],4)*LinOpGrad(par.siz,[],'circular');
        end
    case 'Hess'
        Hess = LinOpHess(par.siz,'circular');
        Reg = CostMixNormSchatt1(Hess.sizeout,par.p)*Hess;
        Reg.setC(par.bounds(1),par.bounds(2));
        Reg.max_iter = par.prox_iter;
        Reg.xtol = par.xtolin;
end

%%
if exist('n','var') && ~isempty(n)
    n_gt.lipp = setInput3D(n,simP,par,par.siz,simP.n0);
    f_res.lipp = (par.k*par.dz)^2*(n_gt.lipp.^2/simP.n0^2 - 1);
    
    figure(10);dispGT3D(par,simP,[],f_res.lipp,par.negative,[],gcf);
else
    f_res.lipp = [];
    n_gt.lipp = zeros(par.siz);
end
if strcmp(par.RegType,'TV')
    Reg.setProxAlgo(par.bounds,par.prox_iter,par.xtolin,...
        myOutputOpti(1,n_gt.lipp,5,par.k,simP.n0,'Lipp',par.negative));
    Reg.optim.ItUpOut = 0;
end

%% Test forward
aroi = 3;
tic
u_gt = gather(H.mapsCell{aroi}*(gpuCpuConverter(f_res.lipp)));
compareFields(u_gt, uM.lipp(:,:,aroi));
drawnow;
toc

%% Run FISTA for Lippmann

for kk = 1:length(regul_set)
    par.regul = regul_set(kk);
    
    fprintf('Starting regul %1.2e, gam %1.2e, Niterin %i\n',par.regul,par.gam,par.Niterin);
    
    OutOp = myOutputOpti(1,f_res.lipp,par.iterVerb,[par.k*par.dz,simP.k],simP.n0,'lipp',par.siz/2);%[135,112,50]);%
    OutOp.doPlot = par.doPlot;OutOp.saveX = par.saveX;
    if isfield(par,'saveAll')
        OutOp.saveAll = par.saveAll;
    end
    
    if isa(Reg,'CostRectangle')
        %F.alpha(end) = par.regul;
    elseif isa(Reg,'CostMultiplication')
        Reg = par.regul*Reg.cost2;
    else
        Reg = par.regul*Reg;
    end
    
    FBS = OptiFBS(F, Reg);
    
    FBS.OutOp = OutOp;
    FBS.fista = true;
    FBS.updateGam = par.redStep;
    FBS.gam = par.gam;
    FBS.ItUpOut = par.ItUpOut;
    FBS.CvOp = TestCvgStepRelative(par.xtol);
    FBS.maxiter = par.Niter;
    FBS.momRestart = par.momRestart;
    FBS.alpha = par.alpha;
    
    input0 = setInput3D(input0,simP,par,par.siz,simP.n0);
    
    FBS.run(gpuCpuConverter(input0));
    
    n_hat.lipp = gather(padarray(sqrt(FBS.xopt/(par.k*par.dz)^2 + 1)*simP.n0,...
        (par.siz - par.siz)/2, simP.n0));
    
    fname = sprintf('%s',simP.ObjectType);
    if strcmp(simP.ObjectType,'mult_bead') && isfield(simP,'nBeads')
        fname = sprintf('%s_%i_low',fname,simP.nBeads);
    end
    if isfield(simP,'dn')
        fname = sprintf('%s_dn_%1.2f',fname,simP.dn(end));
    end
    if isempty(n_gt.lipp)
        curr_snr = nan;
    else
        curr_snr = snr(n_gt.lipp,n_gt.lipp - n_hat.lipp);
    end
    fname = sprintf('%s_lipp_regul_%1.2e_Nx_%i_Nxext_%i_Nz_%i_Nzext_%i_iterin_%i_snr_%1.2f_gam_%1.2e_RegType_%s.mat',...
        fname,regul_set(kk),par.Nx,par.Nxext,par.Nz,par.Nzext,...
        par.Niterin,curr_snr,par.gam,par.RegType);
    fprintf('Saving data as %s...\n',fname);
    OutOp = FBS.OutOp;
    
    curr_nhat = n_hat;
    curr_uM = gather(uM.lipp);
    curr_snr = gather(curr_snr);
    n_hat = struct('lipp',curr_nhat.lipp);
    save(fullfile(fpath,fname),'OutOp','input0','par','n_hat','curr_uM','-v7.3');
    n_hat = curr_nhat;
    clear curr_nhat
end