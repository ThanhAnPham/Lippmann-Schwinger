classdef myOutputOpti < handle
    %% myOutputOpti : class for object OutputOpti objects for tomography
    %  Matlab Inverse Problems Library
    %
    % At each ItUpOut iterations of an optimization algorithm (see Opti generic class),
    % the update method of an OutputOpti object will be called in order to acheive user
    % defined computations:
    %   - compute cost / SNR
    %   - store current iterate / cost value
    %   - plot/display stuffs
    %
    % The present generic class implements a basic update method that:
    %   - display the iteration number
    %   - computes & display the cost (if activated)
    %   - computes & display the SNR if ground truth is provided
    %
    % --Example
    %  VU=OutputOpti(computecost,xtrue,iterVerb)
    % where computecost is a boolean stating if the cost function has to be computed, xtrue is
    % the ground truth and iterVerb corresponds to the interval (in number of iterations) for
    % which the update will only perform the computation (every ItUpOut iterations) but not
    % print/display anithing. Hence iterVerb must be a multiple of ItUpOut Opti parameter.
    %
    % Child classes can derive from this one to define different actions to execute during
    % optimization updates.
    %
    % IMPORTANT: The update method should have an unique imput that is the OPTI object in order to
    % be generic for all Optimization routines. Hence the update method has acces (in reading mode)
    % to all the properties of OPTI objects.
    %
    % -- Properties
    % * |name|        - name of the VerbUpdate class
    % * |computecost| - Boolean, if true the cost function will be computed
    % * |xtrue|       - Ground Truth to compute the error with the solution (if provided)
    % * |evolcost|    - array to save the evolution of the cost function
    % * |evolsnr|     - array to save the evolution of the SNR
    % * |iterVerb|    - message will be displayed every iterVerb iterations (must be a multiple
    %                   of the ItUpOut Opti parameter)
    %
    % -- Methods
    % * |update|  - execute the defined actions
    % * |init|    - initialization (called at the starting of the opti algorithm to initialize
    %               internal variables (counter, array for saving...)
    %
    % See also Opti
    %
    %     Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    properties (SetAccess = protected,GetAccess = public)
        name = 'myOutputOpti'% name of the optimization algorithm
        computecost=false; % Boolean, if true the cost function will be computed
        xtrue;             % Ground Thruth (f)
        ntrue;             % Ground Truth refractive index (n)
        evolcost;          % array saving the evolution of the cost function
        evolsnr;           % array saving the evolution of the error with the groud truth
        evolxopt;          % cell saving the optimization variable xopt
        iternum;           % array saving the iteration number corresponding to evolcost, evolxopt and evolerr entries
        iterVerb=0;        % message will be displayed every iterVerb iterations
        nTimer;
    end
    properties (SetAccess = protected,GetAccess = public)
        normXtrue;         % Norm of the true signal (to compute snr)
        isgt=false;        % Boolean true if Ground Truth is provided
        count;             % internal counter
        normNtrue;         % Norm of the true refractive index (to compute snr)
        k_orig;
        k;
        n0;
        doSNRn = false;
        optiType = '';
        start;
    end
    properties
        doPlot = false;
        negative = false;
        saveX = false;
        pos;
        saveAll = inf;
        gam0 = 1e3;
        niter0 = 0;
    end
    
    methods
        %% Constructor
        function this=myOutputOpti(computecost,xtrue,iterVerb,k,n0,optiType,pos,varargin)
            if nargin==1
                this.computecost=computecost;
            elseif nargin==2
                this.computecost=computecost;
                this.xtrue=xtrue;
            elseif nargin==3
                this.computecost=computecost;
                this.xtrue=xtrue;
                this.iterVerb=iterVerb;
            elseif nargin==5
                this.computecost=computecost;
                this.xtrue=xtrue;
                this.iterVerb=iterVerb;
                this.doSNRn = true;
            elseif nargin==6
                this.computecost=computecost;
                this.xtrue=xtrue;
                this.iterVerb=iterVerb;
                this.doSNRn = true;
                this.optiType = optiType;
            elseif nargin >= 7
                this.computecost=computecost;
                this.xtrue=xtrue;
                this.iterVerb=iterVerb;
                this.doSNRn = true;
                this.optiType = optiType;
                this.pos = pos;
            end
            
            if ~isempty(this.xtrue)
                this.isgt=true;
                this.xtrue=xtrue;
                this.normXtrue=norm(this.xtrue(:));
                
                if this.doSNRn
                    this.k = k(1);
                    this.k_orig = k(end);
                    this.n0 = n0;
                    this.ntrue = this.n0*sqrt(this.xtrue/this.k^2 + 1);
                    this.normNtrue=norm(this.ntrue(:));
                else
                    fprintf('Missing k and n0 inputs...calculating SNR with f\n');
                end
            end
            
            if nargin > 7
                this.negative = varargin{1};
            end
        end
        %% Initialization
        function init(this)
            this.count=1;
            this.evolcost=[];
            this.evolsnr=[];
            this.nTimer = [];
            this.start = tic;
        end
        %% Update method
        function update(this,opti)
            str=sprintf('Iter: %5i',opti.niter);
            if this.computecost
                %if strcmpi(this.optiType(end),'i')
                switch opti.name
                    case 'OptiFBS_hierADMM'
                        this.evolcost(this.count,1) = opti.F.apply(opti.xopt);
                        this.evolcost(this.count,2) = opti.G.apply(opti.xopt);
                    case {'Opti FBS','Opti FBS SMLM','my Opti FBS'}
                        %if isa(opti.F,'CostPartialSummation')
                        %    this.evolcost(this.count,1) = 0;
                        %    for kk = 1:opti.F.Lsub
                        %        ind = opti.F.subset(kk);
                        %        this.evolcost(this.count,1) = this.evolcost(this.count,1)...
                        %            + opti.F.alpha(ind)*opti.F.mapsCell{ind}*(opti.xopt);
                        %    end
                        %    this.evolcost(this.count,1) = this.evolcost(this.count,1)/opti.F.Lsub;
                        %    str = sprintf('%s | %s',str,sprintf('%i ',opti.F.subset));
                        %else
                        this.evolcost(this.count,1) = gather(opti.F.apply(opti.xopt));
                        if isfield(opti.F,'mapsCell') && (isa(opti.F.mapsCell{1}.H2,'OpBPM3D') || isa(opti.F.mapsCell{1}.H2,'OpLipp3D'))
                            for kk = 1:opti.F.Lsub
                                opti.F.mapsCell{opti.F.subset(kk)}.H2.u = [];
                            end
                        end
                        %                        end
                        this.evolcost(this.count,2) = gather(opti.G.apply(opti.xopt));
                    case {'Opti ADMM','OptiADMMReweightedL1','Opti ADMM PU'}
                        for kk = 1:length(opti.Fn)
                            this.evolcost(this.count,kk) = opti.Fn{kk}.apply(opti.Hn{kk}*opti.xopt);
                        end
                        if ~isempty(opti.F0)
                            this.evolcost(this.count,length(opti.Fn) + 1) = opti.F0.apply(opti.xopt);
                        end
                        %this.evolcost(this.count,3) = sum(sum(opti.rho_n(1)*abs(opti.Hnx{end}-opti.yn{end}).^2));
                    case 'Opti Hierarchical ADMM'
                        if strcmp(opti.solver.type,'OptiFBS') %inefficient but lazy
                            this.evolcost(this.count,2) = opti.solver.G*opti.xopt;
                            this.evolcost(this.count,1) = opti.cost*opti.xopt - this.evolcost(this.count,2);
                        else
                            this.evolcost(this.count,1) = opti.cost*opti.xopt;
                        end
                    case 'OptiVMLMB'
                        this.evolcost(this.count,2) = opti.cost.alpha(end)*opti.cost.mapsCell{end}*opti.xopt;
                        this.evolcost(this.count,1) = opti.cost*opti.xopt - this.evolcost(this.count,2);
                    otherwise
                        this.evolcost(this.count,1) = gather(opti.cost*opti.xopt);
                end
                this.evolcost(this.count,:) = real(this.evolcost(this.count,:));
                str = sprintf('%s | Cost: %4.4e',str,this.evolcost(this.count,1));
                if size(this.evolcost,2) > 1
                    for kk = 2:size(this.evolcost,2)
                        str = sprintf('%s + %4.4e',str,this.evolcost(this.count,kk));
                    end
                    str = sprintf('%s = %4.4e',str,sum(this.evolcost(this.count,~isinf(this.evolcost(this.count,:)))));
                end
                %else
                %    cc = opti.cost.apply(opti.xopt);
                %    this.evolcost(this.count)=cc;
                %    str = sprintf('%s | Cost: %4.4e',str,cc);
                %end
            end
            
            if this.isgt
                if this.doSNRn
                    switch this.optiType
                        case {'BPM','BPMi'}
                            snr_out = 20*log10(this.normNtrue/norm(this.ntrue(:)...
                                - (opti.xopt(:) + this.n0)));
                        otherwise
                            snr_out = 20*log10(this.normNtrue/norm(this.ntrue(:)...
                                - this.n0*sqrt(complex(opti.xopt(:)/this.k^2 + 1))));
                    end
                else
                    snr_out=20*log10(this.normXtrue/norm(this.xtrue(:)-opti.xopt(:)));
                end
                str=sprintf('%s | SNR: %4.4e dB',str,snr_out);
                this.evolsnr(this.count)=gather(snr_out);
            end
            if this.saveX==1 || (this.saveX > 1 && mod(opti.niter,this.saveX)==0)
                this.evolxopt{this.count}=gather(opti.xopt);
            end
            this.iternum(this.count)=opti.niter;
            this.nTimer(this.count) = toc(this.start);
            str=sprintf('%s | time : %1.2f s',str,this.nTimer(this.count));
            this.count=this.count+1;
            if (mod(opti.niter,this.iterVerb)==0) || opti.niter==1
                
                fprintf('%s \n',str);
                if this.doPlot% && this.iternum(min(this.count,end))>=50
                    figure(4);clf;
                    if ismatrix(opti.xopt)
                        subplot(121);
                        imagesc(opti.xopt);axis image;colorbar;
                        subplot(122);
                        plot(this.iternum,this.evolcost);axis equal;
                    else
                        if isempty(this.pos)
                            this.pos = ceil(size(opti.xopt)/2);
                        end
                        try
                            figure(gcf);
                            subplot(224);
                            if numel(this.evolcost)==numel(this.iternum)
                                plot(this.iternum,this.evolcost);
                            else
                                plot(this.iternum(1:size(this.evolcost,1)),sum(this.evolcost,2));
                            end
                            %hold on;
                            subplot(221);
                            imagesc(squeeze(opti.xopt(:,:,this.pos(3))));
                            colorbar;axis image;title('XY');caxis(gather([min(opti.xopt(:)), max(opti.xopt(:))]));
                            subplot(222);
                            imagesc(squeeze(opti.xopt(:,this.pos(2),:)));
                            colorbar;axis image;title('XZ');caxis(gather([min(opti.xopt(:)), max(opti.xopt(:))]));
                            subplot(223);
                            imagesc(squeeze(opti.xopt(this.pos(1),:,:))');
                            colorbar;axis image;title('ZY');caxis(gather([min(opti.xopt(:)), max(opti.xopt(:))]));
                            if this.negative
                                colormap(flipud(viridis));
                            end
                            drawnow
                        catch
                            try
                                %hold on;
                                subplot(221);
                                imagesc(squeeze(opti.xopt(:,:,this.pos(3))));
                                colorbar;axis image;title('XY');caxis(gather([min(opti.xopt(:)), max(opti.xopt(:))]));
                                subplot(222);
                                imagesc(squeeze(opti.xopt(:,this.pos(2),:)));
                                colorbar;axis image;title('XZ');caxis(gather([min(opti.xopt(:)), max(opti.xopt(:))]));
                                subplot(223);
                                imagesc(squeeze(opti.xopt(this.pos(1),:,:))');
                                colorbar;axis image;title('ZY');caxis(gather([min(opti.xopt(:)), max(opti.xopt(:))]));
                                if this.negative
                                    colormap(flipud(viridis));
                                end
                                figure(gcf);
                                subplot(224);
                                if numel(this.evolcost)==this.iternum
                                    plot(this.evolcost);
                                else
                                    plot(sum(this.evolcost(:,~isinf(this.evolcost)),2));
                                end
                            end
                        end
                    end
                    drawnow;
                end
            end
            if mod(opti.niter,this.saveAll)==0
                xopt = opti.xopt;
                gam = opti.gam;
                OutOp = this;
                save(sprintf('tmp_save_run_%i%i%i_niter_%i_gam_%1.2e.mat',size(xopt),opti.niter+this.niter0,this.gam0),'xopt','gam','OutOp','-v7.3');
            end
            
        end
    end
end
