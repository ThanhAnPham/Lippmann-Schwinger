classdef OutputOptiLipp < handle
    %% OutputOptiLipp : class for object OutputOpti objects for tomography
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
        name = 'OutputOptiLipp'% name of the optimization algorithm
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
    end
    
    methods
        %% Constructor
        function this=OutputOptiLipp(computecost,xtrue,iterVerb,k,n0,optiType,varargin)
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
            elseif nargin>=6
                this.computecost=computecost;
                this.xtrue=xtrue;
                this.iterVerb=iterVerb;
                this.doSNRn = true;
                this.optiType = optiType;
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
            if nargin > 6
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
                
                switch opti.name
                    case 'Opti FBS'
                        this.evolcost(this.count,1) = opti.F*opti.xopt;
                        this.evolcost(this.count,2) = opti.G*opti.xopt;
                        str = sprintf('%s | Cost: %4.4e + %4.4e= %4.4e',str,this.evolcost(this.count,:),sum(this.evolcost(this.count,:)));
                    otherwise
                        cc = opti.cost.eval(opti.xopt);
                        this.evolcost(this.count)=cc;
                        str = sprintf('%s | Cost: %4.4e',str,cc);
                end
            end
            
            if this.isgt
                if this.doSNRn
                    snr_out = 20*log10(this.normNtrue/norm(this.ntrue(:)...
                        - this.n0*sqrt(opti.xopt(:)/this.k^2 + 1)));
                else
                    snr_out=20*log10(this.normXtrue/norm(this.xtrue(:)-opti.xopt(:)));
                end
                str=sprintf('%s | SNR: %4.4e dB',str,snr_out);
                this.evolsnr(this.count)=snr_out;
            end
            this.evolxopt{this.count}=opti.xopt;
            this.iternum(this.count)=opti.niter;
            this.nTimer(this.count) = toc(this.start);
            str=sprintf('%s | time : %1.2f s',str,this.nTimer(this.count));
            this.count=this.count+1;
            if (mod(opti.niter,this.iterVerb)==0) || opti.niter==1
                
                fprintf('%s \n',str);
                if this.doPlot% && this.iternum(min(this.count,end))>=50
                    try
                        figure(gcf);
                        subplot(121);
                        if numel(this.evolcost)==numel(this.iternum)
                            plot(this.iternum-1,this.evolcost);
                        else
                            plot(this.iternum(1:size(this.evolcost,1))-1,sum(this.evolcost,2));
                        end
                        %hold on;
                        subplot(122);
                        imagesc(opti.xopt);
                        colorbar;axis image;
                        if this.negative
                            colormap(flipud(gray));
                        end
                        drawnow
                    catch
                        try
                            figure(gcf);
                            subplot(121);
                            if numel(this.evolcost)==this.iternum
                                plot(this.evolcost);
                            else
                                plot(sum(this.evolcost,2));
                            end
                            %hold on;
                            subplot(122);
                            imagesc(opti.xopt);
                            colorbar;axis image;
                            if this.negative
                                colormap(flipud(gray));
                            end
                        end
                    end
                end
            end
            
            
        end
    end
end
