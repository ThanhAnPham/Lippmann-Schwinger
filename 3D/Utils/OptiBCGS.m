classdef OptiBCGS < Opti
    % BiConjugate Gradients Stabilized Method optimization algorithm which
    % solves the linear system \\(\\mathrm{Ax=b}\\) by approximating the
    % solution with a vector in a Krylov subspace with minimal residual.
    % Use Matlab built-in function.
    %
    % :param A: square matrix :class:`LinOp`
    % :param b: right-hand term
    %
    % All attributes of parent class :class:`Opti` are inherited.
    %
    %
    % **Example** opti=OptiBCGS(A,b,OutOp)
    %
    % See also :class:`Opti`, :class:`OutputOpti` :class:`Cost`
    
    %%    Copyright (C)
    %     2020 T.-A. Pham thanh-an.pham@epfl.ch
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
    
    % Protected Set and public Read properties
    properties (SetAccess = protected,GetAccess = public)
        A;  % Linear operator
        b;  % right hand side term
        Minv;
        flag0;
        evolcost;
        iternum;
    end
    % Full protected properties
    properties (SetAccess = protected,GetAccess = protected)
        init_siz;
    end
    % Full public
    properties
        xtol = 1e-10;
    end
    
    methods
        %% Constructor
        function this=OptiBCGS(A,b)
            this.name='Opti BiConjugate Gradients Stabilized Method';
            this.init_siz = A.sizeout;
            if A.sizein(1)~=prod(A.sizein)
                A = LinOpShape(A.sizeout,[prod(A.sizeout),1])*A*LinOpShape([prod(A.sizein),1],A.sizein);
            end
            this.A = A;
            this.b=b(:);
            this.cost=CostL2Composition(CostL2([],this.b),A);
            
        end
        %% Set data b
        function set_b(this,b)
            % Set the right-hand side \\(\\mathrm{b}\\)
            b = b(:);
            assert(isequal(this.A.sizeout,size(b)),'A sizeout and size of b must be equal');
            this.b=b;
            this.cost.H1.y=b;
        end
        function set_Minv(this,Minv)
            % Set the preconditioner \\(\\mathrm{Minv}\\) such that \\(\\mathrm{A Minv x=b}\\)
            % A Map Minv must be provided
            
            if Minv.sizein(1)~=prod(Minv.sizein) %must be vector input and output
                Minv = LinOpShape(Minv.sizeout,[prod(Minv.sizeout),1])*Minv*LinOpShape([prod(Minv.sizein),1],Minv.sizein);
            end
            
            assert(isequal(this.A.sizein,Minv.sizeout),'A sizein and Minv sizeout must be equal');
            assert(isequal(size(this.b),Minv.sizein),'size of b and Minv sizein must be equal');
            
            this.Minv = Minv;
        end
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            if isempty(x0) % To restart from current state if wanted
                warning('Starting point x0 set as zero vector');
            else
                this.xopt;
            end
            %warning('Stopping criterion corresponds to the one of built-in function bicgstab');
            %warning('Attribute OutOp is not used');
        end
        
        function run(this,x0)
            % Run the algorithm.
            %
            % :param x0: initial point in \\(\\in X\\), if x0=[] restarts from the current value :attr:`xopt`.
            %
            % **note**: this method does not return anything, the result being stored in public attribute :attr:`xopt`.
            
            this.initialize(x0);
            tstart=tic;
            this.starting_verb();
            if isempty(this.Minv)
                Mprec = this.Minv;
            else
                Mprec = @(x) this.Minv*x;
            end
            [this.xopt,this.flag0,~,this.niter,this.evolcost] = bicgstab(@(x) this.A*x,this.b(:),...
                this.xtol,this.maxiter,Mprec,[],this.xopt(:));
            this.xopt = reshape(this.xopt,this.init_siz);
            this.time=toc(tstart);
            this.iternum = 1:0.5:length(this.evolcost);
            this.ending_verb();
        end
    end
end
