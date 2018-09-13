classdef OpLipp <  Map
    % OpLipp : Lippmann-Schwinger based forward operator
    % 
    % : param G:
    % Please refer to the Map superclass for documentation
    %%    Copyright (C) 2018
    %     T. Pham thanh-an.pham@epfl.ch
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
    
    properties (SetAccess = public,GetAccess = public)
        G;                 % :class:`LinOpFreeSpaceKernel`
        CG;                % :class:`OptiConjGrad`
        MGtild;            % :class: `LinOpGtild` or :class: `LinOpGtild_FFT`
        siz_ext;           % size
        uin;               %
        Niter = 120;       %
        xtol = 1e-8;       %
        u;                 %
        f;                 %
        Muin;              %
        totalfield = true; %
    end
    %% Constructor    
    methods
        function this = OpLipp(G,Gtil,uin,Niter,xtol,totalfield,uem)
            this.name = 'OpLipp';
            this.sizein = G.sizein;
            this.sizeout = Gtil.sizeout;
            this.MGtild = Gtil;
            this.Muin = uem;
            
            if any(size(this.Muin) - this.sizeout)
                this.MGtild = LinOpShape(this.sizeout,size(this.Muin))*Gtild;
                this.sizeout = this.MGtild.sizeout;
            end
            
            this.G = G;
            this.siz_ext = size(uin);
            %Related to a scenario where uin was provided in an extended
            %area that includes the sensors. If this is not the case, Crp
            %is equal to identity.
            Crp = LinOpSelectorPatch(this.siz_ext,...
                1 + (this.siz_ext - this.sizein)/2,(this.siz_ext + this.sizein)/2);
            this.uin = Crp*uin;
            this.Niter = Niter;
            this.xtol = xtol;
            this.totalfield = totalfield;
            this.f = zeros_(this.sizein);
            this.u = Crp*uin;
        end
    end
    methods (Access = protected)
        function [y, utot] = apply_(this,x)
            % Reimplemented from parent class :class:`Map`.
            forw = LinOpLippmannUi(this.G, x);
            
            A = forw'*forw;
            b = forw'*this.uin;
            
            this.CG = OptiConjGrad(A,b);
            this.CG.verbose = false;
            this.CG.OutOp = OutputOpti(1,[],1);
            this.CG.ItUpOut = 0;
            this.CG.maxiter = this.Niter;
            this.CG.CvOp = TestCvgStepRelative(this.xtol);
            
            this.CG.run(this.uin);
            utot = this.CG.xopt;
            
            %fprintf('Niter in forward : %i. ',this.CG.niter);
            
%             if any(isnan(utot(:))) %because of CG not stopping at the beginning
%                 fprintf('NaN in CG (if empty volume, expected)\n');
%                 utot = this.uin;
%             end
            
            if this.totalfield
                y = this.Muin + this.MGtild*(utot.*x);
            else
                y = this.MGtild*(utot.*x);
            end
            
            this.u = utot;
            
        end
        
        function out = applyJacobianT_(this,x,v)
            % Reimplemented from parent class :class:`Map`.
            forw = LinOpLippmannUi(this.G,v);
            GtMtx = this.MGtild'*x;
            A = forw*forw';
            b = forw*(v.*GtMtx);
            
            CGgrad = OptiConjGrad(A,b);
            CGgrad.verbose = false;
            CGgrad.OutOp = OutputOpti(1,[],1);
            CGgrad.ItUpOut = 0;
            CGgrad.maxiter = this.Niter;
            CGgrad.CvOp = TestCvgStepRelative(this.xtol);
            
            CGgrad.run(A'*b);
            
            if ~any(isnan(CGgrad.xopt(:)))
                out = real(conj(this.u).*(this.G'*CGgrad.xopt + GtMtx));
            else
                out = real(conj(this.u).*GtMtx);
            end
            
        end
        
    end
end
