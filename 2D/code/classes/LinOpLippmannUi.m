classdef LinOpLippmannUi < LinOp
    % LinOpLippmannUi: "forward" operator in Lippman-Schwinger model to obtain the incident wavefield u_in
    % this.apply_(x) computes (\Id - G \diag(f))*x, which gives u_in in Lippmann-Schwinger equation
    % :param G: :class:`LinOpFreeSpaceKernel`
    % :param f: :class:`LinOpDiag` operator of the scattering potential 
    %
    % **Example** H = LinOpLippmannUi(G,f)
    %
    % See also :class:`LinOp`, :class:`Map`
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
    
    properties
        G;
    end
    properties (SetAccess = protected,GetAccess = public)
        f;
    end
    %% Constructor    
    methods
        function this = LinOpLippmannUi(G,f)
            this.name = 'LinOpLippmannUi';
            this.sizein = G.sizein;
            this.sizeout = this.sizein;
            this.G = G;
            if isa(f,'LinOpDiag')
                this.f = f;
            elseif ~any(size(f) - this.sizein)
                this.f = LinOpDiag([],f);
            else
                error('f must either be a LinOpDiag or the diagonal of the matrix of size [%i,%i]',this.sizein);
            end
        end
    end
    methods (Access = protected)
        function u_in = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            u_in = x - this.G.apply(this.f*x);
        end
        
        function xout = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            xout = x - this.f'*(this.G'*x);
        end
        function updateF(this,newf)
            if isa(newf,'LinOpDiag')
                this.f = newf;
            elseif ~any(size(newf) - this.sizein)
                this.f = LinOpDiag([],newf);
            else
                error('f must either be a LinOpDiag or the diagonal of the matrix of size [%i,%i]',this.sizein);
            end
        end
        
    end
end