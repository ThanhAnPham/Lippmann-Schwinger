classdef LinOpGtild_FFT < LinOp
    % LinOpGtild_FFT: Convolution with 2D Helmholtz Green function (and its adjoint) with
    % padding and get the measurements at the last plane
    % 
    % :param k: wavenumber 2 * pi / lambda * nb * dx
    % :param N: size of the volume
    % :param Next: size of the padded volume (should include the measurements plane)
    % :param NextOvs: size of the padded padded volume (should be 2 x Next)
    % :param sampleMask: Binary mask for getting the samples
    % All attributes of parent class :class:`Map` are inherited. 
    %
    % **Example** H= LinOpGtild_FFT(kdz,N,Next,NextOvs,sampleMask)
    % See also :class:`Map`
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
        G; %:class: `LinOpFreeSpaceKernel`
        Zpad;
        Crp;
    end
    
    methods
        function this = LinOpGtild_FFT(kdz,N,Next,NextOvs,sampleMask)
            this.name = 'LinOpGtild_FFT';
            this.G = LinOpFreeSpaceKernel(kdz,Next,NextOvs);
            this.sizein = N.*ones(1,2);
            this.Zpad = LinOpSelectorPatch(Next.*ones(1,2),...
                1 + Next/2 - this.sizein/2, Next/2 + this.sizein/2)';
            this.Crp = LinOpSelector(sampleMask);
            this.sizeout = this.Crp.sizeout;
        end
    end
    methods (Access = protected)
        function out = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            out = this.Crp*(this.G*(this.Zpad*x));
        end
        
        function out = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            out = this.Zpad'*(this.G'*(this.Crp'*x));
        end
    end
end