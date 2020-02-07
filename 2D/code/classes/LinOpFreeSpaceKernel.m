classdef LinOpFreeSpaceKernel < LinOp
    % LinOpFreeSpaceKernel: Convolution with 2D Helmholtz Green function (and its adjoint) with
    % padding.
    % 
    % :param k: wavenumber 2 * pi / lambda * nb * dx
    % :param N: size of the volume
    % :param NOvs: size of the padded volume
    % All attributes of parent class :class:`Map` are inherited. 
    %
    % **Example** H= LinOpFreeSpaceKernel(k,N,NOvs)
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
        k;    % Wavenumber $\frac{2 \pi n_b}{\lambda} dz
        mtf;  % Fourier transform of the Green function
        Crp;  % Cropping operator
        NOvs; % Size of the padded image
    end
%% Constructor    
    methods
        function this = LinOpFreeSpaceKernel(k,N,NOvs)
            this.sizein = N.*ones(1,2);
            this.sizeout = this.sizein;
            [rho, maxI] = max(N);
            this.NOvs = NOvs;
            Rscalar = max(NOvs);% R Radius of the OVS signal, rho support of the OVS object
            
            this.Crp = LinOpSelectorPatch(this.NOvs.*ones(1,2),[1,1],this.sizein);
            
            
            this.k = k;
            Vect = (0:Rscalar-1) - Rscalar/2;
            [X,Y] = meshgrid(Vect,Vect);
            R = sqrt(X.^2 + Y.^2);
            
            
            s = 2*pi*R/length(R);
            J0 = @(x) besselj(0,x);J1 = @(x) besselj(1,x);
            H0 = @(x) besselh(0,1,x);H1 = @(x) besselh (1,1,x);
            
            this.mtf = 1 ./( s.^2 - (k)^2 );
            this.mtf = this.mtf .*( 1 + 1i*pi*rho...
                * ( s/2.*H0(k*rho).*J1(rho*s)...
                - k/2*J0(rho*s)*H1(k*rho) )) ;
            
            pos0 = s.^2 ==  k^2;
            this.mtf(pos0) = pi^2/(4*k*Rscalar^2*1i)...
                * ( Rscalar*k/pi*H0(k*rho).*J0(k*rho) ...
                + k*rho*J1(k*rho)*H1(k*rho));
            
            this.mtf = ifftshift(this.mtf);
            
            this.mtf = fftshift(this.mtf);
            switch maxI
                case 1
                    this.mtf = this.mtf(:,1 + end/2 - NOvs(end)/2:end/2 + NOvs(end)/2);
                case 2
                    this.mtf = this.mtf(1 + end/2 - NOvs(1)/2:end/2 + NOvs(1)/2,:);
            end
            this.mtf = ifftshift(this.mtf);
        end
        
    end
    
    methods (Access = protected)
        function out = apply_(this,input)
            % Reimplemented from parent class :class:`LinOp`.
            out = this.Crp*ifft2(fft2(this.Crp'*input) .* this.mtf);
        end
        
        function out = applyAdjoint_(this,input)
            % Reimplemented from parent class :class:`LinOp`.
            out = this.Crp.apply(ifft2(fft2(this.Crp'*input) .* conj(this.mtf)));
        end
        function y = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.NOvs==this.sizein(1)
                y = ifft2(fft2(x).*abs(this.mtf).^2);
            else
                y = this.applyAdjoint(this.apply(x));
            end
        end
        function y = applyHHt_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.NOvs==this.sizein(1)
                y = ifft2(fft2(x).*abs(this.mtf).^2);
            else
                y = this.apply(this.applyAdjoint(x));
            end
        end
    end
    
end