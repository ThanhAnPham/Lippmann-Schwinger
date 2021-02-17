classdef LinOpTiltedPropagator <  LinOp
    %% LinOpTiltedPropagator : Discrete Fresnel operator
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = LinOpTiltedPropagator(lambda, n0, z,dxy,sz)
    % Fresnel transform operator with
    %           lambda   wavelength            [m]
    %         n0      % refractive index of the medium
    %         z       % depth of propagation  [m]
    %         dxy     % pixel size            [m]
    %         if the option 'FeitFleck' is set then it will use the Feit and Fleck model of propagation instead:
    %           M. D. Feit and J. A. Fleck, ?Bean nonparaxiality, filament formaform, and beam breakup in the self-focusing of optical beams,? J. Opt. Soc. Am. B, vol. 5, pp. 633? 640, March 1988.
    %         if the option 'DontUseComplex' is set, complex are repsented as    an extra dimension of size 2 containning Real and imagenary parts of x
    %
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % See also LinOp
    
    %     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
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
    
    properties (Constant=true)
        PUPIL = 3   % Fresnel propagation through an objective of pupil of radius R
        FEITFLECK = 2 % use the Feit and Fleck model of propagation instead:
        % M. D. Feit and J. A. Fleck, ?Bean nonparaxiality, filament formaform, and beam breakup in the self-focusing of optical beams,? J. Opt. Soc. Am. B, vol. 5, pp. 633? 640, March 1988.
        AS = 1 % Use Angular Spectrum method
    end
    properties (SetAccess = protected,GetAccess = public)
        lambda  % wavelenght            [m]
        n0      % refractive index of the medium
        k0      % wavenumber in vaccum  [m-1]
        k       % wavenumber in the medium [m-1]
        z       % depth of propagation  [m]
        dxy     % pixel size            [m]
        Nx      % number of pixels along X
        Ny      % number of pixels along Y
        Nt      % number of angle
        pad     % size with padding  % NOT IMPLEMENTED
        Fu      % Fresnel function along the u axis
        Fv      % Fresnel function along the v axis
        ephi       % Fresnel function
        type = 0; % if true
        theta   % tilt angle
        sizeinC;
        sizeoutC;
        phase;
        mod;
        lp;
        %local_freq;
        mu;
        mv;
    end
    methods
        function this = LinOpTiltedPropagator(lambda, n0, z,dxy, sz, theta, varargin)
            
            this.name ='LinOpTiltedPropagator';
            this.isInvertible=false;
            
            this.lambda = lambda;
            
            assert(isPositiveScalar(n0),'The refractive index n0 should be a positive scalar');
            this.n0 = n0;
            
            this.z = z;
            
            assert(isPositiveScalar(dxy),'The pixel size dxy should be a positive scalar');
            this.dxy = dxy;
            
            assert(issize(sz) && (length(sz)==2),'The input size sz should be a conformable  to size(2D) ');
            
            this.theta = theta;
            this.sizein = sz ;
            this.sizeout = sz;
            this.sizeinC = this.sizein;
            this.sizeoutC = this.sizeout;
            
            this.Nx = sz(1);
            this.Ny = sz(2);
            
            for c=1:length(varargin)
                switch varargin{c}
                    case('Pupil')
                        this.type = this.PUPIL;
                        NA = varargin{c+1};
                    case('FeitFleck')
                        this.type = this.FEITFLECK;
                    case('AS')
                        this.type = this.AS;
                end
            end
            
            Nx_ = this.Nx;
            Ny_ = this.Ny;
            
            
            %  frequency grid
            v = 1./(Nx_ * dxy) *( [0:ceil( Nx_/2)-1, -floor( Nx_/2):-1]' -   n0*dxy *   sin(theta(1)/180*pi)/lambda.*Nx_);
            u = 1./(Ny_ * dxy) * ([0:ceil( Ny_/2)-1, -floor( Ny_/2):-1] -    n0*dxy *   sin(theta(2)/180*pi)/lambda.*Ny_);
            
            
            [mu,mv] =  meshgrid(  u.^2,  v.^2);
            Mesh = mv + mu;
            
            if (max(Mesh(:))>(n0/ lambda)^2)
                mod =  (Mesh<= (n0/ lambda)^2);
                Mesh(~mod)=0.;
            else
                mod = 1;
            end
            switch  this.type
                case  this.PUPIL % pupil
                    ephi_ = mod.*(Mesh<=(NA/ lambda)^2);
                case  this.AS % Angular spectrum
                    ephi_ =  mod.*exp(2i* pi *  z.* real(sqrt((  n0/ lambda)^2- Mesh)));
                case  this.FEITFLECK
                    ephi_=  mod.*exp(-2i* pi * z.*lambda / n0 * Mesh ./ real(1 + sqrt(1 - (lambda/n0)^2 *Mesh)));
                otherwise
                    % separable Fresnel function
                    ephi_ =mod.* exp(-1i* pi *  z.* lambda / n0 .*Mesh);
            end
            this.mod = mod;
            this.phase = 2i* pi * real(sqrt((  n0/ lambda)^2- Mesh));
            this.mu = mu;
            this.mv = mv;
            %this.local_freq = z*sqrt(mu)./(sqrt((  n0/ lambda)^2- Mesh)) <= sz(1)*dxy...
            %    & z*sqrt(mv)./(sqrt((  n0/ lambda)^2- Mesh)) <= sz(2)*dxy;
            this.lp = sinc(dxy*v(:))*sinc(dxy*u(:))';%*dxy^2;
            this.ephi = ephi_;%.*this.lp.*this.local_freq;
        end
    end
    methods (Access = protected)
        function y = apply_(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d]',this.sizein);
            y = ifft2( this.ephi  .*   fft2(x));
        end
        function y = applyAdjoint_(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d]',this.sizeout);
            y =  ifft2(conj(this.ephi) .*  fft2(x));
        end
        
    end
end
