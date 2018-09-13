classdef LinOpGtild < LinOp
    % LinOpGtild: Convolution with 2D Helmholtz Green function (and its adjoint) with
    % a matrix-based convolution of the exact Green function expression
    % only sampled at the measurements positions
    % 
    % :param k: wavenumber 2 * pi / lambda * nb * dx
    % :param N: size of the volume
    % :param SenPos: structure with fields x, z for the measurements positions (origin at the center of the sample)
    % :param dN: step size of the discretized volume
    % :param stored: boolean for storing the matrix or having a lookup table (less precise)
    % :param step: step for the lookup table
    % All attributes of parent class :class:`Map` are inherited. 
    %
    % **Example** H= LinOpGtild(k,N,SenPos,dN,stored,step)
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
        k;              % wavenumber 2*pi*nb/lambda*dx
        dN;             % step size (e.g. nm / pixel), for now MUST be isotropic
        Nsensors;       % number of sensors
        SenPos;         % sensors position wrt sample center (ROI center), Nsensor x Ndims
                        % unit : distance unit (e.g. microns, nm)
                        % will be converted in discrete unit by dividing with dN (e.g. nm / (nm / pixel) = pixel
        N;              % size of the volume
        Gmat;           % Matrix encoding the convolution for the sensors position
        setCard;        % Cardinality of the set (number of pixels)
        setPos;         % Position of ROI pixels
        stored = false; % boolean for storing the matrix; false : uses a lookup table
        ndims;          
        step;           
        kernel;         
        Nelements;     
        minR;
        maxR;
        Rrange;
    end
%% Constructor    
    methods
        function this = LinOpGtild(k,N,SenPos,dN,stored,step)
            this.name = 'LinOpGtild';
            if isstruct(SenPos)
                this.ndims = length(fieldnames(SenPos)) - isfield(SenPos,'Nmeas');
                if this.ndims==2
                    SenPos = [SenPos.z(:), SenPos.x(:)];
                else
                    SenPos = [SenPos.x(:), SenPos.y(:),SenPos.z(:)];%might need to switch
                end
            end
            this.sizein = N.*ones(1,2);
            if isscalar(N)
                this.setCard = N^this.ndims;
            else
                this.setCard = prod(N);
            end
            this.Nsensors = size(SenPos,1);
            this.sizeout = [this.Nsensors,1];
            this.stored = stored;
            this.step = step;
            
            this.dN = dN;
            this.k = k;
            this.SenPos = SenPos;
            this.N = N;
            X = (0:N(1)-1) - N(1)/2 + mod(N(1),2)/2;
            X = repmat(X(:),[N(end),1]);
            Y = (0:N(end)-1) - N(end)/2 + mod(N(end),2)/2;
            Y = repmat(Y,[N(1),1]);
            this.setPos = [Y(:),X(:)];
            if this.stored
                this.Gmat = zeros_(this.Nsensors,prod(this.sizein));
            end
            this.maxR = 0;
            this.minR = inf;
            for kk = 1:this.Nsensors
                R = (repmat(SenPos(kk,:)./dN,[this.setCard,1]) - this.setPos).^2;
                R = sqrt(sum(R,2));
                if any(R(:) < dN)
                    warning('Kernel has a singularity at 0, sensors should not be included in the ROI (or too close)');
                end
                
                max_currR = max(R(:));
                if max_currR > this.maxR
                    this.maxR = max_currR;
                end
                min_currR = min(R(:));
                if min_currR < this.minR
                    this.minR = min_currR;
                end
                
                if this.stored
                    this.Gmat(kk,:) = 1i*besselh(0,1,k*R)/4;
                end
            end
            if ~this.stored
                this.minR = this.minR*this.dN;
                this.maxR = this.maxR*this.dN;
                this.Rrange = this.minR:this.step:this.maxR;
                this.kernel = 1i/4*besselh(0,1,k*this.Rrange/dN);
                this.Nelements = length(this.kernel);
            end
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected)
        function out = apply_(this,input)
            % Reimplemented from parent class :class:`LinOp`.
            if this.stored
                out = this.Gmat*input(:);
            else
                out = zeros_(this.Nsensors,1);
                for kk = 1:this.Nsensors
                     
                    R = (sqrt(sum((repmat(this.SenPos(kk,:)./this.dN,[this.setCard,1])...
                        - this.setPos).^2,2))*this.dN - this.minR)/this.step;% [dX.^2, dY.^2]
                    %R = (sqrt(sum(R,2))*this.dN - this.minR)/this.step;%sqrt(dX.^2 + dY.^2)
                    out(kk,:) = ((this.kernel(min(1 + ceil(R),this.Nelements))...
                        + this.kernel(1 + floor(R)))/2)*input(:);
                end
            end
        end
        
        function out = applyAdjoint_(this,input)
            % Reimplemented from parent class :class:`LinOp`.
            if this.stored
                out = reshape(this.Gmat'*input(:),this.sizein);
            else
                out = zeros_(this.setCard,1);
                for kk = 1:this.setCard
                    R = (sqrt(sum((this.SenPos./this.dN...
                        - repmat(this.setPos(kk,:),[this.Nsensors,1])).^2,2))*this.dN...
                        - this.minR)/this.step;% sqrt((dX.^2 + dY.^2))
                    %R = (sqrt(sum(R,2))*this.dN -this.minR)/this.step;%sqrt(dX.^2 + dY.^2)
                    out(kk,:) = conj((this.kernel(min(1 + ceil(R),this.Nelements))...
                        + this.kernel(1 + floor(R)))/2)*input(:);
                end
            end
        end
        
    end
    
end