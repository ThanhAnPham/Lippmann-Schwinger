classdef LippmannUi3D < LinOp
    % LippmannUi3D : This operator applies (I - G*f)*x. If x=u, the output
    % is equal to u_{in}.
    % 
    % :param G: FreeSpaceKernel3D (Green function)
    % :param f: scattering potential (optional), can be set after construction.
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Note**: This operator is used during LSm forward computation
    %
    % **References**
    %
    % [1] Three-Dimensional Optical Diffraction Tomography with Lippmann-Schwinger Model.,
    %     IEEE Transactions on Computational Imaging (2020), T-a. Pham, E. Soubies, 
    %     A. Ayoub, J. Lim, D. Psaltis, and M. Unser.
    %
    % **Example** H=LippmannUi3D(G,f)
    %
    % See also :class:`Map` :class:`LinOp`
    
    properties (Access = public)
        G;
        f = [];
    end
    methods
        function this = LippmannUi3D(G,varargin)
            this.name = 'LippmannUi3D';
            this.sizein = G.sizein;
            this.sizeout = this.sizein;
            this.G = G;
            if nargin > 1
                this.f = LinOpDiag(size(varargin{1}),varargin{1});
            end
        end
    end
    methods (Access = protected)
        function u_in = apply_(this,x)
            u_in = x - this.G*(this.f*x);
        end
        
        function xout = applyAdjoint_(this,x)
            xout = x - this.f'*(this.G'*(x));
        end
        function updateF(this,newf)
            this.f = LinOpDiag(newf);
        end
        
    end
end