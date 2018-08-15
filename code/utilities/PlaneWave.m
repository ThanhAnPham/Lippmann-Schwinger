classdef PlaneWave < handle
    % 2D Plane wave object. The plane wave is defined by the wave vector. Here
    % we use a unitless wave vector defined in the variable "Frequency".
    % Frequency has the X component is the first dimension and the Y
    % component in the second dimension.
    %
    %T0D0 :
    %   Extend the code to non square case and 3D case.
    %
    %History :
    %   18/02/2016 : creation, Luc Zeng
    %
    properties
        Frequency              %
    end
    
    
%%    
    methods
        function obj = PlaneWave(Frequency)
            obj.Frequency = Frequency;
        end
        function u = Apply(this,Nx,Nz)
            % Propagates the wave over a distance N
            Vectx = (0:Nx-1) - Nx/2;
            Vectz = (0:Nz-1) - Nz/2;
            [Z,X] = meshgrid(Vectz,Vectx);
            u = exp(1i*(this.Frequency(1)*Z + this.Frequency(2)*X));
        end      
    end
    
end