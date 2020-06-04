classdef PlaneWave3D < handle
    % PlaneWave3D : 3D Plane wave object. The plane wave is defined by the wave vector.
    %               We use a unitless wave vector defined in the variable "Frequency".
    %               Frequency has the X component is the first dimension, the Y
    %               component in the second dimension and the Z
    %               component in the third dimension.
    %
    properties
        Frequency              %
    end
    
    
%%    
    methods
        function obj = PlaneWave3D(Frequency)
            obj.Frequency = Frequency;
        end
        function u = Apply(this,Nx,Ny,Nz)
            % Propagates the wave over a distance N
            Vectx = (1:Nx) - Nx/2;
            Vecty = (1:Ny) - Ny/2;
            Vectz = (1:Nz) - Nz/2;
            [X,Y,Z] = meshgrid(Vectx,Vecty,Vectz);
            u = exp(1i*(this.Frequency(1)*X + this.Frequency(2)*Y + this.Frequency(3)*Z));
        end      
    end
    
end