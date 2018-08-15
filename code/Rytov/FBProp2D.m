function [n,F,FF,FF_mask] = FBProp2D(approximation,uM,uin,unwrapped_imag_log_ur,simP,varargin)

lambda0 = simP.lambda0;%0.406;			% Vacuum waelength of the laser (laser diode with approx. 100 micron coherence length)
n0 = simP.n0;					% Background refractive index (immersion oil)
NA = 1.2;%1.2					% Numerical aperture of the imaging objective (illumination and detection is symmetric)
if nargin > 5 && ~isempty(varargin{1})
    thetas = varargin{1};
    Ntheta = length(thetas);
else
    thetas = simP.thetas;
    Ntheta = simP.Ntheta;
end
if nargin > 6
   kz_inc = varargin{2};
   kx_inc = varargin{3};     
end
% Grid definition

Nx = simP.Nx;					% Number of pixels in the grid
Nz = simP.Nz;
if Nx~=Nz
    Nz = Nx;
end
dx = simP.dx;%camera_pixel_size/magnification;	% Effective pixel size in the object plane
dz = simP.dz;

Lx = dx*Nx;					% Width of the computation window
Lz = dz*Nz;

dkx = 2*pi/Lx;				% Pixel size in the Fourier domain
x = (-Nx/2:Nx/2-1)*dx;		% Coordinate vector in the spatial domain
kx = (-Nx/2:Nx/2-1)*dkx;	% Coordinate vector in the Fourier domain

dkz = 2*pi/Lz;

X = x';
Kx = kx';
K = sqrt(Kx.^2);

k0 = 2*pi/lambda0;			% Vaccum k-vector
k = k0*n0;					% k-vector in the background medium of refractive index n0
kmax = k0*NA;				% Largest k-vector transmitted by an optical system of given NA
objective_aperture_ind = find(double(K < kmax));	% Array containing the indices of the pixels within the aperture of radius kmax
Nind = length(objective_aperture_ind);				% Number of pixels within the aperture of radius kmax

% Incident plane waves vectors
if ~exist('kx_inc','var')
    kx_inc = k*sin(thetas);
    kz_inc = sqrt(k.^2 - kx_inc.^2);	% Because the object is immersed in background index n0, we have kx^2 + ky^2 + kz^2 = k^2
end

s0 = [kx_inc(:) kz_inc(:)]/k;		% Normalized k-vectors of the incident waves

sx = Kx(objective_aperture_ind)/k;		% Normalized k-vectors of the detection directions. The minus sign is because of the orientation of the detection plane
sz = real(sqrt(1 - sx.^2));
s = [sx sz];

%****************

% Loop over the view in order to populate the Fourier spectrum F

%propagation_window = exp(-((X(:, :, 1).^2)/(0.8*Lx/2)^2).^10);

%roi_center_x = 512;
%roi_center_y = 512;


FF = zeros(Nx, Nz);							% FF will contain the 3D Fourier transform of the reconstructed scattering potential F
FF_mask = zeros(Nx, Nz);					% Set of the voxels that will be populated in F (for monitoring only)
cumulated_contributions = zeros(Nx, Nz);	% Arrays in which we put the counts in case a voxel gets assigned more than once

for ind_view = 1:Ntheta
    
    %unwrapped_imag_log_ur = unwrapped_imag_log_ur(roi_center_y-Ny_raw/2+1:roi_center_y+Ny_raw/2, roi_center_x-Nx_raw/2+1:roi_center_x+Nx_raw/2);
    %u = u(roi_center_y-Ny_raw/2+1:roi_center_y+Ny_raw/2, roi_center_x-Nx_raw/2+1:roi_center_x+Nx_raw/2);
    %ui = ui(roi_center_y-Ny_raw/2+1:roi_center_y+Ny_raw/2, roi_center_x-Nx_raw/2+1:roi_center_x+Nx_raw/2);
    
    %unwrapped_imag_log_ur = unwrapped_imag_log_ur(1:2:end, 1:2:end);
    %u = u(1:2:end, 1:2:end);
    %ui = ui(1:2:end, 1:2:end);
    u = uM(:,ind_view);
    ui = uin(:,ind_view);
    ur = u./ui;
    us = u - ui;
    
    if strcmpi(approximation, 'Born')
        % Born approximation
        um_0 = us;
    elseif strcmpi(approximation, 'Rytov')
        % Rytov approximation
        imag_log_ur = unwrapped_imag_log_ur(:,ind_view);%imag(log(ur));
        real_log_ur = real(log(ur));
        
        imag_log_ur(isnan(imag_log_ur)) = 0;
        real_log_ur(isnan(real_log_ur)) = 0;
        
        %clf, imagesc(imag_log_ur), caxis([-10 10]), colorbar, drawnow;%, pause;
        %ind_view
        %continue;
        
        um_0 = real_log_ur + 1i*imag_log_ur;
        % We modulate the field by the incident plane wave.
        %Note that the field are initially corrected for the illumination plane wave so that we have to put it back here
        um_0 = um_0 .* exp(1i*(kx_inc(ind_view)*X));
    else
        error('Invalid approximation')
    end
    
    %um_0 = um_0 .* propagation_window;
    far_field = fftshift(fft(fftshift(um_0)));
    far_field = -1i*k*sz/(2*pi*dz).*far_field(objective_aperture_ind);
    
    S = s - repmat(s0(ind_view, :), Nind, 1);		% According to Wolf formula, S is the Fourier component that is probed by wave s under illumination s0, with S = s - s0
    K_nearest_neighbor = round(S./repmat([dkx dkz]/k, Nind, 1));	% We calculate the normalized S vector in pixel units. That will be the surface of the Ewald sphere. We do simple nearest neighbor interpolation
    px = K_nearest_neighbor(:, 1) + round((Nx-1)/2) + 1;                  % More compact notation and set the center of the coordinate system to the corner of the array
    pz = K_nearest_neighbor(:, 2) + round((Nz-1)/2) + 1;
    
    fourier_space_indices = (pz-1)*Nx + px;	% Calculate the indices in the Fourier domain, 
    %sub2ind([Nx,Nz],px,pz);%
    FF(fourier_space_indices) = FF(fourier_space_indices) + far_field;% Assign the far field data to the surface of the Ewald sphere
    
    FF_mask(fourier_space_indices) = 1;	% Set FF_mask to 1 on the Ewald sphere
    cumulated_contributions(fourier_space_indices) = cumulated_contributions(fourier_space_indices) + 1;	% We count an increment of 1 for each voxel that we have just populated
    
end


FF = FF./max(cumulated_contributions, 1);	% We divide the populated Fourier spectrum by the cumulated contributions

F = fftshift(ifftn(fftshift(FF)));	% We reconstruct the scattering potentatial defined as F = k^2/(4pi)(n^2 - n0^2)

n = real(sqrt(4*pi*F/k0^2 + n0^2));	% Retrieve the refractive index from the scattering potential

end

