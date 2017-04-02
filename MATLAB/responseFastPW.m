function psf_t = responseFastPW(x, z, elemSpace, apod_P_f, steerAng, focDepth, f, c, t)
%
% Pulse wave spatial response (on rx or tx) from xdcr array
%
% x: lateral dimension in mm
% z: axial dimension in mm
% elemSpace: element spacing in mm
% apod_P_f: apodization + pulse spectrum vs frequency (rows) and array element (columns)
% steerAng: steering angle (deg)
% focDepth: focus depth (mm)
% P_f: pulse as a function of frequency 
% f: row vector array of frequencies in pulse (MHz)
% c: speed of sound (mm/usec)
% t: measurement time vector in usec (one-way)

tic

% Get K-Space Axes
dx = mean(diff(x)); nxFFT = numel(x); % Spatial Sampling for FFT Axis
kx = (1/(2*dx))*(-1:2/nxFFT:1-2/nxFFT); % FFT Axis for Lateral Spatial Frequency
dz = mean(diff(z)); nzFFT = numel(z); % Spatial Sampling for FFT Axis
kz = (1/(2*dz))*(-1:2/nzFFT:1-2/nzFFT); % FFT Axis for Lateral Spatial Frequency
[Kx_old, F] = meshgrid(kx, f); % Get Temporal + Lateral Spatial Frequency Mesh
[Kx, Kz] = meshgrid(kx, kz); % Get Spatial Frequency Mesh

% X-Coordinate of Apertures
Nelem = size(apod_P_f, 2); 
xpos = (-(Nelem-1)/2:(Nelem-1)/2)*elemSpace; 
apod_P_f_x = interp1(xpos, apod_P_f', x, 'spline', 0);
apod_P_f_x = reshape(apod_P_f_x, [numel(x), numel(f)]);

% Calculating Time Delays
delayIdeal = -(sqrt(x.^2+focDepth^2-2*focDepth*x*sind(steerAng))-focDepth)/c;

% Aperture Function and Its Fourier Transform
a_x = apod_P_f_x .* exp(-1i*2*pi*delayIdeal'*f);
A_Kx_F = fftshift(fft(a_x, [], 1), 1)'; % Fourier Transform Along Lateral Axis

% Convert to 2-D Kspace Image
K = sqrt(Kx.^2 + Kz.^2);
A_Kx_Kz = interp2(Kx_old, F, A_Kx_F, Kx, c*K, 'linear', 0);
A_Kx_Kz = A_Kx_Kz .* (1i./(4*pi*K));
A_Kx_Kz(Kz < 0) = 0; % Get Rid of Backward Going Waves

% Get One-Way Response for Different Time Points
psf_t = zeros([size(A_Kx_Kz), numel(t)]);
for t_idx = 1:numel(t)
    A_Kx_Kz_t = A_Kx_Kz .* exp(-1i*2*pi*c*K*t(t_idx));
    psf_t(:,:,t_idx) = ifft2(ifftshift(A_Kx_Kz_t)); % Invert 2-D Spatial FFT
end

toc

end

