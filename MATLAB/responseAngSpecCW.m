function psf = responseAngSpecCW(x, z, elemSpace, apod_P_f, steerAng, focDepth, f, c)
%
% Continuous wave spatial sensitivity (on rx or tx) from xdcr array
%
% x: lateral dimension in mm
% z: axial dimension in mm
% elemSpace: element spacing in mm
% apod_P_f: apodization + pulse spectrum vs frequency (rows) and array element (columns)
% steerAng: steering angle (deg)
% focDepth: focus depth (mm)
% f: row vector array of frequencies in pulse (MHz)
% c: speed of sound (mm/usec)

disp('CW Response Calculation START'); tic;

% X-Coordinate of Apertures
Nelem = size(apod_P_f, 2); 
xpos = (-(Nelem-1)/2:(Nelem-1)/2)*elemSpace; 
apod_P_f_x = interp1(xpos, apod_P_f', x, 'spline', 0);
apod_P_f_x = reshape(apod_P_f_x, [numel(x), numel(f)]);

% Calculating Time Delays
delayIdeal=(sqrt(x.^2+focDepth^2-2*focDepth*x*sind(steerAng))-focDepth)/c;

% Aperture Function and Its Fourier Transform
a_x = apod_P_f_x .* exp(-1i*2*pi*delayIdeal'*f);
A_Kx = fftshift(fft(a_x, [], 1), 1); % Fourier Transform Along Lateral Axis

% Assembly of Continuous Wave Response By Angular Spectrum
A_Kx_z = repmat(A_Kx, [1, 1, numel(z)]); % Stack in Axial Direction
dx = mean(diff(x)); nxFFT = numel(x); % Spatial Sampling for FFT Axis
kx = (1/(2*dx))*(-1:2/nxFFT:1-2/nxFFT); % FFT Axis for Lateral Spatial Frequency
[F, Kx, Z] = meshgrid(f, kx, z); % Create Grid for Angular Spectrum
Kz = sqrt((F/c).^2 - Kx.^2); % Axial Spatial Frequency
H = exp(1i*2*pi*Kz.*Z); % Propagation Filter in Spatial Frequency Domain
psf_Kx_z = A_Kx_z .* H; % Apply Propagation Filter
psf = ifft(ifftshift(psf_Kx_z, 1), [], 1); % Invert Lateral Spatial FFT

% Reshape For Expected Output
psf = permute(psf, [3, 1, 2]);
toc; disp('CW Response Calculation END');

end
