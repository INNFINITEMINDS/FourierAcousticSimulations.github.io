clear
clc

% Pulse Definition
fc = 3.0; fracBW = 0.5; Nf = round(fracBW*512); 
f = ((-Nf/2:Nf/2-1)/Nf)*4*fc*fracBW+fc; % MHz
f=f(f>0); P_f = exp(-pi*((f-fc)/(fracBW*fc)).^2);

% Aperture Definition
c = 1.54; % mm/usec
lambda = c/fc; elemSpace=0.3*lambda; % mm
Nelem = 128; apod = rectwin(Nelem);
steerAng = 0; % degrees
focDepth = 20; % mm
alpha = 0.0; % (dB/mm/MHz) attenuation of medium

% Complex Apodization as Function of Frequency
apod_P_f = P_f' * apod';

% Simulation Space and Time
tFoc = focDepth/c; t = (0:0.05:2)*tFoc; m = 1; n = 2; 
Nx0 = 256; x = (-(Nx0-1)/2:(Nx0-1)/2)*(elemSpace/m); dov = 1; 
Nu1=round(dov*c*max(t)/(elemSpace/n)); z=((0:Nu1-1))*elemSpace/n;

dB = 60; % Display Dynamic Range (Decibels)

%% Continuous Wave Response

% Continuous Wave Response
f_idx = 100; % Index into frequency vector: Get CW Response at One Frequency
psfCW = responseCW(x, z, elemSpace, apod', steerAng, focDepth, f(f_idx), c, alpha);

% Display Continuous Wave Response at Certain Frequency
psfCW = squeeze(psfCW); psfCWMag = abs(psfCW); 
maxpsfCW = max(abs(psfCW(~isinf(psfCW) & ~isnan(psfCW))));
figure; imagesc(x,z,20*log10(psfCWMag/(maxpsfCW)),[-dB 0]);
zoom on; axis equal; axis xy; axis image;
ylabel('z Axial Distance (mm)');
xlabel('x Azimuthal Distance (mm)');
title(['Frequency = ', num2str(f(f_idx)), ' MHz']);

% K-Space of the One-Way Transmit Response
psfCW(isinf(psfCW) | isnan(psfCW)) = 0;
kspace = fftshift(fft2(real(psfCW)));
dx = mean(diff(x)); nxFFT = numel(x); % Spatial Sampling for FFT Axis
kx = (1/(2*dx))*(-1:2/nxFFT:1-2/nxFFT); % FFT Axis for Lateral Spatial Frequency
dz = mean(diff(z)); nzFFT = numel(z); % Spatial Sampling for FFT Axis
kz = (1/(2*dz))*(-1:2/nzFFT:1-2/nzFFT); % FFT Axis for Lateral Spatial Frequency
figure; imagesc(kx, kz, abs(kspace)); 
axis image; axis xy; title('K-Space of Transmit Response'); 
xlabel('lateral frequency [1/mm]'); ylabel('axial frequency [1/mm]');


%% Transmit K-Space Representation and Pulse Wave Propagation

% Continuous Wave and Pulse Wave Responses
psf_t = responsePW(x, z, elemSpace, steerAng, focDepth, f, apod_P_f, c, t, alpha);

% K-Space of the One-Way Transmit Response
psf_t(isinf(psf_t) | isnan(psf_t)) = 0;
kspace = fftshift(fft2(real(psf_t(:,:,(end+1)/2))));
dx = mean(diff(x)); nxFFT = numel(x); % Spatial Sampling for FFT Axis
kx = (1/(2*dx))*(-1:2/nxFFT:1-2/nxFFT); % FFT Axis for Lateral Spatial Frequency
dz = mean(diff(z)); nzFFT = numel(z); % Spatial Sampling for FFT Axis
kz = (1/(2*dz))*(-1:2/nzFFT:1-2/nzFFT); % FFT Axis for Lateral Spatial Frequency
figure; imagesc(kx, kz, abs(kspace)); axis image; axis xy; 
xlabel('lateral frequency [1/mm]'); ylabel('axial frequency [1/mm]');
title('K-Space of Transmit Response');

% Plotting the Result
maxpsf_t = max(abs(psf_t(~isinf(psf_t) & ~isnan(psf_t)))); 
figure; M = moviein(length(t)); kk = 1; shouldSaveMovie = true;
while(1)
    psf_tMag = abs(psf_t(:,:,kk));
    imagesc(x,z,20*log10(psf_tMag/(maxpsf_t)),[-dB 0]);
    zoom on; axis equal; axis xy; axis image;
    ylabel('z Axial Distance (mm)');
    xlabel('x Azimuthal Distance (mm)');
    M(kk) = getframe(gcf);
    if shouldSaveMovie
        [A,map] = rgb2ind(frame2im(M(kk)),256);
        if kk == 1;
            imwrite(A,map,'FastestMethod','gif', 'Loopcount',inf);
        else
            imwrite(A,map,'FastestMethod','gif','WriteMode','append');
        end
    end
    if kk == length(t)
        kk = 1; shouldSaveMovie = false;
    else 
        kk = kk + 1;
    end
end

%% Two Way (Pulse-Echo) Response and K-Space Representation

% Get the Two-Way Response (Point Spread Function)
apod_tx = sqrt(P_f')*apod'; txSteerAng = steerAng; txFocDepth = focDepth;
apod_rx = sqrt(P_f')*ones(size(apod')); rxSteerAng = steerAng; rxFocDepth = focDepth;
psf = response2WayPW(x, z, elemSpace, apod_tx, txSteerAng, ...
    txFocDepth, apod_rx, rxSteerAng, rxFocDepth, f, c, alpha);
maxpsf_t = max(abs(psf(~isinf(psf) & ~isnan(psf)))); 
figure; imagesc(x,z,20*log10(abs(psf)/(maxpsf_t)),[-dB 0]);
zoom on; axis equal; axis xy; axis image; axis xy; 
ylabel('z Axial Distance (mm)');
xlabel('x Azimuthal Distance (mm)');
title('Two-Way Response Point Spread Function');

% K-Space of the Two-Way Response
psf(isinf(psf) | isnan(psf)) = 0; kspace = fftshift(fft2(real(psf)));
dx = mean(diff(x)); nxFFT = numel(x); % Spatial Sampling for FFT Axis
kx = (1/(2*dx))*(-1:2/nxFFT:1-2/nxFFT); % FFT Axis for Lateral Spatial Frequency
dz = mean(diff(z)); nzFFT = numel(z); % Spatial Sampling for FFT Axis
kz = (1/(2*dz))*(-1:2/nzFFT:1-2/nzFFT); % FFT Axis for Lateral Spatial Frequency
figure; imagesc(kx, kz, abs(kspace)); axis image; axis xy; 
xlabel('lateral frequency [1/mm]'); ylabel('axial frequency [1/mm]');
title('K-Space of Two-Way Response');