function psf_t = responsePW(x, z, elemSpace, steerAng, focDepth, f, apod_P_f, c, t, atten)
%
% Pulse wave spatial response (on rx or tx) from xdcr array
%
% x: lateral dimension in mm
% z: axial dimension in mm
% elemSpace: element spacing in mm
% apod_P_f: apodization + pulse spectrum vs frequency (rows) and array element (columns)
% steerAng: steering angle (deg)
% focDepth: focus depth (mm)
% f: row vector array of frequencies in pulse (MHz)
% c: speed of sound (mm/usec)
% t: measurement time vector in usec (one-way)

disp('Pulsed Wave Calculations Begin Now');

% Propagation of Pulse in Time Domain
delayFactors = exp(-1i*2*pi*f'*t); % Rows = Frequency; Columns = Time
disp('Looping Over Frequencies Begins Now');

% Obtaining Continuous Wave Response
psfCWs = responseCW(x, z, elemSpace, apod_P_f, steerAng, focDepth, f, c, atten);
psfCWs_flat = reshape(psfCWs, [length(z)*length(x), length(f)]);

% Construction of Pulse-Wave Response from Continuous Wave Response
psfPWs = psfCWs_flat*delayFactors;
psf_t = reshape(psfPWs, [length(z),length(x),length(t)]);

end

