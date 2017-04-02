function psf = responseAngSpec2WayPW(x, z, elemSpace, apod_tx_P_f, ...
    txSteerAng, txFocDepth, apod_rx_P_f, rxSteerAng, rxFocDepth, f, c)
%
% Pulse wave spatial response (on rx or tx) from xdcr array
%
% x: lateral dimension in mm
% z: axial dimension in mm
% elemSpace: element spacing in mm
% apod: apodization vector
% steerAng: steering angle (deg)
% focDepth: focus depth (mm)
% P_f: pulse as a function of frequency 
% f: row vector array of frequencies in pulse (MHz)
% c: speed of sound (mm/usec)
% t: measurement time vector in usec (one-way)

disp('Pulsed Wave Calculations Begin Now');

% X-Coordinate of Apertures
Nelem = size(apod_tx_P_f, 2); 
xpos = (-(Nelem-1)/2:(Nelem-1)/2)*elemSpace; 

% Propagation of Pulse in Time Domain
t = 2*rxFocDepth/c;
delayFactors = exp(-1i*2*pi*f'*t); % Rows = Frequency; Columns = Time

disp('Looping Over Frequencies Begins Now');

% Obtaining Continuous Wave Response
psfCWsTx = responseAngSpecCW(x, z, elemSpace, apod_tx_P_f, txSteerAng, txFocDepth, f, c);
psfCWsRx = responseAngSpecCW(x, z, elemSpace, apod_rx_P_f, rxSteerAng, rxFocDepth, f, c);
psfCWs_flat = reshape(psfCWsTx .* psfCWsRx, [length(z)*length(x), length(f)]);

% Construction of Pulse-Wave Response from Continuous Wave Response
psfPWs = psfCWs_flat*delayFactors;
psf = reshape(psfPWs, [length(z),length(x),length(t)]);

end

