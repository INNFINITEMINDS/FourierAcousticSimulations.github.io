function psf = response2WayPW(x, z, elemSpace, apod_tx_P_f, txSteerAng, ...
    txFocDepth, apod_rx_P_f, rxSteerAng, rxFocDepth, f, c, atten)
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

disp('Pulsed Wave Calculations Begin Now');

% X-Coordinate of Apertures
Nelem = size(apod_tx_P_f, 2); 
xpos = (-(Nelem-1)/2:(Nelem-1)/2)*elemSpace; 

% Propagation of Pulse in Time Domain
t = 2*rxFocDepth/c;
delayFactors = exp(-1i*2*pi*f'*t); % Rows = Frequency; Columns = Time

% Obtaining Continuous Wave Response
psfCWsTx = responseCW(x, z, elemSpace, apod_tx_P_f, txSteerAng, txFocDepth, f, c, atten);
psfCWsRx = responseCW(x, z, elemSpace, apod_rx_P_f, rxSteerAng, rxFocDepth, f, c, atten);
psfCWs_flat = reshape(psfCWsTx .* psfCWsRx, [length(z)*length(x), length(f)]);
psfPWs = psfCWs_flat*delayFactors;
psf = reshape(psfPWs, [length(z),length(x)]);

% % Obtaining Continuous Wave Response
% psf = zeros(length(z), length(x));
% disp('Looping Over Frequencies Begins Now');
% for f_idx = 1:numel(f)
%     psfCWTx = responseCW(x, z, elemSpace, apod_tx_P_f(f_idx,:), ...
%         txSteerAng, txFocDepth, f(f_idx), c, atten);
%     psfCWRx = responseCW(x, z, elemSpace, apod_rx_P_f(f_idx,:), ...
%         rxSteerAng, rxFocDepth, f(f_idx), c, atten);
%     psf = psf + delayFactors(f_idx)*(psfCWTx .* psfCWRx);
%     disp(['f = ' num2str(f(f_idx)), ' MHz Completed']);
% end

end

