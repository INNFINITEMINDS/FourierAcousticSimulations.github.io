function psf = responseCW(x, z, elemSpace, apod_P_f, steerAng, focDepth, f, c, atten)
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

% Take Care of Units for Attenuation
att = atten/(20*log10(exp(1))); % dB to Neper

% X-Coordinate of Apertures
Nelem = size(apod_P_f, 2); 
xpos = (-(Nelem-1)/2:(Nelem-1)/2)*elemSpace; 

% Simulation Space Grid
[X, Z] = meshgrid(x, z); 
X = X(:); Z = Z(:); 

% Defining Green's Function for Impulse response:
R = @(X, Z) sqrt(X.^2 + Z.^2);
G = @(X, Z) exp(-att*R(X,Z)*f) .* ...
    exp(1i*2*pi*R(X,Z)*f/c) ./ sqrt(R(X,Z)*ones(size(f)));
% The Green's Function below is exact but takes more time to compute
% G = @(X, Z) exp(-att*R(X,Z)*f) .* besselh(0,1,2*pi*R(X,Z)*f/c);

% Calculating Time Delays
delayIdeal=(sqrt(xpos.^2+focDepth^2-2*focDepth*xpos*sind(steerAng))-focDepth)/c;

% Aperture Function
a = apod_P_f'.*exp(-1i*2*pi*delayIdeal'*f);

% Assembly of Continuous Wave Response
psf = zeros([numel(X), numel(f)]);
for kk = 1:numel(xpos)
    psf = psf + G(X-xpos(kk), Z)*diag(a(kk,:));
    disp(['Element ' num2str(kk) ' Completed'])
end
psf = reshape(psf, [numel(z), numel(x), numel(f)]);

toc; disp('CW Response Calculation END');

end
