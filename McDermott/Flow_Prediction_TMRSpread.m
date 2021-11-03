clear;
clc;

inchPerMeter = 39.37;
PaPerPsi = 6894.75729;

rho_air     = 1.184;

%% INPUTS
% Experimental Boundaries
TMR     = (0.5:0.5:2.0)';
binRange = 30;

Pmax    = 200*PaPerPsi;

massSensitivity     = 0.020;    % kg. Lowest mass that can be measured.
deltaMass           = 0.05;     % [-]. Fraction of the mass flow rate desired to be resolved.
angSector           = 1/10;     % Angle sweep needed for a symmetric slice of the pintle (in the phi direction)

Vavailable      = 0.20;         % m^3. Approximately 55 gallons of water.

load('propellants2238.mat')
load('geom2238.mat')

% Orifice Properties
d   = geom.d_h;
L   = geom.tau_p;
r   = d/10;

% Annulus Properties
d_ann   = 0.0125/inchPerMeter;
L_ann   = 0.25/inchPerMeter;
r_ann   = 0.1/inchPerMeter;

% Water Flow Fluid Properties
pintle = LOX;
annulus = RP1;
properties = fieldnames(LOX);
for n1 = 1:length(properties)
    if ~strcmp(properties{n1},'A')
        pintle.(properties{n1}) = [];
        annulus.(properties{n1}) = [];
    end
end
% Fluid Properties
pintle.rho     = 998;       % Liquid density. kg m^-3
pintle.mu      = 1.002e-3;  % Dynamic Viscosity, kg m^-1 s^-1
pintle.Pvap    = 2339;      % Pa. Vapor pressure, 300K
pintle.sigma    = 0.0717;   % N/m. Surface tension
pintle.Cd       = 0.59;

maxPintle = struct();
maxPintle.mDot  = pintle.Cd*pintle.A*sqrt(2*pintle.rho*Pmax);
maxPintle.v     = maxPintle.mDot/(pintle.A*pintle.rho);

annulus.rho     = 998;       % Liquid density. kg m^-3
annulus.mu      = 1.002e-3;  % Dynamic Viscosity, kg m^-1 s^-1
annulus.Pvap    = 2339;      % Pa. Vapor pressure, 300K
annulus.sigma    = 0.0717;
annulus.Cd      = 0.56;

maxAnnulus = struct();
maxAnnulus.mDot = annulus.Cd*annulus.A*sqrt(2*annulus.rho*Pmax);
maxAnnulus.v    = maxAnnulus.mDot/(annulus.A*annulus.rho);


%% SETUP
% Pressure and Mass Flow Properties
varyPintle      = TMR < ((maxPintle.mDot*maxPintle.v)/(maxAnnulus.mDot*maxAnnulus.v));
varyAnnulus     = ~varyPintle;

% Limited by pintle pressure - maximize annulus flow.
annulus.mDot    = repmat(maxAnnulus.mDot,sum(varyPintle),1);
pintle.mDot     = sqrt(TMR(varyPintle).*(pintle.rho*pintle.A)/(annulus.rho*annulus.A)).*annulus.mDot(varyPintle);

% Limited by annulus pressure - maximize pintle flow.
pintle.mDot     = [pintle.mDot; repmat(maxPintle.mDot,sum(varyAnnulus),1)];
annulus.mDot    = [annulus.mDot; sqrt(1./TMR(varyAnnulus).*(annulus.rho*annulus.A)/(pintle.rho*pintle.A)).*pintle.mDot(varyAnnulus)];

annulus.v   = annulus.mDot/(annulus.rho*annulus.A);
pintle.v    = pintle.mDot/(pintle.rho*pintle.A);

% Approximate pressure drops used for Cd calculation. Actual pressure drop
% calculated after.
dPAnn   = (1/2)*annulus.rho*((annulus.v/annulus.Cd).^2);
dPPintle= (1/2)*pintle.rho*((pintle.v/pintle.Cd).^2);
thetaPred = @(T) 0.6243*atan(3.377*T);
angles =thetaPred(TMR)*(180/pi);
summary = [TMR,dPAnn/PaPerPsi,dPPintle/PaPerPsi,angles,angles-binRange/2,angles+binRange/2];
