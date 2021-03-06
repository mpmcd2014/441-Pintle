clear;
clc;

inchPerMeter = 39.37;
PaPerPsi = 6894.75729;

rho_air     = 1.184;

%% INPUTS
% Experimental Boundaries
TMR     = (0.2:0.2:4)';
NTest   = 10;
massFlowFrac = 1;
reuseFrac = 0.33;

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

annulus.rho     = 998;       % Liquid density. kg m^-3
annulus.mu      = 1.002e-3;  % Dynamic Viscosity, kg m^-1 s^-1
annulus.Pvap    = 2339;      % Pa. Vapor pressure, 300K
annulus.sigma    = 0.0717;


%% SETUP
% Pressure and Mass Flow Properties
varyPintle      = TMR < ((LOX.mDot*LOX.v)/(RP1.mDot*RP1.v));
varyAnnulus     = ~varyPintle;

% Limited by pintle pressure - maximize annulus flow.
annulus.mDot    = repmat(sqrt(annulus.rho*annulus.A*massFlowFrac*RP1.mDot*RP1.v),sum(varyPintle),1);
pintle.mDot     = sqrt(TMR(varyPintle).*(pintle.rho*pintle.A)/(annulus.rho*annulus.A)).*annulus.mDot(varyPintle);

% Limited by annulus pressure - maximize pintle flow.
pintle.mDot     = [pintle.mDot; repmat(sqrt(pintle.rho*pintle.A*massFlowFrac*LOX.mDot*LOX.v),sum(varyAnnulus),1)];
annulus.mDot    = [annulus.mDot; sqrt(1./TMR(varyAnnulus).*(annulus.rho*annulus.A)/(pintle.rho*pintle.A)).*pintle.mDot(varyAnnulus)];

annulus.v   = annulus.mDot/(annulus.rho*annulus.A);
pintle.v    = pintle.mDot/(pintle.rho*pintle.A);

% Approximate pressure drops used for Cd calculation. Actual pressure drop
% calculated after.
dPAnn   = (1/2)*annulus.rho*((annulus.v/0.8).^2);
dPPintle= (1/2)*pintle.rho*((pintle.v/0.8).^2);
Pa      = 101325;
%% LOX: RESULTS 1869, SPREAD3A_F.txt

% Fluid Properties
LOX.mu  = 1.526e-4;                 % Dynamic Viscosity, kg m^-1 s^-1
LOX.Pvap = 37*PaPerPsi;             % Vapor pressure (LOx at 100 K)
LOX.sigma = 0.010762;               % N/m. Surface tension.

rho_chamber     = 4;
% CHAMBER PRESSURE STATE
deltaP = 110*PaPerPsi;          % Approximate pressure drop across pintle holes
Pchamber = 400*PaPerPsi;        % Chamber pressure

P1 = deltaP + Pchamber;         % Upstream pressure
P2 = Pchamber;

Re  = (LOX.rho*LOX.v*L)/LOX.mu;

% DISCHARGE COEFFICIENT
    % Nurick et al.:
    Cd_LOX  = ((1/(.868-.0425*sqrt(L/d)))+((20/Re)*(1+2.25*(L/d)))-((.005*(L/d))/(1+7.5*((log(.0015*Re))^2))))^-1;
    % ANSYS Fluent Method
    [Cd_LOX(2),K_LOX] = orificeCavitation(LOX.rho, LOX.mu, d, L, r, P1, P2, LOX.Pvap,0);

% DIAMETER
D_LOX = dropDiameter(d,LOX.v,LOX.rho,LOX.sigma);

% BREAKUP DISTANCE
We_LOX  = LOX.rho*(LOX.v^2)*d/LOX.sigma;
rb_LOX  = 0.5*d*14.2*((rho_chamber/LOX.rho)^(-2/3))*(We_LOX^(-1/3));

%% RP1: RESULTS
RP1.mu  = 0.0024;
RP1.Pvap = 0;
RP1.sigma = 0.0235;

Re  = (RP1.rho*RP1.v*L)/RP1.mu;

% DISCHARGE COEFFICIENT
    % Nurick et al.:
    Cd_RP1  = ((1/(.868-.0425*sqrt(L/d)))+((20/Re)*(1+2.25*(L/d)))-((.005*(L/d))/(1+7.5*((log(.0015*Re))^2))))^-1;
    % ANSYS Fluent Method
    [Cd_RP1(2),K_RP1] = orificeCavitation(RP1.rho, RP1.mu, d_ann, L_ann, r_ann, P1, P2, RP1.Pvap,0);

% DIAMETER
D_RP1 = dropDiameter(d,RP1.v,RP1.rho,RP1.sigma);

% BREAKUP DISTANCE
We_RP1  = RP1.rho*(RP1.v^2)*d/RP1.sigma;
rb_RP1  = 0.5*d*14.2*((rho_chamber/RP1.rho)^(-2/3))*(We_RP1^(-1/3));
%% H2O - PINTLE
pintle.Re  = (pintle.rho*pintle.v*L)/pintle.mu;

%DISCHARGE COEFFICIENT
    % Nurick et al.:
    pintle.Cd  = ((1/(.868-.0425*sqrt(L/d)))+((20./pintle.Re)*(1+2.25*(L/d)))-((.005*(L/d))./(1+7.5*((log(.0015*pintle.Re)).^2)))).^-1;
    % ANSYS Fluent Method
    [pintle.Cd(:,2),K_H2O] = orificeCavitation(pintle.rho,pintle.mu, d, L, r, dPPintle+Pa, Pa,pintle.Pvap,0);
    
pintle.deltaP   = (1/2)*pintle.rho*((pintle.v./pintle.Cd).^2);
% DIAMETER
pintle.D_H2O   = dropDiameter(d,pintle.v,pintle.rho,pintle.sigma);

% BREAKUP DISTANCE
pintle.We  = pintle.rho*(pintle.v.^2)*d/pintle.sigma;
pintle.rb  = 0.5*d*14.2*((rho_air/pintle.rho)^(-2/3))*(pintle.We.^(-1/3));

%% H2O - ANNULUS
annulus.Re  = (annulus.rho*annulus.v*L)/annulus.mu;

%DISCHARGE COEFFICIENT
    % Nurick et al.:
    annulus.Cd  = ((1/(.868-.0425*sqrt(L/d)))+((20./annulus.Re)*(1+2.25*(L/d)))-((.005*(L/d))./(1+7.5*((log(.0015*annulus.Re)).^2)))).^-1;
    % ANSYS Fluent Method
    [annulus.Cd(:,2),K_H2O] = orificeCavitation(annulus.rho,annulus.mu, d_ann, L_ann, r_ann, dPPintle+Pa, Pa,annulus.Pvap,0);
    
annulus.deltaP   = (1/2)*annulus.rho*((annulus.v./annulus.Cd).^2);
% DIAMETER
annulus.D_H2O   = dropDiameter(d,annulus.v,annulus.rho,annulus.sigma);

% BREAKUP DISTANCE
annulus.We  = annulus.rho*(annulus.v.^2)*d/annulus.sigma;
annulus.rb  = 0.5*d*14.2*((rho_air/annulus.rho)^(-2/3))*(annulus.We.^(-1/3));

%% INSTRUMENTATION REQUIREMENTS
mDotTotal = pintle.mDot + annulus.mDot;

tMin    = massSensitivity*(1/deltaMass)*(1/angSector)*(1/min(mDotTotal));

N   = (annulus.rho*Vavailable)/mean(mDotTotal)/tMin;

V_NOREUSE       = 3*sum(interp1(TMR,mDotTotal,(linspace(TMR(1),TMR(end),NTest))))*tMin/annulus.rho;
V_REUSE         = (1-reuseFrac)*3*sum(interp1(TMR,mDotTotal,(linspace(TMR(1),TMR(end),NTest))))*tMin/annulus.rho;
V_EXTRAREUSE    = (1-reuseFrac)*5*sum(interp1(TMR,mDotTotal,(linspace(TMR(1),TMR(end),NTest))))*tMin/annulus.rho;


%% PLOTTING
close all
figure()
plot(TMR,pintle.We,TMR,annulus.We);
plotformat(gca());
xlabel('TMR')
ylabel('Weber Number')
figure()
plot(TMR,pintle.rb*1000,TMR,annulus.rb*1000);
xlabel('TMR')
ylabel('Radius at Breakup [mm]');
plotformat(gca());
figure()
plot(TMR,pintle.D_H2O*1e6,TMR,annulus.D_H2O*1e6);
plotformat(gca());

function plotformat(a)
a.FontName = 'Times New Roman';
a.FontSize = 12;
end
