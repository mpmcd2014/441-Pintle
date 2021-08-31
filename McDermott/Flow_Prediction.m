clear;
clc;

TMR = 0.2:0.2:4;
Pmax    = 
%% SETUP
inchPerMeter = 39.37;
PaPerPsi = 6894.75729;

load('propellants1869.mat')
load('geom1869.mat')

% Orifice Properties
d   = geom.d_h;
L   = geom.tau_p;
r   = d/10;

% Flow Environment
deltaP = 110*PaPerPsi;          % Approximate pressure drop across pintle holes
Pchamber = 400*PaPerPsi;        % Chamber pressure

%% LOX: RESULTS 1869, SPREAD3A_F.txt
% Fluid Properties
LOX.mu  = 1.526e-4;                 % Dynamic Viscosity, kg m^-1 s^-1
LOX.Pvap = 37*PaPerPsi;             % Vapor pressure (LOx at 100 K)
LOX.sigma = 0.010762;               % N/m. Surface tension.

rho_chamber     = 4;
% CHAMBER PRESSURE STATE
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
%% H2O - PINTLE
pintle = LOX;
properties = fieldnames(LOX);
for n1 = 1:length(properties)
    pintle.(properties{n1}) = [];
end
% Fluid Properties
pintle.rho     = 998;                      % Liquid density. kg m^-3
pintle.mu      = 1.002e-3;                  % Dynamic Viscosity, kg m^-1 s^-1
pintle.Pvap    = 2339;                    % Vapor pressure
pintle.sigma    = 0.0717;

rho_air     = 1.184;

% WATER FLOW PRESSURE STATE
P1  = 101325+deltaP;
P2  = 101325;
Re  = (pintle.rho*pintle.v*L)/pintle.mu;


pintle.v       = LOX.v*sqrt(LOX.rho/rho);

%DISCHARGE COEFFICIENT
    % Nurick et al.:
    pintle.Cd  = ((1/(.868-.0425*sqrt(L/d)))+((20/Re)*(1+2.25*(L/d)))-((.005*(L/d))/(1+7.5*((log(.0015*Re))^2))))^-1;
    % ANSYS Fluent Method
    [pintle.Cd(2),K_H2O] = orificeCavitation(pintle.rho,pintle.mu, d, L, r, P1, P2,pintle.Pvap,0);
    
% DIAMETER
D_H2O   = dropDiameter(pintle.d,pintle.v,pintle.rho,pintle.sigma);

% BREAKUP DISTANCE
We_H2O  = pintle.rho*(pintle.v^2)*d/pintle.sigma;
rb_H2O  = 0.5*d*14.2*((rho_air/pintle.rho)^(-2/3))*(We_H2O^(-1/3));