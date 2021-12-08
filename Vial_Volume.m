%% Title

% Vial Volume Calculator
% For use on AME 441 Pintle Project
% Sean Cornish
% 9/10/2021

clc; clear; close all;

%% Variables

%Pintle Diameter [m]
dP = .0188;

%Water Density [kg/m^3]
rhoW = 997;

%mass flow rate [kg/s]
mDot = 2.06;

%Test runtime [s]
t = 30;

% Spray Cone Angle (fron Horizontal) [deg]
theta = 79;

%Collection radius [m]
%distance from collision to collection
rC = .07;

%collection vial Diameter [m]
dC = .003;

%% Calcluations

%functional spray radius [m]
rFS = rC*cosd(theta);

% Functional Total Radius [m]
rFT = dP/2 + rFS;

%gamma [deg]
% slice of perverbial pie
gamma = 2*asind((dC/2)/rFT);

%Total Volume Flow [m^3]
vT = mDot*t/rhoW;

%volume in vial [m^3]
vV = (gamma/360)*vT;

%Mililiters [mL]
vML = vV*1000000;

fprintf("the maximum volume of water in a vial for this test would be %f mL",vML);



