%{
10/24/2021
Initial Testing Analysis
%}
clc; clear; close all;
%% Read in Data
%%some intializations
%filenamePrefix = 
%filenameExtension = 

filename = 'Data\RPLRadial502_sensors_21-10-22_1824.csv';
tmpData=csvread(filename,1,0);

t =tmpData(:,1); % Time in seconds
p = tmpData(:,7); % Pressure upstream of article in psi
lc_1=tmpData(:,16); % Load Cell (lbf)
lc_2=tmpData(:,17);
lc_3=tmpData(:,18);
lc_4=tmpData(:,19);
clear tmpData;
%% Clean up Data
start = find(p,1,'first');
stop = find(p,1,'last');
start = 10000;
stop = 14000;
t = t(start:stop);
p = p(start:stop);
lc_1 = lc_1(start:stop);
lc_2 = lc_2(start:stop);
lc_3 = lc_3(start:stop);
lc_4 = lc_4(start:stop);
m = (lc_1 + lc_2 +lc_3 + lc_4); %lbf

%% Independent Variables
rho = 998; % kg/m^3. Density of water at 70 deg F.
% Pintle Geometry
orif_d = 0.0625/39.37; % m. Orifice Diameter [in]
n_orif = 20;
A = pi*(orif_d/2)^2*n_orif; % m^2. Pintle Flow area [in^2]

%% Calculate values
p_avg_ovr = mean(p)*6894.75; % Pa.
mFlowRate_ovr = (m(end)-m(1))/(t(end)-t(1))*(1/2.205); % kg/s.
Cd_ovr = mFlowRate_ovr/(A*sqrt(2*rho*p_avg_ovr));


%% Plots
figure
plot(t,p)
hold on
plot(t,m)

