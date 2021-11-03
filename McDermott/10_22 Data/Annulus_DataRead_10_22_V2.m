%{
10/24/2021
Initial Testing Analysis
%}
clc; clear; close all;
%% CONSTANTS
g   = 9.80665;  % m/s^2.
rho = 998; % kg/m^3. Density of water at 70 deg F.
% Pintle Geometry
orif_d = 0.0625/39.37; % m. Orifice Diameter [in]
n_orif = 20;
% A = pi*(orif_d/2)^2*n_orif; % m^2. Pintle Flow area [in^2]

% Annulus Geometry
ann_d = .767/39.37; % m
pint_out_d = 0.742/39.37; % m
A = pi*(ann_d/2)^2 - pi*(pint_out_d/2)^2; % m^2

%% Read in Data
%%some intializations

filenamePrefix = "Data\Annulus\Annulus_10-22_";
filenameExtension = ".csv";
presVals = ["100","125", "150" ,"175"];
data = struct();
start = [8800, 9200, 11584, 14240];
stop = [14500, 14000, 17941, 20679];
expandRange = 2000;

freqThreshold = 0.1;        % Include frequencies that contribute more than 10% of amplitude.
Pvary       = 0.8;          % Fraction of maximum pressure through which pressure can vary during test.

for n1 = 1:length(presVals)
    filename = filenamePrefix + presVals(n1) + filenameExtension;
    tmpData=csvread(filename,1,0);
    fieldName = "psi" + presVals(n1);
    data.(fieldName).t = tmpData(:,1); % Time in seconds
    data.(fieldName).p = tmpData(:,7); % Pressure upstream of article in psi
    data.(fieldName).lc_1 = tmpData(:,16); % Load Cell (lbf)
    data.(fieldName).lc_2 = tmpData(:,17);
    data.(fieldName).lc_3 = tmpData(:,18);
    data.(fieldName).lc_4 = tmpData(:,19);
    data.(fieldName).m = (data.(fieldName).lc_1 + data.(fieldName).lc_2 +data.(fieldName).lc_3 + data.(fieldName).lc_4); %lbf
    clear tmpData;
end

summary = zeros(length(presVals),4);
UPRel = zeros(length(presVals),1);
mDotRel = zeros(length(presVals),1);
CdRel = zeros(length(presVals),1);
for n1 = 1:length(presVals)
%Clean up Data
    fieldName = "psi" + presVals(n1);
    
    data.(fieldName).t = data.(fieldName).t(start(n1)-2*expandRange:stop(n1)+expandRange);
    data.(fieldName).p = data.(fieldName).p(start(n1)-2*expandRange:stop(n1)+expandRange);
    data.(fieldName).m = data.(fieldName).m(start(n1)-2*expandRange:stop(n1)+expandRange);
    
    tmp = FilterData_V1(data.(fieldName).t,[data.(fieldName).p,data.(fieldName).m],freqThreshold);
    data.(fieldName).Psmooth = tmp(:,1);
    data.(fieldName).msmooth = tmp(:,2);

    minInd = find((data.(fieldName).Psmooth-Pvary*max(data.(fieldName).Psmooth))>0,1,'first');
    maxInd = find((data.(fieldName).Psmooth-Pvary*max(data.(fieldName).Psmooth))>0,1,'last');
    
    %Calculate values
    data.(fieldName).p_avg = mean(data.(fieldName).Psmooth(data.(fieldName).Psmooth > Pvary*max(data.(fieldName).Psmooth)))*6894.75; % Pa.
    flowLine = polyfit(data.(fieldName).t(minInd:maxInd),data.(fieldName).msmooth(minInd:maxInd)*(1/2.205),1);
    data.(fieldName).mDot = [flowLine(1) (data.(fieldName).msmooth(maxInd)-data.(fieldName).msmooth(minInd))/(data.(fieldName).t(maxInd)-data.(fieldName).t(minInd))*(1/2.205)]; % kg/s.
    data.(fieldName).Cd = data.(fieldName).mDot/(A*sqrt(2*rho*data.(fieldName).p_avg));
    
    % UNCERTAINTY CALCULATION
    data.(fieldName).D = struct();
    data.(fieldName).D.P = (max(data.(fieldName).Psmooth(minInd:maxInd)) - min(data.(fieldName).Psmooth(minInd:maxInd)))*6894.75/2;
    data.(fieldName).D.mDot = abs(flowLine(1) - (data.(fieldName).msmooth(maxInd)-data.(fieldName).msmooth(minInd))/(data.(fieldName).t(maxInd)-data.(fieldName).t(minInd))*(1/2.205));
    % Relative Uncertainties
    data.(fieldName).D.PRel = data.(fieldName).D.P/data.(fieldName).p_avg;
    data.(fieldName).D.mDotRel = data.(fieldName).D.mDot/data.(fieldName).mDot;
    data.(fieldName).D.CdRel   = sqrt(data.(fieldName).D.mDotRel^2 + (data.(fieldName).D.PRel*0.5)^2);
    
    summary(n1,:) = [data.(fieldName).p_avg/1000,data.(fieldName).mDot(1),data.(fieldName).Cd,data.(fieldName).D.CdRel];
    UPRel(n1) = data.(fieldName).D.PRel;
    mDotRel(n1) = data.(fieldName).D.mDotRel;
    CdRel(n1) = data.(fieldName).D.CdRel;
    
    figure
    plot(data.(fieldName).t,data.(fieldName).msmooth)
end



%% CUSTOM FUNCTIONS
function formatPlot(ax)
ax.FontName = 'Times New Roman';
ax.FontSize = 18;
end