%{
10/24/2021
Initial Testing Analysis
%}
clc; clear; close all;
%% CONSTANTS
g   = 9.80665;  % m/s^2.
rho = 998;      % kg/m^3. Density of water at 70 deg F.
% Pintle Geometry
orif_d = 0.0625/39.37;      % m. Orifice Diameter [in]
n_orif = 20;
A = pi*(orif_d/2)^2*n_orif; % m^2. Pintle Flow area [in^2]

% Annulus Geometry
ann_d = .767/39.37;         % m
pint_out_d = 0.742/39.37;   % m

%% Read in Data
%%some intializations

filenamePrefix = "Data\Pintle\Pintle_10-22_";
filenameExtension = ".csv";
presVals = ["50","100", "125" ,"150"];
nRuns = [3,1,1,1];
data = struct();
start = {[9200, 10316, 12416], [9252], [10390], [10006] };
stop = {[15832, 15711, 14818], [15440], [18208], [13781]};

freqThreshold = 0.1;        % Include frequencies that contribute more than 10% of amplitude.
Pvary       = 0.8;          % Fraction of maximum pressure through which pressure can vary during test.

totRuns = 0;
for n1 = 1:length(presVals)
    for n2 = 1:nRuns(n1)  
    filename = filenamePrefix + presVals(n1) +"-"+ string(n2) + filenameExtension;
    tmpData=csvread(filename,1,0);
    fieldName = "psi" + presVals(n1) + "_" + string(n2);
    data.(fieldName).t = tmpData(:,1); % Time in seconds
    data.(fieldName).p = tmpData(:,7); % Pressure upstream of article in psi
    data.(fieldName).lc_1 = tmpData(:,16); % Load Cell (lbf)
    data.(fieldName).lc_2 = tmpData(:,17);
    data.(fieldName).lc_3 = tmpData(:,18);
    data.(fieldName).lc_4 = tmpData(:,19);
    data.(fieldName).m = (data.(fieldName).lc_1 + data.(fieldName).lc_2 +data.(fieldName).lc_3 + data.(fieldName).lc_4); %lbf
    clear tmpData;
    totRuns = totRuns + 1;
    end
end
 %% Data Analysis
summary = zeros(totRuns,4);
 UPRel = zeros(totRuns,1);
 mDotRel = zeros(totRuns,1);
 CdRel = zeros(totRuns,1);
 currRun = 1;
 for n1 = 1:length(presVals)
     for n2 = 1:nRuns(n1)
        %Clean up Data
        filename = filenamePrefix + presVals(n1) +"-"+ string(n2) + filenameExtension;
        fieldName = "psi" + presVals(n1) + "_" + string(n2);
        
        tmp = FilterData_V1(data.(fieldName).t,[data.(fieldName).p,data.(fieldName).m],freqThreshold);
        data.(fieldName).Psmooth = tmp(:,1);
        data.(fieldName).msmooth = tmp(:,2);
        
        minInd = find((data.(fieldName).Psmooth-Pvary*max(data.(fieldName).Psmooth))>0,1,'first');
        maxInd = find((data.(fieldName).Psmooth-Pvary*max(data.(fieldName).Psmooth))>0,1,'last');
        %Calculate values
        data.(fieldName).p_avg = mean(data.(fieldName).Psmooth(data.(fieldName).Psmooth > Pvary*max(data.(fieldName).Psmooth)))*6894.75; % Pa.
        data.(fieldName).mDot = (data.(fieldName).msmooth(maxInd)-data.(fieldName).msmooth(minInd))/(data.(fieldName).t(maxInd)-data.(fieldName).t(minInd))*(1/2.205); % kg/s.
        data.(fieldName).Cd = data.(fieldName).mDot/(A*sqrt(2*rho*data.(fieldName).p_avg));
        
        % UNCERTAINTY CALCULATION
        data.(fieldName).D = struct();
        data.(fieldName).D.P = (max(data.(fieldName).Psmooth(minInd:maxInd)) - min(data.(fieldName).Psmooth(minInd:maxInd)))*6894.75/2;
        flowLine = polyfit(data.(fieldName).t(minInd:maxInd),data.(fieldName).msmooth(minInd:maxInd)*(1/2.205),1);
        data.(fieldName).D.mDot = abs(flowLine(1) - data.(fieldName).mDot);
        % Relative Uncertainties
        data.(fieldName).D.PRel = data.(fieldName).D.P/data.(fieldName).p_avg;
        data.(fieldName).D.mDotRel = data.(fieldName).D.mDot/data.(fieldName).mDot;
        data.(fieldName).D.CdRel   = sqrt(data.(fieldName).D.mDotRel^2 + (data.(fieldName).D.PRel*0.5)^2);
        
        %Plots
        figure
        plot(data.(fieldName).t(minInd:maxInd),data.(fieldName).Psmooth(minInd:maxInd))
        hold on
        plot(data.(fieldName).t(minInd:maxInd),data.(fieldName).msmooth(minInd:maxInd))
        title(fieldName)
        
        summary(currRun,:) = [data.(fieldName).p_avg/1000,data.(fieldName).mDot,data.(fieldName).Cd,data.(fieldName).D.CdRel];
        UPRel(currRun) = data.(fieldName).D.PRel;
        mDotRel(currRun) = data.(fieldName).D.mDotRel;
        CdRel(currRun) = data.(fieldName).D.CdRel;
        currRun = currRun + 1;
     end
 end
