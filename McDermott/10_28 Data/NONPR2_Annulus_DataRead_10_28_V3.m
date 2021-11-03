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

filenamePrefix = "Data\";
filenameExtension = "_10-28.csv";
presVals = ["50","100","125", "150" ,"175"];
nRuns = [1,1,1,1,1];
data = struct();
start = [8800, 9200, 11584, 14240];
stop = [14500, 14000, 17941, 20679];
expandRange = 2000;

freqThreshold = 0.2;        % Include frequencies that contribute more than 10% of amplitude.
Pvary       = 0.8;          % Fraction of maximum pressure through which pressure can vary during test.

for n1 = 1:length(presVals)
    for n2 = 1:nRuns(n1)
        filename = filenamePrefix +"annulus"+presVals(n1)+"-"+ string(n2)+ filenameExtension;
        tmpData=csvread(filename,1,0);
        fieldName = "psi" + presVals(n1);
        data.(fieldName).t = tmpData(:,1);                  % s. Time
        data.(fieldName).p = tmpData(:,6)*6894.75;          % Pa. Pressure upstream of article.
        data.(fieldName).lc_1 = tmpData(:,16)*(1/2.205);    % kg. Load Cell
        data.(fieldName).lc_2 = tmpData(:,17)*(1/2.205);    % kg. Load Cell
        data.(fieldName).lc_3 = tmpData(:,18)*(1/2.205);    % kg. Load Cell
        data.(fieldName).lc_4 = tmpData(:,19)*(1/2.205);    % kg. Load Cell
        data.(fieldName).m = (data.(fieldName).lc_1 + data.(fieldName).lc_2 +data.(fieldName).lc_3 + data.(fieldName).lc_4);    % kg.
        clear tmpData;
    end
end

summary = zeros(length(presVals),4);
UPRel = zeros(length(presVals),1);
mDotRel = zeros(length(presVals),1);
CdRel = zeros(length(presVals),1);
for n1 = 1:length(presVals)
%Clean up Data
    fieldName = "psi" + presVals(n1);
        
    tmp = FilterData_V1(data.(fieldName).t,[data.(fieldName).p,data.(fieldName).m],freqThreshold);
    data.(fieldName).Psmooth = tmp(:,1);
    data.(fieldName).msmooth = tmp(:,2);

    minInd = find((data.(fieldName).Psmooth-Pvary*max(data.(fieldName).Psmooth))>0,1,'first');
    maxInd = find((data.(fieldName).Psmooth-Pvary*max(data.(fieldName).Psmooth))>0,1,'last');
    
    %Calculate values
    data.(fieldName).p_avg = mean(data.(fieldName).Psmooth(data.(fieldName).Psmooth > Pvary*max(data.(fieldName).Psmooth))); % Pa.
    [data.(fieldName).flowLine,data.(fieldName).S]= polyfit(data.(fieldName).t(minInd:maxInd),data.(fieldName).msmooth(minInd:maxInd),1);
    [~,data.(fieldName).delta] = polyval(data.(fieldName).flowLine,data.(fieldName).t(minInd:maxInd),data.(fieldName).S);
    %data.(fieldName).mDot = [data.(fieldName).flowLine(1) (data.(fieldName).msmooth(maxInd)-data.(fieldName).msmooth(minInd))/(data.(fieldName).t(maxInd)-data.(fieldName).t(minInd))*(1/2.205)]; % kg/s.
    data.(fieldName).mDot = data.(fieldName).flowLine(1);
    data.(fieldName).Cd = data.(fieldName).mDot/(A*sqrt(2*rho*data.(fieldName).p_avg));
    
    % UNCERTAINTY CALCULATION
    data.(fieldName).D = struct();
    data.(fieldName).D.P = (max(data.(fieldName).Psmooth(minInd:maxInd)) - min(data.(fieldName).Psmooth(minInd:maxInd)))/2;
    s1 = ((data.(fieldName).msmooth(end)+data.(fieldName).delta(end)) - (data.(fieldName).msmooth(1)-data.(fieldName).delta(1)))/(data.(fieldName).t(maxInd)-data.(fieldName).t(minInd));
    s2 = ((data.(fieldName).msmooth(end)-data.(fieldName).delta(end)) - (data.(fieldName).msmooth(1)+data.(fieldName).delta(1)))/(data.(fieldName).t(maxInd)-data.(fieldName).t(minInd));
    data.(fieldName).D.mDot = abs(s2-s1)/2;
    %data.(fieldName).D.mDot = abs(data.(fieldName).m(1) - (data.(fieldName).msmooth(maxInd)-data.(fieldName).msmooth(minInd))/(data.(fieldName).t(maxInd)-data.(fieldName).t(minInd)));
    % Relative Uncertainties
    data.(fieldName).D.PRel = data.(fieldName).D.P/data.(fieldName).p_avg;
    data.(fieldName).D.mDotRel = data.(fieldName).D.mDot/data.(fieldName).mDot;
    data.(fieldName).D.CdRel   = sqrt(data.(fieldName).D.mDotRel^2 + (data.(fieldName).D.PRel*0.5)^2);
    
    summary(n1,:) = [data.(fieldName).p_avg/1000,data.(fieldName).mDot,data.(fieldName).Cd,data.(fieldName).D.CdRel];
    UPRel(n1) = data.(fieldName).D.PRel;
    mDotRel(n1) = data.(fieldName).D.mDotRel;
    CdRel(n1) = data.(fieldName).D.CdRel;
    
    figure
    plot(data.(fieldName).t,data.(fieldName).msmooth)
    hold on
    plot(data.(fieldName).t(minInd:maxInd),polyval(data.(fieldName).flowLine,data.(fieldName).t(minInd:maxInd)))
    data.(fieldName).ax = gca();
    data.(fieldName).ax.FontName = 'Times New Roman';
    data.(fieldName).ax.FontSize = 18;
end

%% SINGLED-OUT PLOT FOR PROGRESS REPORT 2
choose = 2;
sample = "psi"+presVals(choose);
ni = find((data.(sample).Psmooth-Pvary*max(data.(sample).Psmooth))>0,1,'first');
xi = find((data.(sample).Psmooth-Pvary*max(data.(sample).Psmooth))>0,1,'last');
% Cleaned Signal
figure
plot(data.(sample).t-min(data.(sample).t),data.(sample).msmooth-min(data.(sample).msmooth),'LineWidth',1.5);
hold on
tsample = [data.(sample).t(ni),data.(sample).t(xi)]-min(data.(sample).t);
plot(tsample,polyval(data.(sample).flowLine,tsample+min(data.(sample).t))-min(data.(sample).msmooth),'LineWidth',1)
b = gca();
b.FontName = 'Times New Roman';
b.FontSize = 18;
ylim(b,[0,ceil((max(data.(sample).msmooth)-min(data.(sample).msmooth)))]);
ylabel(b, 'Water Mass Flowed [kg]');
xlabel(b, 'Time [s]');
legend({'Cleaned Data','Best Fit Line'},'Location','northwest')

% Original Data
figure
plot(data.(sample).t-min(data.(sample).t),data.(sample).m-data.(sample).m(1));
a = gca();
a.FontName = 'Times New Roman';
a.FontSize = 18;
ylim(a,[0,ceil((max(data.(sample).msmooth)-min(data.(sample).msmooth)))]);
ylabel(a, 'Water Mass Flowed [kg]');
xlabel(a, 'Time [s]');


%% CUSTOM FUNCTIONS
function formatPlot(ax)
ax.FontName = 'Times New Roman';
ax.FontSize = 18;
end