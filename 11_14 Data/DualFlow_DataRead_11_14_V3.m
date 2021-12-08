%{
11/15/2021
Initial Testing Analysis
%}
clc; clear; close all;
%% CONSTANTS
g   = 9.80665;  % m/s^2.
rho = 998; % kg/m^3. Density of water at 70 deg F.
% Pintle Geometry
CdP     = 0.60;
orif_d = 0.0625/39.37; % m. Orifice Diameter [in]
n_orif = 20;
AP = pi*(orif_d/2)^2*n_orif; % m^2. Pintle Flow area [in^2]

% Annulus Geometry
CdAnn   = 0.57;
ann_d = .766/39.37; % m
pint_out_d = 0.7415/39.37; % m
AAnn = pi*(ann_d/2)^2 - pi*(pint_out_d/2)^2; % m^2

% ============== 1 CORRESPONDS TO RADIAL/PINTLE ===============
% ============== 2 CORRESPONDS TO AXIAL/ANNULUS ===============

Cd  = [CdP,CdAnn];
A   = [AP,AAnn];
BF  = n_orif*orif_d/(pi*pint_out_d);

% Uncertainty Parameters
t_est   = [1.363,1.383];        % 90% T-estimator for 11 pintle samples and 9 annulus samples.
DCd     = [0.0282,0.0226].*t_est;
dlin    = 0.0001/39.37;   % m. Uncertianty in pintle/annulus dimension measurement.
DA      = [sqrt(2)*(0.001/0.0625),sqrt(((2*dlin*(ann_d/2))/((ann_d/2)^2-(pint_out_d/2)^2))^2 + ((2*dlin*(pint_out_d/2))/((ann_d/2)^2-(pint_out_d/2)^2))^2)];
Drho    = 0.001;

thetaPred = @(T) 0.6243*atan(3.377*T);

printResults = {'TMR','Theta','mTotCalc','mTotMeas'};

%% Read in Data
%%some intializations

filenamePrefix = "Data\";
filenameExtension = "_21-11-14";
testNums = 3:7;
data = struct();
expandRange = 2000;

freqThreshold = 0.1;        % Include frequencies that contribute more than 10% of amplitude.
Pvary       = 0.9;          % Fraction of maximum pressure through which pressure can vary during test.
nTMR      = size(testNums,1);

testNames = {};
for n1 = 1:length(testNums)
    filename = filenamePrefix +"Test"+string(testNums(n1))+filenameExtension;
    tmpData=csvread(filename,23,0);
    fieldName = "Test" + string(testNums(n1));
    data.(fieldName).t = tmpData(:,1);                  % s. Time
    data.(fieldName).res    = repmat(struct('P',[],'Psmooth',[],'minInd',[],'maxInd',[]),2,1);
    data.(fieldName).res(2).Pt = tmpData(:,4)*6894.75;          % psi. Annulus tank pressure.
    data.(fieldName).res(1).Pt = tmpData(:,5)*6894.75;          % psi. Pintle tank pressure.
    data.(fieldName).res(2).P = tmpData(:,6)*6894.75;           % Pa. Pressure upstream of annulus.
    data.(fieldName).res(1).P = tmpData(:,7)*6894.75;          % Pa. Pressure upstream of pintle.
    if sum(data.(fieldName).res(1).P) == 0 || sum(data.(fieldName).res(2).P) == 0
        warning(fieldName+" contains no pressure data.")
    else
        testNames(end+1) = {fieldName};
    end
    data.(fieldName).m = sum(tmpData(:,16:19),2)*(1/2.205); % kg. Mass collected by all four load cells. (data.(fieldName).lc_1 + data.(fieldName).lc_2 +data.(fieldName).lc_3 + data.(fieldName).lc_4);    % kg.
    clear tmpData;
end

%% Process Data
summary = cell(nTMR+1,4);
summary(1,:) = {'tTest','TMR','Pt','P'};
UPRel = zeros(nTMR,1);
mDotRel = zeros(nTMR,1);
CdRel = zeros(nTMR,1);

mTestComp = zeros(length(testNames),3);
for n1 = 1:length(testNames)
%Clean up Data
    fieldName = testNames{n1};
    
    
    for n2 = 1:2
        tmp = FilterData_V1(data.(fieldName).t,[data.(fieldName).res(n2).P,data.(fieldName).m],freqThreshold);
        data.(fieldName).res(n2).Psmooth = tmp(:,1);
        data.(fieldName).msmooth = tmp(:,2);
        data.(fieldName).res(n2).mDot = Cd(n2)*A(n2)*sqrt(2*rho*data.(fieldName).res(n2).P);
        
        data.(fieldName).res(n2).minInd = find((data.(fieldName).res(n2).Psmooth - Pvary*max(data.(fieldName).res(n2).Psmooth))>0,1,'first');
        data.(fieldName).res(n2).maxInd = find((data.(fieldName).res(n2).Psmooth - Pvary*max(data.(fieldName).res(n2).Psmooth))>0,1,'last');
        tmpRange = (data.(fieldName).res(n2).minInd:data.(fieldName).res(n2).maxInd);
        
        data.(fieldName).res(n2).Pavg = mean(data.(fieldName).res(n2).Psmooth(tmpRange));
        data.(fieldName).res(n2).mDotAvg = mean(data.(fieldName).res(n2).mDot(tmpRange));
    end
    dualFlow    = intersect(data.(fieldName).res(1).minInd:data.(fieldName).res(1).maxInd,data.(fieldName).res(2).minInd:data.(fieldName).res(2).maxInd);
    if isempty(dualFlow)
        warning("No intersecting steady data for "+fieldName)
    else
        figure
        plot(data.(fieldName).t(dualFlow),data.(fieldName).res(1).Psmooth(dualFlow)/6894.75);
        hold on
        plot(data.(fieldName).t(dualFlow),data.(fieldName).res(2).Psmooth(dualFlow)/6894.75);
        title(fieldName)

        for n2 = 1:2
                    data.(fieldName).res(n2).PtSummary = [mean(data.(fieldName).res(n2).Pt(tmpRange)),min(data.(fieldName).res(n2).Pt(tmpRange)),max(data.(fieldName).res(n2).Pt(tmpRange))];
        end
        data.(fieldName).tTest  = data.(fieldName).t(dualFlow(end)) - data.(fieldName).t(dualFlow(1));
        data.(fieldName).TMR    = (((data.(fieldName).res(1).mDot(dualFlow)).^2)/(rho*A(1)))./(((data.(fieldName).res(2).mDot(dualFlow)).^2)/(rho*A(2)));
        data.(fieldName).Theta  = thetaPred(data.(fieldName).TMR)*180/pi;

        %data.(fieldName).TMRSummary     = [mean(data.(fieldName).TMR),min(data.(fieldName).TMR),max(data.(fieldName).TMR)];
        %data.(fieldName).ThetaSummary   = [mean(data.(fieldName).Theta),min(data.(fieldName).Theta),max(data.(fieldName).Theta)];

        data.(fieldName).mTotCalc   = trapz(data.(fieldName).t,data.(fieldName).res(1).mDot) + trapz(data.(fieldName).t,data.(fieldName).res(2).mDot);
        data.(fieldName).mTotMeas   = data.(fieldName).m(end) - data.(fieldName).m(1);
        if data.(fieldName).mTotMeas < 0
            try
                data.(fieldName).mTotMeas = data.(fieldName).m(data.(fieldName).res(n2).maxInd+expandRange) - data.(fieldName).m(1);
            catch
                data.(fieldName).mTotMeas = data.(fieldName).m(data.(fieldName).res(n2).maxInd) - data.(fieldName).m(1);
            end
        end

        strLen  = char("%-"+string(15*length(printResults))+"s\n");
        fprintf(strLen,"RESULTS: "+fieldName);
        for n2 = 1:length(printResults)
            fprintf('%15s',printResults{n2});
        end
        fprintf('\n');
        for n2 = 1:length(printResults)
            if length(data.(fieldName).(printResults{n2})) > 1
                fprintf('%15.2f',mean(data.(fieldName).(printResults{n2})));
            else
                fprintf('%15.2f',data.(fieldName).(printResults{n2}));
            end
        end
        fprintf('\n');

        % Produce Summary
        for n2 = 1:size(summary,2)
            if isfield(data.(fieldName),summary{1,n2})
                summary{n1+1,n2} = [mean(data.(fieldName).(summary{1,n2})),min(data.(fieldName).(summary{1,n2})),max(data.(fieldName).(summary{1,n2}))];
            else
                summary{n1+1,n2} = [mean(data.(fieldName).res(1).(summary{1,n2})(dualFlow)),min(data.(fieldName).res(1).(summary{1,n2})(dualFlow)),max(data.(fieldName).res(1).(summary{1,n2})(dualFlow)),...
                    mean(data.(fieldName).res(2).(summary{1,n2})(dualFlow)),min(data.(fieldName).res(2).(summary{1,n2})(dualFlow)),max(data.(fieldName).res(2).(summary{1,n2})(dualFlow))];
            end
        end
    end
    
    mTestComp(n1,1) = (1/10)*(trapz(data.(fieldName).t(dualFlow),data.(fieldName).res(1).mDot(dualFlow)) + trapz(data.(fieldName).t(dualFlow),data.(fieldName).res(2).mDot(dualFlow)));
    mTestComp(n1,2) = (1/10)*(trapz(data.(fieldName).t,data.(fieldName).res(1).mDot + data.(fieldName).res(2).mDot));
    mTestComp(n1,3) = (1/10)*(data.(fieldName).m(dualFlow(end))-data.(fieldName).m(dualFlow(1)));
    
    % UNCERTAINTY CALCULATION
    data.(fieldName).D = struct();
    data.(fieldName).Pavg     = [mean(data.(fieldName).res(1).P(dualFlow)),mean(data.(fieldName).res(2).P(dualFlow))];
    data.(fieldName).D.P = [(max(data.(fieldName).res(1).P(dualFlow)) - min(data.(fieldName).res(1).P(dualFlow)))/2,(max(data.(fieldName).res(2).P(dualFlow)) - min(data.(fieldName).res(2).P(dualFlow)))/2];
    % Relative Uncertainties
    data.(fieldName).D.PRel = data.(fieldName).D.P./data.(fieldName).Pavg;
    
    data.(fieldName).D.TMRRel   = sqrt(4*(DCd(1)^2+DA(1)^2+(Drho/2)^2+(data.(fieldName).D.PRel(1)/2)^2)+(Drho)^2+(DA(1))^2+4*(DCd(2)^2+DA(2)^2+(Drho/2)^2+(data.(fieldName).D.PRel(2)/2)^2)+(Drho)^2+(DA(2))^2);
end