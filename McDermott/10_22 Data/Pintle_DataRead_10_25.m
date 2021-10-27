%{
10/24/2021
Initial Testing Analysis
%}
clc; clear; close all;

%% Independent Variables
rho = 998; % kg/m^3. Density of water at 70 deg F.
% Pintle Geometry
orif_d = 0.0625/39.37; % m. Orifice Diameter [in]
n_orif = 20;
A = pi*(orif_d/2)^2*n_orif; % m^2. Pintle Flow area [in^2]

% Annulus Geometry
ann_d = .767/39.37; % m
pint_out_d = 0.742/39.37; % m
%A = pi*(ann_d/2)^2 - pi*(pint_out_d/2)^2; % m^2

%% Read in Data
%%some intializations

filenamePrefix = "Data\Pintle\Pintle_10-22_";
filenameExtension = ".csv";
presVals = ["50","100", "125" ,"150"];
nRuns = [3,1,1,1];
data = struct();
start = {[9200, 10316, 12416], [9252], [10390], [10006] };
stop = {[15832, 15711, 14818], [15440], [18208], [13781]};
for i = 1:length(presVals)
    for n = 1:nRuns(i)  
    filename = filenamePrefix + presVals(i) +"-"+ string(n) + filenameExtension;
    tmpData=csvread(filename,1,0);
    fieldName = "psi" + presVals(i) + "_" + string(n);
    data.(fieldName).t = tmpData(:,1); % Time in seconds
    data.(fieldName).p = tmpData(:,7); % Pressure upstream of article in psi
    data.(fieldName).lc_1 = tmpData(:,16); % Load Cell (lbf)
    data.(fieldName).lc_2 = tmpData(:,17);
    data.(fieldName).lc_3 = tmpData(:,18);
    data.(fieldName).lc_4 = tmpData(:,19);
    data.(fieldName).m = (data.(fieldName).lc_1 + data.(fieldName).lc_2 +data.(fieldName).lc_3 + data.(fieldName).lc_4); %lbf
    clear tmpData;
    end
end
 %% Data Analysis
 for i = 1:length(presVals)
     for n = 1:nRuns(i)
        %Clean up Data
        filename = filenamePrefix + presVals(i) +"-"+ string(n) + filenameExtension;
        fieldName = "psi" + presVals(i) + "_" + string(n);
        data.(fieldName).t = data.(fieldName).t(start{i}(n):stop{i}(n));
        data.(fieldName).p = data.(fieldName).p(start{i}(n):stop{i}(n));
        data.(fieldName).m = data.(fieldName).m(start{i}(n):stop{i}(n));
        %Calculate values
        data.(fieldName).p_avg = mean(data.(fieldName).p)*6894.75; % Pa.
        data.(fieldName).mFlowRate = (data.(fieldName).m(end)-data.(fieldName).m(1))/(data.(fieldName).t(end)-data.(fieldName).t(1))*(1/2.205); % kg/s.
        data.(fieldName).Cd = data.(fieldName).mFlowRate/(A*sqrt(2*rho*data.(fieldName).p_avg));
        %Plots
        figure
        plot(data.(fieldName).t,data.(fieldName).p)
        hold on
        plot(data.(fieldName).t,data.(fieldName).m)
        title(fieldName)
        
        % UNCERTAINTY CALCULATION
        
        
        figure
        hold on
        plot(data.(fieldName).t, data.(fieldName).lc_1(start{i}(n):stop{i}(n)))
        plot(data.(fieldName).t, data.(fieldName).lc_2(start{i}(n):stop{i}(n)))
        plot(data.(fieldName).t, data.(fieldName).lc_3(start{i}(n):stop{i}(n)))
        plot(data.(fieldName).t, data.(fieldName).lc_4(start{i}(n):stop{i}(n)))
        legend('lc_1','lc_2','lc_3','lc_4')
        title(fieldName)
     end
 end
