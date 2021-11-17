%{
11/17/2021
Testing Visualization V1
%}
clear; close all;

%% READ DATA
filename = "Data summary.xlsx";

pintles     = ["HD0","HD0_5","HD1_0","HD2_0"];

data     = struct();
for n1 = 1:length(pintles)
    raw = xlsread(filename,pintles{n1},"A2:I25");
    data.(pintles(n1)).angles   = flip(raw(~isnan(raw(:,1)),1));
    data.(pintles(n1)).TMR      = raw(1,2:end);
    data.(pintles(n1)).TMR_H    = raw(2,2:end);
    data.(pintles(n1)).TMR_L    = raw(3,2:end);
    
        data.(pintles(n1)).ms   = flip(raw(7:end,2:end));
    TMRLabel = {};
    %{
    for n2 = 1:length(data.(pintles(n1)).TMR)
        TMRLabel(end+1) = {char(string(data.(pintles(n1)).TMR(n2)))};
    end
    %}
    TMRLabel = linspace(floor(min(data.(pintles(n1)).TMR*100))/100-0.01,ceil(max(data.(pintles(n1)).TMR*100))/100+0.01,30);
    data.(pintles(n1)).msPlot   = zeros(length(data.(pintles(n1)).angles),length(TMRLabel));
    for n2 = 1:length(data.(pintles(n1)).TMR)
        plotInd = find((data.(pintles(n1)).TMR(n2) - TMRLabel)<0,1,'first');
        data.(pintles(n1)).msPlot(:,plotInd) = data.(pintles(n1)).ms(:,n2);
    end
    %TMRLabel = string(TMRLabel);
    figure
    c = bar3(data.(pintles(n1)).angles,data.(pintles(n1)).msPlot);
    c(1).Parent.XTickLabel = string(TMRLabel(round(linspace(1,length(TMRLabel),length(c(1).Parent.XTickLabel)))));
    c(1).Parent.YDir = 'normal';
end

%% PLOT DATA
%{
for n1 = 1:length(pintles)
    bar3(
%}