%{
11/20/2021
Testing Visualization V1
%}
clear; close all;

%% READ DATA
filename    = "Data summary.xlsx";
savePlots   = true;

pintles     = ["HD0","HD0_5","HD1_0","HD1_5","HD2_0"];
HDs         = [0,0.5,1.0,1.5,2.0];

compareFun  = struct();
compareFun.Arctan   = @(T) atand(T);
compareFun.Blakely  = @(T) 0.6243*atan(3.377*T)*180/pi;
compareFun.Cheng    = @(T) acosd(1./(1+T));
names = fieldnames(compareFun);

minSize     = 60;
lineOption  = "peaks";
data     = struct();
lin     = {};
linColorMat     = [[0,0.73,0.18];[0.72,0.27,1];[0.3,0.75,0.93]];
%scatColorMat    = [[1,0,0];[0,1,0];[0.215,0.66,0.72];[0.85,0.33,0.10];[0.93,0.69,0.13]];
%scatColorMat    = [[0.8,0.53,0.26];[0.93,0.69,0.13];[0,1,0];[0.215,0.66,0.72];[0.57,0.49,0.66]];
scatColorMat    = [[1,0.6,0];[1,1,0];[0,1,0];[0,1,1];[0.35,0,1]];
fprintf('%-5s%-2s%10s%10s%10s%10s%10s%10s\n','H/D','|','Centroid','','Peaks','','Optical','')
fprintf('%-5s%-2s%10s%10s%10s%10s%10s%10s\n','','|','m','b','m','b','m','b')
for n1 = 1:length(pintles)
    raw = xlsread(filename,pintles{n1},"A2:G28");
    data.(pintles(n1)).angles   = flip(raw(~isnan(raw(:,1)),1));
    data.(pintles(n1)).TMR      = raw(1,2:end);
    data.(pintles(n1)).TMR_H    = raw(2,2:end);
    data.(pintles(n1)).TMR_L    = raw(3,2:end);
    data.(pintles(n1)).optical  = raw(5,2:end);
    data.(pintles(n1)).Doptical = raw(6,2:end);
    data.(pintles(n1)).ms   = flip(raw(8:end-2,2:end));
    data.(pintles(n1)).mExpect  = raw(end,2:end);
    
    data.(pintles(n1)).opticalP         = polyfit(data.(pintles(n1)).TMR,data.(pintles(n1)).optical,1);
    % Plotting
    figure();
    a(n1) = axes;       %#ok<*SAGROW>
    data.(pintles(n1)).totMass  = zeros(1,length(data.(pintles(n1)).TMR));
    data.(pintles(n1)).centroid = zeros(1,length(data.(pintles(n1)).TMR));
    data.(pintles(n1)).peaks    = zeros(1,length(data.(pintles(n1)).TMR));
    for n2 = 1:length(data.(pintles(n1)).TMR)
        binInd = ~isnan(data.(pintles(n1)).ms(:,n2));
        minVal = min(data.(pintles(n1)).ms(binInd,n2));
        [maxVal,maxInd]     = max(data.(pintles(n1)).ms(binInd,n2));
        data.(pintles(n1)).totMass(n2)      = sum(data.(pintles(n1)).ms(binInd,n2));
        validAngles     = data.(pintles(n1)).angles(binInd);
        
        data.(pintles(n1)).centroid(n2)     = sum(data.(pintles(n1)).ms(binInd,n2).*data.(pintles(n1)).angles(binInd))/data.(pintles(n1)).totMass(n2);
        data.(pintles(n1)).centroidP        = polyfit(data.(pintles(n1)).TMR,data.(pintles(n1)).centroid,1);
        data.(pintles(n1)).peaks(n2)        = validAngles(maxInd);
        data.(pintles(n1)).peaksP           = polyfit(data.(pintles(n1)).TMR,data.(pintles(n1)).peaks,1);
        
        sizeVec = 3000*(data.(pintles(n1)).ms(binInd,n2)-minVal)/data.(pintles(n1)).totMass(n2)+minSize;
        colorVec = sizeVec/max(sizeVec)*scatColorMat(n1,:);%,zeros(6,1)];%,zeros(6,1)];
        scat(n1,n2) = scatter(a(n1),repmat(data.(pintles(n1)).TMR(n2),1,6),data.(pintles(n1)).angles(binInd),sizeVec,colorVec,'filled','o');
        scat(n1,n2).LineWidth           = 1.5;
        scat(n1,n2).MarkerFaceAlpha     = 0.75;
        hold on
    end
    TMRBound    = [min(data.(pintles(n1)).TMR), max(data.(pintles(n1)).TMR)];
    switch lineOption
        case "centroid"
            lin{n1} = plot(TMRBound,polyval(data.(pintles(n1)).centroidP,TMRBound));
        case "peaks"
            lin{n1} = plot(TMRBound,polyval(data.(pintles(n1)).peaksP,TMRBound));
        otherwise
            lin{n1} = plot(TMRBound,polyval(data.(pintles(n1)).centroidP,TMRBound),TMRBound,polyval(data.(pintles(n1)).peaksP,TMRBound));
    end
    for n2 = 1:length(lin{n1})
        lin{n1}(n2).Color = [0,0,0];
        lin{n1}(n2).LineWidth = 3;
        lin{n1}(n2).LineStyle = ':';
    end
    fprintf('%-5.1f%-2s%+10.3f%+10.3f%+10.3f%+10.3f%+10.3f%+10.3f\n',HDs(n1),'|',data.(pintles(n1)).centroidP(1),data.(pintles(n1)).centroidP(2),data.(pintles(n1)).peaksP(1),data.(pintles(n1)).peaksP(2),data.(pintles(n1)).opticalP(1),data.(pintles(n1)).opticalP(2))
    
    e(n1) = errorbar(data.(pintles(n1)).TMR,data.(pintles(n1)).optical,data.(pintles(n1)).Doptical);
end

%% COMPARISON PLOT
figure
n1 = 1;

for n2 = 1:length(data.(pintles(n1)).TMR)
    binInd = ~isnan(data.(pintles(n1)).ms(:,n2));
    sizeVec = 3000*(data.(pintles(n1)).ms(binInd,n2)-minVal)/data.(pintles(n1)).totMass(n2)+minSize;
    colorVec = sizeVec/max(sizeVec)*scatColorMat(n1,:);%,zeros(6,1)];%,zeros(6,1)];
    nscat(n2) = scatter(repmat(data.(pintles(n1)).TMR(n2),1,6),data.(pintles(n1)).angles(binInd),sizeVec,colorVec,'filled','o');
    nscat(n2).LineWidth           = 1.5;
    nscat(n2).MarkerFaceAlpha     = 0.45;
    hold on
end
switch lineOption
    case "centroid"
        lin{n1} = plot(TMRBound,polyval(data.(pintles(n1)).centroidP,TMRBound));
    case "peaks"
        lin{n1} = plot(TMRBound,polyval(data.(pintles(n1)).peaksP,TMRBound));
    otherwise
        lin{n1} = plot(TMRBound,polyval(data.(pintles(n1)).centroidP,TMRBound),TMRBound,polyval(data.(pintles(n1)).peaksP,TMRBound));
end
for n2 = 1:length(lin{n1})
    %lin{n1}(n2).Color = linColorMat(n2,:);
    lin{n1}(n2).Color = [0,0,0];
    lin{n1}(n2).LineWidth = 3;
    lin{n1}(n2).LineStyle = ':';
end
f = errorbar(data.(pintles(n1)).TMR,data.(pintles(n1)).optical,data.(pintles(n1)).Doptical);
f.Marker    = 'd';
f.MarkerSize= 10;
f.LineWidth = 1.5;
f.LineStyle = 'none';
f.MarkerFaceColor = [1,0,0];
f.Color     = [1,0,0];
x = linspace(0.4,1.2);
for n2 = 1:length(names)
    c(n2) = plot(x,compareFun.(names{n2})(x),'LineWidth',2.25);
    c(n2).Color = linColorMat(n2,:);
end

b = gca();
b.FontName = 'Times New Roman';
b.FontSize = 18;
xlim(b,[0.3,1.25]);
ylim(b,[10,80]);
xlabel(b,'\itTMR')
ylabel(b,'{\it{\theta}} [\circ]','Interpreter','tex')
title(b,"H/D = "+string(HDs(n1)))


%% Legend: Literature
figure
hold on
for n2 = 1:length(names)
    c(n2) = plot(x,compareFun.(names{n2})(x),'LineWidth',2.25);
    c(n2).Color = linColorMat(n2,:);
end
compLgd = legend(c(1,:),{'Classical','Blakely','Cheng'},'Location','southeast','FontName','Times New Roman','FontSize',24);
title(compLgd,'Literature');

%% Legend: Present work
figure
n1 = 1;
for n2 = 1:length(data.(pintles(n1)).TMR)
    binInd = ~isnan(data.(pintles(n1)).ms(:,n2));
    sizeVec = 3000*(data.(pintles(n1)).ms(binInd,n2)-minVal)/data.(pintles(n1)).totMass(n2)+minSize;
    colorVec = sizeVec/max(sizeVec)*scatColorMat(n1,:);%,zeros(6,1)];%,zeros(6,1)];
    nscat(n2) = scatter(repmat(data.(pintles(n1)).TMR(n2),1,6),data.(pintles(n1)).angles(binInd),sizeVec,colorVec,'filled','o');
    nscat(n2).LineWidth           = 1.5;
    nscat(n2).MarkerFaceAlpha     = 0.45;
    hold on
end
switch lineOption
    case "centroid"
        lin{n1} = plot(TMRBound,polyval(data.(pintles(n1)).centroidP,TMRBound));
    case "peaks"
        lin{n1} = plot(TMRBound,polyval(data.(pintles(n1)).peaksP,TMRBound));
    otherwise
        lin{n1} = plot(TMRBound,polyval(data.(pintles(n1)).centroidP,TMRBound),TMRBound,polyval(data.(pintles(n1)).peaksP,TMRBound));
end
for n2 = 1:length(lin{n1})
    %lin{n1}(n2).Color = linColorMat(n2,:);
    lin{n1}(n2).Color = [0,0,0];
    lin{n1}(n2).LineWidth = 3;
    lin{n1}(n2).LineStyle = ':';
end
f = errorbar(data.(pintles(n1)).TMR,data.(pintles(n1)).optical,data.(pintles(n1)).Doptical);
f.Marker    = 'd';
f.MarkerSize= 10;
f.LineWidth = 1.5;
f.LineStyle = 'none';
f.MarkerFaceColor = [1,0,0];
f.Color     = [1,0,0];
presLgd = legend([nscat(end) lin{1} f],{'Mass Collected','Mass Best Fit','Optical'},'Location','southeast','FontName','Times New Roman','FontSize',24);
title(presLgd,'Present Work');
%}
%% FORMATTING AND I/O
for n1 = 1:length(a)
    a(n1).FontName = 'Times New Roman';
    a(n1).FontSize = 18;
    xlim(a(n1),[0.3,1.25]);
    ylim(a(n1),[10,80]);
    xlabel(a(n1),'\itTMR')
    ylabel(a(n1),'{\it{\theta}} [\circ]','Interpreter','tex')
    if mod(HDs(n1),1) == 0
        title(a(n1),"H/D = "+string(HDs(n1))+".0")
    else
        title(a(n1),"H/D = "+string(HDs(n1)))
    end
    
    e(n1).Marker    = 'd';
    e(n1).MarkerSize= 10;
    e(n1).LineWidth = 1.5;
    e(n1).LineStyle = 'none';
    e(n1).MarkerFaceColor = [1,0,0];
    e(n1).Color     = [1,0,0];
end

if savePlots
    tmp = findobj('Type','figure');
    currPath = "Graph Images\"+string(datetime('now','Format','MM-dd_hh-mm_a'));
    mkdir(currPath)
    for n1 = 1:length(tmp)-3
        %saveas(tmp(end-n1+1),currPath+"\"+pintles(n1)+".png");
        print(tmp(end-n1+1),pintles(n1),'-dsvg')
    end
end