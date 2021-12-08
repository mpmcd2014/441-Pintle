%{
11/20/2021
Testing Visualization V1
%}
clear; close all;

%% READ DATA
filename = "Data summary.xlsx";

pintles     = ["HD0","HD0_5","HD1_0","HD1_5","HD2_0"];
HDs         = [0,0.5,1.0,1.5,2.0];

minSize     = 60;
lineOption  = "";
data     = struct();
lin     = {};
lineColors = [[0,0.45,0.74];[0,0,0]];
fprintf('%-10s%-5s%10s%10s%10s%10s\n','H/D','|','Centroid','','Peaks','')
fprintf('%-10s%-5s%10s%10s%10s%10s\n','','|','m','b','m','b')
for n1 = 1:length(pintles)
    raw = xlsread(filename,pintles{n1},"A2:G27");
    data.(pintles(n1)).angles   = flip(raw(~isnan(raw(:,1)),1));
    data.(pintles(n1)).TMR      = raw(1,2:end);
    data.(pintles(n1)).TMR_H    = raw(2,2:end);
    data.(pintles(n1)).TMR_L    = raw(3,2:end);
    data.(pintles(n1)).ms   = flip(raw(7:end-2,2:end));
    data.(pintles(n1)).mExpect  = raw(end,2:end);
    
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
        colorVec = [zeros(6,1),zeros(6,1),sizeVec/max(sizeVec)];%,zeros(6,1)];%,zeros(6,1)];
        scat(n1,n2) = scatter(a(n1),repmat(data.(pintles(n1)).TMR(n2),1,6),data.(pintles(n1)).angles(binInd),sizeVec,colorVec,'filled','s');%,'s','MarkerEdgeColor','none');
        %data.(pintles(n1)).msPlot(:,plotInd) = data.(pintles(n1)).ms(:,n2);
        hold on
    end
    TMRBound    = [min(data.(pintles(n1)).TMR), max(data.(pintles(n1)).TMR)];
    switch lineOption
        case "centroid"
            lin{n1} = plot(TMRBound,polyval(data.(pintles(n1)).centroidP,TMRBound),'b','LineWidth',1.5);
        case "peaks"
            lin{n1} = plot(TMRBound,polyval(data.(pintles(n1)).peaksP,TMRBound),'b','LineWidth',1.5);
        otherwise
            lin{n1} = plot(TMRBound,polyval(data.(pintles(n1)).centroidP,TMRBound),TMRBound,polyval(data.(pintles(n1)).peaksP,TMRBound));
    end
    for n2 = 1:length(lin{n1})
        lin{n1}(n2).Color = lineColors(n2,:);
        lin{n1}(n2).LineWidth = 1.5;
    end
    
    fprintf('%-10.1f%-5s%+10.3f%+10.3f%+10.3f%+10.3f\n',HDs(n1),'|',data.(pintles(n1)).centroidP(1),data.(pintles(n1)).centroidP(2),data.(pintles(n1)).peaksP(1),data.(pintles(n1)).peaksP(2))
end

for n1 = 1:length(a)
    a(n1).FontName = 'Times New Roman';
    %a(n1).FaceAlpha = 0.5;
    a(n1).FontSize = 14;
    xlabel(a(n1),'\itTMR')
    ylabel(a(n1),'$\theta$','Interpreter','latex')
    title(a(n1),"H/D = "+string(HDs(n1)))
end
%% PLOT DATA
%{
for n1 = 1:length(pintles)
    bar3(
%}