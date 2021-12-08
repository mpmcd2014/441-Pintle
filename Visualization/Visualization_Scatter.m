%{
11/17/2021
Testing Visualization V1
%}
clear; close all;

%% READ DATA
filename = "Data summary.xlsx";

pintles     = ["HD0","HD0_5","HD1_0","HD1_5","HD2_0"];
HDs         = [0,0.5,1.0,1.5,2.0];

minSize     = 60;
data     = struct();
a = repmat(axes,length(pintles),1);
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
    a(n1) = axes;
    data.(pintles(n1)).totMass = zeros(1,length(data.(pintles(n1)).TMR));
    for n2 = 1:length(data.(pintles(n1)).TMR)
        binInd = ~isnan(data.(pintles(n1)).ms(:,n2));
        minVal = min(data.(pintles(n1)).ms(binInd,n2));
        maxVal = max(data.(pintles(n1)).ms(binInd,n2));
        data.(pintles(n1)).totMass(n2) = sum(data.(pintles(n1)).ms(binInd,n2));
        sizeVec = 3000*(data.(pintles(n1)).ms(binInd,n2)-minVal)/data.(pintles(n1)).totMass(n2)+minSize;
        colorVec = [sizeVec/max(sizeVec),zeros(6,1),zeros(6,1)];
        b = scatter(a(n1),repmat(data.(pintles(n1)).TMR(n2),1,6),data.(pintles(n1)).angles(binInd),sizeVec,colorVec,'filled','s');%,'s','MarkerEdgeColor','none');
        %data.(pintles(n1)).msPlot(:,plotInd) = data.(pintles(n1)).ms(:,n2);
        hold on
        fprintf('|%9.2f%9s|',data.(pintles(n1)).TMR(n2),'')
    end
    fprintf('\n');
    for n2 = 1:length(data.(pintles(n1)).TMR)
        fprintf('|%9.2f%9.2f|',data.(pintles(n1)).mExpect(n2),data.(pintles(n1)).totMass(n2));
    end
    fprintf('\n==================================================\n');
    %TMRLabel = string(TMRLabel);
end

for n1 = 1:length(a)
    a(n1).FontName = 'Times New Roman';
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