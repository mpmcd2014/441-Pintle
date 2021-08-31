% 08/18/2021
% Mark McDermott



clear; close all

TMR = linspace(0,4);
curveNames = {};
%% LIQUID-LIQUID
    %% BLAKELY ET AL.
    BlakelyBF  = [0.25, 0.35, 0.45, 0.55];
    a   = [.9262,.8116,.7997,.6243];
    b   = [5.463,3.276,2.108,3.377];

    BlakelyFit = @(A,B,R) A*atand(B*R);

    figure
    hold on
    for n1 = 1:length(BlakelyBF)
        plot(TMR,BlakelyFit(a(n1),b(n1),TMR),'LineStyle','--','LineWidth',0.5)
        curveNames(end+1,1) = {"Blakely, BF = "+string(BlakelyBF(n1))};
    end
    %% ESCHER
    EscherFit = @(R) (105.5-24.5*R).*(R.^0.5);

    plot(TMR(TMR <= 1),EscherFit(TMR(TMR <= 1)),'-k','LineWidth',1.1);
    curveNames(end+1,1) = {"Escher, BF = 0.3-0.8"};

    %% CHENG
    ChengFit = @(R) acosd(1./(1+R));

    plot(TMR,ChengFit(TMR),'-r','LineWidth',1.1);
    curveNames(end+1,1) = {"Cheng, BF = 1"};
%% FORMATTING
a1 = gca();
a1.FontName = 'Times New Roman';
a1.FontSize = 12;
legend(curveNames)
xlabel('Total Momentum Ratio')
ylabel('Spray Cone Angle [deg]')

%% GAS-LIQUID
    %% SON
    %{ 
    "Design Procedure for a ..."

    theta_pt = 0;
    xi  = (90-xi)/90;
    s   = 1.15+1.35*xi;
    p   = 1.3+.9*xi;

    alpha   = 90*xi*exp((s-0.2)/(1+((K/90).^p))-s);
    %}
    %{
    "Effects of Momentum Ratio and weber Number..."
    We = [10,50,100];
    c1 = 38.36;
    c2 = -.096;
    SonFit  = @(R,We) c1*((R./We).^c2);
    for n1 = 1:length(We)
        plot(TMR,SonFit(TMR,We(n1)))
        curveNames(end+1,1) = {"Son, BF=1, We="+string(We)};
    end
    %}

    %% LEE


