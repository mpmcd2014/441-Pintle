%% ANNULUS AREA
clear;

ro  = 0.767/39.37;
ri  = 0.742/39.37;

del     = 0.001;

delArea = (2*ro/(ro^2-ri^2))*sqrt((ro*del)^2+(ri*del)^2);
fprintf('%-30s = %.2f%%\n','Area Relative Uncertainty',delArea*100);

%% TMR
clear;
delM    = 0.1;
delRho  = 0.1;
delAa   = 4.3;
delAr   = 1.6;

delTMR  = sqrt(4*(delM^2)+2*(delRho^2)+delAa^2+delAr^2);
fprintf('%-30s = %.2f%%\n','TMR Relative Uncertainty',delTMR);