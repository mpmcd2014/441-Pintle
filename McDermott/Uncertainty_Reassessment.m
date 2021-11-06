clear;

delm    = 0.03;
delAa   = 0.015;
delAr   = 0.015;
delRho  = 0.001;
delP    = 0.03;

%% PROGRESS REPORT 2
fprintf('Progress Report 2\n');
delCda   = sqrt(delm^2+delAa^2+(delRho/2)^2+(delP/2)^2);
delCdr   = sqrt(delm^2+delAr^2+(delRho/2)^2+(delP/2)^2);
delTMR  = sqrt((4*delm^2 + 9*delAa^2 + 3*delRho^2 + 2*delP^2)+(4*delm^2 + 9*delAr^2 + 3*delRho^2 + 2*delP^2));
fprintf('%15s%10.4f\n','del Cda = ',delCda);
fprintf('%15s%10.4f\n','del Cdr = ',delCdr);
fprintf('%15s%10.4f\n','del TMR = ',delTMR);

%% With t-estimator
fprintf('T-Estimator\n');
delCda  = 0.022*1.782;
delCdr  = 0.028*1.782;
delTMR  = sqrt((delP^2 + (2*delCda)^2 + (delAa^2)) + (delP^2 + (2*delCdr)^2 + (delAr^2)));
fprintf('%-15s%10.4f\n','del Cda = ',delCda);
fprintf('%-15s%10.4f\n','del Cdr = ',delCdr);
fprintf('%-15s%10.4f\n','del TMR = ',delTMR);