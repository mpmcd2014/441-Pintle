clear; close all;
g = 9.8066;

thetaRes    = 1;    % deg

%% SPRAY CHARACTERISTICS
thetaR      = 60;       % deg
H           = 0.5;        % m
v0          = 30;       % m/s
dOrifice    = 0.0625/39.37;     % m
mDot        = 2*0.8*1/12;   % Amount of mass flow, assuming 80% concentrated at this point, and 1/12 a revolution

% DROPLET PROPERTIES
droplet.D   = 50e-6;
droplet.rho = 998;
droplet.sigma = 0.0717;   % N/m. Surface tension


%% SETUP
results = struct('vf',[0,0,0],'x',[],'y',[],'t',[],'thetaMeas',[],'Fmeas',[],'rf',[]);
%% LINEAR
linear = results;
linear.vf   = [v0*sind(thetaR),-v0*cosd(thetaR),v0];
linear.thetaMeas = atand(abs(linear.vf(1)/linear.vf(2)));
linear.rf   = H*tand(linear.thetaMeas);
linear.x    = [0,linear.rf];
linear.y    = [H,0];
linear.t    = sqrt((diff(linear.x)^2) + (diff(linear.y)^2))/v0;

linear.Fsense       = mDot*v0*cosd(linear.thetaMeas);
linear.Fmeas        = linear.Fsense/cosd(linear.thetaMeas);

%% GRAVITY
gravity = results;
gravity.vf  = [v0*sind(thetaR),-sqrt((v0*cosd(thetaR))^2+2*g*H),0];
gravity.vf(3) = rssq(gravity.vf(1:2));
%tf          = (v0*cosd(thetaR)-sqrt(v0^2*cosd(thetaR)^2-2*g*H))/(-g);
gravity.t      = linspace(0,(gravity.vf(2)-(-v0*cosd(thetaR)))/(-g));
gravity.x       = [0, gravity.vf(1)*gravity.t];
gravity.y       = [H, H - v0*cosd(thetaR)*gravity.t - 0.5*g*(gravity.t.^2)];

gravity.rf      = gravity.x(end);
gravity.thetaMeas = atand(gravity.rf/H);

gravity.Fsense  = mDot*gravity.vf(3)*cosd(atand(gravity.vf(1)/gravity.vf(2)));
gravity.Fmeas   = gravity.Fsense/cosd(gravity.thetaMeas);

%% DROPLET WITH DRAG
diameters   = [25e-6,50e-6,100e-6,250e-6,500e-6];
drag = repmat(results,length(diameters),1);
% FLOW PARAMETERS

% AIR PROPERTIES
% Cengel A-15, T = 25 deg C
air.rho = 1.184;                % kg/m^3
air.mu  = 1.849e-5;             % kg/m/s

for n1 = 1:length(diameters)
    droplet.D   = diameters(n1);
    [drag(n1).rf,tf,x,drag(n1).v,drag(n1).t] = SprayBreakupDrag(droplet, 'air', H, dOrifice, v0, thetaR, 1e-4);
    drag(n1).x  = x(:,1);
    drag(n1).y  = x(:,2);
    drag(n1).vf = drag(n1).v(end,:);
    
    drag(n1).thetaDrop  = atand(abs(drag(n1).v(:,1)./drag(n1).v(:,2)));
    drag(n1).thetaMeas  = atand(drag(n1).rf/H);

    drag(n1).Fsense     = mDot*drag(n1).vf(3)*cosd(atand(abs(drag(n1).vf(1)/drag(n1).vf(2))));
    drag(n1).Fmeas      = drag(n1).Fsense/cosd(drag(n1).thetaMeas);
end


%% PLOT COMPARISON
figure()
hold on
plot(linear.x,linear.y,'--k','LineWidth',1)
plot(gravity.x,gravity.y,'.-g','LineWidth',1)

names = {'Linear','Gravity'};
for n1 = 1:length(diameters)
    plot(drag(n1).x,drag(n1).y,'LineWidth',1)
    devPoint = find((thetaR-drag(n1).thetaDrop)>thetaRes,1,'first');
    %scatter(drag(n1).x(devPoint),drag(n1).y(devPoint),'x','MarkerEdgeColor','r');
    names(end+1) = {[char(string(diameters(n1)*1e6)),' um']};
end

legend(names)
a1 = gca();
a1.FontName = 'Times New Roman';
a1.FontSize = 14;
xlabel(a1,'x [m]')

figure()
bar(categorical(names),[linear.Fmeas,gravity.Fmeas,drag(:).Fmeas])
ylabel(a1,'y [m]')