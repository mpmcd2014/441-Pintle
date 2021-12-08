function [rf,tf,x,v,t,F,theta] = SprayDragPrediction_DEBUG(droplet, medium, h0, v0, theta0, dt)
% DRAGPREDICTION Predicted droplet path in Earth gravity in air.


if nargin == 5
    dt = 0.001;
end

iterLim = 1e5;
g = 9.807;                      % m/s^2

% Ambient fluid selections
if ischar(medium) || isstring(medium)
    switch medium
        case 'air'
            medium = struct();
            % AIR PROPERTIES
            % Cengel A-15, T = 25 deg C
            medium.rho = 1.184;                % kg/m^3
            medium.mu  = 1.849e-5;             % kg/m/s
    end
end
% https://pages.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2016.pdf
Cd = @(Re) (24/Re) + ((2.6*(Re/5))/(1+((Re/5)^1.52))) + ((0.411*((Re/2.63e5)^-7.94))/(1+((Re/2.63e5)^-8))) + (.25*(Re/1e6)/(1+(Re/1e6)));

% DROPLET PROPERTIES
m   = droplet.rho.*(4/3)*pi*((droplet.D/2)^3);     % kg. Droplet mass.
A   = (pi/4)*(droplet.D^2);                % m^2. Droplet drag area.

% Initialization
N   = floor(fzero(@(t) -.5*g*t^2 - v0*cosd(theta0)*t + h0,1)/dt + 1);
t   = zeros(N,1);
theta = repmat(theta0,N,1);
x   = repmat([0,h0],N,1);
v   = repmat([v0*sind(theta0),-v0*cosd(theta0),v0],N,1);
F   = zeros(N,1);
Cd_D = zeros(N,1);

n1 = 1;
while (x(n1,2) > 0) && (n1 <= iterLim)
    Re_D    = medium.rho*v(n1,3)*droplet.D/medium.mu;           % Droplet Reynolds number
    Cd_D(n1)    = Cd(Re_D);                         % Droplet drag coefficient
    F(n1)   = .5*Cd_D(n1)*medium.rho*(v(n1,3)^2)*A;      % Drag force on droplet
    
    % Predict velocities
    v(n1+1,1)   = v(n1,1) - (F(n1)/m)*dt*sind(theta(n1));
    v(n1+1,2)   = v(n1,2) + (F(n1)/m)*dt*cos(theta(n1)) - g*dt;
    v(n1+1,3)   = rssq(v(n1+1,1:2));
    
    theta(n1+1)     = atand(abs(v(n1+1,1)/v(n1+1,2)));
    
    % Predict positions
    x(n1+1,1:2)   = x(n1,1:2) + v(n1+1,1:2)*dt;
    t(n1+1)     = t(n1) + dt;
    
    n1 = n1+1;
end

% Remove unused elements
x   = x(1:n1-1,:);
v   = v(1:n1-1,:);
theta = theta(1:n1-1,:);
t   = t(1:n1-1,:);
F = F(1:n1-1,:);

rf  = x(end,1);
tf  = t(end);
    
