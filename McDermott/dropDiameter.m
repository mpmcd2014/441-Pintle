function d_32 = dropDiameter(dH, v, rho, surfTens)
%DROPDIAMETER Estimates diameter of droplets from flow through an orifice.
% dH = hole diameter
% rho = fluid density
% surfTens = surface tension
% d_32 = Sauter mean diameter

lambda = dH./8;

We = (rho.*(v.^2).*lambda)./surfTens;

d_32 = 133.0*lambda.*(We.^-0.74);

end

