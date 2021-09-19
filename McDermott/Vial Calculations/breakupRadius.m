function rb = breakupRadius(rho1,rho2,v,L,surfTens)


We  = (rho1.*(v.^2).*L)./surfTens;

rb = zeros(size(rho1));

l1 = We<=800;
rb(l1) = 0.5*L(l1)*0.167.*We(l1);
rb(~l1) = 0.5*L(~l1)*14.2*((rho2(~l1)./rho1(~l1)).^(-2/3))*(We(~l1).^(-1/3));