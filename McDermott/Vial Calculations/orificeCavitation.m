function [Cd,K] = orificeCavitation(rho, mu, d, L, r, P1, P2, Pvap, output)
%GETCAVITATION Calculates cavitation conditions through an orifice.
%   Calculates cavitation behavior of flow through an orifice using the
%   relations described in the Fluent Theory Guide (FTG), Section 16.12.1,
%   the Plain Orifice Model.

Re = (d.*rho./mu).*sqrt( (2*(P1 - P2))./rho);       % Reynolds number. Fluent Theory Guide, 16-331
K = (P1 - Pvap) ./ (P1 - P2);                       % Cavitation number. FTG 16-332

Cc = 1./(sqrt((1/(0.611^2)) - (11.4*r./d)));        % Coefficient of Contraction. FTG 16-333

Kincep = 1.9*(1-(r./d)).^2 - 1000./Re;                                  % Cavitation number for inception of cavitation. For short nozzles, Kincep ~1.9. FTG 16-337
Kcrit = 1 + (1./( (1+(L./(4*d))).*(1+(2000./Re)).*exp(70*r./d)));       % Kcrit, where flipped nozzle flow occurs (i.e. liquid jet surrounded by downstream gas). FTG 16-338

if K <= Kincep
    if K < Kcrit
        code = 1;
        Cd = 0.611;
    else
        code = 2;
        Cd = Cc.*sqrt(K);
    end
else
    if K < Kcrit
        code = 3;
        Cd = 0.611;
    else
        code = 4;
        Cd = 1./( (1./(0.827-0.0085*(L./d))) + 20*( (1+2.25*(L./d))./Re));      % FTG 16-339,16-340
    end
end

if output > 0
    fprintf('Cavitation Number:\t\t%.2f\n', K);
    switch code
        case 1
            fprintf('The flow is flipped. Potential flow analysis gives theoretical value for Cd.\n')       % FTG 16-342
        case 2
            fprintf('The flow is cavitating.\n')
        case 3
            fprintf('The flow is flipped. Potential flow analysis gives theoretical value for Cd.\n')       % FTG 16-342
        case 4
            fprintf('The flow is single phase.\n')
    end
    fprintf('Cd = %.3f\n', Cd);
end

end

