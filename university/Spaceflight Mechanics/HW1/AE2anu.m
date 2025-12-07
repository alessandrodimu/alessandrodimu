function[anu]=AE2anu(ecc,AE)
% AE2anu
%
% da anomalia eccentrica ad anomalia vera
% IMPORTANTE: controllare bene gli angoli!!!

if isscalar(ecc)==0 || isreal(ecc)==0 || ecc<0
    error('ecc deve essere uno scalare reale maggiore o uguale a 0');
end
if isscalar(AE)==0 || isreal(AE)==0
    error('AE deve essere un valore reale');
end

anu=2*atan(sqrt((1+ecc)/(1-ecc))*tan(AE/2));
end

