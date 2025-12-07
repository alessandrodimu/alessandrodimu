function[XM]=AE2XM(ecc,AE)
% AE2XM
%
% da anomalia eccentrica ad anomalia media

if isscalar(ecc)==0 || isreal(ecc)==0 || ecc<0
    error('ecc deve essere uno scalare reale maggiore o uguale a 0');
end
if isscalar(AE)==0 || isreal(AE)==0
    error('AE deve essere uno scalare reale');
end

XM=AE-ecc*sin(AE);
end