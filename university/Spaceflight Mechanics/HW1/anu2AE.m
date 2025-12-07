function[AE]=anu2AE(ecc,anu)
% anu2AE
%
% da eccentricità e anomalia vera ad anomalia eccentrica

if isscalar(ecc)==0 || isreal(ecc)==0 || ecc<0
    error('ecc deve essere uno scalare reale maggiore o uguale a 0');
end
if isscalar(anu)==0 || isreal(anu)==0
    error('anu deve essere uno scalare reale');
end

% riscrivo l'angolo anu in [0,2*pi)
anu=wrapTo2Pi(anu);
if anu==2*pi
    anu=0;
end

supporto=anu/2;
if anu>=-pi/2 || anu<=pi/2
    AE=2*atan(sqrt((1-ecc)/(1+ecc))*tan(supporto));
else
    AE=2*atan(sqrt((1-ecc)/(1+ecc))*tan(supporto))+pi;
end

% riscrivo l'angolo AE in [0,2*pi)
AE=wrapTo2Pi(AE);
if AE==2*pi
    AE=0;
end

end