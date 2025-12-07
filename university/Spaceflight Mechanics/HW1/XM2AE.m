function[AE]=XM2AE(ecc,XM)
% XM2AE
%
% da anomalia media ad anomalia eccentrica

eps=10^-7;
err=eps+1;
AE0=XM;
while err>eps
    AE=AE0-(AE0-ecc*sin(AE0)-XM)/(1-ecc*cos(AE0));
    err=abs(AE0-AE);
    AE0=AE;
end