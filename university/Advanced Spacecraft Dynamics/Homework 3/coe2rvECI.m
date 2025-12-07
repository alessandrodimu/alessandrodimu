function[rVectECI,vVectECI]=coe2rvECI(COE,amu)
% coe2rvECI
%
% da parametri orbitali a posizione e velocità in ECI

[n,m]=size(COE);
if isvector(COE)==0 || n~=6 || m~=1
    error('COE deve essere un vettore colonna 6x1');
end
if isscalar(amu)==0 || isreal(amu)==0 || amu<0
    error('amu deve essere uno scalare reale positivo')
end
    
% ottengo gli elementi orbitali classici dal vettore COE
a=COE(1);
ecc=COE(2);
AINC=COE(3);
GOM=COE(4);
POM=COE(5);
anu=COE(6);

p=a*(1-ecc^2);
r=p/(1+ecc*cos(anu));
rVectPER=[r*cos(anu);r*sin(anu);0];
vVectROT=sqrt(amu/p)*[ecc*sin(anu);(1+ecc*cos(anu)),;0];
ROT2PER=[cos(anu) -sin(anu) 0;sin(anu) cos(anu) 0;0 0 1];
vVectPER=ROT2PER*vVectROT;
PER2ECI=[cos(POM)*cos(GOM)-sin(POM)*cos(AINC)*sin(GOM) -sin(POM)*cos(GOM)-cos(POM)*cos(AINC)*sin(GOM) sin(AINC)*sin(GOM);...
    cos(POM)*sin(GOM)+sin(POM)*cos(AINC)*cos(GOM) -sin(POM)*sin(GOM)+cos(POM)*cos(AINC)*cos(GOM) -sin(AINC)*cos(GOM);...
    sin(POM)*sin(AINC) cos(POM)*sin(AINC) cos(AINC)];
rVectECI=PER2ECI*rVectPER;
vVectECI=PER2ECI*vVectPER;
end