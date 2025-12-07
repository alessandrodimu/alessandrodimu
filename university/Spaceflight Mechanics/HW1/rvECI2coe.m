function[COE]=rvECI2coe(rVect,vVect,amu)
% rvECI2coe
%
% da posizione e velocità in ECI a parametri orbitali
% IMPORTANTE: ci sono problemi quanto nVers(2)==0 ,quando eVers(3)==0 e
% quando dot(rVect,vVect)==0
%
% provare a riscrivere dopo!!!

[n1,m1]=size(rVect);
if isvector(rVect)==0 || n1~=3 || m1~=1
    error('rVect deve essere un vettore colonna 3x1');
end
[n2,m2]=size(vVect);
if isvector(vVect)==0 || n2~=3 || m2~=1
    error('vVect deve essere un vettore colonna 3x1');
end
if isscalar(amu)==0 || isreal(amu)==0 || amu<0
    error('amu deve essere uno scalare reale positivo')
end

hVect=cross(rVect,vVect);
hVers=hVect/norm(hVect);
h=norm(hVect);

% eccentricità
eVect=cross(vVect,hVect)/amu-rVect/norm(rVect);
eVers=eVect/norm(eVect);
ecc=norm(eVect);

% semiasse maggiore
p=h^2/amu;
a=p/(1-ecc^2);

% inclinazione è definita in [0,pi]
% non ci sono problemi sul segno
AINC=acos(dot(hVers,[0,0,1]));

nVers=cross([0,0,1],hVers)/norm(cross([0,0,1],hVers));

% RAAN è definita in [0,2*pi]
% ci sono problemi sul segno, Immagine(acos) in [0,pi]
% guardiamo la 2a componente del versore n
GOM=acos(dot([1,0,0],nVers));
if nVers(2)>0
    GOM=+GOM;
elseif nVers(2)<0
    GOM=2*pi-GOM;
end

% argomento del perigeo è definito in [0,2*pi];
% ci sono problemi sul segno, Immagine(acos) in [0,pi]
% guardiamo la 3a componente del versore e
POM=acos(dot(nVers,eVers));
if eVers(3)>0
    POM=+POM;
elseif eVers(3)<0
    POM=2*pi-POM;
end

% anomalia vera è definita in [0,2*pi];
% ci sono problemi sul segno, Immagine(acos) in [0,pi]
% guardiamo il prodotto tra rVect e vVect
anu=acos(dot(eVers,rVect)/norm(rVect));
if dot(rVect,vVect)>0
    anu=+anu;
elseif dot(rVect,vVect)<0
    anu=2*pi-anu;
end

% scrivo il vettore COE
COE=[a ; ecc ; AINC ; GOM ; POM ; anu];
end