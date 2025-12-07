function [BN,rapp_AED,rapp_ERP,rapp_EUL,rapp_QUA]=aiuto_es6_1(rapp_0,modo)
% matrice BN = matrice assetto
% rapp_AED = rappresentazione alpha e delta
% rapp_ERP = rappresentazione elementi rotazione principale
% rapp_EUL = rappresentazione angoli di eulero 3-2-1
% rapp_QUA = rappresentazione quaternioni

conv=pi/180;

supporto=[1 2 3 4];
modi_disp={'AED','ERP','EUL','QUA'};
confronto=strcmp(modo,modi_disp);
index=supporto(confronto==1);
if index==1             % modalita AED
    a1=rapp_0(1);
    d1=rapp_0(2);
    a2=rapp_0(3);
    a1=a1*conv;
    d1=d1*conv;
    a2=a2*conv;
    B1(1)=cos(d1)*cos(a1);
    B1(2)=cos(d1)*sin(a1);
    B1(3)=sin(d1);
    %inserire protezione per d1=0
    b21=cos(a2);
    b22=sin(a2);
    b23=-(B1(1)*b21+B1(2)*b22)/B1(3);
    b2=sqrt(b21^2+b22^2+b23^2);
    B2(1)=b21/b2;
    B2(2)=b22/b2;
    B2(3)=b23/b2;
    B3=cross(B1,B2);
    BN=[B1; B2; B3];
elseif index==2         % modalita ERP
    a1=rapp_0(1);
    d1=rapp_0(2);
    phi=rapp_0(3);
    a1=a1*conv;
    d1=d1*conv;
    phi=phi*conv;
    E(1)=cos(d1)*cos(a1);
    E(2)=cos(d1)*sin(a1);
    E(3)=sin(d1);
    BN=cos(phi)*eye(3)+(1-cos(phi))*E'*E-sin(phi)*skew(E);
elseif index==3         % modalita EUL
    phi=rapp_0(1);
    the=rapp_0(2);
    psi=rapp_0(3);
    phi=phi*conv;
    the=the*conv;
    psi=psi*conv;
    sphi=sin(phi);
    cphi=cos(phi);
    sthe=sin(the);
    cthe=cos(the);
    spsi=sin(psi);
    cpsi=cos(psi);
    B1=[cthe*cpsi cthe*spsi -sthe];
    B2=[sphi*sthe*cpsi-cphi*spsi sphi*sthe*spsi+cphi*cpsi sphi*cthe];
    B3=[cphi*sthe*cpsi+sphi*spsi cphi*sthe*spsi-sphi*cpsi cphi*cthe];
    BN=[B1; B2; B3];
elseif index==4         % modalita QUA
    a1=rapp_0(1);
    d1=rapp_0(2);
    phi=rapp_0(3);
    a1=a1*conv;
    d1=d1*conv;
    phi=phi*conv;
    q0=cos(phi/2);
    sph2=sin(phi/2);
    q(1)=sph2*cos(d1)*cos(a1);
    q(2)=sph2*cos(d1)*sin(a1);
    q(3)=sph2*sin(d1);
    BN=(q0*q0-q*q')*eye(3)+2*q'*q-2*q0*skew(q);    
end

% decodifiche
[alf1,del1,alf2,del2,alf3,del3]=decod_1(BN);
[E,Phi]=decod_2(BN);
[phi,the,psi]=decod_3(BN);
[q0,q1,q2,q3]=decod_4(BN);

% conversioni
alf1=rad2deg(alf1);
del1=rad2deg(del1);
alf2=rad2deg(alf2);
del2=rad2deg(del2);
alf3=rad2deg(alf3);
del3=rad2deg(del3);

alfE=rad2deg(atan2(E(2),E(1)));
delE=rad2deg(atan(E(3)/sqrt(E(2)^2+E(1)^2)));
Phi=rad2deg(Phi);

phi=rad2deg(phi);
the=rad2deg(the);
psi=rad2deg(psi);

% risultati
rapp_AED=[alf1,del1,alf2,del2,alf3,del3];
rapp_ERP=[alfE,delE,Phi];
rapp_EUL=[phi,the,psi];
rapp_QUA=[q0,q1,q2,q3];

end
