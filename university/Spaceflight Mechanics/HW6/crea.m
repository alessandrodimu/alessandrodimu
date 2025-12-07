function [BN]= crea 
conv=pi/180;

ferma=1;
while ferma == 1
disp('Modalità di input');
disp('1) Coseni direttori');
disp('2) Elementi della rotazione principale');
disp('3) Angoli di Eulero');
disp('4) Quaternioni');
corpo=input('Scegliere una delle precedenti voci: ');
ferma=0;
switch corpo
    case 1
        a1=input('Longitudine Versore b1:');
        d1=input('Latitudine  Versore b1:');
        a2=input('Longitudine Versore b2:');
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
        
    case 2
        a1=input('Longitudine  Versore  e : ');
        d1=input('Latitudine   Versore  e : ');
        phi=input ('Angolo di Rotazione phi : ');
        a1=a1*conv;
        d1=d1*conv;
        phi=phi*conv;
        E(1)=cos(d1)*cos(a1);
        E(2)=cos(d1)*sin(a1);
        E(3)=sin(d1);
        BN=cos(phi)*eye(3)+(1-cos(phi))*E'*E-sin(phi)*skew(E);

    case 3
        phi=input('Roll   phi : ');
        the=input('Pitch  the : ');
        psi=input('Yaw    psi : ');
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
    
    case 4
        a1=input('Longitudine  Versore  e : ');
        d1=input('Latitudine   Versore  e : ');
        phi=input('Angolo di Rotazione phi : ');
        a1=a1*conv;
        d1=d1*conv;
        phi=phi*conv;
        q0=cos(phi/2);
        sph2=sin(phi/2);
        q(1)=sph2*cos(d1)*cos(a1);
        q(2)=sph2*cos(d1)*sin(a1);
        q(3)=sph2*sin(d1);
        BN=(q0*q0-q*q')*eye(3)+2*q'*q-2*q0*skew(q);
        
otherwise
        disp('Voce non presente');
        ferma=1;
end
end

