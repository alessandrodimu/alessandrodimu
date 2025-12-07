function [angoli_eul,t]=aiuto_es6_3EUL(vettore_eul_ini,vettore_rot,tempo)

global k1 k2 k3 M1 M2 M3
Tol0 = 1e-12; 
Tol1 = 1e-09;
options = odeset('RelTol',Tol0,'AbsTol',Tol1);
conv=pi/180;


phi=vettore_eul_ini(1)*conv;
the=vettore_eul_ini(2)*conv;
psi=vettore_eul_ini(3)*conv;

w1=vettore_rot(1);
w2=vettore_rot(2);
w3=vettore_rot(3);

%Manca l'inerzia del corpo
%k1=(j2-j3)/j1;
%k2=(j3-j1)/j2;
%k3=(j1-j2)/j3;
k1=0;
k2=0;
k3=0;

M1=0;
M2=0;
M3=0;


% integro usando il set di angoli di Eulero 3-2-1
[t,Wa]=ode45(@int_eul0,[0 tempo],[w1 w2 w3 phi the psi],options);
angoli_eul(:,1)=wrapTo180(Wa(:,4)/conv);
angoli_eul(:,2)=wrapTo180(Wa(:,5)/conv);
angoli_eul(:,3)=wrapTo180(Wa(:,6)/conv);
end
