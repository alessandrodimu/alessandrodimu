function [angoli_eul,t]=aiuto_es6_3QUA(vettore_q_ini,vettore_rot,tempo)

global k1 k2 k3 M1 M2 M3
Tol0 = 1e-12; 
Tol1 = 1e-09;
options = odeset('RelTol',Tol0,'AbsTol',Tol1);
conv=pi/180;


q0=vettore_q_ini(1);
q1=vettore_q_ini(2);
q2=vettore_q_ini(3);
q3=vettore_q_ini(4);

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


% integro usando il set di quaternioni
[t,Wa]=ode45(@int_qua0,[0 tempo],[w1 w2 w3 q0 q1 q2 q3],options);
nel=size(Wa);
angoli_eul=zeros(nel(1),3);
for ind=1:1:nel(1)
    n1=[1-2*(Wa(ind,6)^2+Wa(ind,7)^2) ; ...
        2*(Wa(ind,6)*Wa(ind,5)-Wa(ind,4)*Wa(ind,7)) ; ...
        2*(Wa(ind,7)*Wa(ind,5)+Wa(ind,4)*Wa(ind,6))];
    n2=[2*(Wa(ind,5)*Wa(ind,6)+Wa(ind,4)*Wa(ind,7)) ; ...
        1-2*(Wa(ind,5)^2+Wa(ind,7)^2) ; ...
        2*(Wa(ind,7)*Wa(ind,6)-Wa(ind,4)*Wa(ind,5))];
    n3=[2*(Wa(ind,5)*Wa(ind,7)-Wa(ind,4)*Wa(ind,6)) ; ...
        2*(Wa(ind,6)*Wa(ind,7)+Wa(ind,4)*Wa(ind,5)) ; ...
        1-2*(Wa(ind,5)^2+Wa(ind,6)^2)];
    BN=[n1 n2 n3];
    [phi,the,psi]=decod_3(BN);
    angoli_eul(ind,1)=phi/conv;
    angoli_eul(ind,2)=the/conv;
    angoli_eul(ind,3)=psi/conv;
end

