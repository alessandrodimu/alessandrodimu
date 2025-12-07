close all
clear all
clc
global k1 k2 k3 M1 M2 M3
Tol0 = 1e-12; 
Tol1 = 1e-09;
options = odeset('RelTol',Tol0,'AbsTol',Tol1);
conv=pi/180;

[BN]=crea;

[c11,c21,c31,c12,c22,c23]=decod_0(BN);
[alf1,del1,alf2,del2,alf3,del3]=decod_1(BN);
[E,Phi]=decod_2(BN);
[phi,the,psi]=decod_3(BN);
[q0,q1,q2,q3]=decod_4(BN);

w1=E(1);
w2=E(2);
w3=E(3);
W=[w1 ; w2 ; w3];
tempo=4*pi;

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


% integro usando il set di coseni direttori
[t,Wa1]=ode45(@int_cos0,[0 tempo],[w1 w2 w3 c11 c21 c31 c12 c22 c23],options);
%figure (0);
%plot(t,Wa1(:,4),'-',t,Wa1(:,5),'-.',t,Wa1(:,6),'.')
%figure (1);
%plot(t,Wa1(:,7),'-',t,Wa1(:,8),'-.',t,Wa1(:,9),'.')
disp('caso 1');
disp(size(Wa1));
nel=size(Wa1);
Wb1=zeros(nel(1),3);
for ind=1:1:nel(1)
    n1=[Wa1(ind,4); Wa1(ind,5); Wa1(ind,6)];
    n2=[Wa1(ind,7); Wa1(ind,8); Wa1(ind,9)];
    n3=cross(n1,n2);
    BN_1=[n1 n2 n3];
    [phi,the,psi]=decod_3(BN_1);
    Wb1(ind,1)=phi/conv;
    Wb1(ind,2)=the/conv;
    Wb1(ind,3)=psi/conv;
end
figure;
subplot(2,2,1)
plot(t,Wb1(:,1),'-b',t,Wb1(:,2),'-.g',t,Wb1(:,3),'--r')
legend(texlabel('phi'),texlabel('theta'),texlabel('psi'));
title('Set di coseni direttori')
xlabel('time (s)')
ylabel(texlabel('phi theta psi (deg)'))
axis([0 12.6 -180 180])
set(gca,'YTick',[-180 -90 0 90 180],'YTickMode','manual')


% integro usando il set di elementi di rotazione principale
[t,Wa2]=ode45(@int_erp0,[0 tempo],[w1 w2 w3 Phi w1 w2 w3],options);
disp('caso 2');
disp(size(Wa2));
nel=size(Wa2);
Wb2=zeros(nel(1),3);
for ind=1:1:nel(1)
    E_2=[Wa2(ind,5) ; Wa2(ind,6) ; Wa2(ind,7)]/norm([Wa2(ind,5),Wa2(ind,6),Wa2(ind,7)]);
    BN_2=cos(Wa2(ind,4))*eye(3)+(1-cos(Wa2(ind,4)))*E_2*(E_2)'-sin(Wa2(ind,4))*skew(E_2);
    [phi_2,the_2,psi_2]=decod_3(BN_2);
    Wb2(ind,1)=phi_2/conv;
    Wb2(ind,2)=the_2/conv;
    Wb2(ind,3)=psi_2/conv;
end
subplot(2,2,2)
plot(t,Wb2(:,1),'-b',t,Wb2(:,2),'-.g',t,Wb2(:,3),'--r')
legend(texlabel('phi'),texlabel('theta'),texlabel('psi'));
title('Set di elementi di rotazione principale')
xlabel('time (s)')
ylabel(texlabel('phi theta psi (deg)'))
axis([0 12.6 -180 180])
set(gca,'YTick',[-180 -90 0 90 180],'YTickMode','manual')


% integro usando il set di angoli di Eulero 3-2-1
[t,Wa3]=ode45(@int_eul0,[0 tempo],[w1 w2 w3 phi the psi],options);
disp('caso 3');
disp(size(Wa3));
ne1=size(Wa3);
Wb3(:,1)=wrapTo180(Wa3(:,4)/conv);
Wb3(:,2)=wrapTo180(Wa3(:,5)/conv);
Wb3(:,3)=wrapTo180(Wa3(:,6)/conv);
subplot(2,2,3)
plot(t,Wb3(:,1),'-b',t,Wb3(:,2),'-.g',t,Wb3(:,3),'--r')
legend(texlabel('phi'),texlabel('theta'),texlabel('psi'));
title('Set di angoli di Eulero 3-2-1')
xlabel('time (s)')
ylabel(texlabel('phi theta psi (deg)'))
axis([0 12.6 -180 180])
set(gca,'YTick',[-180 -90 0 90 180],'YTickMode','manual')

n1=[];
n2=[];
n3=[];


% integro usando il set di quaternioni
[t,Wa4]=ode45(@int_qua0,[0 tempo],[w1 w2 w3 q0 q1 q2 q3],options);
disp('caso 4');
disp(size(Wa4));
nel=size(Wa4);
Wb4=zeros(nel(1),3);
for ind=1:1:nel(1)
    n1=[1-2*(Wa4(ind,6)^2+Wa4(ind,7)^2) ; ...
        2*(Wa4(ind,6)*Wa4(ind,5)-Wa4(ind,4)*Wa4(ind,7)) ; ...
        2*(Wa4(ind,7)*Wa4(ind,5)+Wa4(ind,4)*Wa4(ind,6))];
    n2=[2*(Wa4(ind,5)*Wa4(ind,6)+Wa4(ind,4)*Wa4(ind,7)) ; ...
        1-2*(Wa4(ind,5)^2+Wa4(ind,7)^2) ; ...
        2*(Wa4(ind,7)*Wa4(ind,6)-Wa4(ind,4)*Wa4(ind,5))];
    n3=[2*(Wa4(ind,5)*Wa4(ind,7)-Wa4(ind,4)*Wa4(ind,6)) ; ...
        2*(Wa4(ind,6)*Wa4(ind,7)+Wa4(ind,4)*Wa4(ind,5)) ; ...
        1-2*(Wa4(ind,5)^2+Wa4(ind,6)^2)];
    BN_4=[n1 n2 n3];
    [phi_4,the_4,psi_4]=decod_3(BN_4);
    Wb4(ind,1)=phi_4/conv;
    Wb4(ind,2)=the_4/conv;
    Wb4(ind,3)=psi_4/conv;
end
subplot(2,2,4)
plot(t,Wb4(:,1),'-b',t,Wb4(:,2),'-.g',t,Wb4(:,3),'--r')
legend(texlabel('phi'),texlabel('theta'),texlabel('psi'));
title('Set di quaternioni')
xlabel('time (s)')
ylabel(texlabel('phi theta psi (deg)'))
axis([0 12.6 -180 180])
set(gca,'YTick',[-180 -90 0 90 180],'YTickMode','manual')


% conversioni per disp
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


% decodifiche
supporto1=[alf1,del1,alf2,del2,alf3,del3];
fprintf('alpha1\t | delta1\t | alpha2\t | delta2\t | alpha3\t | delta3\n');
fprintf('%5.2f\t | %5.2f\t | %5.2f\t | %5.2f\t | %5.2f\t | %5.2f\n',supporto1);
fprintf('\n');

supporto2=[alfE,delE,Phi];
fprintf('alphaE\t | deltaE\t | Phi\n');
fprintf('%5.2f\t | %5.2f\t | %5.2f\n',supporto2);
fprintf('\n');

supporto3=[phi,the,psi];
fprintf('phi\t\t | theta\t | psi\n');
fprintf('%5.2f\t | %5.2f\t | %5.2f\n',supporto3);
fprintf('\n');

supporto4=[q0,q1,q2,q3];
fprintf('q0\t\t | q1\t\t | q2\t\t | q3\n');
fprintf('%5.4f\t | %5.4f\t | %5.4f\t | %5.4f\n',supporto4);
fprintf('\n');