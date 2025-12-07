% ESERCITAZIONE 6
% PARTE 1

% rapp_AED = rappresentazione alpha e delta
% rapp_ERP = rappresentazione elementi rotazione principale
% rapp_EUL = rappresentazione angoli di eulero 3-2-1
% rapp_QUA = rappresentazione quaternioni

% lavoro con MatricecCnew.xls
rapp_AED_0=[40 20 -15];

% SIT_1 = assetto iniziale
[BN_1,rapp_AED_1,rapp_ERP_1,rapp_EUL_1,rapp_QUA_1]=aiuto_es6_1(rapp_AED_0,'AED');
% SIT_2 = assetto dopo rotazione DPHI
[BN_2,rapp_AED_2,rapp_ERP_2,rapp_EUL_2,rapp_QUA_2]=aiuto_es6_1(rapp_ERP_1+[0 0 10],'ERP');
% SIT_3 = assetto dopo rotazione DPHI_EUELERO
[BN_3,rapp_AED_3,rapp_ERP_3,rapp_EUL_3,rapp_QUA_3]=aiuto_es6_1(rapp_EUL_2+[20 0 0],'EUL');
% SIT_4 = assetto dopo rotazione DTHETA_EULERO
[BN_4,rapp_AED_4,rapp_ERP_4,rapp_EUL_4,rapp_QUA_4]=aiuto_es6_1(rapp_EUL_3+[0 30 0],'EUL');
% SIT_5 = assetto dopo rotazione DPSI_EULERO/assetto finale
[BN_5,rapp_AED_5,rapp_ERP_5,rapp_EUL_5,rapp_QUA_5]=aiuto_es6_1(rapp_EUL_4+[0 0 40],'EUL');



% ESERCITAZIONE 6
% PARTE 2

% calcolo velocita di rotazione nel riferimento corpo
alfE_1=deg2rad(rapp_ERP_1(1));
delE_1=deg2rad(rapp_ERP_1(2));
E_1=[cos(alfE_1)*cos(delE_1) sin(alfE_1)*cos(delE_1) sin(delE_1)];
w12=10*pi/180*E_1;

w23=10*pi/180*[1 0 0];

phi_3=deg2rad(rapp_EUL_3(1));
w34=10*pi/180*[0 cos(phi_3) -sin(phi_3)];

phi_4=deg2rad(rapp_EUL_4(1));
the_4=deg2rad(rapp_EUL_4(2));
w45=10*pi/180*[-sin(the_4) sin(phi_4)*cos(the_4) cos(phi_4)*cos(the_4)];



% ESERCITAZIONE 6
% PARTE 3

% integro le equazioni cinematiche usando i quaternioni
[angoli_QUA_12,t_QUA_12]=aiuto_es6_3QUA(rapp_QUA_1,w12,1);
[angoli_QUA_23,t_QUA_23]=aiuto_es6_3QUA(rapp_QUA_2,w23,2);
[angoli_QUA_34,t_QUA_34]=aiuto_es6_3QUA(rapp_QUA_3,w34,3);
[angoli_QUA_45,t_QUA_45]=aiuto_es6_3QUA(rapp_QUA_4,w45,4);

angoli_QUA_tot=[angoli_QUA_12 ; angoli_QUA_23(2:end,1:3) ; ...
    angoli_QUA_34(2:end,1:3) ; angoli_QUA_45(2:end,1:3)];
t_QUA_tot=[t_QUA_12 ; 1+t_QUA_23(2:end) ; 3+t_QUA_34(2:end) ; 6+t_QUA_45(2:end)];

figure(1)
plot(t_QUA_tot,angoli_QUA_tot(:,1),'-b',t_QUA_tot,angoli_QUA_tot(:,2),'-.g',t_QUA_tot,angoli_QUA_tot(:,3),'--r')
title('Andamento di phi, theta e psi (integrazione con quaternioni)');
legend(texlabel('phi'),texlabel('theta'),texlabel('psi'));
xlabel('time (s)')
ylabel(texlabel('phi theta psi (deg)'))
axis([0 10 -180 180])
set(gca,'YTick',[-180 -90 0 90 180],'YTickMode','manual')
grid on
hold on
plot(1*ones(500),linspace(-180,180,500),'--k','HandleVisibility','off');
plot(3*ones(500),linspace(-180,180,500),'--k','HandleVisibility','off');
plot(6*ones(500),linspace(-180,180,500),'--k','HandleVisibility','off');
hold off

% integro le equazioni cinematiche usando gli angoli di Eulero 3-2-1
[angoli_EUL_12,t_EUL_12]=aiuto_es6_3EUL(rapp_EUL_1,w12,1);
[angoli_EUL_23,t_EUL_23]=aiuto_es6_3EUL(rapp_EUL_2,w23,2);
[angoli_EUL_34,t_EUL_34]=aiuto_es6_3EUL(rapp_EUL_3,w34,3);
[angoli_EUL_45,t_EUL_45]=aiuto_es6_3EUL(rapp_EUL_4,w45,4);

angoli_EUL_tot=[angoli_EUL_12 ; angoli_EUL_23(2:end,1:3) ; ...
    angoli_EUL_34(2:end,1:3) ; angoli_EUL_45(2:end,1:3)];
t_EUL_tot=[t_EUL_12 ; 1+t_EUL_23(2:end) ; 3+t_EUL_34(2:end) ; 6+t_EUL_45(2:end)];

figure(2)
plot(t_EUL_tot,angoli_EUL_tot(:,1),'-b',t_EUL_tot,angoli_EUL_tot(:,2),'-.g',t_EUL_tot,angoli_EUL_tot(:,3),'--r')
title('Andamento di phi, theta e psi (integrazione con angoli di Eulero 3-2-1)');
legend(texlabel('phi'),texlabel('theta'),texlabel('psi'));
xlabel('time (s)')
ylabel(texlabel('phi theta psi (deg)'))
axis([0 10 -180 180])
set(gca,'YTick',[-180 -90 0 90 180],'YTickMode','manual')
grid on
hold on
plot(1*ones(500),linspace(-180,180,500),'--k','HandleVisibility','off');
plot(3*ones(500),linspace(-180,180,500),'--k','HandleVisibility','off');
plot(6*ones(500),linspace(-180,180,500),'--k','HandleVisibility','off');
hold off



% ESERCITAZIONE 6
% PARTE 4


% aumento le varie rotazioni del caso precedente
[BN_1p,rapp_AED_1p,rapp_ERP_1p,rapp_EUL_1p,rapp_QUA_1p]=aiuto_es6_1(rapp_AED_0,'AED');
[BN_2p,rapp_AED_2p,rapp_ERP_2p,rapp_EUL_2p,rapp_QUA_2p]=aiuto_es6_1(rapp_ERP_1p+[0 0 100],'ERP');
[BN_3p,rapp_AED_3p,rapp_ERP_3p,rapp_EUL_3p,rapp_QUA_3p]=aiuto_es6_1(rapp_EUL_2p+[200 0 0],'EUL');
[BN_4p,rapp_AED_4p,rapp_ERP_4p,rapp_EUL_4p,rapp_QUA_4p]=aiuto_es6_1(rapp_EUL_3p+[0 300 0],'EUL');
[BN_5p,rapp_AED_5p,rapp_ERP_5p,rapp_EUL_5p,rapp_QUA_5p]=aiuto_es6_1(rapp_EUL_4p+[0 0 400],'EUL');

% aggiorno le varie velocita angolari
w12p=10*w12;
w23p=10*w23;
w34p=10*w34;
w45p=10*w45;

% integro le equazioni cinematiche usando i quaternioni
[angoli_QUA_12p,t_QUA_12p]=aiuto_es6_3QUA(rapp_QUA_1p,w12p,1);
[angoli_QUA_23p,t_QUA_23p]=aiuto_es6_3QUA(rapp_QUA_2p,w23p,2);
[angoli_QUA_34p,t_QUA_34p]=aiuto_es6_3QUA(rapp_QUA_3p,w34p,3);
[angoli_QUA_45p,t_QUA_45p]=aiuto_es6_3QUA(rapp_QUA_4p,w45p,4);

angoli_QUA_totp=[angoli_QUA_12p ; angoli_QUA_23p(2:end,1:3) ; ...
    angoli_QUA_34p(2:end,1:3) ; angoli_QUA_45p(2:end,1:3)];
t_QUA_totp=[t_QUA_12p ; 1+t_QUA_23p(2:end) ; 3+t_QUA_34p(2:end) ; 6+t_QUA_45p(2:end)];

figure(3)
plot(t_QUA_totp,angoli_QUA_totp(:,1),'-b',t_QUA_totp,angoli_QUA_totp(:,2),'-.g',t_QUA_totp,angoli_QUA_totp(:,3),'--r')
title('Andamento di phi, theta e psi (integrazione con quaternioni)');
legend(texlabel('phi'),texlabel('theta'),texlabel('psi'));
xlabel('time (s)')
ylabel(texlabel('phi theta psi (deg)'))
axis([0 10 -180 180])
set(gca,'YTick',[-180 -90 0 90 180],'YTickMode','manual')
grid on
hold on
plot(1*ones(500),linspace(-180,180,500),'--k','HandleVisibility','off');
plot(3*ones(500),linspace(-180,180,500),'--k','HandleVisibility','off');
plot(6*ones(500),linspace(-180,180,500),'--k','HandleVisibility','off');
hold off

% integro le equazioni cinematiche usando gli angoli di Eulero 3-2-1
[angoli_EUL_12p,t_EUL_12p]=aiuto_es6_3EUL(rapp_EUL_1p,w12p,1);
[angoli_EUL_23p,t_EUL_23p]=aiuto_es6_3EUL(rapp_EUL_2p,w23p,2);
[angoli_EUL_34p,t_EUL_34p]=aiuto_es6_3EUL(rapp_EUL_3p,w34p,3);
[angoli_EUL_45p,t_EUL_45p]=aiuto_es6_3EUL(rapp_EUL_4p,w45p,4);

angoli_EUL_totp=[angoli_EUL_12p ; angoli_EUL_23p(2:end,1:3) ; ...
    angoli_EUL_34p(2:end,1:3) ; angoli_EUL_45p(2:end,1:3)];
t_EUL_totp=[t_EUL_12p ; 1+t_EUL_23p(2:end) ; 3+t_EUL_34p(2:end) ; 6+t_EUL_45p(2:end)];

figure(4)
plot(t_EUL_totp,angoli_EUL_totp(:,1),'-b',t_EUL_totp,angoli_EUL_totp(:,2),'-.g',t_EUL_totp,angoli_EUL_totp(:,3),'--r')
title('Andamento di phi, theta e psi (integrazione con angoli di Eulero 3-2-1)');
legend(texlabel('phi'),texlabel('theta'),texlabel('psi'));
xlabel('time (s)')
ylabel(texlabel('phi theta psi (deg)'))
axis([0 10 -180 180])
set(gca,'YTick',[-180 -90 0 90 180],'YTickMode','manual')
grid on
hold on
plot(1*ones(500),linspace(-180,180,500),'--k','HandleVisibility','off');
plot(3*ones(500),linspace(-180,180,500),'--k','HandleVisibility','off');
plot(6*ones(500),linspace(-180,180,500),'--k','HandleVisibility','off');
hold off
