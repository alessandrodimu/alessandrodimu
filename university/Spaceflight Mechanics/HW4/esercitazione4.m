% ESERCITAZIONE 4
% PARTE 1

% ipotesi semplificative considerate
% - orbite dei pianeti circolari e complanari
% - tempi di permanenza nelle sfere di influenza dei pianeti trascurabili

% dati necessari
R_terra=149.5*10^6;
R_marte=227.8*10^6;
r_terra=6371;
r_marte=3397;
Ts_terra=86400;
Ts_marte=88775;
mu_terra=398600;
mu_marte=42828;
mu_sole=1.327*10^11;
V_terra=sqrt(mu_sole/R_terra);
V_marte=sqrt(mu_sole/R_marte);
T_terra=2*pi*sqrt(R_terra^3/mu_sole);
T_marte=2*pi*sqrt(R_marte^3/mu_sole);
n_terra=2*pi/T_terra;
n_marte=2*pi/T_marte;

% trasferta di Hohmann tra Terra e Marte (andata e ritorno)
a_hohmann=(R_terra+R_marte)/2;
e_hohmann=(R_marte-R_terra)/(R_marte+R_terra);
V_departure_1=sqrt(mu_sole/a_hohmann*(1+e_hohmann)/(1-e_hohmann));
Voo_departure=abs(V_departure_1-V_terra);

V_arrival_1=sqrt(mu_sole/a_hohmann*(1-e_hohmann)/(1+e_hohmann));
Voo_arrival=abs(V_arrival_1-V_marte);

TOF_hohmann=pi*sqrt(a_hohmann^3/mu_sole);
T_sinodico=2*pi/abs(n_marte-n_terra);
gamma_1_ini=pi-n_marte*TOF_hohmann;
gamma_1_fin=pi-n_terra*TOF_hohmann;
gamma_1_ini_ritorno=-gamma_1_fin;
t_wait=(-2*gamma_1_fin-2*pi)/(n_marte-n_terra);
TOF_1=t_wait+2*TOF_hohmann;



% ESERCITAZIONE 4
% PARTE 2

% dimostrazione altezza minima 200 km per sforzo propulsivo...
r_leo=linspace(100,500,41)+r_terra;
V_leo=sqrt(mu_terra./r_leo);
V_hyp_per_2=sqrt(Voo_departure^2+2*mu_terra./r_leo);
V0_eq=sqrt(2*(mu_terra/r_terra-mu_terra./(2*r_leo)));

deltaV_2_ideale=V_hyp_per_2-V_leo;
deltaV_2=deltaV_2_ideale+V0_eq;

figure(1)
plot(r_leo-r_terra,deltaV_2,'b-','LineWidth',2);
title('Andamento del costo di immissione');
xlabel('quota z (km)');
ylabel('DV2 (km/s)');
grid on

figure(2)
plot(r_leo-r_terra,V_leo.^2,'r-','LineWidth',2);
title('Andamento della resistenza (circa)');
xlabel('quota z (km)');
ylabel('{V^2}_{LEO} (km/s)');
grid on

% sensibilità
fattore_sensibilita=(V_hyp_per_2/Voo_departure).^2;
figure(3)
plot(r_leo-r_terra,fattore_sensibilita,'LineWidth',2);
title('Andamento del fattore di sensibilita in funzione della quota');
xlabel('quota z (km)');
ylabel('fattore di sensibilita (V0 / Voo)^2');
grid on

% ipotizzo un errore assoluto sulla velocita di iniezione di 0.1 km/s
err_assoluto_V0=0.1;
err_relativo_V0=err_assoluto_V0./V_hyp_per_2;
err_relativo_Voo=fattore_sensibilita.*err_relativo_V0;
err_assoluto_Voo=err_relativo_Voo*Voo_departure;

% scelgo i valori corrispondenti a 200 km dalla superficie terrestre
r_leo=r_leo(11);
V_leo=V_leo(11);
V_hyp_per_2=V_hyp_per_2(11);
V0_eq=V0_eq(11);
deltaV_2_ideale=deltaV_2_ideale(11);
deltaV_2=deltaV_2(11);
err_assoluto_Voo=err_assoluto_Voo(11);

% proprieta geometriche iperbole di partenza
energia_hyp_2=Voo_departure^2/2;
momento_hyp_2=V_hyp_per_2*r_leo;
e_hyp_2=sqrt(1+2*energia_hyp_2*momento_hyp_2^2/mu_terra^2);
beta_hyp_2=acos(1/e_hyp_2);
a_hyp_2=-mu_terra/(2*energia_hyp_2);

% calcolo massa dello spacecraft diretto verso Marte
% la massa utile e l'unica che viene spedita verso Marte
g0=9.81*10^-3;
M_ini_2=4000;
Isp=290;
eff_str=10;

M_fin_2=M_ini_2*exp(-deltaV_2_ideale/(Isp*g0));
M_p_2=M_ini_2-M_fin_2;
M_s_2=M_p_2/(eff_str-1);
M_u_2=M_fin_2-M_s_2;



% ESERCITAZIONE 4
% PARTE 3

% considero l'arrivo su Marte nei seguenti casi
% (A) ---> atterraggio morbido
% (B) ---> orbita circolare di 3597 km
% (C1) --> orbita circolare con periodo pari a 1 giorno solare 
%          marziano raggiunta direttamente
% (C2) --> stessa orbita circolare raggiunta con manovra a due impulsi
% (D) ---> stessa orbita circolare raggiunta con manovra a tre impulsi
%          (con il secondo raggio di 100000 km)
% (E) ---> orbita ellittica con periodo 1 giorno solare marziano e
%          periastro pari a 3597 km

% Punto A
V_hyp_per_3A=sqrt(Voo_arrival^2+2*mu_marte/r_marte);
deltaV_3A=V_hyp_per_3A;

% Punto B
r_3B=3597;
V_hyp_per_3B=sqrt(Voo_arrival^2+2*mu_marte/r_3B);
V_3B=sqrt(mu_marte/r_3B);
deltaV_3B=V_hyp_per_3B-V_3B;

% Punto C1
r_3C1=((Ts_marte/(2*pi))^2*mu_marte)^(1/3);
V_hyp_per_3C1=sqrt(Voo_arrival^2+(2*mu_marte/r_3C1));
V_3C1=sqrt(mu_marte/r_3C1);
deltaV_3C1=V_hyp_per_3C1-V_3C1;

% Punto C2
r_per_3C2=r_3B;
r_apo_3C2=r_3C1;
a_3C2=(r_per_3C2+r_apo_3C2)/2;
e_3C2=(r_apo_3C2-r_per_3C2)/(r_per_3C2+r_apo_3C2);
V_per_3C2=sqrt(mu_marte/a_3C2*(1+e_3C2)/(1-e_3C2));
V_apo_3C2=sqrt(mu_marte/a_3C2*(1-e_3C2)/(1+e_3C2));

deltaV_3C2_1=V_hyp_per_3B-V_per_3C2;
deltaV_3C2_2=V_3C1-V_apo_3C2;
deltaV_3C2=deltaV_3C2_1+deltaV_3C2_2;

% Punto D
r_ini_3D=r_3B;
r_mid_3D=100000;
r_fin_3D=r_3C1;
a_3D_1=(r_mid_3D+r_ini_3D)/2;
e_3D_1=(r_mid_3D-r_ini_3D)/(r_mid_3D+r_ini_3D);
a_3D_2=(r_mid_3D+r_fin_3D)/2;
e_3D_2=(r_mid_3D-r_fin_3D)/(r_mid_3D+r_fin_3D);

V_ell_per_3D_1=sqrt(mu_marte/a_3D_1*(1+e_3D_1)/(1-e_3D_1));
V_ell_apo_3D_1=sqrt(mu_marte/a_3D_1*(1-e_3D_1)/(1+e_3D_1));
V_ell_per_3D_2=sqrt(mu_marte/a_3D_2*(1+e_3D_2)/(1-e_3D_2));
V_ell_apo_3D_2=sqrt(mu_marte/a_3D_2*(1-e_3D_2)/(1+e_3D_2));

deltaV_3D_1=V_hyp_per_3B-V_ell_per_3D_1;
deltaV_3D_2=V_ell_apo_3D_2-V_ell_apo_3D_1;
deltaV_3D_3=V_ell_per_3D_2-V_3C1;
deltaV_3D=deltaV_3D_1+deltaV_3D_2+deltaV_3D_3;

% Punto E
a_3E=r_3C1;
r_per_3E=r_3B;
r_apo_3E=2*a_3E-r_per_3E;
e_3E=(r_apo_3E-r_per_3E)/(r_apo_3E+r_per_3E);
V_ell_per_3E=sqrt(mu_marte/a_3E*(1+e_3E)/(1-e_3E));
deltaV_3E=V_hyp_per_3B-V_ell_per_3E;



% ESERCITAZIONE 4
% PARTE 4

% fly-by di Marte implica un cambio di inclinazione del piano orbitale
r_per_4=200+r_marte;
delta=2*asin(1/(1+r_per_4*Voo_arrival^2/mu_marte));
DI=asin(Voo_arrival*sin(delta)/V_marte);
DI_max=asin(Voo_arrival/V_marte);



% ESERCITAZIONE 4
% PARTE 5

% calcolo nuova composizione della massa dello spacecraft
M_u_5=0.99*M_u_2;
dM_p=0.01*M_u_2/(1+1/(eff_str-1));
dM_s=0.01*M_u_2-dM_p;
M_p_5=M_p_2+dM_p;
M_s_5=M_s_2+dM_s;
M_ini_5=M_u_5+M_p_5+M_s_5;
M_fin_5=M_u_5+M_s_5;

% calcolo finestra di lancio
DeltaV_5=Isp*g0*log(M_ini_5/M_fin_5);
V_hyp_per_5=DeltaV_5+sqrt(mu_terra/r_leo);
Voo_5=sqrt(V_hyp_per_5^2-2*mu_terra/r_leo);
V_departure_5=V_terra+Voo_5;
momento_5=V_departure_5*R_terra;
p_5=momento_5^2/mu_sole;
e_5=p_5/R_terra-1;
a_5=p_5/(1-e_5^2);
ni_5=[acos((p_5-R_marte)/(R_marte*e_5)) ...
    2*pi-acos((p_5-R_marte)/(R_marte*e_5))];
AE_5=2*atan(sqrt((1-e_5)/(1+e_5))*tan(ni_5/2));
AE_5(2)=AE_5(2)+2*pi;
TOF_5=sqrt(a_5^3/mu_sole)*(AE_5-e_5*sin(AE_5));
gamma_5=ni_5-TOF_5*n_marte;
D_GAMMA=gamma_5(1)-gamma_5(2);
DT=D_GAMMA/(n_terra-n_marte)/86400;

% disegno le traiettorie di apertura e chiusura finestra di lancio
anomalia_terra=n_terra*DT*86400/2;

tab_terra=Orbit2D_Ell(R_terra,0,0,0,'cerchio');
tab_marte=Orbit2D_Ell(R_marte,0,0,0,'cerchio');
tab1=Orbit2D_Ell(a_5,e_5,TOF_5(1),0,'sole');
tab2=Orbit2D_Ell(a_5,e_5,TOF_5(2),0,'sole');
tab1=ROT_2D(tab1,anomalia_terra);
tab2=ROT_2D(tab2,-anomalia_terra);

figure(4)
hold on
plot(0,0,'ko','MarkerSize',12);
plot(tab_terra(1,:),tab_terra(2,:),'k');
plot(tab_marte(1,:),tab_marte(2,:),'k');
plot(tab1(1,:),tab1(2,:),'b');
plot(tab1(1,1),tab1(2,1),'bo','MarkerFaceColor','b');
plot(tab1(1,end),tab1(2,end),'bo','MarkerFaceColor','b');
plot(tab2(1,:),tab2(2,:),'r');
plot(tab2(1,1),tab2(2,2),'ro','MarkerFaceColor','r');
plot(tab2(1,end),tab2(2,end),'ro','MarkerFaceColor','r');
hold off
title('Missione Terra-Marte: finestra di lancio');
daspect([1 1 1]);
grid on