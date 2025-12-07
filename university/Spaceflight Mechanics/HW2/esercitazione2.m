% ESERCITAZIONE 2
% PARTE 1

mu=398600;
T_siderale=86164;
T_solare=86400;
LON=-80.6*pi/180;
deltaPhi=28.5*pi/180;

% informazioni su leo di partenza
r_leo=6600;
v_leo=sqrt(mu/r_leo);
T_leo=2*pi*sqrt(r_leo^3/mu);

% informazioni su geo di arrivo
r_geo=((T_siderale/(2*pi))^2*mu)^(1/3);
v_geo=sqrt(mu/r_geo);
T_geo=T_siderale;

% informazioni su gto di trasferimento
a_gto=(r_leo+r_geo)/2;
e_gto=(r_geo-r_leo)/(r_leo+r_geo);
v_gto_per=sqrt(mu/a_gto)*sqrt((1+e_gto)/(1-e_gto));
v_gto_apo=sqrt(mu/a_gto)*sqrt((1-e_gto)/(1+e_gto));
T_gto=2*pi*sqrt(a_gto^3/mu);

% 01/01/2008 @ 00:00:00 --> MJD = 54466
% 31/01/2008 @ 00:00:00 --> MJD = 54496
t_0=54466*T_solare;                     % 01/01/2008 @ 00:00:00
t_1=(54496+10/24)*T_solare;             % 31/01/2008 @ 10:00:00
t_2=t_1+T_leo/4;                        % 31/01/2008 @ 10:22:14
t_3=t_2+T_gto/2;                        % 31/01/2008 @ 15:37:50
t_volo=t_3-t_1;

alpha_green_0=100*pi/180;
alpha_green=alpha_green_0+2*pi/T_siderale*(t_1-t_0);
RAAN=wrapTo2Pi(LON+alpha_green-pi/2);

% Strategia 1.1 (A)
deltaV_1A=v_gto_per-v_leo;
deltaV_2A=(v_geo-v_gto_apo)+(2*v_geo*sin(deltaPhi/2));
deltaV_totA=deltaV_1A+deltaV_2A;

% Strategia 1.2 (B)
deltaV_1B=v_gto_per-v_leo;
deltaV_2B=(2*v_gto_apo*sin(deltaPhi/2))+(v_geo-v_gto_apo);
deltaV_totB=deltaV_1B+deltaV_2B;

% Strategia 1.3 (C)
deltaV_1C=v_gto_per-v_leo;
deltaV_2C=sqrt(v_gto_apo^2+v_geo^2-2*v_gto_apo*v_geo*cos(deltaPhi));
deltaV_totC=deltaV_1C+deltaV_2C;

% calcolo coe delle tre traiettorie
coe_leo=[r_leo ; 0 ; 28.5*pi/180 ; NaN ; RAAN; NaN]; 
coe_gto=[a_gto ; e_gto ; 28.5*pi/180 ; pi ; RAAN; 0];
coe_geo=[r_geo ; 0 ; 0 ; NaN ; RAAN; NaN];


% ESERCITAZIONE 2
% PARTE 2

% Nuova strategia (D)
deltaPhi_1D=linspace(0,2.85*pi/180,100);
deltaPhi_2D=deltaPhi-deltaPhi_1D;
deltaV_1D=sqrt(v_leo^2+v_gto_per^2-2*v_leo*v_gto_per*cos(deltaPhi_1D));
deltaV_2D=sqrt(v_gto_apo^2+v_geo^2-2*v_gto_apo*v_geo*cos(deltaPhi_2D));
deltaV_totD=deltaV_1D+deltaV_2D;
index=find(deltaV_totD==min(deltaV_totD));

% disegno andamento deltaV_totD in funzione di deltaPhi_1
figure(1)
hold on
plot(deltaPhi_1D*180/pi,deltaV_totD);
plot(deltaPhi_1D(index)*180/pi,deltaV_totD(index),'r*');
title('deltaV(deltaPhi_1)')
xlabel('deltaPhi_1 (deg)')
ylabel('deltaV (km/s)')
hold off
grid on

% trovo i valori voluti
deltaPhi_1D=deltaPhi_1D(index);
deltaPhi_2D=deltaPhi_2D(index);
deltaV_1D=deltaV_1D(index);
deltaV_2D=deltaV_2D(index);
deltaV_totD=deltaV_totD(index);


% ESERCITAZIONE 2
% PARTE 3

T_wo=3/4*T_siderale-t_volo;
a_wo=((T_wo/(2*pi))^2*mu)^(1/3);
v_wo_apo=sqrt(2*mu/r_geo-mu/a_wo);
e_wo=r_geo/a_wo-1;
r_wo_per=a_wo*(1-e_wo);

% uso il teorema dei seni per calcolare nuovi deltaPhi
supp1=asin(v_gto_apo/deltaV_2C*sin(deltaPhi));
supp2=pi-asin(v_geo/v_wo_apo*sin(supp1));
deltaPhi_2E=pi-supp1-supp2;
deltaPhi_1E=deltaPhi-deltaPhi_2E;

% nuova strategia (E)
deltaV_1E=v_gto_per-v_leo;
deltaV_2E=sqrt(v_gto_apo^2+v_wo_apo^2-2*v_gto_apo*v_wo_apo*cos(deltaPhi_1E));
deltaV_3E=sqrt(v_wo_apo^2+v_geo^2-2*v_wo_apo*v_geo*cos(deltaPhi_2E));
deltaV_totE=deltaV_1E+deltaV_2E+deltaV_3E;

% uso il teorema dei seni per calcolare nuovi deltaPhi
deltaPhi_1F=deltaPhi_1D;
supp3=asin(v_gto_apo/deltaV_2D*sin(deltaPhi_2D));
supp4=pi-asin(v_geo/v_wo_apo*sin(supp3));
deltaPhi_3F=pi-supp3-supp4;
deltaPhi_2F=deltaPhi_2D-deltaPhi_3F;

% nuova strategia (F)
deltaV_1F=deltaV_1D;
deltaV_2F=sqrt(v_gto_apo^2+v_wo_apo^2-2*v_gto_apo*v_wo_apo*cos(deltaPhi_2F));
deltaV_3F=sqrt(v_wo_apo^2+v_geo^2-2*v_wo_apo*v_geo*cos(deltaPhi_3F));
deltaV_totF=deltaV_1F+deltaV_2F+deltaV_3F;


% ESERCITAZIONE 2
% PARTE 4

% opzionale