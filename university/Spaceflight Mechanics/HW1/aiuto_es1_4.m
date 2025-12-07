function[coe1,coe2,coe3]=aiuto_es1_4(h)
% es1_4_calcolo_molniya
% serve solo come supporto all'esercitazione 1 per semplificare lo script

T_solare=86400;
T_siderale=86164;
mu=398600;

a=(mu*(T_siderale/2)^2/(4*pi^2))^(1/3);
rp=6378+h;
e=1-rp/a;
i=63.4*pi/180;
RAAN_greenwich_0=100*pi/180;
MJD_0=49718*T_solare;
MJD=49732*T_solare;
w_rot_terra=2*pi/T_siderale;
RAAN_greenwich=RAAN_greenwich_0+w_rot_terra*(MJD-MJD_0);
LONG_mosca=37.618*pi/180;
GOM1=RAAN_greenwich+LONG_mosca-pi/2;

GOM2=GOM1-2*pi/3;
GOM3=GOM1+2*pi/3;
POM=270*pi/180;

XM1=pi-4*3600*sqrt(mu/a^3);
XM2=XM1-2*pi/3;
XM3=XM1+2*pi/3;
anu1=AE2anu(e,XM2AE(e,XM1));
anu2=AE2anu(e,XM2AE(e,XM2));
anu3=AE2anu(e,XM2AE(e,XM3));


% ho trovato coe1 al tempo t_A-4hrs = 14/01/1995 @ 20:00:00
coe1=[a;e;i;GOM1;POM;anu1];

% ho trovato coe2 al tempo t_A-4hrs = 14/01/1995 @ 20:00:00
coe2=[a;e;i;GOM2;POM;anu2];

% ho trovato coe3 al tempo t_A-4hrs = 14/01/1995 @ 20:00:00
coe3=[a;e;i;GOM3;POM;anu3];

end