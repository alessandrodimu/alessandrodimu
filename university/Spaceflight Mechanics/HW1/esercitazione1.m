% ESERCITAZIONE 1
% PARTE 1

% calcolo il vettore tspan: disegno un pto sulla traiettoria ogni 10 sec
tspan=linspace(0,24*3600,8641);

% dati iniziali r0, v0 per ogni orbita
r0_GPS=[13426.726354 ; -20645.023752 ; 10108.071628];
v0_GPS=[2.638622 ; 0.329454 ; -2.805018];

r0_Molniya=[2869.944136 ; -1656.963019 ; -6617.757380];
v0_Molniya=[4.814217 ; 8.338469 ; 0.000000];

r0_Geosincrona=[-29260.792529 ; 23865.119693 ; 23865.119693];
v0_Geosincrona=[-2.334835 ; -1.212853 ; -1.212853];

% calcolo matrici di posizione e velocita e disegno la traiettoria per ogni orbita
[rMatrixECI_GPS,vMatrixECI_GPS]=ShellPropagator(r0_GPS,v0_GPS,tspan);
DisegnaTerra(1)
DisegnaTraj3D(1,rMatrixECI_GPS,'b-')
title('Traiettoria orbita GPS')

[rMatrixECI_Molniya,vMatrixECI_Molniya]=ShellPropagator(r0_Molniya,v0_Molniya,tspan);
DisegnaTerra(2)
DisegnaTraj3D(2,rMatrixECI_Molniya,'b-')
title('Traiettoria orbita Molniya')

[rMatrixECI_Geosincrona,vMatrixECI_Geosincrona]=ShellPropagator(r0_Geosincrona,v0_Geosincrona,tspan);
DisegnaTerra(3)
DisegnaTraj3D(3,rMatrixECI_Geosincrona,'b-')
title('Traiettoria orbita Geosincrona')



% ESERCITAZIONE 1
% PARTE 2

% 01/01/1995 @ 00:00:00 --> MJD = 49718
% 15/01/1995 @ 00:00:00 --> MJD = 49732
% 16/01/1995 @ 00:00:00 --> MJD = 49733
T_siderale=86164;
T_solare=86400;
tspan_MJD=linspace(49732*T_solare,49733*T_solare,8641);

% disegno la traccia a terra per ogni orbita
load coastlines coastlon coastlat

figure(4)
plot(coastlon,coastlat);
axis([-180 180 -90 90]);
[LatVect_GPS,LongVect_GPS]=ShellTrace(rMatrixECI_GPS,tspan_MJD);
DisegnaTraccia2D(4,LatVect_GPS,LongVect_GPS,'ro')
title('Traccia a terra orbita GPS')

figure(5)
plot(coastlon,coastlat);
axis([-180 180 -90 90]);
[LatVect_Molniya,LongVect_Molniya]=ShellTrace(rMatrixECI_Molniya,tspan_MJD);
DisegnaTraccia2D(5,LatVect_Molniya,LongVect_Molniya,'ro')
title('Traccia a terra orbita Molniya')

figure(6)
plot(coastlon,coastlat);
axis([-180 180 -90 90]);
[LatVect_Geosincrona,LongVect_Geosincrona]=ShellTrace(rMatrixECI_Geosincrona,tspan_MJD);
DisegnaTraccia2D(6,LatVect_Geosincrona,LongVect_Geosincrona,'ro')
title('Traccia a terra orbita Geosincrona')



% ESERCITAZIONE 1
% PARTE 3

mu=398600;
coe=rvECI2coe(r0_Molniya,v0_Molniya,mu);

% u3=[u1;u2] anu3=[anu1;anu2] AE3=[AE1;AE2] XM3=[XM1;XM2]
u3=[asin(sin(45*pi/180)/sin(coe(3)));pi-asin(sin(45*pi/180)/sin(coe(3)))];
anu3=[u3(1)-coe(5);u3(2)-coe(5)];
AE3=[anu2AE(coe(2),anu3(1));anu2AE(coe(2),anu3(2))];
XM3=[AE2XM(coe(2),AE3(1));AE2XM(coe(2),AE3(2))];
tempo_12=sqrt(coe(1)^3/mu)*(XM3(2)-XM3(1));

% se calcolo il tempo in un arco di 24 ore (il tempo di due periodi orb)
delta_t=2*tempo_12;



% ESERCITAZIONE 1
% PARTE 4

% calcolo il vettore tspan: disegno un pto sulla traiettoria ogni 10 sec
tspan4=linspace((49732-1/6)*T_solare,(49732+1/6)*T_solare,2881);

% calcolo i vari coe delle 3 orbite Molniya
[coe4_1,coe4_2,coe4_3]=aiuto_es1_4(1200);

% disegno le traiettorie in 3D delle tre orbite Molniya
DisegnaTerra(7)

[rVectECI_Molniya4_1,vVectECI_Molniya4_1]=coe2rvECI(coe4_1,mu);
[rMatrixECI_Molniya4_1,vMatrixECI_Molniya4_1]=ShellPropagator(rVectECI_Molniya4_1,vVectECI_Molniya4_1,tspan4);
DisegnaTraj3D(7,rMatrixECI_Molniya4_1,'r-')

[rVectECI_Molniya4_2,vVectECI_Molniya4_2]=coe2rvECI(coe4_2,mu);
[rMatrixECI_Molniya4_2,vMatrixECI_Molniya4_2]=ShellPropagator(rVectECI_Molniya4_2,vVectECI_Molniya4_2,tspan4);
DisegnaTraj3D(7,rMatrixECI_Molniya4_2,'b-')

[rVectECI_Molniya4_3,vVectECI_Molniya4_3]=coe2rvECI(coe4_3,mu);
[rMatrixECI_Molniya4_3,vMatrixECI_Molniya4_3]=ShellPropagator(rVectECI_Molniya4_3,vVectECI_Molniya4_3,tspan4);
DisegnaTraj3D(7,rMatrixECI_Molniya4_3,'g-')

title('Traiettoria costellazione Molniya')


% disegno le tracce a terra delle 3 orbite Molniya
figure(8)
plot(coastlon,coastlat);
axis([-180 180 -90 90]);

[LatVect_Molniya4_1,LongVect_Molniya4_1]=ShellTrace(rMatrixECI_Molniya4_1,tspan4);
DisegnaTraccia2D(8,LatVect_Molniya4_1,LongVect_Molniya4_1,'ro')

[LatVect_Molniya4_2,LongVect_Molniya4_2]=ShellTrace(rMatrixECI_Molniya4_2,tspan4);
DisegnaTraccia2D(8,LatVect_Molniya4_2,LongVect_Molniya4_2,'bo')

[LatVect_Molniya4_3,LongVect_Molniya4_3]=ShellTrace(rMatrixECI_Molniya4_3,tspan4);
DisegnaTraccia2D(8,LatVect_Molniya4_3,LongVect_Molniya4_3,'go')

title('Traccia a terra costellazione Molniya')
