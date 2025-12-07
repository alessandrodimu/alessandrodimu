close all
clear all
clc
global k1 k2 k3 M1 M2 M3

disp('Si assume unitario il modulo del momento angolare');
disp('e l''energia cinetica minima (J1=0.5)');
norma_Hc=1;
Emin=1;
J1=0.5*norma_Hc/Emin;
disp('Inserire l''energia cinetica massima');
Emax=input('Emax(>1 - consigliato 2):');
J3=J1*Emin/Emax;
disp('Definire il valore di J2 inserendo un numero compreso tra i due estremi');
disp('0: E2 = Emin (corpo assialsimmetrico prolato)')
disp('1: E2 = Emax (corpo assialsimmetrico oblato)')
fraz=input('(consigliato .2):');
E2=Emin+fraz*(Emax-Emin);
J2=J1/E2;
J=[J1; J2; J3]; 
disp('Inserendo l''energia del moto con un numero compreso tra i due estremi');
disp('0: Er = Emin (rotazione intorno a b1)')
disp('1: Er = Emax (rotazione intorno a b3)')
fraz=input('(consigliato .4):');
Er=Emin+fraz*(Emax-Emin);
vett_Er=[linspace(Emin,Emax,14) E2];
vett_Er=sort(vett_Er);

% calcolo k1, k2, k3, M1, M2, M3
k1=(J2-J3)/J1;
k2=(J3-J1)/J2;
k3=(J1-J2)/J3;

M1=0;
M2=0;
M3=0;

% calcolo w1 e w3 assumendo w2 = 0
w1=sqrt((1-2*J3*Er)/(J1*(J1-J3)));
w2=0;
w3=sqrt((1-2*J1*Er)/(J3*(J3-J1)));

W=[w1; w2; w3];
Hc=[J1*w1; J2*w2; w3*J3];
ellissoide(J,W,Hc,Er,[1 2]);

% disegno la poloide sull'ellissoide di costante momento
vett_W=zeros(3,length(vett_Er));
for ind_2=1:length(vett_Er)
    w1p=sqrt((1-2*J3*vett_Er(ind_2))/(J1*(J1-J3)));
    w2p=0;
    w3p=sqrt((1-2*J1*vett_Er(ind_2))/(J3*(J3-J1)));
    Wp=[w1p ; w2p ; w3p];
    Hcp=[J1*w1p; J2*w2p; J3*w3p];
    Erp=vett_Er(ind_2);
    vett_W(:,ind_2)=Wp;
end
polhode(J,vett_W,Hcp,Erp,3);


% disegno l'erpoloide sul piano normale ad Hc nel sistema di riferimento
% inerziale e l'ellissoide di costante energia con la rispettiva poloide
herpolhode(J,W,Hc,Er,4)


% studia il moto smorzato con il modello del pozzo di energia
% aggiungendo un'opportuna coppia normale a Hc
moto_frenato(J,W,Hc,Er,[5 6 7 8]);
