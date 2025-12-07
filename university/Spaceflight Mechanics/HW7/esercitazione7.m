% ESERCITAZIONE 7
% PARTE 1

% ha qualche problema nella rappresentazione dell'erpolodia nel piano e nel
% tracciamento della polodia nel caso dello smorzamento
% DA RIVEDERE

% ipotesi considerate
% - J1 > J2 > J3
% - il vettore Hc ha norma 1
% - Emin = 1

% valori iniziali consigliati (caso A)
Emin_A=1;
Emax_A=2;
E2_A=Emin_A+0.2*(Emax_A-Emin_A);
J1=0.5/Emin_A;
J2=J1/E2_A;
J3=J1*Emin_A/Emax_A;
J_A=[J1; J2; J3];
Er_A=Emin_A+0.4*(Emax_A-Emin_A);

% altri valori iniziali (caso B)
Emin_B=1;
Emax_B=3;
E2_B=Emin_B+0.4*(Emax_B-Emin_B);
J1=0.5/Emin_B;
J2=J1/E2_B;
J3=J1*Emin_B/Emax_B;
J_B=[J1; J2; J3];
Er_B=Emin_B+0.7*(Emax_B-Emin_B);


% calcolo w1 e w3 assumendo w2 = 0
w1=sqrt((1-2*J_A(3)*Er_A)/(J_A(1)*(J_A(1)-J_B(3))));
w2=0;
w3=sqrt((1-2*J_A(1)*Er_A)/(J_A(3)*(J_A(3)-J_B(1))));
W_A=[w1; w2; w3];
Hc_A=[J_A(1)*w1; J_A(2)*w2; J_A(3)*w3];

w1=sqrt((1-2*J_B(3)*Er_B)/(J_B(1)*(J_B(1)-J_B(3))));
w2=0;
w3=sqrt((1-2*J_B(1)*Er_B)/(J_B(3)*(J_B(3)-J_B(1))));
W_B=[w1; w2; w3];
Hc_B=[J_B(1)*w1; J_B(2)*w2; J_B(3)*w3];


% calcolo k1, k2, k3, M1, M2, M3
global k1 k2 k3 M1 M2 M3
k1=(J_A(2)-J_A(3))/J_A(1);
k2=(J_A(3)-J_A(1))/J_A(2);
k3=(J_A(1)-J_A(2))/J_A(3);

M1=0;
M2=0;
M3=0;

ellissoide(J_A,W_A,Hc_A,Er_A,[1 2]);

% disegno la poloide sull'ellissoide di costante momento (per il caso A)
vect_Er=[linspace(Emin_A,Emax_A,20) E2_A];
vect_Er=sort(vect_Er);
for ind_1=1:length(vect_Er)
    w1=sqrt((1-2*J_A(3)*vect_Er(ind_1))/(J_A(1)*(J_A(1)-J_A(3))));
    w2=0;
    w3=sqrt((1-2*J_A(1)*vect_Er(ind_1))/(J_A(3)*(J_A(3)-J_A(1))));
    
    W_pol=[w1; w2; w3];
    Hc_pol=[J_A(1)*w1; J_A(2)*w2; J_A(3)*w3];
    Erp=vect_Er(ind_1);
    
    polhode(J_A,W_pol,Hc_pol,Erp,3);
end
title('Poloidi per caso A');



% ESERCITAZIONE 7
% PARTE 2

% disegno l'erpoloide sul piano normale ad Hc nel sistema di riferimento
% inerziale e l'ellissoide di costante energia con la rispettiva poloide
herpolhode(J_A,W_A,Hc_A,Er_A,4);
title('Erpolodia per caso A');



% ESERCITAZIONE 7
% PARTE 3

% smorzamento con il modello del pozzo di energia (applicando una coppia)
% riporto il percorso di omega sull'ellissoide
% caso A --> uso i valori iniziali consigliati 
% caso B --> uso altri valori iniziali
% disegno andamenti di Er, Hc, omega
moto_frenato(J_A,W_A,Hc_A,Er_A,[5 6 7 8]);
figure(5)
title('Andamento del vettore omega per dissipazione (A)')
figure(6)
title('Andamento di w (modulo e componenti) (A)');
figure(7)
title('Andamento di H (modulo e componenti) (A)');
figure(8)
title('Andamento di Er (A)');

moto_frenato(J_B,W_B,Hc_B,Er_B,[9 10 11 12]);
figure(9)
title('Andamento del vettore omega per dissipazione (B)')
figure(10)
title('Andamento di w (modulo e componenti) (B)');
figure(11)
title('Andamento di H (modulo e componenti) (B)');
figure(12)
title('Andamento di Er (B)');

