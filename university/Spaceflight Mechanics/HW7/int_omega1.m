function dw=int_omega1(t,w)
global k1 k2 k3 J Er eps 

dw=zeros(3,1);

Hc=[J(1)*w(1) ; J(2)*w(2) ; J(3)*w(3)];
W=[w(1) ; w(2) ; w(3)];
vers_M=Hc/norm(Hc);
vers_K=cross(Hc,W)/norm(cross(Hc,W));
vers_L=cross(vers_K,vers_M);

% perche c'e quel meno??
Emin=1;
T=-eps*(Er-Emin)*exp(-eps*t)/dot(vers_L,W)*vers_L;

dw(1)=k1*w(2)*w(3)+T(1)/J(1);
dw(2)=k2*w(3)*w(1)+T(2)/J(2);
dw(3)=k3*w(1)*w(2)+T(3)/J(3);

