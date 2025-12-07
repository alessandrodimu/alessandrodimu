function dw=int_omega(t,w)
global k1 k2 k3 M1 M2 M3

dw=zeros(3,1);

dw(1)=k1*w(2)*w(3)+M1;
dw(2)=k2*w(3)*w(1)+M2;
dw(3)=k3*w(1)*w(2)+M3;

