function dw=int_cos0(t,w)
global k1 k2 k3 M1 M2 M3

dw=zeros(9,1);

dw(1)=k1*w(2)*w(3)+M1;
dw(2)=k2*w(3)*w(1)+M2;
dw(3)=k3*w(1)*w(2)+M3;
dw(4)=w(3)*w(5)-w(2)*w(6); %c11
dw(5)=w(1)*w(6)-w(3)*w(4); %c21
dw(6)=w(2)*w(4)-w(1)*w(5); %c31
dw(7)=w(3)*w(8)-w(2)*w(9); %c12
dw(8)=w(1)*w(9)-w(3)*w(7); %c22
dw(9)=w(2)*w(7)-w(1)*w(8); %c32
