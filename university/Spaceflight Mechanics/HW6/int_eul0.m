function dw=int_eul0(t,w)
global k1 k2 k3 M1 M2 M3

dw=zeros(6,1);

dw(1)=k1*w(2)*w(3)+M1;
dw(2)=k2*w(3)*w(1)+M2;
dw(3)=k3*w(1)*w(2)+M3;
dw(4)=1/cos(w(5))*(cos(w(5))*w(1)+sin(w(4))*sin(w(5))*w(2)+cos(w(4))*sin(w(5))*w(3));     %phi
dw(5)=cos(w(4))*w(2)-sin(w(4))*w(3);                                                      %theta
dw(6)=1/cos(w(5))*(sin(w(4))*w(2)+cos(w(4))*w(3));                                        %psi

