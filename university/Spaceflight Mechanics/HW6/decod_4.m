function [q0,q1,q2,q3]=decod_4(BN)

[E,phi]= decod_2(BN);
q0=cos(phi/2);
q1=E(1)*sin(phi/2);
q2=E(2)*sin(phi/2);
q3=E(3)*sin(phi/2);

end

