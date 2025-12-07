function dw=model_prop(t,w,omega)
% Purpose: Integrate quaternion [q1 ; q2 ; q3 ; q4].
% Input: w = [q1 q2 q3 q4]
% Output dw = [dq1 dq2 dq3 dq4]

dw=zeros(4+9,1);

quat=w(1:4);
OMEGA=[-skew(omega) omega ; -omega' 0];
dw(1:4)=0.5*OMEGA*quat;

phi=w(5:13);
PHI=reshape(phi,3,3);
dphi=-skew(omega)*PHI;
dphi=reshape(dphi,1,9);
dw(5:13)=dphi;

end