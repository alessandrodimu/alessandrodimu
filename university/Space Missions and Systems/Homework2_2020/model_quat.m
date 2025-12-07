function dw=model_quat(t,w,omega_est)
% Purpose: Integrate quaternion [q1 ; q2 ; q3 ; q4].
% Input: w = [q1 q2 q3 q4]
% Output dw = [dq1 dq2 dq3 dq4]

dw=zeros(4,1);

quat=w(1:4);
OMEGA=[-skew(omega_est) omega_est ; -omega_est' 0];
dw(1:4)=0.5*OMEGA*quat;

end