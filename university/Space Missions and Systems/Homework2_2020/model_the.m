function dw=model_the(t,w,omega_est)
% Purpose: Integrate matrix THETA.
% Input: w = [T11 T12 T13 T21 T22 T23 T31 T32 T33]
% Output dw = [dT11 dT12 dT13 dT21 dT22 dT23 dT31 dT32 dT33]

dw=zeros(9,1);

theta=w(1:9);
THETA=reshape(theta,3,3);
dTHETA=-skew(omega_est)*THETA;
dtheta=reshape(dTHETA,1,9);

dw(1:9)=dtheta;

end

