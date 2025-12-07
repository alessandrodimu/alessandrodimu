function dw=model_psi(t,w,omega_est)
% Purpose: Integrate matrix PSI.
% Input: w = [P11 P12 P13 P21 P22 P23 P31 P32 P33]
% Output dw = [dP11 dP12 dP13 dP21 dP22 dP23 dP31 dP32 dP33]

dw=zeros(9,1);

psi=w(1:9);
PSI=reshape(psi,3,3);
dPSI=-skew(omega_est)*PSI-0.5*eye(3);
dpsi=reshape(dPSI,1,9);

dw(1:9)=dpsi;

end