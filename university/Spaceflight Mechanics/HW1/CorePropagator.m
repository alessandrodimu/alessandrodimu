function[r1,v1]=CorePropagator(r0,v0,dt)
% CorePropagator
%
% da r0, v0, dt a r1, v1

mu=398600;

coe=rvECI2coe(r0,v0,mu);         % ricavo elementi orbitali classici
AE0=anu2AE(coe(2),coe(6));       % calcolo anomalia eccentrica al tempo t0
XM0=AE2XM(coe(2),AE0);           % calcolo anomalia media al tempo t0
XM1=XM0+sqrt(mu/coe(1)^3)*dt;    % calcola anomalia media al tempo t1
AE1=XM2AE(coe(2),XM1);           % calcolo anomalia eccentrica al tempo t1
anu1=AE2anu(coe(2),AE1);         % calcolo anomalia vera al tempo t1
coe(6)=anu1;                     % aggiorno gli elementi orbitali
[r1,v1]=coe2rvECI(coe,mu);       % calcolo posizione e velocità al tempo t1

end