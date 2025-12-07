function[r1,v1]=CorePropagator(r0,v0,dt,mu)
% CorePropagator
%
% da r0, v0, dt a r1, v1

% ricavo elementi orbitali classici
coe=rvECI2coe(r0,v0,mu);

ecc=coe(2);
anu0=coe(6);

% calcolo anomalia eccentrica al tempo t0
anu0=wrapTo2Pi(anu0);
if anu0==2*pi
    anu0=0;
end

supp1=anu0/2;
if anu0>=-pi/2 || anu0<=pi/2
    AE0=2*atan(sqrt((1-ecc)/(1+ecc))*tan(supp1));
else
    AE0=2*atan(sqrt((1-ecc)/(1+ecc))*tan(supp1))+pi;
end

% calcolo anomalia media al tempo t0
XM0=AE0-ecc*sin(AE0);

% calcola anomalia media al tempo t1
XM1=XM0+sqrt(mu/coe(1)^3)*dt;

% calcolo anomalia eccentrica al tempo t1
eps=10^-7;
err=eps+1;
supp2=XM1;
while err>eps
    AE1=supp2-(supp2-ecc*sin(supp2)-XM1)/(1-ecc*cos(supp2));
    err=abs(supp2-AE1);
    supp2=AE1;
end

% calcolo anomalia vera al tempo t1
anu1=2*atan(sqrt((1+ecc)/(1-ecc))*tan(AE1/2));

% calcolo posizione e velocità al tempo t1
coe(6)=anu1;

[r1,v1]=coe2rvECI(coe,mu); 
end