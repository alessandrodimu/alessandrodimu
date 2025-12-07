function dw=hamiltonian_eqs(t,w)

global mu beta

r=w(1);
theta=w(2);
vr=w(3);
vtheta=w(4);
K=w(5);
pr=w(6);
pvr=w(7);
pvtheta=w(8);

dw=zeros(5+3,1);

alpha=atan((-3*pvr-sqrt(9*pvr^2+8*pvtheta^2))/(4*pvtheta));


% state equations
dw(1)=vr;
dw(2)=vtheta/r;
dw(3)=vtheta^2/r-mu/r^2+beta/r^2*cos(alpha)^3;
dw(4)=-vr*vtheta/r+beta/r^2*cos(alpha)^2*sin(alpha);
dw(5)=-1;

% costate equations
dw(6)=pvr*(vtheta^2/r^2-2*mu/r^3+2*beta/r^3*cos(alpha)^3)+...
    pvtheta*(-vr*vtheta/r^2+2*beta/r^3*cos(alpha)^2*sin(alpha));
dw(7)=-pr+pvtheta*vtheta/r;
dw(8)=-2*pvr*vtheta/r+pvtheta*vr/r;

end

