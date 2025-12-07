function dw=state_eqs(t,w)

global mu beta

r=w(1);
theta=w(2);
vr=w(3);
vtheta=w(4);
alpha=w(5);

dw=zeros(5,1);


% state equations
dw(1)=vr;
dw(2)=vtheta/r;
dw(3)=vtheta^2/r-mu/r^2+beta/r^2*cos(alpha)^3;
dw(4)=-vr*vtheta/r+beta/r^2*cos(alpha)^2*sin(alpha);
dw(5)=0;

end

