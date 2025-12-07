function dw=dyn_and_kin(t,w)

r=w(1);
vr=w(3);
gamma_r=w(4);

dw=zeros(4,1);

global m w_e mu Cd S R_e ro_0 beta

h=r-R_e;                        % in m
ro=ro_0*exp(-beta*h);           % in kg/m^3

D=0.5*Cd*ro*S*vr^2;             % in N
a_D=D/m;                        % in m/s^2

dw(1)=vr*sin(gamma_r);                                                      % in m/s
dw(2)=vr*cos(gamma_r)/r;                                                    % in rad/s
dw(3)=(-mu/r^2)*sin(gamma_r)-a_D+w_e^2*r*sin(gamma_r);                      % in m/s^2
dw(4)=(vr/r-mu/(r^2*vr))*cos(gamma_r)+(w_e^2*r/vr)*cos(gamma_r)+2*w_e;      % in rad/s

end

