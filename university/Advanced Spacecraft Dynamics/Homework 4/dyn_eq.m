function dw = dyn_eq(t,w)
% DYN_EQ contains the dynamics equations for the equinoctial elements that
% describe the motion of the spacecraft considering also the perturbation
% effects due to J2

x=w(1:7);

dw=zeros(7,1);

global mu R_e J2 DU K p_d e_d i_d u_T_max c


% estimate G matrix
eta=1+x(2)*cos(x(6))+x(3)*sin(x(6));
G=sqrt(x(1)/mu)*[0 2*x(1)/eta 0;...
   sin(x(6)) ((eta+1)*cos(x(6))+x(2))/eta -x(3)*(x(4)*sin(x(6))-x(5)*cos(x(6)))/eta;...
   -cos(x(6)) ((eta+1)*sin(x(6))+x(3))/eta x(2)*(x(4)*sin(x(6))-x(5)*cos(x(6)))/eta;...
   0 0 cos(x(6))*(1+x(4)^2+x(5)^2)/(2*eta);...
   0 0 sin(x(6))*(1+x(4)^2+x(5)^2)/(2*eta)];


% estimate the perturbation accelerations due to J2
r=x(1)/eta;

stheta=sin(x(6))*x(4)/sqrt(x(4)^2+x(5)^2)-cos(x(6))*x(5)/sqrt(x(4)^2+x(5)^2);
ctheta=cos(x(6))*x(4)/sqrt(x(4)^2+x(5)^2)+sin(x(6))*x(5)/sqrt(x(4)^2+x(5)^2);
si=2*sqrt(x(4)^2+x(5)^2)/(1+x(4)^2+x(5)^2);
ci=(1-(x(4)^2+x(5)^2))/(1+x(4)^2+x(5)^2);

a_P_r=(3*mu*R_e^2*J2/r^4)*(3*stheta^2*si^2-1)/2;
a_P_theta=-(3*mu*R_e^2*J2/r^4)*si^2*stheta*ctheta;
a_P_h=-(3*mu*R_e^2*J2/r^4)*si*ci*stheta;
a_P=[a_P_r ; a_P_theta ; a_P_h]./(DU^2);


% estimate the constraint vector
psi=[x(1)-p_d;x(2)^2+x(3)^2-e_d^2;x(4)^2+x(5)^2-tan(i_d/2)^2];
dpsi=[1 0 0 0 0;...
      0 2*x(2) 2*x(3) 0 0;...
      0 0 0 2*x(4) 2*x(5)];


% estimate the control vector
b=G'*dpsi'*K*psi;
u_T=-u_T_max*(x(7)*(b+a_P))/(max(u_T_max,norm(x(7)*(b+a_P))));


% estimate the thrust acceleration
a_T=u_T/x(7);
a=a_P+a_T;


% integration (both dynamical eq + lyapunov)
dw(1:5)=G*a;
dw(6)=sqrt(mu/x(1)^3)*eta^2+...
      sqrt(x(1)/mu)*(x(4)*sin(x(6))-x(5)*cos(x(6)))/eta*a(3);
dw(7)=-norm(u_T)/c;

