function dw=dyn_kin_and_adjoints(t,w)

global mu c u_T_max

r=w(1);
csi=w(2);
vr=w(3);
vt=w(4);
lambda_1=w(5);
lambda_3=w(6);
lambda_4=w(7);

dw=zeros(7,1);

% estimation of optimal control
u1=u_T_max;
s_u2=-lambda_3/sqrt(lambda_3^2+lambda_4^2);
c_u2=-lambda_4/sqrt(lambda_3^2+lambda_4^2);

% estimation of thrust acceleration
a_Tr=s_u2*u1*c/(c-u1*t);
a_Tt=c_u2*u1*c/(c-u1*t);

% dynamics and kinematics equations (state equations)
dw(1)=vr;
dw(2)=vt/r;
dw(3)=-mu/r^2+vt^2/r+a_Tr;
dw(4)=-vr*vt/r+a_Tt;

% adjoint equations (costate equations)
dw(5)=lambda_3*(vt^2/r^2-2*mu/r^3)-lambda_4*vr*vt/r^2;
dw(6)=-lambda_1+lambda_4*vt/r;
dw(7)=-2*lambda_3*vt/r+lambda_4*vr/r;


end

