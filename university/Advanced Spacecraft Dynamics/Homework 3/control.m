function dw=control(t,w)
% CONTROL contains the feeback loop (control algorithm + actuator dynamics
% + kinematics) provided for attitude tracking using an array of RWs

global in r0 v0 q0e_0 gainA_inv gainB Jc A mu

q0e=w(1);
qe=w(2:4);
omega=w(5:7);
omega_s=w(8:11);

dw=zeros(11,1);


% orbit propagation
[r,v]=CorePropagator(r0,v0,t,mu);
v_r=dot(v,r/norm(r));
h=cross(r,v);

% rotation matrix from C to B
R_CtoB=(q0e^2-qe'*qe)*eye(3)+2*(qe*qe')-2*q0e*skew(qe);

% desired angular velocity and acceleration
anu_dot=norm(h)/norm(r)^2;
omega_C=[0 ; -anu_dot ; 0];                     % in the frame C
anu_ddot=-2*norm(h)*v_r/norm(r)^3;
omega_c_dot=[0 ; -anu_ddot ; 0];                % in the frame C

% error angular velocity
omega_d=omega-R_CtoB*omega_C;

qe_til=skew(qe);
omega_til=skew(omega);

% commanded torque
Mc=[0;0;0];
Tc=omega_til*Jc*omega-Mc+Jc*R_CtoB*omega_c_dot-Jc*gainA_inv*gainB*omega_d...
   -sign(q0e_0)*Jc*gainA_inv*qe;

% actual torque
switch in
    case 1
        Ta=Tc;
    case 2
        Ta=-omega_til*A*omega_s+Tc;
end

% odes
dw(1)=-0.5*omega_d'*qe;
dw(2:4)=0.5*(q0e*eye(3)+qe_til)*omega_d;
dw(5:7)=Jc\(-omega_til*Jc*omega+Ta+Mc);
dw(8:11)=-A'*((A*A')\Tc);


end