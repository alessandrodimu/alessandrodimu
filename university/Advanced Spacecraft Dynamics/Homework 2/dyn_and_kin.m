function dw=dyn_and_kin(t,w)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

omega_1=w(1:3);
sigma=w(4:7);
v_P1=w(8:10);
q01=w(11);
q1=w(12:14);
theta=w(15:18);
r_P1=w(19:21);

dw=zeros(10+11,1);

global m_b1 m_bb J_b1 J_bb k1 k2 N l R mu

% A) evaluate the joint partials Gamma_Gk and Gamma_Gk_dot
Gamma_G1=[1;0;0];
Gamma_G2=[1;0;0];
Gamma_G3=[1;0;0];
Gamma_G4=[1;0;0];
% Gamma_G1_dot = [0;0;0] and also for all the other joints


% B) for each B_i retrieve from the variables u and x the following:
% omega_1, v_P1, R_NtoB1 and r_P1

% omega_1, sigma, v_P1 and r_P1 were already calculated
R_Nto1=(q01^2-q1'*q1)*eye(3)+2*(q1*q1')-2*q01*skew(q1);


% C) find all relative rotation matrices
B1=(R_Nto1*N')';
B2=([ 1 0 0 ;...
       0 cos(pi+theta(1)) sin(pi+theta(1)) ;...
       0 -sin(pi+theta(1)) cos(pi+theta(1)) ]*B1')';
B3=([ 1 0 0 ;...
       0 cos(theta(2)) sin(theta(2)) ;...
       0 -sin(theta(2)) cos(theta(2)) ]*B1')';
B4=([ 1 0 0 ;...
       0 cos(pi/2+theta(3)) sin(pi/2+theta(3)) ;...
       0 -sin(pi/2+theta(3)) cos(pi/2+theta(3)) ]*B1')';
B5=([ 1 0 0 ;...
       0 cos(3*pi/2+theta(4)) sin(3*pi/2+theta(4)) ;...
       0 -sin(3*pi/2+theta(4)) cos(3*pi/2+theta(4)) ]*B1')';


% D) find the path vectors to insert in matrices OM and V
r_2toG1=N'*B2*[ 0 ; -l/2 ; 0 ];
r_3toG2=N'*B3*[ 0 ; -l/2 ; 0 ];
r_4toG3=N'*B4*[ 0 ; -l/2 ; 0 ];
r_5toG4=N'*B5*[ 0 ; -l/2 ; 0 ];

r_1toG1=N'*B1*[ 0 ; -R ; 0 ];
r_1toG2=N'*B1*[ 0 ; R ; 0 ];
r_1toG3=N'*B1*[ 0 ; 0 ; R ];
r_1toG4=N'*B1*[ 0 ; 0 ; -R ];

r_2to1=r_2toG1-r_1toG1;
r_3to1=r_3toG2-r_1toG2;
r_4to1=r_4toG3-r_1toG3;
r_5to1=r_5toG4-r_1toG4;

% r_P2=r_P1-r_2to1;
% r_P3=r_P1-r_3to1;
% r_P4=r_P1-r_4to1;
% r_P5=r_P1-r_5to1;

r21_til=skew(r_2to1);
r31_til=skew(r_3to1);
r41_til=skew(r_4to1);
r51_til=skew(r_5to1);
r2G1_til=skew(r_2toG1);
r3G2_til=skew(r_3toG2);
r4G3_til=skew(r_4toG3);
r5G4_til=skew(r_5toG4);


% E) evaluate OM and V and find omega_i and v_Pi for all the bodies
OM=[eye(3) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3) ;...
    B2'*B1 Gamma_G1 zeros(3,1) zeros(3,1) zeros(3,1) zeros(3) ;...
    B3'*B1 zeros(3,1) Gamma_G2 zeros(3,1) zeros(3,1) zeros(3) ;...
    B4'*B1 zeros(3,1) zeros(3,1) Gamma_G3 zeros(3,1) zeros(3) ;...
    B5'*B1 zeros(3,1) zeros(3,1) zeros(3,1) Gamma_G4 zeros(3) ];

V=[zeros(3) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) eye(3) ;...
   r21_til*N'*B1 r2G1_til*N'*B2*Gamma_G1 zeros(3,1) zeros(3,1) zeros(3,1) eye(3) ;...
   r31_til*N'*B1 zeros(3,1) r3G2_til*N'*B3*Gamma_G2 zeros(3,1) zeros(3,1) eye(3) ;...
   r41_til*N'*B1 zeros(3,1) zeros(3,1) r4G3_til*N'*B4*Gamma_G3 zeros(3,1) eye(3) ;...
   r51_til*N'*B1 zeros(3,1) zeros(3,1) zeros(3,1) r5G4_til*N'*B5*Gamma_G4 eye(3) ];

omega=OM*[ omega_1 ; sigma ; v_P1 ];
omega_1=omega(1:3);
omega_2=omega(4:6);
omega_3=omega(7:9);
omega_4=omega(10:12);
omega_5=omega(13:15);

omega1_til=skew(omega_1);
omega2_til=skew(omega_2);
omega3_til=skew(omega_3);
omega4_til=skew(omega_4);
omega5_til=skew(omega_5);

% v_P=V*[ omega_1 ; sigma ; v_P1 ];
% v_P1=v_P(1:3);
% v_P2=v_P(4:6);
% v_P3=v_P(7:9);
% v_P4=v_P(10:12);
% v_P5=v_P(13:15);


% F) evaluate the remainder accelerations alpha_R and a_R
alpha_1R=zeros(3,1);
alpha_2R=omega2_til*Gamma_G1*sigma(1);
alpha_3R=omega3_til*Gamma_G2*sigma(2);
alpha_4R=omega4_til*Gamma_G3*sigma(3);
alpha_5R=omega5_til*Gamma_G4*sigma(4);
alpha_R=[alpha_1R;alpha_2R;alpha_3R;alpha_4R;alpha_5R];

a_1R=zeros(3,1);
a_2R=-skew(N'*B2*omega2_til*Gamma_G1*sigma(1))*r_2toG1+...
    skew(N'*B1*omega_1)*skew(N'*B1*omega_1)*r_1toG1-...
    skew(N'*B2*omega_2)*skew(N'*B2*omega_2)*r_2toG1;
a_3R=-skew(N'*B3*omega3_til*Gamma_G2*sigma(2))*r_3toG2+...
    skew(N'*B1*omega_1)*skew(N'*B1*omega_1)*r_1toG2-...
    skew(N'*B3*omega_3)*skew(N'*B3*omega_3)*r_3toG2;
a_4R=-skew(N'*B4*omega4_til*Gamma_G3*sigma(3))*r_4toG3+...
    skew(N'*B1*omega_1)*skew(N'*B1*omega_1)*r_1toG3-...
    skew(N'*B4*omega_4)*skew(N'*B4*omega_4)*r_4toG3;
a_5R=-skew(N'*B5*omega5_til*Gamma_G4*sigma(4))*r_5toG4+...
    skew(N'*B1*omega_1)*skew(N'*B1*omega_1)*r_1toG4-...
    skew(N'*B5*omega_5)*skew(N'*B5*omega_5)*r_5toG4;
a_R=[a_1R;a_2R;a_3R;a_4R;a_5R];


% G) evaluate the inertia forces and torques for each B_i
M=[m_b1*eye(3) zeros(3) zeros(3) zeros(3) zeros(3) ;...
   zeros(3) m_bb*eye(3) zeros(3) zeros(3) zeros(3) ;...
   zeros(3) zeros(3) m_bb*eye(3) zeros(3) zeros(3) ;...
   zeros(3) zeros(3) zeros(3) m_bb*eye(3) zeros(3) ;...
   zeros(3) zeros(3) zeros(3) zeros(3) m_bb*eye(3) ];

J=[J_b1 zeros(3) zeros(3) zeros(3) zeros(3) ;...
   zeros(3) J_bb zeros(3) zeros(3) zeros(3) ;...
   zeros(3) zeros(3) J_bb zeros(3) zeros(3) ;...
   zeros(3) zeros(3) zeros(3) J_bb zeros(3) ;...
   zeros(3) zeros(3) zeros(3) zeros(3) J_bb ];

omJom=[omega1_til*J_b1*omega_1 ;...
       omega2_til*J_bb*omega_2 ;...
       omega3_til*J_bb*omega_3 ;...
       omega4_til*J_bb*omega_4 ;...
       omega5_til*J_bb*omega_5 ];

F_iner=-M*a_R;
T_iner=-J*alpha_R-omJom;


% H) evaluate the active forces and torques for each B_i
F_1=-m_b1*mu*r_P1/norm(r_P1)^3;
F_2=-m_bb*mu*r_P1/norm(r_P1)^3;
F_3=-m_bb*mu*r_P1/norm(r_P1)^3;
F_4=-m_bb*mu*r_P1/norm(r_P1)^3;
F_5=-m_bb*mu*r_P1/norm(r_P1)^3;
F=[ F_1 ; F_2 ; F_3 ; F_4 ; F_5 ];

T_2=B2'*(-k1*theta(1)-k2*sigma(1))*B2(:,1);
T_3=B3'*(-k1*theta(2)-k2*sigma(2))*B3(:,1);
T_4=B4'*(-k1*theta(3)-k2*sigma(3))*B4(:,1);
T_5=B5'*(-k1*theta(4)-k2*sigma(4))*B5(:,1);
T_1=-(T_2+T_3+T_4+T_5);
T=[ T_1 ; T_2 ; T_3 ; T_4 ; T_5 ];


% I) evaluate the terms in Kane's equation
A=OM'*J*OM+V'*M*V;
B=OM'*(T+T_iner)+V'*(F+F_iner);


% J) solve for u_dot and x_dot
dw(1:10)=A\B;
dw(11)=-0.5*w(1:3)'*w(12:14);
dw(12:14)=0.5*(w(11)*eye(3)+skew(w(12:14)))*w(1:3);
dw(15:18)=w(4:7);
dw(19:21)=w(8:10);

end

