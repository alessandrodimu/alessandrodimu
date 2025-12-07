% HOMEWORK 1

%% WIPEOUT
clear
close all
clc


%% DATA
global m_b m_w m_d cd kd kw b_w I_s I_t in

in=input('Please select the case: (1-4)\n\n');
if isscalar(in)==0 || isreal(in)==0 || in<1 || in>4 %#ok<COMPNOT>
    error('Expected input must be between 1 and 4')
end

% body
m_b=100;                        % in kg

% wheel (x3)
m_w=5;                          % in kg
R_w=1;                          % in m
b_w=2;                          % in m
I_s=(m_w*R_w^2)/2;              % in kg^2*m
I_t=(m_w*R_w^2)/4;              % in kg^2*m

% damper
m_d=10;                         % in kg
kd=3.5;
cd=30;

% constants
kw=0.1;                         % in kg^2*m/sec
switch in
    case {1,2}
        t_f=5*3600;             % in s
    case {3,4}
        t_f=1*3600;             % in s
end

% frames
b1=[1;0;0];
b2=[0;1;0];
b3=[0;0;1];
B=[b1 b2 b3];

% initial conditions
v_P_0=[0 ; 0 ; 0];              % in m/s
w_0=deg2rad([36 ; 3 ; 3]);      % in rad/s
xi_0=0;                         % in m
xi_dot_0=0;                     % in m/s
switch in
    case 1
        w_S_0=deg2rad([0 ; 0 ; 0]);        % in rad/s 
    case 2
        w_S_0=deg2rad([1080 ; 0 ; 0]);     % in rad/s
    case 3
        w_S_0=deg2rad([1080 ; 60 ; -30]);  % in rad/s
    case 4
        w_S_0=deg2rad([360 ; 60 ; -30]);   % in rad/s
end
eta_0=[v_P_0 ; w_0 ; xi_dot_0 ; w_S_0 ];



% mechanical energy estimation at t0
M=m_b+3*m_w+m_d;

S_BP=0;
S_WP=m_w*(b_w*b1)+m_w*(b_w*b2)+m_w*(b_w*b3);
S_DP=m_d*xi_0*b1;
S_P=S_BP+S_WP+S_DP;

J_BP=diag([350,300,400]);
J_WP=diag([I_s,I_t,I_t])+m_w*(norm(b_w*b1)^2*eye(3,3)-(b_w*b1)*(b_w*b1)')+...
     diag([I_t,I_s,I_t])+m_w*(norm(b_w*b2)^2*eye(3,3)-(b_w*b2)*(b_w*b2)')+...
     diag([I_t,I_t,I_s])+m_w*(norm(b_w*b3)^2*eye(3,3)-(b_w*b3)*(b_w*b3)');
J_DP=m_d*(norm(xi_0*b1)^2*eye(3,3)-(xi_0*b1)*(xi_0*b1'));
J_P=J_BP+J_WP+J_DP;

S_Ptil=skew(S_P);

MT_0=[M*eye(3,3) -S_Ptil m_d*b1 zeros(3,1) zeros(3,1) zeros(3,1);...
      S_Ptil J_P zeros(3,1) I_s*b1 I_s*b2 I_s*b3;...
      m_d*b1' zeros(1,3) m_d 0 0 0;...
      zeros(1,3) I_s*b1' 0 I_s 0 0;...
      zeros(1,3) I_s*b2' 0 0 I_s 0;...
      zeros(1,3) I_s*b3' 0 0 0 I_s];

T_0=0.5*eta_0'*MT_0*eta_0;
V_0=0.5*kd*xi_0^2;
E_0=T_0+V_0;


%% INTEGRATION
tspan=[0 t_f];
options=odeset('RelTol',1.e-11,'AbsTol',1.e-8);
X_0=[eta_0; xi_0 ; E_0];

tic
[t,W]=ode45(@dyn_eq,tspan,X_0,options);
toc
vp_table=W(:,1:3)';
omega_table=W(:,4:6)';
xi_dot_table=W(:,7)';
omega_S_table=W(:,8:10)';
xi_table=W(:,11)';
E_table=W(:,12)';
t=t/3600;


%% PLOT
switch in
    case 1
        str='(case 1)';
    case 2
        str='(case 2)';
    case 3
        str='(case 3)';
    case 4
        str='(case 4)';
end
figure(1)
hold on
plot(t,vp_table(1,:),'r-')
plot(t,vp_table(2,:),'g-')
plot(t,vp_table(3,:),'b-')
legend('$$v_{P1}$$','$$v_{P2}$$','$$v_{P3}$$','Interpreter','latex','Fontsize',12)
title(sprintf('$$\\vec{v}_P (t)$$ %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\vec{v}_P \ (m/sec)$$','Interpreter','latex','Fontsize',12)

figure(2)
hold on
plot(t,omega_table(1,:),'r-')
plot(t,omega_table(2,:),'g-')
plot(t,omega_table(3,:),'b-')
legend('$$\omega_1$$','$$\omega_2$$','$$\omega_3$$','Interpreter','latex','Fontsize',12)
title(sprintf('$$\\vec{\\omega} (t)$$ %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\vec{\omega} \ (rad/sec)$$','Interpreter','latex','Fontsize',12)

figure(3)
hold on
plot(t,omega_S_table(1,:),'r-')
plot(t,omega_S_table(2,:),'g-')
plot(t,omega_S_table(3,:),'b-')
legend('$$\omega_{S,1}$$','$$\omega_{S,2}$$','$$\omega_{S,3}$$','Interpreter','latex','Fontsize',12)
title(sprintf('$$\\vec{\\omega}_S (t)$$ %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\vec{\omega}_{S} \ (rad/sec)$$','Interpreter','latex','Fontsize',12)

figure(4)
plot(t,xi_dot_table,'k-')
title(sprintf('$$\\dot{\\xi} (t)$$ %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\dot{\xi} \ (m/sec)$$','Interpreter','latex','Fontsize',12)

figure(5)
plot(t,xi_table,'k-')
title(sprintf('$$\\xi (t)$$ %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\xi \ (m)$$','Interpreter','latex','Fontsize',12)

figure(6)
plot(t,E_table,'k-')
title(sprintf('$$E (t)$$ %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$E$$','Interpreter','latex','Fontsize',12)