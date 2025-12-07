% HOMEWORK 2

%% WIPEOUT
clear
close all
clc


%% DATA
global m_b1 m_bb J_b1 J_bb k1 k2 N l R mu

in=input('Please select the case: (1-2)\n\n');
if isscalar(in)==0 || isreal(in)==0 || in<1 || in>2 %#ok<COMPNOT>
    error('Expected input must be either 1 or 2')
end

% body 1 (root body)
m_b1=2000;                          % in kg
R=1;                                % in m
L=4;                                % in m
I_s=(m_b1*R^2)/2;                   % in kg^2*m
I_t=(m_b1*R^2)/4+(m_b1*L^2)/12;     % in kg^2*m
J_b1=diag([I_s,I_t,I_t]);

% body 2 through 5
m_bb=100;               % in kg
r_a=0.05;               % in m
switch in
    case 1
        l=2;            % in m
    case 2
        l=6;            % in m
end
I_s2=(m_bb*r_a^2)/2;                    % in kg^2*m
I_t2=(m_bb*r_a^2)/4+(m_bb*l^2)/12;      % in kg^2*m
J_bb=diag([I_t2,I_s2,I_t2]);

% constants
mu=398600.4418*10^9;    % in m^3/s^2
k1=100;                 % in N*m
k2=50;                  % in N*m*sec
switch in
    case 1
        t_f=5*3600;      % in sec
    case 2
        t_f=0.5*3600;    % in sec
end

% frames
n1=[1;0;0] ; n2=[0;1;0]; n3=[0;0;1];
N=[n1 n2 n3];

% initial conditions
R_0=7000*10^3;           % in m

v_P1_0=[0 ; sqrt(mu/R_0) ; 0];
r_P1_0=[R_0 ; 0 ; 0];
omega_1_0=[1 ; 0.1 ; 0.1];
theta_0=deg2rad([-30 ; 30 ; -30 ; 10]);
sigma_0=[-0.2 ; 0.2 ; 0.1 ; -0.3];
R_Nto1=eye(3);
B1_0=R_Nto1*N;
[q01_0,q1_0]=quat_finder(R_Nto1);


%% INTEGRATION
tspan=[0 t_f];
options=odeset('RelTol',1.e-10,'AbsTol',1.e-7);
w0=[ omega_1_0 ; sigma_0 ; v_P1_0 ; q01_0 ; q1_0 ; theta_0 ; r_P1_0 ];

tic
[t,W]=ode45(@dyn_and_kin,tspan,w0,options);
toc
omega_table=W(:,1:3)';
sigma_table=rad2deg(W(:,4:7)');
vp_table=(W(:,8:10)/1000)';
theta_table=rad2deg(W(:,15:18)');
r_P1_table=(W(:,19:21)/1000)';
t=t/3600;


%% PLOT
switch in
    case 1
        str='(case 1)';
    case 2
        str='(case 2)';
end

figure(1)
hold on
plot(t,omega_table(1,:),'r-')
plot(t,omega_table(2,:),'g-')
plot(t,omega_table(3,:),'b-')
legend('$$\omega_{1,1}$$','$$\omega_{1,2}$$','$$\omega_{1,3}$$','Interpreter','latex','Fontsize',12)
title(sprintf('$$\\vec{\\omega}_{B1} (t)$$ %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\vec{\omega}_{B1} \ (rad/sec)$$','Interpreter','latex','Fontsize',12)

figure(2)
hold on
plot(t,vp_table(1,:),'r-')
plot(t,vp_table(2,:),'g-')
plot(t,vp_table(3,:),'b-')
legend('$$v_{P1,1}$$','$$v_{P1,2}$$','$$v_{P1,3}$$','Interpreter','latex','Fontsize',12)
title(sprintf('$$\\vec{v}_{P1} (t)$$ %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\vec{v}_{P1} \ (km/sec)$$','Interpreter','latex','Fontsize',12)

figure(3)
subplot(2,2,1)
plot(t,sigma_table(1,:),'b-')
title('$$\sigma_1 (t)$$','Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\sigma_1 \ (deg/sec)$$','Interpreter','latex','Fontsize',12)
subplot(2,2,2)
plot(t,sigma_table(2,:),'b-')
title('$$\sigma_2 (t)$$','Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\sigma_2 \ (deg/sec)$$','Interpreter','latex','Fontsize',12)
subplot(2,2,3)
plot(t,sigma_table(3,:),'b-')
title('$$\sigma_3 (t)$$','Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\sigma_3 \ (deg/sec)$$','Interpreter','latex','Fontsize',12)
subplot(2,2,4)
plot(t,sigma_table(4,:),'b-')
title('$$\sigma_4 (t)$$','Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\sigma_4 \ (deg/sec)$$','Interpreter','latex','Fontsize',12)
sgtitle(sprintf('$$\\vec{\\sigma} (t)$$ %s',str),'Interpreter','latex','Fontsize',12)

figure(4)
subplot(2,2,1)
plot(t,theta_table(1,:),'b-')
title('$$\theta_1 (t)$$','Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\theta_1 \ (deg)$$','Interpreter','latex','Fontsize',12)
subplot(2,2,2)
plot(t,theta_table(2,:),'b-')
title('$$\theta_2 (t)$$','Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\theta_2 \ (deg)$$','Interpreter','latex','Fontsize',12)
subplot(2,2,3)
plot(t,theta_table(3,:),'b-')
title('$$\theta_3 (t)$$','Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\theta_3 \ (deg)$$','Interpreter','latex','Fontsize',12)
subplot(2,2,4)
plot(t,theta_table(4,:),'b-')
title('$$\theta_4 (t)$$','Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\theta_4 \ (deg)$$','Interpreter','latex','Fontsize',12)
sgtitle(sprintf('$$\\vec{\\theta} (t)$$ %s',str),'Interpreter','latex','Fontsize',12)

figure(5)
plot3(r_P1_table(1,:),r_P1_table(2,:),r_P1_table(3,:),'r-','LineWidth',1)
grid on
axis equal
