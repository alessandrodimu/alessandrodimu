% HOMEWORK 3A

%% WIPEOUT
clear
close all
clc


%% DATA
global in r0 v0 q0e_0 gainA_inv gainB Jc A mu

in=input('Please select the case: (1-2)\n\n');
if isscalar(in)==0 || isreal(in)==0 || in<1 || in>2 %#ok<COMPNOT>
    error('Expected input must be either 1 or 2')
end

% spacecraft
Jc=[1200 100 -200;...
    100 2200 300;...
    -200 300 3100];     % in kgm^2
  
% wheels
I_s=80;                 % in kgm^2
I_t=40;                 % in kgm^2

a1=[1;0;0];
a2=[0;1;0];
a3=[0;0;1];
a4=[1;1;1]/norm([1;1;1]);
A=I_s*[a1,a2,a3,a4];

% orbit information
a=26554;                % in km
e=0.73;
i=deg2rad(63.4);        % in rad
GOM=deg2rad(0);         % in rad
pom=deg2rad(-90);       % in rad

% constants
mu=398600.4;            % in km^3/s^2
omega_n=0.04;           % in rad/s
zita=1;
t_f=600;                % in s

% frames
N=[1 0 0;0 1 0;0 0 1];

% estimation of gains A^(-1) and B
c1=2*omega_n^2;
c2=zita/omega_n;
gainA_inv=c1*eye(3,3);
gainB=c2*eye(3,3);


% initial conditions
omega_0=[-0.05 ; 0.03 ; 0.05];       % in rad/s
omega_s_0=[0.5 ; 0.5 ; -0.5 ; -0.5]; % in rad/s
q_0=[0.1 ; -0.2 ; 0.5];
q0_0=sqrt(1-norm(q_0)^2);
anu_0=deg2rad(90);                   % in rad

coe_0=[a;e;i;GOM;pom;anu_0];
[r0,v0]=coe2rvECI(coe_0,mu);         % r0 and v0 in inertial frame
h=cross(r0,v0);

% estimation of commanded attitute at t0
b_c3_0=-r0/norm(r0);
b_c2_0=-h/norm(h);
b_c1_0=cross(b_c2_0,b_c3_0);

C_0=[b_c1_0 ,b_c2_0,b_c3_0];
R_NtoC=C_0'*N;
[q0c_0,qc_0]=quat_finder(R_NtoC);

% estimation of error quaternions at t0 (from C to B)
qc_0_til=skew(qc_0);

q0e_0=q0c_0*q0_0+qc_0'*q_0;
qe_0=-qc_0*q0_0+q0c_0*q_0-qc_0_til*q_0;


%% INTEGRATION
tspan=[0 t_f];
options=odeset('RelTol',1.e-12,'AbsTol',1.e-9);
w0=[ q0e_0 ; qe_0 ; omega_0 ; omega_s_0 ];

tic
[t,W]=ode45(@control,tspan,w0,options);
toc

q0e_table=W(:,1)';
qe_table=W(:,2:4)';
omega_table=W(:,5:7)';
omega_s_table=W(:,8:11)';


%% ESTIMATION OF TABLES FOR SUBSEQUENT PLOTS
Tc_table=zeros(3,length(t));
Ta_table=zeros(3,length(t));
omega_d_table=zeros(3,length(t));


% orbit propagation
for i=1:length(t)
    [r,v]=CorePropagator(r0,v0,t(i),mu);
   
    v_r=dot(v,r/norm(r));
    h=cross(r,v);

    R_CtoB=(q0e_table(:,i)^2-qe_table(:,i)'*qe_table(:,i))*eye(3,3)...
           +2*(qe_table(:,i)*qe_table(:,i)')-2*q0e_table(:,i)*skew(qe_table(:,i));
    
    anu_dot=norm(h)/norm(r)^2;
    omega_c=[0 ; -anu_dot ; 0];
    anu_ddot=-2*norm(h)*v_r/norm(r)^3;
    omega_c_dot=[0 ; -anu_ddot ; 0];

    % error angular velocity
    omega_d_table(:,i)=omega_table(:,i)-R_CtoB*omega_c;
    
    % commanded torque
    omega_til=skew(omega_table(:,i));
    Mc=[0;0;0];
    Tc_table(:,i)=omega_til*Jc*omega_table(:,i)-Mc+Jc*R_CtoB*omega_c_dot...
                  -Jc*gainA_inv*gainB*omega_d_table(:,i)...
                  -sign(q0e_0)*Jc*gainA_inv*qe_table(:,i);
    
    % actual torque
    switch in
        case 1
            Ta_table(:,i)=Tc_table(:,i);
        case 2
            Ta_table(:,i)=-omega_til*A*omega_s_table(:,i)+Tc_table(:,i);
    end
end


%% PLOT
switch in
    case 1
        str='(case 1)';
    case 2
        str='(case 2)';
end

figure(1)
hold on
plot(t,q0e_table,'k');
plot(t,qe_table(1,:),'r');
plot(t,qe_table(2,:),'g');
plot(t,qe_table(3,:),'b');
legend('$$q_{e0}$$','$$q_{e1}$$','$$q_{e2}$$','$$q_{e3}$$','Interpreter','latex','Fontsize',12);
title(sprintf('Error quaternion %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (sec)$$','Interpreter','latex','Fontsize',12)

figure(2)
hold on
plot(t,omega_d_table(1,:),'r');
plot(t,omega_d_table(2,:),'g');
plot(t,omega_d_table(3,:),'b');
legend('$$\omega_{D1}$$','$$\omega_{D2}$$','$$\omega_{D3}$$','Interpreter','latex','Fontsize',12);
title(sprintf('Error angular velocity %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (sec)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\vec{\omega}_D \ (rad/sec)$$','Interpreter','latex','Fontsize',12);

figure(3)
hold on
plot(t,Tc_table(1,:),'r');
plot(t,Tc_table(2,:),'g');
plot(t,Tc_table(3,:),'b');
legend('$$T_{c1}$$','$$T_{c2}$$','$$T_{c3}$$','Interpreter','latex','Fontsize',12);
title(sprintf('Commanded torque %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (sec)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\vec{T}_c \ (Nm) $$','Interpreter','latex','Fontsize',12);
ylim([-20 10])

switch in
    case 2
        figure(4)
        hold on
        plot(t,omega_s_table(1,:),'r');
        plot(t,omega_s_table(2,:),'g');
        plot(t,omega_s_table(3,:),'b');
        plot(t,omega_s_table(4,:),'k');
        legend('$$\omega_{S1}$$','$$\omega_{S2}$$','$$\omega_{S3}$$','$$\omega_{S4}$$','Interpreter','latex','Fontsize',12);
        title(sprintf('Wheels angular velocity %s',str),'Interpreter','latex','Fontsize',12)
        xlabel('$$t \ (sec)$$','Interpreter','latex','Fontsize',12)
        ylabel('$$\vec{\omega}_S \ (rad/sec) $$','Interpreter','latex','Fontsize',12);

        figure(5)
        hold on
        plot(t,Ta_table(1,:),'r');
        plot(t,Ta_table(2,:),'g');
        plot(t,Ta_table(3,:),'b');
        legend('$$T_{a1}$$','$$T_{a2}$$','$$T_{a3}$$','Interpreter','latex','Fontsize',12);
        title(sprintf('Actual torque %s',str),'Interpreter','latex','Fontsize',12)
        xlabel('$$t \ (sec)$$','Interpreter','latex','Fontsize',12)
        ylabel('$$\vec{T}_a \ (Nm) $$','Interpreter','latex','Fontsize',12);
        ylim([-20 10])
end

