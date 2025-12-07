% HOMEWORK 4

%% WIPEOUT
clear
close all
clc



%% DATA
global mu R_e J2 DU K p_d e_d i_d u_T_max c

in=input('Please select the case: (1-2)\n\n');
if isscalar(in)==0 || isreal(in)==0 || in<1 || in>2 %#ok<COMPNOT>
    error('Expected input must be either 1 or 2')
end

% constants
mu=398600.4;                    % in km^3/s^2
R_e=6378.136;                   % in km
J2=1.082*10^(-3);
DU=R_e;                         % in km
TU=sqrt(DU^3/mu);               % in s
mu=mu/(DU^3/TU^2);              % dimensionless (before in km^3/s^2)
t_f=30*86400/TU;                % dimensionless (before in s)

% initial orbit
switch in
    case 1
        p_0=18000/DU;           % dimensionless (before in km)
        e_0=0.12;
        i_0=deg2rad(45);        % in rad
        GOM_0=deg2rad(30);      % in rad
        pom_0=deg2rad(0);       % in rad
        anu_0=deg2rad(0);       % in rad
    case 2
        p_0=8000/DU;            % dimensionless (before in km)
        e_0=0.12;
        i_0=deg2rad(45);        % in rad
        GOM_0=deg2rad(30);      % in rad
        pom_0=deg2rad(0);       % in rad
        anu_0=deg2rad(0);       % in rad
end

% final orbit
switch in
    case 1
        p_d=20000/DU;           % dimensionless (before in km)
        e_d=0;                  
        i_d=deg2rad(50);        % in rad
    case 2
        p_d=10000/DU;           % dimensionless (before in km)
        e_d=0;                  
        i_d=deg2rad(50);        % in rad
end


% thruster specs
g0=(9.8065*10^(-3))/(DU/TU^2);  % dimensionless (before in m/s^2)
u_T_max=10^(-4)*g0;             % dimensionless (before in m/s^2)
c=30/(DU/TU);                   % dimensionless (before in km/s)

% gains
k1=1;
k2=10^4;
k3=1;
K=diag([k1 k2 k3]);



%% INITIALIZATION

% from orbital parameters to equinoctial elements
l_0=e_0*cos(GOM_0+pom_0);
m_0=e_0*sin(GOM_0+pom_0);
n_0=tan(i_0/2)*cos(GOM_0);
s_0=tan(i_0/2)*sin(GOM_0);
q_0=GOM_0+pom_0+anu_0;



%% INTEGRATION
tspan=[0 t_f];
options=odeset('RelTol',1.e-10,'AbsTol',1.e-8);
x_0=[p_0;l_0;m_0;n_0;s_0;q_0;1];

tic
[t,W]=ode45(@dyn_eq,tspan,x_0,options);
toc


p_table=W(:,1)';

x2_table=W(:,2)';
x3_table=W(:,3)';
e_table=sqrt(x2_table.^2+x3_table.^2);

x4_table=W(:,4)';
x5_table=W(:,5)';
i_table=2*atan(sqrt(x4_table.^2+x5_table.^2));

x7_table=W(:,7)';

V_table=zeros(1,length(t));
for i=1:length(t)
    psi=[p_table(i)-p_d;...
         e_table(i)^2-e_d^2;...
         tan(i_table(i)/2)^2-tan(i_d/2)^2];
    V_table(i)=0.5*psi'*K*psi; 
end



%% PLOT
switch in
    case 1
        str='(case 1)';
    case 2
        str='(case 2)';
end

% conversion for plots
t=t*TU/86400;
p_table=p_table*DU;
i_table=rad2deg(i_table);


figure(1)
plot(t,p_table,'b');
title(sprintf('Semilatus rectum %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (days)$$','Interpreter','latex','Fontsize',12)
ylabel('$$p \ (km)$$','Interpreter','latex','Fontsize',12);

figure(2)
plot(t,e_table,'b');
title(sprintf('Eccentricity %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (days)$$','Interpreter','latex','Fontsize',12)
ylabel('$$e$$','Interpreter','latex','Fontsize',12);

figure(3)
plot(t,i_table,'b');
title(sprintf('Inclination %s',str),'Interpreter','latex','Fontsize',12)
ylim([44 51])
xlabel('$$t \ (days)$$','Interpreter','latex','Fontsize',12)
ylabel('$$i \ (deg)$$','Interpreter','latex','Fontsize',12);

figure(4)
plot(t,x7_table,'b');
title(sprintf('Mass ratio %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (days)$$','Interpreter','latex','Fontsize',12)
ylabel('$$x_{7}$$','Interpreter','latex','Fontsize',12);

figure(5)
plot(t,V_table,'b');
title(sprintf('Lyapunov function %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (days)$$','Interpreter','latex','Fontsize',12)
ylabel('$$V$$','Interpreter','latex','Fontsize',12);

