% HOMEWORK 5

%% WIPEOUT

clear
close all
clc



%% DATA

global mu c u_T_max r_ini v_ini r_fin v_fin

% constants
mu=398600.4;                        % in km^3/s^2
R_e=6378.136;                       % in km
DU=R_e;                             % in km
TU=sqrt(DU^3/mu);                   % in s
mu=mu/(DU^3/TU^2);                  % adimensional


% orbit values
r_ini=10000/DU;                     % adimensional
v_ini=sqrt(mu/r_ini);               % adimensional
r_fin=42164/DU;                     % adimensional
v_fin=sqrt(mu/r_fin);               % adimensional


% propulsion parameters
g0=9.8065*10^(-3)/(DU/TU^2);        % adimensional
c=30/(DU/TU);                       % adimensional
u_T_max=0.01*g0;                    % adimensional



%% INITIALIZATION

% initial and final orbit
r_0=r_ini;                  % adimensional
csi_0=0;                    % adimensional
vr_0=0;                     % adimensional
vt_0=v_ini;                 % adimensional

r_f=r_fin;                  % adimensional
csi_f=[];                   % we don't care about it
vr_f=0;                     % adimensional
vt_f=v_fin;                 % adimensional


% first guess values
lambda_10_g=-1;
lambda_30_g=0.2;
lambda_40_g=-1;
t_f_g=40;                   % adimensional



%% OPTIMIZATION

options_opt=optimset('Display','off','PlotFcns','optimplotfval','TolFun',...
            1e-10,'TolX',1e-10,'MaxFunEvals',800);
x0_g=[lambda_10_g;lambda_30_g;lambda_40_g;t_f_g];

x0_opt=fminsearch('J_function',x0_g,options_opt);

lambda_10=x0_opt(1);
lambda_30=x0_opt(2);
lambda_40=x0_opt(3);
t_f=x0_opt(4);



%% INTEGRATION

tspan=[0 t_f];
options_int=odeset('RelTol',1e-10,'AbsTol',1e-8);
x0=[r_0; csi_0; vr_0 ; vt_0 ; lambda_10 ; lambda_30 ; lambda_40];

[t,W]=ode45(@dyn_kin_and_adjoints,tspan,x0,options_int);

r_table=W(:,1)';
csi_table=W(:,2)';
vr_table=W(:,3)';
vt_table=W(:,4)';
lambda_1_table=W(:,5)';
lambda_2_table=zeros(1,length(t));
lambda_3_table=W(:,6)';
lambda_4_table=W(:,7)';

sinalpha=-lambda_3_table./sqrt(lambda_3_table.^2+lambda_4_table.^2);
cosalpha=-lambda_4_table./sqrt(lambda_3_table.^2+lambda_4_table.^2);
alpha_table=atan2(sinalpha,cosalpha);

t=t*TU/3600;               % in hours

H=lambda_1_table(end)*vr_table(end)+lambda_2_table(end)*vt_table(end)/r_table(end)+...
  lambda_3_table(end)*(-mu/r_table(end)^2+vt_table(end)^2/r_table(end)+u_T_max*c*sin(alpha_table(end))/(c-u_T_max*t(end)))+...
  lambda_4_table(end)*(-vr_table(end)*vt_table(end)/r_table(end)+u_T_max*c*cos(alpha_table(end))/(c-u_T_max*t(end)));


% conversion for subsequent plots
r_table=r_table*DU;                     % in km
vr_table=vr_table*(DU/TU);              % in km/s
vt_table=vt_table*(DU/TU);              % in km/s
alpha_table=rad2deg(alpha_table);       % in deg

r_ini=r_ini*DU;                         % in km
r_fin=r_fin*DU;                         % in km
v_ini=v_ini*DU/TU;                      % in km/s
v_fin=v_fin*DU/TU;                      % in km/s


% print values at the orbit injection
delta_r=abs(r_table(end)-r_fin);
delta_vr=abs(vr_table(end));
delta_vt=abs(vt_table(end)-v_fin);

fprintf('\nErrors at the orbit injection\n\n');
fprintf('delta_r = %e km \n',delta_r);
fprintf('delta_vr = %e km/s \n',delta_vr);
fprintf('delta_vt = %e km/s \n\n',delta_vt);
fprintf('t_f = %3.1f hours\n',t(end));
fprintf('H = %5.4f\n',H);



%% PLOT

figure(2)
plot(t,r_table,'b')
title('Radius','Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$r \ (km)$$','Interpreter','latex','Fontsize',12)
xlim([0 t(end)])

figure(3)
hold on
plot(t,vr_table,'b')
plot(t,vt_table,'r')
title('Velocity components','Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$v \ (km/s)$$','Interpreter','latex','Fontsize',12)
legend('$$v_r$$','$$v_t$$','Interpreter','latex','Fontsize',12)
xlim([0 t(end)])

figure(4)
plot(t,alpha_table,'b')
title('Thrust pointing angle','Interpreter','latex','Fontsize',12)
xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\alpha \ (deg)$$','Interpreter','latex','Fontsize',12)
xlim([0 t(end)])

% figure(5)
% hold on
% plot(t,lambda_1_table,'r')
% plot(t,lambda_2_table,'g')
% plot(t,lambda_3_table,'b')
% plot(t,lambda_4_table,'k')
% title('Lagrangian multipliers','Interpreter','latex','Fontsize',12)
% xlabel('$$t \ (hrs)$$','Interpreter','latex','Fontsize',12)
% ylabel('$$\lambda$$','Interpreter','latex','Fontsize',12)
% legend('$$\lambda_1$$','$$\lambda_2$$','$$\lambda_3$$','$$\lambda_4$$','Interpreter','latex','Fontsize',12)
% xlim([0 t(end)])

figure(6)
[x_ini_table,y_ini_table]=pol2cart(linspace(0,2*pi,1000),r_ini*ones(1,1000));
[x_fin_table,y_fin_table]=pol2cart(linspace(0,2*pi,1000),r_fin*ones(1,1000));
[x_table,y_table]=pol2cart(csi_table,r_table);
hold on
plot(x_ini_table,y_ini_table,'k')
plot(x_fin_table,y_fin_table,'k')
plot(x_table,y_table,'b')
title('Transfer trajectory','Interpreter','latex','Fontsize',12)
xlabel('$$x \ (km)$$','Interpreter','latex','Fontsize',12)
ylabel('$$y \ (km)$$','Interpreter','latex','Fontsize',12)
legend('Initial orbit','Final orbit','Transfer trajectory','Interpreter','latex')
axis('equal')

