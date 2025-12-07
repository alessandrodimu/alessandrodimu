% HOMEWORK 6

%% WIPEOUT

clear
close all
clc



%% DATA

in=input('Please select the case: (1-2)\n\n');
if isscalar(in)==0 || isreal(in)==0 || in<1 || in>2 %#ok<COMPNOT>
    error('Expected input must be either 1 or 2')
end


global m w_e mu Cd S R_e ro_0 beta C v_EI

% Soyuz capsule
m=2400;                     % in kg
S=3.8;                      % in m^2
Cd=1.341;                   % adimensional
Rc=2.235;                   % in m

% constants
mu=398600.4*10^9;           % in m^3/s^2
R_e=6378.136*10^3;          % in m
w_e=2*pi/86164;             % in rad/s
ro_0=1.225;                 % in kg/m^3
beta=0.1378*10^(-3);        % in m^(-1)
g0=9.8065;                  % in m/s^2
kp=1.762*10^(-4);           % in kg^(1/2)*m^(-1)

% initial conditions
h_0=400*10^3;               % in m
e_0=0;                      % adimensional
i_0=deg2rad(0);             % in rad

% conditions at EI
h_EI=120*10^3;              % in m

% tangential braking (considered impulsive)
switch in
    case 1
        dv=0.15*10^3;       % in m/s
    case 2
        dv=0.3*10^3;        % in m/s
end



%% ESTIMATION
% step A
r_0=R_e+h_0;                                            % in m
v_0=sqrt(mu/r_0)-dv;                                    % in m/s
r_EI=R_e+h_EI;                                          % in m

% conservation of specific mechanical energy
sme_0=v_0^2/2-mu/r_0;                                   % in m^2/s^2
v_iini=sqrt(2*(sme_0+mu/r_EI));                         % in m/s     

% conservation of specific angular momentum
sam_0=r_0*v_0;                                          % in m^2/s
gamma_iini=-acos(sam_0/(r_EI*v_iini));                  % in rad


% step B
a_t=-mu/(2*sme_0);                                      % in m
r_per=2*a_t-r_0;                                        % in m
e_t=(r_0-r_per)/(r_0+r_per);                            % adimensional
p_t=sam_0^2/mu;                                         % in m

anu_iini=-acos((p_t-r_EI)/(r_EI*e_t));                  % in rad

% position and velocity vectors in the orbital rotating frame
v_iini=[v_iini*sin(gamma_iini) ; v_iini*cos(gamma_iini) ; 0];   % in m/s
r_iini=r_EI*[1;0;0];                                            % in m

v_rini=v_iini-cross([0;0;w_e],r_iini);                          % in m/s
gamma_rini=atan(v_rini(1)/v_rini(2));                           % in rad


% step C
v_EI=norm(v_rini);                      % in m/s
gamma_EI=gamma_rini;                    % in rad
t_f=10*60;                              % in s

tspan=[0 t_f];
options=odeset('RelTol',1.e-12,'AbsTol',1.e-9);
x_0=[r_EI;pi/2;v_EI;gamma_EI];

[t,W]=ode45(@dyn_and_kin,tspan,x_0,options);

ind_impact=(W(:,1)>R_e);
W=W(ind_impact,:);
t=t(ind_impact);
t_f=t(end);
fprintf('\nThe capsule impacted the ground at %g sec\n', t_f);

r_table=W(:,1)';                    % in m
lambdag_table=W(:,2)';              % in rad/s
vr_table=W(:,3)';                   % in m/s
gammar_table=W(:,4)';               % in rad/s
h_table=r_table-R_e;                % in m


% step D
aD_table=zeros(1,length(t));
qA_table=zeros(1,length(t));
qS_table=zeros(1,length(t));
pD_table=zeros(1,length(t));

for i=1:length(t)
    ro=ro_0*exp(-beta*h_table(i));                  % in kg/m^3
    D=0.5*Cd*ro*S*vr_table(i)^2;                    % in N
    a_D=D/m;                                        % in m/s^2
    
    r_dot=vr_table(i)*sin(gammar_table(i));                 % in m/s
    vr_dot=(-mu/r_table(i)^2)*sin(gammar_table(i))-...      % in m/s^2
        a_D+w_e^2*r_table(i)*sin(gammar_table(i));
    gammar_dot=(vr_table(i)/r_table(i)...                   % in rad/s
        -mu/(r_table(i)^2*vr_table(i)))*cos(gammar_table(i))+...
        (w_e^2*r_table(i)/vr_table(i))*cos(gammar_table(i))+2*w_e;
    
    aD_table(i)=-0.5/sqrt(vr_table(i)^2+w_e^2*r_table(i)^2+2*w_e*r_table(i)*vr_table(i)*cos(gammar_table(i)))*...
        (2*vr_table(i)*vr_dot+2*w_e^2*r_table(i)*r_dot+...
        2*w_e*r_dot*vr_table(i)*cos(gammar_table(i))+2*w_e*r_table(i)*vr_dot*cos(gammar_table(i))-...
        2*w_e*r_table(i)*vr_table(i)*sin(gammar_table(i))*gammar_dot);      % in m/s^2
    qA_table(i)=0.25*ro*vr_table(i)^3;                                      % in W/m^2
    qS_table(i)=kp*(sqrt(ro)*vr_table(i)^3/sqrt(Rc));                       % in W/m^2
    pD_table(i)=0.5*ro*vr_table(i)^2;                                       % in Pa
end    

% conversion for subsequent plots and max values (and times and heights)
aD_max=max(aD_table);                       % in m/s^2
t_aD_max=t(aD_table==aD_max);               % in s
h_aD_max=h_table(aD_table==aD_max)/10^3;    % in km

qA_table=qA_table/10^3;                     % in kW/m^2
qA_max=max(qA_table);                       % in kW/m^2
t_qA_max=t(qA_table==qA_max);               % in s
h_qA_max=h_table(qA_table==qA_max)/10^3;    % in km

qS_table=qS_table/10^3;                     % in kW/m^2
qS_max=max(qS_table);                       % in kW/m^2
t_qS_max=t(qS_table==qS_max);               % in s
h_qS_max=h_table(qS_table==qS_max)/10^3;    % in km

pD_max=max(pD_table);                       % in Pa
t_pD_max=t(pD_table==pD_max);               % in s
h_pD_max=h_table(pD_table==pD_max)/10^3;    % in km


% step e
B=Cd*S/m;                                   % in m^2/kg
C=B*ro_0/(2*beta*sin(gamma_EI));            % adimensional

h_AE=h_table;                               % in m 
vr_AE=v_EI*exp(C*exp(-beta*h_AE));          % in m/s
gamma_AE=gamma_EI*ones(1,length(h_AE));     % in rad

aD_AE=0.5*B*ro_0*v_EI^2*exp(-beta*h_AE+2*C*exp(-beta*h_AE));                % in m/s^2
qA_AE=0.25*ro_0*v_EI^3*exp(-beta*h_AE+3*C*exp(-beta*h_AE));                 % in W/m^2
qS_AE=kp*sqrt(ro_0/Rc)*v_EI^3*exp(-0.5*beta*h_AE+3*C*exp(-beta*h_AE));      % in W/m^2
pD_AE=0.5*ro_0*exp(-beta*h_AE).*vr_AE.^2;                                   % in Pa
 
% conversion for subsequent plots
h_AE=h_AE/10^3;                     % in km
vr_AE=vr_AE/10^3;                   % in km/s
gamma_AE=rad2deg(gamma_AE);         % in rad

aD_AE_max=max(aD_AE);                       % in m/s^2
t_aD_AE_max=t(aD_AE==aD_AE_max);            % in s
h_aD_AE_max=h_AE(aD_AE==aD_AE_max);         % in km

qA_AE=qA_AE/10^3;                           % in kW/m^2
qA_AE_max=max(qA_AE);                       % in kW/m^2
t_qA_AE_max=t(qA_AE==qA_AE_max);            % in s
h_qA_AE_max=h_AE(qA_AE==qA_AE_max);         % in km

qS_AE=qS_AE/10^3;                           % in kW/m^2
qS_AE_max=max(qS_AE);                       % in kW/m^2
t_qS_AE_max=t(qS_AE==qS_AE_max);            % in s
h_qS_AE_max=h_AE(qS_AE==qS_AE_max);         % in km

pD_AE_max=max(pD_AE);                       % in Pa
t_pD_AE_max=t(pD_AE==pD_AE_max);            % in s
h_pD_AE_max=h_AE(pD_AE==pD_AE_max);         % in km



%% PRINTING VALUES

% conversion for printing and plots
r_table=r_table/10^3;                   % in km
h_table=h_table/10^3;                   % in km
vr_table=vr_table/10^3;                 % in km/s
gammar_table=rad2deg(gammar_table);     % in deg

% printing max values
fprintf('\nNumerical integration\n\n');
fprintf('\t\t\t\t |\t max value\t |\t h (km)\t |\t t (sec)\n');
fprintf('----------------------------------------------------------\n');
fprintf('aD_max (m/s^2)\t |\t %7.2f\t |\t %4.2f\t |\t %5.2f\n',[aD_max h_aD_max t_aD_max]);
fprintf('qA_max (MW/m^2)\t |\t %7.2f\t |\t %4.2f\t |\t %5.2f\n',[qA_max/10^3 h_qA_max t_qA_max]);
fprintf('qS_max (kW/m^2)\t |\t %7.2f\t |\t %4.2f\t |\t %5.2f\n',[qS_max h_qS_max t_qS_max]);
fprintf('pD_max (kPa)\t |\t %7.2f\t |\t %4.2f\t |\t %5.2f\n\n',[pD_max/10^3 h_pD_max t_pD_max]);

fprintf('\nAllen-Eggers solution\n\n');
fprintf('\t\t\t\t |\t max value\t |\t h (km)\t |\t t (sec)\n');
fprintf('----------------------------------------------------------\n');
fprintf('aD_max (m/s^2)\t |\t %7.2f\t |\t %4.2f\t |\t %5.2f\n',[aD_AE_max h_aD_AE_max t_aD_AE_max]);
fprintf('qA_max (MW/m^2)\t |\t %7.2f\t |\t %4.2f\t |\t %5.2f\n',[qA_AE_max/10^3 h_qA_AE_max t_qA_AE_max]);
fprintf('qS_max (kW/m^2)\t |\t %7.2f\t |\t %4.2f\t |\t %5.2f\n',[qS_AE_max h_qS_AE_max t_qS_AE_max]);
fprintf('pD_max (kPa)\t |\t %7.2f\t |\t %4.2f\t |\t %5.2f\n\n',[pD_AE_max/10^3 h_pD_AE_max t_pD_AE_max]);



%% PLOT
switch in
    case 1
        str='(case 1)';
    case 2
        str='(case 2)';
end

figure(1)
plot(t/60,h_table,'b');
title(sprintf('Altitude %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (min)$$','Interpreter','latex','Fontsize',12)
ylabel('$$h \ (km)$$','Interpreter','latex','Fontsize',12)
xlim([0 t_f/60])

figure(2)
plot(t/60,vr_table,'b');
title(sprintf('Relative velocity %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (min)$$','Interpreter','latex','Fontsize',12)
ylabel('$$v_{R} \ (km/s)$$','Interpreter','latex','Fontsize',12)
xlim([0 t_f/60])

figure(3)
plot(t/60,gammar_table,'b');
title(sprintf('Flight path angle %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$t \ (min)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\gamma_{R} \ (deg)$$','Interpreter','latex','Fontsize',12)
xlim([0 t_f/60])

figure(4)
[xe_table,ye_table]=pol2cart(linspace(0,2*pi,1000),R_e/10^3*ones(1,1000));
[x_table,y_table]=pol2cart(lambdag_table,r_table);
plot(xe_table,ye_table,'k')
hold on
plot(x_table,y_table,'b')
title(sprintf('Reentry trajectory %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$x \ (km)$$','Interpreter','latex','Fontsize',12)
ylabel('$$y \ (km)$$','Interpreter','latex','Fontsize',12)
axis equal
xlim([1.1*min(x_table) 100])
switch in
    case 1
        ylim([5500 7000])   
    case 2
        ylim([5800 6700])
end


figure(5)
plot(h_table,gammar_table,'b');
hold on
plot(h_AE,gamma_AE,'r');
title(sprintf('Flight path angle %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$h \ (km)$$','Interpreter','latex','Fontsize',12)
ylabel('$$\gamma_R \ (g_0)$$','Interpreter','latex','Fontsize',12)
legend('Numerical integration','Allen-Eggers solution')

figure(6)
plot(h_table,vr_table,'b');
hold on
plot(h_AE,vr_AE,'r');
title(sprintf('Relative velocity %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$h \ (km)$$','Interpreter','latex','Fontsize',12)
ylabel('$$v_R \ (km/s)$$','Interpreter','latex','Fontsize',12)
legend('Numerical integration','Allen-Eggers solution')

figure(7)
plot(h_table,aD_table/g0,'b');
hold on
plot(h_AE,aD_AE/g0,'r');
plot(h_aD_max,aD_max/g0,'bo','MarkerFaceColor','b');
plot(h_aD_AE_max,aD_AE_max/g0,'ro','MarkerFaceColor','r');
title(sprintf('Deceleration %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$h \ (km)$$','Interpreter','latex','Fontsize',12)
ylabel('$$a_D \ (g0)$$','Interpreter','latex','Fontsize',12)
legend('Numerical integration','Allen-Eggers solution')

figure(8)
plot(h_table,qA_table,'b');
hold on
plot(h_AE,qA_AE,'r');
plot(h_qA_max,qA_max,'bo','MarkerFaceColor','b');
plot(h_qA_AE_max,qA_AE_max,'ro','MarkerFaceColor','r');
title(sprintf('Average heat rate per surface unit %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$h \ (km)$$','Interpreter','latex','Fontsize',12)
ylabel('$$q_A \ (kW/m^2)$$','Interpreter','latex','Fontsize',12)
legend('Numerical integration','Allen-Eggers solution')

figure(9)
plot(h_table,qS_table,'b');
hold on
plot(h_AE,qS_AE,'r');
plot(h_qS_max,qS_max,'bo','MarkerFaceColor','b');
plot(h_qS_AE_max,qS_AE_max,'ro','MarkerFaceColor','r');
title(sprintf('Heat rate per surface unit at stagnation point %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$h \ (km)$$','Interpreter','latex','Fontsize',12)
ylabel('$$q_S \ (kW/m^2)$$','Interpreter','latex','Fontsize',12)
legend('Numerical integration','Allen-Eggers solution')

figure(10)
plot(h_table,pD_table,'b');
hold on
plot(h_AE,pD_AE,'r');
plot(h_pD_max,pD_max,'bo','MarkerFaceColor','b');
plot(h_pD_AE_max,pD_AE_max,'ro','MarkerFaceColor','r');
title(sprintf('Dynamic pressure %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$h \ (km)$$','Interpreter','latex','Fontsize',12)
ylabel('$$p_D \ (Pa)$$','Interpreter','latex','Fontsize',12)
legend('Numerical integration','Allen-Eggers solution')


