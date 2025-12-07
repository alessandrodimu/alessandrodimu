% HOMEWORK 2

%% WIPEOUT
clear
close all
clc



%% DATA
global mu L2

in=input('Please select the system: (1-2)\n 1) Earth-Moon system\n 2) Sun-Earth system\n\n');
if isscalar(in)==0 || isreal(in)==0 || in<1 || in>2 %#ok<COMPNOT>
    error('Expected input must be either 1 or 2')
end

if in==2
    in2=input('\nPlease select the Lissajous orbit amplitude: (1-2)\n 1) 35000 km\n 2) 350000 km\n\n');
    if isscalar(in2)==0 || isreal(in2)==0 || in2<1 || in2>2 %#ok<COMPNOT>
        error('Expected input must be either 1 or 2')
    end
end

% constants
G=6.674*10^(-20);                       % in km^3/kg*s^2
mE=5.972*10^24;                         % in kg
mM=7.342*10^22;                         % in kg
mS=1.9885*10^30;                        % in kg

switch in
    case 1
        mu=mM/(mE+mM);                  % adimensional
        D=384400;                       % in km
        w_syn=sqrt(G*(mE+mM)/D^3);      % in rad/s
        TU=1/w_syn;                     % adimensional
    case 2
        mu=mE/(mS+mE);                  % adimensional
        D=149600000;                    % in km
        w_syn=sqrt(G*(mS+mE)/D^3);      % in rad/s
        TU=1/w_syn;                     % adimensional
end



%% FINDING L2 POSITION
% we find the position of L2 from imposing dU/dX=0 with y=0
% then we consider the case (X+mu)>1, and we change reference system
% accordingly, thus giving us a 5th degree polynomial
coeff=[1,3-mu,3-2*mu,-mu,-2*mu,-mu];
XL=roots(coeff);
L2=XL(3)+1-mu;                          % only the third root is real



%% LINEARIZED SYSTEM (xdot=Ax+Bu)
% x=6x1 state vector
% x(1,2,3) are position values, x(4,5,6) are velocity values
% u=3x1 control vector

% parameters
sigma=(1-mu)/(L2+mu)^3+mu/(L2+mu-1)^3;
Uxx=2*sigma+1;
Uyy=1-sigma;
Uzz=-sigma;

% matrices A and B
A=zeros(6,6);
A(1:3,1:3)=zeros(3,3);
A(1:3,4:6)=eye(3,3);
A(4:6,1:3)=diag([Uxx,Uyy,Uzz])*eye(3,3);
A(4:6,4:6)=[0 2 0; -2 0 0; 0 0 0];

B=zeros(6,3);
B(4:6,:)=eye(3,3);



%% REFERENCE SOLUTION
% in-plane and out-of-plane frequencies
lambda=eig(A);
w_xy=norm(lambda(3));
w_z=norm(lambda(5));

% parameter
Am=(w_xy^2+Uxx)/(2*w_xy);

% initial conditions for the Lissajous traj
switch in
    case 1
        x0=0;
        y0=3500/D;
        z0=0;
        vx0=w_xy*y0/Am;
        vy0=-Am*w_xy*x0;
        vz0=-y0*w_z;
        X0=[x0 y0 z0 vx0 vy0 vz0]';
        
    case 2
        x0=0;
        switch in2
            case 1
                y0=35000/D;
            case 2
                y0=350000/D;
        end
        z0=0;
        vx0=w_xy*y0/Am;
        vy0=-Am*w_xy*x0;
        vz0=-y0*w_z;
        X0=[x0 y0 z0 vx0 vy0 vz0]';
end

% integration time is 1 year, each step is 3600s
dt=60*60/TU;
tspan=(0:dt*TU:365*86400)/TU;

% reference solution (LJ traj)
xr=y0/Am*sin(w_xy*tspan);
yr=y0*cos(w_xy*tspan);
zr=-y0*sin(w_z*tspan);
vxr=y0*w_xy/Am*cos(w_xy*tspan);
vyr=-y0*w_xy*sin(w_xy*tspan);
vzr=-y0*w_z*cos(w_z*tspan);

Xr=[xr; yr ; zr ; vxr ; vyr ; vzr];



%% LQR MATRICES
% Q=6x6 weight matrix on the state
% R=3x3 weight matrix on the control

Q=zeros(6,6);
Q(1:3,1:3)=1e3*eye(3,3);
Q(4:6,4:6)=eye(3,3);
R=1e-3*eye(3,3);

% solving Riccati algebraic equation and finding optimal coeff for control
S=icare(A,B,Q,R);
K=R\B'*S;



%% NON LINEAR SYSTEM
% initialization
Xnl=zeros(6,length(tspan));
Xnl(:,1)=X0;
U=zeros(3,length(tspan));
deltaV=0;

X=X0;
options=odeset('RelTol',1e-10,'AbsTol',1e-10);

for i = 1:length(tspan)-1
    U(:,i)=-K*(X-Xr(:,i));
    deltaV=deltaV+(1000*D*dt/TU)*norm(U(:,i));
    [t,W]=ode45(@state_eqs,[tspan(i) tspan(i+1)],X,options,U(:,i));
    Xnl(:,i+1)=W(end,:);
    X=Xnl(:,i+1);
end

U(:,end)=-K*(Xnl(:,end)-Xr(:,end));

% printing some values
fprintf('\nAmplitude = %d km',y0*D);
fprintf('\nDeltaV =%7.2f m/s^2 \n\n',deltaV);



%% PLOTS
% transfer from m2-centric ref system to synodic system
Xnl(1,:)=Xnl(1,:)+L2;
Xr(1,:)=Xr(1,:)+L2;

tf=365;
switch in
    case 1
        str='(E-M system)';
    case 2
        str='(S-E system)';
end

figure(1)
hold on
grid on
axis('equal')
view(-30,30)
p1=plot3(D*Xr(1,:),D*Xr(2,:),D*Xr(3,:),'r');
p2=plot3(D*Xnl(1,:),D*Xnl(2,:),D*Xnl(3,:),'b');
plot3(D*L2,0,0,'ok')
title(sprintf('Lissajous orbit %s',str),'Interpreter','latex','Fontsize',12)
xlabel('$$X \ (km)$$','Interpreter','latex','Fontsize',12)
ylabel('$$Y \ (km)$$','Interpreter','latex','Fontsize',12)
zlabel('$$Z \ (km)$$','Interpreter','latex','Fontsize',12)
legend([p1 p2],{'Reference trajectory','Controlled trajectory'},'Interpreter','latex')

figure(2)
hold on
subplot(3,1,1)
plot(TU/86400*tspan,D*Xnl(1,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$x \ (km)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
subplot(3,1,2)
plot(TU/86400*tspan,D*Xnl(2,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$y \ (km)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
subplot(3,1,3)
plot(TU/86400*tspan,D*Xnl(3,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$z \ (km)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
sgtitle(sprintf('Position %s',str),'Interpreter','latex','Fontsize',12)

figure(3)
hold on
subplot(3,1,1)
plot(TU/86400*tspan,1000*D/TU*Xnl(4,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$v_x \ (m/s)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
subplot(3,1,2)
plot(TU/86400*tspan,1000*D/TU*Xnl(5,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$v_y \ (m/s)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
subplot(3,1,3)
plot(TU/86400*tspan,1000*D/TU*Xnl(6,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$v_z \ (m/s)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
sgtitle(sprintf('Velocity %s',str),'Interpreter','latex','Fontsize',12)

figure(4)
a=(1000*D/TU^2)*U;
hold on
subplot(3,1,1)
plot(TU/86400*tspan,a(1,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$u_x \ (m/s^2)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
subplot(3,1,2)
plot(TU/86400*tspan,a(2,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$u_y \ (m/s^2)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
subplot(3,1,3)
plot(TU/86400*tspan,a(3,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$u_z \ (m/s^2)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
sgtitle(sprintf('Control acceleration %s',str),'Interpreter','latex','Fontsize',12)

figure(5)
err=Xnl-Xr;
hold on
subplot(3,1,1)
plot(TU/86400*tspan,D*err(1,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$\Delta x \ (km)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
subplot(3,1,2)
plot(TU/86400*tspan,D*err(2,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$\Delta y \ (km)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
subplot(3,1,3)
plot(TU/86400*tspan,D*err(3,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$\Delta z \ (km)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
sgtitle(sprintf('Error on position %s',str),'Interpreter','latex','Fontsize',12)

figure(6)
hold on
subplot(3,1,1)
plot(TU/86400*tspan,D*1000/TU*err(4,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$\Delta v_x \ (m/s)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
subplot(3,1,2)
plot(TU/86400*tspan,D*1000/TU*err(5,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$\Delta v_y \ (m/s)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
subplot(3,1,3)
plot(TU/86400*tspan,D*1000/TU*err(6,:),'b')
xlabel('$$t \ (d)$$','FontSize',10,'Interpreter','latex')
ylabel('$$\Delta v_z \ (m/s)$$','FontSize',10,'Interpreter','latex')
xlim([0 tf])
sgtitle(sprintf('Error on velocity %s',str),'Interpreter','latex','Fontsize',12)


