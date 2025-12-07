function J=J_function_3arc(x_guess)
% function that accepts the unknown quantities (alpha1,t1,alpha2,t2,alpha3,t3) 
% and returns the value of magnitude of the boundary condition violations

global RE RM mu

alpha1=x_guess(1);
t1=x_guess(2);
alpha2=x_guess(3);
t2=x_guess(4);
alpha3=x_guess(5);
t3=x_guess(6);



%% INITIALIZATION

% initial values of the state
r_0=RE;
theta_0=0;
vr_0=0;
vtheta_0=sqrt(mu/r_0);



%% INTEGRATION

options_int = odeset('RelTol',1e-12,'AbsTol',1e-12);

% first arc
CI1=[r_0,theta_0,vr_0,vtheta_0,alpha1];
[~,W1]=ode45('state_eqs',[0 t1],CI1,options_int);

% first arc
CI2=[W1(end,1:4),alpha2];
[~,W2]=ode45('state_eqs',[t1 (t1+t2)],CI2,options_int);

% second arc
CI3=[W2(end,1:4),alpha3];
[~,W3]=ode45('state_eqs',[(t1+t2) (t1+t2+t3)],CI3,options_int);

r3=W3(:,1);
vr3=W3(:,3);
vtheta3=W3(:,4);



%% ESTIMATION OF COST FUNCTION

% conditions for the final values of the state
r_f=RM;
vr_f=0;
vtheta_f=sqrt(mu/r_f);

% estimation of the error
delta_r=r3(end)-r_f;
delta_vr=vr3(end)-vr_f;
delta_vtheta=vtheta3(end)-vtheta_f;

% weights
c1=1;
c2=1e6;
c3=1e6;

% cost function
J=abs((t1+t2+t3)/(24*3600))+c1*abs(delta_r)+c2*abs(delta_vr)+c3*abs(delta_vtheta);


end

