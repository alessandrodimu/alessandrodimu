function J=J_function_cont(x_guess)
% function that accepts the unknown quantities (pr_0,pvr_0,pvtheta_0,tf) 
% and returns the value of magnitude of the boundary condition violations

global RE RM mu

pr_0=x_guess(1);
pvr_0=x_guess(2);
pvtheta_0=x_guess(3);
t_f=x_guess(4);



%% INITIALIZATION

% initial values of the state
r_0=RE;
theta_0=0;
vr_0=0;
vtheta_0=sqrt(mu/r_0);
K_0=0;

CI=[r_0,theta_0,vr_0,vtheta_0,K_0,pr_0,pvr_0,pvtheta_0];



%% INTEGRATION

tspan=[0 t_f];
options_int = odeset('RelTol',1e-12,'AbsTol',1e-12);

[~,W]=ode45('hamiltonian_eqs',tspan,CI,options_int);

r=W(:,1);
vr=W(:,3);
vtheta=W(:,4);



%% ESTIMATION OF COST FUNCTION

% conditions for the final values of the state
r_f=RM;
vr_f=0;
vtheta_f=sqrt(mu/r_f);

% estimation of the error
delta_r=r(end)-r_f;
delta_vr=vr(end)-vr_f;
delta_vtheta=vtheta(end)-vtheta_f;

% weights
c1=1;
c2=1e6;
c3=1e6;

% cost function
J=abs(t_f/(24*3600))+c1*abs(delta_r)+c2*abs(delta_vr)+c3*abs(delta_vtheta);


end

