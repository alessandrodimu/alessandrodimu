function J=J_function(x_guess)
% function that accepts the unknown quantities (lam10,lam30,lam40,tf) and
% returns the value of magnitude of the boundary condition violations

global r_ini v_ini r_fin v_fin

lambda_10=x_guess(1);
lambda_30=x_guess(2);
lambda_40=x_guess(3);
t_f=x_guess(4);



%% INITIALIZATION

r_0=r_ini;
csi_0=0;
vr_0=0;
vt_0=v_ini;

x0=[r_0,csi_0,vr_0,vt_0,lambda_10,lambda_30,lambda_40];



%% INTEGRATION

tspan=[0 t_f];
options_int = odeset('RelTol',1e-10,'AbsTol',1e-8);

[~,W]=ode45(@dyn_kin_and_adjoints,tspan,x0,options_int);

r_f=W(end,1);
csi_f=W(end,2);
vr_f=W(end,3);
vt_f=W(end,4);



%% ESTIMATION OF COST FUNCTION

psi=[r_0-r_ini ; csi_0 ; vr_0 ; vt_0-v_ini ;...
     r_f-r_fin ; vr_f ; vt_f-v_fin];

J=norm(psi);


end

