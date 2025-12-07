function [Xprop,Pprop]=prop2(Xback,Pback,omega_meas,Q,t_ini,t_fin,bias)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% options for the integration/propagation
Tol0=1e-10; 
Tol1=1e-10;
options=odeset('RelTol',Tol0,'AbsTol',Tol1);

q_est=Xback(1:4);
b_est=Xback(5:7);

omega_est=omega_meas-b_est;

% quaternion calculation
ic_quat=q_est;
[~,W0]=ode45(@(t,w) model_quat(t,w,omega_est),[t_ini t_fin],ic_quat,options);
qprop=W0(end,1:4)';

% matrix theta calculation
ic_theta=reshape(eye(3),1,9);
[~,W1]=ode45(@(t,w) model_the(t,w,omega_est),[t_ini t_fin],ic_theta,options);
THETA=reshape(W1(end,:),3,3);

% matrix psi calculation
ic_psi=reshape(zeros(3),1,9);
[~,W2]=ode45(@(t,w) model_psi(t,w,omega_est),[t_ini t_fin],ic_psi,options);
PSI=reshape(W2(end,:),3,3);

% matrix phi calculation
PHI=[THETA PSI ; zeros(3) eye(3)];

% output of propagation
Pprop=PHI*Pback*PHI'+Q;
Xprop=[qprop ; bias];