function [qprop,Pprop]=prop(q_ini,P_ini,OMEGA,indexspan)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% options for the integration/propagation
Tol0=1e-10; 
Tol1=1e-10;
options=odeset('RelTol',Tol0,'AbsTol',Tol1);

sigma_u=5*10^-9;
sigma_v=10^-7;
supp1=reshape(eye(3),1,9);

qprop=zeros(length(indexspan),4);
Pprop=zeros(3,3,length(indexspan));

for i=indexspan
    omega=OMEGA(i-1,:)';
    [~,W]=ode45(@(t,w) model_prop(t,w,omega),[0 0.1],[q_ini supp1],options);
    qprop(i,:)=W(end,1:4);
    phi=reshape(W(end,5:13),3,3);
    Pprop(:,:,i)=phi*P_ini*phi';
    q_ini=qprop(i,:);
    P_ini=Pprop(:,:,i);
end
qprop=qprop(indexspan,:);
Pprop=Pprop(:,:,indexspan);

