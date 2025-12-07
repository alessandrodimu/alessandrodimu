function [X,P]=EKF(Xback,Pback,qopt,R,omega_meas,Q,t_ini,t_fin)
% EKF algorithm

% options for the integration/propagation
Tol0=1e-10; 
Tol1=1e-10;
options=odeset('RelTol',Tol0,'AbsTol',Tol1);

% state vector inizialization
q_est=Xback(1:4);
b_est=Xback(5:7);

P_tilde_back=Pback;
omega_est=omega_meas-b_est;

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

% covariance matrix updated
P_est=PHI*P_tilde_back*PHI'+Q;

% state updated
ic_quat=q_est;
[~,W3]=ode45(@(t,w) model_quat(t,w,omega_est),[t_ini t_fin],ic_quat,options);
q_est=W3(end,:)';

z=quat_mult(qopt',q_est');

% sensitivity matrix H
H=[eye(3) zeros(3)];

% kalman gain matrix K
K=P_est*H'*inv(H*P_est*H'+R);

% innovation vector
delta_x=K*z(1:3);
delta_q=[delta_x(1:3) ; 1];
delta_q=delta_q/norm(delta_q);
delta_b=delta_x(4:6);

% state update
q_est(1:3)=-q_est(1:3);
q_est=quat_mult(delta_q',q_est');
b_est=b_est+delta_b;

X=[q_est ; b_est];

% covariance matrix update
P=(eye(6)-K*H)*P_est;

end

