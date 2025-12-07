%% HOMEWORK 2

% load variables
load data_vs2.mat


% index intervals
tspan=Qtrue(:,5)';
index=(1:4000);

conv=(3600*180)/pi;


%% TRIAD ALGORITHM

[qdet1,Ptt1]=triad(OBS,    [2 4],   (1:1000));
[qdet2,Ptt2]=prop(qdet1(end,:),Ptt1(:,:,end),OMEGA_MEAS,(1001:2000));
[qdet3,Ptt3]=triad(OBS,    [1 3],  (2001:3000));
[qdet4,Ptt4]=triad(OBS,    [2 4],   (3001:4000));

% concatenation of values
Qdet=[qdet1 ; qdet2 ; qdet3 ; qdet4];
Ptt_det=cat(3,Ptt1,Ptt2,Ptt3,Ptt4);

% clear useless junk
clear qdet1 qdet2 qdet3 qdet4 ...
      Ptt1 Ptt2 Ptt3 Ptt4

% estimation of the delta_q matrix
dq_det=zeros(length(Qtrue),4);
for i=1:length(Qtrue)
    dq_det(i,:)=quat_mult(Qtrue(i,1:4),Qdet(i,:));
end

delta_phi=2*conv*dq_det(:,1);
delta_the=2*conv*dq_det(:,2);
delta_psi=2*conv*dq_det(:,3);

% estimation of the sigma matrix
sigma_det=zeros(length(Qtrue),3);
for i=1:length(Qtrue)
    sigma_det(i,:)=sqrt(diag(Ptt_det(:,:,i)));
end

sigma_phi=sigma_det(:,1)*conv;
sigma_the=sigma_det(:,2)*conv;
sigma_psi=sigma_det(:,3)*conv;


% plot quaternions for the attitude
figure(1)
plot(tspan,Qdet(:,1),'b-');
hold on
plot(tspan,Qdet(:,2),'g-');
plot(tspan,Qdet(:,3),'r-');
plot(tspan,Qdet(:,4),'c-');
hold off
title('TRIAD - Quaternions')
xlabel('Time (s)')
ylabel('Quaternions')
legend('q_1','q_2','q_3','q_4')


% plot q1, q2, q3 errors and uncertainties
figure(2)
sgtitle('TRIAD - Quaternions (errors and uncertainties)')
subplot(3,1,1)
plot(tspan,dq_det(:,1),'b.');
hold on
plot(tspan,0.5*sigma_phi/conv,'r-');
plot(tspan,-0.5*sigma_phi/conv,'y-');
hold off
title('q_1 error')
ylabel('q_1 error')

subplot(3,1,2)
plot(tspan,dq_det(:,2),'b.');
hold on
plot(tspan,0.5*sigma_the/conv,'r-');
plot(tspan,-0.5*sigma_the/conv,'y-');
hold off
title('q_2 error')
ylabel('q_2 error')

subplot(3,1,3)
plot(tspan,dq_det(:,3),'b.');
hold on
plot(tspan,0.5*sigma_psi/conv,'r-');
plot(tspan,-0.5*sigma_psi/conv,'y-');
hold off
title('q_3 error')
ylabel('q_3 error')


% plot roll, pitch, yaw angles errors and uncertainties
figure(3)
sgtitle('TRIAD - Euler angles (errors and uncertainties)')
subplot(3,1,1)
plot(tspan,delta_phi,'b.');
hold on
plot(tspan,sigma_phi,'r-');
plot(tspan,-sigma_phi,'y-');
hold off
title('\Phi error')
ylabel('\Phi error (arcsec)')
axis([0 400 -60 60])

subplot(3,1,2)
plot(tspan,delta_the,'b.');
hold on
plot(tspan,sigma_the,'r-');
plot(tspan,-sigma_the,'y-');
hold off
title('\Theta error')
ylabel('\Theta error (arcsec)')
axis([0 400 -60 60])

subplot(3,1,3)
plot(tspan,delta_psi,'b.');
hold on
plot(tspan,sigma_psi,'r-');
plot(tspan,-sigma_psi,'y-');
hold off
title('\Psi error')
ylabel('\Psi error (arcsec)')
axis([0 400 -60 60])




%% QUEST ALGORITHM

[qopt1,Ptt1]=quest(OBS,    (1:10),     (1:1000));
[qopt2,Ptt2]=prop(qopt1(end,:),Ptt1(:,:,end),OMEGA_MEAS,(1001:2000));
[qopt3,Ptt3]=quest(OBS,    [1:3 10],   (2001:3000));
[qopt4,Ptt4]=quest(OBS,    [1:8 10],   (3001:3074));
[qopt5,Ptt5]=quest(OBS,    (1:10),     (3075:3659));
[qopt6,Ptt6]=quest(OBS,    [1:8 10],   (3660:3675));
[qopt7,Ptt7]=quest(OBS,    (1:10),     (3676:4000));

% concatenation of values
Qopt=[qopt1 ; qopt2 ; qopt3 ; qopt4 ; qopt5 ; qopt6 ; qopt7];
Ptt_opt=cat(3,Ptt1,Ptt2,Ptt3,Ptt4,Ptt5,Ptt6,Ptt7);
Popt=Ptt_opt/4;

% clear useless junk
clear qopt1 qopt2 qopt3 qopt4 qopt5 qopt6 qopt7 ...
      Ptt1 Ptt2 Ptt3 Ptt4 Ptt5 Ptt6 Ptt7

% estimation of the delta_q matrix
dq_opt=zeros(length(Qtrue),4);
for i=1:length(Qtrue)
    dq_opt(i,:)=quat_mult(Qtrue(i,1:4),Qopt(i,:));
end

delta_phi=2*conv*dq_opt(:,1);
delta_the=2*conv*dq_opt(:,2);
delta_psi=2*conv*dq_opt(:,3);

% estimation of the sigma matrix
sigma_opt=zeros(length(Qtrue),3);
for i=1:length(Qtrue)
    sigma_opt(i,:)=sqrt(diag(Ptt_opt(:,:,i)));
end

sigma_phi=sigma_opt(:,1)*conv;
sigma_the=sigma_opt(:,2)*conv;
sigma_psi=sigma_opt(:,3)*conv;


% plot quaternions for the attitude
figure(4)
plot(tspan,Qopt(:,1),'b-');
hold on
plot(tspan,Qopt(:,2),'g-');
plot(tspan,Qopt(:,3),'r-');
plot(tspan,Qopt(:,4),'c-');
hold off
title('QUEST - Quaternions')
xlabel('Time (s)')
ylabel('Quaternions')
legend('q_1','q_2','q_3','q_4')


% plot q1, q2, q3 errors and uncertainties
figure(5)
sgtitle('QUEST - Quaternions (errors and uncertainties)')
subplot(3,1,1)
plot(tspan,dq_opt(:,1),'b.');
hold on
plot(tspan,0.5*sigma_phi/conv,'r-');
plot(tspan,-0.5*sigma_phi/conv,'y-');
hold off
title('q_1 error')
ylabel('q_1 error')

subplot(3,1,2)
plot(tspan,dq_opt(:,2),'b.');
hold on
plot(tspan,0.5*sigma_the/conv,'r-');
plot(tspan,-0.5*sigma_the/conv,'y-');
hold off
title('q_2 error')
ylabel('q_2 error')

subplot(3,1,3)
plot(tspan,dq_opt(:,3),'b.');
hold on
plot(tspan,0.5*sigma_psi/conv,'r-');
plot(tspan,-0.5*sigma_psi/conv,'y-');
hold off
title('q_3 error')
ylabel('q_3 error')


% plot roll, pitch, yaw angles errors and uncertainties
figure(6)
sgtitle('QUEST - Euler angles (errors and uncertainties)')
subplot(3,1,1)
plot(tspan,delta_phi,'b.');
hold on
plot(tspan,sigma_phi,'r-');
plot(tspan,-sigma_phi,'y-');
hold off
title('\Phi error')
ylabel('\Phi error (arcsec)')
axis([0 400 -60 60])

subplot(3,1,2)
plot(tspan,delta_the,'b.');
hold on
plot(tspan,sigma_the,'r-');
plot(tspan,-sigma_the,'y-');
hold off
title('\Theta error')
ylabel('\Theta error (arcsec)')
axis([0 400 -60 60])

subplot(3,1,3)
plot(tspan,delta_psi,'b.');
hold on
plot(tspan,sigma_psi,'r-');
plot(tspan,-sigma_psi,'y-');
hold off
title('\Psi error')
ylabel('\Psi error (arcsec)')
axis([0 400 -60 60])




%% EKF algorithm

q_est=Qopt(1,:)';
b_est=[0; 0 ; 0];
X=[q_est ; b_est];
Xback=X;

dt0=Qtrue(2,5)-Qtrue(1,5);
sigma_u=5*10^-9;
sigma_v=10^-7;

Pback=[Popt(:,:,1) zeros(3) ; zeros(3) (sigma_u*dt0)^2*eye(3)];

% save covariance, quaternion, bias in matrix forma
cov_save(1,:)=diag(Pback);
Qekf(1,:)=q_est';
best_save(1,:)=b_est';

% estimation of Q
Q11=sigma_v^2*eye(3);
Q22=(sigma_u*dt0)^2*eye(3);
Q=[Q11 zeros(3) ; zeros(3) Q22];

dt=0;
for k=2:length(Qtrue)
    t_ini=Qtrue(k-1,5);
    t_fin=Qtrue(k,5);
    
    omega_meas=OMEGA_MEAS(k-1,:)';
    
    
    if k>1000 && k<=2000

        dt=dt+0.1;
        bias=Xback(5:7);
        
        % re-estimation of Q
        Q11=sigma_v^2*eye(3);
        Q22=(sigma_u*dt)^2*eye(3);
        Qspec=[Q11 zeros(3) ; zeros(3) Q22];

        [X,P]=prop2(Xback,Pback,omega_meas,Qspec,t_ini,t_fin,bias);
        
    else
        if k==2001
            Xback(5:7)=[0 ; 0 ; 0];
        end
        
        % values from QUEST are fed to the EKF
        qopt=Qopt(k,:)';
        R=Popt(:,:,k);
    
        % proper EKF algorithm
        [X,P]=EKF(Xback,Pback,qopt,R,omega_meas,Q,t_ini,t_fin);
    
    end
    % update values for next iteration
    Xback=X;
    Pback=P;
    
    % save values in matrix for graphic
    cov_save(k,:)=diag(P);
    Qekf(k,:)=X(1:4)';
    best_save(k,:)=X(5:7)';
end
cov_tt=cov_save*4;
OMEGA_EST=OMEGA_MEAS-best_save;

dq_ekf=zeros(length(Qtrue),4);
for i=1:length(Qtrue)
    dq_ekf(i,:)=quat_mult(Qtrue(i,1:4),Qekf(i,:));
end

delta_phi=2*conv*dq_ekf(:,1);
delta_the=2*conv*dq_ekf(:,2);
delta_psi=2*conv*dq_ekf(:,3);

sigma_ekf=sqrt(cov_tt(:,1:3));
sigma_phi=sigma_ekf(:,1)*conv;
sigma_the=sigma_ekf(:,2)*conv;
sigma_psi=sigma_ekf(:,3)*conv;

% estimate angular velocity errors
domega=OMEGA_TRUE-OMEGA_EST;


% plot quaternions for the attitude
figure(7)
plot(tspan,Qekf(:,1),'b-');
hold on
plot(tspan,Qekf(:,2),'g-');
plot(tspan,Qekf(:,3),'r-');
plot(tspan,Qekf(:,4),'c-');
hold off
title('EKF - Quaternions')
xlabel('Time (s)')
ylabel('Quaternions')
legend('q_1','q_2','q_3','q_4')


% plot q1, q2, q3 errors and uncertainties
figure(8)
sgtitle('EKF - Quaternions (errors and uncertainties)')
subplot(3,1,1)
plot(tspan,dq_ekf(:,1),'b.');
hold on
plot(tspan,0.5*sigma_phi/conv,'r-');
plot(tspan,-0.5*sigma_phi/conv,'y-');
hold off
title('q_1 error')
ylabel('q_1 error')

subplot(3,1,2)
plot(tspan,dq_ekf(:,2),'b.');
hold on
plot(tspan,0.5*sigma_the/conv,'r-');
plot(tspan,-0.5*sigma_the/conv,'y-');
hold off
title('q_2 error')
ylabel('q_2 error')

subplot(3,1,3)
plot(tspan,dq_ekf(:,3),'b.');
hold on
plot(tspan,0.5*sigma_psi/conv,'r-');
plot(tspan,-0.5*sigma_psi/conv,'y-');
hold off
title('q_3 error')
ylabel('q_3 error')


% plot roll, pitch, yaw angles errors and uncertainties
figure(9)
sgtitle('EKF - Euler angles (errors and uncertainties)')
subplot(3,1,1)
plot(tspan,delta_phi,'b.');
hold on
plot(tspan,sigma_phi,'r-');
plot(tspan,-sigma_phi,'y-');
hold off
title('\Phi error')
ylabel('\Phi error (arcsec)')
axis([0 400 -60 60])

subplot(3,1,2)
plot(tspan,delta_the,'b.');
hold on
plot(tspan,sigma_the,'r-');
plot(tspan,-sigma_the,'y-');
hold off
title('\Theta error')
ylabel('\Theta error (arcsec)')
axis([0 400 -60 60])

subplot(3,1,3)
plot(tspan,delta_psi,'b.');
hold on
plot(tspan,sigma_psi,'r-');
plot(tspan,-sigma_psi,'y-');
hold off
title('\Psi error')
ylabel('\Psi error (arcsec)')
axis([0 400 -60 60])


% plot angular velocities errors
figure(10)
sgtitle('Angular velocity error')
subplot(3,1,1)
plot(tspan,domega(:,1)*conv,'b-');
title('\omega_1 error')
ylabel('\omega_1 error (arcsec/s)')

subplot(3,1,2)
plot(tspan,domega(:,2)*conv,'g-');
title('\omega_2 error')
ylabel('\omega_2 error (arcsec/s)')

subplot(3,1,3)
plot(tspan,domega(:,3)*conv,'r-');
title('\omega_3 error')
ylabel('\omega_3 error (arcsec/s)')


% plot estimated bias
figure(11)
sgtitle('Estimated bias \beta')
subplot(3,1,1)
plot(tspan,best_save(:,1)*conv,'b-');
title('\beta_1')
ylabel('\beta_1 (arcsec/s)')

subplot(3,1,2)
plot(tspan,best_save(:,2)*conv,'g-');
title('\beta_2')
ylabel('\beta_2 (arcsec/s)')

subplot(3,1,3)
plot(tspan,best_save(:,3)*conv,'r-');
title('\beta_2')
ylabel('\beta_3 (arcsec/s)')

