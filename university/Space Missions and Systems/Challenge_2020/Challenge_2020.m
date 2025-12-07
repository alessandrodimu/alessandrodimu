clear all
close all
clc
format long e

%% parte 1
% Observables Loading

obs_file_1=load('optical_ganymede.txt');

t_gan=obs_file_1(:,1);
g1_obs=obs_file_1(:,2);
g2_obs=obs_file_1(:,3);

N_points=length(t_gan);


% Check observables

figure()
% Plot G1 and G2
hold on
plot(t_gan,g1_obs,'k-')
plot(t_gan,g2_obs,'k--')
hold off

title('Ganymede observables')
legend('$\hat{G}_{x}$','$\hat{G}_{y}$','Interpreter','latex');
xlabel('Time [s]')
ylabel('X,Y components [km]')


% Problem constants

global GM_jup GM_eur ni0_eur ni0_gan a_eur a_gan w_eur w_gan

ni0_eur=deg2rad(-44);
ni0_gan=deg2rad(-13);
a_eur=670900;
a_gan=1070400;
GM_jup=126686534.20;
GM_eur=3202.74;

w_eur=sqrt(GM_jup/a_eur^3);
w_gan=sqrt(GM_jup/a_gan^3);


% A priori information

x0=1180486;
y0=0.6;
vx0=0;
vy0=8.330;
GM0=9887;

oldstate0=[x0;y0;vx0;vy0;GM0];
state0=[x0;y0;vx0;vy0;GM0];
dev0=[0;0;0;0;0];  %lavoro con questa 

Pcov0=diag([(1)^2,(1)^2,(0.001)^2,(0.001)^2,(1)^2]);
%L0=inv(Pcov0);  % lavoro con questa
L0=zeros(5,5);
L0(1,1)=1/(1)^2;
L0(2,2)=1/(1)^2;
L0(3,3)=1/(0.001)^2;
L0(4,4)=1/(0.001)^2;
L0(5,5)=1/(1)^2; 


% Trajectory propagation options

Tol0=1e-13; 
Tol1=1e-13;
options=odeset('RelTol',Tol0,'AbsTol',Tol1);

supp1=reshape(eye(5),1,25);
tspan=[0,t_gan'];

iterations=6;


% Run the filter

W1_x_save=zeros(N_points,iterations);
W1_y_save=zeros(N_points,iterations);
W1_vx_save=zeros(N_points,iterations);
W1_vy_save=zeros(N_points,iterations);
W1_GM_gan_save=zeros(N_points,iterations);

% Select weights
w_g1=1000;
w_g2=1000;

for counter=1:iterations
    
    [t_gan,W1]=ode113('Model_and_transition_ch',tspan,[state0' supp1],options);
    
    W1=W1(2:end,:);
    t_gan=t_gan(2:end);
    
    W1_x_save(:,counter)=W1(:,1);
    W1_y_save(:,counter)=W1(:,2);
    W1_vx_save(:,counter)=W1(:,3);
    W1_vy_save(:,counter)=W1(:,4);
    W1_GM_gan_save(:,counter)=W1(:,5);
    
    % initialization of parameters I will need later on
    x_eur_comp=zeros(N_points,1);
    y_eur_comp=zeros(N_points,1);
    x_gan_comp=zeros(N_points,1);
    y_gan_comp=zeros(N_points,1);
    
    direzione_rel_x=zeros(N_points,1);
    direzione_rel_y=zeros(N_points,1);
    
    y_g1=zeros(N_points,1);
    y_g2=zeros(N_points,1);
    
    H=[];
        
    for i=1:N_points
        
        x=W1(i,1);
        y=W1(i,2);
        vx=W1(i,3);
        vy=W1(i,4);
        GM_gan=W1(i,5);
        X=[x y vx vy GM_gan];
        
        t=t_gan(i);
        
        % Estimated computed values
        x_eur_comp(i)=a_eur*cos(ni0_eur+w_eur*t);
        y_eur_comp(i)=a_eur*sin(ni0_eur+w_eur*t);
        x_gan_comp(i)=a_gan*cos(ni0_gan+w_gan*t);
        y_gan_comp(i)=a_gan*sin(ni0_gan+w_gan*t);
        
        supp2=sqrt((x_gan_comp(i)-x)^2+(y_gan_comp(i)-y)^2);
        direzione_rel_x(i)=(x_gan_comp(i)-x)/supp2;
        direzione_rel_y(i)=(y_gan_comp(i)-y)/supp2;
        
        % Build H matrix   
        Htilde=H_tilde_ch(X,t);
        phi=reshape(W1(i,6:30),5,5);
        H=[H ; Htilde*phi];
        
        % Compute residuals
        y_g1(i)=g1_obs(i)-direzione_rel_x(i);
        y_g2(i)=g2_obs(i)-direzione_rel_y(i);
        
    end
    
    % Residual statistics
    MEAN_g1=mean(y_g1);
    MEAN_g2=mean(y_g2);
    RMS_g1=sqrt(y_g1'*y_g1/length(y_g1));
    RMS_g2=sqrt(y_g2'*y_g2/length(y_g2));
    SOS_g1=y_g1'*w_g1*y_g1;
    SOS_g2=y_g2'*w_g2*y_g2;
    
    % Update weights, if necessary
    if counter>3
        w_g1=1/(RMS_g1)^2;
        w_g2=1/(RMS_g2)^2;
    end
    
    % Build weighting matrix
    for i=1:N_points
        if i>=26 && i<=266
            W(2*i-1,2*i-1)=w_g1;
            W(2*i,2*i)=w_g2;
        else
            W(2*i-1,2*i-1)=w_g1;
            W(2*i,2*i)=w_g2;
        end
    end
    
    % Build residual vector
    y=[];
    for i=1:N_points
        y=[y; y_g1(i)];
        y=[y; y_g2(i)];
    end
    
    % Compute L and N matrices
    L=L0+H'*W*H;
    N=L0*dev0+H'*W*y;
    
    % Compute solution
    dev_est=L\N;
    Pcov_est=inv(L);
    
    % Update Solution
    state_update=state0+dev_est;
    state0=state_update;
    dev0=dev0-dev_est;   % shift to be wrt the new state
    
    % Draw Figures
    figure()
    subplot(2,1,1)
    string1=sprintf('G1 residuals: Iter %i\n Mean: %1.2e RMS: %1.3e SOS: %1.3e',[counter,MEAN_g1,RMS_g1,SOS_g1]);
    plot(t_gan,y_g1,'b+')
    xlabel('Time [s]')
    ylabel('[km]')
    title(string1)
    
    subplot(2,1,2)
    string2=sprintf('G2 residuals: Iter %i\n Mean: %1.2e RMS: %1.3e SOS: %1.3e',[counter,MEAN_g2,RMS_g2,SOS_g2]);
    plot(t_gan,y_g2,'r+')
    xlabel('Time [s]')
    ylabel('[km/s]')
    title(string2)
    
end

L=[];
N=[];

%% parte 2

tspan2=linspace(0,309600,5000);
[tspan2,W2]=ode113('Model_and_transition_ch',tspan2,[state0' supp1],options);

x_eur=a_eur*cos(ni0_eur+w_eur*tspan2);
y_eur=a_eur*sin(ni0_eur+w_eur*tspan2);
x_sc=W2(:,1);
y_sc=W2(:,2);

modE=sqrt((x_eur-x_sc).^2+(y_eur-y_sc).^2);
sup=min(modE);
ind=find(modE==sup);

phi=reshape(W2(ind,6:30),5,5);
Pcov_minE=phi*Pcov_est*phi';

figure(iterations+2)
plot(tspan2,modE,'k-');
xlabel('Time [s]')
ylabel('Distance to Europa [km]')
title('Approach to Europa after first observation')

grad=[(x_sc(ind)-x_eur(ind))/sup ; (y_sc(ind)-y_eur(ind))/sup];
sigma_modE=sqrt(grad'*Pcov_minE(1:2,1:2)*grad);


%% parte 3

tspan3=linspace(0,115200,5000);

% before observation
[tspan3,W3a]=ode113('Model_and_transition_ch',tspan3,[oldstate0' supp1],options);

x_gan=a_gan*cos(ni0_gan+w_gan*tspan3);
y_gan=a_gan*sin(ni0_gan+w_gan*tspan3);
x_sc3=W3a(:,1);
y_sc3=W3a(:,2);

modGa=sqrt((x_gan-x_sc3).^2+(y_gan-y_sc3).^2);
sup3a=min(modGa);
ind3=find(modGa==sup3a);

phi=reshape(W3a(ind3,6:30),5,5);
Pcov_minGa=phi*Pcov0*phi';

figure(iterations+3)
plot(tspan3,modGa,'k-');
xlabel('Time [s]')
ylabel('Distance to Ganymede [km]')
title('Approach to Ganymede before any observation')

grad=[(x_sc3(ind3)-x_gan(ind3))/sup3a ; (y_sc3(ind3)-y_gan(ind3))/sup3a];
sigma_modGa=sqrt(grad'*Pcov_minGa(1:2,1:2)*grad);

[Va,Da]=eig(Pcov_minGa(1:2,1:2));
axis_a=sqrt(Da(1,1));
axis_b=sqrt(Da(2,2));
eig_v1=Va(:,1);
eig_v2=Va(:,2);

if eig_v1(2)>0
    theta=acos(eig_v1(1));
elseif eig_v1(2)<0
    theta=2*pi-acos(eig_v1(1));
end
figure(iterations+4)
hold on
draw_ellipse(axis_a,axis_b,theta,iterations+4);
draw_ellipse(2*axis_a,2*axis_b,theta,iterations+4);
draw_ellipse(3*axis_a,3*axis_b,theta,iterations+4);
hold off
xlabel('[km]')
ylabel('[km]')
title('Probability ellipsoid at Ganymede CA (B.O.)')

XCA_gan_before=[x_sc3(ind3) ; y_sc3(ind3) ; x_gan(ind3) ; y_gan(ind3) ; axis_a ; axis_b ; theta];


% after observation
[tspan3,W3b]=ode113('Model_and_transition_ch',tspan3,[state0' supp1],options);

x_gan=a_gan*cos(ni0_gan+w_gan*tspan3);
y_gan=a_gan*sin(ni0_gan+w_gan*tspan3);
x_sc3=W3b(:,1);
y_sc3=W3b(:,2);

modGb=sqrt((x_gan-x_sc3).^2+(y_gan-y_sc3).^2);
sup3b=min(modGb);
ind3=find(modGb==sup3b);

phi=reshape(W3b(ind3,6:30),5,5);
Pcov_minGb=phi*Pcov_est*phi';

figure(iterations+5)
plot(tspan3,modGb,'k-');
xlabel('Time [s]')
ylabel('Distance to Ganymede [km]')
title('Approach to Ganymede after first observation')

grad=[(x_sc3(ind3)-x_gan(ind3))/sup3b ; (y_sc3(ind3)-y_gan(ind3))/sup3b];
sigma_modGb=sqrt(grad'*Pcov_minGb(1:2,1:2)*grad);

[Vb,Db]=eig(Pcov_minGb(1:2,1:2));
axis_a=sqrt(Db(1,1));
axis_b=sqrt(Db(2,2));
eig_v1=Vb(:,1);
eig_v2=Vb(:,2);

if eig_v1(2)>0
    theta=acos(eig_v1(1));
elseif eig_v1(2)<0
    theta=2*pi-acos(eig_v1(1));
end
figure(iterations+6)
hold on
draw_ellipse(axis_a,axis_b,theta,iterations+6);
draw_ellipse(2*axis_a,2*axis_b,theta,iterations+6);
draw_ellipse(3*axis_a,3*axis_b,theta,iterations+6);
hold off
xlabel('[km]')
ylabel('[km]')
title('Probability ellipsoid at Ganymede CA (A.O.)')


%% parte 4
% Observables Loading

obs_file_2=load('optical_europa.txt');

t_eur=obs_file_2(:,1);
e1_obs=obs_file_2(:,2);
e2_obs=obs_file_2(:,3);

N_points4=length(t_eur);


% Check observables

figure()
% Plot E1 and E2
hold on
plot(t_eur,e1_obs,'k-')
plot(t_eur,e2_obs,'k--')
hold off

title('Europa observables')
legend('$\hat{E}_{x}$','$\hat{E}_{y}$','Interpreter','latex');
xlabel('Time [s]')
ylabel('X,Y components [km]')

% A priori information

state4=state0;
dev4=[0;0;0;0;0];  %lavoro con questa 

Pcov4=Pcov_est;
L4=inv(Pcov4);  % lavoro con questa


% Trajectory propagation options

tspan4=[0,t_eur'];



% Run the filter

W4_x_save=zeros(N_points4,iterations);
W4_y_save=zeros(N_points4,iterations);
W4_vx_save=zeros(N_points4,iterations);
W4_vy_save=zeros(N_points4,iterations);
W4_GM_gan_save=zeros(N_points4,iterations);

% Select weights
w_e1=1000;
w_e2=1000;

for counter=1:iterations
    
    [t_eur,W4]=ode113('Model_and_transition_ch',tspan4,[state4' supp1],options);
    
    W4=W4(2:end,:);
    t_eur=t_eur(2:end);
    
    W4_x_save(:,counter)=W4(:,1);
    W4_y_save(:,counter)=W4(:,2);
    W4_vx_save(:,counter)=W4(:,3);
    W4_vy_save(:,counter)=W4(:,4);
    W4_GM_gan_save(:,counter)=W4(:,5);
    
    % initialization of parameters I will need later on
    x_eur_comp=zeros(N_points4,1);
    y_eur_comp=zeros(N_points4,1);
    x_gan_comp=zeros(N_points4,1);
    y_gan_comp=zeros(N_points4,1);
    
    direzione_rel_x=zeros(N_points4,1);
    direzione_rel_y=zeros(N_points4,1);
    
    y_e1=zeros(N_points4,1);
    y_e2=zeros(N_points4,1);
    
    H4=[];
        
    for i=1:N_points4
        
        x=W4(i,1);
        y=W4(i,2);
        vx=W4(i,3);
        vy=W4(i,4);
        GM_gan=W4(i,5);
        X=[x y vx vy GM_gan];
        
        t=t_eur(i);
        
        % Estimated computed values
        x_eur_comp(i)=a_eur*cos(ni0_eur+w_eur*t);
        y_eur_comp(i)=a_eur*sin(ni0_eur+w_eur*t);
        x_gan_comp(i)=a_gan*cos(ni0_gan+w_gan*t);
        y_gan_comp(i)=a_gan*sin(ni0_gan+w_gan*t);
        
        supp2=sqrt((x_eur_comp(i)-x)^2+(y_eur_comp(i)-y)^2);
        direzione_rel_x(i)=(x_eur_comp(i)-x)/supp2;
        direzione_rel_y(i)=(y_eur_comp(i)-y)/supp2;
        
        % Build H matrix   
        Htilde4=H_tilde_ch4(X,t);
        phi=reshape(W4(i,6:30),5,5);
        H4=[H4 ; Htilde4*phi];
        
        % Compute residuals
        y_e1(i)=e1_obs(i)-direzione_rel_x(i);
        y_e2(i)=e2_obs(i)-direzione_rel_y(i);
        
    end
    
    % Residual statistics
    MEAN_e1=mean(y_e1);
    MEAN_e2=mean(y_e2);
    RMS_e1=sqrt(y_e1'*y_e1/length(y_e1));
    RMS_e2=sqrt(y_e2'*y_e2/length(y_e2));
    SOS_e1=y_e1'*w_e1*y_e1;
    SOS_e2=y_e2'*w_e2*y_e2;
    
    % Update weights, if necessary
    if counter>3
        w_e1=1/(RMS_e1)^2;
        w_e2=1/(RMS_e2)^2;
    end
    
    % Build weighting matrix
    for i=1:N_points4
        if i>=26 && i<=266
            W(2*i-1,2*i-1)=w_e1;
            W(2*i,2*i)=w_e2;
        else
            W(2*i-1,2*i-1)=w_e1;
            W(2*i,2*i)=w_e2;
        end
    end
    
    % Build residual vector
    y4=[];
    for i=1:N_points4
        y4=[y4; y_e1(i)];
        y4=[y4; y_e2(i)];
    end
    
    % Compute L and N matrices
    L=L4+H4'*W*H4;
    N=L4*dev4+H4'*W*y4;
    
    % Compute solution
    dev_est4=L\N;
    Pcov_est4=inv(L);
    
    % Update Solution
    state_update4=state4+dev_est4;
    state4=state_update4;
    dev4=dev4-dev_est4;   % shift to be wrt the new state
    
    % Draw Figures
    figure()
    subplot(2,1,1)
    string1=sprintf('E1 residuals: Iter %i\n Mean: %1.2e RMS: %1.3e SOS: %1.3e',[counter,MEAN_e1,RMS_e1,SOS_e1]);
    plot(t_eur,y_e1,'b+')
    xlabel('Time [s]')
    ylabel('[km]')
    title(string1)
    
    subplot(2,1,2)
    string2=sprintf('E2 residuals: Iter %i\n Mean: %1.2e RMS: %1.3e SOS: %1.3e',[counter,MEAN_e2,RMS_e2,SOS_e2]);
    plot(t_eur,y_e2,'r+')
    xlabel('Time [s]')
    ylabel('[km/s]')
    title(string2)
    
end


%% parte 5

tspan2=linspace(0,309600,5000);
[tspan2,W5]=ode113('Model_and_transition_ch',tspan2,[state4' supp1],options);

x_eur=a_eur*cos(ni0_eur+w_eur*tspan2);
y_eur=a_eur*sin(ni0_eur+w_eur*tspan2);
x_sc=W5(:,1);
y_sc=W5(:,2);

modE5=sqrt((x_eur-x_sc).^2+(y_eur-y_sc).^2);
sup5=min(modE5);
ind=find(modE5==sup5);

phi=reshape(W5(ind,6:30),5,5);
Pcov_minE5=phi*Pcov_est4*phi';

figure(2*iterations+8)
plot(tspan2,modE5,'k-');
xlabel('Time [s]')
ylabel('Distance to Europa [km]')
title('Approach to Europa after second observation')

grad=[(x_sc(ind)-x_eur(ind))/sup5 ; (y_sc(ind)-y_eur(ind))/sup5];
sigma_modE5=sqrt(grad'*Pcov_minE5(1:2,1:2)*grad);
Pcov_minE5(2,1)=Pcov_minE5(1,2);

[Vc,Dc]=eig(Pcov_minE5(1:2,1:2));
axis_a=sqrt(Dc(1,1));
axis_b=sqrt(Dc(2,2));
eig_v1=Vc(:,1);
eig_v2=Vc(:,2);

if eig_v1(2)>0
    theta=acos(eig_v1(1));
elseif eig_v1(2)<0
    theta=2*pi-acos(eig_v1(1));
end

XCA_eur_after=[x_sc(ind); y_sc(ind) ; x_eur(ind) ; y_eur(ind) ; axis_a ; axis_b ; theta];


%% parte 6

R_gan=2634.1;
R_eur=1560.8;


% XCA_gan_after
[tspan3,W6]=ode113('Model_and_transition_ch',tspan3,[state4' supp1],options);

x_gan=a_gan*cos(ni0_gan+w_gan*tspan3);
y_gan=a_gan*sin(ni0_gan+w_gan*tspan3);
x_sc=W6(:,1);
y_sc=W6(:,2);

modG=sqrt((x_gan-x_sc).^2+(y_gan-y_sc).^2);
sup6a=min(modG);
ind=find(modG==sup6a);

phi=reshape(W6(ind,6:30),5,5);
Pcov_minG6=phi*Pcov_est4*phi';

grad=[(x_sc(ind)-x_gan(ind))/sup6a ; (y_sc(ind)-y_gan(ind))/sup6a];
sigma_modG6=sqrt(grad'*Pcov_minG6(1:2,1:2)*grad);
Pcov_minG6(2,1)=Pcov_minG6(1,2);

[Vd,Dd]=eig(Pcov_minG6(1:2,1:2));
axis_a=sqrt(Dd(1,1));
axis_b=sqrt(Dd(2,2));
eig_v1=Vd(:,1);
eig_v2=Vd(:,2);

if eig_v1(2)>0
    theta=acos(eig_v1(1));
elseif eig_v1(2)<0
    theta=2*pi-acos(eig_v1(1));
end

XCA_gan_after=[x_sc(ind) ; y_sc(ind) ; x_gan(ind) ; y_gan(ind) ; axis_a ; axis_b ; theta];


% XCA_eur_before
[tspan2,W6]=ode113('Model_and_transition_ch',tspan2,[oldstate0' supp1],options);

x_eur=a_eur*cos(ni0_eur+w_eur*tspan2);
y_eur=a_eur*sin(ni0_eur+w_eur*tspan2);
x_sc=W6(:,1);
y_sc=W6(:,2);

modE=sqrt((x_eur-x_sc).^2+(y_eur-y_sc).^2);
sup6b=min(modE);
ind=find(modE==sup6b);

phi=reshape(W6(ind,6:30),5,5);
Pcov_minE6=phi*Pcov0*phi';

grad=[(x_sc(ind)-x_eur(ind))/sup6b ; (y_sc(ind)-y_eur(ind))/sup6b];
sigma_modE6=sqrt(grad'*Pcov_minE6(1:2,1:2)*grad);

[Ve,De]=eig(Pcov_minE6(1:2,1:2));
axis_a=sqrt(De(1,1));
axis_b=sqrt(De(2,2));
eig_v1=Ve(:,1);
eig_v2=Ve(:,2);

if eig_v1(2)>0
    theta=acos(eig_v1(1));
elseif eig_v1(2)<0
    theta=2*pi-acos(eig_v1(1));
end

XCA_eur_before=[x_sc(ind); y_sc(ind) ; x_eur(ind) ; y_eur(ind) ; axis_a ; axis_b ; theta];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% XCA_gan_before    calcolato nella parte 3
% XCA_eur_before    appena calcolato

% XCA_gan_after     appena calcolato
% XCA_eur_after     calcolato nella parte 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2*iterations+9)
prob1=draw_ellipse6(XCA_gan_before,R_gan,Pcov_minGa(1:2,1:2),2*iterations+9);
xlabel('[km]')
ylabel('[km]')
title('Ganymede CA before observation campaign')

figure(2*iterations+10)
prob2=draw_ellipse6(XCA_gan_after,R_gan,Pcov_minG6(1:2,1:2),2*iterations+10);
xlabel('[km]')
ylabel('[km]')
title('Ganymede CA after observation campaign')

figure(2*iterations+11)
prob3=draw_ellipse6(XCA_eur_before,R_eur,Pcov_minE6(1:2,1:2),2*iterations+11);
xlabel('[km]')
ylabel('[km]')
title('Europa CA before observation campaign')

figure(2*iterations+12)
prob4=draw_ellipse6(XCA_eur_after,R_eur,Pcov_minE5(1:2,1:2),2*iterations+12);
xlabel('[km]')
ylabel('[km]')
title('Europa CA after observation campaign')
