% HOMEWORK 1

%% WIPEOUT
clear
close all
clc



%% DATA
global RE RM mu beta 

in=input('Please select the mode: (1-2)\n 1) Continuous arc\n 2) 3 arcs\n\n');
if isscalar(in)==0 || isreal(in)==0 || in<1 || in>2 %#ok<COMPNOT>
    error('Expected input must be either 1 or 2')
end

% constants
RE=1.49598023*10^8;                 % in km
RM=2.27939200*10^8;                 % in km
c=299792.458;                       % in km/s
mu=1.327124402*10^11;               % in km^3/s^2
WE=1361;                            % in W/m^2
sigma=20*10^3;                      % in kg/km^2
beta=2*WE*RE^2/(c*sigma);           % in km^3/s^2



%% INITIALIZATION

% initial and final state (fixed)
r_0=RE;                             % in km
theta_0=0;                          % in rad
vr_0=0;                             % in km/s
vtheta_0=sqrt(mu/r_0);              % in km/s
K_0=0;                              % adimensional

r_f=RM;                             % in km
theta_f=[];                         % we don't care about it
vr_f=0;                             % in km/s
vtheta_f=sqrt(mu/r_f);              % in km/s
K_f=[];                             % we don't care about it

switch in
    case 1      % Continuous arc
        
        % first guess values for adjoints and time (variable)
        pr_0_g=-5.337648476397696e+04;
        ptheta_0=0;
        pvr_0_g=-1.953801489538171e+10;
        pvtheta_0_g=-8.827947513727990e+10;
        pk_0=1;
        t_f_g=700*(24*3600);
        
        
    case 2      % 3 arcs
        
        % first guess values for alpha and time (variable)
        alpha1_g=1.091537665072149;
        t1_g=2.369216999978725e+02*(24*3600);
        alpha2_g=0.536875320685295;
        t2_g=2.772055010000000e+02*(24*3600);
        alpha3_g=0.870867111400536;
        t3_g=2.610262008155346e+02*(24*3600);
        
        
end




%% OPTIMIZATION

switch in
    case 1      % Continuous arc
        
        x0_g=[ pr_0_g ; pvr_0_g ; pvtheta_0_g ; t_f_g ];
        
        options_opt=optimset('Display','off','PlotFcns','optimplotfval','TolFun',...
            1e-12,'TolX',1e-12,'MaxFunEvals',800);
        x0_opt=fminsearch('J_function_cont',x0_g,options_opt);
        
        pr_0=x0_opt(1);
        pvr_0=x0_opt(2);
        pvtheta_0=x0_opt(3);
        t_f=x0_opt(4);
        
        
    case 2      % 3 arcs
        
        x0_g=[ alpha1_g ; t1_g ; alpha2_g ; t2_g ; alpha3_g ; t3_g ];
        
        options_opt=optimset('Display','off','PlotFcns','optimplotfval','TolFun',...
            1e-12,'TolX',1e-12,'MaxFunEvals',800);
        x0_opt=fminsearch('J_function_3arc',x0_g,options_opt);
        
        alpha1=x0_opt(1);
        t1=x0_opt(2);
        alpha2=x0_opt(3);
        t2=x0_opt(4);
        alpha3=x0_opt(5);
        t3=x0_opt(6);
        
        
end



%% INTEGRATION

switch in
    case 1      % Continuous arc
        
        x0=[r_0; theta_0; vr_0 ; vtheta_0 ; K_0; pr_0 ; pvr_0 ; pvtheta_0];
        
        options_int=odeset('RelTol',1e-12,'AbsTol',1e-12);
        [t,W]=ode45('hamiltonian_eqs',[0 t_f],x0,options_int);

        r_table=W(:,1)';
        theta_table=W(:,2)';
        vr_table=W(:,3)';
        vtheta_table=W(:,4)';
        K_table=W(:,5);
        pr_table=W(:,6)';
        ptheta_table=zeros(1,length(t));
        pvr_table=W(:,7)';
        pvtheta_table=W(:,8)';
        pK_table=ones(1,length(t));

        alpha_table=rad2deg(atan((-3*pvr_table-sqrt(9*pvr_table.^2+8*pvtheta_table.^2))./(4*pvtheta_table)));

        t=t/(24*3600);              % in days

        H_table=pr_table.*vr_table+ptheta_table.*(vtheta_table./r_table)+...
            pvr_table.*((vtheta_table.^2)./r_table-mu./(r_table.^2)+beta./(r_table.^2).*(cos(alpha_table).^3))+...
            pvtheta_table.*((-vr_table.*vtheta_table)./r_table+beta./(r_table.^2).*(cos(alpha_table).^2).*sin(alpha_table))+...
            pK_table*(-1);

        % values at the orbit injection
        delta_r=abs(r_table(end)-r_f);
        delta_vr=abs(vr_table(end)-vr_f);
        delta_vtheta=abs(vtheta_table(end)-vtheta_f);
        
        % print important values
        fprintf('\nErrors at the orbit injection\n\n');
        fprintf('delta_r = %e km \n',delta_r);
        fprintf('delta_vr = %e km/s \n',delta_vr);
        fprintf('delta_vt = %e km/s \n\n',delta_vtheta);
        fprintf('t_f = %3.1f days\n',t(end));
        fprintf('H = %5.4f\n',H_table(end));


    case 2      % 3 arcs
        
        options_int=odeset('RelTol',1e-12,'AbsTol',1e-12);
        
        % first arc
        x0_1=[ r_0 ; theta_0 ; vr_0 ; vtheta_0 ; alpha1 ];
        [t1_table,W1]=ode45('state_eqs',[0 t1],x0_1,options_int);
        
        r1_table=W1(:,1)';
        theta1_table=W1(:,2)';
        vr1_table=W1(:,3)';
        vtheta1_table=W1(:,4)';
        alpha1_table=rad2deg(W1(:,5))';
        t1_table=t1_table/(24*3600)';
        
        % second arc
        x0_2=[ W1(end,1:4)' ; alpha2 ];
        [t2_table,W2]=ode45('state_eqs',[t1 (t1+t2)],x0_2,options_int);
        
        r2_table=W2(:,1)';
        theta2_table=W2(:,2)';
        vr2_table=W2(:,3)';
        vtheta2_table=W2(:,4)';
        alpha2_table=rad2deg(W2(:,5))';
        t2_table=t2_table/(24*3600)';
        
        % third arc
        x0_3=[ W2(end,1:4)' ; alpha3 ];
        [t3_table,W3]=ode45('state_eqs',[(t1+t2) (t1+t2+t3)],x0_3,options_int);
        
        r3_table=W3(:,1)';
        theta3_table=W3(:,2)';
        vr3_table=W3(:,3)';
        vtheta3_table=W3(:,4)';
        alpha3_table=rad2deg(W3(:,5))';
        t3_table=t3_table/(24*3600)';
        
        % values at the orbit injection
        delta_r=abs(r3_table(end)-r_f);
        delta_vr=abs(vr3_table(end)-vr_f);
        delta_vtheta=abs(vtheta3_table(end)-vtheta_f);

        
        % print important values
        fprintf('\nErrors at the orbit injection\n\n');
        fprintf('delta_r = %e km \n',delta_r);
        fprintf('delta_vr = %e km/s \n',delta_vr);
        fprintf('delta_vt = %e km/s \n\n',delta_vtheta);
        fprintf('t_f = %3.1f days\n',t3_table(end));
        
        
end



%% PLOTS

switch in
    case 1      % Continuous arc
        
        % figure(2)
        % plot(t,r_table,'b')
        % title('Radius','Interpreter','latex','Fontsize',12)
        % xlabel('$$t \ (days)$$','Interpreter','latex','Fontsize',12)
        % ylabel('$$r \ (km)$$','Interpreter','latex','Fontsize',12)
        % xlim([0 t(end)])
        % 
        % figure(3)
        % hold on
        % plot(t,vr_table,'b')
        % plot(t,vtheta_table,'r')
        % title('Velocity components','Interpreter','latex','Fontsize',12)
        % xlabel('$$t \ (days)$$','Interpreter','latex','Fontsize',12)
        % ylabel('$$v \ (km/s)$$','Interpreter','latex','Fontsize',12)
        % legend('$$v_r$$','$$v_{\theta}$$','Interpreter','latex','Fontsize',12)
        % xlim([0 t(end)])

        figure(4)
        plot(t,alpha_table,'b')
        title('Thrust pointing angle','Interpreter','latex','Fontsize',12)
        xlabel('$$t \ (days)$$','Interpreter','latex','Fontsize',12)
        ylabel('$$\alpha \ (deg)$$','Interpreter','latex','Fontsize',12)
        xlim([0 t(end)])

        % figure(5)
        % hold on
        % plot(t,pr_table,'r')
        % plot(t,ptheta_table,'g')
        % plot(t,pvr_table,'b')
        % plot(t,pvtheta_table,'k')
        % plot(t,pK_table,'m')
        % title('Lagrangian multipliers','Interpreter','latex','Fontsize',12)
        % xlabel('$$t \ (days)$$','Interpreter','latex','Fontsize',12)
        % ylabel('$$p$$','Interpreter','latex','Fontsize',12)
        % legend('$$p_r$$','$$p_{\theta}$$','$$p_{v_r}$$','$$p_{v_{\theta}}$$','$$p_K$$','Interpreter','latex','Fontsize',12)
        % xlim([0 t(end)])

        figure(6)
        [xE_table,yE_table]=pol2cart(linspace(0,2*pi,1000),RE*ones(1,1000));
        [xM_table,yM_table]=pol2cart(linspace(0,2*pi,1000),RM*ones(1,1000));
        [x_table,y_table]=pol2cart(theta_table,r_table);
        hold on
        plot(xE_table,yE_table,'k')
        plot(xM_table,yM_table,'k')
        plot(x_table,y_table,'b')
        title('Transfer trajectory','Interpreter','latex','Fontsize',12)
        xlabel('$$x \ (km)$$','Interpreter','latex','Fontsize',12)
        ylabel('$$y \ (km)$$','Interpreter','latex','Fontsize',12)
        legend('Initial orbit','Final orbit','Transfer trajectory','Interpreter','latex','Location','southeast')
        axis('equal')
        
        
    case 2      % 3 arcs
        
        figure(4)
        plot([t1_table' t2_table' t3_table'],[alpha1_table alpha2_table alpha3_table],'b')
        title('Thrust pointing angle','Interpreter','latex','Fontsize',12)
        xlabel('$$t \ (days)$$','Interpreter','latex','Fontsize',12)
        ylabel('$$\alpha \ (deg)$$','Interpreter','latex','Fontsize',12)
        xlim([0 t3_table(end)])
        
        figure(6)
        [xE_table,yE_table]=pol2cart(linspace(0,2*pi,1000),RE*ones(1,1000));
        [xM_table,yM_table]=pol2cart(linspace(0,2*pi,1000),RM*ones(1,1000));
        [x1_table,y1_table]=pol2cart(theta1_table,r1_table);
        [x2_table,y2_table]=pol2cart(theta2_table,r2_table);
        [x3_table,y3_table]=pol2cart(theta3_table,r3_table);
        hold on
        plot(xE_table,yE_table,'k')
        plot(xM_table,yM_table,'k')
        plot(x1_table,y1_table,'b')
        plot(x2_table,y2_table,'r')
        plot(x3_table,y3_table,'g')
        title('Transfer trajectory','Interpreter','latex','Fontsize',12)
        xlabel('$$x \ (km)$$','Interpreter','latex','Fontsize',12)
        ylabel('$$y \ (km)$$','Interpreter','latex','Fontsize',12)
        legend('Initial orbit','Final orbit','Transfer arc 1','Transfer arc 2','Transfer arc 3','Interpreter','latex','Location','southeast')
        axis('equal')
        
        
end
