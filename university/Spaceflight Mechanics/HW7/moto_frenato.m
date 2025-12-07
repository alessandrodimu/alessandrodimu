function moto_frenato(J_0,W_0,Hc_0,Er_0,fig_ind)
global J Er eps 

Tol0 = 1e-12; 
Tol1 = 1e-09;
options = odeset('RelTol',Tol0,'AbsTol',Tol1);

tempo=100;
eps=5/tempo;
J=J_0;
Er=Er_0;

% trovo i valori iniziali di integrazione
w10=W_0(1);
w20=W_0(2);
w30=W_0(3);

[t,W_forced]=ode45(@int_omega1,[0 tempo],[w10 w20 w30],options);
nel=size(t);
w1=zeros(nel(1),1);
w2=zeros(nel(1),1);
w3=zeros(nel(1),1);
mod_w=zeros(nel(1),1);
h1=zeros(nel(1),1);
h2=zeros(nel(1),1);
h3=zeros(nel(1),1);
Hc=zeros(nel(1),1);
Er=zeros(nel(1),1);
for i=1:nel(1)
    w1(i,1)=W_forced(i,1);
    w2(i,1)=W_forced(i,2);
    w3(i,1)=W_forced(i,3);
    mod_w(i,1)=sqrt(w1(i,1)^2+w2(i,1)^2+w3(i,1)^2);
    h1(i,1)=J(1)*w1(i,1);
    h2(i,1)=J(2)*w2(i,1);
    h3(i,1)=J(3)*w3(i,1);
    Hc(i,1)=sqrt(h1(i,1)^2+h2(i,1)^2+h3(i,1)^2);
    Er(i,1)=0.5*(J(1)*w1(i,1)^2+J(2)*w2(i,1)^2+J(3)*w3(i,1)^2);
end

[X,Y,Z]=sphere(200);
X=X/J(1);
Y=Y/J(2);
Z=Z/J(3);

figure(fig_ind(1));
hold on
surf(X,Y,Z,'FaceColor',[1 1 0.47],'EdgeColor','none');
title('Andamento del vettore omega per dissipazione')
camlight
lighting gouraud
alpha (0.5)
view (164,24)
daspect([1 1 1])
axis off
grid off

b=1./J;
b=b+.5;
p11=[-b(1) b(1)];
p12=[0 0];
p13=[0 0];
plot3(p11,p12,p13,'k','linewidth',2);
text (b(1)+.5,0,0,'$\hat{b_{1}}$','Interpreter','latex','fontsize',12);
text (-b(1)-.5,0,0,'$\hat{b_{1}}$','Interpreter','latex','fontsize',12);
p21=[0 0];
p22=[-b(2) b(2)];
p23=[0 0];
plot3(p21,p22,p23,'k','linewidth',2);
text (0,b(2)+.5,0,'$\hat{b_{2}}$','Interpreter','latex','fontsize',12);
text (0,-b(2)-.5,0,'$\hat{b_{2}}$','Interpreter','latex','fontsize',12);
p31=[0 0];
p32=[0 0];
p33=[-b(3) b(3)];
plot3(p31,p32,p33,'k','linewidth',2);
text (0,0,b(3)+.5,'$\hat{b_{3}}$','Interpreter','latex','fontsize',12);
text (0,0,-b(3)-.5,'$\hat{b_{3}}$','Interpreter','latex','fontsize',12);

h=animatedline(w1(1),w2(1),w3(1),'Color','r');
for ind=1:10:nel(1)
    addpoints(h,w1(ind,1),w2(ind,1),w3(ind,1));
    drawnow;
end
hold off

figure(fig_ind(2))
plot(t,w1(:,1),'b-',t,w2(:,1),'g-',t,w3(:,1),'m-',t,mod_w(:,1),'r-')
title('Andamento di w (modulo e componenti)');
legend('$\omega_{1}$','$\omega_{2}$','$\omega_{3}$','$\omega$','Interpreter','latex');
xlabel('time (s)')
ylabel(texlabel('vel angolare (rad/s)'))
xlim([0 tempo])
grid on

figure(fig_ind(3))
plot(t,h1(:,1),'b-',t,h2(:,1),'g-',t,h3(:,1),'m-',t,Hc(:,1),'r-')
title('Andamento di H (modulo e componenti)');
legend(texlabel('H_1'),texlabel('H_2'),texlabel('H_3'),texlabel('H_C'));
xlabel('time (s)')
ylabel(texlabel('momento angolare'))
axis([0 tempo -1 1.5])
grid on

figure(fig_ind(4))
plot(t,Er(:,1),'b-')
title('Andamento di Er');
xlabel('time (s)')
ylabel(texlabel('energia rotazionale'))
xlim([0 tempo])
grid on

