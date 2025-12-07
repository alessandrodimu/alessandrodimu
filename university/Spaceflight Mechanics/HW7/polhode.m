function polhode(J,W,Hc,Er,fig_ind)
Tol0 = 1e-12; 
Tol1 = 1e-09;
options = odeset('RelTol',Tol0,'AbsTol',Tol1);

% bisogna trovare un tempo buono per la rappresentazione
tempo=50*pi;

% disegno l'ellissoide di costante momento
[X,Y,Z]=sphere(200);
X=X/J(1);
Y=Y/J(2);
Z=Z/J(3);

figure(fig_ind)
hold on
surf(X,Y,Z,'FaceColor',[1 1 0.47],'EdgeColor','none')
view (164,24)
daspect([1 1 1])
axis off
grid off
title('Poloide e separatrici')

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

supporto=size(W);
for i=1:supporto(2)
    % trovo i valori iniziali di integrazione
    w1=W(1,i);
    w2=W(2,i);
    w3=W(3,i);
    
    % integro w1 w2 w3 per ottenere la poloide
    [~,W1]=ode45(@int_omega,[0 tempo],[w1 w2 w3],options);
    [~,W2]=ode45(@int_omega,[0 tempo],[-w1 w2 -w3],options);
    
    plot3(W1(:,1),W1(:,2),W1(:,3),'r-')
    plot3(W2(:,1),W2(:,2),W2(:,3),'r-')
    
    W1=[];
    W2=[];
end
hold off

end

