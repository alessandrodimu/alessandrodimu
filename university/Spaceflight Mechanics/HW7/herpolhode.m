function herpolhode(J,W,Hc,Er,fig_ind)

Tol0 = 1e-12; 
Tol1 = 1e-09;
options = odeset('RelTol',Tol0,'AbsTol',Tol1);

tempo=2;

figure(fig_ind)
hold on


% disegno il piano su cui giace l'erpoloide
d=2*Er/norm(Hc);
p1=[-6 6 6 -6];
p2=[6 6 -6 -6];
p3=[d d d d];

patch('XData',p1,'YData',p2,'ZData',p3,'FaceColor',[135 206 250]/255,'FaceAlpha',.5)
title('Erpoloide e poloide')
daspect([1 1 1])
axis off
grid off
view(234,14)


% trovo i valori iniziali di integrazione
w1=W(1);
w2=W(2);
w3=W(3);


% calcolo assetto iniziale
n3=Hc/norm(Hc);
n2=cross(Hc,W)/norm(cross(Hc,W));
n1=cross(n2,n3);
BN=[n1 n2 n3];
[c11,c21,c31,c12,c22,c23]=decod_0(BN);


% disegno l'erpoloide sul piano
[~,Wa]=ode45(@int_cos0,[0 tempo],[w1 w2 w3 c11 c21 c31 c12 c22 c23],options);
nel=size(Wa);
Werp_inerz=zeros(nel(1),3);
for ind_1=1:1:nel(1)
    n1=[Wa(ind_1,4); Wa(ind_1,5); Wa(ind_1,6)];
    n2=[Wa(ind_1,7); Wa(ind_1,8); Wa(ind_1,9)];
    n3=cross(n1,n2);
    BN=[n1 n2 n3];
    Werp=[Wa(ind_1,1); Wa(ind_1,2); Wa(ind_1,3)];
    Werp_inerz(ind_1,:)=BN'*Werp;
end
plot3(Werp_inerz(:,1),Werp_inerz(:,2),Werp_inerz(:,3),'b-');


% disegno l'ellissoide di costante energia ruotato
[X,Y,Z]=sphere(200);
X=X*sqrt(2*Er/J(1));
Y=Y*sqrt(2*Er/J(2));
Z=Z*sqrt(2*Er/J(3));

BN_fissato=BN;
for ind_2=1:201
    for ind_3=1:201
        v=[X(ind_2,ind_3) ; Y(ind_2,ind_3)  ;Z(ind_2,ind_3)];
        v_inerz=BN_fissato'*v;
        X(ind_2,ind_3)=v_inerz(1);
        Y(ind_2,ind_3)=v_inerz(2);
        Z(ind_2,ind_3)=v_inerz(3);
    end
end
surf(X,Y,Z,'FaceColor',[1 1 0.47],'EdgeColor','none')
camlight 
lighting gouraud
alpha (0.5)


% disegno le direzioni principali del body
b=sqrt(2*Er./J)+0.5;
PP=0.4*[1 0 0;0 1 0;0 0 1];
PP_rot=BN_fissato'*PP;
P1=[b(1) ; 0 ; 0];
P1=BN_fissato'*P1;
plot3([-P1(1) P1(1)],[-P1(2) P1(2)],[-P1(3) P1(3)],'k','linewidth',2);
P1=P1+PP_rot(:,1);
text (P1(1),P1(2),P1(3),'$\hat{b_{1}}$','Interpreter','latex','fontsize',12);
text (-P1(1),-P1(2),-P1(3),'$\hat{b_{1}}$','Interpreter','latex','fontsize',12);

P2=[0 ; b(2) ; 0];
P2=BN_fissato'*P2;
plot3([-P2(1) P2(1)],[-P2(2) P2(2)],[-P2(3) P2(3)],'k','linewidth',2);
P2=P2+PP_rot(:,2);
text (P2(1),P2(2),P2(3),'$\hat{b_{2}}$','Interpreter','latex','fontsize',12);
text (-P2(1),-P2(2),-P2(3),'$\hat{b_{2}}$','Interpreter','latex','fontsize',12);

P3=[0 ; 0 ; b(3)];
P3=BN_fissato'*P3;
plot3([-P3(1) P3(1)],[-P3(2) P3(2)],[-P3(3) P3(3)],'k','linewidth',2);
P3=P3+PP_rot(:,3);
text (P3(1),P3(2),P3(3),'$\hat{b_{3}}$','Interpreter','latex','fontsize',12);
text (-P3(1),-P3(2),-P3(3),'$\hat{b_{3}}$','Interpreter','latex','fontsize',12);

quiver3(0,0,0,2*Werp_inerz(end,1),2*Werp_inerz(end,2),2*Werp_inerz(end,3),...
    'Color',[255 105 180]/255,'LineWidth',2);
text(2*Werp_inerz(end,1),2*Werp_inerz(end,2),2*Werp_inerz(end,3),...
    '$\vec{\omega}$','Interpreter','latex','fontsize',15);


% disegno la poloide ruotata
[~,W1]=ode45(@int_omega,[0 tempo],[w1 w2 w3],options);
ne2=size(W1);
W1_inerz=zeros(ne2(1),3);
for ind_4=1:ne2(1)
    W1_body=[W1(ind_4,1) ; W1(ind_4,2) ; W1(ind_4,3)];
    W1_inerz(ind_4,:)=BN_fissato'*W1_body;
end
plot3(W1_inerz(:,1),W1_inerz(:,2),W1_inerz(:,3),'r-');
hold off


end

