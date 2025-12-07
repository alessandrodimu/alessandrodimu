function ellissoide (J,W,Hc,Er,fig_ind) 

%vettore Hc
a=sqrt(2*J*Er);

b=zeros(1,3);
for i=1:3
    if(a(i) > 1)
        b(i)=a(i);
    else
        b(i)=1;
    end
end
b=b+.2;
bmax=max(b);

figure(fig_ind(1));
title('Sfera di costante momento ed ellissoide di costante energia');
axis off;
grid off;
hold on;
p11=[-b(1) b(1)];
p12=[0 0];
p13=[0 0];
plot3(p11,p12,p13,'k','linewidth',2);
text (b(1)+.2,0,0,'$\hat{b_{1}}$','Interpreter','latex','fontsize',12);
text (-b(1)-.2,0,0,'$\hat{b_{1}}$','Interpreter','latex','fontsize',12);
p21=[0 0];
p22=[-b(2) b(2)];
p23=[0 0];
plot3(p21,p22,p23,'k','linewidth',2);
text (0,b(2)+.2,0,'$\hat{b_{2}}$','Interpreter','latex','fontsize',12);
text (0,-b(2)-.2,0,'$\hat{b_{2}}$','Interpreter','latex','fontsize',12);
p31=[0 0];
p32=[0 0];
p33=[-b(3) b(3)];
plot3(p31,p32,p33,'k','linewidth',2);
text (0,0,b(3)+.2,'$\hat{b_{3}}$','Interpreter','latex','fontsize',12);
text (0,0,-b(3)-.2,'$\hat{b_{3}}$','Interpreter','latex','fontsize',12);
h1=[0 bmax*Hc(1)];
h2=[0 bmax*Hc(2)];
h3=[0 bmax*Hc(3)];
plot3(h1,h2,h3,'k','linewidth',3);
text (1.3*h1(2),1.3*h2(2),1.3*h3(2),'$\vec{H_{c}}$','Interpreter','latex','fontsize',15);
bmax=bmax/sqrt(W(1)^2+W(2)^2+W(3)^2);
h1=[0 bmax*W(1)];
h2=[0 bmax*W(2)];
h3=[0 bmax*W(3)];
plot3(h1,h2,h3,'k','linewidth',3);
text (1.2*h1(2),1.2*h2(2),1.2*h3(2),'$\vec{\omega}$','Interpreter','latex','fontsize',15);

[X,Y,Z] = sphere(500);
surf (X,Y,Z,'FaceColor',[1 1 0.47],'EdgeColor','none');
camlight 
lighting gouraud
alpha (0.5)

[X,Y,Z] = sphere(500);
X=X*a(1);
Y=Y*a(2);
Z=Z*a(3);
surf (X,Y,Z,'FaceColor','blue','EdgeColor','none');
camlight 
lighting gouraud
%lighting phong
alpha (0.5)
axis equal;
view (164,24);

%vettore omega
a=sqrt(2*Er./J);

for i=1:3
    if(a(i) > 1/J(i))
        b(i)=a(i);
    else
        b(i)=1/J(i);
    end
end
b=b+.5;
bmax=max(b);
    
figure(fig_ind(2));
title('Ellissoidi di costante momento e di costante energia');
axis off;
grid off;
hold on;
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
h1=[0 bmax*Hc(1)];
h2=[0 bmax*Hc(2)];
h3=[0 bmax*Hc(3)];
plot3(h1,h2,h3,'k','linewidth',3);
text (1.3*h1(2),1.3*h2(2),1.3*h3(2),'$\vec{H_{c}}$','Interpreter','latex','fontsize',15);
bmax=bmax/sqrt(W(1)^2+W(2)^2+W(3)^2);
h1=[0 bmax*W(1)];
h2=[0 bmax*W(2)];
h3=[0 bmax*W(3)];
plot3(h1,h2,h3,'k','linewidth',3);
text (1.2*h1(2),1.2*h2(2),1.2*h3(2),'$\vec{\omega}$','Interpreter','latex','fontsize',15);

[X,Y,Z] = sphere(500);
X=X/J(1);
Y=Y/J(2);
Z=Z/J(3);
surf (X,Y,Z,'FaceColor',[1 1 0.47],'EdgeColor','none');
camlight 
lighting gouraud
alpha (0.5)

[X,Y,Z] = sphere(500);
X=X*a(1);
Y=Y*a(2);
Z=Z*a(3);
surf (X,Y,Z,'FaceColor','blue','EdgeColor','none');
camlight 
lighting gouraud
%lighting phong
alpha (0.5)

axis equal;
view (164,24);






