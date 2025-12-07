function[P]=draw_ellipse6(XCA,Rmoon,Pcov,i)
% draw_ellipse will draw the rotated ellipse starting from its parameters
%   INPUT: a       semi-major axis
%          b       semi-minor axis
%          theta   rotation angle around k3

xsc=XCA(1);
ysc=XCA(2);
xmoon=XCA(3);
ymoon=XCA(4);
a=XCA(5);
b=XCA(6);
theta=XCA(7);

phi=linspace(0,2*pi,100);

% disegno ellisse
x_ell=a*cos(phi);
y_ell=b*sin(phi);
X_ell1=[x_ell ; y_ell];
X_ell2=2*X_ell1;
X_ell3=3*X_ell1;

rot=[cos(theta) sin(theta); -sin(theta) cos(theta)];
X_ell1=rot*X_ell1;
X_ell1=X_ell1 + [xsc; ysc]*ones(1,length(phi));
X_ell2=rot*X_ell2;
X_ell2=X_ell2 + [xsc; ysc]*ones(1,length(phi));
X_ell3=rot*X_ell3;
X_ell3=X_ell3 + [xsc; ysc]*ones(1,length(phi));


% disegno pianeta
x_circ=Rmoon*cos(phi);
y_circ=Rmoon*sin(phi);
X_circ=[x_circ ; y_circ] + [xmoon ; ymoon]*ones(1,length(phi));

% disegno quadrato
dirSC=[xsc-xmoon ; ysc-ymoon];
dirSC=dirSC/norm(dirSC);

halfside=(2*Rmoon+400)/2;
halfdiag=sqrt(2)*halfside;

%X_lim=halfdiag*[cos(phi) ; sin(phi)] + [xmoon ; ymoon]*ones(1,length(phi));

%punto1_sq=halfdiag*dirSC;
%X_sq=zeros(2,5);
%X_sq(:,1)=punto1_sq;
%for k=1:4
%    rot2=[cos(k*pi/2) sin(k*pi/2); -sin(k*pi/2) cos(k*pi/2)];
%    X_sq(:,k+1)=rot2*punto1_sq;
%end
%X_sq=X_sq + [xmoon ; ymoon]*ones(1,5);

x_sq2=[halfside halfside -halfside -halfside halfside];
y_sq2=[halfside -halfside -halfside halfside halfside];
X_sq2=[x_sq2 ; y_sq2] + [xmoon ; ymoon]*ones(1,5);

figure(i)
hold on
quiver(xmoon,ymoon,3*Rmoon*dirSC(1),3*Rmoon*dirSC(2),'m-')
plot(xsc,ysc,'r*')
plot(xmoon,ymoon,'r*')

plot(X_ell1(1,:),X_ell1(2,:),'k-')
plot(X_ell2(1,:),X_ell2(2,:),'k-')
plot(X_ell3(1,:),X_ell3(2,:),'k-')

plot(X_circ(1,:),X_circ(2,:),'k-')
%plot(X_sq(1,:),X_sq(2,:),'k-')
plot(X_sq2(1,:),X_sq2(2,:),'r-')
%plot(X_lim(1,:),X_lim(2,:),'k--')
hold off
axis equal

mean=[xsc ysc];
lower_limit=[X_sq2(1,3) X_sq2(2,3)];
upper_limit=[X_sq2(1,1) X_sq2(2,1)];
P=mvncdf(lower_limit,upper_limit,mean,Pcov);


end
