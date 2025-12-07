function dw=Model_and_transition_ch(t,w)
% Purpose: Integrate a simple orbital problem. Forces: Gravity from 3
% different bodies
% Input: w = [X Y U V GM_gan]
% Output dw = [dX dY dU dV dGM_gan]


global GM_jup GM_eur ni0_eur ni0_gan a_eur a_gan w_eur w_gan

x_eur=a_eur*cos(ni0_eur+w_eur*t);
y_eur=a_eur*sin(ni0_eur+w_eur*t);
x_gan=a_gan*cos(ni0_gan+w_gan*t);
y_gan=a_gan*sin(ni0_gan+w_gan*t);

dw=zeros(5+25,1);

% State integration
X=w(1:5);
x=w(1);
y=w(2);
dx=w(3);
dy=w(4);
GM_gan=w(5);

r=sqrt(x^2+y^2);
E=sqrt((x_eur-x)^2+(y_eur-y)^2);
G=sqrt((x_gan-x)^2+(y_gan-y)^2);

dw(1)=dx;
dw(2)=dy;
dw(3)=(-GM_jup/r^3-GM_eur/E^3-GM_gan/G^3)*x+(GM_eur*x_eur)/E^3+(GM_gan*x_gan)/G^3;
dw(4)=(-GM_jup/r^3-GM_eur/E^3-GM_gan/G^3)*y+(GM_eur*y_eur)/E^3+(GM_gan*y_gan)/G^3;
dw(5)=0;


% STM integration
phi=w(6:30);
PHI=reshape(phi,5,5);
A=A_matrix_ch(X,t);
dphi=A*PHI;
dphi=reshape(dphi,1,25);

dw(6:30)=dphi;

end
