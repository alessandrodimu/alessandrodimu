function [A]=A_matrix_ch(X,t)
% The input X is: X=[x,y,u,v,GM_gan]


global GM_jup GM_eur ni0_eur ni0_gan a_eur a_gan w_eur w_gan

x_eur=a_eur*cos(ni0_eur+w_eur*t);
y_eur=a_eur*sin(ni0_eur+w_eur*t);
x_gan=a_gan*cos(ni0_gan+w_gan*t);
y_gan=a_gan*sin(ni0_gan+w_gan*t);


A=zeros(5,5);

% For simplicity we rename variables
x=X(1);
y=X(2);
u=X(3);
v=X(4);
GM_gan=X(5);

r=sqrt(x^2+y^2);
E=sqrt((x_eur-x)^2+(y_eur-y)^2);
G=sqrt((x_gan-x)^2+(y_gan-y)^2);

% First two row
A(1,3)=1;
A(2,4)=1;

% Third row
A(3,1)=(-GM_jup/r^3+(3*GM_jup*x^2)/r^5) + (-GM_eur/E^3+(3*GM_eur*(x-x_eur)^2)/E^5) + (-GM_gan/G^3+(3*GM_gan*(x-x_gan)^2)/G^5);
A(3,2)=((3*GM_jup*x*y)/r^5) + ((3*GM_eur*(x-x_eur)*(y-y_eur))/E^5) + ((3*GM_gan*(x-x_gan)*(y-y_gan))/G^5);
A(3,5)=(x_gan-x)/G^3;

% Fourth row
A(4,1)=((3*GM_jup*x*y)/r^5) + ((3*GM_eur*(x-x_eur)*(y-y_eur))/E^5) + ((3*GM_gan*(x-x_gan)*(y-y_gan))/G^5);
A(4,2)=(-GM_jup/r^3+(3*GM_jup*y^2)/r^5) + (-GM_eur/E^3+(3*GM_eur*(y-y_eur)^2)/E^5) + (-GM_gan/G^3+(3*GM_gan*(y-y_gan)^2)/G^5);
A(4,5)=(y_gan-y)/G^3;

end