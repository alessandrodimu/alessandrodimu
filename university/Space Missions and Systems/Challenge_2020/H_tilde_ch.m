function Ht=H_tilde_ch(X,t)

global ni0_gan a_gan w_gan

x_gan=a_gan*cos(ni0_gan+w_gan*t);
y_gan=a_gan*sin(ni0_gan+w_gan*t);

x=X(1);
y=X(2);

G=sqrt((x_gan-x)^2+(y_gan-y)^2);


Ht=zeros(2,5);

% First row
Ht(1,1)=(x_gan-x)^2/G^3 - 1/G;
Ht(1,2)=(x_gan-x)*(y_gan-y)/G^3;

% Second row
Ht(2,1)=(x_gan-x)*(y_gan-y)/G^3;
Ht(2,2)=(y_gan-y)^2/G^3 - 1/G;

end
