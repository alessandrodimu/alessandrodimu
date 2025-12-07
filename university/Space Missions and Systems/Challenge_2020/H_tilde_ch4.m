function Ht=H_tilde_ch4(X,t)

global ni0_eur a_eur w_eur

x_eur=a_eur*cos(ni0_eur+w_eur*t);
y_eur=a_eur*sin(ni0_eur+w_eur*t);

x=X(1);
y=X(2);

E=sqrt((x_eur-x)^2+(y_eur-y)^2);


Ht=zeros(2,5);

% First row
Ht(1,1)=(x_eur-x)^2/E^3 - 1/E;
Ht(1,2)=(x_eur-x)*(y_eur-y)/E^3;

% Second row
Ht(2,1)=(x_eur-x)*(y_eur-y)/E^3;
Ht(2,2)=(y_eur-y)^2/E^3 - 1/E;

end
