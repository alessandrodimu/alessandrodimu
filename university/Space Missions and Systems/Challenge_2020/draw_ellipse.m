function[]=draw_ellipse(a,b,theta,i)
% draw_ellipse will draw the rotated ellipse starting from its parameters
%   INPUT: a       semi-major axis
%          b       semi-minor axis
%          theta   rotation angle around k3

phi=linspace(0,2*pi,100);

x=a*cos(phi);
y=b*sin(phi);
X=[x ; y];

rot=[cos(theta) sin(theta); -sin(theta) cos(theta)];
X=rot*X;

figure(i)
plot(X(1,:),X(2,:),'k-')

end

