function[]=DisegnaTerra(ind_figura)
% DisegnaTerra

figure(ind_figura)
hold on

k=5;
n=2^k-1;
theta=pi*(-n:2:n)/n;
phi=(pi/2)*(-n:2:n)'/n;
X=6371*cos(phi)*sin(theta);
Y=6371*cos(phi)*cos(theta);
Z=6371*sin(phi)*ones(size(theta));
surf(X,-Y,-Z);
map=imread('mappa.jpg');
h=findobj('Type','surface');
set(h,'CData',map,'FaceColor','texturemap','edgecolor','none');
daspect([1 1 1]);
grid on

hold off
view(3)
end

