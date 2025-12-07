function[]=DisegnaTraj3D(ind_figura,rMatrixECI,color_fig)
% DisegnaTraj3D
%   Detailed explanation goes here

figure(ind_figura)
hold on

% disegno traiettoria
X=rMatrixECI(1,:);
Y=rMatrixECI(2,:);
Z=rMatrixECI(3,:);
plot3(X,Y,Z,color_fig);

% dettagli figura
title('Traiettoria in 3D')
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
daspect([1 1 1]);
ax=gca;
ax.XAxis.Exponent=0;
ax.YAxis.Exponent=0;
ax.ZAxis.Exponent=0;

hold off
view(3)

%save pdf
end

