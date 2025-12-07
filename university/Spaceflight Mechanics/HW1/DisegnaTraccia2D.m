function[]=DisegnaTraccia2D(ind_fig,LatVect,LongVect,color_fig)
% DisegnaTraccia2D

[n1,m1]=size(LatVect);
if isvector(LatVect)==0 || n1~=1
    error('LatVect deve essere un vettore riga 1xn');
end
[n2,m2]=size(LongVect);
if isvector(LongVect)==0 || n2~=1
    error('LongVect deve essere un vettore riga 1xn');
end
if m1~=m2
    error('LatVect e LongVect devono avere lo stesso numero di colonne');
end

figure(ind_fig)
hold on

plot(LongVect*180/pi,LatVect*180/pi,color_fig,'MarkerSize',0.5);

title('Traccia a terra')
xlabel('Longitudine')
ylabel('Latitudine')
grid on

hold off

% save pdf
end