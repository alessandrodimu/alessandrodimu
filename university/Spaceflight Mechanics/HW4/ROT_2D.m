function[tabella_posizioni_2D]=ROT_2D(tabella_posizioni_2D,angolo_rot)
% ROT_2D
%
% Effettua una rotazione vettoriale in 2D
ROT=[cos(angolo_rot) sin(angolo_rot);-sin(angolo_rot) cos(angolo_rot)];
tabella_posizioni_2D=ROT*tabella_posizioni_2D;
end