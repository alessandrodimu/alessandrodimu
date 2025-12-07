function[r,LAT,LONG]=rect2polar(rVectECEF)
% rect2polar

[n,m]=size(rVectECEF);
if isvector(rVectECEF)==0 || n~=3 || m~=1
    error('rVect deve essere un vettore colonna 3x1');
end

r=norm(rVectECEF);
LAT=asin(rVectECEF(3)/r);
LONG=atan2(rVectECEF(2),rVectECEF(1));
end