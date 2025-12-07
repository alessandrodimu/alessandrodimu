function[LatVect,LongVect]=ShellTrace(rMatrixECI,tspan_MJD)
% ShellTrace

[n1,m1]=size(rMatrixECI);
if ismatrix(rMatrixECI)==0 || n1~=3
    error('rMatrixECI deve essere una matrice 3xn');
end
[n2,m2]=size(tspan_MJD);
if isvector(tspan_MJD)==0 || n2~=1
    error('tspan_MJD deve essere un vettore riga 1xn');
end
if m1~=m2
    error('rMatrixECI e tspan_MJD devono avere numero di colonne uguale');
end

LatVect=zeros(1,m2);
LongVect=zeros(1,m2);
for i=1:m2
    D_MJD=tspan_MJD(i);  %-tspan_MJD(1);
    rVectECI=rMatrixECI(:,i);
    [Lat,Long]=CoreTrace(D_MJD,rVectECI);
    LatVect(i)=Lat;
    LongVect(i)=Long;
end