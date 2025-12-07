function[rMatrixECI,vMatrixECI]=ShellPropagator(r0,v0,tspan)
% ShellPropagator
%   Detailed explanation goes here

% controllo su tspan
[n,~]=size(tspan);
if isvector(tspan)==0 || n~=1
    error('tspan deve essere un vettore riga 1xn')
end

rMatrixECI=zeros(3,length(tspan));
vMatrixECI=zeros(3,length(tspan));
for i=1:length(tspan)
    Dt=tspan(i)-tspan(1);
    [supp1,supp2]=CorePropagator(r0,v0,Dt);
    rMatrixECI(:,i)=supp1;
    vMatrixECI(:,i)=supp2;
end

