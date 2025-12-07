function [Ek] = InvertedKeplero(M,e)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
E0=0;
err1=1;
err2=1;
f=@(x) (M-x+e*sin(x));
df=@(x) (-1+e*cos(x));
while err1>=10^(-4) && err2>=10^(-5)
    Ek=E0-f(E0)/df(E0);
    err1=abs(E0-Ek);
    err2=abs(f(Ek));
    E0=Ek;
end