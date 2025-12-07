function [a1,d1,a2,d2,a3,d3]=decod_1(BN) 

a1=atan2(BN(1,2),BN(1,1));
b1=sqrt(BN(1,1)^2+BN(1,2)^2);
if b1==0
    if BN(1,3)>0
        d1=pi/2;
    elseif BN(1,3)<0
        d1=-pi/2;
    end
else
    d1=atan(BN(1,3)/b1);
end
a2=atan2(BN(2,2),BN(2,1));
b2=sqrt(BN(2,1)^2+BN(2,2)^2);
if b2==0
    if BN(2,3)>0
        d2=pi/2;
    elseif BN(2,3)<0
        d2=-pi/2;
    end
else
    d2=atan(BN(2,3)/b2);
end
a3=atan2(BN(3,2),BN(3,1));
b3=sqrt(BN(3,1)^2+BN(3,2)^2);
if b3==0
    if BN(3,3)>0
        d3=pi/2;
    elseif BN(3,3)<0
        d3=-pi/2;
    end
else
    d3=atan(BN(3,3)/b3);
end

