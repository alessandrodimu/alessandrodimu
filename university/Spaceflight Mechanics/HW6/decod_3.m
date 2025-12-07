function [phi,the,psi]=decod_3(BN) 

the=asin(-BN(1,3));
acthe=abs(cos(the));
if(acthe > .000001)
    phi=atan2(BN(2,3),BN(3,3));
    psi=atan2(BN(1,2),BN(1,1));
else
    phi=0;
    if(BN(1,3) > 0) %theta=-pi/2
        psi=atan2(-BN(3,2),-BN(3,1));
    else
        psi=atan2(BN(3,2),BN(3,1));
    end
end

