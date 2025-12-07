function [E,phi]=decod_2(BN) 

cphi=(BN(1,1)+BN(2,2)+BN(3,3)-1)/2;
if(cphi < .999999)
    if(cphi > -.999999)
        phi=acos(cphi);
        sphi2=2*sin(phi);
        E(1)=(BN(2,3)-BN(3,2))/sphi2;
        E(2)=(BN(3,1)-BN(1,3))/sphi2;
        E(3)=(BN(1,2)-BN(2,1))/sphi2;
    else
        E(1)=sqrt((BN(1,1)+1)/2);
        imx=1;
        ia=2;
        ib=3;
        E(2)=sqrt((BN(2,2)+1)/2);
        if(E(2) > E(1))
            imx=2;
            ia=1;
        end
        E(3)=sqrt((BN(3,3)+1)/2);
        if(E(3) > E(imx))
            ib=imx;
            imx=3;
        end
        if(BN(imx,ia) < 0)
            E(ia)=-E(ia);
        end
        if(BN(imx,ib) < 0)
            E(ib)=-E(ib);
        end
    end
else
    phi=0;
    E=[1 0 0];
end

