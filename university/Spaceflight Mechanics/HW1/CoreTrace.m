function[Lat,Long]=CoreTrace(MJD,rVectECI)
% CoreTrace

rVectECEF=ECI2ECEF(MJD,rVectECI);
[~,Lat,Long]=rect2polar(rVectECEF);
end

