function[rVectECEF]=ECI2ECEF(MJD,rVectECI)
% ECI2ECEF

if isscalar(MJD)==0 || isreal(MJD)==0
    error('MJD deve essere uno scalare reale');
end
[n,m]=size(rVectECI);
if isvector(rVectECI)==0 || n~=3 || m~=1
    error('rVectECI deve essere un vettore colonna 3x1');
end

T_solare=86400;
MJD_0=49718*T_solare;            % MJD_0=1/1/1995 @ 00:00:00
T_siderale=86164;
w_rot_terra=2*pi/T_siderale;
RAAN_greenwich_0=100*pi/180;        % e' un alpha in realta'
RAAN_greenwich=RAAN_greenwich_0+w_rot_terra*(MJD-MJD_0);

ECI2ECEF=[cos(RAAN_greenwich) sin(RAAN_greenwich) 0;...
    -sin(RAAN_greenwich) cos(RAAN_greenwich) 0;0 0 1];
rVectECEF=ECI2ECEF*rVectECI;
end

