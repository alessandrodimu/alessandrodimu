function [lat,lon,u,vect_e]=los_algorithm(hh,vect_s,vect_g)

% wgs84 ellipsoid
wgs84=wgs84Ellipsoid('m');
RE=wgs84.SemimajorAxis; % Earth equatorial radius WGS84
RP=wgs84.SemiminorAxis; % Earth polar radius WGS84


% algorithm that searches for the intersection point with ellipsoid
A=(RP+hh)^2*(vect_g(1)^2+vect_g(2)^2)+(RE+hh)^2*vect_g(3)^2;
B=2*((RP+hh)^2*(vect_s(1)*vect_g(1)+vect_s(2)*vect_g(2))+(RE+hh)^2*(vect_s(3)*vect_g(3)));
C=(RP+hh)^2*(vect_s(1)^2+vect_s(2)^2)+(RE+hh)^2*(vect_s(3)^2-(RP+hh)^2);

% distance (u) from satellite and intersect point vector (e)
u=(-B-sqrt(B^2-4*A*C))/(2*A);
vect_e=vect_s+u*vect_g;

if isreal(vect_e)==1
    % geodetic coordinates on said ellipsoid
    lon=atan2d(vect_e(2),vect_e(1));
    geoc_lat=atan2d(vect_e(3),sqrt(vect_e(1)^2+vect_e(2)^2));
    geod_lat=atand((RE)^2/(RP)^2*tand(geoc_lat));
    lat=geod_lat;
else
    lon=NaN;
    lat=NaN;

end

