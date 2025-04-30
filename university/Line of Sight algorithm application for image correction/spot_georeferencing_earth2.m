clear all
clc
close all hidden
wmclose('all')
format long g

%% earth model (wgs84)
wgs84=wgs84Ellipsoid('m');
RE=wgs84.SemimajorAxis; % Earth equatorial radius WGS84
RP=wgs84.SemiminorAxis; % Earth polar radius WGS84

%% satellite scenario
% time propagation
startTime=datetime(2024,7,16,10,15,54);
stopTime=startTime+hours(1);
sampleTime=1;
sc=satelliteScenario(startTime,stopTime,sampleTime);

% satellite info
sat=satellite(sc,'leoSatellite_S2.tle',Name='SAT');
sat.Visual3DModel='smallsat.glb';
sat.Visual3DModelScale=1e5;
ax=coordinateAxes(sat,'Scale',2);


% ground station info
lat_gs=dms2degrees([42 24 17]);
lon_gs=dms2degrees([12 51 26]);
gs=groundStation(sc,lat_gs,lon_gs,Name='Rieti GS');

% pointing sat at nadir
pointAt(sat,'nadir');

% coverage of satellite
camSensor1=conicalSensor(sat,MaxViewAngle=60);
ac=access(camSensor1,gs);
t=accessIntervals(ac);

% figure for satellite scenario
% v1=satelliteScenarioViewer(sc,'CameraReferenceFrame','ECEF','Name','SSV_ECEF');
% v2=satelliteScenarioViewer(sc,'CameraReferenceFrame','Inertial','Name','SSV_ECI');
fieldOfView(camSensor1,'LineColor','g');
groundTrack(sat);


%% satellite state
[pos_eci,vel_eci]=states(sat,'CoordinateFrame','inertial');
[pos_ecef,vel_ecef]=states(sat,'CoordinateFrame','ecef');
[pos_lla,vel_lla]=states(sat,'CoordinateFrame','geographic');
[ncoord,npos]=size(pos_eci);

% maximum coverage achievable at given heigth
orb_rad=zeros(1,npos);
Earth_radius=zeros(1,npos);
for i=1:npos
    orb_rad(i)=norm(pos_eci(:,i));
    Earth_radius(i)=geocradius(pos_lla(1,i),'WGS84');
end
height=orb_rad-Earth_radius;
heightm=mean(height);
rho=asind(RE/(RE+heightm));
camSensor2=conicalSensor(sat,MaxViewAngle=2*rho);
fieldOfView(camSensor2,'LineColor','r');



%% los vector

nlos=3;
w=zeros(3,nlos);w(3,:)=linspace(-36,36,nlos);

m=zeros(3,nlos);
for i=1:nlos
    m(:,i)=[cosd(w(1,1))*sind(w(2,i))*cosd(w(3,i))+sind(w(1,i))*sind(w(3,i));...
            sind(w(1,i))*sind(w(2,i))*cosd(w(3,i))-cosd(w(1,i))*sind(w(3,i));...
            cosd(w(2,i))*cosd(w(3,i))];
end

rpy=zeros(3,npos);      % roll, pitch, yaw axis rotations along the orbital motion,
                        % should change following the gs, but alas is nadir
                        % orientation here!


g=zeros(3,nlos,npos);
for i=1:npos
    g(:,:,i)=PrepAlgor(pos_ecef(:,1),vel_ecef(:,1),m,rpy(:,1));
end



%% we are interested in one time instance, t=0
pos_ecef=pos_ecef(:,1);
pos_lla=pos_lla(:,1);
g=g(:,:,1);




%% load DEM

% intersection with Earth ellipsoid (h=0)
int_table0=cell(nlos,1);
geod_lat=zeros(1,nlos);lon=zeros(1,nlos);

for i=1:nlos
    % algorithm that searches for the intersection point with Earth
    % ellipsoid (h=0)
    [geod_lat(i),lon(i)]=los_algorithm(0,pos_ecef,g(:,i));
    int_table0{i}=[geod_lat(i) lon(i) 0];
end


% search for DEM tiles to load
DEMtiles=unique([floor(lon);floor(geod_lat)]','rows');


figure(1);
plot(lon,geod_lat,'r');
grid minor
axis equal
hold on
for i=1:length(DEMtiles)
    rectangle('Position',[DEMtiles(i,1) DEMtiles(i,2) 1 1]);
end

xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
title('');


% load merged dem
load('mergedDEM.mat');

lat_vector=linspace(RA.LatitudeLimits(2),RA.LatitudeLimits(1),RA.RasterSize(1));
lon_vector=linspace(RA.LongitudeLimits(1),RA.LongitudeLimits(2),RA.RasterSize(2));

% generation of grid
ZZ=A;
[XX,YY]=meshgrid(lon_vector,lat_vector);






%% new part (to adjust for number of possible intersections)
% stats over dem square (1x1 deg)
minA=min(A,[],'all');
maxA=max(A,[],'all');
% meanA=mean(A,'all');
% medianA=median(A,'all');
% stats_dem=[minA maxA meanA medianA];



hh0_mult=(minA:50:maxA);
n_mult=length(hh0_mult);

hh1_mult=zeros(nlos,n_mult);
hh0=zeros(nlos,1);
for i=1:nlos
    for ii=1:n_mult
        % algorithm that searches for the intersection point with ellipsoid
        [geod_lat,lon,~,~]=los_algorithm(hh0_mult(ii),pos_ecef,g(:,i));
    
        % interpolation of raster data, and newfound altitude
        hh1_mult(i,ii)=interp2(XX,YY,ZZ,lon,geod_lat);
    end
    
    % number of intersections
    yn=hh1_mult(i,:)>hh0_mult;
    supp1=abs(diff(yn));
    numb_inters=numel(find(supp1==1));
    
    ind_hh0=find(yn==1,1,'last');
    hh0(i)=hh0_mult(ind_hh0);
end


%% larger ellipsoid (wgs84+hh)
eps=0.1;  %0.1 meter
int_table=cell(nlos,50);
hhf=zeros(size(hh0));
klast=zeros(size(hh0));

for i=1:nlos
dist=1e4;
hh=hh0(i);
k=0;

    while dist>eps
        k=k+1;
    
        % algorithm that searches for the intersection point with ellipsoid
        [geod_lat,lon,~,~]=los_algorithm(hh,pos_ecef,g(:,i));
        int_table{i,k}=[geod_lat lon hh];
        
        if k>1
            dist=distance(int_table{i,k-1}(1),int_table{i,k-1}(2),int_table{i,k}(1),int_table{i,k}(2),wgs84);
        end
            
        % interpolation of raster data, and newfound altitude
        hh=interp2(XX,YY,ZZ,lon,geod_lat);
    
        % exit condition from while loop
        if dist<eps || k>49
            fprintf('dist = %7.4f m\n',dist)
            break
        end
    end
    klast(i)=k;
    hhf(i)=int_table{i,k}(3);
end
hhf=int_table(3,end);






% % map profile height for figures
% dist1=distance(int_tables(1,1),int_tables(2,1),int_tables(1,2),int_tables(2,2),wgs84);
% dist2=distance(int_tables(1,1),int_tables(2,1),int_tables(1,3),int_tables(2,3),wgs84);
% if dist1>dist2
%     limit_lat=[int_tables(1,2) int_tables(1,1)];
%     limit_lon=[int_tables(2,2) int_tables(2,1)];
% else
%     limit_lat=[int_tables(1,3) int_tables(1,1)];
%     limit_lon=[int_tables(2,3) int_tables(2,1)];
% end
% [zq,distq,latq,lonq]=mapprofile(AA,RA,limit_lat,limit_lon);
% 
% los_xy=zeros(1,size(int_table,2));
% los_z=int_table(3,:);
% for k=1:size(int_table,2)
%     los_xy(k)=distance(limit_lat(1),limit_lon(1),int_table(1,k),int_table(2,k),wgs84);
% end




%% some figures
% dem (1x1 deg square tile)
figure(2)
surf(XX,YY,ZZ,'EdgeColor','none')
title('DEM combined tiles (3 arcsec)')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
zlabel('Altitude wrt WGS84 (m)')
demcmap(ZZ)
daspect([1 1 7000])
hold on
for i=1:nlos
    scatter3(int_table{i,klast(i)}(2),int_table{i,klast(i)}(1),int_table{i,klast(i)}(3),'o','MarkerFaceColor','r','MarkerEdgeColor','r');
end
hold off



%% statistics about analysis
dlat=zeros(1,nlos);
dlon=zeros(1,nlos);
ddis=zeros(1,nlos);

for i=1:nlos
    dlat(i)=abs(int_table0{i}(1)-int_table{i,klast(i)}(1));
    dlon(i)=abs(int_table0{i}(2)-int_table{i,klast(i)}(2));
    ddis(i)=distance(int_table0{i}(1),int_table0{i}(2),int_table{i,klast(i)}(1),int_table{i,klast(i)}(2),wgs84);
end

dlat_avg=mean(dlat);
dlon_avg=mean(dlon);
ddis_avg=mean(ddis);

fprintf('dlat_avg = %7.4e arcsec\n',dlat_avg*3600);
fprintf('dlon_avg = %7.4e arcsec\n',dlon_avg*3600);
fprintf('distance_avg = %8.4f meters\n',ddis_avg);



% confronta dati finali con sito web:
% https://www.dcode.fr/earth-elevation
% https://www.gpsvisualizer.com/elevation










function g=PrepAlgor(pos_ecef,vel_ecef,m,rpy)
    
    % rotation matrix between orbital and ecef frames
    c1=pos_ecef/norm(pos_ecef);
    c2=vel_ecef/norm(vel_ecef);
    c3=cross(c2,c1);
    c4=cross(c2,c3);
    
    D=[c2 c3 c4];
    
    
    % rotation matrix between body and orbital frames
    SY=sind(rpy(3));CY=cosd(rpy(3));
    SP=sind(rpy(2));CP=cosd(rpy(2));
    SR=sind(rpy(1));CR=cosd(rpy(1));
    
    M=[CY*CP CY*SP*SR-SY*CR CY*SP*CR+SY*SR;...
       SY*CP SY*SP*SR+CY*CR SY*SP*CR-CY*SR;...
      -SP    CP*SR          CP*CR];
    
    
    % LOS vector in ecef frame
    g=D*M*m;

end




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

end
