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



%% algorithm
c1=zeros(3,npos);c2=zeros(3,npos);c3=zeros(3,npos);c4=zeros(3,npos);
for i=1:npos
    c1(:,i)=pos_ecef(:,i)/norm(pos_ecef(:,i));
    c2(:,i)=vel_ecef(:,i)/norm(vel_ecef(:,i));
    c3(:,i)=cross(c2(:,i),c1(:,i));
    c4(:,i)=cross(c2(:,i),c3(:,i));
end

D=zeros(3,3,npos);
for i=1:npos
    D(:,:,i)=[c2(:,i) c3(:,i) c4(:,i)];
end
clear c1 c2 c3 c4


% w è il vettore mounting, non cambia nel tempo perchè il come è montato
% sul satellite non cambia
fprintf('Please select the target of the LOS:\n');
obs_mode=input('1) L''Aquila, Italy \n2) Bursa, Turkey \n3) Brig, Switzerland \n4) Chamonix, France \n5) manual\n');
switch obs_mode
    case 1
        w=[0;-1.65;12.05];
    case 2
        w=[0;-14.6;55.7];
    case 3
        w=[0;-23;-23];
    case 4
        w=[0;-17.7;-27.4];
    case 5
        w1=input('rotation about the yaw axis w1 [deg]: ');
        w2=input('rotation about the pitch axis w2 [deg]: ');
        w3=input('rotation about the roll axis w3 [deg]: ');
        w=[w1;w2;w3];
end
w=repmat(w,1,npos);
m=zeros(3,npos);
for i=1:npos
    m(:,i)=[cosd(w(1,i))*sind(w(2,i))*cosd(w(3,i))+sind(w(1,i))*sind(w(3,i));...
            sind(w(1,i))*sind(w(2,i))*cosd(w(3,i))-cosd(w(1,i))*sind(w(3,i));...
            cosd(w(2,i))*cosd(w(3,i))];
end

% assetto nel tempo del satellite rispetto al frame orbitale, dovrebbe
% cambiare nel tempo nella modalità in cui il sat deve sempre puntare la gs
roll=zeros(1,npos);
pitch=zeros(1,npos);
yaw=zeros(1,npos);

SY=sind(yaw);CY=cosd(yaw);
SP=sind(pitch);CP=cosd(pitch);
SR=sind(roll);CR=cosd(roll);

M=zeros(3,3,npos);r=zeros(3,npos);g=zeros(3,npos);
for i=1:npos
    M(:,:,i)=[CY(i)*CP(i) CY(i)*SP(i)*SR(i)-SY(i)*CR(i) CY(i)*SP(i)*CR(i)+SY(i)*SR(i);...
              SY(i)*CP(i) SY(i)*SP(i)*SR(i)+CY(i)*CR(i) SY(i)*SP(i)*CR(i)-CY(i)*SR(i);...
              -SP(i) CP(i)*SR(i) CP(i)*CR(i)];
    r(:,i)=M(:,:,i)*m(:,i);
    g(:,i)=D(:,:,i)*r(:,i);
end
clear SY SP SR CY CP CR



%% load dem
% trovare un modo alternativo per scaricare il dem adeguato (non tutto il mondo) in maniera
% automatica
% prendere da qua i dem:
% https://prism-dem-open.copernicus.eu/pd-desk-open-access/prismDownload/COP-DEM_GLO-30-DGED__2023_1/Copernicus_DSM_10_N42_00_E013_00.tar
% cambia i dati del link con floor(lat) e floor(lon) del punto di
% intersezione
switch obs_mode
    case 1
        [AA,RA]=readgeoraster('../dem_files/Copernicus_DSM_10_N42_00_E013_00_DEM.tif','OutputType','double');
        hh0=2500;
    case 2
        [AA,RA]=readgeoraster('../dem_files/Copernicus_DSM_10_N40_00_E029_00_DEM.tif','OutputType','double');
        hh0=1200;
    case 3
        [AA,RA]=readgeoraster('../dem_files/Copernicus_DSM_10_N46_00_E007_00_DEM.tif','OutputType','double');
        hh0=3000;
    case 4
        [AA,RA]=readgeoraster('../dem_files/Copernicus_DSM_10_N45_00_E006_00_DEM.tif','OutputType','double');
        hh0=3500;
    case 5
        lat_obs=input('Approx latitude to observe [deg]: ');
        lon_obs=input('Approx longitude to observe [deg]: ');
        lat_obs=floor(lat_obs);
        lon_obs=floor(lon_obs);

        str_dem=strcat('../dem_files/Copernicus_DSM_10_N',num2str(lat_obs),'_00_E',num2str(lon_obs,'%03.f'),'_00_DEM.tif');
        [AA,RA]=readgeoraster(str_dem,'OutputType','double');
        hh0=1000;
end
lat_vector=linspace(RA.LatitudeLimits(2),RA.LatitudeLimits(1),RA.RasterSize(1));
lon_vector=linspace(RA.LongitudeLimits(1),RA.LongitudeLimits(2),RA.RasterSize(2));

% generation of grid
ZZ=AA;
[XX,YY]=meshgrid(lon_vector,lat_vector);



%% we are interested in one time instance, t=0
pos_ecef=pos_ecef(:,1);
g=g(:,1);
pos_lla=pos_lla(:,1);


%% new part (to adjust for number of possible intersections)
% stats over dem square (1x1 deg)
minAA=min(AA,[],'all');
maxAA=max(AA,[],'all');
meanAA=mean(AA,'all');
medianAA=median(AA,'all');
stats_dem=[minAA maxAA meanAA medianAA];



hh0_mult=(minAA:20:maxAA);
n_mult=length(hh0_mult);

hh1_mult=zeros(1,n_mult);
for i1=1:n_mult
    % algorithm that searches for the intersection point with ellipsoid
    [geod_lat,lon,~,~]=los_algorithm(hh0_mult(i1),pos_ecef,g);

    % interpolation of raster data, and newfound altitude
    hh1_mult(i1)=interp2(XX,YY,ZZ,lon,geod_lat);
end

% number of intersections
yn=hh1_mult>hh0_mult;
supp1=abs(diff(yn));
numb_inters=numel(find(supp1==1));

ind_hh0=find(yn==1,1,'last');
hh0=hh0_mult(ind_hh0);


%% larger ellipsoid (wgs84+hh)
eps=0.1;  %0.1 meter
dist=1e4;
hh=hh0;
k=0;
int_table=zeros(4,1);
while dist>eps
    k=k+1;

    % algorithm that searches for the intersection point with ellipsoid
    [geod_lat,lon,~,~]=los_algorithm(hh,pos_ecef,g);
    int_table(1:3,k)=[geod_lat lon hh]';
    
    if k>1
        dist=distance(int_table(1,k-1),int_table(2,k-1),int_table(1,k),int_table(2,k),wgs84);
    end
        
    % interpolation of raster data, and newfound altitude
    hh=interp2(XX,YY,ZZ,lon,geod_lat);
    int_table(4,k)=hh;

    % exit condition from while loop
    if dist<eps || k>49
        fprintf('dist = %7.4f m\n',dist)
        break
    end
end
hhf=int_table(3,end);
clear A B C us geoc_lat



%% only wgs84
% algorithm that searches for the intersection point with wgs84
[geod_lat,lon,~,~]=los_algorithm(0,pos_ecef,g);
int_tables=[geod_lat lon 0 interp2(XX,YY,ZZ,lon,geod_lat)]';


%% only wgs84+hh(+-500m)
hhs=[min(hh0,hhf)-500 max(hh0,hhf)+500];
for k=1:length(hhs)
    hh=hhs(k);
    [geod_lat,lon,~,~]=los_algorithm(hh,pos_ecef,g);
    int_tables(:,1+k)=[geod_lat lon hh interp2(XX,YY,ZZ,lon,geod_lat)]';
end

% map profile height for figures
dist1=distance(int_tables(1,1),int_tables(2,1),int_tables(1,2),int_tables(2,2),wgs84);
dist2=distance(int_tables(1,1),int_tables(2,1),int_tables(1,3),int_tables(2,3),wgs84);
if dist1>dist2
    limit_lat=[int_tables(1,2) int_tables(1,1)];
    limit_lon=[int_tables(2,2) int_tables(2,1)];
else
    limit_lat=[int_tables(1,3) int_tables(1,1)];
    limit_lon=[int_tables(2,3) int_tables(2,1)];
end
[zq,distq,latq,lonq]=mapprofile(AA,RA,limit_lat,limit_lon);

los_xy=zeros(1,size(int_table,2));
los_z=int_table(3,:);
for k=1:size(int_table,2)
    los_xy(k)=distance(limit_lat(1),limit_lon(1),int_table(1,k),int_table(2,k),wgs84);
end




%% some figures
% webmap
v3=webmap('World Street Map');
wmmarker(pos_lla(1),pos_lla(2),'FeatureName','SAT_SSP','Color','r');
wmmarker(int_tables(1,1),int_tables(2,1),'FeatureName','WGS84','Color','g');
wmmarker(int_tables(1,2),int_tables(2,2),'FeatureName','hh=-500m','Color','g');
wmmarker(int_tables(1,3),int_tables(2,3),'FeatureName','hh=+500m','Color','g');
for k=1:size(int_table,2)
    str_name=['h_',num2str(k-1)];
    wmmarker(int_table(1,k),int_table(2,k),'FeatureName',str_name,'Color','b');
end


% dem (1x1 deg square tile)
figure(1)
surf(XX,YY,ZZ,'EdgeColor','none')
title('DEM square tile (1 arcsec)')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
zlabel('Altitude wrt WGS84 (m)')
demcmap(ZZ)
daspect([1 1 7000])
hold on
contour3(XX,YY,ZZ,'LevelList',0,'Color','k')
scatter3(int_table(2,end),int_table(1,end),int_table(3,end),'o','MarkerFaceColor','b','MarkerEdgeColor','b');
hold off


% webmap 2
figure(2)
geoplot(latq,lonq,'b','LineWidth',1.5)
hold on
for k=1:size(int_tables,2)
    geoscatter(int_tables(1,k),int_tables(2,k),'o','MarkerFaceColor','g','MarkerEdgeColor','g')
end
for k=1:size(int_table,2)
    geoscatter(int_table(1,k),int_table(2,k),'o','MarkerFaceColor','b','MarkerEdgeColor','b')
end
geobasemap satellite
title('Intersection with DEM')


% profile height
figure(3)
plot(distq,zq,'LineWidth',2)
hold on
% los intersections
for k=1:size(int_table,2)
    plot(distq,int_table(3,k)*ones(length(distq)),'k--')
end
% los vector
ax=gca;
ylim=ax.YLim;
m=(los_z(2)-los_z(1))/(los_xy(2)-los_xy(1));
q=los_z(1)-m*los_xy(1);
f= @(y) y/m-q/m;
line([f(ylim(1)) f(ylim(2))],[ylim(1) ylim(2)],'Color','r','LineWidth',2)
for k=1:size(int_table,2)
    scatter(los_xy(k),los_z(k),120,'o','MarkerFaceColor','k','MarkerEdgeColor','k')
    text(los_xy(k),los_z(k),string(k),'Color','white','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8)
end
hold off
xlim([0 distq(end)])
xlabel('Range (m)')
ylabel('Elevation (m)')
title('Profile height of the intersection between LOS and DEM')


int_tables
int_table

dlat=abs(int_tables(1,1)-int_table(1,end));
dlon=abs(int_tables(2,1)-int_table(2,end));
ddis=distance(int_tables(1,1),int_tables(2,1),int_table(1,end),int_table(2,end),wgs84);
fprintf('dlat = %7.4e arcsec\n',dlat*3600);
fprintf('dlon = %7.4e arcsec\n',dlon*3600);
fprintf('distance = %8.4f meters\n',ddis);

% confronta dati finali con sito web:
% https://www.dcode.fr/earth-elevation
% https://www.gpsvisualizer.com/elevation

% calcolare errore dall'articolo
% manca il discorso sui pixel, e dove cade il los nella matrice (l,m)