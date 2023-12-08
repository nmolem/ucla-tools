

gname = '/avatar/nmolem/SPLASH/splash_grd.nc';

mask1= ncread(gname,'mask_rho');
lon  = ncread(gname,'lon_rho');
lat  = ncread(gname,'lat_rho');

lon1 = lon;
lat1 = lat;


ol1_lon = [lon(1,1) lon(1,end) lon(end,end) lon(end,1) lon(1,1)];
ol1_lat = [lat(1,1) lat(1,end) lat(end,end) lat(end,1) lat(1,1)];

gname = 'roms_grd.nc'

lon  = ncread(gname,'lon_rho');
lat  = ncread(gname,'lat_rho');


ol2_lon = [lon(1,1) lon(1,end) lon(end,end) lon(end,1) lon(1,1)];
ol2_lat = [lat(1,1) lat(1,end) lat(end,end) lat(end,1) lat(1,1)];


close all
pcolor(lon1,lat1,mask1);shading interp
hold on
plot(ol1_lon,ol1_lat);
plot(ol2_lon,ol2_lat);
hold off
xlim([-90.3 -89.0])
ylim([ 28.6  29.35])
