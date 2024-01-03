function make_grid(nx,ny,lon,lat,pn,pm,hraw,angle,xsize,ysize,rot,tra_lon,tra_lat,lone,late);
% 
% This is part of Easy Grid
%  (c) 2008, Jeroen Molemaker, UCLA
%
grdname = 'roms_grd.nc';
ROMS_title = ['ROMS grid by Easy Grid. Settings:', ...
    ' nx: '   ,num2str(nx),    ' ny: '   , num2str(ny), ...
    ' xsize: ', num2str(xsize/1e3),' ysize: ',num2str(ysize/1e3), ...
    ' rotate: ',num2str(rot),' Lon: ',num2str(tra_lon),' Lat: ',num2str(tra_lat) ];
%
nxp= nx+2;
nyp= ny+2;
%
% Create the grid file
%
if exist(grdname)
 !rm roms_grd.nc
end
create_grid(nxp,nyp,grdname,ROMS_title)
%
%
f=4*pi*sin(lat)/(24*3600);
%fmax = max(max(f));
%fmin = min(min(f));
%
% Compute the mask
%
mask = 0*hraw + 1;
mask(hraw > 0) = 0;
%
% Fill the grid file
%
ncwrite(grdname,'pm',pm');
ncwrite(grdname,'pn',pn');
ncwrite(grdname,'angle',angle')
ncwrite(grdname,'hraw',hraw')
ncwrite(grdname,'f',f')
ncwrite(grdname,'mask_rho',mask')
ncwrite(grdname,'lon_rho',lon'*180/pi) %% (degrees)
ncwrite(grdname,'lat_rho',lat'*180/pi)  %% (degrees)
ncwrite(grdname,'spherical','T')
ncwrite(grdname,'tra_lon',tra_lon)
ncwrite(grdname,'tra_lat',tra_lat)
ncwrite(grdname,'rotate',rot)
%

