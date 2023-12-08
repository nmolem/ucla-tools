

 nc = netcdf('/batavia/SMCC/smcc_grd_mod.nc');
 lonp= nc{'lon_rho'}(:);
 latp= nc{'lat_rho'}(:);
%lonp = lonp-360;
%nc.title(:);
 hp = nc{'h'}(:);
 close(nc)

 tra_lon =-125.30;
 tra_lat =  37.25;
 rotate  =  24.0;
 xy_flip = 1;


 [lonp,latp] = grd2mesh(lonp,latp,tra_lon,tra_lat,rotate,xy_flip,1);

 nc = netcdf('/batavia/TANEY/taney_grd.nc');
 lonc= nc{'lon_rho'}(:);
 latc= nc{'lat_rho'}(:);
 hc  = nc{'h'}(:);
 close(nc)

 [lonc,latc] = grd2mesh(lonc,latc,tra_lon,tra_lat,rotate,xy_flip,0);


 i0 = find(lonp(1,:)<min(min(lonc)),1,'last');
 i1 = find(lonp(1,:)>max(max(lonc)),1,'first');
 j0 = find(latp(:,1)<min(min(latc)),1,'last');
 j1 = find(latp(:,1)>max(max(latc)),1,'first');


 if 0
  nc = netcdf('/batavia/TANEY/taney_grd.nc');
  lonc= nc{'lon_rho'}(:);
  latc= nc{'lat_rho'}(:);
  close(nc)
  nc = netcdf('/batavia/SMCC/smcc_grd_mod.nc');
  lonp= nc{'lon_rho'}(:);
  latp= nc{'lat_rho'}(:);
  close(nc)
 end

%if xy_flip
% lonp = lonp(i0:i1,j0:j1);
% latp = latp(i0:i1,j0:j1);
% hp   = hp(i0:i1,j0:j1);
%else
  lonp = lonp(j0:j1,i0:i1);
  latp = latp(j0:j1,i0:i1);
  hp = flipdim(hp',1);
  hp   = hp(j0:j1,i0:i1);
%end



 if xy_flip
  hc = flipdim(hc',1);
 end

 return

 hc = interp2(lonp,latp,hp,lonc,latc,'cubic');
 hc0 = hc;
 if xy_flip
  hc = flipdim(hc,1)' ;
 end
 [hx,hy] = gradient(hc);
 dh = sqrt(hx.^2 + hy.^2);

 nc = netcdf('/batavia/TANEY/taney_grd.nc');
 lonc= nc{'lon_rho'}(:) -360;
 latc= nc{'lat_rho'}(:);
 close(nc)
 return




