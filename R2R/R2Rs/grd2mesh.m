function [lon,lat] = grd2mesh(lon,lat,tra_lon,tra_lat,rotate,xy_flip,parent);

 d2r = pi/180;
 r2d = 180/pi;

 lon= lon*d2r;
 lat= lat*d2r;
 lon = lon - tra_lon*pi/180;
 [lon,lat] = tra_sphere(lon,lat,-tra_lat);
 [lon,lat] = rot_sphere(lon,lat,-rotate);
  if xy_flip
   lon = flipdim(lon',1);
   lat = flipdim(lat',1);
   [lon,lat] = rot_sphere(lon,lat,-90);
  end


 if parent
  lon1 = mean(lon,1);
  lat1 = mean(lat,2);
  [lon2,lat2] = meshgrid(lon1,lat1);
  if max(max(abs(lon2-lon)))>1e-6 
    err = (lon2 - lon)./lon;
    imagesc(err);set(gca,'ydir','normal');colorbar
    error 'mesgrid and original dont match; consult Jeroen.'
  end
  lon= lon2*r2d;
  lat= lat2*r2d;
 else
  lon= lon*r2d;
  lat= lat*r2d;
 end

