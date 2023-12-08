  function ang = calc_ang(lon,lat);

   lonu = 0.5*(lon(2:end,:)+lon(1:end-1,:))*pi/180;
   latu = 0.5*(lat(2:end,:)+lat(1:end-1,:))*pi/180;

 %% Compute angles of local grid positive x-axis relative to east
   dellon = lonu(2:end,:)-lonu(1:end-1,:);
   dellat = latu(2:end,:)-latu(1:end-1,:);
   dellon(dellon> pi) = dellon(dellon>pi) - 2*pi;
   dellon(dellon<-pi) = dellon(dellon<-pi)+ 2*pi;
   dellon = dellon.*cos(0.5*(latu(2:end,:)+latu(1:end-1,:)));

   ang_s = atan(dellat./(dellon+1e-16));
   ang_s(dellon<0 & dellat< 0) = ang_s(dellon<0 & dellat< 0) - pi;
   ang_s(dellon<0 & dellat>=0) = ang_s(dellon<0 & dellat>=0) + pi;
   ang_s(ang_s > pi) = ang_s(ang_s > pi) - pi;
   ang_s(ang_s <-pi) = ang_s(ang_s <-pi) + pi;

   ang = lon;
   ang(2:end-1,:) = ang_s;
   ang(1,:)   = ang(2,:);
   ang(end,:) = ang(end-1,:);


