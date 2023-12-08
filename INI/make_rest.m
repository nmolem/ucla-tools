

nc = netcdf('/batavia/SMCAL/smcal_grd_mod.nc');
h =nc{'h'}(:);
[z,Cs] = zlevs4(h,0*h,6.5,0,120,42,'r','new2008');
close(nc)

[ny,nx] = size(h);


nd = netcdf('/batavia/SMCAL/smcal_ini.nc','w');
i0 = 400;
j0 = 550;
t = nd{'temp'}(:);
s = nd{'salt'}(:);
tp = t(:,j0,i0);
sp = s(:,j0,i0);
zp = z(:,j0,i0);


 for j = 1:ny
  for i = 1:nx
    t(:,j,i) = interp1(zp,tp,squeeze(z(:,j,i)),'spline','extrap');
    s(:,j,i) = interp1(zp,tp,squeeze(z(:,j,i)),'spline','extrap');
  end
  round(j/ny*100)
 end

nd{'temp'}(:) = t;
nd{'salt'}(:) = s;
 



