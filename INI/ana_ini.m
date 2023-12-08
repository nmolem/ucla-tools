
  gridfile = '/batavia/MOUNT/mount0_grd.nc';
  inifile  = '/batavia/MOUNT/mount0_ini.nc';
  gridfile = '/batavia/nmolem/PTEST/ADJUST/testp32_grd.nc';
  inifile  = '/batavia/nmolem/PTEST/ADJUST/testp32a_ini.nc';
  gridfile = '/batavia/nmolem/PTEST/ADJUST/shelf48_grd.nc';
  inifile  = '/batavia/nmolem/PTEST/ADJUST/shelf48_ini.nc';
  gridfile = '/batavia/MOUNT/mount128_grd.nc';
  inifile  = '/batavia/MOUNT/mount128_ini.nc';

  profile =  '/batavia/TANEY/profile';
  load(profile)

  nc = netcdf(gridfile);
  lon = nc{'lon_rho'}(:);
  lat = nc{'lat_rho'}(:);
  h   = nc{'h'}(:);
  close(nc)
  [ny,nx] = size(h);

  [x,y] = ll2utm(lon,lat);

   % Child and parent s-coord parameters into chdscd and parscd

   N = 32;
   theta_s = 6.5;
   theta_b = 3.5;
   hc      = 400;
   chdscd.theta_s = 6.5;
   chdscd.theta_b = 3.5;
   chdscd.hc      = 400;
   if ~exist(inifile)
     r2r_create_ini(inifile, gridfile, N, chdscd, 'clobber')
   end

   [zw,Cs_w] = zlevs3(h,0*h,theta_s,theta_b,hc,N,'w','new2008');
   [zr,Cs_r] = zlevs3(h,0*h,theta_s,theta_b,hc,N,'r','new2008');

   jdx = [1:ny];
   idx = [1:nx];
   [x,y] = meshgrid(idx,jdx);
   y0 = min(min(y));
   y1 = max(max(y));


   nc = netcdf(inifile,'w');
   nc{'ocean_time'}(:) = 0;
   nc{'zeta'}(:) = 0.2*(y-y0)/(y1-y0) - 0.1;
   nc{'zeta'}(:) = 0.0;
   nc{'u'}(:) = 0.0;
   nc{'v'}(:) = 0.0;
%  nc{'temp'}(:) =  3.5;
%  nc{'salt'}(:) = 34.0;
   for i = 1:nx
    for j = 1:ny
%     nc{'temp'}(:,j,i) = interp1(zi,t_av,zr(:,j,i));
%     nc{'salt'}(:,j,i) = interp1(zi,s_av,zr(:,j,i));
      nc{'temp'}(1,:,j,i) = interp1(zi,t_av,zr(:,j,i));
      nc{'salt'}(1,:,j,i) = interp1(zi,s_av,zr(:,j,i));
    end
   end
   nc{'salt'}(1,:,:,:) = 33 + 0*zr;
   nc{'temp'}(1,:,:,:) = 20 + 18*(zr/5000).^3; % + 0*tanh(5e-5*(xr-50000));


   nc.Cs_w = Cs_w;
   nc.Cs_r = Cs_r;

   close(nc)




 %%  f v_z =  g*b_x



