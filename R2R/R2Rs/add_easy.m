
%grdfile = '/batavia/SMCC/smcc_grd_mod.nc';
grdfile = '/batavia/NWATL/nwatl_grd.nc';

nc = netcdf(grdfile,'w');

nc.title(:)

 rotate  =  35 
 tra_lon =- 66.0 
 tra_lat =  28.4 
 xy_flip = 0 


nc{'tra_lon'} = ncdouble('one');
nc{'tra_lon'}.long_name = 'easy grid longitudinal translation of base grid';
nc{'tra_lon'}.units = 'degrees';

nc{'tra_lat'} = ncdouble('one');
nc{'tra_lat'}.long_name = 'easy grid latitudinal translation of base grid';
nc{'tra_lat'}.units = 'degrees';

nc{'rotate'} = ncdouble('one');
nc{'rotate'}.long_name = 'easy grid latitudinal translation of base grid';
nc{'rotate'}.units = 'degrees';

nc{'xy_flip'} = ncint('one');
nc{'xy_flip'}.long_name = 'easy grid xy flip of base grid';
nc{'xy_flip'}.units = 'True/False (0/1)';

nc{'tra_lon'}(:) = tra_lon;
nc{'tra_lat'}(:) = tra_lat;
nc{'rotate'}(:)  = rotate;
nc{'xy_flip'}(:) = xy_flip;

close(nc)


