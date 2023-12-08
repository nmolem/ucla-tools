
cname = 'USW4/usw4_grd_mod.nc';
pname = 'PACBIG/pacbig_grd.nc';
cname = 'sample_child_grd.nc';
pname = 'sample_grd_riv.nc';

bcname = 'onlinebc.nc';

OBC = [1,1,0,1]; % east,south,west,north

lonc = ncread(cname,'lon_rho') - 360;
latc = ncread(cname,'lat_rho');
mskc = ncread(cname,'mask_rho');

lonp = ncread(pname,'lon_rho');
latp = ncread(pname,'lat_rho');

[nxc,nyc] = size(lonc);
[nxp,nyp] = size(lonp);

ip = [0:nxp-1];
jp = [0:nyp-1];
[ip,jp] = meshgrid(ip,jp);

ip = ip';
jp = jp';



if OBC(1) % East boundary
 lonbc = lonc(end-1:end,:);
 latbc = latc(end-1:end,:);
 mskbc = mskc(end-1:end,:);


  ic = griddata(lonp,latp,ip,lonbc,latbc).*mskbc;
  jc = griddata(lonp,latp,jp,lonbc,latbc).*mskbc;
  ic(isnan(ic)) = -1;
  jc(isnan(jc)) = -1;


 nccreate(bcname,'i_east','dimensions',{'two',2,'ny',nyc},'datatype','single');
 ncwriteatt(bcname,'i_east','long_name','i-coordinate of east child boundary');
 ncwriteatt(bcname,'i_east','units','nondimensional');

 nccreate(bcname,'j_east','dimensions',{'two',2,'ny',nyc},'datatype','single');
 ncwriteatt(bcname,'j_east','long_name','j-coordinate of east child boundary');
 ncwriteatt(bcname,'j_east','units','nondimensional');

 ncwrite(bcname,'i_east',ic);
 ncwrite(bcname,'j_east',jc);

end
if OBC(2) % South boundary
  lonbc = lonc(:,1:2);
  latbc = latc(:,1:2);
  mskbc = mskc(:,1:2);

  ic = griddata(lonp,latp,ip,lonbc,latbc).*mskbc;
  jc = griddata(lonp,latp,jp,lonbc,latbc).*mskbc;
  ic(isnan(ic)) = -1;
  jc(isnan(jc)) = -1;

 nccreate(bcname,'i_south','dimensions',{'nx',nxc,'two',2},'datatype','single');
 ncwriteatt(bcname,'i_south','long_name','i-coordinate of east child boundary');
 ncwriteatt(bcname,'i_south','units','nondimensional');

 nccreate(bcname,'j_south','dimensions',{'nx',nxc,'two',2},'datatype','single');
 ncwriteatt(bcname,'j_south','long_name','j-coordinate of east child boundary');
 ncwriteatt(bcname,'j_south','units','nondimensional');

 ncwrite(bcname,'i_south',ic);
 ncwrite(bcname,'j_south',jc);

end
    
if OBC(3) % West boundary
 lonbc = lonc(1:2,:);
 latbc = latc(1:2,:);
 mskbc = mskc(1:2,:);

  ic = griddata(lonp,latp,ip,lonbc,latbc).*mskbc;
  jc = griddata(lonp,latp,jp,lonbc,latbc).*mskbc;
  ic(isnan(ic)) = -1;
  jc(isnan(jc)) = -1;

 nccreate(bcname,'i_west','dimensions',{'two',2,'ny',nyc},'datatype','single');
 ncwriteatt(bcname,'i_west','long_name','i-coordinate of east child boundary');
 ncwriteatt(bcname,'i_west','units','nondimensional');

 nccreate(bcname,'j_west','dimensions',{'two',2,'ny',nyc},'datatype','single');
 ncwriteatt(bcname,'j_west','long_name','j-coordinate of east child boundary');
 ncwriteatt(bcname,'j_west','units','nondimensional');

 ncwrite(bcname,'i_west',ic);
 ncwrite(bcname,'j_west',jc);

end
    
if OBC(4) % North boundary
  lonbc = lonc(:,end-1:end);
  latbc = latc(:,end-1:end);
  mskbc = mskc(:,end-1:end);

  ic = griddata(lonp,latp,ip,lonbc,latbc).*mskbc;
  jc = griddata(lonp,latp,jp,lonbc,latbc).*mskbc;

  ic(isnan(ic)) = -1;
  jc(isnan(jc)) = -1;

  nccreate(bcname,'i_north','dimensions',{'nx',nxc,'two',2},'datatype','single');
  ncwriteatt(bcname,'i_north','long_name','i-coordinate of east child boundary');
  ncwriteatt(bcname,'i_north','units','nondimensional');

  nccreate(bcname,'j_north','dimensions',{'nx',nxc,'two',2},'datatype','single');
  ncwriteatt(bcname,'j_north','long_name','j-coordinate of east child boundary');
  ncwriteatt(bcname,'j_north','units','nondimensional');


 
  ncwrite(bcname,'i_north',ic);
  ncwrite(bcname,'j_north',jc);

end
    
    
