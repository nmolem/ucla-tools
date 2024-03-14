
gridfile = '/batavia/NWATL/nwat_grd.nc';
gridfile = '/batavia/nmolem/PROJ/ROMS/CL/CL3D/mont_grd.nc';
gridfile = '/batavia/NATL10/natl10_grd.nc';
gridfile = 'pacmed_grd.nc';
gridfile = '/avatar/nmolem/SPLASH/splash_grd.nc';
gridfile = '/avatar/nmolem/LBIGHT/lbight_grd.nc';
gridfile = '/avatar/nmolem/WESTC/westc_grd.nc';
gridfile = '/avatar/nmolem/LUZON/luzon_grd.nc';
gridfile = '/avatar/nmolem/ERA5_grd.nc';
gridfile = '/avatar/nmolem/NEPAC/nepac_grd_mod.nc';
gridfile = '/avatar/nmolem/GREEN/green_grd.nc';
gridfile = '/avatar/nmolem/TWAKE/twake_grd.nc';
gridfile = '/avatar/nmolem/palau_bathy.nc';
gridfile = '/avatar/nmolem/SFBAY/sfbay_grd.nc';

 mask = ncread(gridfile,'mask_rho');
%mask = ncread(gridfile,'mask');
[ny,nx] = size(mask);

reg = bwlabel(mask,4);

lint =  0; %% size of largest region
lreg = 0; %% number of largest region
nreg = max(max(reg)); %% number of regions
for i = 1:nreg
  int = sum(sum(reg==i));
  if int>lint
    lreg = i;
    lint = int;
  end
end


for ireg = 1:nreg
  if ireg~=lreg
    %% check before setting to zero
    int = sum(sum(reg==ireg));
    if int>nx*ny/10
      disp(['region: ' num2str(ireg) 'is large.'])
    else
     mask(reg==ireg) = 0;
    end
  end
end

if 0
  subplot(1,2,1)
  imagesc(reg);axis xy;colorbar
  subplot(1,2,2)
  imagesc(mask);axis xy;colorbar
end
  ncwrite(gridfile,'mask_rho',mask);


