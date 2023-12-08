% fill_frc_rivers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Put in the daily river runoff data in the 'new' format. 
%
%  2022, Jeroen Molemaker (UCLA)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
%

if 1
 clear all
 close all
 wrk_dir = '/paracas/nmolem/NEPAC/';
 grdname   = [wrk_dir, 'nepac_grd.nc'];
 frcname   = [wrk_dir, 'nepac_riv.nc'];
 wrk_dir = '/paracas/nmolem/LUZON/';
 grdname   = [wrk_dir, 'luzon_grd.nc'];
 frcname   = [wrk_dir, 'luzon_riv.nc'];
 wrk_dir = '/paracas/nmolem/GREEN/';
 grdname   = [wrk_dir, 'green_grd.nc'];
 frcname   = [wrk_dir, 'green_riv.nc'];
 wrk_dir = '/paracas/nmolem/WESTC/';
 grdname   = [wrk_dir, 'westc_grd.nc'];
 frcname   = [wrk_dir, 'westc_riv.nc'];
 wrk_dir = '/paracas/nmolem/SMODE/';
 grdname   = [wrk_dir, 'smode_grd.nc'];
 frcname   = [wrk_dir, 'smode_riv.nc'];
%wrk_dir = '/paracas/nmolem/SFBAY/';
%grdname   = [wrk_dir, 'sfbay_grd.nc'];
%frcname   = [wrk_dir, 'sfbay_riv.nc'];
end
%
%  runoff climatology file names:
%
runoff_dir = '/paracas/DATASETS/';
runoff_data  = [runoff_dir,'coastal-stns-Vol-monthly.updated-May2019.nc'];

mask = ncread(grdname,'mask_rho');
[nx ny] = size(mask);

% create entry for river distribution in the grid file if needed
 novar = 0;
 ncid = netcdf.open(grdname);
 try
   id = netcdf.inqVarID(ncid,'river_flux'); 
 catch exception
   novar = 1; 
 end
 netcdf.close(ncid);

 if novar
   nccreate(grdname,'river_flux','Dimensions',{'xi_rho',nx,'eta_rho',ny},'datatype','single');
   ncwriteatt(grdname,'river_flux','long_name','river volume flux partition');
   ncwriteatt(grdname,'river_flux','units','none');
   rflx = 0*mask;
   ncwrite(grdname,'river_flux',rflx);
 end

%
%%%%%%%%%%%%%%%%%%% END USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% disp(' ')
% disp(' Read in the grid_data')
 lon  = ncread(grdname,'lon_rho');
 lat  = ncread(grdname,'lat_rho');
 dxi  = ncread(grdname,'pm');
 dyi  = ncread(grdname,'pn');
 mask = ncread(grdname,'mask_rho');
 h    = ncread(grdname,'hraw');
 dx  = 1./dxi;
 dy  = 1./dyi;
 lon = mod(lon,360);

 [nx, ny] = size(lon);
 

 if 1 % reset river flux if needed
   rflx = 0*mask;
   ncwrite(grdname,'river_flux',rflx);
 end

% find the appropriate location of the river in the grid.

if 1
 lon_frc = ncread(runoff_data,'lon_mou');
 lat_frc = ncread(runoff_data,'lat_mou');
 tim_frc = ncread(runoff_data,'time');
 flx_frc = ncread(runoff_data,'FLOW');
 m2s_frc = ncread(runoff_data,'ratio_m2s');
 riv_frc = ncread(runoff_data,'riv_name');
 vol_frc = ncread(runoff_data,'vol_stn');
 [nr,nt] = size(flx_frc);
 lon_frc = mod(lon_frc,360);
 else
  % create your own rivers
  nt = 360; % once per day
  nr =   1; % 1 is plenty
  lon_frc(1) = -120.632;
  lat_frc(1) =  34.55;
  tim_frc = [0.5:1:359.5]; % time in days
  flx_frc = 2e3*ones(360,1); % m3/s
 end



 d2r = pi/180;
 iloc = 0*lon_frc;
 jloc = 0*lon_frc;
 mindist = 0*lon_frc;
 for i = 1:length(lon_frc)
   dist = gc_dist(lon*d2r,lat*d2r,lon_frc(i)*d2r,lat_frc(i)*d2r);
   mindist(i) = min(dist(:));
   if mindist(i) < mean(dx(:))
     [iloc(i),jloc(i)] = find(dist==mindist(i));
   end
 end


 % toss the rivers outside the grid
 out = mindist>mean(dx(:));
 iloc(out) = [];
 jloc(out) = [];
 mindist(out) = [];
 nriv = length(iloc) 
 riv_flx = zeros(nt,nriv);
 riv_lon = zeros(1,nriv);
 riv_lat = zeros(1,nriv);
 riv_vol = zeros(1,nriv);
 iriv = 0;
 for i = 1:nr
   if ~out(i)
     iriv = iriv+1;
     riv_flx(:,iriv) = flx_frc(i,:).*m2s_frc(i) ;
     riv_lon(iriv)   = lon_frc(i);
     riv_lat(iriv)   = lat_frc(i);
     riv_nam(:,iriv) = riv_frc(:,i) ;
     riv_vol(iriv)   = vol_frc(i).*m2s_frc(i) ;
   end
 end
 riv_tim = tim_frc;

% toss the NaN
 out=find(isnan(riv_vol)==0) ;
 nriv = length(out)
 riv_flx = riv_flx(:,out) ;
 riv_lon = riv_lon(out) ;
 riv_lat = riv_lat(out) ;
 riv_nam = riv_nam(:,out) ;
 riv_vol = riv_vol(out) ;

% Edit the 2d field with river fluxes


 for iriv = 1:nriv

   if iriv==1
    jloc(iriv) = jloc(iriv) +200;
   end
   	 
   disp(['River : ' num2str(iriv) '/' num2str(nriv) ' --> ' riv_nam(:,iriv)'])
   disp(['(lon,lat)=(' num2str(riv_lon(iriv)) ',' num2str(riv_lat(iriv)) ')'])
   disp(['Mean Flux = ' num2str(riv_vol(iriv)*1e9/(365*24*3600)) ' m3/s'])
   Cml = sum(riv_vol(1:iriv))/sum(riv_vol) ;
   disp(['Cumulative Flux = ' num2str(Cml*100) '%'])
   disp('***************')

   i0 = max(iloc(iriv)-400,1);
   i1 = min(iloc(iriv)+400,nx);
   j0 = max(jloc(iriv)-400,1);
   j1 = min(jloc(iriv)+400,ny);

   % edit the mask if needed

   % edit the river flux
   river_click

 end

 % Check the 2d field with river fluxes and partition the flux
 mask = ncread(grdname,'mask_rho');
 rflx = ncread(grdname,'river_flux');
 for iriv = 1:nriv
   iriv
   err_rflx = rflx>0&mask>0;
   if sum(err_rflx(:))
     rflx(err_rflx) = 0;
     disp('some fluxes were wrong and were removed')
   end
   faces = mask(1:end-2,2:end-1)+mask(3:end,2:end-1)+mask(2:end-1,1:end-2)+mask(2:end-1,3:end);
   rflx_sm = rflx(2:end-1,2:end-1);
   err2_rflx = rflx_sm>0&faces==0;
   if sum(err2_rflx(:))
     rflx_sm(err2_rflx) = 0;
     disp('some fluxes were wrong and were removed')
     rflx(2:end-1,2:end-1) = rflx_sm;
   end

   sum_flx = sum(rflx(:)>=iriv&rflx(:)<iriv+1) 
   if sum_flx>0
     rflx(rflx>=iriv&rflx<iriv+1) = iriv + 1./sum_flx;
   end
 end
 ncwrite(grdname,'river_flux',rflx);


 % Done with the river partitions, now write the volumes and tracer values

 % time , make seasonal cycle 
 riv_tim = datenum( round(floor(riv_tim*1e-2)) , round((riv_tim*1e-2-floor(riv_tim*1e-2))*1e2) , 15 ) ;
 riv_flx2=riv_flx;
 for t=1:12
     mnth=find(str2num(datestr(riv_tim,'mm'))==t) ;
     for m=1:length(mnth)
     riv_flx2(mnth(m),:) = nanmean(riv_flx(mnth,:),1) ;
     end
 end
 date_orig = datenum(2000,01,01) ;
 riv_tim = riv_tim - date_orig ;
 indfrst = find(riv_tim<0,1,'last') ;
 riv_tim = riv_tim(indfrst:end) ; 
 riv_flx = riv_flx2(indfrst:end,:) ;
 ntimes = length(riv_tim) ;

 nt = 2; % Temperature and Salinity only

 riv_trc = zeros(nriv,nt,ntimes);
 for iriv = 1:nriv
%  riv_trc(iriv,1,:) = 5 + 13*(cos(tim_frc*2*pi/365 + 0.85*pi)+1); % idealized annual cycle of temp
   riv_trc(iriv,1,:) = 17;
   riv_trc(iriv,2,:) = 1.0; % Salinity
 end

% create entry for river volume and tracer data in the forcing file if needed
 nofile = 0;
 try
   ncid = netcdf.open(frcname);
%  id = netcdf.inqVarID(ncid,'river_volume'); 
 catch exception
   nofile = 1; 
 end

 if nofile
   ncid = netcdf.create(frcname,'NC_NOCLOBBER');
   did = netcdf.defDim(ncid,'xi_rho',nx);
   did = netcdf.defDim(ncid,'eta_rho',ny);
   netcdf.close(ncid)

   nccreate(frcname,'river_volume','Dimensions',{'nriver',nriv,'river_time',ntimes},'datatype','single');
   ncwriteatt(frcname,'river_volume','long_name','River volume flux');
   ncwriteatt(frcname,'river_volume','units','m^3/s');

   nccreate(frcname,'river_tracer','Dimensions',{'nriver',nriv,'ntracers',nt,'river_time',ntimes},'datatype','single');
   ncwriteatt(frcname,'river_tracer','long_name','River tracer data');
   ncwriteatt(frcname,'river_tracer','units','variable');

   nccreate(frcname,'river_time','Dimensions',{'river_time',ntimes},'datatype','single');
   ncwriteatt(frcname,'river_time','long_name','river data time');
   ncwriteatt(frcname,'river_time','units','yearday');
 else
   netcdf.close(ncid);
 end

 ncwrite(frcname,'river_volume',riv_flx');
 ncwrite(frcname,'river_tracer',riv_trc);
 ncwrite(frcname,'river_time',riv_tim');

% Now, do the same for pipes

return

 pip_trc = zeros(npip,nt,ntimes);
 for iriv = 1:npip
%  riv_trc(iriv,1,:) = 5 + 13*(cos(tim_frc*2*pi/365 + 0.85*pi)+1); % idealized annual cycle of temp
   pip_trc(iriv,1,:) = 17;
   pip_trc(iriv,2,:) = 1.0; % Salinity
 end

 if 1 % create entry for river volume and tracer data in the forcing file if needed
   nccreate(frcname,'pipe_volume','Dimensions',{'nriver',nriv,'river_time',ntimes},'datatype','single');
   ncwriteatt(frcname,'pipe_volume','long_name','River volume flux');
   ncwriteatt(frcname,'pipe_volume','units','m^3/s');

   nccreate(frcname,'pipe_tracer','Dimensions',{'nriver',nriv,'ntracers',nt,'river_time',ntimes},'datatype','single');
   ncwriteatt(frcname,'pipe_tracer','long_name','River tracer data');
   ncwriteatt(frcname,'pipe_tracer','units','variable');

   nccreate(frcname,'pipe_time','Dimensions',{'river_time',ntimes},'datatype','single');
   ncwriteatt(frcname,'pipe_time','long_name','river data time');
   ncwriteatt(frcname,'pipe_time','units','yearday');
 end

 ncwrite(frcname,'pipe_volume',pip_flx');
 ncwrite(frcname,'pipe_tracer',pip_trc);
 ncwrite(frcname,'pipe_time',pip_tim');








