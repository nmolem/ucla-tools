% fill_frc_hrrr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make a ROMS bulk forcing file using HRRR hourly data
%
%  2023, Jeroen Molemaker, Pierre Damien, UCLA
%
%  Possible future work
% - Extend lon_frc (and data) when grids straddle the prime meridian
% - Investigate over which time period the fluxes are integrated (centered or shifted)
% - Add the effect of surface pressure on Qair and other things
% - Maybe force the model with Tair and Dew point instead of Humidity
%
%%%%%%%%%%%%%%%%%%%%% USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
%
%  frc climatology file names:

frc_dir = '/paracas/DATASETS/HRRR/2016/';

% Set a date range for the forcing file
start_date = datenum(2016,09,01);
end_date   = datenum(2016,12,30);

%grdname = '/paracas/nmolem/PACHUG/pachug_grd.nc';
%disname = '/paracas/nmolem/PACHUG/pachug_cdist.mat';
%rivname = '/paracas/nmolem/PACHUG/pachug_swf.nc';
%root_name='/paracas/nmolem/PACHUG/pachug';

%grdname  = '/paracas/nmolem/NWPAC/nwpac_grd.nc';
%disname  = '/paracas/nmolem/NWPAC/nwpac_cdist.mat';
%root_name= '/paracas/nmolem/NWPAC/nwpac';

%grdname   ='/paracas/nmolem/NORMAR/normar_grd.nc';
%disname   ='/paracas/nmolem/NORMAR/normar_cdist.mat';
%root_name ='/paracas/nmolem/NORMAR/normar';

%grdname   ='/paracas/nmolem/GREEN/green_grd.nc';
%disname   ='/paracas/nmolem/GREEN/green_cdist.mat';
%root_name ='/paracas/nmolem/GREEN/green';

%grdname  = '/paracas/nmolem/NEPAC/nepac_grd.nc';
%disname  = '/paracas/nmolem/NEPAC/nepac_cdist.mat';
%root_name='/paracas/nmolem/NEPAC/nepac';

%grdname  = '/paracas/nmolem/LUZON/luzon_grd.nc';
%disname  = '/paracas/nmolem/LUZON/luzon_cdist.mat';
%root_name= '/paracas/nmolem/LUZON/luzon';

%grdname  = '/paracas/nmolem/ANGAUR/angaur_grd.nc';
%disname  = '/paracas/nmolem/ANGAUR/angaur_cdist.mat';
%root_name= '/paracas/nmolem/ANGAUR/angaur';

grdname  = '/paracas/nmolem/WESTC/westc_grd.nc';
disname  = '/paracas/nmolem/WESTC/westc_cdist.mat';
root_name= '/paracas/nmolem/WESTC/westc';

%grdname  = '/paracas/nmolem/SFBAY/sfbay_grd.nc';
%disname  = '/paracas/nmolem/SFBAY/sfbay_cdist.mat';
%root_name= '/paracas/nmolem/SFBAY/sfbay';

%grdname  = '/paracas/nmolem/SMODE/smode_grd.nc';
%disname  = '/paracas/nmolem/SMODE/smode_cdist.mat';
%root_name= '/paracas/nmolem/SMODE/smode';

%grdname  = '/paracas/nmolem/HRTEST/hrtest_grd.nc';
%root_name= '/paracas/nmolem/HRTEST/hrtest';

swcorrname = [frc_dir 'SSR_correction.nc'];


coarse_frc   = 0; % forcing files at half the resolution of the grid
add_rivers   = 0; % Adds river runoff as additional precipitation

%
%%%%%%%%%%%%%%%%%%% END USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%

%
   hrrrlist = dir([frc_dir 'HRRR_2*']);
   nfiles = length(hrrrlist);

   % find the right files
   stimes = zeros(nfiles,1);
   for i=1:nfiles
     datname = [frc_dir hrrrlist(i).name];
     stime(i) = double(ncread(datname,'time',[1],[1]))/(24*3600) + datenum(1970,1,1);
   end
   t0 = find(stime<start_date,1,'last');
   t1 = find(stime>end_date,1,'first');

   t0 = t0+1
   % trim list of filenames
   display('taking full set')
%  hrrrlist = hrrrlist(t0:t1);
   nfiles = length(hrrrlist);

   % Read in data coordinates and set data trim.
   disp(' ')
   disp(' Read in the target grid')
   if coarse_frc
%    if ~exist(grdname,'file')
       disp('adding coarse lon/lat to grid')
       add_coarse_grd(grdname);
%    end
     disp('Reading coarse grid')
     lon = ncread(grdname,'lon_coarse');
     lat = ncread(grdname,'lat_coarse');
     ang = ncread(grdname,'angle_coarse');
   else
     disp('Reading full grid')
     lon = ncread(grdname,'lon_rho');
     lat = ncread(grdname,'lat_rho');
     ang = ncread(grdname,'angle');
   end
   [nx,ny] = size(lon);

   addpath ~/matlab/m_map
   m_proj('lambert','lon',[220 300],'lat',[20 55],'clo',262.5,'par',[20.588 55.],'ell','sphere')
   mult = 1.0472044568;

   lon = mod(lon,360); %% HRRR forcing data is in the 0 to 360 longitude range


   disp(' ')
   disp(' Read in the data grid')

   datname = [frc_dir hrrrlist(1).name];
   lon_frc = ncread(datname,'longitude');
   lat_frc = ncread(datname,'latitude');
   ang_frc = calc_ang(lon_frc,lat_frc);
   [xf,yf] = m_ll2xy(lon_frc,lat_frc);
   [xg,yg] = m_ll2xy(lon,lat);

   xg = (xg-xf(1,1))*mult;
   yg = (yg-yf(1,1))*mult;

   grd.x = xg;
   grd.y = yg;

   xf = (xf-xf(1,1))*mult;
   yf = (yf-yf(1,1))*mult;

   x0 = min(xg(:));
   x1 = max(xg(:));
   y0 = min(yg(:));
   y1 = max(yg(:));

   x_frc = ncread(datname,'x');
   y_frc = ncread(datname,'y');

   if x0<min(x_frc)
    i0 = 1;
   else
    i0 = find(x_frc<x0,1,'last');
   end
   if x1>max(x_frc)
     i1 = length(x_frc);
   else
     i1 = find(x_frc>x1,1,'first');
   end

   j0 = find(y_frc<y0,1,'last');
   j1 = find(y_frc>y1,1,'first');
   fnx = i1-i0+1;
   fny = j1-j0+1;


   x_frc = x_frc(i0:i1);
   y_frc = y_frc(j0:j1);

   mask = 1 - ncread([frc_dir 'hrrr_land.nc'],'LAND_surface',[i0 j0 1],[fnx fny 1]); 
%  mask = 1 + 0*lon_frc;


   if 0 % turn off for now
   % prepare for short wave radiation correction
   corr_time = ncread(swcorrname,'time');
   swr_mult = zeros(nx,ny,12);
   for i=1:12
     corr_swr  = ncread(swcorrname,'ssr_corr',[i0 j0 i],[fnx fny 1]); 
     corr_swr = fliplr(corr_swr);
     corr_swr(mask<1) = nan;
     corr_swr = inpaint_nans(corr_swr,2);
     swr_mult(:,:,i) = interp2(lon_frc,lat_frc,corr_swr',lon,lat,'makima');
   end
   end

   if add_rivers
     % prepare for the addition of river runoff to rain data
     riv_time = ncread(rivname,'swf_time');
     riv_flux = ncread(rivname,'swf_flux');
   end

   % find the difference between the ROMS grid angles and the HRRR grid angles
   % and compute rotation factors

   ang_frc = interp2(x_frc,y_frc,ang_frc(i0:i1,j0:j1)',xg,yg,'linear');
   ang = ang-ang_frc;
%  return

   cosa = cos(ang);
   sina = sin(ang);



   % Loop over data files
   for i = 1:nfiles
     datname = [frc_dir hrrrlist(i).name];
     disp([' Processing: ' hrrrlist(i).name])

     dat_time = ncread(datname,'time');
     nrecord = length(dat_time);

     date_num = double(dat_time(1)/(24*3600)) + datenum(1970,1,1);
     label = datestr(date_num,'YYYYmmdd')

     frcname = [root_name '_frc.' label '.nc'] 

     data.datname = datname;
     data.x = x_frc;
     data.y = y_frc;
     data.i0 = i0;
     data.j0 = j0;
     data.fnx = fnx;
     data.fny = fny;
     data.mask = mask;

     if exist(frcname)
       delete(frcname)
     end
     create_frc_bulk(grdname,frcname,coarse_frc);
     ncwriteatt(frcname,'/','Data Source','HRRR (3 km nominal res)');


     for irec = 1:nrecord
        disp(['Record: ' num2str(irec)])

	% ---- time -----
	time = ncread(datname,'time',[irec],[1]);

	% translate to days since 2000,1,1
%       days = double(dat_time(1)/24) + datenum(1900,1,1) - datenum(2000,1,1) 
        days = double(time)/24. + datenum(1900,1,1) - datenum(2000,1,1);
	ncwrite(frcname,'time',days,[irec]);

	% ---- 10 meter winds -----
	u = get_frc_hrrr(data,grd,'UGRD_10maboveground',irec,'makima');
	v = get_frc_hrrr(data,grd,'VGRD_10maboveground',irec,'makima');
	
	% rotate from HRRR grid angles to ROMS grid angles 
        ugrid = cosa.*u + sina.*v;
        vgrid =-sina.*u + cosa.*v;

	ncwrite(frcname,'uwnd',ugrid,[1 1 irec]);
	ncwrite(frcname,'vwnd',vgrid,[1 1 irec]);

	% ---- Incoming Radiation -----
	swr = get_frc_hrrr(data,grd,'DSWRF_surface',irec,'linear'); % downward_shortwave_flux [W/m2]
	lwr = get_frc_hrrr(data,grd,'DLWRF_surface',irec,'linear'); % downward_longwave_flux  [W/m2]

	if 0
	% ---- Correction to swr using  COREv2 dataset -----
	% This is a multiplicative correcting, trying to account for 
	% for errors in the ERA5 cloud cover
	% temporal interpolation in a monthly climatology
	yr_day = mod(days,365.25);
	if yr_day>=corr_time(12)
           id0 = 12;
           id1 =  1;
	   cf  = (yr_day - corr_time(12))/(corr_time(1)+365.25-corr_time(12));
        elseif yr_day<=corr_time(1)
           id0 = 12;
           id1 =  1;
	   cf  = (yr_day - corr_time(12)+365.25)/(corr_time(1)+365.25-corr_time(12));
        else
           id0 = find(yr_day>corr_time,1,'last');
	   id1 = id0+1;
	   cf  = (yr_day - corr_time(id0))/(corr_time(id1)-corr_time(id0));
	end

	swrd= swr;
	swr_cr = cf*swr_mult(:,:,id0) + (1-cf)*swr_mult(:,:,id1);
	swr= swr.*swr_cr;
	end

	ncwrite(frcname,'swrad',swr,[1 1 irec]);
	ncwrite(frcname,'lwrad',lwr,[1 1 irec]);

	% ---- Absolute humidity -----

	qair = get_frc_hrrr(data,grd,'SPFH_2maboveground',irec,'linear'); % Specific Humidity [kg/kg]
	ncwrite(frcname,'qair',qair,[1 1 irec]);

	t2m = get_frc_hrrr(data,grd,'TMP_2maboveground',irec,'linear'); % 2 meter air temp [K]
	t2m = t2m - 273.15; % convert to degrees C
	ncwrite(frcname,'Tair',t2m, [1 1 irec]);

	% ---- Rain ------
%       rain = get_frc_hrrr(data,grd,'PRATE_surface',irec,'linear'); % precipitation [kg/m^2/s]
	rain = get_frc_hrrr(data,grd,'APCP_surface',irec,'linear'); % 1 hour precipitation accumulation [kg/m^2]

	% Translate to cm/day (1e-3*100*24*3600)
	% Translate to cm/day (1e-3*100*24)
	if 0 % output in m/s
          rain = rain/3600;
        else % output in cm/day
%         rain = rain*0.1*24*3600;
	  rain = rain*0.1*24;
        end

	% ---- Correction to rain using Dai and Trenberth data  -----
	% Temporal interpolation in a monthly climatology
        if add_rivers
          yr_day = mod(days,365.25);
          if yr_day>=riv_time(12)
            id0 = 12;
            id1 =  1;
	    cf  = (yr_day -  riv_time(12))/( riv_time(1)+365.25- riv_time(12));
          elseif yr_day<= riv_time(1)
            id0 = 12;
            id1 =  1;
	   cf  = (yr_day -  riv_time(12)+365.25)/( riv_time(1)+365.25- riv_time(12));
          else
            id0 = find(yr_day>riv_time,1,'last');
	    id1 = id0+1;
	    cf  = (yr_day -  riv_time(id0))/( riv_time(id1)- riv_time(id0));
	  end
          rain = rain + cf*riv_flux(:,:,id0) + (1-cf)*riv_flux(:,:,id1);
        end

	ncwrite(frcname,'rain',rain, [1 1 irec]);

     end
   end % 


