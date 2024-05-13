%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make a ROMS bulk forcing file using HRRR hourly data
%  Called by: 'make_frc'
%
%  2020-2024, Jeroen Molemaker, UCLA
%
%  Possible future work
% - Add the effect of surface pressure on Qair and other things
% - Maybe force the model with Tair and Dew point instead of Humidity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   frclist = dir([frc_dir frc_source '_Y*']);
   nfiles = length(frclist);

   % HRRR time is in seconds since 1970

   % find the right files
   stime = zeros(nfiles,1);
   for i=1:nfiles
     datname = [frc_dir frclist(i).name];
     stime(i) = double(ncread(datname,'time',[1],[1]))/(24*3600) + datenum(1970,1,1);
   end
   t0 = find(stime<=start_date,1,'last');
   t1 = find(stime>end_date,1,'first');
   if end_date>=stime(end)
    t1 = nfiles
   end
%  t0 = t0+1


   % trim list of filenames
   frclist = frclist(t0:t1);
   nfiles = length(frclist);

   % Read in data coordinates and set data trim.
   disp(' ')
   if coarse_frc
     add_coarse_grd(grdname);

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

   datname = [frc_dir frclist(1).name];
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
     datname = [frc_dir frclist(i).name];
     disp([' Processing: ' frclist(i).name])

     dat_time = ncread(datname,'time');
     nrecord = length(dat_time);

     date_num = double(dat_time(1)/(24*3600)) + datenum(1970,1,1);
     label = datestr(date_num,'YYYYmm')

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

       if mod(irec,10)==0
        disp(['Record: ' num2str(irec)])
       end

       % ---- time -----
       time = ncread(datname,'time',[irec],[1]);


	% translate to days since 2000,1,1
        days = double(time)/24./3600. + datenum(1970,1,1) - datenum(2000,1,1);
	ncwrite(frcname,'time',days,[irec]);
        ncwrite(frcname,'rad_time',days-0.5/24,[irec]);



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


