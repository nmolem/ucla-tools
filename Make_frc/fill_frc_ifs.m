
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make a ROMS bulk forcing file using ECNWF hourly data
%  Called by: 'make_frc'
%
%  2020-2024, Jeroen Molemaker, Pierre Damien, UCLA
%
%  Possible future work
% - Add the effect of surface pressure on Qair and other things
% - Maybe force the model with Tair and Dew point instead of Humidity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   frclist = dir([frc_dir frc_source '_Y*']);
   nfiles = length(frclist);

   % find the right files
   stime = zeros(nfiles,1);
   for i=1:nfiles
     datname = [frc_dir frclist(i).name];
     stime(i) = double(ncread(datname,'time',[1],[1]))/24 + datenum(1900,1,1);
   end
   t0 = find(stime<=start_date,1,'last');
   t1 = find(stime>end_date,1,'first');

   if stime(t1)-1./24 >= end_date
    t1 = t1-1;
   end

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

   lon0 = min(lon(:));
   lon1 = max(lon(:));
   lat0 = min(lat(:));
   lat1 = max(lat(:));

   % Figure out if the grid is straddling the dateline.

   % Assume that straddles are in the [-180,180] range
   if lon0<0 & lon1>0
     grd.straddle = 1;
   else
     grd.straddle = 0;
     lon = mod(lon,360); %% Atmospheric data is in the 0 to 360 longitude range
   end

   cosa = cos(ang);
   sina = sin(ang);

   grd.lon = lon;
   grd.lat = lat;

   disp(' ')
   disp(' Read in the data grid')

   datname = [frc_dir frclist(1).name];
   lon_frc = ncread(datname,'longitude');
   lat_frc = ncread(datname,'latitude');

   if grd.straddle
     % lon0 is negative, lon1 positive
     i0 = find(lon_frc-360<lon0,1,'last');
     i1 = find(lon_frc    >lon1,1,'first');
     lon_frc = [lon_frc(i0:end)'-360 lon_frc(1:i1)']';
   else
     i0 = find(lon_frc<lon0,1,'last');
     i1 = find(lon_frc>lon1,1,'first');
     lon_frc = lon_frc(i0:i1);
   end

   % Atmospheric data is written in upside down order, latitude is decreasing
   j0 = find(lat_frc>lat1,1,'last');
   j1 = find(lat_frc<lat0,1,'first');
   lat_frc = flipud(lat_frc(j0:j1));

   fnx = length(lon_frc);
   fny = length(lat_frc);

   if grd.straddle
     mask = [ ncread(maskname,'mask',[i0 j0],[inf fny])' ...
              ncread(maskname,'mask',[ 1 j0],[ i1 fny])' ]';
   else
     mask =   ncread(maskname,'mask',[i0 j0],[fnx fny]);
   end
   mask = fliplr(mask); % To deal with upside down data

   if wind_dropoff
     if ~exist(disname)
       cdist = comp_cdist(grdname,disname,coarse_frc);
     else
       load(disname)
     end
     cdist = cdist/1e3;
     mult = 1-0.4*exp(-0.08*cdist);
   end

   % Prepare for short and long wave radiation correction
   if rad_corr
     corr_time = ncread(rcorname,'time'); % time in year days
     swr_mult = zeros(nx,ny,12);
     lwr_mult = zeros(nx,ny,12);
     for i=1:12
       if grd.straddle
         corr_swr = [ ncread(rcorname,'ssr_corr',[i0 j0 i],[inf fny 1])' ...
                      ncread(rcorname,'ssr_corr',[ 1 j0 i],[ i1 fny 1])' ]';
       else
         corr_swr =   ncread(rcorname,'ssr_corr',[i0 j0 i],[fnx fny 1]);
       end
       corr_swr = fliplr(corr_swr);
       corr_swr(mask<1) = nan;
       corr_swr = inpaint_nans(corr_swr,2);
       swr_mult(:,:,i) = interp2(lon_frc,lat_frc,corr_swr',lon,lat,'makima');

       if grd.straddle
         corr_lwr = [ ncread(rcorname,'strd_corr',[i0 j0 i],[inf fny 1])' ...
                      ncread(rcorname,'strd_corr',[ 1 j0 i],[ i1 fny 1])' ]';
       else
         corr_lwr =   ncread(rcorname,'strd_corr',[i0 j0 i],[fnx fny 1]);
       end
       corr_lwr = fliplr(corr_lwr);
       corr_lwr(mask<1) = nan;
       corr_lwr = inpaint_nans(corr_lwr,2);
       lwr_mult(:,:,i) = interp2(lon_frc,lat_frc,corr_lwr',lon,lat,'makima');
     end
   end % rad_corr

   if add_rivers
     % prepare for the addition of river runoff to rain data
     riv_time = ncread(rivname,'swf_time');
     riv_flux = ncread(rivname,'swf_flux');
   end

   % Loop over data files
   for i = 1:nfiles
     datname = [frc_dir frclist(i).name];
     disp([' Processing: ' frclist(i).name])

     dat_time = ncread(datname,'time');
     nrecord = length(dat_time);

     date_num = double(dat_time(1)/24) + datenum(1900,1,1);
     label = datestr(date_num,'YYYYmm')

     frcname = [root_name '_frc.' label '.nc'] 

     data.datname = datname;
     data.lon = lon_frc;
     data.lat = lat_frc;
     data.i0 = i0;
     data.i1 = i1;
     data.j0 = j0;
     data.fnx = fnx;
     data.fny = fny;
     data.mask = mask;

     if exist(frcname)
       delete(frcname)
     end
     create_frc_bulk(grdname,frcname,coarse_frc);
     ncwriteatt(frcname,'/','Data Source',dsatt);


     for irec = 1:nrecord
       disp(['Record: ' num2str(irec)])

       % ---- time -----
       time = ncread(datname,'time',[irec],[1]);

       % time for the fluxes (radiation and rain) is centered between
       % the hourly times of the outputs. This is because the raw output
       % of the IFS model is the integrated flux between the current and the
       % previous time (also see the division by 3600 below)

       % translate to days since 2000,1,1
       days = double(time)/24. + datenum(1900,1,1) - datenum(2000,1,1);
       ncwrite(frcname,'time',days,[irec]);
       ncwrite(frcname,'rad_time',days-0.5/24,[irec]);

       % ---- 10 meter winds -----
       u = get_frc_era(data,grd,'u10',irec,'makima');
       v = get_frc_era(data,grd,'v10',irec,'makima');
	
       % rotate to grid angles 
       ugrid = cosa.*u + sina.*v;
       vgrid =-sina.*u + cosa.*v;

       if wind_dropoff
         ugrid = ugrid.*mult;
         vgrid = vgrid.*mult;
       end

       ncwrite(frcname,'uwnd',ugrid,[1 1 irec]);
       ncwrite(frcname,'vwnd',vgrid,[1 1 irec]);

       % ---- Incoming Radiation -----
       swr = get_frc_era(data,grd,'ssr',irec,'linear');  % downward_shortwave_flux [J/m2]
       lwr = get_frc_era(data,grd,'strd',irec,'linear'); % downward_longwave_flux [J/m2]

	% Translate to fluxes from values that are integrated over 1 hour
	% Governing times are offset by 30 minutes (rad_time)
	swr = swr/3600;
	lwr = lwr/3600;

	if rad_corr
          % Multiplicative correction to net shortwave and incoming
	  % long wave radiation. This correction is trying to account for 
	  % for errors in the model's cloud coveing a monthly satellite
	  % climatology

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
          swr_cr = cf*swr_mult(:,:,id0) + (1-cf)*swr_mult(:,:,id1);
          swr= swr.*swr_cr;
          lwr_cr = cf*lwr_mult(:,:,id0) + (1-cf)*lwr_mult(:,:,id1);
          lwr= lwr.*lwr_cr;
        end % rad_corr

	ncwrite(frcname,'swrad',swr,[1 1 irec]);
	ncwrite(frcname,'lwrad',lwr,[1 1 irec]);

	% ---- Absolute humidity -----

	t2m = get_frc_era(data,grd,'t2m',irec,'linear'); % 2 meter air temp [K]
	d2m = get_frc_era(data,grd,'d2m',irec,'linear'); % 2 meter air dew point temp [K]
	t2m = t2m - 273.15; % K to C
	d2m = d2m - 273.15; % K to C
        Qair=(exp((17.625*d2m)./(243.04+d2m))./exp((17.625*t2m)./(243.04+t2m)));  % Relative humidity fraction

        % Relative to absolute humidity assuming constant pressure

        patm=1010.0; 

        cff=(1.0007+3.46e-6*patm).*6.1121.*exp(17.502*t2m./(240.97+t2m));
        cff=cff.* Qair;
        qair =0.62197.*(cff./(patm-0.378.*cff)) ;

	ncwrite(frcname,'qair',qair,[1 1 irec]);
	ncwrite(frcname,'Tair',t2m, [1 1 irec]);

	% ---- Rain ------
	rain = get_frc_era(data,grd,'tp',irec,'linear'); % precipitation [m]

	% Translate to fluxes from values integrated over 1 hour
	if 0 % output in m/s
          rain = rain/3600;
        else % output in cm/day
	  rain = rain*100*24;
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


