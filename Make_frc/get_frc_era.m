	function fld = get_frc_era(data,grd,varname,irec,method)
	%

	datname  = data.datname;
	lon_frc = data.lon;
	lat_frc = data.lat;
	i0 = data.i0;
	j0 = data.j0;
	fnx = data.fnx;
	fny = data.fny;
        mask= data.mask;

	lon = grd.lon;
	lat = grd.lat;
	straddle = grd.straddle;

	if straddle
          i1 = data.i1;
          frc = [ ncread(datname,varname,[i0 j0 irec],[inf fny 1])' ...
                  ncread(datname,varname,[ 1 j0 irec],[ i1 fny 1])' ]';
	else
          frc = ncread(datname,varname,[i0 j0 irec],[fnx fny 1]); 
	end
  	frc = fliplr(frc);

        frc(mask<1) = nan;
        frc = inpaint_nans(frc,2);
        fld = interp2(lon_frc,lat_frc,frc',lon,lat,method);

	return
