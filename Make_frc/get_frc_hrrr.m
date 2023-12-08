	function fld = get_frc_hrrr(data,grd,varname,irec,method)

	datname  = data.datname;
	x_frc = data.x;
	y_frc = data.y;
	i0 = data.i0;
	j0 = data.j0;
	fnx = data.fnx;
	fny = data.fny;
        mask= data.mask;

	x = grd.x;
	y = grd.y;

        frc = ncread(datname,varname,[i0 j0 irec],[fnx fny 1]); 

        frc(mask<1) = nan;
        frc = inpaint_nans(frc,2);
        fld = interp2(x_frc,y_frc,frc',x,y,method,9999);

	return
