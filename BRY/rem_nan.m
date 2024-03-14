
  % Use the outer boundary data from a history file
  % to produce boundary forcing file
  %
  % We rely on a strict naming convection for varnames:
  % u_south,v_south,temp_south, etc, etc

  fname = '/zulu/nmolem/SFB/sfb_rst.nc';

  inf = ncinfo(fname);

  nvar = length(inf.Variables)

  for i = 1:nvar
    vname = inf.Variables(i).Name;
    f = ncread(fname,vname);
    nf = sum(isnan(f(:)));
    if nf>0
      disp(['fixing ' vname])
    end
    f(isnan(f)) = 0;
    ncwrite(fname,vname,f)
  end
