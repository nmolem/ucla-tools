
  % Vertical interpolation for lateral boundary condition values
  % We rely on a strict naming convection for varnames:
  % u_south,v_south,temp_south, etc, etc

  gname = '/zulu/nmolem/GWAKE/gwake_grd.nc';
  fname = '/zulu/nmolem/GWAKE/gwake_bry.201408.nc';

  oname = [fname(1:end-3) '_int.nc'];

  ts= 6;
  tb= 6;
  hc= 250;
  nz= 100;

  tso= 6;
  tbo= 6;
  hco= 250;
  nzo= 200;

  h = ncread(gname,'h');

  bdimz = 's_rho'; % add something when dealing with 's_w' vars.

  inf = ncinfo(fname);

  nvar = length(inf.Variables)

  for i = 1:nvar
    tic
    vname = inf.Variables(i).Name
    f = ncread(fname,vname);
    lname  =  ncreadatt(fname,vname,'long_name');
    units  =  ncreadatt(fname,vname,'units');


    if ndims(f)<3 % surface var or time
      if size(f,2)==1 % time
        nccreate(oname,vname,'dimensions',{'time',0},'datatype','double');
      else
        [h_bnd,bdimx] = get_hbnd(vname,h);
        nx = size(f,1);
        nccreate(oname,vname,'dimensions',{bdimx,nx,'time',0},'datatype','double');
      end
      ncwrite(oname,vname,f)
      ncwriteatt(oname,vname,'long_name',lname)
      ncwriteatt(oname,vname,'units',units)
      continue
    end
    nx = size(f,1);
    nt = size(f,3);

    % Get appropriate topography at the boundary
    [h_bnd,bdimx] = get_hbnd(vname,h);
    z_bnd  = zlevs3(h_bnd,0*h_bnd,ts,tb,hc,nz,'r','new2008');
    zo_bnd = zlevs3(h_bnd,0*h_bnd,tso,tbo,hco,nzo,'r','new2008');
    z_bnd  = squeeze(z_bnd)';
    zo_bnd = squeeze(zo_bnd)';

    nccreate(oname,vname,'dimensions',{bdimx,nx,bdimz,nzo,'time',0},'datatype','double');
    ncwriteatt(oname,vname,'long_name',lname)
    ncwriteatt(oname,vname,'units',units)

    toc

    % Vertical interpolation
    fo = zeros(nx,nzo,nt);

    for i = 1:nx
      fint = interp1(z_bnd(i,:),squeeze(f(i,:,:)),zo_bnd(i,:),'linear','extrap');
      fo(i,:,:) = fint;
    end
    ncwrite(oname,vname,fo)
    toc
  end
