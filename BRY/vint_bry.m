
  % Vertical interpolation for lateral boundary condition values
  % We rely on a strict naming convection for varnames:
  % u_south,v_south,temp_south, etc, etc

  gname = '/paracas/nmolem/GREEN/green_grd.nc';
  fname = '/paracas/nmolem/GREEN/green_bry.20140821002229.nc';
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

    if ndims(f)<3 % surface var or time
      if ndims(f)<2 % time
	      F always has 2 or more dims
        nccreate(oname,vname,'dimensions',{'time',0},'datatype','double');
      else
        [h_bnd,bdimx] = get_hbnd(vname,h);
        nx = size(f,1);
        nccreate(oname,vname,'dimensions',{bdimx,nx,'time',0},'datatype','double');
      end
      ncwrite(oname,vname,f)
      continue
    end
    nx = size(f,1);
    nt = size(f,3);

    % Get appropriate topography at the boundary
%   [h_bnd,bdimx] = get_hbnd(vname,h);
    h_bnd = h_bnd';
    z_bnd = zlevs3(h_bnd,0*h_bnd,ts,tb,hc,nz,'r','new2008')';
    zo_bnd = zlevs3(h_bnd,0*h_bnd,tso,tbo,hco,nzo,'r','new2008')';

    toc

    % Vertical interpolation
    fo = zeros(nx,nzo,nt);

    for i = 1:nx
      fint = interp1(z_bnd(i,:),squeeze(f(i,:,:)),zo_bnd(i,:),'linear','extrap');
      fo(i,:,:) = fint;
    end
    nccreate(oname,vname,'dimensions',{bdimx,nx,bdimz,nzo,'time',0},'datatype','double');
    ncwrite(oname,vname,fo)
    toc
  end
