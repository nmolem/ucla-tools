
  % Use the outer boundary data from a history file
  % to produce boundary forcing file
  %
  % We rely on a strict naming convection for varnames:
  % u_south,v_south,temp_south, etc, etc

% gname = '/zulu/nmolem/GWAKE/gwake_grd.nc';
  fname = '/zulu/nmolem/SFB/sfb_his.nc';

  oname = [fname(1:end-7) '_bry.nc'];

  bdimz = 's_rho'; % add something when dealing with 's_w' vars.

  inf = ncinfo(fname);

  nvar = length(inf.Variables)

  for i = 11:nvar-2
    tic
    vname = inf.Variables(i).Name
    vsz = inf.Variables(i).Size;
    rank = length(inf.Variables(i).Size);

    lname  =  ncreadatt(fname,vname,'long_name');
    units  =  ncreadatt(fname,vname,'units');


    if rank==1
      % time
      f = ncread(fname,vname);
      ovname = 'bry_time'
      nccreate(oname,ovname,'dimensions',{'time',0},'datatype','double');
      ncwrite(oname,ovname,f/(3600*24))
      ncwriteatt(oname,ovname,'long_name',lname)
      ncwriteatt(oname,ovname,'units','days')
      continue
    elseif rank==3 % surface data
      vdm1 = inf.Variables(i).Dimensions(1).Name;
      vdm2 = inf.Variables(i).Dimensions(2).Name;
      % south
        ovname = [vname '_south'];
        nccreate(oname,ovname,'dimensions',{vdm1, vsz(1),'time',0,},'datatype','double');
        ncwriteatt(oname,ovname,'long_name',lname)
        ncwriteatt(oname,ovname,'units','units')
	f = squeeze(ncread(fname,vname,[1 1 1],[vsz(1) 1 vsz(3)]));
	f(isnan(f)) = 0;
        ncwrite(oname,ovname,f)
      % north
        ovname = [vname '_north'];
        nccreate(oname,ovname,'dimensions',{vdm1, vsz(1),'time',0,},'datatype','double');
        ncwriteatt(oname,ovname,'long_name',lname)
        ncwriteatt(oname,ovname,'units','units')
	f = squeeze(ncread(fname,vname,[1 vsz(2) 1],[vsz(1) 1 vsz(3)]));
	f(isnan(f)) = 0;
        ncwrite(oname,ovname,f)
      % west
        ovname = [vname '_west'];
        nccreate(oname,ovname,'dimensions',{vdm2, vsz(2),'time',0,},'datatype','double');
        ncwriteatt(oname,ovname,'long_name',lname)
        ncwriteatt(oname,ovname,'units','units')
	f = squeeze(ncread(fname,vname,[1 1 1],[1 vsz(2) vsz(3)]));
	f(isnan(f)) = 0;
        ncwrite(oname,ovname,f)
      % east
        ovname = [vname '_east'];
        nccreate(oname,ovname,'dimensions',{vdm2, vsz(2),'time',0,},'datatype','double');
        ncwriteatt(oname,ovname,'long_name',lname)
        ncwriteatt(oname,ovname,'units','units')
	f = squeeze(ncread(fname,vname,[vsz(1) 1 1],[1 vsz(2) vsz(3)]));
	f(isnan(f)) = 0;
        ncwrite(oname,ovname,f)
    elseif rank==4 % volume data
      vdm1 = inf.Variables(i).Dimensions(1).Name;
      vdm2 = inf.Variables(i).Dimensions(2).Name;
      vdm3 = inf.Variables(i).Dimensions(3).Name;
      % south
        ovname = [vname '_south'];
        nccreate(oname,ovname,'dimensions',{vdm1, vsz(1),vdm3, vsz(3),'time',0,},'datatype','double');
        ncwriteatt(oname,ovname,'long_name',lname)
        ncwriteatt(oname,ovname,'units','units')
	f = squeeze(ncread(fname,vname,[1 1 1 1],[vsz(1) 1 vsz(3) vsz(4)]));
	f(isnan(f)) = 0;
        ncwrite(oname,ovname,f)
      % north
        ovname = [vname '_north'];
        nccreate(oname,ovname,'dimensions',{vdm1, vsz(1),vdm3, vsz(3),'time',0,},'datatype','double');
        ncwriteatt(oname,ovname,'long_name',lname)
        ncwriteatt(oname,ovname,'units','units')
	f = squeeze(ncread(fname,vname,[1 vsz(2) 1 1],[vsz(1) 1 vsz(3) vsz(4)]));
	f(isnan(f)) = 0;
        ncwrite(oname,ovname,f)
      % west
        ovname = [vname '_west'];
        nccreate(oname,ovname,'dimensions',{vdm2, vsz(2),vdm3, vsz(3),'time',0,},'datatype','double');
        ncwriteatt(oname,ovname,'long_name',lname)
        ncwriteatt(oname,ovname,'units','units')
	f = squeeze(ncread(fname,vname,[1 1 1 1],[1 vsz(2) vsz(3) vsz(4)]));
	f(isnan(f)) = 0;
        ncwrite(oname,ovname,f)
      % east
        ovname = [vname '_east'];
        nccreate(oname,ovname,'dimensions',{vdm2, vsz(2),vdm3, vsz(3),'time',0,},'datatype','double');
        ncwriteatt(oname,ovname,'long_name',lname)
        ncwriteatt(oname,ovname,'units','units')
	f = squeeze(ncread(fname,vname,[vsz(1) 1 1 1],[1 vsz(2) vsz(3) vsz(4)]));
	f(isnan(f)) = 0;
        ncwrite(oname,ovname,f)
    end
  end

