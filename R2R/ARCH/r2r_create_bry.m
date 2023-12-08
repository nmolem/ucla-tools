function r2r_create_bry(bryname,grdname,obcflag,param,cycle);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   function l2r_create_bry(bryname,grdname,obcflag,...
%                          chdscd,cycle)
%
%   Input:
%
%   bryname      Netcdf climatology file name (character string)
%   grdname      Netcdf grid file name (character string)
%   obcflag      open boundary flag (1=open, [S E N W])
%   chdscd       S-coordinate parameters (object)
%   cycle        Length (days) for cycling the climatology (real)

%
%
% get S-coordinate parameters
%
theta_b = param.theta_b;
theta_s = param.theta_s;
hc      = param.hc;
N       = param.N;
%
%
%  Read the grid file and check the topography
%
nc      = netcdf(grdname, 'nowrite');
h       = nc{'h'}(:);
maskr   = nc{'mask_rho'}(:);
Lp      = length(nc('xi_rho'));
Mp      = length(nc('eta_rho'));
status = close(nc);
hmin   = min(min(h(maskr==1)));

L  = Lp - 1;
M  = Mp - 1;
Np = N  + 1;

%
%  Create the boundary file
%
type    = 'BOUNDARY file';
history = 'ROMS';
nc      = netcdf(bryname,'clobber');
result  = redef(nc);
%
%  Create dimensions
%
time = [1:5];
nc('xi_u')     = L;
nc('xi_rho')   = Lp;
nc('eta_v')    = M;
nc('eta_rho')  = Mp;
nc('s_rho')    = N;
nc('bry_time') = 0;
nc('one')      = 1;
%
  [nxu,ny] = size(u);
  [nx,nyv] = size(v);
  nccreate(bryname,'ocean_time','Dimensions', {'time', inf});
  nccreate(bryname,'u','Dimensions', {'xi_u',nxu,'eta_rho',ny, 'time', inf});
  nccreate(bryname,'v','Dimensions', {'xi_rho',nx,'eta_v',nyv, 'time', inf});
  nccreate(bryname,'temp','Dimensions', {'xi_rho',nx,'eta_rho',ny, 'time', inf});
  nccreate(bryname,'salt','Dimensions', {'xi_rho',nx,'eta_rho',ny, 'time', inf});
  nccreate(bryname,'zeta','Dimensions', {'xi_rho',nx,'eta_rho',ny, 'time', inf});

%  Create variables and attributes
%
nccreate(bryname,'theta_s','Dimensions',{'one',1});
nc{'theta_s'} = ncdouble('one') ;
nc{'theta_s'}.long_name = 'S-coordinate surface control parameter';
nc{'theta_s'}.units = 'nondimensional';
%
nc{'theta_b'} = ncdouble('one') ;
nc{'theta_b'}.long_name = 'S-coordinate bottom control parameter';
nc{'theta_b'}.units = 'nondimensional';
%
nc{'hc'} = ncdouble('one') ;
nc{'hc'}.long_name = 'S-coordinate parameter, critical depth';
nc{'hc'}.units = 'meter';
%
nc{'bry_time'} = ncdouble('bry_time') ;
nc{'bry_time'}.long_name = 'time for boundary data';
nc{'bry_time'}.units = 'day';
%
if cycle>0
 nc{'bry_time'}.cycle_length = cycle;
end
%
if obcflag(1)==1  %%   Southern boundary
%
  nc{'temp_south'} = ncfloat('bry_time','s_rho','xi_rho') ;
  nc{'temp_south'}.long_name = 'southern boundary potential temperature';
  nc{'temp_south'}.units = 'Celsius';
%
  nc{'salt_south'} = ncfloat('bry_time','s_rho','xi_rho') ;
  nc{'salt_south'}.long_name = 'southern boundary salinity';
  nc{'salt_south'}.units = 'PSU';
%
  nc{'u_south'} = ncfloat('bry_time','s_rho','xi_u') ;
  nc{'u_south'}.long_name = 'southern boundary u-momentum component';
  nc{'u_south'}.units = 'meter second-1';
%
  nc{'v_south'} = ncfloat('bry_time','s_rho','xi_rho') ;
  nc{'v_south'}.long_name = 'southern boundary v-momentum component';
  nc{'v_south'}.units = 'meter second-1';
%
  nc{'ubar_south'} = ncfloat('bry_time','xi_u') ;
  nc{'ubar_south'}.long_name = 'southern boundary vertically integrated u-momentum component';
  nc{'ubar_south'}.units = 'meter second-1';
%
  nc{'vbar_south'} = ncfloat('bry_time','xi_rho') ;
  nc{'vbar_south'}.long_name = 'southern boundary vertically integrated v-momentum component';
  nc{'vbar_south'}.units = 'meter second-1';
%
  nc{'zeta_south'} = ncfloat('bry_time','xi_rho') ;
  nc{'zeta_south'}.long_name = 'southern boundary sea surface height';
  nc{'zeta_south'}.units = 'meter';
end
%
%
if obcflag(2)==1  %%   Eastern boundary
%
  nc{'temp_east'} = ncfloat('bry_time','s_rho','eta_rho') ;
  nc{'temp_east'}.long_name = 'eastern boundary potential temperature';
  nc{'temp_east'}.units = 'Celsius';
%
  nc{'salt_east'} = ncfloat('bry_time','s_rho','eta_rho') ;
  nc{'salt_east'}.long_name = 'eastern boundary salinity';
  nc{'salt_east'}.units = 'PSU';
%
  nc{'u_east'} = ncfloat('bry_time','s_rho','eta_rho') ;
  nc{'u_east'}.long_name = 'eastern boundary u-momentum component';
  nc{'u_east'}.units = 'meter second-1';
%
  nc{'v_east'} = ncfloat('bry_time','s_rho','eta_v') ;
  nc{'v_east'}.long_name = 'eastern boundary v-momentum component';
  nc{'v_east'}.units = 'meter second-1';
%
  nc{'ubar_east'} = ncfloat('bry_time','eta_rho') ;
  nc{'ubar_east'}.long_name = 'eastern boundary vertically integrated u-momentum component';
  nc{'ubar_east'}.units = 'meter second-1';
%
  nc{'vbar_east'} = ncfloat('bry_time','eta_v') ;
  nc{'vbar_east'}.long_name = 'eastern boundary vertically integrated v-momentum component';
  nc{'vbar_east'}.units = 'meter second-1';
%
  nc{'zeta_east'} = ncfloat('bry_time','eta_rho') ;
  nc{'zeta_east'}.long_name = 'eastern boundary sea surface height';
  nc{'zeta_east'}.units = 'meter';
end
%
if obcflag(3)==1  %%   Northern boundary
%
  nc{'temp_north'} = ncfloat('bry_time','s_rho','xi_rho') ;
  nc{'temp_north'}.long_name = 'northern boundary potential temperature';
  nc{'temp_north'}.units = 'Celsius';
%
  nc{'salt_north'} = ncfloat('bry_time','s_rho','xi_rho') ;
  nc{'salt_north'}.long_name = 'northern boundary salinity';
  nc{'salt_north'}.units = 'PSU';
%
  nc{'u_north'} = ncfloat('bry_time','s_rho','xi_u') ;
  nc{'u_north'}.long_name = 'northern boundary u-momentum component';
  nc{'u_north'}.units = 'meter second-1';
%
  nc{'v_north'} = ncfloat('bry_time','s_rho','xi_rho') ;
  nc{'v_north'}.long_name = 'northern boundary v-momentum component';
  nc{'v_north'}.units = 'meter second-1';
%
  nc{'ubar_north'} = ncfloat('bry_time','xi_u') ;
  nc{'ubar_north'}.long_name = 'northern boundary vertically integrated u-momentum component';
  nc{'ubar_north'}.units = 'meter second-1';
%
  nc{'vbar_north'} = ncfloat('bry_time','xi_rho') ;
  nc{'vbar_north'}.long_name = 'northern boundary vertically integrated v-momentum component';
  nc{'vbar_north'}.units = 'meter second-1';
%
  nc{'zeta_north'} = ncfloat('bry_time','xi_rho') ;
  nc{'zeta_north'}.long_name = 'northern boundary sea surface height';
  nc{'zeta_north'}.units = 'meter';
end
%
if obcflag(4)==1  %%   Western boundary
%
  nc{'temp_west'} = ncfloat('bry_time','s_rho','eta_rho') ;
  nc{'temp_west'}.long_name = 'western boundary potential temperature';
  nc{'temp_west'}.units = 'Celsius';
%
  nc{'salt_west'} = ncfloat('bry_time','s_rho','eta_rho') ;
  nc{'salt_west'}.long_name = 'western boundary salinity';
  nc{'salt_west'}.units = 'PSU';
%
  nc{'u_west'} = ncfloat('bry_time','s_rho','eta_rho') ;
  nc{'u_west'}.long_name = 'western boundary u-momentum component';
  nc{'u_west'}.units = 'meter second-1';
%
  nc{'v_west'} = ncfloat('bry_time','s_rho','eta_v') ;
  nc{'v_west'}.long_name = 'western boundary v-momentum component';
  nc{'v_west'}.units = 'meter second-1';
%
  nc{'ubar_west'} = ncfloat('bry_time','eta_rho') ;
  nc{'ubar_west'}.long_name = 'western boundary vertically integrated u-momentum component';
  nc{'ubar_west'}.units = 'meter second-1';
%
  nc{'vbar_west'} = ncfloat('bry_time','eta_v') ;
  nc{'vbar_west'}.long_name = 'western boundary vertically integrated v-momentum component';
  nc{'vbar_west'}.units = 'meter second-1';
%
  nc{'zeta_west'} = ncfloat('bry_time','eta_rho') ;
  nc{'zeta_west'}.long_name = 'western boundary sea surface height';
  nc{'zeta_west'}.units = 'meter';
end
%
%
% Create global attributes
%
nc.title     = ['Boundary file for' grdname];
nc.date      = date;
nc.grd_file  = grdname;
nc.type      = type;
nc.history   = history;
%
% Leave define mode
%
result = endef(nc);
%
% Write variables
%
nc{'theta_s'}(:)      = theta_s;
nc{'theta_b'}(:)      = theta_b;
nc{'Tcline'}(:)       = hc;
nc{'hc'}(:)           = hc;
close(nc)
