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
%  Create the boundary file with global attributes
%
 ncwriteatt(bryname,'/','Title',['Boundary forcing for' grdname]);
 ncwriteatt(bryname,'/','Date',date);
%
%  Create variables and attributes
%
nccreate(bryname,'theta_s','Dimensions',{'one',1});
ncwriteatt(bryname,'theta_s','Long name','S-coordinate surface control parameter');
ncwriteatt(bryname,'theta_s','Units','Nondimensional');
%
nccreate(bryname,'theta_b','Dimensions',{'one',1});
ncwriteatt(bryname,'theta_b','Long name','S-coordinate bottom control parameter');
ncwriteatt(bryname,'theta_b','Units','Nondimensional');
%
nccreate(bryname,'hc','Dimensions',{'one',1});
ncwriteatt(bryname,'hc','Long name','S-coordinate critical depth control parameter');
ncwriteatt(bryname,'hc','Units','Meter');
%
nccreate(bryname,'bry_time','Dimensions',{'one',1});
ncwriteatt(bryname,'bry_time','Long name','Since January 1, 00:00 am');
ncwriteatt(bryname,'bry_time','Units','Day');
if cycle>0
 nc{'bry_time'}.cycle_length = cycle;
 ncwriteatt(bryname,'bry_time','cycle_length',cycle);
end
%
if obcflag(1)==1  %%   Southern boundary
%
  nccreate(bryname,'temp_south','Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'temp_south','Long name','Southern boundary potential temperature');
  ncwriteatt(bryname,'temp_south','Units','Celsius');
%
  nccreate(bryname,'salt_south','Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'salt_south','Long name','Southern boundary salinity');
  ncwriteatt(bryname,'salt_south','Units','PSU');
%
  nccreate(bryname,'u_south','Dimensions',{'xi_u',L,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'u_south','Long name','Southern boundary u-momentum');
  ncwriteatt(bryname,'u_south','Units','meter/second');
%
  nccreate(bryname,'v_south','Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'v_south','Long name','Southern boundary v-momentum');
  ncwriteatt(bryname,'v_south','Units','meter/second');
%
  nccreate(bryname,'ubar_south','Dimensions',{'xi_u',L,'bry_time', inf});
  ncwriteatt(bryname,'ubar_south','Long name','Southern boundary vertical mean u-momentum');
  ncwriteatt(bryname,'ubar_south','Units','meter/second');
%
  nccreate(bryname,'vbar_south','Dimensions',{'xi_rho',Lp,'bry_time', inf});
  ncwriteatt(bryname,'vbar_south','Long name','Southern boundary vertical mean v-momentum');
  ncwriteatt(bryname,'vbar_south','Units','meter/second');
%
  nccreate(bryname,'zeta_south','Dimensions',{'xi_rho',Lp,'bry_time', inf});
  ncwriteatt(bryname,'zeta_south','Long name','Southern boundary sea surface height');
  ncwriteatt(bryname,'zeta_south','Units','meter');
end
%
%
if obcflag(2)==1  %%   Eastern boundary
%
  nccreate(bryname,'temp_east','Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'temp_east','Long name','Eastern boundary potential temperature');
  ncwriteatt(bryname,'temp_east','Units','Celsius');
%
  nccreate(bryname,'salt_east','Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'salt_east','Long name','Eastern boundary salinity');
  ncwriteatt(bryname,'salt_east','Units','PSU');
%
  nccreate(bryname,'u_east','Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'u_east','Long name','Eastern boundary u-momentum');
  ncwriteatt(bryname,'u_east','Units','meter/second');
%
  nccreate(bryname,'v_east','Dimensions',{'eta_v',M,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'v_east','Long name','Eastern boundary v-momentum');
  ncwriteatt(bryname,'v_east','Units','meter/second');
%
  nccreate(bryname,'ubar_east','Dimensions',{'eta_rho',Mp,'bry_time', inf});
  ncwriteatt(bryname,'ubar_east','Long name','Eastern boundary vertical mean u-momentum');
  ncwriteatt(bryname,'ubar_east','Units','meter/second');
%
  nccreate(bryname,'vbar_east','Dimensions',{'eta_v',M,'bry_time', inf});
  ncwriteatt(bryname,'vbar_east','Long name','Eastern boundary vertical mean v-momentum');
  ncwriteatt(bryname,'vbar_east','Units','meter/second');
%
  nccreate(bryname,'zeta_east','Dimensions',{'eta_rho',Mp,'bry_time', inf});
  ncwriteatt(bryname,'zeta_east','Long name','Eastern boundary sea surface height');
  ncwriteatt(bryname,'zeta_east','Units','meter');
end
%
if obcflag(3)==1  %%   Northern boundary
  nccreate(bryname,'temp_north','Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'temp_north','Long name','Northern boundary potential temperature');
  ncwriteatt(bryname,'temp_north','Units','Celsius');
%
  nccreate(bryname,'salt_north','Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'salt_north','Long name','Northern boundary salinity');
  ncwriteatt(bryname,'salt_north','Units','PSU');
%
  nccreate(bryname,'u_north','Dimensions',{'xi_u',L,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'u_north','Long name','Northern boundary u-momentum');
  ncwriteatt(bryname,'u_north','Units','meter/second');
%
  nccreate(bryname,'v_north','Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'v_north','Long name','Northern boundary v-momentum');
  ncwriteatt(bryname,'v_north','Units','meter/second');
%
  nccreate(bryname,'ubar_north','Dimensions',{'xi_u',L,'bry_time', inf});
  ncwriteatt(bryname,'ubar_north','Long name','Northern boundary vertical mean u-momentum');
  ncwriteatt(bryname,'ubar_north','Units','meter/second');
%
  nccreate(bryname,'vbar_north','Dimensions',{'xi_rho',Lp,'bry_time', inf});
  ncwriteatt(bryname,'vbar_north','Long name','Northern boundary vertical mean v-momentum');
  ncwriteatt(bryname,'vbar_north','Units','meter/second');
%
  nccreate(bryname,'zeta_north','Dimensions',{'xi_rho',Lp,'bry_time', inf});
  ncwriteatt(bryname,'zeta_north','Long name','Northern boundary sea surface height');
  ncwriteatt(bryname,'zeta_north','Units','meter');
end
%
if obcflag(4)==1  %%   Western boundary
%
  nccreate(bryname,'temp_west','Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'temp_west','Long name','Western boundary potential temperature');
  ncwriteatt(bryname,'temp_west','Units','Celsius');
%
  nccreate(bryname,'salt_west','Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'salt_west','Long name','Western boundary salinity');
  ncwriteatt(bryname,'salt_west','Units','PSU');
%
  nccreate(bryname,'u_west','Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'u_west','Long name','Western boundary u-momentum');
  ncwriteatt(bryname,'u_west','Units','meter/second');
%
  nccreate(bryname,'v_west','Dimensions',{'eta_v',M,'s_rho',N,'bry_time', inf});
  ncwriteatt(bryname,'v_west','Long name','Western boundary v-momentum');
  ncwriteatt(bryname,'v_west','Units','meter/second');
%
  nccreate(bryname,'ubar_west','Dimensions',{'eta_rho',Mp,'bry_time', inf});
  ncwriteatt(bryname,'ubar_west','Long name','Western boundary vertical mean u-momentum');
  ncwriteatt(bryname,'ubar_west','Units','meter/second');
%
  nccreate(bryname,'vbar_west','Dimensions',{'eta_v',M,'bry_time', inf});
  ncwriteatt(bryname,'vbar_west','Long name','Western boundary vertical mean v-momentum');
  ncwriteatt(bryname,'vbar_west','Units','meter/second');
%
  nccreate(bryname,'zeta_west','Dimensions',{'eta_rho',Mp,'bry_time', inf});
  ncwriteatt(bryname,'zeta_west','Long name','Western boundary sea surface height');
  ncwriteatt(bryname,'zeta_west','Units','meter');
end
%
% Write variables
%
ncwrite('theta_s'}(:)      = theta_s;
ncwrite('theta_b'}(:)      = theta_b;
ncwrite('hc'}(:)           = hc;
close(nc)
