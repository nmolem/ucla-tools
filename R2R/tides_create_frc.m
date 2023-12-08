function r2r_create_bry(frcfile,grdfile,obcflag,param,cycle);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   function tides_create_frc(frcname,grdname,
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
h       = ncread(grdfile,'h')';
maskr   = ncread(grdfile,'mask_rho')';
[Mp,Lp] = size(h);
%hmin   = min(min(h(maskr==1)));

L  = Lp - 1;
M  = Mp - 1;
Np = N  + 1;

%
%  Create the boundary file
%
%
%  Create dimensions
%c('xi_u')     = L;
%c('xi_rho')   = Lp;
%c('eta_v')    = M;
%c('eta_rho')  = Mp;
%c('s_rho')    = N;
%c('bry_time') = 0;
%c('one')      = 1;
%
%  Create variables and attributes
%
nccreate(frcfile,'Frequencies','dimensions',{'one',1},'datatype','single');
ncwriteatt(frcfile,'theta_s','long_name','S-coordinate surface control parameter');
ncwriteatt(frcfile,'theta_s','units','nondimensional');
%
nccreate(frcfile,'theta_b','dimensions',{'one',1},'datatype','single');
ncwriteatt(frcfile,'theta_b','long_name','S-coordinate bottom control parameter');
ncwriteatt(frcfile,'theta_b','units','nondimensional');
%
nccreate(frcfile,'hc','dimensions',{'one',1},'datatype','single');
ncwriteatt(frcfile,'hc','long_name','S-coordinate parameter critical depth');
ncwriteatt(frcfile,'hc','units','meter');
%
nccreate(frcfile,'Tcline','dimensions',{'one',1},'datatype','single');
ncwriteatt(frcfile,'Tcline','long_name','S-coordinate surface/bottom layer width');
ncwriteatt(frcfile,'Tcline','units','meter');
%
nccreate(frcfile,'bry_time','dimensions',{'bry_time',0},'datatype','single');
ncwriteatt(frcfile,'bry_time','long_name','Time for boundary data');
ncwriteatt(frcfile,'bry_time','units','day');
%
if cycle>0
  ncwriteatt(frcfile,'bry_time','cycle_length',cycle);
end
%
if obcflag(1)==1  %%   Southern boundary
%
  nccreate(frcfile,'temp_south','dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'temp_south','long_name','Southern boundary potential temperature');
  ncwriteatt(frcfile,'temp_south','units','Celsius');
%
  nccreate(frcfile,'salt_south','dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'salt_south','long_name','Southern boundary salinity');
  ncwriteatt(frcfile,'salt_south','units','PSU');
%
  nccreate(frcfile,'u_south','dimensions',{'xi_u',L,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'u_south','long_name','Southern u-momentum component');
  ncwriteatt(frcfile,'u_south','units','meter second -1');
%
  nccreate(frcfile,'v_south','dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'v_south','long_name','Southern v-momentum component');
  ncwriteatt(frcfile,'v_south','units','meter second -1');
%
  nccreate(frcfile,'w_south','dimensions',{'xi_rho',Lp,'s_w',Np,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'w_south','long_name','Southern w-flux component');
  ncwriteatt(frcfile,'w_south','units','meter second -1');
%
  nccreate(frcfile,'ubar_south','dimensions',{'xi_u',L,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'ubar_south','long_name','Southern integrated u-momentum component');
  ncwriteatt(frcfile,'ubar_south','units','meter second -1');
%
  nccreate(frcfile,'vbar_south','dimensions',{'xi_rho',Lp,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'vbar_south','long_name','Southern integrated v-momentum component');
  ncwriteatt(frcfile,'vbar_south','units','meter second -1');
%
  nccreate(frcfile,'zeta_south','dimensions',{'xi_rho',Lp,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'zeta_south','long_name','Southern boundary sea surface height');
  ncwriteatt(frcfile,'zeta_south','units','meter');
end
%
%
if obcflag(2)==1  %%   Eastern boundary
%
  nccreate(frcfile,'temp_east','dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'temp_east','long_name','Eastern boundary potential temperature');
  ncwriteatt(frcfile,'temp_east','units','Celsius');
%
  nccreate(frcfile,'salt_east','dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'salt_east','long_name','Eastern boundary salinity');
  ncwriteatt(frcfile,'salt_east','units','PSU');
%
  nccreate(frcfile,'u_east','dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'u_east','long_name','Eastern u-momentum component');
  ncwriteatt(frcfile,'u_east','units','meter second -1');
%
  nccreate(frcfile,'v_east','dimensions',{'eta_v',M,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'v_east','long_name','Eastern v-momentum component');
  ncwriteatt(frcfile,'v_east','units','meter second -1');
%
  nccreate(frcfile,'w_east','dimensions',{'eta_rho',Mp,'s_w',Np,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'w_east','long_name','Eastern w-flux component');
  ncwriteatt(frcfile,'w_east','units','meter second -1');
%
  nccreate(frcfile,'ubar_east','dimensions',{'eta_rho',Mp,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'ubar_east','long_name','Eastern integrated u-momentum component');
  ncwriteatt(frcfile,'ubar_east','units','meter second -1');
%
  nccreate(frcfile,'vbar_east','dimensions',{'eta_v',M,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'vbar_east','long_name','Eastern integrated v-momentum component');
  ncwriteatt(frcfile,'vbar_east','units','meter second -1');
%
  nccreate(frcfile,'zeta_east','dimensions',{'eta_rho',Mp,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'zeta_east','long_name','Eastern boundary sea surface height');
  ncwriteatt(frcfile,'zeta_east','units','meter');
end
%
if obcflag(3)==1  %%   Northern boundary
  nccreate(frcfile,'temp_north','dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'temp_north','long_name','Northern boundary potential temperature');
  ncwriteatt(frcfile,'temp_north','units','Celsius');
%
  nccreate(frcfile,'salt_north','dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'salt_north','long_name','Northern boundary salinity');
  ncwriteatt(frcfile,'salt_north','units','PSU');
%
  nccreate(frcfile,'u_north','dimensions',{'xi_u',L,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'u_north','long_name','Northern u-momentum component');
  ncwriteatt(frcfile,'u_north','units','meter second -1');
%
  nccreate(frcfile,'v_north','dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'v_north','long_name','Northern v-momentum component');
  ncwriteatt(frcfile,'v_north','units','meter second -1');
%
  nccreate(frcfile,'w_north','dimensions',{'xi_rho',Lp,'s_w',Np,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'w_north','long_name','Northern w-flux component');
  ncwriteatt(frcfile,'w_north','units','meter second -1');
%
  nccreate(frcfile,'ubar_north','dimensions',{'xi_u',L,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'ubar_north','long_name','Northern integrated u-momentum component');
  ncwriteatt(frcfile,'ubar_north','units','meter second -1');
%
  nccreate(frcfile,'vbar_north','dimensions',{'xi_rho',Lp,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'vbar_north','long_name','Northern integrated v-momentum component');
  ncwriteatt(frcfile,'vbar_north','units','meter second -1');
%
  nccreate(frcfile,'zeta_north','dimensions',{'xi_rho',Lp,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'zeta_north','long_name','Northern boundary sea surface height');
  ncwriteatt(frcfile,'zeta_north','units','meter');
end
%
if obcflag(4)==1  %%   Western boundary
  nccreate(frcfile,'temp_west','dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'temp_west','long_name','Western boundary potential temperature');
  ncwriteatt(frcfile,'temp_west','units','Celsius');
%
  nccreate(frcfile,'salt_west','dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'salt_west','long_name','Western boundary salinity');
  ncwriteatt(frcfile,'salt_west','units','PSU');
%
  nccreate(frcfile,'u_west','dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'u_west','long_name','Western u-momentum component');
  ncwriteatt(frcfile,'u_west','units','meter second -1');
%
  nccreate(frcfile,'v_west','dimensions',{'eta_v',M,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'v_west','long_name','Western v-momentum component');
  ncwriteatt(frcfile,'v_west','units','meter second -1');
%
  nccreate(frcfile,'w_west','dimensions',{'eta_rho',Mp,'s_w',Np,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'w_west','long_name','Western w-flux component');
  ncwriteatt(frcfile,'w_west','units','meter second -1');
%
  nccreate(frcfile,'ubar_west','dimensions',{'eta_rho',Mp,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'ubar_west','long_name','Western integrated u-momentum component');
  ncwriteatt(frcfile,'ubar_west','units','meter second -1');
%
  nccreate(frcfile,'vbar_west','dimensions',{'eta_v',M,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'vbar_west','long_name','Western integrated v-momentum component');
  ncwriteatt(frcfile,'vbar_west','units','meter second -1');
%
  nccreate(frcfile,'zeta_west','dimensions',{'eta_rho',Mp,'bry_time',0},'datatype','single');
  ncwriteatt(frcfile,'zeta_west','long_name','Western boundary sea surface height');
  ncwriteatt(frcfile,'zeta_west','units','meter');
end
%
%
%  Write global attributes
%
ncwriteatt(frcfile,'/','Title',['Boundary file for' grdfile]);
ncwriteatt(frcfile,'/','Date',date);
%
% Write variables
%
ncwrite(frcfile,'theta_s',theta_s);
ncwrite(frcfile,'theta_b',theta_b);
ncwrite(frcfile,'Tcline',hc);
ncwrite(frcfile,'hc',hc);

