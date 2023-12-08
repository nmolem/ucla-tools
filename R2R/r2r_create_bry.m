function r2r_create_bry(bryfile,grdfile,obcflag,param,cycle);
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
%
%  Create variables and attributes
%
nccreate(bryfile,'theta_s','dimensions',{'one',1},'datatype','single');
ncwriteatt(bryfile,'theta_s','long_name','S-coordinate surface control parameter');
ncwriteatt(bryfile,'theta_s','units','nondimensional');
%
nccreate(bryfile,'theta_b','dimensions',{'one',1},'datatype','single');
ncwriteatt(bryfile,'theta_b','long_name','S-coordinate bottom control parameter');
ncwriteatt(bryfile,'theta_b','units','nondimensional');
%
nccreate(bryfile,'hc','dimensions',{'one',1},'datatype','single');
ncwriteatt(bryfile,'hc','long_name','S-coordinate parameter critical depth');
ncwriteatt(bryfile,'hc','units','meter');
%
nccreate(bryfile,'Tcline','dimensions',{'one',1},'datatype','single');
ncwriteatt(bryfile,'Tcline','long_name','S-coordinate surface/bottom layer width');
ncwriteatt(bryfile,'Tcline','units','meter');
%
nccreate(bryfile,'bry_time','dimensions',{'bry_time',0},'datatype','single');
ncwriteatt(bryfile,'bry_time','long_name','Time for boundary data');
ncwriteatt(bryfile,'bry_time','units','day');
%
if cycle>0
  ncwriteatt(bryfile,'bry_time','cycle_length',cycle);
end
%
if obcflag(1)==1  %%   Southern boundary
%
  nccreate(bryfile,'temp_south','dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'temp_south','long_name','Southern boundary potential temperature');
  ncwriteatt(bryfile,'temp_south','units','Celsius');
%
  nccreate(bryfile,'salt_south','dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'salt_south','long_name','Southern boundary salinity');
  ncwriteatt(bryfile,'salt_south','units','PSU');
%
  nccreate(bryfile,'u_south','dimensions',{'xi_u',L,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'u_south','long_name','Southern u-momentum component');
  ncwriteatt(bryfile,'u_south','units','meter second -1');
%
  nccreate(bryfile,'v_south','dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'v_south','long_name','Southern v-momentum component');
  ncwriteatt(bryfile,'v_south','units','meter second -1');
%
  nccreate(bryfile,'w_south','dimensions',{'xi_rho',Lp,'s_w',Np,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'w_south','long_name','Southern w-flux component');
  ncwriteatt(bryfile,'w_south','units','meter second -1');
%
  nccreate(bryfile,'ubar_south','dimensions',{'xi_u',L,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'ubar_south','long_name','Southern integrated u-momentum component');
  ncwriteatt(bryfile,'ubar_south','units','meter second -1');
%
  nccreate(bryfile,'vbar_south','dimensions',{'xi_rho',Lp,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'vbar_south','long_name','Southern integrated v-momentum component');
  ncwriteatt(bryfile,'vbar_south','units','meter second -1');
%
  nccreate(bryfile,'zeta_south','dimensions',{'xi_rho',Lp,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'zeta_south','long_name','Southern boundary sea surface height');
  ncwriteatt(bryfile,'zeta_south','units','meter');
end
%
%
if obcflag(2)==1  %%   Eastern boundary
%
  nccreate(bryfile,'temp_east','dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'temp_east','long_name','Eastern boundary potential temperature');
  ncwriteatt(bryfile,'temp_east','units','Celsius');
%
  nccreate(bryfile,'salt_east','dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'salt_east','long_name','Eastern boundary salinity');
  ncwriteatt(bryfile,'salt_east','units','PSU');
%
  nccreate(bryfile,'u_east','dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'u_east','long_name','Eastern u-momentum component');
  ncwriteatt(bryfile,'u_east','units','meter second -1');
%
  nccreate(bryfile,'v_east','dimensions',{'eta_v',M,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'v_east','long_name','Eastern v-momentum component');
  ncwriteatt(bryfile,'v_east','units','meter second -1');
%
  nccreate(bryfile,'w_east','dimensions',{'eta_rho',Mp,'s_w',Np,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'w_east','long_name','Eastern w-flux component');
  ncwriteatt(bryfile,'w_east','units','meter second -1');
%
  nccreate(bryfile,'ubar_east','dimensions',{'eta_rho',Mp,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'ubar_east','long_name','Eastern integrated u-momentum component');
  ncwriteatt(bryfile,'ubar_east','units','meter second -1');
%
  nccreate(bryfile,'vbar_east','dimensions',{'eta_v',M,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'vbar_east','long_name','Eastern integrated v-momentum component');
  ncwriteatt(bryfile,'vbar_east','units','meter second -1');
%
  nccreate(bryfile,'zeta_east','dimensions',{'eta_rho',Mp,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'zeta_east','long_name','Eastern boundary sea surface height');
  ncwriteatt(bryfile,'zeta_east','units','meter');
end
%
if obcflag(3)==1  %%   Northern boundary
  nccreate(bryfile,'temp_north','dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'temp_north','long_name','Northern boundary potential temperature');
  ncwriteatt(bryfile,'temp_north','units','Celsius');
%
  nccreate(bryfile,'salt_north','dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'salt_north','long_name','Northern boundary salinity');
  ncwriteatt(bryfile,'salt_north','units','PSU');
%
  nccreate(bryfile,'u_north','dimensions',{'xi_u',L,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'u_north','long_name','Northern u-momentum component');
  ncwriteatt(bryfile,'u_north','units','meter second -1');
%
  nccreate(bryfile,'v_north','dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'v_north','long_name','Northern v-momentum component');
  ncwriteatt(bryfile,'v_north','units','meter second -1');
%
  nccreate(bryfile,'w_north','dimensions',{'xi_rho',Lp,'s_w',Np,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'w_north','long_name','Northern w-flux component');
  ncwriteatt(bryfile,'w_north','units','meter second -1');
%
  nccreate(bryfile,'ubar_north','dimensions',{'xi_u',L,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'ubar_north','long_name','Northern integrated u-momentum component');
  ncwriteatt(bryfile,'ubar_north','units','meter second -1');
%
  nccreate(bryfile,'vbar_north','dimensions',{'xi_rho',Lp,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'vbar_north','long_name','Northern integrated v-momentum component');
  ncwriteatt(bryfile,'vbar_north','units','meter second -1');
%
  nccreate(bryfile,'zeta_north','dimensions',{'xi_rho',Lp,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'zeta_north','long_name','Northern boundary sea surface height');
  ncwriteatt(bryfile,'zeta_north','units','meter');
end
%
if obcflag(4)==1  %%   Western boundary
  nccreate(bryfile,'temp_west','dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'temp_west','long_name','Western boundary potential temperature');
  ncwriteatt(bryfile,'temp_west','units','Celsius');
%
  nccreate(bryfile,'salt_west','dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'salt_west','long_name','Western boundary salinity');
  ncwriteatt(bryfile,'salt_west','units','PSU');
%
  nccreate(bryfile,'u_west','dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'u_west','long_name','Western u-momentum component');
  ncwriteatt(bryfile,'u_west','units','meter second -1');
%
  nccreate(bryfile,'v_west','dimensions',{'eta_v',M,'s_rho',N,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'v_west','long_name','Western v-momentum component');
  ncwriteatt(bryfile,'v_west','units','meter second -1');
%
  nccreate(bryfile,'w_west','dimensions',{'eta_rho',Mp,'s_w',Np,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'w_west','long_name','Western w-flux component');
  ncwriteatt(bryfile,'w_west','units','meter second -1');
%
  nccreate(bryfile,'ubar_west','dimensions',{'eta_rho',Mp,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'ubar_west','long_name','Western integrated u-momentum component');
  ncwriteatt(bryfile,'ubar_west','units','meter second -1');
%
  nccreate(bryfile,'vbar_west','dimensions',{'eta_v',M,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'vbar_west','long_name','Western integrated v-momentum component');
  ncwriteatt(bryfile,'vbar_west','units','meter second -1');
%
  nccreate(bryfile,'zeta_west','dimensions',{'eta_rho',Mp,'bry_time',0},'datatype','single');
  ncwriteatt(bryfile,'zeta_west','long_name','Western boundary sea surface height');
  ncwriteatt(bryfile,'zeta_west','units','meter');
end
%
%
%  Write global attributes
%
ncwriteatt(bryfile,'/','Title',['Boundary file for' grdfile]);
ncwriteatt(bryfile,'/','Date',date);
%
% Write variables
%
ncwrite(bryfile,'theta_s',theta_s);
ncwrite(bryfile,'theta_b',theta_b);
ncwrite(bryfile,'Tcline',hc);
ncwrite(bryfile,'hc',hc);

