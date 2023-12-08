import scipy.io as sio
import numpy as np
from Scientific.IO.NetCDF import NetCDFFile as netcdf
from get_coef import *

#print  get_h_coef.__doc__
#print  findll.__doc__
#print  interp_mat_2d.__doc__

parentgrid = '/batavia/GULFS/gulfs_grd.nc'
childgrid  = '/batavia/WEDDY/weddy_grd.nc'

nzp = 40;
nzc = 40;

nc = netcdf(parentgrid,'r')
lonp = np.asfortranarray(nc.variables['lon_rho'][:])
latp = np.asfortranarray(nc.variables['lat_rho'][:])
nc.close()

nc = netcdf(childgrid,'r')
lonc = np.asfortranarray(nc.variables['lon_rho'][:])
latc = np.asfortranarray(nc.variables['lat_rho'][:])
nc.close()

aidx,aval = get_h_coef(lonp,latp,lonc,latc)

j0 = np.min(aidx[:,:,0])
j1 = np.max(aidx[:,:,0])
i0 = np.min(aidx[:,:,1])
i1 = np.max(aidx[:,:,1])

nc = netcdf(childgrid,'r')
lonc = np.asfortranarray(nc.variables['lon_rho'][:])
latc = np.asfortranarray(nc.variables['lat_rho'][:])
nc.close()

zp = set_deth(h,thetasp,thetabp,hcp,nzp,'r')
zc = set_deth(h,thetasp,thetabp,hcp,nzp,'r')

fi = interp_mat_2d(lonp,aidx,aval)
sio.savemat('fi_data.mat', mdict={'fi': fi})





#findll(lonp,latp,x0,y0,i,j)

