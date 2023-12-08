import scipy.io as sio
import numpy as np
from Scientific.IO.NetCDF import NetCDFFile as netcdf
from forder import forder

print forder.__doc__

path = '/batavia/WEDDY/HIS/ZLEV/weddy_z.'
path2= '/batavia/WEDDY/HIS/weddy_his.'

for fr in range(600,635):

   fnum = (fr-1)/4
   fnum = fnum*4
   ninFile = fr - fnum -1   # zero based indexing

   ncfile  = path + '{0:04d}'.format(fr) + '.nc'


   nc = netcdf(ncfile,'r')
   a = np.asfortranarray(nc.variables['sz'])
   (nt,nz,ny,nx) = a.shape
   print a.shape
   a.resize( (nz,ny,nx) )
   (nz,ny,nx) = a.shape
   dum = a[1,1,1]
   print a[1,1,1]
   dum = forder(a,nx,ny,nz)
   print a[1,1,1]
   print dum
   break

