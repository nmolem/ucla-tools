import scipy.io as sio
import numpy as np
from Scientific.IO.NetCDF import NetCDFFile as netcdf
from forder import forder

path = '/batavia/WEDDY/HIS/ZLEV/weddy_z.'
path2= '/batavia/WEDDY/HIS/weddy_his.'

for fr in range(600,635):

   fnum = (fr-1)/4
   fnum = fnum*4
   ninFile = fr - fnum -1   # zero based indexing

   matfile = path + '{0:04d}'.format(fr) + '.mat'
   ncfile  = path + '{0:04d}'.format(fr) + '.nc'
   ncfile2 = path2+ '{0:04d}'.format(fnum) + '.nc'

   nc2= netcdf(ncfile2,'r')
   time = nc2.variables['ocean_time']
   otime = time[ninFile]
   nc2.close()

#  a = asfortranarray(nc.variable['var'])
#  print otime

   print(matfile)
#  matdata = sio.loadmat(matfile)
#  (nz,ny,nx) = matdata['sz'].shape

#  a = matdata['sz']
#  b =a.reshape( (nz,ny,nx), order='F')


   nc = netcdf(ncfile,'r')
#  a = nc.variables['sz'][:]
   a = np.asfortranarray(nc.variables['sz'])
   (nt,nz,ny,nx) = a.shape
   print a.shape
   a.resize( (nz,ny,nx) )
   (nz,ny,nx) = a.shape
   print a[1,1,1]
   forder(a,nx,ny,nz)
   print a[1,1,1]
   break

   nc.createDimension('time',1)
   nc.createDimension('depth',nz)
   nc.createDimension('ny',ny)
   nc.createDimension('nx',nx)

   nc_otime = nc.createVariable('ocean_time','d',('time',))
   nc_otime[:] = otime
#
   nc_z1 = nc.createVariable('depth','d',('depth',))
   nc_z1[:] = matdata['z1'].reshape(nz,)
#
   nc_sh = nc.createVariable('ssh','d',('time','ny','nx'))
   nc_sh[:,:] = matdata['ssh']
#
   nc_sz = nc.createVariable('sz','d',('time','depth','ny','nx'))
   nc_sz[:,:,:] = matdata['sz']
#
   nc_tz = nc.createVariable('tz','d',('time','depth','ny','nx'))
   nc_tz[:,:,:] = matdata['tz']
#
   nc_rz = nc.createVariable('rz','d',('time','depth','ny','nx'))
   nc_rz[:,:,:] = matdata['rz']
#
   nc_rz = nc.createVariable('uz','d',('time','depth','ny','nx'))
   nc_rz[:,:,:] = matdata['uz']
#
   nc_rz = nc.createVariable('vz','d',('time','depth','ny','nx'))
   nc_rz[:,:,:] = matdata['vz']
#
   nc_rz = nc.createVariable('wz','d',('time','depth','ny','nx'))
   nc_rz[:,:,:] = matdata['wz']
#
   nc.close()
