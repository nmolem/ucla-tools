!----------------------------------------------------------------------------------------------
      subroutine gnomonic(lon,lat,xg,yg,nx,ny)
!
!     Gnomonic projection is distance preserving for
!     interpolation purposes.
!
!     Input:  lon,lat in degrees.
!
      implicit none

!     import/export
      integer(kind=8)              ,intent(in) :: nx,ny
      real(kind=8),dimension(ny,nx),intent(in) :: lon,lat
      real(kind=8),dimension(ny,nx),intent(out):: xg,yg
!     local
      real(kind=8),parameter :: pi = 3.141592653589793
      real(kind=8)           :: dlon,dlat,lon0,lat0
      real(kind=8),dimension(:,:),allocatable :: lonr,latr,cosc

       dlon = maxval(lon)-minval(lon)
       dlat = maxval(lat)-minval(lat)
       if (  (dlon.gt.100).or.(dlat.gt.100)  ) then
        print *, 'This area is too large for gnomonic projections!!'
        xg = lon 
        yg = lat 
        return
       endif

       allocate(lonr(ny,nx))
       allocate(latr(ny,nx))
       allocate(cosc(ny,nx))
       lonr = lon*pi/180 
       latr = lat*pi/180 
       lon0 = sum(lonr)/(nx*ny)
       lat0 = sum(latr)/(nx*ny)

       cosc = sin(lat0)*sin(latr) + cos(lat0)*cos(latr)*cos(lonr-lon0)

       xg = cos(latr)*sin(lon-lon0)/cosc
       yg = (cos(lat0)*sin(latr) - sin(lat0)*cos(latr)*cos(lonr-lon0) )/cosc

       deallocate(lonr)
       deallocate(latr)
       deallocate(cosc)

       end
