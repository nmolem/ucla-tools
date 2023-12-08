      subroutine r2r_bry_subgrid(lonp,latp,lonc,latc,nxp,nyp,nxc,nyc,obcflag)
      implicit none
*     import/export
      real(kind=8),dimension(nxp,nyp) :: lonp,latp
      real(kind=8),dimension(nxc,nyc) :: lonc,latc
      integer(kind=8),dimension(4)    :: obcflag
*
      do bnd = 1,4
        if (.not.obcflag(bnd) )
          continue
        else
          allocate(blon(2,nyc))
          allocate(blat(2,nyc))
          blon = lonc(1:2,1:end)
          blat = latc(1:2,1:end)
          deallocate(blon)
          deallocate(blat)
        endif
      enddo
      


