      subroutine forder(a,dum,nx,ny,nz)
      implicit none
*     IMPORT/EXPORT
      integer(kind=8),                    intent(in)   :: nx,ny,nz
      real(kind=8)   ,dimension(nz,ny,nx),intent(inout):: a
      real(kind=8)                       ,intent(out):: dum
cf2py intent(in)    nx
cf2py intent(in)    ny
cf2py intent(in)    nz
cf2py intent(inout) a
cf2py intent(out) dum
*     LOCAL
      integer i,j,k,cnt,n
      integer(kind=8),allocatable,dimension(:,:,:) :: b
*
      dum = 666.
      n = 3
      allocate(b(n,n,n))
      cnt = 0
      do  k = 1,n
       do  j = 1,n
        do i = 1,n
          b(i,j,k) = cnt
          cnt = cnt+1
        enddo 
       enddo
      enddo 

      a(1:3,1:3,1:3) = b
      dum =  666.
      deallocate(b)
      print *,dum
*
      end
