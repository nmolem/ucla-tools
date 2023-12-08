       subroutine get_h_coef(xp,yp,xc,yc,Aidx,Aval,nxp,nyp,nxc,nyc)
!
       implicit none 
!      import/export
       real(kind=8)   ,dimension(nyp,nxp)    ,intent(in)  :: xp,yp
       real(kind=8)   ,dimension(nyc,nxc)    ,intent(in)  :: xc,yc
       real(kind=8)   ,dimension(nyc,nxc,2,2),intent(out) :: Aval
       integer(kind=8),dimension(nyc,nxc,2)  ,intent(out) :: Aidx
       integer(kind=8)                       ,intent(in)  :: nxp,nyp,nxc,nyc

!f2py intent(in)  xp,yp
!f2py intent(in)  xc,yc
!f2py intent(out) Aidx
!f2py intent(out) Aval
!f2py intent(in)  nxp,nyp
!f2py intent(in)  nxc,nyc
!      local
       integer(kind=8)               :: ip,jp,i,j
       real(kind=8)   ,dimension(2,2):: w
       real(kind=8)                  :: eps
       integer(kind=8)               :: i0,i1,j0,j1
!
       eps = 1e-6
       ip = 1
       jp = 1
       do i = 1,nxc
         do j = 1,nyc
           call findll(xp,yp,xc(j,i),yc(j,i),ip,jp,nxp,nyp)

           call bilinear(xp(jp:jp+1,ip:ip+1),yp(jp:jp+1,ip:ip+1),xc(j,i),yc(j,i),w)

           Aval(j,i,:,:) = w
           Aidx(j,i,1) = jp
           Aidx(j,i,2) = ip
         enddo
       enddo

       end
!----------------------------------------------------------------------------------------------
       subroutine findll(x,y,x0,y0,i,j,nx,ny)
!
!      Finds the index (i,j) of the lower left corner of the 
!      enclosing box around x0,y0. Incoming i,j are used as initial guess
!
       implicit none 
!      import/export
       real(kind=8),dimension(ny,nx),intent(in)   :: x,y
       real(kind=8)                 ,intent(in)   :: x0,y0
       integer(kind=8)              ,intent(inout):: i,j
       integer(kind=8)              ,intent(in)   :: nx,ny
!f2py intent(in)    x,y
!f2py intent(in)    x0,y0
!f2py intent(inout) i,j
!f2py intent(in)    nx,ny
!      local
       integer(kind=8)  :: done,idir,outx,outy,it
       real(kind=8)     :: det,orient,new_orient
!
       done = 0 
       idir = 1 
       outx = 0
       outy = 0

       do it = 1, nx+ny !! set to something related to nx,ny
         if (done.gt.1) exit
         if (idir) then !! move in the i-direction
           orient = (y0 - y(j,i))*(x(j+1,i) - x(j,i))-(x0 - x(j,i))*(y(j+1,i) - y(j,i)) 
!          print *,orient, i,j,x(j,i),y(j,i)
           if (orient.le.0) then !! try a line to the right
             new_orient = (y0 - y(j,i+1))*(x(j+1,i+1) - x(j,i+1)) - (x0 - x(j,i+1))*(y(j+1,i+1) - y(j,i+1))
             if (new_orient.le.0) then !! We're still on the right-hand side
               if (i+1.eq.nx) then
!                print *, 'already at the right end'
                 done = done+1
                 idir = 0
                 outx = 1
               else
                 i = i+1 
                 done = 0
                 outx = 0
               endif
             else
!              print *, 'We have found an i'
!              print *,new_orient, i,j,x(j,i+1),y(j,i+1)
               done = done+1
               idir = 0;  !! change direction
             endif
           else       !! We're on the left of the line
             if (i.eq.1) then
!              print *, 'already at the left end'
               done = done+1
               outx = 1
               idir = 0
             else
               i = i-1 
               done = 0
               outx = 0
             endif
           endif
         else  !! move in the j-direction
           orient = (y0 - y(j,i))*(x(j,i+1) - x(j,i)) - (x0 - x(j,i))*(y(j,i+1) - y(j,i))
           if (orient.ge.0) then !! We're 'above', try a line higher
             new_orient = (y0 - y(j+1,i))*(x(j+1,i+1) - x(j+1,i)) - (x0 - x(j+1,i))*(y(j+1,i+1) - y(j+1,i))
             if (new_orient.ge.0) then  !! We're still 'above' the local line at j+1
               if (j+1.eq.ny) then
!                print *, 'already at the top end'
                 done = done+1
                 outy = 1
                 idir = 1
               else
                 j = j+1;
                 done = 0
                 outy = 0
               endif
             else     !! The line at j+1 is above us; change direction of search
!              print *, 'we are in place'
               done = done+1
               idir = 1;  !! change direction
             endif
           else  !! We're 'below'
             if (j.eq.1) then
!              print *, 'already at the bottom end'
               done = done+1
               outy = 1
               idir = 1
             else
               j = j-1;
               done = 0
               outy = 0
             endif
           endif
         endif
       enddo
!
       if (outx.gt.0.or.outy.gt.0) then
          print *, 'outside parent' 
          stop
       endif
!
       end
!----------------------------------------------------------------------------------------------
       subroutine interp_mat_2d(fp,fi,Aidx,Aval,nxp,nyp,nxc,nyc)
!
!      uses the pre-computed indices and weight to interpolated fp to fi
!
       implicit none 
!      import/export
       real(kind=8)   ,dimension(nyp,nxp)    ,intent(in)  :: fp
       real(kind=8)   ,dimension(nyc,nxc)    ,intent(out) :: fi
       real(kind=8)   ,dimension(nyc,nxc,2,2),intent(in)  :: Aval
       integer(kind=8),dimension(nyc,nxc,2  ),intent(in)  :: Aidx
       integer(kind=8)                       ,intent(in)  :: nxp,nyp
       integer(kind=8)                       ,intent(in)  :: nxc,nyc
!f2py intent(in)    fp
!f2py intent(out)   fi
!f2py intent(in)    Aval,Aidx
!f2py intent(in)    nxp,nyp
!f2py intent(in)    nxc,nyc
!      local
       integer(kind=8)  :: i,j,il,jl
!
       do i = 1, nxc 
         do j = 1, nyc
           jl = Aidx(j,i,1)
           il = Aidx(j,i,2)
           fi(j,i) = sum( fp(jl:jl+1,il:il+1)*Aval(j,i,:,:) )
         enddo
       enddo

       end
!----------------------------------------------------------------------------------------------
       subroutine bilinear(x,y,xc,yc,w)
       implicit none
!      import/export
       real(kind=8),dimension(4),intent(in) :: x,y
       real(kind=8)             ,intent(in) :: xc,yc
       real(kind=8),dimension(4),intent(out):: w
!      local
       real(kind=8) :: A,B,C,s,t
!
!      p1:  (jp  ,ip  )
!      p2:  (jp  ,ip  )
!      p3:  (jp  ,ip  )
!      p4:  (jp+1,ip+1)

       A = (x(1)-xc)*(y(1)-y(3)) - (y(1)-yc)*(x(1)-x(3))
       B = 0.5*( (x(1)-xc)*(y(2)-y(4)) - (y(1)-yc)*(x(2)-x(4))+(x(2)-xc)*(y(1)-y(3)) - (y(2)-yc)*(x(1)-x(3)) )
       C = (x(2)-xc)*(y(2)-y(4)) - (y(2)-yc)*(x(2)-x(4))

!      print *, 'Bilinear'
!      print *,(A-B + sqrt(B*B-A*C) )/(A-2*B + C)
!      print *,(A-B - sqrt(B*B-A*C) )/(A-2*B + C)

       s = (A-B + sqrt(B*B-A*C) )/(A-2*B + C)
       if ( (s.gt.1).or.(s.lt.0) ) then
         s = (A-B - sqrt(B*B-A*C) )/(A-2*B + C)
       endif

       t = ( (1-s)*(x(1)-xc) + s*(x(2)-xc) ) / ( (1-s)*(x(1)-x(3)) + s*(x(2)-x(4)) )

       w(1) = (1-s)*(1-t);
       w(2) =    s*(1-t);
       w(3) = (1-s)*t;
       w(4) =    s*t;

       end
!----------------------------------------------------------------------------------------------
       subroutine interp_mat_3d(fp,fi,Aidx,Aval,nxp,nyp,nzp,nxc,nyc,nzc)
!
!      uses the pre-computed indices and weight to interpolated fp to fi
!
       implicit none 
!      import/export
       real(kind=8)   ,dimension(nzp,nyp,nxp)  ,intent(in)  :: fp
       real(kind=8)   ,dimension(nzc,nyc,nxc)  ,intent(out) :: fi
       real(kind=8)   ,dimension(nyc,nxc,2,2)  ,intent(in)  :: Ahval
       integer(kind=8),dimension(nyc,nxc,2  )  ,intent(in)  :: Ahidx
       real(kind=8)   ,dimension(nzc,nyc,nxc,2),intent(in)  :: Avval
       integer(kind=8),dimension(nzc,nyc,nxc)  ,intent(in)  :: Avidx
       integer(kind=8)                         ,intent(in)  :: nxp,nyp,nzp
       integer(kind=8)                         ,intent(in)  :: nxc,nyc,nzc
!f2py intent(in)    fp
!f2py intent(out)   fi
!f2py intent(in)    Avval,Avidx
!f2py intent(in)    Ahval,Ahidx
!f2py intent(in)    nzp,nyp,nxp
!f2py intent(in)    nzc,nyc,nxc
!      local
       integer(kind=8)  :: i,j,k,il,jl,kl
!
       do i = 1, nxc 
         do j = 1, nyc
           jl = Aidx(j,i,1)
           il = Aidx(j,i,2)
           do k = 1, nzc
             kl = Avidx(k,j,i)
             fi(k,j,i) = sum( fp(kl  ,jl:jl+1,il:il+1)*Aval(j,i,:,:))*Avval(k,j,i,1) +  &
                         sum( fp(kl+1,jl:jl+1,il:il+1)*Aval(j,i,:,:))*Avval(k,j,i,2)
           enddo
         enddo
       enddo

       end
!----------------------------------------------------------------------------------------------
