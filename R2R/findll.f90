       subroutine findll(x,y,x0,y0,i,j,nx,ny)
!
!      Finds the index (i,j) of the lower left corner of the 
!      enclosing box around x0,y0. Incoming i0,j0 are used as initial guess
!
       implicit none 
!      import/export
       real(kind=8),dimension(ny,nx),intent(in)   :: x,y
       real(kind=8)                 ,intent(in)   :: x0,y0
       integer(kind=8)              ,intent(inout):: i,j
       integer(kind=8)              ,intent(in)   :: nx,ny
cf2py intent(in)    x,y
cf2py intent(in)    x0,y0
cf2py intent(inout) i,j
cf2py intent(in)    nx,ny
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
*          print *,orient, i,j,x(j,i),y(j,i)
           if (orient.le.0) then !! try a line to the right
             new_orient = (y0 - y(j,i+1))*(x(j+1,i+1) - x(j,i+1)) - (x0 - x(j,i+1))*(y(j+1,i+1) - y(j,i+1))
             if (new_orient.le.0) then !! We're still on the right-hand side
               if (i+1.eq.nx) then
                 print *, 'already at the right end'
                 exit
                 done = done+1
                 idir = 0
                 outx = 1
               else
                 i = i+1 
                 done = 0
                 outx = 0
               endif
             else
               print *, 'We have found an i'
*              print *,new_orient, i,j,x(j,i+1),y(j,i+1)
               done = done+1
               idir = 0;  !! change direction
             endif
           else       !! We're on the left of the line
             if (i.eq.1) then
               print *, 'already at the left end'
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
                 print *, 'already at the top end'
                 done = done+1
                 outy = 1
                 idir = 1
               else
                 j = j+1;
                 done = 0
                 outy = 0
               endif
             else     !! The line at j+1 is above us; change direction of search
               print *, 'we are in place'
               done = done+1
               idir = 1;  !! change direction
             endif
           else  !! We're 'below'
             if (j.eq.1) then
               print *, 'already at the bottom end'
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

       print *, 'done: ',done
       print *, 'outx: ', outx
       print *, 'outy: ', outy
       print *, 'i,j:  ', i,j
!
       end



