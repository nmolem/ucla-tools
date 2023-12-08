       subroutine get_v_coef(zp,zc,Avidx,Avval)
       implicit none 
!      import/export
       real(kind=8)   ,dimension(nzp,nyp,nxp)    ,intent(in)  :: zp
       real(kind=8)   ,dimension(nzc,nyc,nxc)    ,intent(in)  :: zc
       real(kind=8)   ,dimension(nzc,nyc,nxc,2,2),intent(out) :: Avval
       integer(kind=8),dimension(nzc,nyc,nxc,2)  ,intent(out) :: Avidx
       integer(kind=8)                           ,intent(in)  :: nzp,nxp,nyp
       integer(kind=8)                           ,intent(in)  :: nzc,nxc,nyc
!f2py intent(in)  zp,zc
!f2py intent(in)  xc,yc
!f2py intent(out) Avidx
!f2py intent(out) Avval
!f2py intent(in)  nxp,nyp,nzp
!f2py intent(in)  nxc,nyc,nzc
!      local
       integer(kind=8) :: i,j,k,kp
       real(kind=8)   ,dimension(nzp,nyp,nxp)    ,intent(in)  :: zp

       allocate(zpi(nzp,nyc,nxc))
       do k = 1,nzp
         call interp_mat_2d(zp(k,:,:),zpi(k,:,:),Aidx,Aval,nxp,nyp,nxc,nyc)
       enddo

       do i = 1,nxc
         do j = 1,nyc
           kp = 1
           do k = 1,nzc
             do while (zpi(kp).lt.zc(k).and.kp.lt.nzp) 
               ip = ip+1
             enddo
             if (kp.eq.1.or.zpi(ip).lt.zc(k))
               !! zc is below zpi(1) or zc is above zpi(nzp)
               Avidx(k,j,i) = ip
               Avval(k,j,i,1) = 1.0
               Avval(k,j,i,2) = 0.0
             else
               alp = (zc(ic)-zp(ip-1))/(zp(ip)-zp(ip-1));
               Avidx(k,j,i) = ip-1
               Avval(k,j,i,1) = 1.0-alp
               Avval(k,j,i,2) = alp
             endif
           enddo
         enddo
       enddo
       deallocate(zpi)
*          
    end

