!==================================================================

      SUBROUTINE initpop 
      use var_inc
      implicit none

      real, dimension(lx,ly,lz):: usqr, G
      integer ip

      usqr = ux*ux + uy*uy + uz*uz
      usqr = 1.5*usqr   

! initialise density for incompressible formulation
!   
!     rho = rho0
      rho = 0.0

      f(0,:,:,:) = ww0*(rho - usqr) 

      do ip = 1,6
        G = (cix(ip)*ux + ciy(ip)*uy + ciz(ip)*uz)
        f(ip,:,:,:) = ww1*(rho + 3.0*G + 4.5*G*G - usqr) 
      end do

      do ip = 7,npop-1
        G = (cix(ip)*ux + ciy(ip)*uy + ciz(ip)*uz)
        f(ip,:,:,:) = ww2*(rho + 3.0*G + 4.5*G*G - usqr) 
      end do

      END SUBROUTINE initpop
!==================================================================

      SUBROUTINE initrand(inseed)
      use mpi
      use var_inc
      implicit none
 
      integer(kind = 4) irand
      integer i, myseed, inseed  
      integer, dimension(nproc):: iseed 
 
      if(myid == 0)then
        call srand(inseed)
        do i = 1,nproc
          iseed(i) = irand()
        end do 
      end if
 
      call MPI_BCAST(iseed,nproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
      myseed = iseed(myid + 1)
      call srand(myseed) 
     
      END SUBROUTINE initrand
!==================================================================

      SUBROUTINE initvel
      use mpi
      use var_inc
      implicit none

      integer i,j,k,jj,kk
      real yplus,y9,z9,A9,u9,alpha,beta9,gamma1,gamma2,phi1,phi2
      real ccc1, ccc9, ccc10, cc

!     A9 = 0.1*sqrt(2.0)
      A9 = 0.0
      alpha = 1.0
      beta9 = 1.0
      gamma1 = 2.0
      gamma2 = 1.0
      phi1 = 0.111
      phi2 = 0.555
      cc = 60.0

      ccc1 = -real(ny)/pi2/alpha/ystar*A9*ustar/cc/cc

! mean flow
      ux = 0.0
      uz = 0.0
      do i=1,nxh
      yplus = (real(i)-0.5)/ystar

      if(yplus.lt.10.8)then
      uy(i,:,:) = yplus*ustar
      uy(nx+1-i,:,:) = yplus*ustar
      else
      u9 = alog(yplus)/0.41 + 5.0
      u9 = u9*ustar
      uy(i,:,:) = u9
      uy(nx+1-i,:,:) = u9
      end if
      end do 


!     add some perturbation
      do k=1,lz
      kk = k + indz*lz
      z9 = pi2*(real(kk)-0.5)/real(nz)
      do j=1,ly
      jj = j + indy*ly
      y9 = pi2*(real(jj)-0.5)/real(ny)
      do i=1,nxh
      yplus = (real(i)-0.5)/ystar
      ccc9 = exp(-yplus/cc)
      u9 =  ccc1*yplus*ccc9*sin(alpha*y9 + beta9*z9)
      uy(i,j,k) = uy(i,j,k) + u9
      uy(nx+1-i,j,k) = uy(nx+1-i,j,k) + u9

      ccc10 = A9*ustar*(1.-ccc9-yplus/cc*ccc9)

      u9 = ccc10*cos(alpha*y9 + beta9*z9)
      ux(i,j,k) = ux(i,j,k) + u9
      ux(nx+1-i,j,k) = ux(nx+1-i,j,k) + u9

!     u9 = 0.5*ccc10/sqrt(2.0)*(sin(gamma1*y9+pi2*phi1)+sin(gamma2*y9+pi2*phi2) )
!     u9 = ccc10*sin(gamma1*y9+pi2*phi1)
!     uz(i,j,k) = uz(i,j,k) +  u9
!     uz(nx+1-i,j,k) = uz(nx+1-i,j,k) + u9

      end do
      end do
      end do

      RETURN
      END SUBROUTINE initvel    




