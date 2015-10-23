!=============================================================================
!@subroutine pipe_links
!@desc Determines which node links cross the pipe wall boundary and the distance
!      the boundary is from the node. This data is then used later to execute
!	   interpolation bounceback. This should only need to be called once.
!=============================================================================
	subroutine pipe_links
	use mpi
	use var_inc
	implicit none

	integer i, j, k, imove, jmove, kmove, ip
	real xx0, zz0, xx1, zz1, rr0, rr01, rr1, rr11, pradm2
	real aa, bb, cc, alpha, alpha0, alpha1, alphaz
		
	pradm2 = prad - 2.d0
	alphaz = dfloat(indz*lz)
!        nlink = 0
	do i = 1, lx
	do k = 1, lz
	  xx0 = dfloat(i) - 0.5d0 - pxcenter
          zz0 = dfloat(k) - 0.5d0 + alphaz - pzcenter

	  rr0 = xx0*xx0 + zz0*zz0
	  rr01 = rr0 - (prad*1.d0)**2

	  if(rr01 .ge. 0)then
	     ibnodes(i,:,k) = 1

	     else if(rr0 .ge. pradm2**2)then
	  !Determine which node links are crossing the solid boundary
	   do ip = 1,npop-1
	      imove = i + cix(ip) 
	      kmove = k + ciz(ip)

	      xx1 = dfloat(imove) - 0.5d0 - pxcenter
	      zz1 = dfloat(kmove) - 0.5d0 + alphaz - pzcenter
	
	      rr1 = xx1**2 + zz1**2  
              rr11 = rr1 - (prad*1.d0)**2  

		if(rr11 .ge. 0)then
		  do j = 1, ly
		  nlink = nlink + 1
	          !Save link information
	          xlink(nlink) = i
	          ylink(nlink) = j
	          zlink(nlink) = k

	          iplink(nlink) = ip
	          mlink(nlink) = -1

  !Calculate the position of the solid boundary on the link
  ![x0+alpha(x1-x0)]^2 + [z0+alpha(z1-z0)]^2 = prad^2
  !Quadratic formula
	        aa = rr0 + rr1 - 2.d0*(xx0*xx1 + zz0*zz1)
	        bb = 2.d0*(xx0*(xx1-xx0) + zz0*(zz1-zz0))    
	        cc = rr0 - prad**2 
!                if(bb*bb-4.d0*aa*cc .lt. 0) write(*,*)'myid,i,j,k,ip = ',myid,i,j,k,ip
	        alpha1 = dsqrt(bb*bb - 4.d0*aa*cc)    
	        alpha = (-bb + alpha1)/(2.d0*aa)

	        alink(nlink) = alpha
!                if(myid.eq.0) write(*,*)i,j,k,ip,alpha
	                      !Check to see if we will need data from neighboring MPI tasks for interpolation
	                      call parse_MPI_links(ip, i, j, k, alpha, nlink) !@file partlib.f90
                      enddo
					endif

		          enddo
		       endif

			enddo
		enddo


		end subroutine pipe_links
