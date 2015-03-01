!===========================================================================
! this is to monitor the mean and maximum flow velocity and particle
! velocity
      subroutine diag
      use mpi
      use var_inc
      implicit none

!********THIS IS CHANGED*******************************
!      integer, dimension(2):: idwp, idomgp

      integer, dimension(1):: idwp, idomgp

      real ttt, vmax, wpmax, omgpmax,volf
      real, dimension(lx,ly,lz):: vel
      real, dimension(npart):: wpmag, omgpmag
      real, dimension(nproc):: vmax0,vmax0t
      integer, dimension(nproc):: im,jm,km,im0,jm0,km0
      real umean,vmean,wmean,urms,vrms,wrms,urmst,vrmst,wrmst
      real umeant,vmeant,wmeant
      integer i,j,k,imout,jmout,kmout,nfluid0,nfluid

      character (len = 120):: fnm
!!!!!!!!!!!
      umean = 0.0
      vmean = 0.0
      wmean = 0.0

      nfluid0 = count(ibnodes < 0)
      umean = sum (ux,MASK = (ibnodes < 0) )
      vmean = sum (uy,MASK = (ibnodes < 0) )
      wmean = sum (uz,MASK = (ibnodes < 0) )
      urms = sum (ux*ux,MASK = (ibnodes < 0) )
      vrms = sum (uy*uy,MASK = (ibnodes < 0) )
      wrms = sum (uz*uz,MASK = (ibnodes < 0) )

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      CALL MPI_ALLREDUCE(nfluid0,nfluid,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(umean,umeant,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(vmean,vmeant,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(wmean,wmeant,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(urms,urmst,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(vrms,vrmst,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(wrms,wrmst,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)

!!!!!!!!!!!

      vel = sqrt(ux*ux + uy*uy + uz*uz)
      where(ibnodes > 0) vel = 0.0

!     vmax0 = maxval(vel)
      im = 0
      jm = 0
      km = 0
      vmax0 = 0.0
      do k=1,lz
      do j=1,ly
      do i=1,lx
      if(vel(i,j,k).gt.vmax0(myid) ) then
      vmax0(myid) = vel(i,j,k)
      im(myid)=i
      jm(myid)=j +  indy*ly
      km(myid)=k +  indz*lz
      end if
      end do
      end do
      end do
!
!     call MPI_REDUCE(vmax0,vmax,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(im,im0,nproc,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(jm,jm0,nproc,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(km,km0,nproc,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(vmax0,vmax0t,nproc,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!     if(myid == 0)then
!       wpmax = 0.0
!       omgpmax = 0.0
!       idwp = 0
!       idomgp = 0

!       if(ipart)then
!         wpmag = sqrt(wp(1,:)**2 + wp(2,:)**2 + wp(3,:)**2)
!         wpmax = maxval(wpmag)
!         idwp = maxloc(wpmag)

!         omgpmag = sqrt(omgp(1,:)**2 + omgp(2,:)**2 + omgp(3,:)**2)
!         omgpmax = maxval(omgpmag)
!         idomgp = maxloc(omgpmag)
!       end if

!     end if

          rhoerr = maxval(rho,MASK = (ibnodes < 0))
          call MPI_ALLREDUCE(rhoerr,rhomax,1,MPI_REAL8,MPI_MAX, &
                             MPI_COMM_WORLD,ierr)
          rhoerr = minval(rho,MASK = (ibnodes < 0))
          call MPI_ALLREDUCE(rhoerr,rhomin,1,MPI_REAL8,MPI_MIN, &
                             MPI_COMM_WORLD,ierr)
          if(myid == 0 ) write(*,*)istep, rhomax, rhomin

      if(myid == 0)then
       umeant = umeant/real(nfluid)
       vmeant = vmeant/real(nfluid)
       wmeant = wmeant/real(nfluid)
       urmst = sqrt(urmst/real(nfluid) - umeant**2)
       vrmst = sqrt(vrmst/real(nfluid) - vmeant**2)
       wrmst = sqrt(wrmst/real(nfluid) - wmeant**2)

       umeant = umeant/ustar
       vmeant = vmeant/ustar
       wmeant = wmeant/ustar
       urmst = urmst/ustar
       vrmst = vrmst/ustar
       wrmst = wrmst/ustar

       volf = 1. - real(nfluid)/real(nx*ny*nz)

        fnm = trim(dirdiag)//'diag.dat'

        open(26, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted', position = 'append')

        ttt = real(istep)

        vmax = 0.0
        do i = 1, nproc
        if(vmax0t(i).gt.vmax)then
        vmax = vmax0t(i)
        imout = im0(i)
        jmout = jm0(i)
        kmout = km0(i)
        end if

        end do


       write(26,260) ttt,vmax,imout,jmout,kmout,umeant,vmeant,wmeant,urmst,vrmst,wrmst,volf, &
                rhomax,rhomin
       write(*,*)'ttt,vmax,imout,jmout,kmout,umean,vmean,wmean,umeans,vmeans,wmeans=', &
         ttt,vmax,imout,jmout,kmout,umeant,vmeant,wmeant,urmst,vrmst,wrmst,volf, &
                rhomax,rhomin
        close(26)

      end if

260   format(2x, 2(1pe16.6),3I6,9(1pe16.6) )

      end subroutine diag

