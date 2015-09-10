!=============================================================================
!@subroutine para
!@desc Initializes parameters, variables, and constants
!=============================================================================      
      subroutine para
      use var_inc
      implicit none

      logical dirExist
      integer i,j,k
      real vmax
!=======================================================
! General parameters
!=======================================================
      ! true: if a entirely new run is desired
      ! false: wish to load an old flow
      newrun = .false.
      ! true: if one wants to load an old pre-relaxed flow
      ! false: wish to load an previous simulaiton
      newinitflow = .false. 

      ! Initial timestep one want to load flow/ particle data from
      ! It can be any step # at which the "endrunflow" & "endrunpart" data are saved
      istpload = 2000
      ! Number of time steps in main loop
      nsteps = 2000

      s = 1 !Global switch used for 2-array algorithm

      istep0 = 0
      istep00 = 1 
      istat = 0
      imovie = 0 

      ! wall clock time limit, must be the same as in RUN file, unit in [min]
      time_lmt = 3600.0
      ! wall clock time buffer, for saving the data, [min]
      time_buff = 10.0
      ! wall clock time upper bound, for computation, [s]
      time_bond = (time_lmt - time_buff)*60.d0
!=======================================================
! Physical constants and parameters
!=======================================================
      pi = 4.0*atan(1.0) 
      pi2 = 2.0*pi
      
      ! specify the force magnitude
!     visc = 0.0032
!     Rstar = 200.0
!     visc = 0.0036
      visc = 0.0036
      Rstar = 180.0
      ustar = 2.0*Rstar*visc/real(nx)
      force_in_y = 2.*rho0*ustar*ustar/real(nx)
      ystar = visc/ustar
      force_mag = 1.0

      ! not used
      kpeak = 4         ! It does not matter. Starting from stagnant flow
      u0in = 0.0503     ! It does not matter. Starting from stagnant flow

      iseedf = 232300

      ! get rid of vscale
      vscale = 1.0
      escale = 1.0/(vscale*vscale)
      dscale = pi2/real(ny)/vscale**3
      tscale = pi2*vscale/real(ny)

!     nek = nx/3
      nek = int(nx7/2 - 1.5)
 
      tau = 3.0*visc + 0.5
!=======================================================
! Initialize MRT related constants
!=======================================================
      MRTtype = 1

      s9 = 1.0/tau
      s13 = s9

      select case(MRTtype) 
      case(1) ! Linear analysis data for stability 
        s1 = 1.19
        s2 = 1.4
        s4 = 1.2
        s10 = 1.4
        s16 = 1.98

        omegepsl = 0.0
        omegepslj = -475.0/63.0
        omegxx = 0.0
      case(2) ! To recover LBGK  
        s1 = s9
        s2 = s9
        s4 = s9
        s10 = s9
        s16 = s9

        omegepsl = 3.0
        omegepslj = -11.0/2.0
        omegxx = -1.0/2.0
      case(3) ! To reduce dissipation
        s1 = 1.8
        s2 = s1
        s4 = s9
        s10 = s1
        s16 = s1

        omegepsl = 3.0
        omegepslj = -11.0/2.0
        omegxx = -1.0/2.0
      end select

      coef1 = -2.0/3.0
      coef2 = -11.0
      coef3 = 8.0
      coef4 = -4.0
      coef5 = 2.0

      coef3i = 1.0/coef3
      coef4i = 1.0/coef4

      val1 = 19.0
      val2 = 2394.0
      val3 = 252.0
      val4 = 10.0
      val5 = 40.0
      val6 = 36.0
      val7 = 72.0
      val8 = 12.0
      val9 = 24.0

      val1i = 1.0/val1 
      val2i = 1.0/val2
      val3i = 1.0/val3
      val4i = 1.0/val4
      val5i = 1.0/val5
      val6i = 1.0/val6
      val7i = 1.0/val7
      val8i = 1.0/val8
      val9i = 1.0/val9

      ww0 = 1.0/3.0
      ww1 = 1.0/18.0
      ww2 = 1.0/36.0

! The order must be exactly as in D'Humieres et al (2002), (A1)
!      Phil. Trans. R. Soc. Lond. A (2002) 360, 437-451
! index 0:  ( 0  0  0)
! index 1:  ( 1  0  0)
! index 2:  (-1  0  0)
! index 3:  ( 0  1  0)
! index 4:  ( 0 -1  0)
! index 5:  ( 0  0  1)
! index 6:  ( 0  0 -1)
!
! index 7:  (+1 +1  0)
! index 8:  (-1 +1  0)
! index 9:  (+1 -1  0)
! index 10: (-1 -1  0)
!
! index 11: (+1  0 +1)
! index 12: (-1  0 +1)
! index 13: (+1  0 -1)
! index 14: (-1  0 -1)
!
! index 15: ( 0 +1 +1)
! index 16: ( 0 -1 +1)
! index 17: ( 0 +1 -1)
! index 18: ( 0 -1 -1)

      cix = (/0, 1, -1, 0, 0, 0, 0,1,-1, 1,-1,1,-1, 1,-1,0, 0, 0, 0 /) 
      ciy = (/0, 0, 0, 1, -1, 0, 0,1, 1,-1,-1,0, 0, 0, 0,1,-1, 1,-1 /) 
      ciz = (/0, 0, 0, 0, 0,  1,-1,0, 0, 0, 0,1, 1,-1,-1,1, 1,-1,-1/) 
      
      !Arrray used for finding the opposite velocity
      ipopp=(/0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15/) 

!=======================================================
! Create MPI topology
!=======================================================
      nprocY = 20 !MPI topology width
      nprocZ = nproc/nprocY !MPI topology height

      ly = ny/nprocY         !local division of dist for procs in y dir
      lz = nz/nprocZ         !local division of dist for procs in z dir

      if((mod(myid,nprocY)+1).eq.nprocY) then
        lyext=2
      else
        lyext=0
      endif
      lly = ly+lyext

      !Determine MPI neighbor Ids
      indy = mod(myid,nprocY)
      indz = int(myid/nprocY)

      mzp = mod(indz+1,nprocZ) * nprocY + indy  !top
      mzm = mod(indz + nprocZ - 1,nprocZ) * nprocY + indy !bottom

      myp = indz*nprocY + mod(indy+1,nprocY) !right
      mym = indz*nprocY + mod(indy+nprocY-1,nprocY) !left


      mypzp = mod(indz+1,nprocZ)*nprocY + mod(indy+1,nprocY) !top-right corner
      mypzm = mod(indz+nprocZ-1,nprocZ)*nprocY + mod(indy+1,nprocY) !bottom-right corner
      mymzp = mod(indz+1,nprocZ)*nprocY + mod(indy+nprocY-1,nprocY) !top-left corner
      mymzm = mod(indz+nprocZ-1,nprocZ)*nprocY + mod(indy+nprocY-1,nprocY) !bottom-left corner
	
      if(mzp == mzm .OR. myp == mym)then
        if(myid==0)then
          write(*,*)'WARNING: MPI topology is too small!'
          write(*,*)'Particle Laden code will not work!'
        endif
      endif

      rhoepsl = 1.e-05
!=======================================================
! Declare reading and writing directories
!=======================================================
      dircntdflow0 = trim('/glade/scratch/ngeneva/Particle_Laden_Channel_Flow_Code/')
      dircntdpart0 = trim('/glade/scratch/ngeneva/Particle_Laden_Channel_Flow_Code/')

      dirgenr = '/glade/scratch/ngeneva/Particle_Laden_Channel_Flow_Code/'
      dirdiag = trim(dirgenr)//'diag/'
      dirstat = trim(dirgenr)//'stat/'
      dirprobe = trim(dirgenr)//'probe/'
      dirinitflow = trim(dirgenr)//'initflow/'
      dircntdflow = trim(dirgenr)//'cntdflow/'
      dirflowout = trim(dirgenr)//'flowout/'
      dirmoviedata = trim(dirgenr)//'moviedata/'

      ! saving and loading directories relevant to particle
      dircntdpart = trim(dirgenr)//'cntdpart/'
      dirpartout = trim(dirgenr)//'partout/'

      !Benchmarking directories
      dirbench = trim(dirgenr)//'benchmark/'
      dirbenchmatlab = trim(dirbench)//'matlab/'
      dirbenchbead = trim(dirbench)//'bead/'
      dirbenchflow = trim(dirbench)//'flow/'

      !Make sure that directories exist with makedir
      !@file para.f90
      if(myid==0)then
        call makedir(dirdiag)
        call makedir(dirstat)
        call makedir(dirprobe)
        call makedir(dirinitflow)
        call makedir(dircntdflow)
        call makedir(dirflowout)
        call makedir(dirmoviedata)
        call makedir(dircntdpart)
        call makedir(dirpartout)

        call makedir(dirbench)
        call makedir(dirbenchbead)
        call makedir(dirbenchflow)
        call makedir(dirbenchmatlab)
      endif
!=======================================================
! Particle related parameters
!=======================================================
      ipart = .true.
      released = .FALSE.

      !Flag indicating node is inside solid particle
      !Used in the direct request system in beads_collision
      IBNODES_TRUE = 12345678

      if(ipart)then
        volp = 4.0/3.0*pi*rad**3
        amp = rhopart*volp
        aip = 0.4*amp*rad**2 ! moment of inertia for solid sphere

        ! rhog is used to calculate the (body force - buoyancy)
        ! Define g by setting Stokes settling velocity / ustar
        ! ws_normalized =  Stokes settling velocity / ustar
        ! (1) For with gravity case
!       ws_normalized = 1.0
!       g_lbm = 9.0/2.0*ws_normalized/(rhopart/rho0 - 1.0)  
!       g_lbm = g_lbm/rad**2*visc*ustar
!       rhog = (rhopart/rho0 - 1.0)*rho0*g_lbm
        ! (2) For no gravity case
        rhog = 0.0
        g_lbm = ustar*ustar/ystar

        ! stiff coefficient from Feng & Michaelides (2005) JCP 202, pp 20-51, eqn(28)
        stf0 = 0.025
        stf1 = 0.002
        stf0_w = 0.025
        stf1_w = 0.002
!       stf0 = 0.01
!       stf1 = 0.001
!       stf0_w = 0.01
!       stf1_w = 0.001
        
        !Random seed for inital particle placement
        iseedp = 12345
        !Max number of solid particles that can exist in a single local MPI domain (NEEDS FIX)
        msize = int(10.0*real(npart)/real(nproc))
        ! max number of node that will require filling
        ! umax is determined from poiseuille flow in the y direction
        vmax = ((force_in_y*force_mag*nx**2)/(8*visc))
        maxbfill = vmax*msize*(pi*rad**2)
        if(myid == 0)write(*,*) maxbfill
        ! maximum number of boundary links per local domain
        ! 5 = avg number of vertexes a plane crosses in D3Q19
        maxlink = 5*msize*(4*pi*rad**2)
 
        wwp = (/ww1, ww1, ww1, ww1, ww1, ww1, ww2, ww2, ww2,           &
                ww2, ww2, ww2, ww2, ww2, ww2, ww2, ww2, ww2/) 

      end if
!=======================================================
! Computing wave numbers and other values for the forcing
! Not used in particle laden
!=======================================================
      do i = 1, nx, 2
       kxr(i)   = int(i/2)
       kxr(i+1) = kxr(i)
      end do

      do j = 1, ny
       if ( j.lt.ny/2+2 ) then
        kyr(j) = j - 1
       else
        kyr(j) = -(ny+1-j)
       endif
      end do

      do k = 1, nz
       if ( k.lt.nz/2+2 ) then
        kzr(k) = k - 1
       else
        kzr(k) = -(nz+1-k)
       endif
      end do

      iseedf = -iseedf
      ivf(:)  = 0
      iyf     = 0
      b1r = 0.0
      b2r = 0.0
      b3r = 0.0
      
      end subroutine para
!=============================================================================
!@subroutine allocarray
!@desc Allocates all globally used arrays in the var_inc module
!=============================================================================  
      subroutine allocarray
      use var_inc
      implicit none

      allocate (f(0:npop-1,lx,ly,lz,2))
      allocate (rho(lx,ly,lz))
      allocate (rhop(lx,ly,lz))
      allocate (ux(lx,ly,lz))
      allocate (uy(lx,ly,lz))
      allocate (uz(lx,ly,lz))
      allocate (ox(lx,ly,lz))
      allocate (oy(lx,ly,lz))
      allocate (oz(lx,ly,lz))
      allocate (kx(lx+2,ly+lyext,lz))
      allocate (ky(lx+2,ly+lyext,lz))
      allocate (kz(lx+2,ly+lyext,lz))
      allocate (k2(lx+2,ly+lyext,lz))
      allocate (ik2(lx+2,ly+lyext,lz))
      allocate (vx(lx+2,ly+lyext,lz))
      allocate (vy(lx+2,ly+lyext,lz))
      allocate (vz(lx+2,ly+lyext,lz))
      allocate (wx(lx+2,ly+lyext,lz))
      allocate (wy(lx+2,ly+lyext,lz))
      allocate (wz(lx+2,ly+lyext,lz))
      allocate (ibnodes(0:lx+1,0:ly+1,0:lz+1))

      allocate(force_realx(lx,ly,lz))
      allocate(force_realy(lx,ly,lz))
      allocate(force_realz(lx,ly,lz))

      ibnodes = -1

      if(ipart)then !If particles are being used
      
      allocate (ovlpflag(npart))
      allocate (ovlpflagt(npart))
      allocate (ipglb(msize))
      allocate (ypglb(3,npart))
      allocate (ypglbp(3,npart))
      allocate (yp(3,msize))
      allocate (thetap(3,msize))
      allocate (wp(3,npart))
      allocate (wpp(3,npart))
      allocate (omgp(3,npart))
      allocate (omgpp(3,npart))
      allocate (dwdt(3,msize))
      allocate (domgdt(3,msize))
      allocate (fHIp(3,npart))
      allocate (flubp(3,npart))
      allocate (forcep(3,msize))
      allocate (forcepp(3,msize))
      allocate (torqp(3,npart))
      allocate (torqpp(3,msize))

      allocate (xlink(maxlink))
      allocate (ylink(maxlink))
      allocate (zlink(maxlink))
      allocate (iplink(maxlink))
      allocate (mlink(maxlink))
      allocate (alink(maxlink))
      allocate (iblinks(0:2,2,maxlink))
      allocate (ipf_mym(2*19*lx*(lz + 4)))
      allocate (ipf_myp(2*19*lx*(lz + 4)))
      allocate (ipf_mzm(2*19*lx*ly))
      allocate (ipf_mzp(2*19*lx*ly))

      allocate (ibnodes0(0:lx+1,0:ly+1,0:lz+1))
      allocate (isnodes(lx,ly,lz))
      allocate (isnodes0(lx,ly,lz))

      allocate (xbfill(maxbfill))
      allocate (ybfill(maxbfill))
      allocate (zbfill(maxbfill))
      allocate (idbfill(maxbfill))
      allocate (fillMPIrequest(0:nproc,8))
      allocate (localReqData(8))

      allocate(fillRecvYm(0:npop-1,lx,0:lz+1))
      allocate(fillRecvYp(0:npop-1,lx,0:lz+1))
      allocate(fillRecvZm(0:npop-1,lx,ly))
      allocate(fillRecvZp(0:npop-1,lx,ly))

      isnodes = -1

      end if

!bench marking
      allocate (collision_MRT_bnch(nsteps))
      allocate (streaming_bnch(nsteps))
      allocate (macrovar_bnch(nsteps))

      allocate (beads_collision_bnch(nsteps))
      allocate (beads_lubforce_bnch(nsteps))
      allocate (beads_move_bnch(nsteps))
      allocate (beads_redistribute_bnch(nsteps))
      allocate (beads_links_bnch(nsteps))
      allocate (beads_filling_bnch(nsteps))

      end subroutine allocarray
!=============================================================================
!@subroutine constructMPItypes
!@desc Creates any custom MPI datatypes used
!=============================================================================  
      subroutine constructMPItypes
      use mpi
      use var_inc
      implicit none

      call MPI_TYPE_CONTIGUOUS(4, MPI_INTEGER, MPI_IPF_NODE, ierr)
      call MPI_TYPE_COMMIT(MPI_IPF_NODE, ierr)

      end subroutine constructMPItypes
!=============================================================================
!@subroutine makedir
!@desc Makes sure that the given directory exists, if it does not it attempts
!      to create it.
!@param dirPath = character array which holds the directory path that needs
!                 to be checked
!=============================================================================  
      subroutine makedir (dirPath)
      use var_inc

      implicit none
      logical dirExist
      character(len=120):: dirPath !Needs to be the same as declared length!

      !inquire( file=trim(DirPath)//'/.', exist=dirExist )  ! Works with gfortran
      inquire(directory = trim(dirPath), exist = dirExist)  ! Works with ifort (yellowstone)
      if(.NOT.dirExist)then
        write(*,*) trim(dirPath)//' not found. Creating...'
        call system('mkdir -p '// trim(dirPath)) !Execute system command to create directory
      endif

      end subroutine makedir
