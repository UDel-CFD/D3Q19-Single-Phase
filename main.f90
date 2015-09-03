!=============================================================================
!       MPI F90 Code for Simulating 3D Decaying Homogeneous Isotropic 
!       Turbulence (DHIT) with finite-size freely-moving particles 
!       embedded in the cube.
!
!       Using LBM D3Q19 MRT Model.
!
!       MPI version is designed to operate on NCAR's Yellowstone
!
!       The DHIT code was originally written by Dr. Lian-Ping Wang, June 2000.
!
!       The particle part of the code was implemented by Hui Gao, Jan. 2010.
!       The original DHIT code was re-structured as well, using Dr. Yan Peng's
!       code as a reference.
!
!           x   y   z
!       {   0   0   0   }      rest
!
!       {   +   0   0   }      dir  1
!       {   -   0   0   }      dir  2
!       {   0   +   0   }      dir  3
!       {   0   -   0   }      dir  4
!       {   0   0   +   }      dir  5
!       {   0   0   -   }      dir  6
!
!       {   +   +   0   }      dir  7
!       {   -   +   0   }      dir  8
!       {   +   -   0   }      dir  9
!       {   -   -   0   }      dir  10
!       {   +   0   +   }      dir  11
!       {   -   0   +   }      dir  12
!       {   +   0   -   }      dir  13
!       {   -   0   -   }      dir  14
!       {   0   +   +   }      dir  15
!       {   0   -   +   }      dir  16
!       {   0   +   -   }      dir  17
!       {   0   -   -   }      dir  18
!=============================================================================
      PROGRAM main
      use mpi
      use var_inc
      implicit none

      integer:: i,j,k

      !Initialize MPI library and get rank/ topology size
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)      

      !Initialize all variables used
      !@file para.f90
      call para
      !Allocate all global arrays used in var_inc module
      !@file para.f90
      call allocarray
      !Construct custom MPI data types
      !@file para.f90
      call constructMPItypes

!=======================================================
! FLOW INITIALIZATION/ LOADING
!=======================================================
      IF(newrun)THEN
        if(newinitflow)then  

          !call initrand(iflowseed)
          ! Initialize velocity feild for the flow
          !@file intial.f90
          call initvel
          ! Calculate forcing values used in collision_MRT
          !@file collision.f90
          call FORCING

          ux = 0.0
          uy = 0.0
          uz = 0.0

          !Initialize particle populations on each node
          !@file partlib.f90
          call initpop

          istep = 0

          ! Pre-relaxation of density field after initial forcing
          do
            rhop = rho
            !Update density values
            !@file collision.f90
            call rhoupdat

            call collision_MRT !@file collision.f90

            !Determine the maximum density change
            rhoerr = maxval(abs(rho - rhop))        
            call MPI_ALLREDUCE(rhoerr,rhoerrmax,1,MPI_REAL8, &
                        MPI_MAX,MPI_COMM_WORLD,ierr)
            if(myid == 0 .and. mod(istep,1) == 0) write(*,*)istep, rhoerrmax

            !If density change is smaller than tolerance or loop has executed enough, break
            if(rhoerrmax <= rhoepsl .or. istep > 15000)then
              if(myid == 0) write(*,*)'final relaxation => ',istep, rhoerrmax
              exit
            end if
            istep = istep + 1 
          end do

          ! Check the maximum density fluctuations relative to rho0 = 1.0
          rhoerr = maxval(rho)        
          call MPI_ALLREDUCE(rhoerr, rhomax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
          rhoerr = minval(rho)        
          call MPI_ALLREDUCE(rhoerr, rhomin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
          if(myid == 0 ) write(*,*)istep, rhomax, rhomin


          ! save initial pre-relaxed flow      
          call saveinitflow !@file saveload.f90
          call macrovar !@file collison.f90
          istep0 = 0
          istep = istep0
          !For testing only
          call statistc !@file saveload.f90
          call statistc2 !@file saveload.f90

        else ! Load initial pre-relaxed flow
          call loadinitflow !@file saveload.f90
          call macrovar !@file collison.f90
          istep0=0
          istep = istep0
          !For testing only
          call statistc !@file saveload.f90
          call statistc2 !@file saveload.f90
        end if  

      ELSE !Load an existing flow
        if(myid == 0)write(*,*)'Loading flow...'
        call loadcntdflow !@file saveload.f90

        if(ipart .and. istpload > irelease)then
          call loadcntdpart  !@file saveload.f90
          call beads_links !@file partlib.f90
          released = .TRUE.
        end if

        call macrovar !@file collision.f90
      END IF
!=======================================================
! END OF FLOW INITIALIZATION/ LOADING
!=======================================================

      ! Calculate forcing values used in collision_MRT
      ! Forcing is independent of time so only call once here
      !@file collision.f90
      call FORCING
!     call FORCINGP
      time_start = MPI_WTIME()

!=======================================================
! MAIN LOOP
!=======================================================
      do istep = istep0+1,istep0+nsteps 

        if(myid.eq.0 .and. mod(istep,50).eq.0)then
        	write(*,*) 'istep=',istep
        endif

        ! Release partilces only after proper skewness (~ -0.5) has been established
        ! Initialise particle center position and velocity
        if(ipart .and. istep == irelease)then
          istep00 = 0
          if(myid.eq.0) write(*,*) 'Releasing beads'
          
          call initpartpos !@file partlib.f90
          call initpartvel !@file partlib.f90
          call initpartomg !@file partlib.f90
          call beads_links !@file partlib.f90

          istep00 = 1
          released = .TRUE.
        end if

        !Executes Collision and Propagation of the Fluid
        !@file collision.f90
        bnchstart = MPI_WTIME()
        call collision_MRT
        collision_MRT_bnch(istep-istpload) = MPI_WTIME() - bnchstart


        if(ipart .and. istep >= irelease)then !If we have solid particles in our flow
          !Calculates interpolation bounce back for the fluid off the solid particles
          !@file partlib.f90
          bnchstart = MPI_WTIME()
          call beads_collision
          beads_collision_bnch(istep-istpload) = MPI_WTIME() - bnchstart

          !Determine the lubrication forces acting on the solid particles
          !@file partlib.f90
          bnchstart = MPI_WTIME()
          call beads_lubforce
          beads_lubforce_bnch(istep-istpload) = MPI_WTIME() - bnchstart

          !Update the position of the solid particles
          !@file partlib.f90
          bnchstart = MPI_WTIME()
          call beads_move
          beads_move_bnch(istep-istpload) = MPI_WTIME() - bnchstart

          !Redistribute solid particles to their respective MPI task based on their global position
          !@file partlib.f90
          bnchstart = MPI_WTIME()
          call beads_redistribute
          beads_redistribute_bnch(istep-istpload) = MPI_WTIME() - bnchstart

          !Determine which lattice links cross the fluid/ solid particle boundary,
          !the position the boundary is located on the link, which nodes are inside a solid particle,
          !and what data we will need from neighboring MPI tasks for interpolation bounce back
          !@file partlib.f90
          bnchstart = MPI_WTIME()
          call beads_links
          beads_links_bnch(istep-istpload) = MPI_WTIME() - bnchstart

          !Repopulate nodes that previously existed inside a solid particle
          !@file partlib.f90
          bnchstart = MPI_WTIME()
          call beads_filling
          beads_filling_bnch(istep-istpload) = MPI_WTIME() - bnchstart
        end if

        !Calculate macroscopic variables
        !@file collision.f90
        bnchstart = MPI_WTIME()
        call macrovar
        macrovar_bnch(istep-istpload) = MPI_WTIME() - bnchstart

        !Calculate and output respective data/diagnostics
        if(mod(istep,ndiag) == 0)  call diag !@file saveload.f90
        if(mod(istep,nstat) == 0)  call statistc !@file saveload.f90
        if(mod(istep,nstat) == 0)  call statistc2 !@file saveload.f90
        if(mod(istep,nflowout) == 0) call outputflow !@file saveload.f90

        if(ipart .and. istep >= irelease .and. mod(istep,npartout) == 0)call outputpart !@file saveload.f90

!       if(ipart .and. istep >= irelease .and. mod(istep,nmovieout) == 0) then
!          call moviedata
!          call sijstat03
!          go to 101
!       end if

!       if(mod(istep,nstat) == 0) call rmsstat
!       if(mod(istep,nsij) == 0) call sijstat03   
!       if(mod(istep,nsij) == 0) call sijstat 

        !Check current wall-clock time
        if(mod(istep,ntime) == 0)then
           !Determine current runtime and collect all runtimes from MPI tasks
           time_end = MPI_WTIME()
           time_diff = time_end - time_start

           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
           call MPI_ALLREDUCE(time_diff, time_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
          ! If our current runtime exceeds our specified limit, stop and save the run
          if(time_max > time_bond) exit
        endif

      end do
!=======================================================
! END OF MAIN LOOP
!=======================================================

      !Get main loop wall-clock time 
      time_end = MPI_WTIME()
      time_diff = time_end - time_start

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(time_diff,time_max,1,MPI_REAL8, MPI_MAX,MPI_COMM_WORLD,ierr)

!Probe each processor for final flow comparing
       call probe
       !call outputvort
! save data for continue run
      !call savecntdflow
!save param variables
!      call input_outputf(2)
!save bead positions
      !if(ipart .and. istep > irelease) call savecntdpart    

!Record Benchmarks
      call benchflow
      call benchbead
      call benchmatlab

      call benchtotal

101   continue

      if(myid.eq.0)then
      write(*,*)'time_max = ',time_max
!     write(50,*)' ',time_stream_max_array
      endif

      call MPI_FINALIZE(ierr)

      END PROGRAM main
