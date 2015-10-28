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

            call collision_MRT !@file collision.f90

            call rhoupdat
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
          if(myid.eq.0) write(*,*)'start to save initial flow'     
          call saveinitflow !@file saveload.f90
          if(myid.eq.0) write(*,*)'save successfully'
!          call macrovar !@file collison.f90
          istep0 = 0
          istep = istep0
          !For testing only
!          call statistc !@file saveload.f90
!          call statistc2 !@file saveload.f90

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
!          if(myid.eq.0)write(*,*)'Initial beads_link passed'
          released = .TRUE.
        end if

        call macrovar !@file collision.f90
      END IF

      if(myid.eq.0)write(*,*)'Loading successful'

!=======================================================
! END OF FLOW INITIALIZATION/ LOADING
!=======================================================

      ! Calculate forcing values used in collision_MRT
      ! Forcing is independent of time so only call once here
      !@file collision.f90
      call FORCING
!     call FORCINGP
      time_start = MPI_WTIME()

      !Create Pipe boundary
      !@file pipelib.f90
      call pipe_links
!=======================================================
      ! Simplified version of time measurement
      ! Should be commented out in regular runs

       time_collcum = 0.d0
       time_streacum = 0.d0
       time_macrocum = 0.d0
       time_bcolcum = 0.d0
       time_blubcum = 0.d0
       time_bmovcum = 0.d0
       time_bredcum = 0.d0
       time_blincum = 0.d0
       time_bfilcum = 0.d0


!=======================================================
! MAIN LOOP
!=======================================================
       if(myid.eq.0) write(*,*)'Enter main loop'
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

        time_start = MPI_WTIME()
        !Executes Collision and Propagation of the Fluid
        !@file collision.f90
        !bnchstart = MPI_WTIME()
        time1 = MPI_WTIME()
        call collision_MRT
        time2 = MPI_WTIME()
        time_coll = time2 - time1
        !collision_MRT_bnch(istep-istpload) = MPI_WTIME() - bnchstart
        
!        if (myid.eq.0)write(*,*)'Passed collision and streaming'
        if(ipart .and. istep >= irelease)then !If we have solid particles in our flow
          !Calculates interpolation bounce back for the fluid off the solid particles
          !@file partlib.f90
          !bnchstart = MPI_WTIME()
          time1 = MPI_WTIME()
          call beads_collision
          time2 = MPI_WTIME()
          time_bcol = time2 - time1
          !beads_collision_bnch(istep-istpload) = MPI_WTIME() - bnchstart

          !Determine the lubrication forces acting on the solid particles
          !@file partlib.f90
          !bnchstart = MPI_WTIME()
          time1 = MPI_WTIME()
          call beads_lubforce
          time2 = MPI_WTIME()
          time_blub = time2 - time1
          !beads_lubforce_bnch(istep-istpload) = MPI_WTIME() - bnchstart

          !Update the position of the solid particles
          !@file partlib.f90
          !bnchstart = MPI_WTIME()
          time1 = MPI_WTIME()
          call beads_move
          time2 = MPI_WTIME()
          time_bmov = time2 - time1
          !beads_move_bnch(istep-istpload) = MPI_WTIME() - bnchstart

          !Redistribute solid particles to their respective MPI task based on their global position
          !@file partlib.f90
          !bnchstart = MPI_WTIME()
          time1 = MPI_WTIME()
          call beads_redistribute
          time2 = MPI_WTIME()
          time_bred = time2 - time1
          !beads_redistribute_bnch(istep-istpload) = MPI_WTIME() - bnchstart

          !Determine which lattice links cross the fluid/ solid particle boundary,
          !the position the boundary is located on the link, which nodes are inside a solid particle,
          !and what data we will need from neighboring MPI tasks for interpolation bounce back
          !@file partlib.f90
          !bnchstart = MPI_WTIME()
          time1 = MPI_WTIME()
          call beads_links
          time2 = MPI_WTIME()
          time_blin = time2 - time1
          !beads_links_bnch(istep-istpload) = MPI_WTIME() - bnchstart

          !Repopulate nodes that previously existed inside a solid particle
          !@file partlib.f90
          !bnchstart = MPI_WTIME()
          time1 = MPI_WTIME()
          call beads_filling
          time2 = MPI_WTIME()
          time_bfil = time2 - time1
          !beads_filling_bnch(istep-istpload) = MPI_WTIME() - bnchstart
!          if(myid.eq.0)write(*,*)'Passed beads_filling'
        end if

          call beads_collision
        !Calculate macroscopic variables
        !@file collision.f90
        !bnchstart = MPI_WTIME()
        time1 = MPI_WTIME()
        call macrovar
        time2 = MPI_WTIME()
        time_macro = time2 - time1
        time_end = MPI_WTIME()
        time_diff = time_end - time_start
        !macrovar_bnch(istep-istpload) = MPI_WTIME() - bnchstart

        if(mod(istep,100) == 0) then
          call avedensity
          rhomeang = rhomeang/dfloat(nfluidtotal)
          call cordensity
        endif


! Get the maximum time in all processors
        call MPI_ALLREDUCE(time_diff,time_diffmax,1,MPI_REAL8,   &
               MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(time_coll,time_collmax,1,MPI_REAL8,   &
               MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(time_strea,time_streamax,1,MPI_REAL8,   &
               MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(time_bcol,time_bcolmax,1,MPI_REAL8,   &
               MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(time_blub,time_blubmax,1,MPI_REAL8,   &
               MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(time_bmov,time_bmovmax,1,MPI_REAL8,   &
               MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(time_bred,time_bredmax,1,MPI_REAL8,   &
               MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(time_blin,time_blinmax,1,MPI_REAL8,   &
               MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(time_bfil,time_bfilmax,1,MPI_REAL8,   &
               MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(time_macro,time_macromax,1,MPI_REAL8,   &
               MPI_MAX,MPI_COMM_WORLD,ierr)

        if (myid .eq. 0) then
         time_diffcum = time_diffcum + time_diffmax
         time_collcum = time_collcum + time_collmax
         time_streacum = time_streacum + time_streamax
         time_bcolcum = time_bcolcum + time_bcolmax
         time_blubcum = time_blubcum + time_blubmax
         time_bmovcum = time_bmovcum + time_bmovmax
         time_bredcum = time_bredcum + time_bredmax
         time_blincum = time_blincum + time_blinmax
         time_bfilcum = time_bfilcum + time_bfilmax
         time_macrocum = time_macrocum + time_macromax
        end if


        !diag calculates and outputs respective data/diagnostics
        if(mod(istep,ndiag) == 0)  call diag !@file saveload.f90
        !partstatis is to output the position, velocity
        !and force of each particle
!        if(mod(istep,200) == 0) call partstatis
        ! statistc4 outputs the mean, rms velocity&vorticity of particle phase, the vorticity calculation is called here
!        if(mod(istep,200) == 0)  call statistc4
        !  moviedata2 generates 2D vorticity contours for 2D visualization (only use when necessary)
        ! if(mod(istep,200) == 0) call moviedata2
        ! statistc2 calculates the mean and rms velocity profiles of fluid phase
        if(mod(istep,nstat) == 0)  call statistc2 !@file saveload.f90

        ! statistc3 calculates the mean and rms vorticity profiles of fluid phase
!        if(mod(istep,nstat) == 0) call statistc3

        if(mod(istep,nflowout) == 0) call outputuy !@file saveload.f90

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
        ! if(mod(istep,ntime) == 0)then
           !Determine current runtime and collect all runtimes from MPI tasks
        !   time_end = MPI_WTIME()
        !   time_diff = time_end - time_start

         !  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         !  call MPI_ALLREDUCE(time_diff, time_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
          ! If our current runtime exceeds our specified limit, stop and save the run
          if(time_max > time_bond) exit
        ! endif !Line 356

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
       !call probe
       !call outputvort
! save data for continue run
      call savecntdflow
!save param variables
!      call input_outputf(2)
!save bead positions
      if(ipart .and. istep > irelease) call savecntdpart    

!Record Benchmarks
!      call benchflow
!      call benchbead
!      call benchmatlab

!      call benchtotal

101   continue

      if(myid.eq.0)then
      write(*,*)'time_max = ',time_max
!     write(50,*)' ',time_stream_max_array
      write(*,*)'total time per step = ',time_diffcum/float(nsteps)
      write(*,*)'collision per step = ',time_collcum/float(nsteps)
      write(*,*)'streaming per step = ',time_streacum/float(nsteps)
      write(*,*)'bead_col per step = ',time_bcolcum/float(nsteps)
      write(*,*)'bead_lub per step = ',time_blubcum/float(nsteps)
      write(*,*)'bead_mov per step = ',time_bmovcum/float(nsteps)
      write(*,*)'bead_red per step = ',time_bredcum/float(nsteps)
      write(*,*)'bead_lin per step = ',time_blincum/float(nsteps)
      write(*,*)'bead_fil per step = ',time_bfilcum/float(nsteps)
      write(*,*)'marcovar per step = ',time_macrocum/float(nsteps)
      endif

      call MPI_FINALIZE(ierr)

      END PROGRAM main