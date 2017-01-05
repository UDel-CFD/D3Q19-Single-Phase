!=============================================================================
!University of Delaware LBM D3Q19 Single Phase Simulation 
!Copyright (C) 2017 Lian-Ping Wang

!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.

!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

      !Liscense Text
      if(myid==0)then
        write(*,*)'================'
        write(*,'(2X, A79)')'University of Delaware LBM D3Q19 Channel Flow Copyright (C) 2017 Lian-Ping Wang'
        write(*,'(2X, A81)')'This program comes with ABSOLUTELY NO WARRANTY; See LICENSE.txt file for details.'
        write(*,'(2X, A86)')'This is free software, and you are welcome to redistribute it under certain conditions'
        write(*,*)'================'
      endif

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
          istep0=0
          istep = istep0
          !For testing only
          call statistc !@file saveload.f90
          call statistc2 !@file saveload.f90
        end if  

      ELSE !Load an existing flow
        if(myid == 0)write(*,*)'Loading flow...'
        call loadcntdflow !@file saveload.f90
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
      
      !Call macrovar to initialize velocity variables
      call macrovar !@file collision.f90
      time_start = MPI_WTIME()

!=======================================================
! MAIN LOOP
!=======================================================
      do istep = istep0+1,istep0+nsteps 

        if(myid.eq.0 .and. mod(istep,50).eq.0)then
        	write(*,*) 'istep=',istep
        endif

        !Update or shut off perturb forcing, not used in particle laden
        !if(istep .gt. npforcing)then
        !else if(istep .lt. npforcing)then
        ! call FORCINGP
        !else if(istep .eq. npforcing)then
        ! call FORCING
        !end if

        !@file collision.f90
        call collision_MRT
        
        !Calculate macroscopic variables
        !@file collision.f90
        call macrovar

        if(ipart .and. mod(istep,100) == 0)then
          !Remove average density to correct mass loss from interpolation
          !@file collision.f90
          call avedensity
        endif


        !diag calculates and outputs respective data/diagnostics
        if(mod(istep,ndiag) == 0)  call diag !@file saveload.f90
        !partstatis is to output the position, velocity and force of each particle
        !if(ipart .and. mod(istep,200) == 0) call partstatis
        ! statistc4 outputs the mean, rms velocity&vorticity of particle phase, the vorticity calculation is called here
        !if(mod(istep,200) == 0)  call statistc4
        !  moviedata2 generates 2D vorticity contours for 2D visualization (only use when necessary)
        ! if(mod(istep,200) == 0) call moviedata2

        ! statistc2 calculates the mean and rms velocity profiles of fluid phase
        !if(mod(istep,nstat) == 0)  call statistc2 !@file saveload.f90
        ! statistc3 calculates the mean and rms vorticity profiles of fluid phase
        !if(mod(istep,nstat) == 0) call statistc3

        if(mod(istep,nflowout) == 0) call outputflow !@file saveload.f90

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
         endif !Line 356

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

101   continue

      if(myid.eq.0)then
        write(*,*)'time_max = ',time_max
      endif

      call MPI_FINALIZE(ierr)

      END PROGRAM main
