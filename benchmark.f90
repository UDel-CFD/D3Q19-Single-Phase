!==========================================================================================
!==========================================================================================
     subroutine benchflow
     use var_inc

     implicit none
     integer iprc1, iprc2, iprc3, i
     character (len = 150):: fnm

     iprc1 = myid / 100
     iprc2 = mod(myid,100) / 10
     iprc3 = mod(myid,10)

     fnm = trim(dirbenchflow)//'bmflow.' &
     				//char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

     open(60, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

     write(60,500) 'collision_MRT//streaming//macrovar'
     do i = 1, nsteps
     	write(60,600) collision_MRT_bnch(i), streaming_bnch(i), macrovar_bnch(i)
     enddo
     close(60)

500   format(2x,1A34)
600   format(2x,3(1f16.10))
     end subroutine benchflow
!==========================================================================================
!==========================================================================================
     subroutine benchbead
     use var_inc

     implicit none
     integer iprc1, iprc2, iprc3, i
     character (len = 150):: fnm

     iprc1 = myid / 100
     iprc2 = mod(myid,100) / 10
     iprc3 = mod(myid,10)

     fnm = trim(dirbenchbead)//'bmbead.' &
     				//char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

     open(60, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

     write(60,500) 'beads_collision//beads_lubforce//beads_move' &
        				//'//beads_redistribute//beads_links//beads_filling'
     do i = 1, nsteps
     	write(60,600) beads_collision_bnch(i), beads_lubforce_bnch(i), &
                beads_move_bnch(i), beads_redistribute_bnch(i), beads_links_bnch(i),&
                beads_filling_bnch(i)
     enddo
     close(60)

500   format(2x,1A91)
600   format(2x,6(1f16.10))
     end subroutine benchbead
!==========================================================================================
!==========================================================================================
     subroutine benchmatlab
     use var_inc

     implicit none
     integer iprc1, iprc2, iprc3, i
     character (len = 150):: fnm

     iprc1 = myid / 100
     iprc2 = mod(myid,100) / 10
     iprc3 = mod(myid,10)

     fnm = trim(dirbenchmatlab)//'bmmatlab.' &
     				//char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

     open(60, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

     do i = 1, nsteps
     	write(60,600) collision_MRT_bnch(i), beads_collision_bnch(i), streaming_bnch(i), &
 	    			beads_lubforce_bnch(i), beads_move_bnch(i), beads_redistribute_bnch(i), &
 	    			beads_links_bnch(i), beads_filling_bnch(i), macrovar_bnch(i)
     enddo
     close(60)

600   format(2x,9(1f16.10))
     end subroutine benchmatlab
!==========================================================================================
!==========================================================================================
     subroutine benchtotal
     use mpi
     use var_inc

     implicit none
     integer id
     character (len = 120):: fnm
     character (len = 25):: filenum
     real, dimension(nproc) :: totalbnch

     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     CALL MPI_GATHER(time_diff, 1, MPI_REAL8, totalbnch, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

     write (filenum, "(I10)") nsteps
     if(myid == 0)then !print ut probe data to file
        fnm = trim(dirbench)//'runtimes'//trim(adjustl(filenum))//'.dat'

        open(60, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        do id = 1,nproc
          write(60,600) id-1, totalbnch(id)
        enddo
      endif
      close(60)

500   format(2x,1A29)
600   format(2x,1i3,1f16.10)
     end subroutine benchtotal
!==========================================================================================
!==========================================================================================
     subroutine benchbeads_collision
     use var_inc

     implicit none
     integer iprc1, iprc2, iprc3, i
     character (len = 150):: fnm
     character (len = 10):: filenum

     write (filenum, "(I10)") myid

     fnm = trim(dirbench)//'beads_collision/bench'//trim(adjustl(filenum))//'.dat'

     open(60, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

     do i = 1, nsteps
        !write(60,600) beads_collision_ex2(i), beads_collision_loop(i), beads_collision_allre(i)
     enddo
     close(60)

600   format(2x,9(1f16.10))
     end subroutine benchbeads_collision
!==========================================================================================
!==========================================================================================