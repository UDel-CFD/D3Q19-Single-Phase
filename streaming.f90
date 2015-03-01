!==================================================================
      subroutine streaming
      use mpi
      use var_inc
      implicit none
 
      real, dimension(lx,ly)     :: tmpz1,tmpz2
      real, dimension(lx,lz)     :: tmpy1,tmpy2
      real, dimension(0:npop-1,ly,lz)     :: tmpxL,tmpxR

!     Save two layers to be used after streaming
      tmpxL = f(:,1,:,:)
      tmpxR = f(:,lx,:,:)

! x+ move Dir 1,7,9,11,13  

      time1in = MPI_WTIME()
      f(1,:,:,:) = cshift(f(1,:,:,:), shift = -1, dim = 1)
      f(7,:,:,:) = cshift(f(7,:,:,:), shift = -1, dim = 1)
      f(9,:,:,:) = cshift(f(9,:,:,:), shift = -1, dim = 1)
      f(11,:,:,:) = cshift(f(11,:,:,:), shift = -1, dim = 1)
      f(13,:,:,:) = cshift(f(13,:,:,:), shift = -1, dim = 1)
      time2in = MPI_WTIME()
      time_stXp = time_stXp + (time2in-time1in)
    
! x- move Dir 2,8,10,12,14

      time1in = MPI_WTIME()
      f(2,:,:,:) = cshift(f(2,:,:,:), shift = 1, dim = 1)
      f(8,:,:,:) = cshift(f(8,:,:,:), shift = 1, dim = 1)
      f(10,:,:,:) = cshift(f(10,:,:,:), shift = 1, dim = 1)
      f(12,:,:,:) = cshift(f(12,:,:,:), shift = 1, dim = 1)
      f(14,:,:,:) = cshift(f(14,:,:,:), shift = 1, dim = 1)
      time2in = MPI_WTIME()
      time_stXm = time_stXm + (time2in-time1in)

! y+ move Dir 3,7,8,15,17  

      time1in = MPI_WTIME()

      tmpy1 = f(3,:,ly,:)
      call exchng0y(tmpy1,tmpy2,3,1)
      f(3,:,:,:) = cshift(f(3,:,:,:), shift = -1, dim = 2)
      f(3,:,1,:) = tmpy2

      tmpy1 = f(7,:,ly,:)
      call exchng0y(tmpy1,tmpy2,7,1)
      f(7,:,:,:) = cshift(f(7,:,:,:), shift = -1, dim = 2)
      f(7,:,1,:) = tmpy2

      tmpy1 = f(8,:,ly,:)
      call exchng0y(tmpy1,tmpy2,8,1)
      f(8,:,:,:) = cshift(f(8,:,:,:), shift = -1, dim = 2)
      f(8,:,1,:) = tmpy2

      tmpy1 = f(15,:,ly,:)
      call exchng0y(tmpy1,tmpy2,15,1)
      f(15,:,:,:) = cshift(f(15,:,:,:), shift = -1, dim = 2)
      f(15,:,1,:) = tmpy2

      tmpy1 = f(17,:,ly,:)
      call exchng0y(tmpy1,tmpy2,17,1)
      f(17,:,:,:) = cshift(f(17,:,:,:), shift = -1, dim = 2)
      f(17,:,1,:) = tmpy2

      time2in = MPI_WTIME()
      time_stYp = time_stYp + (time2in-time1in)

! y- move Dir 4,9,10,16,18  

      time1in = MPI_WTIME()

      tmpy1 = f(4,:,1,:)
      call exchng0y(tmpy1,tmpy2,4,-1)
      f(4,:,:,:) = cshift(f(4,:,:,:), shift = 1, dim = 2)
      f(4,:,ly,:) = tmpy2


      tmpy1 = f(9,:,1,:)
      call exchng0y(tmpy1,tmpy2,9,-1)
      f(9,:,:,:) = cshift(f(9,:,:,:), shift = 1, dim = 2)
      f(9,:,ly,:) = tmpy2


      tmpy1 = f(10,:,1,:)
      call exchng0y(tmpy1,tmpy2,10,-1)
      f(10,:,:,:) = cshift(f(10,:,:,:), shift = 1, dim = 2)
      f(10,:,ly,:) = tmpy2


      tmpy1 = f(16,:,1,:)
      call exchng0y(tmpy1,tmpy2,16,-1)
      f(16,:,:,:) = cshift(f(16,:,:,:), shift = 1, dim = 2)
      f(16,:,ly,:) = tmpy2

      tmpy1 = f(18,:,1,:)
      call exchng0y(tmpy1,tmpy2,18,-1)
      f(18,:,:,:) = cshift(f(18,:,:,:), shift = 1, dim = 2)
      f(18,:,ly,:) = tmpy2

      time2in = MPI_WTIME()
      time_stYm = time_stYm + (time2in-time1in)

! z+ move Dir 5,11,12,15,16   

      time1in = MPI_WTIME()

      tmpz1 = f(5,:,:,lz)   
      call exchng0z(tmpz1,tmpz2,5,1)
      f(5,:,:,:) = cshift(f(5,:,:,:), shift = -1, dim = 3)
      f(5,:,:,1) = tmpz2


      tmpz1 = f(11,:,:,lz)   
      call exchng0z(tmpz1,tmpz2,11,1)
      f(11,:,:,:) = cshift(f(11,:,:,:), shift = -1, dim = 3)
      f(11,:,:,1) = tmpz2


      tmpz1 = f(12,:,:,lz)   
      call exchng0z(tmpz1,tmpz2,12,1)
      f(12,:,:,:) = cshift(f(12,:,:,:), shift = -1, dim = 3)
      f(12,:,:,1) = tmpz2


      tmpz1 = f(15,:,:,lz)   
      call exchng0z(tmpz1,tmpz2,15,1)
      f(15,:,:,:) = cshift(f(15,:,:,:), shift = -1, dim = 3)
      f(15,:,:,1) = tmpz2


      tmpz1 = f(16,:,:,lz)   
      call exchng0z(tmpz1,tmpz2,16,1)
      f(16,:,:,:) = cshift(f(16,:,:,:), shift = -1, dim = 3)
      f(16,:,:,1) = tmpz2

      time2in = MPI_WTIME()
      time_stZp = time_stZp + (time2in-time1in)

! z- move Dir: 6,13,14,17,18   

      time1in = MPI_WTIME()

      tmpz1 = f(6,:,:,1)   
      call exchng0z(tmpz1,tmpz2,6,-1)
      f(6,:,:,:) = cshift(f(6,:,:,:), shift = 1, dim = 3)
      f(6,:,:,lz) = tmpz2


      tmpz1 = f(13,:,:,1)   
      call exchng0z(tmpz1,tmpz2,13,-1)
      f(13,:,:,:) = cshift(f(13,:,:,:), shift = 1, dim = 3)
      f(13,:,:,lz) = tmpz2


      tmpz1 = f(14,:,:,1)   
      call exchng0z(tmpz1,tmpz2,14,-1)
      f(14,:,:,:) = cshift(f(14,:,:,:), shift = 1, dim = 3)
      f(14,:,:,lz) = tmpz2


      tmpz1 = f(17,:,:,1)   
      call exchng0z(tmpz1,tmpz2,17,-1)
      f(17,:,:,:) = cshift(f(17,:,:,:), shift = 1, dim = 3)
      f(17,:,:,lz) = tmpz2


      tmpz1 = f(18,:,:,1)   
      call exchng0z(tmpz1,tmpz2,18,-1)
      f(18,:,:,:) = cshift(f(18,:,:,:), shift = 1, dim = 3)
      f(18,:,:,lz) = tmpz2

      time2in = MPI_WTIME()
      time_stZm = time_stZm + (time2in-time1in)
! middle-link bounce back
! x+ move Dir 1,7,9,11,13  
      f( 1,1,:,:)=tmpxL( 2,:,:)
      f( 7,1,:,:)=tmpxL(10,:,:)
      f( 9,1,:,:)=tmpxL( 8,:,:)
      f(11,1,:,:)=tmpxL(14,:,:)
      f(13,1,:,:)=tmpxL(12,:,:)
! x- move Dir 2,8,10,12,14
      f( 2,lx,:,:)=tmpxR( 1,:,:)
      f(10,lx,:,:)=tmpxR( 7,:,:)
      f( 8,lx,:,:)=tmpxR( 9,:,:)
      f(14,lx,:,:)=tmpxR(11,:,:)
      f(12,lx,:,:)=tmpxR(13,:,:)

      end subroutine streaming
!==================================================================

      subroutine exchng0z(tmpz1,tmpz2,tag,zdir)
      use mpi
      use var_inc
      implicit none
      
      real, dimension(lx,ly)     :: tmpz1,tmpz2

      integer tag, zdir, ilen, ips, ipr
      integer status_array(MPI_STATUS_SIZE,2), req(2)

      ilen = lx*ly

      if(zdir > 0)then
        ips = mzp
        ipr = mzm
      else
        ips = mzm
        ipr = mzp
      end if

      call MPI_ISEND(tmpz1,ilen,MPI_REAL8,ips,tag,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmpz2,ilen,MPI_REAL8,ipr,tag,MPI_COMM_WORLD,req(2),ierr)
      call MPI_WAITALL(2,req,status_array,ierr)


      end subroutine exchng0z 
!==================================================================

      subroutine exchng0y(tmpy1,tmpy2,tag,ydir)
      use mpi
      use var_inc
      implicit none

      real, dimension(lx,lz)     :: tmpy1,tmpy2

      integer tag, ydir, ilen, ips, ipr
      integer status_array(MPI_STATUS_SIZE,2), req(2)

      ilen = lx*lz

      if(ydir > 0)then
        ips = myp
        ipr = mym
      else
        ips = mym
        ipr = myp
      end if

      call MPI_ISEND(tmpy1,ilen,MPI_REAL8,ips,tag,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmpy2,ilen,MPI_REAL8,ipr,tag,MPI_COMM_WORLD,req(2),ierr)
      call MPI_WAITALL(2,req,status_array,ierr)


      end subroutine exchng0y
!==================================================================



      
