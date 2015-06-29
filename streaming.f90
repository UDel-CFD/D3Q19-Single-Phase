!==================================================================
      subroutine streaming
      use mpi
      use var_inc
      implicit none
 
      real, dimension(lx,ly)     :: tmpz1,tmpz2
      real, dimension(lx,lz)     :: tmpy1,tmpy2
      real, dimension(ly,lz,0:npop-1)     :: tmpxL,tmpxR
      integer i, t1, t2

!     Save two layers to be used after streaming
      tmpxL = f(1,:,:,:)
      tmpxR = f(lx,:,:,:)

! x+ move Dir 1,7,9,11,13  

      time1in = MPI_WTIME()
      f(:,:,:,1) = cshift(f(:,:,:,1), shift = -1, dim = 1)
      f(:,:,:,7) = cshift(f(:,:,:,7), shift = -1, dim = 1)
      f(:,:,:,9) = cshift(f(:,:,:,9), shift = -1, dim = 1)
      f(:,:,:,11) = cshift(f(:,:,:,11), shift = -1, dim = 1)
      f(:,:,:,13) = cshift(f(:,:,:,13), shift = -1, dim = 1)
      time2in = MPI_WTIME()
      time_stXp = time_stXp + (time2in-time1in)
    
! x- move Dir 2,8,10,12,14

      time1in = MPI_WTIME()
      f(:,:,:,2) = cshift(f(:,:,:,2), shift = 1, dim = 1)
      f(:,:,:,8) = cshift(f(:,:,:,8), shift = 1, dim = 1)
      f(:,:,:,10) = cshift(f(:,:,:,10), shift = 1, dim = 1)
      f(:,:,:,12) = cshift(f(:,:,:,12), shift = 1, dim = 1)
      f(:,:,:,14) = cshift(f(:,:,:,14), shift = 1, dim = 1)
      time2in = MPI_WTIME()
      time_stXm = time_stXm + (time2in-time1in)

! y+ move Dir 3,7,8,15,17  
      time1in = MPI_WTIME()

      tmpy1 = f(:,ly,:,3)
      t1 = MPI_WTIME()
      call exchng0y(tmpy1,tmpy2,3,1)
      f(:,:,:,3) = cshift(f(:,:,:,3), shift = -1, dim = 2)
      f(:,1,:,3) = tmpy2(:,:)

      tmpy1 = f(:,ly,:,7)
      call exchng0y(tmpy1,tmpy2,7,1)
      f(:,:,:,7) = cshift(f(:,:,:,7), shift = -1, dim = 2)
      f(:,1,:,7) = tmpy2(:,:)

      tmpy1 = f(:,ly,:,8)
      call exchng0y(tmpy1,tmpy2,8,1)
      f(:,:,:,8) = cshift(f(:,:,:,8), shift = -1, dim = 2)
      f(:,1,:,8) = tmpy2(:,:)

      tmpy1 = f(:,ly,:,15)
      call exchng0y(tmpy1,tmpy2,15,1)
      f(:,:,:,15) = cshift(f(:,:,:,15), shift = -1, dim = 2)
      f(:,1,:,15) = tmpy2(:,:)

      tmpy1 = f(:,ly,:,17)
      call exchng0y(tmpy1,tmpy2,17,1)
      f(:,:,:,17) = cshift(f(:,:,:,17), shift = -1, dim = 2)
      f(:,1,:,17) = tmpy2(:,:)

      time2in = MPI_WTIME()
      time_stYp = time_stYp + (time2in-time1in)

! y- move Dir 4,9,10,16,18  
      time1in = MPI_WTIME()

      tmpy1 = f(:,1,:,4)
      call exchng0y(tmpy1,tmpy2,4,-1)
      f(:,:,:,4) = cshift(f(:,:,:,4), shift = 1, dim = 2)
      f(:,ly,:,4) = tmpy2(:,:)

      tmpy1 = f(:,1,:,9)
      call exchng0y(tmpy1,tmpy2,9,-1)
      f(:,:,:,9) = cshift(f(:,:,:,9), shift = 1, dim = 2)
      f(:,ly,:,9) = tmpy2(:,:)

      tmpy1 = f(:,1,:,10)
      call exchng0y(tmpy1,tmpy2,10,-1)
      f(:,:,:,10) = cshift(f(:,:,:,10), shift = 1, dim = 2)
      f(:,ly,:,10) = tmpy2(:,:)


      tmpy1 = f(:,1,:,16)
      call exchng0y(tmpy1,tmpy2,16,-1)
      f(:,:,:,16) = cshift(f(:,:,:,16), shift = 1, dim = 2)
      f(:,ly,:,16) = tmpy2(:,:)

      tmpy1 = f(:,1,:,18)
      call exchng0y(tmpy1,tmpy2,18,-1)
      f(:,:,:,18) = cshift(f(:,:,:,18), shift = 1, dim = 2)
      f(:,ly,:,18) = tmpy2(:,:)

      time2in = MPI_WTIME()
      time_stYm = time_stYm + (time2in-time1in)

! z+ move Dir 5,11,12,15,16   

      time1in = MPI_WTIME()

      tmpz1 = f(:,:,lz,5)   
      call exchng0z(tmpz1,tmpz2,5,1)
      f(:,:,:,5) = cshift(f(:,:,:,5), shift = -1, dim = 3)
      f(:,:,1,5) = tmpz2


      tmpz1 = f(:,:,lz,11)   
      call exchng0z(tmpz1,tmpz2,11,1)
      f(:,:,:,11) = cshift(f(:,:,:,11), shift = -1, dim = 3)
      f(:,:,1,11) = tmpz2


      tmpz1 = f(:,:,lz,12)   
      call exchng0z(tmpz1,tmpz2,12,1)
      f(:,:,:,12) = cshift(f(:,:,:,12), shift = -1, dim = 3)
      f(:,:,1,12) = tmpz2


      tmpz1 = f(:,:,lz,15)   
      call exchng0z(tmpz1,tmpz2,15,1)
      f(:,:,:,15) = cshift(f(:,:,:,15), shift = -1, dim = 3)
      f(:,:,1,15) = tmpz2


      tmpz1 = f(:,:,lz,16)   
      call exchng0z(tmpz1,tmpz2,16,1)
      f(:,:,:,16) = cshift(f(:,:,:,16), shift = -1, dim = 3)
      f(:,:,1,16) = tmpz2

      time2in = MPI_WTIME()
      time_stZp = time_stZp + (time2in-time1in)

! z- move Dir: 6,13,14,17,18   

      time1in = MPI_WTIME()

      tmpz1 = f(:,:,1,6)   
      call exchng0z(tmpz1,tmpz2,6,-1)
      f(:,:,:,6) = cshift(f(:,:,:,6), shift = 1, dim = 3)
      f(:,:,lz,6) = tmpz2


      tmpz1 = f(:,:,1,13)   
      call exchng0z(tmpz1,tmpz2,13,-1)
      f(:,:,:,13) = cshift(f(:,:,:,13), shift = 1, dim = 3)
      f(:,:,lz,13) = tmpz2


      tmpz1 = f(:,:,1,14)   
      call exchng0z(tmpz1,tmpz2,14,-1)
      f(:,:,:,14) = cshift(f(:,:,:,14), shift = 1, dim = 3)
      f(:,:,lz,14) = tmpz2


      tmpz1 = f(:,:,1,17)   
      call exchng0z(tmpz1,tmpz2,17,-1)
      f(:,:,:,17) = cshift(f(:,:,:,17), shift = 1, dim = 3)
      f(:,:,lz,17) = tmpz2


      tmpz1 = f(:,:,1,18)   
      call exchng0z(tmpz1,tmpz2,18,-1)
      f(:,:,:,18) = cshift(f(:,:,:,18), shift = 1, dim = 3)
      f(:,:,lz,18) = tmpz2

      time2in = MPI_WTIME()
      time_stZm = time_stZm + (time2in-time1in)
! middle-link bounce back
! x+ move Dir 1,7,9,11,13  
      f(1,:,:,1) =tmpxL(:,:,2)
      f(1,:,:,7) =tmpxL(:,:,10)
      f(1,:,:,9) =tmpxL(:,:,8)
      f(1,:,:,11)=tmpxL(:,:,14)
      f(1,:,:,13)=tmpxL(:,:,12)
! x- move Dir 2,8,10,12,14
      f(lx,:,:,2) =tmpxR(:,:,1)
      f(lx,:,:,10)=tmpxR(:,:,7)
      f(lx,:,:,8) =tmpxR(:,:,9)
      f(lx,:,:,14)=tmpxR(:,:,11)
      f(lx,:,:,12)=tmpxR(:,:,13)

      if(ipart)then
        !Inject fluid solid particle collisions into fluid domain
        do i=1, ifsc_inject
          f(fsc_inject(i)%x,fsc_inject(i)%y,fsc_inject(i)%z,fsc_inject(i)%ip) = fsc_inject(i)%dist
        enddo
        deallocate(fsc_inject)
      endif

      end subroutine streaming
!==================================================================
      subroutine streamingp
      use mpi
      use var_inc
      implicit none
 
      real, dimension(5,lx,ly)     :: tmpzpS,tmpzpR,tmpzmS,tmpzmR
      real, dimension(5,lx,lz)     :: tmpypS,tmpypR,tmpymS,tmpymR
      real, dimension(0:npop-1,ly,lz)     :: tmpxL,tmpxR
      integer status(MPI_STATUS_SIZE), status_array(MPI_STATUS_SIZE,4), reqs(4), reqr(2)
      integer i, index, ylen, zlen

      ylen = 5*lx*lz
      zlen = 5*lx*ly
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

! Create y+ stream buffer Dir 3,7,8,15,17
      tmpypS(1,:,:) = f(3,:,ly,:)
      tmpypS(2,:,:) = f(7,:,ly,:)
      tmpypS(3,:,:) = f(8,:,ly,:)
      tmpypS(4,:,:) = f(15,:,ly,:)
      tmpypS(5,:,:) = f(17,:,ly,:)
! Create Y- stream buffer Dir 4,9,10,16,18
      tmpymS(1,:,:) = f(4,:,1,:)
      tmpymS(2,:,:) = f(9,:,1,:)
      tmpymS(3,:,:) = f(10,:,1,:)
      tmpymS(4,:,:) = f(16,:,1,:)
      tmpymS(5,:,:) = f(18,:,1,:)
! Start MPI calls
      call MPI_IRECV(tmpypR,ylen,MPI_REAL8,myp,1,MPI_COMM_WORLD,reqr(1),ierr)
      call MPI_IRECV(tmpymR,ylen,MPI_REAL8,mym,0,MPI_COMM_WORLD,reqr(2),ierr)
      call MPI_ISEND(tmpypS,ylen,MPI_REAL8,myp,0,MPI_COMM_WORLD,reqs(1),ierr)
      call MPI_ISEND(tmpymS,ylen,MPI_REAL8,mym,1,MPI_COMM_WORLD,reqs(2),ierr)
! Stream the Y+ inside
      f(3,:,:,:) = cshift(f(3,:,:,:), shift = -1, dim = 2)
      f(7,:,:,:) = cshift(f(7,:,:,:), shift = -1, dim = 2)
      f(8,:,:,:) = cshift(f(8,:,:,:), shift = -1, dim = 2)
      f(15,:,:,:) = cshift(f(15,:,:,:), shift = -1, dim = 2)
      f(17,:,:,:) = cshift(f(17,:,:,:), shift = -1, dim = 2)
! Stream the Y- inside
      f(4,:,:,:) = cshift(f(4,:,:,:), shift = 1, dim = 2)
      f(9,:,:,:) = cshift(f(9,:,:,:), shift = 1, dim = 2)
      f(10,:,:,:) = cshift(f(10,:,:,:), shift = 1, dim = 2)
      f(16,:,:,:) = cshift(f(16,:,:,:), shift = 1, dim = 2)
      f(18,:,:,:) = cshift(f(18,:,:,:), shift = 1, dim = 2)
! Wait for MPI responses then stream the corresponding edge
      do i =1, 2
         call MPI_WAITANY(2, reqr, index, status, ierr)
         !Stream edge (Y+ direction)
         if(status(MPI_SOURCE) == mym)then
            f(3,:,1,:) = tmpymR(1,:,:)
            f(7,:,1,:) = tmpymR(2,:,:)
            f(8,:,1,:) = tmpymR(3,:,:)
            f(15,:,1,:) = tmpymR(4,:,:)
            f(17,:,1,:) = tmpymR(5,:,:)
         else !Stream edge (Y- direction)
            f(4,:,ly,:) = tmpypR(1,:,:)
            f(9,:,ly,:) = tmpypR(2,:,:)
            f(10,:,ly,:) = tmpypR(3,:,:)
            f(16,:,ly,:) = tmpypR(4,:,:)
            f(18,:,ly,:) = tmpypR(5,:,:)
         endif
      enddo
! Create z+ stream buffer Dir 5,11,12,15,16
      tmpzpS(1,:,:) = f(5,:,:,lz)
      tmpzpS(2,:,:) = f(11,:,:,lz)
      tmpzpS(3,:,:) = f(12,:,:,lz)
      tmpzpS(4,:,:) = f(15,:,:,lz)
      tmpzpS(5,:,:) = f(16,:,:,lz)
! Create z- stream buffer Dir 6,13,14,17,18
      tmpzmS(1,:,:) = f(6,:,:,1)
      tmpzmS(2,:,:) = f(13,:,:,1)
      tmpzmS(3,:,:) = f(14,:,:,1)
      tmpzmS(4,:,:) = f(17,:,:,1)
      tmpzmS(5,:,:) = f(18,:,:,1)
! Start MPI calls
      call MPI_IRECV(tmpzpR,zlen,MPI_REAL8,mzp,1,MPI_COMM_WORLD,reqr(1),ierr)
      call MPI_IRECV(tmpzmR,zlen,MPI_REAL8,mzm,0,MPI_COMM_WORLD,reqr(2),ierr)
      call MPI_ISEND(tmpzpS,zlen,MPI_REAL8,mzp,0,MPI_COMM_WORLD,reqs(3),ierr)
      call MPI_ISEND(tmpzmS,zlen,MPI_REAL8,mzm,1,MPI_COMM_WORLD,reqs(4),ierr)
! Stream the z+ inside
      f(5,:,:,:) = cshift(f(5,:,:,:), shift = -1, dim = 3)
      f(11,:,:,:) = cshift(f(11,:,:,:), shift = -1, dim = 3)
      f(12,:,:,:) = cshift(f(12,:,:,:), shift = -1, dim = 3)
      f(15,:,:,:) = cshift(f(15,:,:,:), shift = -1, dim = 3)
      f(16,:,:,:) = cshift(f(16,:,:,:), shift = -1, dim = 3)
! Stream the z- inside
      f(6,:,:,:) = cshift(f(6,:,:,:), shift = 1, dim = 3)
      f(13,:,:,:) = cshift(f(13,:,:,:), shift = 1, dim = 3)
      f(14,:,:,:) = cshift(f(14,:,:,:), shift = 1, dim = 3)
      f(17,:,:,:) = cshift(f(17,:,:,:), shift = 1, dim = 3)
      f(18,:,:,:) = cshift(f(18,:,:,:), shift = 1, dim = 3)
! Wait for MPI responses then stream the corresponding edge
      do i =1, 2
         call MPI_WAITANY(2, reqr, index, status, ierr)
         !Stream side edge (z+ direction)
         if(status(MPI_SOURCE) == mzm)then
            f(5,:,:,1) = tmpzmR(1,:,:)
            f(11,:,:,1) = tmpzmR(2,:,:)
            f(12,:,:,1) = tmpzmR(3,:,:)
            f(15,:,:,1) = tmpzmR(4,:,:)
            f(16,:,:,1) = tmpzmR(5,:,:)
         else !Stream side edge (z- direction)
            f(6,:,:,lz) = tmpzpR(1,:,:)
            f(13,:,:,lz) = tmpzpR(2,:,:)
            f(14,:,:,lz) = tmpzpR(3,:,:)
            f(17,:,:,lz) = tmpzpR(4,:,:)
            f(18,:,:,lz) = tmpzpR(5,:,:)
         endif
     enddo
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

      if(ipart == .TRUE.)then
        !Inject fluid solid particle collisions into fluid domain
        do i=1, ifsc_inject
          f(fsc_inject(i)%ip,fsc_inject(i)%x,fsc_inject(i)%y,fsc_inject(i)%z) = fsc_inject(i)%dist
        enddo
        deallocate(fsc_inject)
      endif
! Wait for any sends that did not complete
      call MPI_WAITALL(4,reqs,status_array,ierr)
      end subroutine streamingp
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
      subroutine exchngy(tmpy1,tmpy2,tag,ydir)
      use mpi
      use var_inc
      implicit none

      real, dimension(5,lx,lz)     :: tmpy1,tmpy2

      integer tag, ydir, ilen, ips, ipr
      integer status_array(MPI_STATUS_SIZE,2), req(2)

      ilen = 5*lx*lz

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


      end subroutine exchngy
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



      
