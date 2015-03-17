      subroutine initpartpos

      use mpi
      use var_inc
      implicit none

      integer ip,cnt, idf,idfy,idfz, i, j, k
      integer npsflag, npsflagt, npst, ovlpcnt
      real ranpp
      real radpgap

      radpgap = rad + mingap_w

      iseedp = -iseedp
      iyp = 0
      ivp = 0
      ovlpflagt = 1

      DO

      if(myid == 0)then
       do ip = 1,npart
         if(ovlpflagt(ip) > 0)then
           ypglb(1,ip) = radpgap +(real(nx)-2.0*radpgap)*ranpp(iseedp,iyp,ivp)
           ypglb(2,ip) = real(ny)*ranpp(iseedp,iyp,ivp)
           ypglb(3,ip) = real(nz)*ranpp(iseedp,iyp,ivp)
         end if
       end do
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ypglb,3*npart,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      nps = 0
      ipglb = 0
      yp = 0.0

      do ip = 1,npart
        idfy = floor(ypglb(2,ip)/real(ny)*real(nprocY))
        idfz = floor(ypglb(3,ip)/real(nz)*real(nprocZ))
        idf =  idfz * nprocY + idfy 
        if(idf == myid)then
          nps = nps + 1
          yp(1,nps) = ypglb(1,ip)
          yp(2,nps) = ypglb(2,ip)
          yp(3,nps) = ypglb(3,ip)
          ipglb(nps) = ip
        end if
      end do

      npsflag = 0
      if(nps > msize)then
        write(*,*) 'fault: nps = ', nps, ' > msize',msize,'nproc,npart',nproc,npart
       npsflag = 1
      end if
      call MPI_ALLREDUCE(npsflag,npsflagt,1,MPI_INTEGER,               &
                         MPI_SUM,MPI_COMM_WORLD,ierr)
      if(npsflagt > 0) stop

      ovlpflag = 0
      ovlpflagt = 0
      call ovlpchk
      call MPI_ALLREDUCE(ovlpflag,ovlpflagt,npart,MPI_INTEGER,         &
                         MPI_SUM,MPI_COMM_WORLD,ierr)
      ovlpcnt = count(ovlpflagt /= 0)
      if(myid == 0) write(*,*) 'ovlpcnt = ', ovlpcnt
      if(ovlpcnt == 0) then
        if(myid==0)write(*,*) 'Particle position initialization done.'
!        write(*,*) 'nps = ',nps
        call MPI_ALLREDUCE(nps,npst,1,MPI_INTEGER,MPI_SUM,             &
                           MPI_COMM_WORLD,ierr)
        if(npst == npart)then
!          write(*,*) 'npst == npart = ', npst
        else
          write(*,*) 'fault: npst != npart'
          stop
        end if
        exit
      end if

      END DO

      end subroutine initpartpos
!===========================================================================

      Function ranpp(idum,iy1,iv1)
      implicit none

!     Minimal random number generator of Park and Miller with
!     Bays-Durham shuffle and added safegards, see Numerical Recipes

      integer idum, IA,IM,IQ,IR,NTAB,NDIV
      real ranpp,AM,EPS,RNMX
      parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,    &
                 NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      integer j,k,iv1(NTAB),iy1

      if(idum.le.0 .or. iy1.eq.0) then
        idum=max(-idum,1)
        do j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if(idum.lt.0) idum=idum+IM
          if(j.le.NTAB) iv1(j)=idum
        enddo
        iy1=iv1(1)
      endif

      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if(idum.lt.0) idum=idum+IM
      j=1+iy1/NDIV
      iy1=iv1(j)
      iv1(j)=idum
      ranpp=min(AM*iy1,RNMX)

      return
      end
!===========================================================================

      subroutine ovlpchk
      use var_inc
      implicit none

      integer ip, id, ii, iflag
      real xc, yc, zc, xb, yb, zb, dist, hgap

      do ip = 1,nps
! particles with center in local processor
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)

        if(id < npart)then

        iflag = 0
        do ii = id+1,npart
          xb = ypglb(1,ii)
          yb = ypglb(2,ii)
          zb = ypglb(3,ii)

! the algorithm below is designed by my wife, Dr. Hui Li
!       if((xb - xc) > real(nxh)) xb = xb - real(nx)
!       if((xb - xc) < -real(nxh)) xb = xb + real(nx)

          if((yb - yc) > real(nyh)) yb = yb - real(ny)
          if((yb - yc) < -real(nyh)) yb = yb + real(ny)

          if((zb - zc) > real(nzh)) zb = zb - real(nz)
          if((zb - zc) < -real(nzh)) zb = zb + real(nz)

          dist = sqrt((xc - xb)**2 + (yc - yb)**2 + (zc -zb)**2)

          hgap = dist - 2.0*rad

          if(istep00 == 0) hgap = hgap - mingap

          if(hgap < 0.0)then
            iflag = 1
             if(istep .gt. irelease)then
              write(*,'(A20,F14.8)')'In ovlpchk hgap < 0:',hgap
              write(*,'(A6,I10)')'Istep ',istep
              write(*,'(A44,3I4)')'Ip(Direction), Id(Bead 1 Id), Ii(Bead 2 Id):',ip,id,ii
              write(*,'(A41,2I4)')'MyId(Processor), nps(Beads in Processor):',myid,nps
              write(*,'(A15)')'X Y Z Positions'
              write(*,'(A7,3(2X,F15.8))')'Bead 1:',xc,yc,zc
              write(*,'(A7,3(2X,F15.8))')'Bead 2:',ypglb(1,ii),ypglb(2,ii),ypglb(3,ii)
              write(*,'(A10)')'=========='
             endif
            exit
          end if
        end do

        if(iflag > 0) ovlpflag(id) = 1

        end if
      end do

      end subroutine ovlpchk
!===========================================================================
! this memory-optimized subroutine is the same as the original version of 
! "initpartvel00" above, except that the large arrays "uxt, uyt, uzt" with
! dimension of (lx,ly,0:lz+1) are no longer employed as in the original version.

      subroutine initpartvel
      use mpi
      use var_inc
      implicit none

      integer ip, id, ilen, idir, iflag   
      real xc, yc, zc, v0      

      real, dimension(3,npart):: wp0
      real,allocatable,dimension(:,:):: tmpL, tmpR
      real,allocatable,dimension(:,:):: tmpU, tmpD


      allocate (tmpL(lx,ly))   
      allocate (tmpR(lx,ly))   
      allocate (tmpU(lx,0:lz+1))
      allocate (tmpD(lx,0:lz+1))

      wp0 = 0.0

! "iflag = 1", indicates velocity interpolation
      iflag = 1

! interpolate ux       

      call exchng1(ux(:,:,1),tmpR,ux(:,:,lz),tmpL,ux(:,1,:),tmpU,ux(:,ly,:),tmpD)

      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)

        idir =  1 

        call trilinear(tmpL,tmpR,tmpU,tmpD,xc,yc,zc,idir,iflag,v0,id)


        wp0(1,id) = v0
      end do

! interpolate uy    

      call exchng1(uy(:,:,1),tmpR,uy(:,:,lz),tmpL,uy(:,1,:),tmpU,uy(:,ly,:),tmpD)


      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)

        idir = 2

        call trilinear(tmpL,tmpR,tmpU,tmpD,xc,yc,zc,idir,iflag,v0,id)

        wp0(2,id) = v0
      end do

! interpolate uz  

      call exchng1(uz(:,:,1),tmpR,uz(:,:,lz),tmpL,uz(:,1,:),tmpU,uz(:,ly,:),tmpD)


      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)

        idir = 3  

        call trilinear(tmpL,tmpR,tmpU,tmpD,xc,yc,zc,idir,iflag,v0,id)

        wp0(3,id) = v0
      end do

      deallocate (tmpL)   
      deallocate (tmpR)   
      deallocate (tmpU)
      deallocate (tmpD)
      
! collect info. for wp
      ilen = 3*npart
      call MPI_ALLREDUCE(wp0,wp,ilen,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

      end subroutine initpartvel
!===========================================================================

      subroutine exchng1(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8)
      use mpi
      use var_inc
      implicit none

!                 tmp1   tmp2   tmp3     tmp4   tmp5    tmp6    tmp7    tmp8
!      exchng1(ux(:,:,1),tmpR,ux(:,:,lz),tmpL,ux(:,1,:),tmpU,ux(:,ly,:),tmpD)


      integer ilenz, ileny
      real, dimension(lx,ly)  :: tmp1, tmp2, tmp3, tmp4
      real, dimension(lx,lz)  :: tmp5, tmp7
      real, dimension(lx,0:lz+1):: tmp5l, tmp6, tmp7l, tmp8  

      integer status_array(MPI_STATUS_SIZE,4), req(4)
     
      ilenz = lx*ly
      ileny = lx*(lz+2)

      call MPI_IRECV(tmp2,ilenz,MPI_REAL8,mzp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp4,ilenz,MPI_REAL8,mzm,1,MPI_COMM_WORLD,req(2),ierr)

      call MPI_ISEND(tmp1,ilenz,MPI_REAL8,mzm,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp3,ilenz,MPI_REAL8,mzp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      tmp5l(:,0) = tmp4(:,1)
      tmp5l(:,1:lz) = tmp5(:,:)
      tmp5l(:,lz+1) = tmp2(:,1)

      tmp7l(:,0) = tmp4(:,ly)
      tmp7l(:,1:lz) = tmp7(:,:)
      tmp7l(:,lz+1) = tmp2(:,ly)

      call MPI_IRECV(tmp6,ileny,MPI_REAL8,myp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp8,ileny,MPI_REAL8,mym,1,MPI_COMM_WORLD,req(2),ierr)

      call MPI_ISEND(tmp5l,ileny,MPI_REAL8,mym,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp7l,ileny,MPI_REAL8,myp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      end subroutine exchng1
!===========================================================================

!===========================================================================
! this is the 2nd version of trilinear interpolation subroutine used 
! in accordance with "initpartvel" and "initpartomg".

      subroutine trilinear(tmpL,tmpR,tmpU,tmpD,xc,yc,zc,idir,iflag,v0,id) 
      use var_inc
      implicit none

!              tmp1      tmp2 tmp3       tmp4 tmp5      tmp6 tmp7       tmp8
!      exchng1(ux(:,:,1),tmpR,ux(:,:,lz),tmpL,ux(:,1,:),tmpU,ux(:,ly,:),tmpD)

      integer idir, iflag, i, j, k, ip1, jp1, kp1, id

      real xc, yc, zc, dx, dy, dz   
      real tx, ty, tz, tx0,  ty0, tz0
      real v0, v1, v2, v3, v4, v5, v6, v7, v8

      real, dimension(lx,ly):: tmpL, tmpR   
      real, dimension(lx,0:lz+1):: tmpU, tmpD

      dx = 1.0
      dy = 1.0
      dz = 1.0

      i = int((xc + 0.5) / dx)
      j = int((yc - real(indy*ly) + 0.5) / dy)
      k = int((zc - real(indz*lz) + 0.5) / dz)

      ip1 = i + 1
      jp1 = j + 1
      kp1 = k + 1

      tx = (xc + 0.5 - real(i)) / dx
      ty = (yc - real(indy*ly) + 0.5 - real(j)) / dy
      tz = (zc - real(indz*lz) + 0.5 - real(k)) / dz

      tx0 = 1.0 - tx
      ty0 = 1.0 - ty
      tz0 = 1.0 - tz

! periodicity
     ! if(i == 0) i = nx
     ! if(ip1 > nx) ip1 = 1

      IF(iflag == 1)THEN  
! translational velocity interpolation  

        if(k ==  0)then

         if( j == 0)then        !case 3L
          v1 = tmpD(i,0)
          v2 = tmpD(ip1,0)
          v3 = tmpL(ip1,1)
          v4 = tmpL(i,1)
          v5 = tmpD(i,1)
          v6 = tmpD(ip1,1)
          if(idir == 1)then
           v7 = ux(ip1,jp1,kp1)
           v8 = ux(i,jp1,kp1)
          endif
          if(idir == 2)then
           v7 = uy(ip1,jp1,kp1)
           v8 = uy(i,jp1,kp1)
          endif
          if(idir == 3)then
           v7 = uz(ip1,jp1,kp1)
           v8 = uz(i,jp1,kp1)
          endif

         else if(jp1 > ly)then  !case 2L

          v1 = tmpL(i,ly)
          v2 = tmpL(ip1,ly)
          v3 = tmpU(ip1,0)
          v4 = tmpU(i,0)
          if(idir == 1)then
           v5 = ux(i,j,kp1)
           v6 = ux(ip1,j,kp1)
          endif
          if(idir == 2)then
           v5 = uy(i,j,kp1)
           v6 = uy(ip1,j,kp1)
          endif
          if(idir == 3)then
           v5 = uz(i,j,kp1)
           v6 = uz(ip1,j,kp1)
          endif
          v7 = tmpU(ip1,1)
          v8 = tmpU(i,1)

         else                   !case 1L

          v1 = tmpL(i,j)
          v2 = tmpL(ip1,j)
          v3 = tmpL(ip1,jp1)
          v4 = tmpL(i,jp1)
          if(idir == 1)then
           v5 = ux(i,j,kp1)
           v6 = ux(ip1,j,kp1)
           v7 = ux(ip1,jp1,kp1)
           v8 = ux(i,jp1,kp1)
          endif
          if(idir == 2)then
           v5 = uy(i,j,kp1)
           v6 = uy(ip1,j,kp1)
           v7 = uy(ip1,jp1,kp1)
           v8 = uy(i,jp1,kp1)
          endif
          if(idir == 3)then
           v5 = uz(i,j,kp1)
           v6 = uz(ip1,j,kp1)
           v7 = uz(ip1,jp1,kp1)
           v8 = uz(i,jp1,kp1)
          endif
         endif

        else if(kp1 > lz)then 

          if(j == 0)then        !case 3R
           v1 = tmpD(i,lz)
           v2 = tmpD(ip1,lz)
           if(idir == 1)then
            v3 = ux(ip1,jp1,k)
            v4 = ux(i,jp1,k)
           endif
           if(idir == 2)then
            v3 = uy(ip1,jp1,k)
            v4 = uy(i,jp1,k)
           endif
           if(idir == 3)then
            v3 = uz(ip1,jp1,k)
            v4 = uz(i,jp1,k)
           endif
           v5 = tmpD(i,lz+1)
           v6 = tmpD(ip1,lz+1)
           v7 = tmpR(ip1,1)
           v8 = tmpR(i,1)

          else if(jp1 > ly)then !case 2R
           if(idir == 1)then
            v1 = ux(i,j,k)
            v2 = ux(ip1,j,k)
           endif
           if(idir == 2)then
            v1 = uy(i,j,k)
            v2 = uy(ip1,j,k)
           endif
           if(idir == 3)then
            v1 = uz(i,j,k)
            v2 = uz(ip1,j,k)
           endif
           v3 = tmpU(ip1,lz)
           v4 = tmpU(i,lz)
           v5 = tmpR(i,ly)
           v6 = tmpR(ip1,ly)
           v7 = tmpU(ip1,lz+1)
           v8 = tmpU(i,lz+1) 

          else                  !case 1R
           if(idir == 1)then
            v1 = ux(i,j,k)
            v2 = ux(ip1,j,k)
            v3 = ux(ip1,jp1,k)
            v4 = ux(i,jp1,k)
           endif
           if(idir == 2)then
            v1 = uy(i,j,k)
            v2 = uy(ip1,j,k)
            v3 = uy(ip1,jp1,k)
            v4 = uy(i,jp1,k)
           endif
           if(idir == 3)then
            v1 = uz(i,j,k)
            v2 = uz(ip1,j,k)
            v3 = uz(ip1,jp1,k)
            v4 = uz(i,jp1,k)
           endif
           v5 = tmpR(i,j)
           v6 = tmpR(ip1,j)
           v7 = tmpR(ip1,jp1)
           v8 = tmpR(i,jp1)
          endif

        else

          if(jp1 > ly)then          !case 5

           if(idir == 1)then
            v1 = ux(i,j,k)
            v2 = ux(ip1,j,k)
            v5 = ux(i,j,kp1)
            v6 = ux(ip1,j,kp1)
           endif
           if(idir == 2)then
            v1 = uy(i,j,k)
            v2 = uy(ip1,j,k)
            v5 = uy(i,j,kp1)
            v6 = uy(ip1,j,kp1)
           endif
           if(idir == 3)then
            v1 = uz(i,j,k)
            v2 = uz(ip1,j,k)
            v5 = uz(i,j,kp1)
            v6 = uz(ip1,j,kp1)
           endif
           v3 = tmpU(ip1,k)
           v4 = tmpU(i,k)
           v7 = tmpU(ip1,kp1)
           v8 = tmpU(i,kp1)

          else if(j==0)then      !case 4 

           v1 = tmpD(i,k)
           v2 = tmpD(ip1,k)
           v5 = tmpD(i,kp1)
           v6 = tmpD(ip1,kp1)
           if(idir == 1)then
            v3 = ux(ip1,jp1,k)
            v4 = ux(i,jp1,k)
            v7 = ux(ip1,jp1,kp1)
            v8 = ux(i,jp1,kp1)
           endif
           if(idir == 2)then
            v3 = uy(ip1,jp1,k)
            v4 = uy(i,jp1,k)
            v7 = uy(ip1,jp1,kp1)
            v8 = uy(i,jp1,kp1)
           endif
           if(idir == 3)then
            v3 = uz(ip1,jp1,k)
            v4 = uz(i,jp1,k)
            v7 = uz(ip1,jp1,kp1)
            v8 = uz(i,jp1,kp1)
           endif

          else !case for center region
           if(idir == 1)then
            v1 = ux(i,j,k)
            v2 = ux(ip1,j,k)
            v3 = ux(ip1,jp1,k)
            v4 = ux(i,jp1,k)
            v5 = ux(i,j,kp1)
            v6 = ux(ip1,j,kp1)
            v7 = ux(ip1,jp1,kp1)
            v8 = ux(i,jp1,kp1)
           endif
           if(idir == 2)then
            v1 = uy(i,j,k)
            v2 = uy(ip1,j,k)
            v3 = uy(ip1,jp1,k)
            v4 = uy(i,jp1,k)
            v5 = uy(i,j,kp1)
            v6 = uy(ip1,j,kp1)
            v7 = uy(ip1,jp1,kp1)
            v8 = uy(i,jp1,kp1)
           endif
           if(idir == 3)then
            v1 = uz(i,j,k)
            v2 = uz(ip1,j,k)
            v3 = uz(ip1,jp1,k)
            v4 = uz(i,jp1,k)
            v5 = uz(i,j,kp1)
            v6 = uz(ip1,j,kp1)
            v7 = uz(ip1,jp1,kp1)
            v8 = uz(i,jp1,kp1)
           endif

          endif

        endif

      ELSEIF(iflag == 2)THEN
!! voritcity interpolation 
        if(k ==  0)then

         if( j == 0)then        !case 3L
          v1 = tmpD(i,0)
          v2 = tmpD(ip1,0)
          v3 = tmpL(ip1,1)
          v4 = tmpL(i,1)
          v5 = tmpD(i,1)
          v6 = tmpD(ip1,1)
          if(idir == 1)then
           v7 = ox(ip1,jp1,kp1)
           v8 = ox(i,jp1,kp1)
          endif
          if(idir == 2)then
           v7 = oy(ip1,jp1,kp1)
           v8 = oy(i,jp1,kp1)
          endif
          if(idir == 3)then
           v7 = oz(ip1,jp1,kp1)
           v8 = oz(i,jp1,kp1)
          endif


         else if(jp1 > ly)then  !case 2L

          v1 = tmpL(i,ly)
          v2 = tmpL(ip1,ly)
          v3 = tmpU(ip1,0)
          v4 = tmpU(i,0)
          if(idir == 1)then
           v5 = ox(i,j,kp1)
           v6 = ox(ip1,j,kp1)
          endif
          if(idir == 2)then
           v5 = oy(i,j,kp1)
           v6 = oy(ip1,j,kp1)
          endif
          if(idir == 3)then
           v5 = oz(i,j,kp1)
           v6 = oz(ip1,j,kp1)
          endif
          v7 = tmpU(ip1,1)
          v8 = tmpU(i,1)

         else                   !case 1L


          v1 = tmpL(i,j)
          v2 = tmpL(ip1,j)
          v3 = tmpL(ip1,jp1)
          v4 = tmpL(i,jp1)
          if(idir == 1)then
           v5 = ox(i,j,kp1)
           v6 = ox(ip1,j,kp1)
           v7 = ox(ip1,jp1,kp1)
           v8 = ox(i,jp1,kp1)
          endif
          if(idir == 2)then
           v5 = oy(i,j,kp1)
           v6 = oy(ip1,j,kp1)
           v7 = oy(ip1,jp1,kp1)
           v8 = oy(i,jp1,kp1)
          endif
          if(idir == 3)then
           v5 = oz(i,j,kp1)
           v6 = oz(ip1,j,kp1)
           v7 = oz(ip1,jp1,kp1)
           v8 = oz(i,jp1,kp1)
          endif
         endif

        else if(kp1 > lz)then

          if(j == 0)then        !case 3R
           v1 = tmpD(i,lz)
           v2 = tmpD(ip1,lz)
           if(idir == 1)then
            v3 = ox(ip1,jp1,k)
            v4 = ox(i,jp1,k)
           endif
           if(idir == 2)then
            v3 = oy(ip1,jp1,k)
            v4 = oy(i,jp1,k)
           endif
           if(idir == 3)then
            v3 = oz(ip1,jp1,k)
            v4 = oz(i,jp1,k)
           endif
           v5 = tmpD(i,lz+1)
           v6 = tmpD(ip1,lz+1)
           v7 = tmpR(ip1,1)
           v8 = tmpR(i,1)


          else if(jp1 > ly)then !case 2R
           if(idir == 1)then
            v1 = ox(i,j,k)
            v2 = ox(ip1,j,k)
           endif
           if(idir == 2)then
            v1 = oy(i,j,k)
            v2 = oy(ip1,j,k)
           endif
           if(idir == 3)then
            v1 = oz(i,j,k)
            v2 = oz(ip1,j,k)
           endif
           v3 = tmpU(ip1,lz)
           v4 = tmpU(i,lz)
           v5 = tmpR(i,ly)
           v6 = tmpR(ip1,ly)
           v7 = tmpU(ip1,lz+1)
           v8 = tmpU(i,lz+1)

          else                  !case 1R
           if(idir == 1)then
            v1 = ox(i,j,k)
            v2 = ox(ip1,j,k)
            v3 = ox(ip1,jp1,k)
            v4 = ox(i,jp1,k)
           endif
           if(idir == 2)then
            v1 = oy(i,j,k)
            v2 = oy(ip1,j,k)
            v3 = oy(ip1,jp1,k)
            v4 = oy(i,jp1,k)
           endif
           if(idir == 3)then
            v1 = oz(i,j,k)
            v2 = oz(ip1,j,k)
            v3 = oz(ip1,jp1,k)
            v4 = oz(i,jp1,k)
           endif
           v5 = tmpR(i,j)
           v6 = tmpR(ip1,j)
           v7 = tmpR(ip1,jp1)
           v8 = tmpR(i,jp1)
          endif

        else

          if(jp1 > ly)then          !case 5

           if(idir == 1)then
            v1 = ox(i,j,k)
            v2 = ox(ip1,j,k)
            v5 = ox(i,j,kp1)
            v6 = ox(ip1,j,kp1)
           endif
           if(idir == 2)then
            v1 = oy(i,j,k)
            v2 = oy(ip1,j,k)
            v5 = oy(i,j,kp1)
            v6 = oy(ip1,j,kp1)
           endif
           if(idir == 3)then
            v1 = oz(i,j,k)
            v2 = oz(ip1,j,k)
            v5 = oz(i,j,kp1)
            v6 = oz(ip1,j,kp1)
           endif
           v3 = tmpU(ip1,k)
           v4 = tmpU(i,k)
           v7 = tmpU(ip1,kp1)
           v8 = tmpU(i,kp1)

         else if(j==0) then         !case 4
           v1 = tmpD(i,k)
           v2 = tmpD(ip1,k)
           v5 = tmpD(i,kp1)
           v6 = tmpD(ip1,kp1)
           if(idir == 1)then
            v3 = ox(ip1,jp1,k)
            v4 = ox(i,jp1,k)
            v7 = ox(ip1,jp1,kp1)
            v8 = ox(i,jp1,kp1)
           endif
           if(idir == 2)then
            v3 = oy(ip1,jp1,k)
            v4 = oy(i,jp1,k)
            v7 = oy(ip1,jp1,kp1)
            v8 = oy(i,jp1,kp1)
           endif
           if(idir == 3)then
            v3 = oz(ip1,jp1,k)
            v4 = oz(i,jp1,k)
            v7 = oz(ip1,jp1,kp1)
            v8 = oz(i,jp1,kp1)
           endif

          else                       !case for center region

           if(idir == 1)then
            v1 = ox(i,j,k)
            v2 = ox(ip1,j,k)
            v3 = ox(ip1,jp1,k)
            v4 = ox(i,jp1,k)
            v5 = ox(i,j,kp1)
            v6 = ox(ip1,j,kp1)
            v7 = ox(ip1,jp1,kp1)
            v8 = ox(i,jp1,kp1)
           endif
           if(idir == 2)then
            v1 = oy(i,j,k)
            v2 = oy(ip1,j,k)
            v3 = oy(ip1,jp1,k)
            v4 = oy(i,jp1,k)
            v5 = oy(i,j,kp1)
            v6 = oy(ip1,j,kp1)
            v7 = oy(ip1,jp1,kp1)
            v8 = oy(i,jp1,kp1)
           endif
           if(idir == 3)then
            v1 = oz(i,j,k)
            v2 = oz(ip1,j,k)
            v3 = oz(ip1,jp1,k)
            v4 = oz(i,jp1,k)
            v5 = oz(i,j,kp1)
            v6 = oz(ip1,j,kp1)
            v7 = oz(ip1,jp1,kp1)
            v8 = oz(i,jp1,kp1)
           endif

          endif

        endif

      END IF 


!      if(id==1)then
!      write(*,*)'id,idir,iflag,v1,v2,v3,v4,v5,v6,v7,v8 =',id,idir,iflag,v1,v2,v3,v4,v5,v6,v7,v8
!      endif

      v0 = tx0*ty0*tz0*v1 + tx*ty0*tz0*v2 + tx*ty*tz0*v3               &
         + tx0*ty*tz0*v4 + tx0*ty0*tz*v5 +tx*ty0*tz*v6                 &
         + tx*ty*tz*v7 + tx0*ty*tz*v8

!      if(idir == 1 .and. iflag ==2)then
!      write(*,*)'id,idir,iflag,v0 = ',id,idir,iflag,v0
!      endif

      end subroutine trilinear
!===========================================================================
! this memory-optimized subroutine is the same as the original version of 
! "initpartomg00" above, except that the large arrays "oxt, oyt, ozt" with
! dimension of (lx,ly,0:lz+1) are no longer employed as in the original version.

      subroutine initpartomg
      use mpi
      use var_inc
      implicit none

      integer ip, id, ilen, idir, iflag  
      real xc, yc, zc, v0  

      real, dimension(3,npart):: omgp0
      real,allocatable,dimension(:,:):: tmpL, tmpR
      real,allocatable,dimension(:,:):: tmpU, tmpD


! prepare vorticity field ox, oy, and oz
      call vortcalc

      allocate (tmpL(lx,ly))
      allocate (tmpR(lx,ly))
      allocate (tmpU(lx,0:lz+1))
      allocate (tmpD(lx,0:lz+1))

      omgp0 = 0.0

! "iflag = 2", indicates vorticity interpolation
      iflag = 2  

! interpolate ox

      call exchng1(ox(:,:,1),tmpR,ox(:,:,lz),tmpL,ox(:,1,:),tmpU,ox(:,ly,:),tmpD)

      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)

        idir = 1
        call trilinear(tmpL,tmpR,tmpU,tmpD,xc,yc,zc,idir,iflag,v0,id)   

! vorticity is twice the rotation rate of a small fluid line segment
! oriented along a principal direction of rate-of-strain tensor
        omgp0(1,id) = 0.5*v0
      end do

! interpolate oy   

      call exchng1(oy(:,:,1),tmpR,oy(:,:,lz),tmpL,oy(:,1,:),tmpU,oy(:,ly,:),tmpD)

      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)

        idir = 2    
        call trilinear(tmpL,tmpR,tmpU,tmpD,xc,yc,zc,idir,iflag,v0,id)   

        omgp0(2,id) = 0.5*v0
      end do

! interpolate oz       

      call exchng1(oz(:,:,1),tmpR,oz(:,:,lz),tmpL,oz(:,1,:),tmpU,oz(:,ly,:),tmpD)

      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)

        idir = 3  
        call trilinear(tmpL,tmpR,tmpU,tmpD,xc,yc,zc,idir,iflag,v0,id)

        omgp0(3,id) = 0.5*v0
      end do

      deallocate (tmpL)
      deallocate (tmpR)
      deallocate (tmpU)
      deallocate (tmpD)

! collect info. for omgp
      ilen = 3*npart
      call MPI_ALLREDUCE(omgp0,omgp,ilen,MPI_REAL8,MPI_SUM,             &
                         MPI_COMM_WORLD,ierr)

!      if(myid==0)then
!      write(*,*)'omgp',omgp
!      endif

110   continue    
! to compare with Elghobashi's data
      omgp = 0.0

      end subroutine initpartomg
!===========================================================================

      subroutine beads_links
      use var_inc
      implicit none

      integer ip,ippp, i, j, k, ipop, imove, jmove, kmove
      real xc, yc, zc, aa, bb, cc    
      real alpha0, alpha1, alpha     
      real xx0, yy0, zz0, rr0, rr01      
      real xx1, yy1, zz1, rr1, rr11,xxtt   
      real xpnt, ypnt, zpnt, radp2 

      ilink = -1
      mlink = -1
      ibnodes = -1
      isnodes = -1 
      alink = 0.0
      radp2 = rad + 2.d0
 
      do ip = 1,npart
        xc = ypglb(1,ip)
        yc = ypglb(2,ip)
        zc = ypglb(3,ip)
 
        do k = 1,lz
          zpnt = dfloat(k) - 0.5d0 + dfloat(indz*lz)

! use the nearest particle center instead of the real center
          if((zc - zpnt) > dfloat(nzh)) zc = zc - dfloat(nz)
          if((zc - zpnt) < -dfloat(nzh)) zc = zc + dfloat(nz)

          zz0 = zpnt - zc
        IF(abs(zz0) <= radp2)THEN

        do j = 1,ly 
          ypnt = dfloat(j) - 0.5d0 + dfloat(indy*ly)

! use the nearest particle center instead of the real center
          if((yc - ypnt) > dfloat(nyh)) yc = yc - dfloat(ny)
          if((yc - ypnt) < -dfloat(nyh)) yc = yc + dfloat(ny)

          yy0 = ypnt - yc
! option 1
!       IF(abs(yy0) <= radp2)THEN
! option 2
        IF(sqrt(yy0*yy0+zz0*zz0) <= radp2)THEN

        do i = 1,lx 
          xpnt = dfloat(i) - 0.5d0    

! use the nearest particle center instead of the real center
          !if((xc - xpnt) > dfloat(nxh)) xc = xc - dfloat(nx)
          !if((xc - xpnt) < -dfloat(nxh)) xc = xc + dfloat(nx)



          xx0 = xpnt - xc
! option 1
!       IF(abs(xx0) <= radp2)THEN
! option 2
          rr0 = xx0**2 + yy0**2 + zz0**2 
         IF(sqrt(rr0) <= radp2)THEN
          rr01 = rr0 - (rad*1.d0)**2 

          if(rr01 <= 0.d0)then
            ibnodes(i,j,k) = 1 
            isnodes(i,j,k) = ip
          else
            do ipop = 1,npop-1
              imove = i + cix(ipop) 
              jmove = j + ciy(ipop) 
              kmove = k + ciz(ipop) 

              xx1 = dfloat(imove) - 0.5d0 - xc
              yy1 = dfloat(jmove) - 0.5d0 - yc + dfloat(indy*ly)
              zz1 = dfloat(kmove) - 0.5d0 - zc + dfloat(indz*lz)
!              rr1 = sqrt(xx1*xx1 + yy1*yy1 + zz1*zz1)
              rr1 = xx1**2 + yy1**2 + zz1**2  
              rr11 = rr1 - (rad*1.d0)**2  

!          if(imove < 1 .or. imove>nx)then
!          ilink(ipop,i,j,k)=2
!          end if

              if(rr11 <= 0.d0)then
! note (i,j,k) is the flow-domain point, and (imove,jmove,kmove) the wall node
! alpha is the percentage of link from wall measured from (i,j,k)
! alpha = percentage of distance from the flow-domain node to actual wall
! [x1+(x0-x1)*alpha]^2 + [y1+(y0-y1)*alpha]^2 + [z1+(z0-z1)*alpha]^2 = rad^2
! then alpha = 1 - alpha

              ilink(ipop,i,j,k) = 1
              mlink(ipop,i,j,k) = ip

              aa = rr0 + rr1 - 2.d0*(xx0*xx1 + yy0*yy1 + zz0*zz1)
              bb = xx1*(xx0-xx1) + yy1*(yy0-yy1) + zz1*(zz0-zz1)    
              cc = rr11  
              alpha0 = bb/aa 
              alpha1 = sqrt(alpha0**2 - cc/aa)    
              alpha = -alpha0 + alpha1

              alpha = 1.d0 - alpha
              alink(ipop,i,j,k) = alpha

              if(alpha < 0.d0 .or. alpha > 1.d0)then 
                write(*,*) 'fault: alpha = ', alpha, ' @ ilink(',      &
                           ipop, ', ', i, ', ', j, ', ', k, ')' 
!                stop
              end if

              end if

            end do
          end if

        END IF
        end do 
        END IF
        end do 
        END IF
        end do 
      end do 

!      do j = 1,32 
!       do k = 5,8  
!        do i = 1,lx 

!         if(istep==630 .and. myid == 1)then
!         write(90,*)' ',ibnodes(i,j,k)
!          write(90,*)(ilink(ippp,i,j,k),ippp=1,18)
!         endif

!        enddo
!       enddo
!      enddo

      end subroutine beads_links
!===========================================================================
! this is a modified version of the subroutine "beads_collision" above.
! the goal is to save some run time memory by avoiding employing large
! arrays, such as f9(0:npop-1,lx,ly,-1:lz+2), and ibnodes9(lx,ly,-1:lz+2) 

      subroutine beads_collision
      use mpi 
      use var_inc
      implicit none

      integer id, i, j, k, ii, jj, kk,  ipop, ipp, ix, iy, iz, ilen 
      integer im1, im2, ip1, jm1, jm2, jp1, km1, km2, kp1   
      integer ibm1, ibm2
      real alpha, xc, yc, zc, xx0, yy0, zz0
      real uwx, uwy, uwz, uwpro, ff1, ff2, ff3 
      real w1, w2, w3, omg1, omg2, omg3 
      real c1, c2, c3, dff, dxmom, dymom, dzmom 
      real xpnt, ypnt, zpnt, f9   

       character (len = 100):: fnm2
!      character (len = 100):: fnm_a,fnm_a1,fnm_a2,fnm_a3,fnm_a4,fnm_a5,fnm_a6,fnm_a7,fnm_a8
!      character (len = 100):: fnm_b,fnm_b1,fnm_b2,fnm_b3,fnm_b4,fnm_b5,fnm_b6,fnm_b7,fnm_b8
!      character (len = 100):: fnm20,fnm21,fnm25,fnm26,fnm27,fnm28

      real,dimension(lx,ly,lz):: f9print,alphaprint,ff1print

      real,allocatable,dimension(:,:,:,:):: tmpfL, tmpfR 
      real,allocatable,dimension(:,:,:,:):: tmpfU, tmpfD

      real,allocatable,dimension(:,:,:):: tmp3, tmp4    
      real,allocatable,dimension(:,:,:):: tmp5, tmp6

      integer,allocatable,dimension(:,:,:):: tmpiL, tmpiR
      integer,allocatable,dimension(:,:,:):: tmpiU, tmpiD

      integer,allocatable,dimension(:,:,:):: tmp3i, tmp4i
      integer,allocatable,dimension(:,:,:):: tmp5i, tmp6i

      integer,allocatable,dimension(:,:,:):: if9L, if9R  
      integer,allocatable,dimension(:,:,:):: if9U, if9D

      real, dimension(3,npart):: fHIp0, torqp0 

!      fnm = '/ptmp/lwang/debug2D/ibnodes0.2D.dat'
!      fnm1 = '/ptmp/lwang/debug2D/ibnodes1.2D.dat'
!      fnm5 = '/ptmp/lwang/debug2D/ibnodes5.2D.dat'
!      fnm6 = '/ptmp/lwang/debug2D/ibnodes6.2D.dat'
!      fnm7 = '/ptmp/lwang/debug2D/ibnodes7.2D.dat'
!      fnm8 = '/ptmp/lwang/debug2D/ibnodes8.2D.dat'

!      fnm_a = '/ptmp/lwang/debug2D/fcheck.dat'
!      fnm_a1 = '/ptmp/lwang/debug2D/fcheck1.dat'
!      fnm_a2 = '/ptmp/lwang/debug2D/fcheck2.dat'
!      fnm_a3 = '/ptmp/lwang/debug2D/fcheck3.dat'
!      fnm_a4 = '/ptmp/lwang/debug2D/fcheck4.dat'
!      fnm_a5 = '/ptmp/lwang/debug2D/fcheck5.dat'
!      fnm_a6 = '/ptmp/lwang/debug2D/fcheck6.dat'
!      fnm_a7 = '/ptmp/lwang/debug2D/fcheck7.dat'
!      fnm_a8 = '/ptmp/lwang/debug2D/fcheck8.dat'

!      fnm_b = '/ptmp/lwang/debug2D/ibcheck.dat'
!      fnm_b1 = '/ptmp/lwang/debug2D/ibcheck1.dat'
!      fnm_b2 = '/ptmp/lwang/debug2D/ibcheck2.dat'
!      fnm_b3 = '/ptmp/lwang/debug2D/ibcheck3.dat'
!      fnm_b4 = '/ptmp/lwang/debug2D/ibcheck4.dat'
!      fnm_b5 = '/ptmp/lwang/debug2D/ibcheck5.dat'
!      fnm_b6 = '/ptmp/lwang/debug2D/ibcheck6.dat'
!      fnm_b7 = '/ptmp/lwang/debug2D/ibcheck7.dat'
!      fnm_b8 = '/ptmp/lwang/debug2D/ibcheck8.dat'

!      fnm20 = '/ptmp/lwang/debug2D/ibnodes20.2D.dat'
!      fnm21 = '/ptmp/lwang/debug2D/ibnodes21.2D.dat'
!      fnm25 = '/ptmp/lwang/debug2D/ibnodes25.2D.dat'
!      fnm26 = '/ptmp/lwang/debug2D/ibnodes26.2D.dat'
!      fnm27 = '/ptmp/lwang/debug2D/ibnodes27.2D.dat'
!      fnm28 = '/ptmp/lwang/debug2D/ibnodes28.2D.dat'


!      open(10, file = trim(fnm), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(11, file = trim(fnm1), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(15, file = trim(fnm5), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(16, file = trim(fnm6), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(17, file = trim(fnm7), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(18, file = trim(fnm8), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')

!      open(20, file = trim(fnm_a), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(21, file = trim(fnm_a1), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(22, file = trim(fnm_a2), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(23, file = trim(fnm_a3), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(24, file = trim(fnm_a4), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(25, file = trim(fnm_a5), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(26, file = trim(fnm_a6), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(27, file = trim(fnm_a7), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(28, file = trim(fnm_a8), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')

!      open(30, file = trim(fnm_b), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(31, file = trim(fnm_b1), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(32, file = trim(fnm_b2), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(33, file = trim(fnm_b3), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(34, file = trim(fnm_b4), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(35, file = trim(fnm_b5), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(36, file = trim(fnm_b6), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(37, file = trim(fnm_b7), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(38, file = trim(fnm_b8), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
 
!      open(20, file = trim(fnm20), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(21, file = trim(fnm21), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(25, file = trim(fnm25), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(26, file = trim(fnm26), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(27, file = trim(fnm27), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(28, file = trim(fnm28), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')


      allocate(tmpfL(0:npop-1,lx,ly,-1:0))
      allocate(tmpfR(0:npop-1,lx,ly,lz+1:lz+2))
      allocate(tmpfU(0:npop-1,lx,ly+1:ly+2,-1:lz+2))
      allocate(tmpfD(0:npop-1,lx,-1:0,-1:lz+2))


!      if(istep.eq.630)then
!      if(myid.eq.56) write(21,*)' before L',((j,k,f(1,5,j,k),j=1,ly),k=lz-1,lz)
!      if(myid.eq.8) write(22,*)' before R',((j,k,f(1,5,j,k),j=1,ly),k=1,2)
!      if(myid.eq.1) write(23,*)' before U',((j,k,f(1,5,j,k),j=1,2),k=1,lz)
!      if(myid.eq.7) write(24,*)' before D',((j,k,f(1,5,j,k),j=ly-1,ly),k=1,lz)

!      if(myid.eq.57) write(25,*)' before UL',((j,k,f(1,5,j,k),j=1,2),k=lz-1,lz)
!      if(myid.eq.63) write(26,*)'before DL',((j,k,f(1,5,j,k),j=ly-1,ly),k=lz-1,lz)
!      if(myid.eq.9) write(27,*)' before UR',((j,k,f(1,5,j,k),j=1,2),k=1,2)
!      if(myid.eq.15) write(28,*)' before DR',((j,k,f(1,5,j,k),j=ly-1,ly),k=1,2)
!      if(myid.eq.56) close(21)
!      if(myid.eq.8) close(22)
!      if(myid.eq.1) close(23)
!      if(myid.eq.7) close(24)
!      if(myid.eq.57) close(25)
!      if(myid.eq.63) close(26)
!      if(myid.eq.9) close(27)
!      if(myid.eq.15) close(28)
!      end if

   
      call exchng2(f(:,:,:,1:2),tmpfR,f(:,:,:,lz-1:lz),tmpfL,    &
                   f(:,:,1:2,:),tmpfU,f(:,:,ly-1:ly,:),tmpfD)

!      if(istep.eq.630)then
!      if(myid.eq.0) then
!      write(20,*)'L',((j,k,tmpfL(1,5,j,k),j=1,ly),k=-1,0)
!      write(20,*)'R',((j,k,tmpfR(1,5,j,k),j=1,ly),k=lz+1,lz+2)
!      write(20,*)'U',((j,k,tmpfU(1,5,j,k),j=ly+1,ly+2),k=1,lz)
!      write(20,*)'D',((j,k,tmpfD(1,5,j,k),j=-1,0),k=1,lz)

!      write(20,*)'UL',((j,k,tmpfU(1,5,j,k),j=ly+1,ly+2),k=-1,0)
!      write(20,*)'DL',((j,k,tmpfD(1,5,j,k),j=-1,0),k=-1,0)
!      write(20,*)'UR',((j,k,tmpfU(1,5,j,k),j=ly+1,ly+2),k=lz+1,lz+2)
!      write(20,*)'DR',((j,k,tmpfD(1,5,j,k),j=-1,0),k=lz+1,lz+2)
!       close(20)
!      end if
!      end if

      allocate (tmpiL(lx,ly,-1:0))
      allocate (tmpiR(lx,ly,lz+1:lz+2))
      allocate (tmpiU(lx,ly+1:ly+2,-1:lz+2))
      allocate (tmpiD(lx,-1:0,-1:lz+2))

!       do ii= 1,nx
!       do jj= 1,ly 
!       do kk= 1,2 
!       if(istep==630 .and. myid==53 .and. ii==1 .and. jj==ly-1 .and. kk==25)write(*,*)'Writing out here (1)' 
!       if(istep==631 .and. myid==35 .and. ii==1 .and. jj==ly-1 .and. kk==25)write(*,*)'Writing out here (1)'
!       if(istep==630 .and. myid==53)write(10,210) ibnodes(ii,jj,kk),ii,jj,kk
!       if(istep==631 .and. myid==35)write(20,210) ibnodes(ii,jj,kk),ii,jj,kk
!       enddo
!       enddo
!       enddo

!       do ii = 1,nx
!       do jj = 1,ly 
!       do kk = lz+1,lz+2 
!       if(istep==630 .and. myid==36 .and. ii==1 .and. jj==-1 .and. kk==25)write(*,*)'Writing out here (2)'
!       if(istep==631 .and. myid==36 .and. ii==1 .and. jj==-1 .and. kk==25)write(*,*)'Writing out here (2)'
!       if(istep==630 .and. myid==52)write(11,210)tmpiR(ii,jj,kk),ii,jj,kk
!       if(istep==631 .and. myid==36)write(21,210)tmpiD(ii,jj,kk),ii,jj,kk
!       enddo
!       enddo
!       enddo

!      close(10)
!      close(11)


!     if(myid==0)write(*,*)'I am in collision  (1)'
!     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!      if(istep.eq.630)then
!      if(myid.eq.56) write(31,*)' before L',((j,k,ibnodes(5,j,k),j=1,ly),k=lz-1,lz)
!      if(myid.eq.8) write(32,*)' before R',((j,k,ibnodes(5,j,k),j=1,ly),k=1,2)
!      if(myid.eq.1) write(33,*)' before U',((j,k,ibnodes(5,j,k),j=1,2),k=1,lz)
!      if(myid.eq.7) write(34,*)' before D',((j,k,ibnodes(5,j,k),j=ly-1,ly),k=1,lz)

!      if(myid.eq.57) write(35,*)' before UL',((j,k,ibnodes(5,j,k),j=1,2),k=lz-1,lz)
!      if(myid.eq.63) write(36,*)' before DL',((j,k,ibnodes(5,j,k),j=ly-1,ly),k=lz-1,lz)
!      if(myid.eq.9) write(37,*)' before UR',((j,k,ibnodes(5,j,k),j=1,2),k=1,2)
!      if(myid.eq.15) write(38,*)' before DR',((j,k,ibnodes(5,j,k),j=ly-1,ly),k=1,2)

!!     if(myid.eq.56) close(31)
!!     if(myid.eq.8) close(32)
!!     if(myid.eq.1) close(33)
!!     if(myid.eq.7) close(34)
!!     if(myid.eq.57) close(35)
!!     if(myid.eq.63) close(36)
!!     if(myid.eq.9) close(37)
!!     if(myid.eq.15) close(38)
!      end if

      call exchng2iNew(ibnodes(:,:,1:2),tmpiR,ibnodes(:,:,lz-1:lz),tmpiL, &
                       ibnodes(:,1:2,:),tmpiU,ibnodes(:,ly-1:ly,:),tmpiD)

!      if(istep.eq.630)then
!      if(myid.eq.0) then
!      write(30,*)'L',((j,k,tmpiL(5,j,k),j=1,ly),k=-1,0)
!      write(30,*)'R',((j,k,tmpiR(5,j,k),j=1,ly),k=lz+1,lz+2)
!      write(30,*)'U',((j,k,tmpiU(5,j,k),j=ly+1,ly+2),k=1,lz)
!      write(30,*)'D',((j,k,tmpiD(5,j,k),j=-1,0),k=1,lz)

!      write(30,*)'UL',((j,k,tmpiU(5,j,k),j=ly+1,ly+2),k=-1,0)
!      write(30,*)'DL',((j,k,tmpiD(5,j,k),j=-1,0),k=-1,0)
!      write(30,*)'UR',((j,k,tmpiU(5,j,k),j=ly+1,ly+2),k=lz+1,lz+2)
!      write(30,*)'DR',((j,k,tmpiD(5,j,k),j=-1,0),k=lz+1,lz+2)
!       close(30)
!      end if
!      end if

!       do ii = 1,nx
!       do jj = 1,ly
!       do kk = 1,2
!       if(istep==630 .and. myid==35  .and. ii==1 .and. jj==ly-1 .and. kk==25)write(*,*)'Writing out here (7)'
!       if(istep==631 .and. myid==35  .and. ii==1 .and. jj==ly-1 .and. kk==25)write(*,*)'Writing out here (7)'
!       if(istep==630 .and. myid==53)write(15,210)ibnodes(ii,jj,kk),ii,jj,kk
!       if(istep==631 .and. myid==35)write(25,210)ibnodes(ii,jj,kk),ii,jj,kk
!       enddo
!       enddo
!       enddo


!       do ii= 1,nx
!       do jj= 1,ly 
!       do kk = lz+1,lz+2 
!       if(istep==630 .and. myid==52 .and. ii==1 .and. jj==-1 .and. kk==25)write(*,*)'Writing out here (6)'
!       if(istep==631 .and. myid==36 .and. ii==1 .and. jj==-1 .and. kk==25)write(*,*)'Writing out here (6)'
!       if(istep==630 .and. myid==52)write(16,210)tmpiR(ii,jj,kk),ii,jj,kk
!       if(istep==631 .and. myid==36)write(26,210)tmpiD(ii,jj,kk),ii,jj,kk
!       enddo
!       enddo
!       enddo


!      close(15)
!      close(16)
!      close(25)
!      close(26)


      fHIp0 = 0.0
      torqp0 = 0.0

      allocate (if9L(0:npop-1,lx,ly)) 
      allocate (if9R(0:npop-1,lx,ly)) 
      allocate (if9U(0:npop-1,lx,0:lz+1))
      allocate (if9D(0:npop-1,lx,0:lz+1))

      if9L = 0
      if9R = 0
      if9U = 0
      if9D = 0

      do k = 1,lz
      do j = 1,ly 
      do i = 1,lx 

      do ipop = 1,npop-1
        if(ilink(ipop,i,j,k) > 0)then

        id = mlink(ipop,i,j,k) 
        alpha = alink(ipop,i,j,k) 

        ipp = ipopp(ipop)       

        xc = ypglb(1,id) 
        yc = ypglb(2,id) 
        zc = ypglb(3,id) 

        xpnt = real(i) - 0.5
        ypnt = real(j) - 0.5 + real(indy*ly)
        zpnt = real(k) - 0.5 + real(indz*lz)

! use the nearest particle center instead of the real center
!        if((xc - xpnt) > real(nxh)) xc = xc - real(nx)
!        if((xc - xpnt) < -real(nxh)) xc = xc + real(nx)

        if((yc - ypnt) > real(nyh)) yc = yc - real(ny)
        if((yc - ypnt) < -real(nyh)) yc = yc + real(ny)

        if((zc - zpnt) > real(nzh)) zc = zc - real(nz)
        if((zc - zpnt) < -real(nzh)) zc = zc + real(nz)

        w1 = wp(1,id)
        w2 = wp(2,id)
        w3 = wp(3,id)

        omg1 = omgp(1,id)
        omg2 = omgp(2,id)
        omg3 = omgp(3,id)

        ix = cix(ipop)
        iy = ciy(ipop)
        iz = ciz(ipop)

        im1 = i - ix
        jm1 = j - iy
        km1 = k - iz
        im2 = i - 2*ix
        jm2 = j - 2*iy
        km2 = k - 2*iz
        ip1 = i + ix
        jp1 = j + iy
        kp1 = k + iz

! periodicity
!        if(im1 < 1) im1 = im1 + nx 
!        if(im1 > nx) im1 = im1 - nx 

!        if(im2 < 1) im2 = im2 + nx  
!        if(im2 > nx) im2 = im2 - nx  

!        if(ip1 < 1) ip1 = ip1 + nx  
!        if(ip1 > nx) ip1 = ip1 - nx  

        xx0 = xpnt + real(ix)*alpha - xc 
        yy0 = ypnt + real(iy)*alpha - yc 
        zz0 = zpnt + real(iz)*alpha - zc 

        ff1 = f(ipop,i,j,k)

        uwx = w1 + omg2*zz0 - omg3*yy0
        uwy = w2 + omg3*xx0 - omg1*zz0
        uwz = w3 + omg1*yy0 - omg2*xx0

        uwpro = uwx*real(ix) + uwy*real(iy) + uwz*real(iz)


! bounce-back scheme on moving particle surfaces
        IF(alpha <= 0.5)THEN

!--------
        if(im1 > nx) then
          ibm1 = 2
        else if(im1 < 1) then
          ibm1 = 2
          else
          if(jm1 > ly) then
          ff2 = tmpfU(ipop,im1,jm1,km1)
          ibm1 = tmpiU(im1,jm1,km1)
          else if (jm1 < 1) then
          ff2 = tmpfD(ipop,im1,jm1,km1)
          ibm1 = tmpiD(im1,jm1,km1)
          else
            if(km1 > lz) then
            ff2 = tmpfR(ipop,im1,jm1,km1)
            ibm1 = tmpiR(im1,jm1,km1)
            else if(km1 < 1 ) then
            ff2 = tmpfL(ipop,im1,jm1,km1)
            ibm1 = tmpiL(im1,jm1,km1)
            else
              ff2 = f(ipop,im1,jm1,km1)
              ibm1 = ibnodes(im1,jm1,km1)
              end if
            end if
          end if
! ------
         if(im2 > nx) then
          ibm2 = 2
         else if(im2 < 1) then
          ibm2 = 2
          else
          if(jm2 > ly) then
          ff3 = tmpfU(ipop,im2,jm2,km2)
          ibm2 = tmpiU(im2,jm2,km2)
          else if (jm2 < 1) then
          ff3 = tmpfD(ipop,im2,jm2,km2)
          ibm2 = tmpiD(im2,jm2,km2)
          else
            if(km2 > lz) then
            ff3 = tmpfR(ipop,im2,jm2,km2)
            ibm2 = tmpiR(im2,jm2,km2)
            else if(km2 < 1 ) then
            ff3 = tmpfL(ipop,im2,jm2,km2)
            ibm2 = tmpiL(im2,jm2,km2)
            else
            ff3 = f(ipop,im2,jm2,km2)
            ibm2 = ibnodes(im2,jm2,km2)
            end if
           end if
          end if
!------------
          if(ibm1 > 0)then
! no interpolation, use simple bounce back
            f9 = ff1 - 6.0*wwp(ipop)*uwpro

          else if(ibm2 > 0)then
! use 2-point interpolation
            f9 = 2.0*alpha*(ff1 - ff2) + ff2 - 6.0*wwp(ipop)*uwpro

          else
! use 3-point interpolation scheme of Lallemand and Luo (2003) JCP
            c1 = alpha*(1.0 + 2.0*alpha)
            c2 = 1.0 - 4.0*alpha*alpha
            c3 = -alpha*(1.0 - 2.0*alpha)
            f9 = c1*ff1 + c2*ff2 + c3*ff3 - 6.0*wwp(ipop)*uwpro
          end if

!          if(id==460 .and. ipop==8 .and. i==197 .and. j==1 .and. k==28 .and. istep==630)then
!           write(*,*)'im1,jm1,km1,im2,jm2,km2,ff1,ff2,ff3,ibm1,ibm2,alpha,wwp(ipop),uwpro',  &
!                      im1,jm1,km1,im2,jm2,km2,ff1,ff2,ff3,ibm1,ibm2,alpha,wwp(ipop),uwpro
!          endif

        ELSE 

          ff2 = f(ipp,i,j,k)
          
          if(im1 > nx) then
          ibm1 = 2
          else if(im1 < 1) then
          ibm1 = 2
          else
          if(jm1 > ly) then
          ff3 = tmpfU(ipp,im1,jm1,km1)
          ibm1 = tmpiU(im1,jm1,km1)
          else if (jm1 < 1) then
          ff3 = tmpfD(ipp,im1,jm1,km1)
          ibm1 = tmpiD(im1,jm1,km1)
          else
            if(km1 > lz) then
            ff3 = tmpfR(ipp,im1,jm1,km1)
            ibm1 = tmpiR(im1,jm1,km1)
            else if(km1 < 1 ) then
            ff3 = tmpfL(ipp,im1,jm1,km1)
            ibm1 = tmpiL(im1,jm1,km1)
            else
            ff3 = f(ipp,im1,jm1,km1)
            ibm1 = ibnodes(im1,jm1,km1)
          end if
            end if
          end if

          if(ibm1 > 0)then
! use 2-point interpolation
            c1 = 0.5 / alpha
            c2 = 1.0 - c1
            f9 = c1*ff1 + c2*ff2 - 6.0*wwp(ipop)*c1*uwpro

          else
! use 3-point interpolation scheme of Lallemand and Luo (2003) JCP
            c1 = 1.0 / alpha / (2.0*alpha + 1.0)
            c2 = (2.0*alpha - 1.0) / alpha
            c3 = 1.0 - c1 - c2
            f9 = c1*ff1 + c2*ff2 + c3*ff3 - 6.0*wwp(ipop)*c1*uwpro
          end if

       END IF

!       if(id==6393 .and. ipop==10 .and. i==132 .and. j==32 .and. k==26 .and. istep==630)then
!       write(*,*)'im1,jm1,km1,im2,jm2,km2,ff1,ff2,ff3,f9,ibm1,ibm2,', &
!                 'alpha,wwp(ipop),uwpro',                             &
!                  im1,jm1,km1,im2,jm2,km2,ff1,ff2,ff3,f9,ibm1,ibm2,   &
!                  alpha,wwp(ipop),uwpro
!       endif

!       do ii= 1,nx
!       do jj= -1,0
!       do kk = 25,28 
!       if(istep==630 .and. myid==36)write(*,*)'Writing out here (8)'
!       if(istep==630 .and. myid==36)write(17,210)tmpiD(ii,jj,kk),ii,jj,kk
!       enddo
!       enddo
!       enddo

!       do ii = 1,nx
!       do jj = ly-1,ly
!       do kk = 25,28 
!       if(istep==630 .and. myid==35)write(*,*)'Writing out here (9)'
!       if(istep==630 .and. myid==35)write(18,210)ibnodes(ii,jj,kk),ii,jj,kk
!       enddo
!       enddo
!       enddo

!       close(17)
!       close(18)

!       if(istep==630 .and. id == 6395)then
!        write(*,*)'id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k',  &
!                   id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k
!       endif 

!       if(istep==630 .and. id == 460)then
!        write(*,*)'id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k',  &
!                   id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k
!       endif

!       if(istep==631 .and. id == 1)then
!        write(*,*)'id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k',  &
!                   id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k
!       endif

!       if(istep==631 .and. id == 2)then
!        write(*,*)'id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k',  &
!                   id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k
!       endif


!       if(istep==630 .and. id == 632)then
!        write(*,*)'id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k',  &
!                   id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k
!       endif

!       if(istep==630 .and. id == 1169)then
!        write(*,*)'id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k',  &
!                   id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k
!       endif

!       if(istep==630 .and. id == 1324)then
!        write(*,*)'id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k',  &
!                   id,ipop,myid,ix,iy,iz,f9,ff1,alpha,i,j,k
!       endif


          if(jp1 > ly) then
          tmpfU(ipp,ip1,jp1,kp1) = f9
          if9U(ipp,ip1,kp1) = 1
          else if (jp1 < 1) then
          tmpfD(ipp,ip1,jp1,kp1) = f9
          if9D(ipp,ip1,kp1) = 1
          else
            if(kp1 > lz) then
            tmpfR(ipp,ip1,jp1,kp1) = f9
            if9R(ipp,ip1,jp1) = 1
            else if(kp1 < 1 ) then
            tmpfL(ipp,ip1,jp1,kp1) = f9
            if9L(ipp,ip1,jp1) = 1
            else
            f(ipp,ip1,jp1,kp1) = f9
            end if
          end if

! compute force and torque acting on particles
        dff = ff1 + f9
        dxmom = dff*real(ix)
        dymom = dff*real(iy)
        dzmom = dff*real(iz)

        fHIp0(1,id) = fHIp0(1,id) + dxmom
        fHIp0(2,id) = fHIp0(2,id) + dymom
        fHIp0(3,id) = fHIp0(3,id) + dzmom

        torqp0(1,id) = torqp0(1,id) + dzmom*yy0 - dymom*zz0
        torqp0(2,id) = torqp0(2,id) + dxmom*zz0 - dzmom*xx0 
        torqp0(3,id) = torqp0(3,id) + dymom*xx0 - dxmom*yy0
 
        end if 
      end do 
      end do 
      end do 
      end do 


!       do ii= 1,nx
!       do jj= -1,0
!       do kk = 25,28
!       if(istep==630 .and. myid==36 .and. ii==1 .and. jj==-1 .and. kk==25)write(*,*)'Writing out here (8)'
!       if(istep==631 .and. myid==36 .and. ii==1 .and. jj==-1 .and. kk==25)write(*,*)'Writing out here (8)'
!       if(istep==630 .and. myid==36)write(17,210)tmpiD(ii,jj,kk),ii,jj,kk
!       if(istep==631 .and. myid==36)write(27,210)tmpiD(ii,jj,kk),ii,jj,kk
!       enddo
!       enddo
!       enddo

!       do ii = 1,nx
!       do jj = ly-1,ly
!       do kk = 25,28
!       if(istep==630 .and. myid==35 .and. ii==1 .and. jj==ly-1 .and. kk==25)write(*,*)'Writing out here (9)'
!       if(istep==631 .and. myid==35 .and. ii==1 .and. jj==ly-1 .and. kk==25)write(*,*)'Writing out here (9)'
!       if(istep==630 .and. myid==35)write(18,210)ibnodes(ii,jj,kk),ii,jj,kk
!       if(istep==631 .and. myid==35)write(28,210)ibnodes(ii,jj,kk),ii,jj,kk
!       enddo
!       enddo
!       enddo

!       close(17)
!       close(18)
!       close(27)
!       close(28)


      deallocate(tmpiL)
      deallocate(tmpiR)
      deallocate(tmpiU)
      deallocate(tmpiD)

! We combine the two communications and the order is in reverse inside
! Also the re-assignment of f is done inside the suroutine
! So exchng3i is not needed any more.

      call exchng3(tmpfL(:,:,:,0),tmpfR(:,:,:,lz+1), &
              tmpfD(:,:,0,0:lz+1),tmpfU(:,:,ly+1,0:lz+1), &
             if9L,if9R,if9D,if9U)

      deallocate(tmpfL)   
      deallocate(tmpfR)  
      deallocate(tmpfU)
      deallocate(tmpfD)

      deallocate(if9U)
      deallocate(if9D)
      deallocate(if9L)   
      deallocate(if9R) 
    
! collect info. for fHIp, and torqp
      ilen = 3*npart

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(fHIp0,fHIp,ilen,MPI_REAL8,MPI_SUM,             &
                         MPI_COMM_WORLD,ierr)      
      call MPI_ALLREDUCE(torqp0,torqp,ilen,MPI_REAL8,MPI_SUM,           &
                         MPI_COMM_WORLD,ierr)      

!      if(istep==630)write(*,*)'fHIp(:,292) = ',fHIp(:,292)

!       do i = 1,3
!        do j = 1,npart
!         if(fHIp(i,j) .gt. 100)then
!          if(myid==0)write(*,*)'istep, fHIp(i,j) = ',istep,fHIp(i,j),i,j 
!         endif
!        enddo
!       enddo


!210   format(2x,4i5)
      end subroutine beads_collision
!===========================================================================
      subroutine exchng2(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8)
      use mpi
      use var_inc
      implicit none

!                   tmp1         tmp2  tmp3             tmp4  tmp5         tmp6  tmp7             tmp8
!      call exchng2(f(:,:,:,1:2),tmpfR,f(:,:,:,lz-1:lz),tmpfL,f(:,:,1:2,:),tmpfU,f(:,:,ly-1:ly,:),tmpfD)

      integer ileny, ilenz
      real, dimension(0:npop-1,lx,ly,2):: tmp1, tmp2, tmp3, tmp4
      real, dimension(0:npop-1,lx,2,lz):: tmp5, tmp7
      real, dimension(0:npop-1,lx,2,-1:lz+2):: tmp5l, tmp6, tmp7l, tmp8

      integer error,status_array(MPI_STATUS_SIZE,4), req(4)

      ilenz = npop*lx*ly*2
      ileny = npop*lx*(lz+4)*2

      call MPI_IRECV(tmp2,ilenz,MPI_REAL8,mzp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp4,ilenz,MPI_REAL8,mzm,1,MPI_COMM_WORLD,req(2),ierr)

      call MPI_ISEND(tmp1,ilenz,MPI_REAL8,mzm,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp3,ilenz,MPI_REAL8,mzp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

      tmp5l(:,:,1,-1) = tmp4(:,:,1,1)
      tmp5l(:,:,1,0) = tmp4(:,:,1,2)
      tmp5l(:,:,2,-1) = tmp4(:,:,2,1)
      tmp5l(:,:,2,0) = tmp4(:,:,2,2)

      tmp5l(:,:,:,1:lz) = tmp5(:,:,:,:)

      tmp5l(:,:,1,lz+1) = tmp2(:,:,1,1)
      tmp5l(:,:,1,lz+2) = tmp2(:,:,1,2)
      tmp5l(:,:,2,lz+1) = tmp2(:,:,2,1)
      tmp5l(:,:,2,lz+2) = tmp2(:,:,2,2)

      tmp7l(:,:,1,-1) = tmp4(:,:,ly-1,1)
      tmp7l(:,:,1,0) = tmp4(:,:,ly-1,2)
      tmp7l(:,:,2,-1) = tmp4(:,:,ly,1)
      tmp7l(:,:,2,0) = tmp4(:,:,ly,2)

      tmp7l(:,:,:,1:lz) = tmp7(:,:,:,:)

      tmp7l(:,:,1,lz+1) = tmp2(:,:,ly-1,1)
      tmp7l(:,:,1,lz+2) = tmp2(:,:,ly-1,2)
      tmp7l(:,:,2,lz+1) = tmp2(:,:,ly,1)
      tmp7l(:,:,2,lz+2) = tmp2(:,:,ly,2)

      call MPI_IRECV(tmp6,ileny,MPI_REAL8,myp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp8,ileny,MPI_REAL8,mym,1,MPI_COMM_WORLD,req(2),ierr)

      call MPI_ISEND(tmp5l,ileny,MPI_REAL8,mym,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp7l,ileny,MPI_REAL8,myp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

      end subroutine exchng2
!===========================================================================

      subroutine exchng2i(tmp1i,tmp2i,tmp3i,tmp4i,tmp5i,tmp6i,tmp7i,tmp8i)
      use mpi
      use var_inc
      implicit none

!                    tmp1i            tmp2i tmp3i                tmp4i
!                    tmp5i            tmp6i tmp7i                tmp8i
!      call exchng2i(ibnodes(:,:,1:2),tmpiR,ibnodes(:,:,lz-1:lz),tmpiL, &
!                    ibnodes(:,1:2,:),tmpiU,ibnodes(:,ly-1:ly,:),tmpiD)

      integer ilenz, ileny, ii, jj, kk
      integer, dimension(lx,ly,2):: tmp1i, tmp2i, tmp3i, tmp4i
      integer, dimension(lx,2,lz):: tmp5i, tmp7i
      integer, dimension(lx,2,-1:lz+2):: tmp5il, tmp6i, tmp7il, tmp8i
!      character (len = 100):: fnm2,fnm3,fnm4
!      character (len = 100):: fnm22,fnm23,fnm24

      integer error,status_array(MPI_STATUS_SIZE,4), req(4)
      integer err1,err2,err3,err4,err5,err6,err7,err8,err9,err10,err11,err12

!      fnm2 = '/ptmp/lwang/debug2D/ibnodes2.2D.dat'
!      fnm3 = '/ptmp/lwang/debug2D/ibnodes3.2D.dat'
!      fnm4 = '/ptmp/lwang/debug2D/ibnodes4.2D.dat'
     
!      fnm22 = '/ptmp/lwang/debug2D/ibnodes22.2D.dat'
!      fnm23 = '/ptmp/lwang/debug2D/ibnodes23.2D.dat'
!      fnm24 = '/ptmp/lwang/debug2D/ibnodes24.2D.dat'

!      open(12, file = trim(fnm2), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(13, file = trim(fnm3), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(14, file = trim(fnm4), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')

!      open(22, file = trim(fnm22), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(23, file = trim(fnm23), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(24, file = trim(fnm24), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')

      ilenz = lx*ly*2
      ileny = lx*(lz+4)*2

!      call MPI_IRECV(tmp2i,ilenz,MPI_INTEGER,mzp,2,MPI_COMM_WORLD,req(1),err1)
!      call MPI_ISEND(tmp1i,ilenz,MPI_INTEGER,mzm,2,MPI_COMM_WORLD,req(2),err1)
!      call MPI_WAITALL(2,req,status_array,err2)
!      call MPI_IRECV(tmp4i,ilenz,MPI_INTEGER,mzm,3,MPI_COMM_WORLD,req(1),err4)
!      call MPI_ISEND(tmp3i,ilenz,MPI_INTEGER,mzp,3,MPI_COMM_WORLD,req(2),err4)
!      call MPI_WAITALL(2,req,status_array,err5)

!Graeme's idea for blocking recv
!     tmp1i            tmp2i tmp3i                tmp4i
!     ibnodes(:,:,1:2) tmpiR ibnodes(:,:,lz-1:lz) tmpiL

      if(myid==0)write(*,*)'I am at second MPI_ISEND (1)'
      call MPI_ISEND(tmp3i,ilenz,MPI_INTEGER,mzp,3,MPI_COMM_WORLD,req(1),err2)

      if(myid==0)write(*,*)'I am at first MPI_ISEND (1)'
      call MPI_ISEND(tmp1i,ilenz,MPI_INTEGER,mzm,2,MPI_COMM_WORLD,req(2),err1)
      
      if(myid==0)write(*,*)'I am at first MPI_RECV (1)'
      call MPI_RECV(tmp2i,ilenz,MPI_INTEGER,mzp,2,MPI_COMM_WORLD,req(3),err3)
      
      if(myid==0)write(*,*)'I am at second MPI_RECV (1)'
      call MPI_RECV(tmp4i,ilenz,MPI_INTEGER,mzm,3,MPI_COMM_WORLD,req(4),err4)

      if(myid==0)write(*,*)'I am at the end of L/R send in exchan21 (1)'
      return

      call MPI_BARRIER(MPI_COMM_WORLD,err5)

      tmp5il(:,1,-1) = tmp4i(:,1,1)
      tmp5il(:,1,0) = tmp4i(:,1,2)
      tmp5il(:,2,-1) = tmp4i(:,2,1)
      tmp5il(:,2,0) = tmp4i(:,2,2)

      tmp5il(:,:,1:lz) = tmp5i(:,:,:)

      tmp5il(:,1,lz+1) = tmp2i(:,1,1)
      tmp5il(:,1,lz+2) = tmp2i(:,1,2)
      tmp5il(:,2,lz+1) = tmp2i(:,2,1)
      tmp5il(:,2,lz+2) = tmp2i(:,2,2)

      tmp7il(:,1,-1) = tmp4i(:,ly-1,1)
      tmp7il(:,1,0) = tmp4i(:,ly-1,2)
      tmp7il(:,2,-1) = tmp4i(:,ly,1)
      tmp7il(:,2,0) = tmp4i(:,ly,2)

      tmp7il(:,:,1:lz) = tmp7i(:,:,:)

      tmp7il(:,1,lz+1) = tmp2i(:,ly-1,1)
      tmp7il(:,1,lz+2) = tmp2i(:,ly-1,2)
      tmp7il(:,2,lz+1) = tmp2i(:,ly,1)
      tmp7il(:,2,lz+2) = tmp2i(:,ly,2)

!      do ii = 1,nx
!      do jj = 1,2
!      do kk = 25,28 
!       if(istep==630 .and. myid==35 .and. ii==1 .and. jj==1 .and. kk==25)write(*,*)'Writing out here (3)'
!       if(istep==631 .and. myid==35 .and. ii==1 .and. jj==1 .and. kk==25)write(*,*)'Writing out here (3)'
!       if(istep==630 .and. myid==35)write(12,211)tmp7i(ii,jj,kk),ii,jj,kk
!       if(istep==631 .and. myid==35)write(22,211)tmp7i(ii,jj,kk),ii,jj,kk
!      enddo
!      enddo
!      enddo

!      do ii = 1,nx
!      do jj = 1,2
!      do kk = 25,28 
!       if(istep==630 .and. myid==35 .and. ii==1 .and. jj==1 .and. kk==25)write(*,*)'Writing out here (4)'
!       if(istep==631 .and. myid==35 .and. ii==1 .and. jj==1 .and. kk==25)write(*,*)'Writing out here (4)'
!       if(istep==630 .and. myid==35)write(13,211)tmp7il(ii,jj,kk),ii,jj,kk
!       if(istep==631 .and. myid==35)write(23,211)tmp7il(ii,jj,kk),ii,jj,kk
!      enddo
!      enddo
!      enddo

!      close(12)
!      close(13)
!      close(23)
!      close(24)

!      call MPI_IRECV(tmp6i,ileny,MPI_INTEGER,myp,4,MPI_COMM_WORLD,req(1),err7)
!      call MPI_ISEND(tmp5il,ileny,MPI_INTEGER,mym,4,MPI_COMM_WORLD,req(2),err7)
!      call MPI_WAITALL(2,req,status_array,err8)
!      call MPI_IRECV(tmp8i,ileny,MPI_INTEGER,mym,5,MPI_COMM_WORLD,req(1),err10)
!      call MPI_ISEND(tmp7il,ileny,MPI_INTEGER,myp,5,MPI_COMM_WORLD,req(2),err10)
!      call MPI_WAITALL(2,req,status_array,err11)

      call MPI_BARRIER(MPI_COMM_WORLD,err6)

      if(myid==0)write(*,*)'I am in exchng2i(2),ileny,size(tmp5il)',  &
                       'size(tmp7il)',ileny,size(tmp5il),size(tmp7il)

      call MPI_BARRIER(MPI_COMM_WORLD,err11)

!Graeme's idea for blocking recv
      call MPI_ISEND(tmp6i,ileny,MPI_INTEGER,myp,4,MPI_COMM_WORLD,req(1),err7)
      call MPI_ISEND(tmp8i,ileny,MPI_INTEGER,mym,5,MPI_COMM_WORLD,req(2),err8)
      call MPI_RECV(tmp5il,ileny,MPI_INTEGER,myp,5,MPI_COMM_WORLD,req(3),err9)
      call MPI_RECV(tmp7il,ileny,MPI_INTEGER,mym,4,MPI_COMM_WORLD,req(4),err10)

      call MPI_BARRIER(MPI_COMM_WORLD,err12)

      if(myid==0)write(*,*)'I am in exchng2i (3)'

!      do ii = 1,nx
!      do jj = 1,2
!      do kk = 25,28 
!       if(istep==630 .and. myid==36  .and. ii==1 .and. jj==1 .and. kk==25)write(*,*)'Writing out here (5)'
!       if(istep==631 .and. myid==36  .and. ii==1 .and. jj==1 .and. kk==25)write(*,*)'Writing out here (5)'
!       if(istep==630 .and. myid==36)write(14,211)tmp8i(ii,jj,kk)
!       if(istep==631 .and. myid==36)write(24,211)tmp8i(ii,jj,kk)
!      enddo
!      enddo
!      enddo

!      close(14)
!      close(24)

211   format(2x,4i5)
      end subroutine exchng2i
!===========================================================================
      subroutine exchng2iNew(tmpSendL,tmpRecvR,tmpSendR,tmpRecvL,tmpSendD,tmpRecvUl,tmpSendU,tmpRecvDl)
      use mpi
      use var_inc
      implicit none

!      call exchng2iNew(ibnodes(:,:,1:2),tmpiR,ibnodes(:,:,lz-1:lz),tmpiL, &
!                       ibnodes(:,1:2,:),tmpiU,ibnodes(:,ly-1:ly,:),tmpiD)

!      allocate (tmpiL(lx,ly,-1:0))
!      allocate (tmpiR(lx,ly,lz+1:lz+2))
!      allocate (tmpiU(lx,ly+1:ly+2,-1:lz+2))
!      allocate (tmpiD(lx,-1:0,-1:lz+2))

      integer ileny, ilenz, i, j, ii,jj,kk,k
      integer,dimension(lx,ly,2):: tmpSendL,tmpSendR,tmpRecvL,tmpRecvR
      integer,dimension(lx,2,lz):: tmpSendU,tmpSendD
      integer,dimension(lx,2,-1:lz+2):: tmpRecvUl,tmpRecvDl,tmpSendUl,tmpSendDl
!     integer err1,err2,err3,err4,err5,err6,err7,err8
!     integer reqSendL,reqSendR,reqSendU,reqSendD,reqRecvL,reqRecvR,reqRecvU,reqRecvD
      integer status_array(MPI_STATUS_SIZE,4),req(4)

!      character (len = 100):: fnm2,fnm3,fnm4
!      character (len = 100):: fnm22,fnm23,fnm24

!      fnm2 = '/ptmp/lwang/debug2D/ibnodes2.2D.dat'
!      fnm3 = '/ptmp/lwang/debug2D/ibnodes3.2D.dat'
!      fnm4 = '/ptmp/lwang/debug2D/ibnodes4.2D.dat'

!      fnm22 = '/ptmp/lwang/debug2D/ibnodes22.2D.dat'
!      fnm23 = '/ptmp/lwang/debug2D/ibnodes23.2D.dat'
!      fnm24 = '/ptmp/lwang/debug2D/ibnodes24.2D.dat'

!     open(12, file = trim(fnm2), status = 'unknown',                 &
!                form = 'formatted', position = 'append')
!     open(13, file = trim(fnm3), status = 'unknown',                 &
!                form = 'formatted', position = 'append')
!     open(14, file = trim(fnm4), status = 'unknown',                 &
!                form = 'formatted', position = 'append')

!      open(22, file = trim(fnm22), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(23, file = trim(fnm23), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')
!      open(24, file = trim(fnm24), status = 'unknown',                 &
!                 form = 'formatted', position = 'append')

      ilenz = lx*ly*2
      ileny = lx*(lz+4)*2

!       do ii = 1,nx
!       do jj = 1,ly 
!       do kk = lz+1,lz+2 
!       if(istep==630 .and. myid==38 .and. ii==1 .and. jj==224 .and. kk==1)write(*,*)'Writing out here (3)'
!       if(istep==631 .and. myid==35 .and. ii==1 .and. jj==1 .and. kk==25)write(*,*)'Writing out here (3)'
!       if(istep==630 .and. myid==52)write(12,211)tmpRecvR(ii,jj,kk),ii,jj,kk
!       if(istep==631 .and. myid==35)write(22,211)tmpSendU(ii,jj,kk),ii,jj,kk
!       enddo
!       enddo
!       enddo

!      close(12)

      call MPI_IRECV(tmpRecvR,ilenz,MPI_INTEGER,mzp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmpRecvL,ilenz,MPI_INTEGER,mzm,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmpSendL,ilenz,MPI_INTEGER,mzm,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmpSendR,ilenz,MPI_INTEGER,mzp,1,MPI_COMM_WORLD,req(4),ierr)
      call MPI_WAITALL(4,req,status_array,ierr)

!       do ii = 1,nx
!       do jj = 1,ly
!       do kk = lz+1,lz+2 
!       if(istep==630 .and. myid==35 .and. ii==1 .and. jj==1 .and. kk==25)write(*,*)'Writing out here (4)'
!       if(istep==631 .and. myid==35 .and. ii==1 .and. jj==1 .and. kk==25)write(*,*)'Writing out here (4)'
!       if(istep==630 .and. myid==52)write(13,211)tmpRecvR(ii,jj,kk),ii,jj,kk
!       if(istep==631 .and. myid==35)write(23,211)tmpSendUl(ii,jj,kk),ii,jj,kk
!       enddo
!       enddo
!       enddo

!      close(13)

!     do j=1,2
!      do i = 1,2
!       tmpSendUl(:,i,j-2) = tmpRecvL(:,ly-2+i,j)
!      enddo
!     enddo

!     do j=1,2
!      do i=1,2
!       tmpSendUl(:,i,lz+j) = tmpRecvR(:,ly-2+i,j)
!      enddo
!     enddo

        tmpSendUl(:,1:2,-1:0) = tmpRecvL(:,ly-1:ly,1:2)
        tmpSendUl(:,1:2,lz+1:lz+2) = tmpRecvR(:,ly-1:ly,1:2)
        tmpSendUl(:,:,1:lz) = tmpSendU(:,:,:)

!     do j=1,2
!      do i=1,2
!       tmpSendDl(:,i,j-2) = tmpRecvL(:,i,j)
!      enddo
!     enddo

!     do j=1,2
!      do i=1,2
!       tmpSendDl(:,i,lz+j) = tmpRecvR(:,i,j)
!      enddo
!     enddo

        tmpSendDl(:,1:2,-1:0) = tmpRecvL(:,1:2,1:2)
        tmpSendDl(:,1:2,lz+1:lz+2) = tmpRecvR(:,1:2,1:2)
        tmpSendDl(:,:,1:lz) = tmpSendD(:,:,:)


!      if(istep.eq.630)then
!      if(myid.eq.1) write(33,*)' before U',((j,k,tmpSendDl(5,j,k),j=1,2),k=1,lz)
!      if(myid.eq.7) write(34,*)' before D',((j,k,tmpSendUl(5,j,k),j=1,2),k=1,lz)

!      if(myid.eq.57) write(35,*)' before UL',((j,k,tmpSendDl(5,j,k),j=1,2),k=lz-1,lz)
!      if(myid.eq.63) write(36,*)' before DL',((j,k,tmpSendUl(5,j,k),j=1,2),k=lz-1,lz)
!      if(myid.eq.9) write(37,*)' before UR',((j,k,tmpSendDl(5,j,k),j=1,2),k=1,2)
!      if(myid.eq.15) write(38,*)' before DR',((j,k,tmpSendUl(5,j,k),j=1,2),k=1,2)

!      if(myid.eq.56) close(31)
!      if(myid.eq.8) close(32)
!      if(myid.eq.1) close(33)
!      if(myid.eq.7) close(34)
!      if(myid.eq.57) close(35)
!      if(myid.eq.63) close(36)
!      if(myid.eq.9) close(37)
!      if(myid.eq.15) close(38)
!      end if

      call MPI_IRECV(tmpRecvUl,ileny,MPI_INTEGER,myp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmpRecvDl,ileny,MPI_INTEGER,mym,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmpSendDl,ileny,MPI_INTEGER,mym,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmpSendUl,ileny,MPI_INTEGER,myp,1,MPI_COMM_WORLD,req(4),ierr)
      call MPI_WAITALL(4,req,status_array,ierr)

!     do ii = 1,nx
!     do jj = 1,2
!     do kk = 25,28
!      if(istep==630 .and. myid==36  .and. ii==1 .and. jj==1 .and. kk==25)write(*,*)'Writing out here (5)'
!       if(istep==631 .and. myid==36  .and. ii==1 .and. jj==1 .and. kk==25)write(*,*)'Writing out here (5)'
!      if(istep==630 .and. myid==36)write(14,211)tmpRecvDl(ii,jj,kk)
!       if(istep==631 .and. myid==36)write(24,211)tmpRecvDl(ii,jj,kk)
!     enddo
!     enddo
!     enddo

!     close(14)
!      close(24)

!211   format(2x,4i5)
      end subroutine exchng2iNew      
!===========================================================================      
      subroutine exchng3(tmp1,tmp3,tmp5,tmp7,&
            tmp1i,tmp3i,tmp5i,tmp7i)
      use mpi
      use var_inc
      implicit none

!                  1              2    3                 4
!     call exchng3(tmpfL(:,:,:,0),tmp3,tmpfR(:,:,:,lz+1),tmp4, &
!                  5              6    7                 8
!                  tmpfD(:,:,0,:),tmp5,tmpfU(:,:,ly+1,:),tmp6)


      integer error, ilenz, ileny
      real, dimension(0:npop-1,lx,ly):: tmp1, tmp2, tmp3, tmp4
      real, dimension(0:npop-1,lx,lz):: tmp6, tmp8
      real, dimension(0:npop-1,lx,0:lz+1):: tmp5, tmp7, tmp6e, tmp8e
      integer, dimension(0:npop-1,lx,ly):: tmp1i, tmp2i, tmp3i, tmp4i
      integer, dimension(0:npop-1,lx,lz):: tmp6i, tmp8i
      integer, dimension(0:npop-1,lx,0:lz+1):: tmp5i, tmp7i, tmp6ie, tmp8ie

      integer status_array(MPI_STATUS_SIZE,4), req(4)

      ilenz = npop*lx*ly
      ileny = npop*lx*(lz+2)
!Change: this needs to be in reverse order, UD first, then LR
      call MPI_IRECV(tmp6e,ileny,MPI_REAL8,myp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp8e,ileny,MPI_REAL8,mym,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp5,ileny,MPI_REAL8,mym,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp7,ileny,MPI_REAL8,myp,1,MPI_COMM_WORLD,req(4),ierr)
      call MPI_WAITALL(4,req,status_array,ierr)

      call MPI_IRECV(tmp6ie,ileny,MPI_INTEGER,myp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp8ie,ileny,MPI_INTEGER,mym,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp5i,ileny,MPI_INTEGER,mym,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp7i,ileny,MPI_INTEGER,myp,1,MPI_COMM_WORLD,req(4),ierr)
      call MPI_WAITALL(4,req,status_array,ierr)

      tmp6 = tmp6e(:,:,1:lz)
      tmp6i = tmp6ie(:,:,1:lz)
      tmp8 = tmp8e(:,:,1:lz)
      tmp8i = tmp8ie(:,:,1:lz)
! Now corners
      where(tmp6ie(:,:,0) == 1) tmp1(:,:,ly) = tmp6e(:,:,0) 
      where(tmp6ie(:,:,lz+1) == 1) tmp3(:,:,ly) = tmp6e(:,:,lz+1) 
      where(tmp8ie(:,:,0) == 1) tmp1(:,:,1) = tmp8e(:,:,0) 
      where(tmp8ie(:,:,lz+1) == 1) tmp3(:,:,1) = tmp8e(:,:,lz+1) 

      where(tmp6ie(:,:,0) == 1) tmp1i(:,:,ly) = tmp6ie(:,:,0) 
      where(tmp6ie(:,:,lz+1) == 1) tmp3i(:,:,ly) = tmp6ie(:,:,lz+1) 
      where(tmp8ie(:,:,0) == 1) tmp1i(:,:,1) = tmp8ie(:,:,0) 
      where(tmp8ie(:,:,lz+1) == 1) tmp3i(:,:,1) = tmp8ie(:,:,lz+1) 

! in: tmp1 and tmp3;    out:tmp2,tmp4
      call MPI_IRECV(tmp2,ilenz,MPI_REAL8,mzp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp4,ilenz,MPI_REAL8,mzm,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp1,ilenz,MPI_REAL8,mzm,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp3,ilenz,MPI_REAL8,mzp,1,MPI_COMM_WORLD,req(4),ierr)
      call MPI_WAITALL(4,req,status_array,ierr)

      call MPI_IRECV(tmp2i,ilenz,MPI_INTEGER,mzp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp4i,ilenz,MPI_INTEGER,mzm,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp1i,ilenz,MPI_INTEGER,mzm,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp3i,ilenz,MPI_INTEGER,mzp,1,MPI_COMM_WORLD,req(4),ierr)
      call MPI_WAITALL(4,req,status_array,ierr)

      where(tmp2i == 1) f(:,:,:,lz) = tmp2  
      where(tmp4i == 1) f(:,:,:,1) = tmp4    
      where(tmp6i == 1) f(:,:,ly,:) = tmp6
      where(tmp8i == 1) f(:,:,1,:) = tmp8

      end subroutine exchng3
!===========================================================================

      subroutine exchng3i(tmp1i,tmp2i,tmp3i,tmp4i,tmp5i,tmp6i,tmp7i,tmp8i)
      use mpi
      use var_inc
      implicit none

!                    1    2     3    4     5    6     7    8
!      call exchng3i(if9L,tmp3i,if9R,tmp4i,if9D,tmp5i,if9U,tmp6)

      integer ileny, ilenz
      integer, dimension(0:npop-1,lx,ly):: tmp1i, tmp2i, tmp3i, tmp4i
      integer, dimension(0:npop-1,lx,lz):: tmp5i, tmp7i
      integer, dimension(0:npop-1,lx,0:lz+1):: tmp5il, tmp6i, tmp7il, tmp8i

      integer error, status_array(MPI_STATUS_SIZE,4), req(4)

      ilenz = npop*lx*ly
      ileny = npop*lx*(lz+2) 

      call MPI_IRECV(tmp2i,ilenz,MPI_INTEGER,mzp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp4i,ilenz,MPI_INTEGER,mzm,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp1i,ilenz,MPI_INTEGER,mzm,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp3i,ilenz,MPI_INTEGER,mzp,1,MPI_COMM_WORLD,req(4),ierr)
      call MPI_WAITALL(4,req,status_array,ierr)

      tmp5il(:,:,0) = tmp2i(:,:,1)
      tmp5il(:,:,1:lz) = tmp5i(:,:,:)
      tmp5il(:,:,lz+1) = tmp4i(:,:,1)

      tmp7il(:,:,0) = tmp2i(:,:,ly)
      tmp7il(:,:,1:lz) = tmp7i(:,:,:)
      tmp7il(:,:,lz+1) = tmp4i(:,:,ly)



      end subroutine exchng3i
!===========================================================================

      subroutine beads_lubforce
      use mpi 
      use var_inc
      implicit none

      integer ii, jj, id, ilen,ip
      real dist0, dxij, dyij, dzij, dist, hgap
      real xc, yc, zc, xb, yb, zb, xd 
      real fp, fpc, fpd, dist1
      real, dimension(3,npart)  :: flubp0

      dist0 = 2.0*rad + mingap
!     fpc = volp*rhog
      fpc = volp*rho0*g_lbm

      flubp0 = 0.0

! First particle-particle lubrication add-on
      do ii = 1,nps
        xc = yp(1,ii)
        yc = yp(2,ii)
        zc = yp(3,ii)
        id = ipglb(ii)
                

        if(id < npart)then

        do jj = id+1,npart
          xb = ypglb(1,jj)
          yb = ypglb(2,jj)
          zb = ypglb(3,jj)

! use the nearest particle center instead of the real center
         ! if((xb - xc) > real(nxh)) xb = xb - real(nx)
         ! if((xb - xc) < -real(nxh)) xb = xb + real(nx)

          if((yb - yc) > real(nyh)) yb = yb - real(ny)
          if((yb - yc) < -real(nyh)) yb = yb + real(ny)

          if((zb - zc) > real(nzh)) zb = zb - real(nz)
          if((zb - zc) < -real(nzh)) zb = zb + real(nz)


          dxij = xc - xb
          dyij = yc - yb
          dzij = zc - zb
          dist = sqrt(dxij*dxij + dyij*dyij + dzij*dzij)
          hgap = dist - 2.0*rad

          if(dist < dist0)then
! adopt from Feng & Michaelides (2005) JCP 202, pp 20-51, eqn(28)
            dist1 = 0.0
            if(hgap < 0.0) dist1 = -hgap

            fp = (((dist0 - dist) / mingap)**2 / stf0                  &
               + dist1 / mingap / stf1)*fpc
          fpd = fp / dist 
          
          flubp0(1,id) = flubp0(1,id) + fpd*dxij          
          flubp0(2,id) = flubp0(2,id) + fpd*dyij          
          flubp0(3,id) = flubp0(3,id) + fpd*dzij          

          flubp0(1,jj) = flubp0(1,jj) - fpd*dxij          
          flubp0(2,jj) = flubp0(2,jj) - fpd*dyij          
          flubp0(3,jj) = flubp0(3,jj) - fpd*dzij          
!         write(*,*)'id,jj,fp=',id,jj,fp
          end if
          
        end do
        end if
      end do
      
! Next particle-wall lubrication add-on
      dist0 = mingap_w + rad

      do ii = 1,nps
      id = ipglb(ii)

      xd = yp(1,ii)

! bottom wall first
      dist = xd
! the truncation distance of lubforce from walls
      hgap = dist - rad
      if (dist < dist0) then 
!when the distance between particle and wall is less than mingap_w 
! the wall start to give a lubforce to particles
        dist1 = 0.0
        if(hgap < 0.0)  dist1 = -hgap
            fp = (((dist0 - dist) / mingap_w)**2 / stf0_w            &
               + dist1 / mingap_w / stf1_w)*fpc
        flubp0(1,id) = flubp0(1,id) + fp
!       write(*,*)'id,fp=',id,fp
       end if

! top wall
      dist = real(nx) - xd
! the truncation distance of lubforce from walls
      hgap = dist - rad
      if (dist < dist0) then
!when the distance between particle and wall is less than mingap_w &
! the wall start to give a lubforce to particles
        dist1 = 0.0
        if(hgap < 0.0) dist1 = -hgap
            fp = (((dist0 - dist) / mingap_w)**2 / stf0_w                  &
               + dist1 / mingap_w / stf1_w)*fpc
          flubp0(1,id) = flubp0(1,id) - fp
!         write(*,*)'id,fp=',id,fp
       end if
 
      end do

! collect info. for flubp
      ilen = 3*npart

      call MPI_ALLREDUCE(flubp0,flubp,ilen,MPI_REAL8,MPI_SUM,           &
                         MPI_COMM_WORLD,ierr)      

!      do ii=1,npart
!      if(myid==0)write(*,*)'istep,i,flubp(:,ii)',istep,ii,flubp(:,ii)
!      enddo

      end subroutine beads_lubforce
!===========================================================================

      subroutine beads_move
      use mpi 
      use var_inc
      implicit none

      integer ip, id, ilen
      real dt
      real, dimension(3):: dwdtn, domgdtn
      real, dimension(3,npart):: wp0, omgp0

      dt = 1.0

      wp0 = 0.0
      omgp0 = 0.0

      if(istep == irelease) thetap = 0.0

! save as previous values for wp and omgp 
      wpp = wp
      omgpp = omgp

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do ip = 1,nps
        id = ipglb(ip)


! update particle translational velocity
        forcep(:,ip) = fHIp(:,id) + flubp(:,id)
        


!       if(myid==0 .and. ip==1)then
!       write(*,*)'ip',ip,'fHIp(:,id)',fHIp(:,id),'flubp(:,id)',flubp(:,id),'forcep(:,ip)',forcep(:,ip)
!       endif


! if gravity is desired in "-x" dir
!       forcep(1,ip) = forcep(1,ip) - volp*rhog
! Add the pressure-gradient force for the particle in "+y" direction
        forcep(2,ip) = forcep(2,ip) + volp*force_in_y


        if(istep == irelease) forcepp(:,ip) = forcep(:,ip)
        dwdtn = 0.5*(forcep(:,ip) + forcepp(:,ip)) / amp

        if(istep == irelease) dwdt(:,ip) = dwdtn
        wp(:,id) = wpp(:,id) + 0.5*dt*(dwdt(:,ip) + dwdtn)


! update particle position
        yp(:,ip) = yp(:,ip) + 0.5*dt*(wp(:,id) + wpp(:,id))        
 
! update particle angular velocity
        if(istep == irelease) torqpp(:,ip) = torqp(:,id)
        domgdtn = 0.5*(torqp(:,id) + torqpp(:,ip)) / aip

        if(istep == irelease) domgdt(:,ip) = domgdtn
        omgp(:,id) = omgpp(:,id) + 0.5*dt*(domgdt(:,ip) + domgdtn)


! update particle angular displacement
        thetap(:,ip) = thetap(:,ip) + 0.5*dt*(omgp(:,id) + omgpp(:,id))

! update d(wp)/dt, d(omgp)/dt
!        dwdt(:,ip) = (wp(:,id) - wpp(:,id)) / dt  
        dwdt(:,ip) = dwdtn 

!        domgdt(:,ip) = (omgp(:,id) - omgpp(:,id)) / dt
        domgdt(:,ip) = domgdtn  

! save as previous values
        forcepp(:,ip) = forcep(:,ip)

        torqpp(:,ip) = torqp(:,id)

! update wp0 and mgp0 
        wp0(:,id) = wp(:,id)

        omgp0(:,id) = omgp(:,id)

      end do

! save as previous values for ibnodes and isnodes
      ibnodes0 = ibnodes
      isnodes0 = isnodes

! update info. for wp and omgp
      ilen = 3*npart


!      write(*,*)'myid,wp0',myid,wp0
!      write(*,*)'myid,omgp0',myid,omgp0

      call MPI_ALLREDUCE(wp0,wp,ilen,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(omgp0,omgp,ilen,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!       if(myid==0)then
!       write(*,*)'wp =',wp
!       write(*,*)'omgp = ',omgp
!       endif

      end subroutine beads_move
!===========================================================================

      Subroutine beads_redistribute
      use mpi 
      use var_inc
      implicit none

      integer ip, n2l, n2r, n2d, n2dl, n2dr, n2u, n2ul, n2ur, ns, ntmp0 
      integer npsflag, npsflagt, npst, id, ovlpcnt, ckpw1, ckpw2
      integer, dimension(msize):: id2l, id2r, id2u, id2d, ids
      integer, dimension(msize):: id2ul, id2ur, id2dl, id2dr
      integer, dimension(msize):: tmp2li, tmp2ri, tmp2ui, tmp2di, tmpsi, tmp0i
      integer, dimension(msize):: tmp2dli, tmp2dri, tmp2uli, tmp2uri

      real xc, yc, zc
      real, dimension(3,npart):: yp0 
      real, dimension(21,msize):: tmp2l, tmp2r, tmp2d, tmp2u, tmps, tmp0
      real, dimension(21,msize):: tmp2dl, tmp2dr, tmp2ul, tmp2ur
    
      ckpw1 = 0
      ckpw2 = 0
      n2l = 0
      n2r = 0
      n2u = 0
      n2d = 0
      n2ul = 0
      n2ur = 0
      n2dl = 0
      n2dr = 0
      ns = 0

      tmp2li = 0
      tmp2ri = 0
      tmp2ui = 0.0
      tmp2uli = 0.0
      tmp2uri = 0.0
      tmp2di = 0.0
      tmp2dli = 0.0
      tmp2dri = 0.0
      tmpsi = 0

      tmp2l = 0.0
      tmp2r = 0.0
      tmp2d = 0.0
      tmp2dl = 0.0
      tmp2dr = 0.0
      tmp2u = 0.0
      tmp2ul = 0.0
      tmp2ur = 0.0
      tmps = 0.0

      do ip = 1,nps
        yc = yp(2,ip)
        zc = yp(3,ip)
       
 
        if(zc < real(indz*lz))then

          if(yc < real(indy*ly))then             !1
           n2dl = n2dl + 1
           id2dl(n2dl) = ip
          else if(yc >= real((indy+1)*ly))then   !7
           n2ul = n2ul + 1
           id2ul(n2ul) = ip
          else                                   !8
           n2l = n2l + 1
           id2l(n2l) = ip
          endif

        else if(zc >= real((indz+1)*lz))then

          if(yc < real(indy*ly))then             !3
           n2dr = n2dr + 1
           id2dr(n2dr) = ip
          else if(yc >= real((indy+1)*ly))then   !5
           n2ur = n2ur + 1
           id2ur(n2ur) = ip
          else                                   !4
           n2r = n2r + 1
           id2r(n2r) = ip
          endif

        else

          if(yc < real(indy*ly))then             !2
           n2d = n2d + 1
           id2d(n2d) = ip
          else if(yc >= real((indy+1)*ly))then   !6
           n2u = n2u + 1
           id2u(n2u) = ip
          else                                   !stay
           ns = ns + 1
           ids(ns) = ip
          endif

        end if

      end do

      if(n2l > 0) call preswap(n2l,id2l,tmp2l,tmp2li)      !8

      if(n2r > 0) call preswap(n2r,id2r,tmp2r,tmp2ri)      !4

      if(n2u > 0) call preswap(n2u,id2u,tmp2u,tmp2ui)      !6

      if(n2d > 0) call preswap(n2d,id2d,tmp2d,tmp2di)      !2

      if(n2ul > 0) call preswap(n2ul,id2ul,tmp2ul,tmp2uli) !7

      if(n2ur > 0) call preswap(n2ur,id2ur,tmp2ur,tmp2uri) !5

      if(n2dl > 0) call preswap(n2dl,id2dl,tmp2dl,tmp2dli) !1

      if(n2dr > 0) call preswap(n2dr,id2dr,tmp2dr,tmp2dri) !3

      if(ns > 0)  call preswap(ns,ids,tmps,tmpsi)          !stay


      nps = ns
      
      if(nps > 0)then   
        ipglb(1:nps) = tmpsi(1:nps)
        yp(:,1:nps) = tmps(1:3,1:nps)
        thetap(:,1:nps) = tmps(4:6,1:nps)
        dwdt(:,1:nps) = tmps(7:9,1:nps)
        domgdt(:,1:nps) = tmps(10:12,1:nps)
        forcep(:,1:nps) = tmps(13:15,1:nps)
        forcepp(:,1:nps) = tmps(16:18,1:nps)
        torqpp(:,1:nps) = tmps(19:21,1:nps)
      end if 
 
! send tmp2l and tmp2li to the left slab and receive tmp0 and tmp0i from right
      call exchng4(n2l,ntmp0,tmp2li,tmp2l,tmp0i,tmp0,0,-1)

! concatenate tmp0 and tmp0i received from right slab, with nps updated
      call concat(ntmp0,tmp0i,tmp0)

! send tmp2r and tmp2ri to the right slab and recive tmp0 and tmp0i from left
      call exchng4(n2r,ntmp0,tmp2ri,tmp2r,tmp0i,tmp0,0,1)

! concatenate tmp0 and tmp0i received from left slab, with nps updated
      call concat(ntmp0,tmp0i,tmp0)

! send tmp2u and tmp2ui to the up slab and receive tmp0 and tmp0i from down
      call exchng4(n2u,ntmp0,tmp2ui,tmp2u,tmp0i,tmp0,1,0)

! concatenate tmp0 and tmp0i received from down slab, with nps updated
      call concat(ntmp0,tmp0i,tmp0)

! send tmp2d and tmp2di to the down slab and receive tmp0 and tmp0i from up
      call exchng4(n2d,ntmp0,tmp2di,tmp2d,tmp0i,tmp0,-1,0)

! concatenate tmp0 and tmp0i received from up slab, with nps updated
      call concat(ntmp0,tmp0i,tmp0)




! send tmp2ul and tmp2uli to the upper left slab and receive tmp0 and tmp0i from down right 
      call exchng4(n2ul,ntmp0,tmp2uli,tmp2ul,tmp0i,tmp0,1,-1)

! concatenate tmp0 and tmp0i received from ??? slab, with nps updated
      call concat(ntmp0,tmp0i,tmp0)

! send tmp2ur and tmp2uri to the upper right slab and receive tmp0 and tmp0i from down left 
      call exchng4(n2ur,ntmp0,tmp2uri,tmp2ur,tmp0i,tmp0,1,1)

! concatenate tmp0 and tmp0i received from ??? slab, with nps updated
      call concat(ntmp0,tmp0i,tmp0)

! send tmp2dl and tmp2dli to the down left slab and receive tmp0 and tmp0i from up right 
      call exchng4(n2dl,ntmp0,tmp2dli,tmp2dl,tmp0i,tmp0,-1,-1)

! concatenate tmp0 and tmp0i received from ??? slab, with nps updated
      call concat(ntmp0,tmp0i,tmp0)

! send tmp2dr and tmp2dri to the down right slab and receive tmp0 and tmp0i from up left 
      call exchng4(n2dr,ntmp0,tmp2dri,tmp2dr,tmp0i,tmp0,-1,1)

! concatenate tmp0 and tmp0i received from ??? slab, with nps updated
      call concat(ntmp0,tmp0i,tmp0)
      

! check nps
      npsflag = 0
      if(nps > msize)then
        write(*,*) 'fault: nps = ', nps, ' > msize'
         npsflag = 1 
      end if

      call MPI_ALLREDUCE(npsflag,npsflagt,1,MPI_INTEGER,               &
                         MPI_SUM,MPI_COMM_WORLD,ierr)

      if(npsflagt > 0) stop

      call MPI_ALLREDUCE(nps,npst,1,MPI_INTEGER,MPI_SUM,               &
                         MPI_COMM_WORLD,ierr)
      if(npst /= npart)then
        if(myid == 0) write(*,*) 'fault: npst != npart'
        stop
      end if
        
! check if particles overlap with wall
      do ip = 1,nps
        xc = yp(1,ip)        
        yc = yp(2,ip)        
        zc = yp(3,ip)        

        if( (xc+rad) > real(nx)) ckpw2 = ckpw2 + 1
               
        if( (xc-rad) < 0.0) ckpw1 = ckpw1 + 1

        if(yc >= real(ny)) yp(2,ip) = yp(2,ip) - real(ny)
        if(yc < 0.0) yp(2,ip) = yp(2,ip) + real(ny)

        if(zc >= real(nz)) yp(3,ip) = yp(3,ip) - real(nz)
        if(zc < 0.0) yp(3,ip) = yp(3,ip) + real(nz)
      end do

!      open (97, file = trim('part_wall'), status = 'unknown',   &
!      form = 'formatted')
!      write(97,*) istep, ckpw1, ckpw2
!      close (97)

! update yp0 
      yp0 = 0.0

      do ip = 1,nps
        id = ipglb(ip)
        yp0(:,id) = yp(:,ip) 
      end do

! update ypglb  
      ypglbp = ypglb
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(yp0,ypglb,3*npart,MPI_REAL8,MPI_SUM,           &
                         MPI_COMM_WORLD,ierr)      


! overlap check
      ovlpflag = 0
      ovlpflagt = 0

      call ovlpchk

      call MPI_ALLREDUCE(ovlpflag,ovlpflagt,npart,MPI_INTEGER,         &
                         MPI_SUM,MPI_COMM_WORLD,ierr)

      ovlpcnt = count(ovlpflagt /= 0)

      if(ovlpcnt > 0 )then
        if(myid==0) write(*,*)'in beads_redis ',ovlpcnt,                & 
          ' particles overlap at istep = ', istep
        call MPI_FINALIZE(ierr)
        stop
      end if

      end subroutine beads_redistribute
!===========================================================================

      subroutine preswap(n2t,id2t,tmp2t,tmp2ti) 
      use var_inc
      implicit none

      integer ii, ip, n2t 
      integer, dimension(msize):: id2t, tmp2ti
      real, dimension(21,msize):: tmp2t

      do ii = 1,n2t
        ip = id2t(ii)

        tmp2ti(ii) = ipglb(ip)
        tmp2t(1:3,ii) = yp(:,ip)
        tmp2t(4:6,ii) = thetap(:,ip)
        tmp2t(7:9,ii) = dwdt(:,ip)
        tmp2t(10:12,ii) = domgdt(:,ip)
        tmp2t(13:15,ii) = forcep(:,ip)
        tmp2t(16:18,ii) = forcepp(:,ip)
        tmp2t(19:21,ii) = torqpp(:,ip)
      end do

      end subroutine preswap 
!===========================================================================

      subroutine exchng4(n2t,ntmp0,tmp2ti,tmp2t,tmp0i,tmp0,ydir,zdir) 
      use mpi
      use var_inc
      implicit none

      integer n2t, ntmp0, ydir, zdir, ips, ipr  
      integer, dimension(msize):: tmp2ti, tmp0i
      real, dimension(21,msize):: tmp2t, tmp0

      integer status(MPI_STATUS_SIZE)
      integer status_array2(MPI_STATUS_SIZE,2), req2(2)
      integer status_array4(MPI_STATUS_SIZE,4), req4(4)


      if(zdir > 0)then

        if(ydir > 0)then
          ips = mypzp
          ipr = mymzm
        else if(ydir < 0)then
          ips = mymzp
          ipr = mypzm
        else
          ips = mzp
          ipr = mzm
        endif

      else if(zdir < 0)then

        if(ydir > 0)then
          ips = mypzm
          ipr = mymzp
        else if(ydir < 0)then
          ips = mymzm
          ipr = mypzp
        else
          ips = mzm
          ipr = mzp
        endif

      else

        if(ydir > 0)then
          ips = myp
          ipr = mym
        else if(ydir < 0)then
          ips = mym
          ipr = myp
        endif

      endif

! the MPI_SENDRECV() below does not pass compiling. CISL stuff suggested 
! use an integer dummy variable of value 1 to replace the
! hardcoding value of 1 in the sendcount, and then the code compiles 
! successfully. Also, the sendtag and recvtag must match identically, 
! or the code  would deadlock.
!      call MPI_SENDRECV(n2t,1,MPI_INTEGER,ips,0,                       &
!                        ntmp0,1,MPI_INTEGER,ipr,0,                     &
!                        MPI_COMM_WORLD,status,ierr)

      call MPI_ISEND(n2t,1,MPI_INTEGER,ips,0,MPI_COMM_WORLD,req2(1),ierr)
      call MPI_IRECV(ntmp0,1,MPI_INTEGER,ipr,0,MPI_COMM_WORLD,req2(2),ierr)
      call MPI_WAITALL(2,req2,status_array2,ierr)

      if(n2t > 0)then
        call MPI_ISEND(tmp2ti,n2t,MPI_INTEGER,ips,1,                   &
                       MPI_COMM_WORLD,req4(1),ierr)
        call MPI_ISEND(tmp2t,21*n2t,MPI_REAL8,ips,2,                    &
                       MPI_COMM_WORLD,req4(2),ierr)
      end if

      if(ntmp0 > 0)then
        call MPI_IRECV(tmp0i,ntmp0,MPI_INTEGER,ipr,1,                  &
                       MPI_COMM_WORLD,req4(3),ierr)
        call MPI_IRECV(tmp0,21*ntmp0,MPI_REAL8,ipr,2,                   &
                       MPI_COMM_WORLD,req4(4),ierr)
      end if

      if(n2t > 0 .and. ntmp0 > 0)                                      &
        call MPI_WAITALL(4,req4,status_array4,ierr)
      if(n2t > 0 .and. ntmp0 == 0)then
        call MPI_WAIT(req4(1),status,ierr)
        call MPI_WAIT(req4(2),status,ierr)
      end if
      if(n2t == 0 .and. ntmp0 > 0)then
        call MPI_WAIT(req4(3),status,ierr)
        call MPI_WAIT(req4(4),status,ierr)
      end if

      end subroutine exchng4 
!===========================================================================
      subroutine concat(ntmp0,tmp0i,tmp0)
      use var_inc
      implicit none

      integer ntmp0, nps1, nps2 
      integer, dimension(msize):: tmp0i
      real, dimension(21,msize):: tmp0

      if(ntmp0 > 0)then
        nps1 = nps + 1
        nps2 = nps + ntmp0
        nps = nps2

        ipglb(nps1:nps2) = tmp0i(1:ntmp0)
        yp(:,nps1:nps2) = tmp0(1:3,1:ntmp0)
        thetap(:,nps1:nps2) = tmp0(4:6,1:ntmp0)
        dwdt(:,nps1:nps2) = tmp0(7:9,1:ntmp0)
        domgdt(:,nps1:nps2) = tmp0(10:12,1:ntmp0)
        forcep(:,nps1:nps2) = tmp0(13:15,1:ntmp0)
        forcepp(:,nps1:nps2) = tmp0(16:18,1:ntmp0)
        torqpp(:,nps1:nps2) = tmp0(19:21,1:ntmp0)
      end if

      end subroutine concat
!===========================================================================

      subroutine beads_filling
      use mpi
      use var_inc
      implicit none

      integer id, ix, iy, iz, ipop, ipmx, ii, nghb
      integer i, j, k, ip1, jp1, kp1
      integer ibp1,ib0p1  

      real xc, yc, zc, xpnt, ypnt, zpnt
      real w1, w2, w3, omg1, omg2, omg3
      real aa, bb, cc, ddt, ddt0, ddt1  
      real xp1, yp1, zp1, xp2, yp2, zp2
      real xx0, yy0, zz0, prod0, prod
      real rho9, u9, v9, w9

      real, dimension(0:npop-1):: f9, feq9 
      real,allocatable,dimension(:,:,:):: tmpfL, tmpfR    
      real,allocatable,dimension(:,:,:):: tmpfU, tmpfD
      
      integer,allocatable,dimension(:,:):: tmpiL, tmpiR
      integer,allocatable,dimension(:,:):: tmpiU, tmpiD 

      integer,allocatable,dimension(:,:):: tmpiL0, tmpiR0   
      integer,allocatable,dimension(:,:):: tmpiU0, tmpiD0

      allocate(tmpfL(0:npop-1,lx,ly))   
      allocate(tmpfR(0:npop-1,lx,ly))   
      allocate(tmpfU(0:npop-1,lx,0:lz+1))
      allocate(tmpfD(0:npop-1,lx,0:lz+1))

      allocate(tmpiL(lx,ly))
      allocate(tmpiR(lx,ly))
      allocate(tmpiU(lx,0:lz+1))
      allocate(tmpiD(lx,0:lz+1))

      allocate(tmpiL0(lx,ly))
      allocate(tmpiR0(lx,ly))
      allocate(tmpiU0(lx,0:lz+1))
      allocate(tmpiD0(lx,0:lz+1))


      call exchng5(f(:,:,:,1),tmpfR,f(:,:,:,lz),tmpfL,          &
                   f(:,:,1,:),tmpfU,f(:,:,ly,:),tmpfD)

      call exchng5i(ibnodes(:,:,1),tmpiR,ibnodes(:,:,lz),tmpiL, &
                    ibnodes(:,1,:),tmpiU,ibnodes(:,ly,:),tmpiD)


      call exchng5i(ibnodes0(:,:,1),tmpiR0,ibnodes0(:,:,lz),tmpiL0,  &
                    ibnodes0(:,1,:),tmpiU0,ibnodes0(:,ly,:),tmpiD0)


      do k = 1,lz
      do j = 1,ly
      do i = 1,lx
      IF(ibnodes0(i,j,k) > 0 .and. ibnodes(i,j,k) < 0)THEN
!     write(*,*)'Start istep,i,j,k='
      id = isnodes0(i,j,k)

      xc = ypglbp(1,id)*1.d0      
      yc = ypglbp(2,id)*1.d0    
      zc = ypglbp(3,id)*1.d0   

      xpnt = dfloat(i) - 0.5d0 
      ypnt = dfloat(j) - 0.5d0 + dfloat(indy*ly)   
      zpnt = dfloat(k) - 0.5d0 + dfloat(indz*lz)

! use the nearest particle center instead of the real center
!      if((xc - xpnt) > dfloat(nxh)) xc = xc - dfloat(nx)
!      if((xc - xpnt) < -dfloat(nxh)) xc = xc + dfloat(nx)

      if((yc - ypnt) > dfloat(nyh)) yc = yc - dfloat(ny)
      if((yc - ypnt) < -dfloat(nyh)) yc = yc + dfloat(ny)

      if((zc - zpnt) > dfloat(nzh)) zc = zc - dfloat(nz)
      if((zc - zpnt) < -dfloat(nzh)) zc = zc + dfloat(nz)

      w1 = -0.5d0*(wp(1,id)*1.d0 + wpp(1,id)*1.d0)
      w2 = -0.5d0*(wp(2,id)*1.d0 + wpp(2,id)*1.d0)
      w3 = -0.5d0*(wp(3,id)*1.d0 + wpp(3,id)*1.d0)
      omg1 = -0.5d0*(omgp(1,id)*1.d0 + omgpp(1,id)*1.d0)
      omg2 = -0.5d0*(omgp(2,id)*1.d0 + omgpp(2,id)*1.d0)
      omg3 = -0.5d0*(omgp(3,id)*1.d0 + omgpp(3,id)*1.d0)

      aa = w1*w1 + w2*w2 + w3*w3
      bb = (xpnt - xc)*w1 + (ypnt - yc)*w2 + (zpnt -zc)*w3 
      cc = (xpnt - xc)**2 + (ypnt - yc)**2 + (zpnt -zc)**2 - (rad*1.d0)**2 

      ddt0 = bb/aa 
      ddt1 = sqrt(ddt0**2 - cc/aa) 
      ddt = -ddt0 + ddt1  
  
!     if(ddt < 0.d0 .or. ddt > 1.d0)then
!        write(*,*) 'fault: ddt = ', ddt, ' @ ibnodes0(', i,            &
!                   ', ', j, ', ', k, ')'
!        write(*,*) 'ddtt = ', ddtt
!        write(*,*) 'ddt0 = ', ddt0, ' ddt1 = ', ddt1, ' aa = ', aa, ' bb = ', bb, ' cc = ', cc  
!        stop
!     end if
      if(ddt < 0.d0) ddt = 0.d0
      if(ddt > 1.d0) ddt = 1.d0

      xp1 = xpnt + w1*ddt
      yp1 = ypnt + w2*ddt
      zp1 = zpnt + w3*ddt

! (xp2, yp2, zp2) is the point on the particle surface. It is THROUGH this
! point the previous solid node (xpnt, ypnt, zpnt) moves to fluid region.
      xp2 = xp1 + (omg2*(zp1-zc) - omg3*(yp1-yc))*ddt
      yp2 = yp1 + (omg3*(xp1-xc) - omg1*(zp1-zc))*ddt
      zp2 = zp1 + (omg1*(yp1-yc) - omg2*(xp1-xc))*ddt

      xx0 = xp2 - xc
      yy0 = yp2 - yc
      zz0 = zp2 - zc

! Lallemand and Luo, JCP 184, 2003, pp.414
! identify ipmx, the discrete velocity direction which maximizes the
! quantity n^(hat) dot e_alpha, where n^(hat) is the out-normal vector
! of the wall at the point (xp2, yp2, zp2).
      prod0 = -100.0
      do ipop = 1,npop-1
        ix = cix(ipop)
        iy = ciy(ipop)
        iz = ciz(ipop)
        prod = real(ix)*xx0 + real(iy)*yy0 + real(iz)*zz0
        if(ipop <= 6) prod = prod*sqrt(2.0)
        if(prod > prod0)then
          ipmx = ipop
          prod0 = prod
        end if
      end do


! Caiazzo A. progress in CFD, vol. 8, 2008
! equilibrium + non-equilibrium refill
! first obtain the non-equilibrium part by copying from the neighbouring
! fluid node along the ipmx direction
      ix = cix(ipmx)
      iy = ciy(ipmx)
      iz = ciz(ipmx)

      ip1 = i + ix
      jp1 = j + iy
      kp1 = k + iz

! periodicity
!     if(ip1 < 1) ip1 = ip1 + lx
!     if(ip1 > lx) ip1 = ip1 - lx
          
          if(ip1 < 1) then
          ibp1 = 2
          ib0p1 = 2
          else if (ip1 > lx) then
          ibp1 = 2
          ib0p1 = 2
          else
          if(jp1 > ly) then
          ibp1 = tmpiU(ip1,kp1)
          ib0p1 = tmpiU0(ip1,kp1)
          else if (jp1 < 1) then
          ibp1 = tmpiD(ip1,kp1)
          ib0p1 = tmpiD0(ip1,kp1)
          else
            if(kp1 > lz) then
            ibp1 = tmpiR(ip1,jp1)
            ib0p1 = tmpiR0(ip1,jp1)
            else if(kp1 < 1 ) then
            ibp1 = tmpiL(ip1,jp1)
            ib0p1 = tmpiL0(ip1,jp1)
            else
            ibp1 = ibnodes(ip1,jp1,kp1)
            ib0p1 = ibnodes0(ip1,jp1,kp1)
            end if
          end if
        end if
!!!!!
        IF(ibp1 < 0 .and. ib0p1 < 0)THEN
        rho9 = 0.0 
        u9 = 0.0
        v9 = 0.0
        w9 = 0.0
      
          if(jp1 > ly) then
          f9 = tmpfU(:,ip1,kp1)
          else if (jp1 < 1) then
          f9 = tmpfD(:,ip1,kp1)
          else
            if(kp1 > lz) then
            f9 = tmpfR(:,ip1,jp1)
            else if(kp1 < 1 ) then
            f9 = tmpfL(:,ip1,jp1)
            else
            f9 = f(:,ip1,jp1,kp1)
            end if
          end if
 
        do ipop = 0,npop-1
          rho9 = rho9 + f9(ipop)  
          u9 = u9 + real(cix(ipop))*f9(ipop)    
          v9 = v9 + real(ciy(ipop))*f9(ipop)   
          w9 = w9 + real(ciz(ipop))*f9(ipop)   
        end do
        call feqpnt(u9,v9,w9,rho9,feq9)

! note: below LHS = f(:,i,j,k), NOT f9(:,i,j,k)
        f(:,i,j,k) = f9 - feq9
      ELSE
        f(:,i,j,k) = 0.0
      END IF

! now calculate the equilibrium part
! first obtain the local mean density
      nghb = 0
      rho9 = 0.0

      do ipop = 1,npop-1
        ix = cix(ipop)
        iy = ciy(ipop)
        iz = ciz(ipop)

        ip1 = i + ix
        jp1 = j + iy
        kp1 = k + iz

! periodicity
      !  if(ip1 < 1) ip1 = ip1 + lx
      !  if(ip1 > lx) ip1 = ip1 - lx

        if(ip1 > lx) then
          ibp1 = 2
          ib0p1 = 2
        else if (ip1 < 1) then
          ibp1 = 2
          ib0p1 = 2
        else
          if(jp1 > ly) then
          ibp1 = tmpiU(ip1,kp1)
          ib0p1 = tmpiU0(ip1,kp1)
          else if (jp1 < 1) then
          ibp1 = tmpiD(ip1,kp1)
          ib0p1 = tmpiD0(ip1,kp1)
          else
            if(kp1 > lz) then
            ibp1 = tmpiR(ip1,jp1)
            ib0p1 = tmpiR0(ip1,jp1)
            else if(kp1 < 1 ) then
            ibp1 = tmpiL(ip1,jp1)
            ib0p1 = tmpiL0(ip1,jp1)
            else
            ibp1 = ibnodes(ip1,jp1,kp1)
            ib0p1 = ibnodes0(ip1,jp1,kp1)
            end if
          end if
        end if

        IF(ibp1 < 0 .and. ib0p1 < 0)THEN
          nghb = nghb + 1


          if(jp1 > ly) then
          f9 = tmpfU(:,ip1,kp1)
          else if (jp1 < 1) then
          f9 = tmpfD(:,ip1,kp1)
          else
            if(kp1 > lz) then
            f9 = tmpfR(:,ip1,jp1)
            else if(kp1 < 1 ) then
            f9 = tmpfL(:,ip1,jp1)
            else
            f9 = f(:,ip1,jp1,kp1)
            end if
          end if

          do ii = 0,npop-1
            rho9 = rho9 + f9(ii)
          end do
        end if
      end do

      if(nghb > 0)then
! use the locally averaged density if available
        rho9 = rho9 / real(nghb)
      ELSE
! otherwise use the global average density rho0=1.0
!       rho9 = rho0
        rho9 = 0.0d0
      END IF


! then calculate the previous solid node's velocity
      xx0 = xpnt - xc
      yy0 = ypnt - yc
      zz0 = zpnt - zc

      w1 = wp(1,id)
      w2 = wp(2,id)
      w3 = wp(3,id)

      omg1 = omgp(1,id)
      omg2 = omgp(2,id)
      omg3 = omgp(3,id)

      u9 = w1 + omg2*zz0 - omg3*yy0
      v9 = w2 + omg3*xx0 - omg1*zz0
      w9 = w3 + omg1*yy0 - omg2*xx0

      call feqpnt(u9,v9,w9,rho9,feq9)

! equilibrium + non-equilibrium refill
      f(:,i,j,k) = f(:,i,j,k) + feq9

!     write(*,*)'istep,i,j,k,nghb,f(:,i,j,k)=',istep,i,j,k,nghb,f(:,i,j,k)
      END IF
      end do
      end do
      end do

      deallocate(tmpfL)
      deallocate(tmpfR)      
      deallocate(tmpfU)
      deallocate(tmpfD)

      deallocate (tmpiL)      
      deallocate (tmpiR)
      deallocate (tmpiU)
      deallocate (tmpiD)

      deallocate (tmpiL0)      
      deallocate (tmpiR0)
      deallocate (tmpiU0)
      deallocate (tmpiD0)   

      end subroutine beads_filling
!===========================================================================

      subroutine exchng5(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8)
      use mpi
      use var_inc
      implicit none

!                     1          2       3        4
!                     5          6       7        8
!      call exchng5(f(:,:,:,1),tmpfR,f(:,:,:,lz),tmpfL,          &
!                   f(:,:,1,:),tmpfU,f(:,:,ly,:),tmpfD)


      integer ileny, ilenz
      real, dimension(0:npop-1,lx,ly)    :: tmp1, tmp2, tmp3, tmp4
      real, dimension(0:npop-1,lx,lz)    :: tmp5, tmp7
      real, dimension(0:npop-1,lx,0:lz+1):: tmp5l, tmp6, tmp7l, tmp8


      integer status_array(MPI_STATUS_SIZE,4), req(4)

      ilenz = npop*lx*ly
      ileny = npop*lx*(lz+2)

      call MPI_IRECV(tmp2,ilenz,MPI_REAL8,mzp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp4,ilenz,MPI_REAL8,mzm,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp1,ilenz,MPI_REAL8,mzm,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp3,ilenz,MPI_REAL8,mzp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      tmp5l(:,:,0) = tmp4(:,:,1)
      tmp5l(:,:,1:lz) = tmp5(:,:,:)
      tmp5l(:,:,lz+1) = tmp2(:,:,1)

      tmp7l(:,:,0) = tmp4(:,:,ly)
      tmp7l(:,:,1:lz) = tmp7(:,:,:)
      tmp7l(:,:,lz+1) = tmp2(:,:,ly)

      call MPI_IRECV(tmp6,ileny,MPI_REAL8,myp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp8,ileny,MPI_REAL8,mym,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp5l,ileny,MPI_REAL8,mym,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp7l,ileny,MPI_REAL8,myp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      end subroutine exchng5
!===========================================================================

      subroutine exchng5i(tmp1i,tmp2i,tmp3i,tmp4i,tmp5i,tmp6i,tmp7i,tmp8i)
      use mpi
      use var_inc
      implicit none

!                         1           2         3          4
!                         5           6         7          8
!      call exchng5i(ibnodes(:,:,1),tmpiR,ibnodes(:,:,lz),tmpiL, &
!                    ibnodes(:,1,:),tmpiU,ibnodes(:,ly,:),tmpiD)


!                         1             2         3             4 
!                         5             6         7             8
!      call exchng5i(ibnodes0(:,:,1),tmpiR0,ibnodes0(:,:,lz),tmpiL0,  &
!                    ibnodes0(:,1,:),tmpiU0,ibnodes0(:,ly,:),tmpiD0)



      integer ilenz,ileny
      integer, dimension(lx,ly):: tmp1i, tmp2i, tmp3i, tmp4i  
      integer, dimension(lx,lz)    :: tmp5i, tmp7i
      integer, dimension(lx,0:lz+1):: tmp5il, tmp6i, tmp7il, tmp8i

      integer status_array(MPI_STATUS_SIZE,4), req(4)

      ilenz = lx*ly
      ileny = lx*(lz+2)

      call MPI_IRECV(tmp2i,ilenz,MPI_INTEGER,mzp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp4i,ilenz,MPI_INTEGER,mzm,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp1i,ilenz,MPI_INTEGER,mzm,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp3i,ilenz,MPI_INTEGER,mzp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      tmp5il(:,0) = tmp4i(:,1)
      tmp5il(:,1:lz) = tmp5i(:,:)
      tmp5il(:,lz+1) = tmp2i(:,1)

      tmp7il(:,0) = tmp4i(:,ly)
      tmp7il(:,1:lz) = tmp7i(:,:)
      tmp7il(:,lz+1) = tmp2i(:,ly)

      call MPI_IRECV(tmp6i,ileny,MPI_INTEGER,myp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp8i,ileny,MPI_INTEGER,mym,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp5il,ileny,MPI_INTEGER,mym,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp7il,ileny,MPI_INTEGER,myp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      end subroutine exchng5i
!===========================================================================

      subroutine feqpnt(u9,v9,w9,rho9,feq9)
      use mpi
      use var_inc
      implicit none

      integer ip
      real u9, v9, w9, rho9
      real feq9(0:npop-1)
      real usqr, G

      usqr = u9*u9 + v9*v9 + w9*w9
      usqr = 1.5*usqr

      feq9(0) = ww0*(rho9 - usqr)

      do ip = 1,6
        G = (cix(ip)*u9 + ciy(ip)*v9 + ciz(ip)*w9)
        feq9(ip) = ww1*(rho9 + 3.0*G + 4.5*G*G - usqr)
      end do
 
      do ip = 7,npop-1
        G = (cix(ip)*u9 + ciy(ip)*v9 + ciz(ip)*w9)
        feq9(ip) = ww2*(rho9 + 3.0*G + 4.5*G*G - usqr)
      end do

      end subroutine feqpnt 
!===========================================================================


