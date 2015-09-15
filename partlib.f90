!=============================================================================
!@subroutine initpartpos
!@desc Randomly spawns solid particles into the flow such that the amount of
!      solid particles per local MPI task does not exceed the set maximum and
!      all solid particles do not overlap.
!=============================================================================
      subroutine initpartpos
      use mpi
      use var_inc
      implicit none

      integer ip,cnt, idf,idfy,idfz, i, j, k
      integer npsflag, npsflagt, npst, ovlpcnt
      real ranpp
      real radpgap

      radpgap = rad + mingap_w

      iseedp = -iseedp !Random seed
      iyp = 0
      ivp = 0
      ovlpflagt = 1 !Solid particle overlap flags

      DO
        if(myid == 0)then
        !For each solid particle that needs a new location generate new global position
         do ip = 1,npart
           if(ovlpflagt(ip) > 0)then
             ypglb(1,ip) = radpgap +(real(nx)-2.0*radpgap)*ranpp(iseedp,iyp,ivp)
             ypglb(2,ip) = real(ny)*ranpp(iseedp,iyp,ivp)
             ypglb(3,ip) = real(nz)*ranpp(iseedp,iyp,ivp)
           end if
         end do
        end if

        !Distribute all global solid particle positions across all MPI tasks
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ypglb,3*npart,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

        nps = 0
        ipglb = 0
        yp = 0.0

        !Determine which solid particle centers reside within the MPI tasks local domain
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

        !Check to see if the number of solid particles in the local domain exceeds set limit
        npsflag = 0
        if(nps > msize)then
          write(*,*) 'fault: nps = ', nps, ' > msize',msize,'nproc,npart',nproc,npart
         npsflag = 1
        end if
        call MPI_ALLREDUCE(npsflag,npsflagt,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        if(npsflagt > 0) stop

        !Now check for particle overlap
        ovlpflag = 0
        ovlpflagt = 0
        call ovlpchk
        call MPI_ALLREDUCE(ovlpflag,ovlpflagt,npart,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        ovlpcnt = count(ovlpflagt /= 0)
        if(myid == 0) write(*,*) 'ovlpcnt = ', ovlpcnt

        !If there's no particle overlap, we are done with initialization
        if(ovlpcnt == 0) then
          if(myid==0)write(*,*) 'Particle position initialization done.'
          call MPI_ALLREDUCE(nps,npst,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
          if(npst .ne. npart)then
            write(*,*) 'fault: npst != npart'
            stop
          end if
          exit !Break out of loop
        end if
        !Other wise respawn particles that do overlap to new random locations
      END DO

      end subroutine initpartpos
!=============================================================================
!@subroutine ranpp
!@desc Minimal random number generator of Park and Miller with Bays-Durham 
!      shuffle and added safegards, see Numerical Recipes.
!@param idum = integer
!@param iy1 = integer
!@param iv1 = integer
!=============================================================================
      Function ranpp(idum,iy1,iv1)
      implicit none

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
!=============================================================================
!@subroutine ovlpchk
!@desc Checks for solid particle overlap
!=============================================================================
      subroutine ovlpchk
      use var_inc
      implicit none

      integer ip, id, ii, iflag
      real xc, yc, zc, xb, yb, zb, dist, hgap

      !Loop through particles with center in local domain
      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)

        if(id < npart)then
          iflag = 0
          !Parse through other particles with a id greater than current
          do ii = id+1,npart
            xb = ypglb(1,ii)
            yb = ypglb(2,ii)
            zb = ypglb(3,ii)

            !Get closest center with periodic boundaries
            !if((xb - xc) > real(nxh)) xb = xb - real(nx)
            !if((xb - xc) < -real(nxh)) xb = xb + real(nx)

            if((yb - yc) > real(nyh)) yb = yb - real(ny)
            if((yb - yc) < -real(nyh)) yb = yb + real(ny)

            if((zb - zc) > real(nzh)) zb = zb - real(nz)
            if((zb - zc) < -real(nzh)) zb = zb + real(nz)

            dist = sqrt((xc - xb)**2 + (yc - yb)**2 + (zc -zb)**2)
            hgap = dist - 2.0*rad

            if(istep00 == 0) hgap = hgap - mingap

            !If we have particle overlap something when wrong
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
!=============================================================================
!@subroutine initpartvel
!@desc Initialize solid particle linear velocity through interpolation of the
!      surrounding fluid
!=============================================================================
      subroutine initpartvel
      use mpi
      use var_inc
      implicit none

      integer ip, id, ilen, idir, iflag   
      real xc, yc, zc, v0      

      real, dimension(3,npart):: wp0
      real,allocatable,dimension(:,:):: tmpYm, tmpYp, tmpZm, tmpZp

      allocate (tmpYm(lx,0:lz+1))   
      allocate (tmpYp(lx,0:lz+1))
      allocate (tmpZm(lx,ly))
      allocate (tmpZp(lx,ly))

      !Set all solid particle velocities to 0
      wp0 = 0.0

      !"iflag = 1", indicates velocity interpolation
      iflag = 1

      !Interpolate solid particle velocity in the x direction
      !Get neighboring task velocity data      
      call initpart_exchng(ux(:,1,:),tmpYm,ux(:,ly,:),tmpYp,ux(:,:,1),tmpZm,ux(:,:,lz),tmpZp)
      !Parse through solid particles
      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)
        idir =  1 
        !Interpolated velocity v0
        call trilinear(tmpYm,tmpYp,tmpZm,tmpZp,xc,yc,zc,idir,iflag,v0,id)
        wp0(1,id) = v0
      end do

      !Interpolate solid particle velocity in the y direction
      !Get neighboring task velocity data  
      call initpart_exchng(uy(:,1,:),tmpYm,uy(:,ly,:),tmpYp,uy(:,:,1),tmpZm,uy(:,:,lz),tmpZp)
      !Parse through solid particles
      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)
        idir = 2
        !Interpolated velocity v0
        call trilinear(tmpYm,tmpYp,tmpZm,tmpZp,xc,yc,zc,idir,iflag,v0,id)
        wp0(2,id) = v0
      end do

      !Interpolate solid particle velocity in the z direction
      !Get neighboring task velocity data  
      call initpart_exchng(uz(:,1,:),tmpYm,uz(:,ly,:),tmpYp,uz(:,:,1),tmpZm,uz(:,:,lz),tmpZp)
      !Parse through solid particles
      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)
        idir = 3  
        !Interpolated velocity v0
        call trilinear(tmpYm,tmpYp,tmpZm,tmpZp,xc,yc,zc,idir,iflag,v0,id)
        wp0(3,id) = v0
      end do

      deallocate (tmpYm)   
      deallocate (tmpYp)   
      deallocate (tmpZp)
      deallocate (tmpZm)
      
      !Collect velocity of all solid particles across all tasks
      ilen = 3*npart
      call MPI_ALLREDUCE(wp0,wp,ilen,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

      end subroutine initpartvel
!=============================================================================
!@subroutine initpartomg
!@desc Initialize solid particle rotational velocity through interpolation of
!      the vorticity of the surrounding fluid
!=============================================================================
      subroutine initpartomg
      use mpi
      use var_inc
      implicit none

      integer ip, id, ilen, idir, iflag  
      real xc, yc, zc, v0  

      real, dimension(3,npart):: omgp0
      real,allocatable,dimension(:,:):: tmpYm, tmpYp, tmpZm, tmpZp

      allocate (tmpYm(lx,0:lz+1))   
      allocate (tmpYp(lx,0:lz+1))
      allocate (tmpZm(lx,ly))
      allocate (tmpZp(lx,ly))
            omgp0 = 0.0

      !Prepare vorticity field ox, oy, and oz
      call vortcalc

      !"iflag = 2", indicates vorticity interpolation
      iflag = 2  

      !Interpolate solid particle vorticity in the x direction
      !Get neighboring task vorticity data  
      call initpart_exchng(ox(:,1,:),tmpYm,ox(:,ly,:),tmpYp,ox(:,:,1),tmpZm,ox(:,:,lz),tmpZp)
      !Parse through solid particles
      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)
        idir = 1
        !Interpolated vorticity v0
        call trilinear(tmpYm,tmpYp,tmpZm,tmpZp,xc,yc,zc,idir,iflag,v0,id)   
        ! vorticity is twice the rotation rate of a small fluid line segment
        ! oriented along a principal direction of rate-of-strain tensor
        omgp0(1,id) = 0.5*v0
      end do

      !Interpolate solid particle vorticity in the y direction
      !Get neighboring task vorticity data  
      call initpart_exchng(oy(:,1,:),tmpYm,oy(:,ly,:),tmpYp,oy(:,:,1),tmpZm,oy(:,:,lz),tmpZp)
      !Parse through solid particles
      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)
        idir = 2
        !Interpolated vorticity v0
        call trilinear(tmpYm,tmpYp,tmpZm,tmpZp,xc,yc,zc,idir,iflag,v0,id)   
        omgp0(2,id) = 0.5*v0
      end do

      !Interpolate solid particle vorticity in the z direction
      !Get neighboring task vorticity data  
      call initpart_exchng(oz(:,1,:),tmpYm,oz(:,ly,:),tmpYp,oz(:,:,1),tmpZm,oz(:,:,lz),tmpZp)
      !Parse through solid particles
      do ip = 1,nps
        xc = yp(1,ip)
        yc = yp(2,ip)
        zc = yp(3,ip)
        id = ipglb(ip)
        idir = 3
        !Interpolated vorticity v0
        call trilinear(tmpYm,tmpYp,tmpZm,tmpZp,xc,yc,zc,idir,iflag,v0,id)
        omgp0(3,id) = 0.5*v0
      end do

      deallocate (tmpYm)   
      deallocate (tmpYp)   
      deallocate (tmpZp)
      deallocate (tmpZm)

      !Collect vorticity of all solid particles across all tasks
      ilen = 3*npart
      call MPI_ALLREDUCE(omgp0,omgp,ilen,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

      omgp = 0.0 !Not sure why this is here

      end subroutine initpartomg
!=============================================================================
!@subroutine initpart_exchng
!@desc MPI exchange to send/receive velocity or vorticity data with neighboring
!      tasks. Used for initializing solid particles
!@param tmpYmSi,tmpYpSi = real array; local data to be sent to neighboring Y tasks
!@param tmpYmR, tmpYpR = real array; stores data recieved from Y tasks
!@param tmpZmS,tmpZpS = real array; local data to be sent to neighboring Z tasks
!@param tmpZmR, tmpZpR = real array; stores data recieved from Z tasks
!=============================================================================
      subroutine initpart_exchng(tmpYmSi,tmpYmR,tmpYpSi,tmpYpR,tmpZmS,tmpZmR,tmpZpS,tmpZpR)
      use mpi
      use var_inc
      implicit none

      integer ilenz, ileny
      real, dimension(lx,ly)  :: tmpZmS, tmpZmR, tmpZpS, tmpZpR
      real, dimension(lx,lz)  :: tmpYmSi, tmpYpSi
      real, dimension(lx,0:lz+1):: tmpYmS, tmpYmR, tmpYpS, tmpYpR  

      integer status_array(MPI_STATUS_SIZE,4), req(4)
     
      ilenz = lx*ly
      ileny = lx*(lz+2)
      !Send data to neighboring MPI tasks in the Z direction
      call MPI_IRECV(tmpZpR,ilenz,MPI_REAL8,mzp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmpZmR,ilenz,MPI_REAL8,mzm,1,MPI_COMM_WORLD,req(2),ierr)

      call MPI_ISEND(tmpZmS,ilenz,MPI_REAL8,mzm,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmpZpS,ilenz,MPI_REAL8,mzp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

      !Combine data cornor data from Z neighbor with Y send buffer
      tmpYmS(:,0) = tmpZmR(:,1)
      tmpYmS(:,1:lz) = tmpYmSi(:,:)
      tmpYmS(:,lz+1) = tmpZpR(:,1)

      tmpYpS(:,0) = tmpZmR(:,ly)
      tmpYpS(:,1:lz) = tmpYpSi(:,:)
      tmpYpS(:,lz+1) = tmpZpR(:,ly)

      !Send data to neighboring MPI tasks in the Y direction
      call MPI_IRECV(tmpYpR,ileny,MPI_REAL8,myp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmpYmR,ileny,MPI_REAL8,mym,1,MPI_COMM_WORLD,req(2),ierr)

      call MPI_ISEND(tmpYmS,ileny,MPI_REAL8,mym,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmpYpS,ileny,MPI_REAL8,myp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

      end subroutine initpart_exchng
!=============================================================================
!@subroutine trilinear
!@desc Uses trilinear interpolation to interpolates intial solidparticle 
!      velocity or rotation. Used by initpartvel and initpartomg
!@param tmpYm,tmpYp,tmpZm,tmpZp = real array; velocity or vorticity
!      data from neighboring MPI tasks
!@param xc,yc,zc = integer; solid particle's center local position
!@param idir = integer, deminsion under consideration. 
!       idir = 1 : x, idir = 2 : y, idir = 3 : z
!@param iflag = integer; flag indicating velocity or vorticity
!       iflag = 1 : velocity, iflag = 2 : vorticity
!@param v0 = real; interpolated value
!=============================================================================
      subroutine trilinear(tmpYm,tmpYp,tmpZm,tmpZp,xc,yc,zc,idir,iflag,v0) 
      use var_inc
      implicit none

      integer idir, iflag, i, j, k, n, it, jt, kt
      real xc, yc, zc, dx, dy, dz, v0
      real tx, ty, tz, tx0,  ty0, tz0

      real, dimension(lx,ly):: tmpZm,tmpZp   
      real, dimension(lx,0:lz+1):: tmpYm, tmpYp 
      integer, dimension(8):: trix, triy, triz
      real, dimension(8):: v

      !Hold indices used to get interpolation points
      trix = (/0,1,1,0,0,1,1,0/)
      triy = (/0,0,1,1,0,0,1,1/)
      triz = (/0,0,0,0,1,1,1,1/)

      dx = 1.0
      dy = 1.0
      dz = 1.0

      !Get the corner position of the box used to interpolate
      i = int((xc + 0.5) / dx)
      j = int((yc - real(indy*ly) + 0.5) / dy)
      k = int((zc - real(indz*lz) + 0.5) / dz)

      !Calculate constants used for interpolation
      tx = (xc + 0.5 - real(i)) / dx
      ty = (yc - real(indy*ly) + 0.5 - real(j)) / dy
      tz = (zc - real(indz*lz) + 0.5 - real(k)) / dz

      tx0 = 1.0 - tx
      ty0 = 1.0 - ty
      tz0 = 1.0 - tz

      !Retreive the 8 interpolation values needed
      do n = 1, 8
        !Get location of interpolation point
        it = i + trix(n)
        jt = j + triy(n)
        kt = k + triz(n)

        !If node point is outside our local domain
        if(jt < 1)then
          v(n) = tmpYm(it,kt)
        elseif(jt > ly)then
          v(n) = tmpYp(it,kt)
        elseif(kt < 1)then
          v(n) = tmpZm(it,jt)
        elseif(kt > lz)then
          v(n) = tmpZp(it,jt)
        !If node point is in our local domain
        else
          !Interpolate velocity
          if(iflag == 1)then
            if(idir == 1)then
             v(n) = ux(it,jt,kt)
            elseif(idir == 2)then
             v(n) = uy(it,jt,kt)
            elseif(idir == 3)then
             v(n) = uz(it,jt,kt)
            else
              write(*,'(A39,I3)')'Trilinear (partlib.f90): Invalid idir: ', idir
              write(*,'(A43)')'idir = 1 : x, idir = 2 : y, idir = 3 : z'
            endif
          !Interpolate vorticity
          elseif(iflag == 2)then
            if(idir == 1)then
             v(n) = ox(it,jt,kt)
            elseif(idir == 2)then
             v(n) = oy(it,jt,kt)
            elseif(idir == 3)then
             v(n) = oz(it,jt,kt)
            else
              write(*,'(A39,I3)')'Trilinear (partlib.f90): Invalid idir: ', idir
              write(*,'(A43)')'idir = 1 : x, idir = 2 : y, idir = 3 : z'
            endif
          else
            write(*,'(A40,I3)')'Trilinear (partlib.f90): Invalid iflag: ', iflag
            write(*,'(A43)')'iflag = 1 : velocity, iflag = 2 : vorticity'
          endif
        endif
      enddo

      !Interpolate
      v0 = tx0*ty0*tz0*v(1) + tx*ty0*tz0*v(2) + tx*ty*tz0*v(3)  &
         + tx0*ty*tz0*v(4) + tx0*ty0*tz*v(5) +tx*ty0*tz*v(6)    &
         + tx*ty*tz*v(7) + tx0*ty*tz*v(8)

      end subroutine trilinear
!=============================================================================
!@subroutine beads_links
!@desc beads_links does a variety of functions. The main purpose of beads_links
!      is to determine which velocity links cross a solid particle boundary.
!      This subroutine also updates which nodes are inside a solid particle,
!      which nodes will require filling, which MPI neighbors we require data for
!      filling, and which distribution data we need from neighbors for inter-
!      polation bounceback.
!=============================================================================
      subroutine beads_links
      use var_inc
      implicit none

      integer ip,ippp, i, j, k, id, imove, jmove, kmove
      real xc, yc, zc, aa, bb, cc    
      real alpha0, alpha1, alpha     
      real xx0, yy0, zz0, rr0, rr01      
      real xx1, yy1, zz1, rr1, rr11,xxtt   
      real xpnt, ypnt, zpnt, radp2

      integer nfbeads
      logical bfilled
      real alphay, alphaz
      real,dimension(2*msize):: idfbeads, xfbeads, yfbeads, zfbeads

      !Reset indexers and flag array for interpolation bounce back requests
      ipf_mymc = 0; ipf_mypc = 0; ipf_mzmc = 0; ipf_mzpc = 0
      iblinks(:,:,:) = 0

      !Reset fill request data
      localReqData = .FALSE.
      radp2 = rad + 2.d0
      nfbeads = 0
      nlink = 0
      nbfill = 0
      alphay = dfloat(indy*ly)+ 0.5d0
      alphaz = dfloat(indz*lz)+ 0.5d0
      
      !Execute a global search for particles that are touching our local domain (course search)
      do id=1,npart
        yc = ypglb(2,id)
        zc = ypglb(3,id)

        !Use the nearest particle center instead of the real center
        if((yc - alphay) > dfloat(nyh)) yc = yc - dfloat(ny)
        if((yc - alphay) < -dfloat(nyh)) yc = yc + dfloat(ny)
        if((zc - alphaz) > dfloat(nzh)) zc = zc - dfloat(nz)
        if((zc - alphaz) < -dfloat(nzh)) zc = zc + dfloat(nz)

        !If a solid particle is touching our local domain save its data
        if((alphay+ly+radp2)-yc > 0 .AND. (alphay-radp2)-yc < 0 &
            .AND. (alphaz+lz+radp2)-zc > 0 .AND. (alphaz-radp2)-zc < 0)then
          nfbeads = nfbeads + 1
          idfbeads(nfbeads) = id
          xfbeads(nfbeads) = ypglb(1,id)
          yfbeads(nfbeads) = ypglb(2,id)
          zfbeads(nfbeads) = ypglb(3,id)
        endif
      enddo

      !Interate across our local fluid domain (fine search)
      do i=1, lx
        do j=1, ly
          do k=1, lz
            xpnt = dfloat(i) - 0.5d0
            ypnt = dfloat(j) - 1.d0 + alphay
            zpnt = dfloat(k) - 1.d0 + alphaz
            bfilled = .FALSE.

            !Loop through particles saved from the course search
            do 111 id=1, nfbeads
              !Check to see if this node is close to the solid particle
              xc = xfbeads(id)
              xx0 = xpnt - xc
              if(xx0**2 > radp2**2)GO TO 111

              yc = yfbeads(id)
              if((yc - ypnt) > dfloat(nyh)) yc = yc - dfloat(ny)
              if((yc - ypnt) < -dfloat(nyh)) yc = yc + dfloat(ny)
              yy0 = ypnt - yc
              if(xx0**2+yy0**2 > radp2**2)GO TO 111

              zc = zfbeads(id)
              if((zc - zpnt) > dfloat(nzh)) zc = zc - dfloat(nz)
              if((zc - zpnt) < -dfloat(nzh)) zc = zc + dfloat(nz)
              zz0 = zpnt - zc
              rr0 = xx0**2 + yy0**2 + zz0**2

              if(sqrt(rr0) <= radp2)then
                rr01 = rr0 - (rad*1.d0)**2
                !If we are inside the solid particle update ibnodes and isnodes
                if(rr01 <= 0.d0)then
                  ibnodes(i,j,k) = 1 
                  isnodes(i,j,k) = idfbeads(id)
                  bfilled = .TRUE.
                !If not check links for solid particle boundary
                else
                  !Determine which node links are crossing the solid boundary
                  do ip = 1,npop-1
                    imove = i + cix(ip) 
                    jmove = j + ciy(ip)
                    kmove = k + ciz(ip)

                    xx1 = dfloat(imove) - 0.5d0 - xc
                    yy1 = dfloat(jmove) - 0.5d0 - yc + dfloat(indy*ly)
                    zz1 = dfloat(kmove) - 0.5d0 - zc + dfloat(indz*lz)

                    rr1 = xx1**2 + yy1**2 + zz1**2  
                    rr11 = rr1 - (rad*1.d0)**2  

                    !If the neighboring node in the given direction is inside the solid particle
                    !we know that this link crosses the solid boundary
                    if(rr11 <= 0.d0)then
                      nlink = nlink + 1
                      if(nlink.ge.maxlink)then
                        write(*,'(A44,I5)') 'Number of links exceeded array size, maxlink: ',maxlink
                        stop
                      endif

                      !Save link information
                      xlink(nlink) = i
                      ylink(nlink) = j
                      zlink(nlink) = k

                      iplink(nlink) = ip
                      mlink(nlink) = idfbeads(id)

                      !Calculate the position of the solid boundary on the link
                      aa = rr0 + rr1 - 2.d0*(xx0*xx1 + yy0*yy1 + zz0*zz1)
                      bb = xx1*(xx0-xx1) + yy1*(yy0-yy1) + zz1*(zz0-zz1)    
                      cc = rr11  
                      alpha0 = bb/aa 
                      alpha1 = sqrt(alpha0**2 - cc/aa)    
                      alpha = -alpha0 + alpha1

                      alpha = 1.d0 - alpha      
                      alink(nlink) = alpha

                      if(alpha < 0.d0 .or. alpha > 1.d0)then 
                        write(*,*) 'fault: alpha = ', alpha, ' @ ilink(',      &
                                   ip, ', ', i, ', ', j, ', ', k, ')' 
                      end if
                      !Check to see if we will need data from neighboring MPI tasks for interpolation
                      call parse_MPI_links(ip, i, j, k, alpha, nlink)
                      
                    end if
                  end do
                end if
              endif

111          continue !npart
              !If a node is not inside a bead and ibnodes is > 1, we know it needs to be filled
              if(.NOT. bfilled .AND. ibnodes(i,j,k) > 0)then
                nbfill = nbfill + 1
                if(nbfill.ge.maxbfill)then
                  write(*,'(A44,I5,I4)') 'Number of fill node exceeded array size, maxbfill: ',maxbfill, myid
                  stop
                endif

                !Safe filling data
                xbfill(nbfill) = i
                ybfill(nbfill) = j
                zbfill(nbfill) = k
                idbfill(nbfill) = isnodes(i,j,k)
                !Correct ibnodes and isnodes
                ibnodes(i,j,k) = -1
                isnodes(i,j,k) = -1

                !Determine if we need to request any data for filling
                if(k == 1)then
                  localReqData(2) = .TRUE.
                  if(j == 1) localReqData(7) = .TRUE.
                  if(j == ly) localReqData(8) = .TRUE.
                else if(k == lz)then
                  localReqData(1) = .TRUE.
                  if(j == 1) localReqData(5) = .TRUE.
                  if(j == ly) localReqData(6) = .TRUE.
                else
                endif
                if(j == 1) localReqData(3) = .TRUE.
                if(j == ly) localReqData(4) = .TRUE.
              endif !If needs filling
          enddo !z
        enddo !y
      enddo !x

      !Update ghost in ibnodes nodes
      !This is needed in beads filling and post streaming interpolation
      ibnodes(:,0:ly+1,0)=-1    
      ibnodes(:,0:ly+1,lz+1)=-1
      do k=0,lz+1,lz+1
        do i = 1, lx
          do j = 0, ly+1
            xpnt = dfloat(i) - 0.5d0
            ypnt = dfloat(j) - 1.d0 + alphay
            zpnt = dfloat(k) - 1.d0 + alphaz

            do 114 id=1, nfbeads
              xc = xfbeads(id)
              xx0 = xpnt - xc
              if(xx0**2 > radp2**2)GO TO 114

              yc = yfbeads(id)
              if((yc - ypnt) > dfloat(nyh)) yc = yc - dfloat(ny)
              if((yc - ypnt) < -dfloat(nyh)) yc = yc + dfloat(ny)
              yy0 = ypnt - yc
              if(xx0**2+yy0**2 > radp2**2)GO TO 114

              zc = zfbeads(id)
              if((zc - zpnt) > dfloat(nzh)) zc = zc - dfloat(nz)
              if((zc - zpnt) < -dfloat(nzh)) zc = zc + dfloat(nz)
              zz0 = zpnt - zc
              rr0 = xx0**2 + yy0**2 + zz0**2
              
              if(sqrt(rr0) <= radp2)then
                rr01 = rr0 - (rad*1.d0)**2
                if(rr01 <= 0.d0)ibnodes(i,j,k) = 1 
              end if                     
114          continue !npart
          enddo !y
        enddo !x
      enddo !z

      ibnodes(:,0,1:lz)=-1    
      ibnodes(:,ly+1,1:lz)=-1   
      do j=0,ly+1,ly+1
        do i = 1, lx
          do k = 1, lz
            xpnt = dfloat(i) - 0.5d0
            ypnt = dfloat(j) - 1.d0 + alphay
            zpnt = dfloat(k) - 1.d0 + alphaz

            do 115 id=1, nfbeads
              xc = xfbeads(id)
              xx0 = xpnt - xc
              if(xx0**2 > radp2**2)GO TO 115

              yc = yfbeads(id)
              if((yc - ypnt) > dfloat(nyh)) yc = yc - dfloat(ny)
              if((yc - ypnt) < -dfloat(nyh)) yc = yc + dfloat(ny)
              yy0 = ypnt - yc
              if(xx0**2+yy0**2 > radp2**2)GO TO 115

              zc = zfbeads(id)
              if((zc - zpnt) > dfloat(nzh)) zc = zc - dfloat(nz)
              if((zc - zpnt) < -dfloat(nzh)) zc = zc + dfloat(nz)
              zz0 = zpnt - zc
              rr0 = xx0**2 + yy0**2 + zz0**2
              
              if(sqrt(rr0) <= radp2)then
                rr01 = rr0 - (rad*1.d0)**2
                if(rr01 <= 0.d0)ibnodes(i,j,k) = 1 
              end if                     
115         continue !npart
          enddo !y
        enddo !x
      enddo !z

      end subroutine beads_links
!=============================================================================
!@subroutine parse_MPI_links
!@desc Checks to see if we require any data from neighboring MPI tasks for 
!      interpolation bounce back for a given node link that crosses a 
!      solid boundary. If so add data to appropriate request arrays.
!@param ipi = integer; direction of the link under consideration
!@param i,j,k = integer; position node with the link under consideration
!@param alpha = real; distance solid boundary is from the node 0>alpha>1
!@param n = integer; unique index of given link
!=============================================================================
      subroutine parse_MPI_links(ipi,i,j,k,alpha,n)
      use var_inc

      integer ipi, i, j, k, ip, ip2, ix, iy, ik, n
      integer im1, jm1, km1, im2, jm2, km2, ip1, jp1, kp1
      real alpha
      logical im2f

      im2f = .TRUE.

      ix = cix(ipi)
      iy = ciy(ipi)
      iz = ciz(ipi)

      ip1 = i + ix
      jp1 = j + iy
      kp1 = k + iz

      im1 = i - ix
      jm1 = j - iy
      km1 = k - iz

      im2 = i - 2*ix
      jm2 = j - 2*iy
      km2 = k - 2*iz
      
      !Determine which velocity direction is needed
      if(alpha > 0.5)then
        ip = ipopp(ipi)
        ip2 = ipopp(ipi)
        !Mid-link bounce back handling on x boundaries
        if(im1 < 1 .or. im1 > lx)then
          goto 116
        elseif(im2 < 1 .or. im2 > lx)then
          !Adjust for bounce back representation
          ip2 = ipi
          im2 = i - ix
          jm2 = j - iy
          km2 = k - iz
        endif
      else
        ip = ipi
        ip2 = ipi
        !Mid-link bounce back handling on x boundaries
        if(im1 < 1 .or. im1 > lx)then
          goto 116
        elseif(im2 < 1 .or. im2 > lx)then
          im2f = .FALSE. !No need to this if its out of bounds
        endif
      endif

      !Determine if we need to request data from neighboring tasks
      !If the distribution needed is outside our local domain 4 things occur
      !1. Indexers for storing request data are incremented (ipf_mymc, ipf_mypc, ipf_mzmc, ipf_mzpc)
      !2. Request data of the distribution we need is stored with custom data type ipf_node (interpolation fluid node)
      !3. A flag used to determine which array we get interpolation data from is stored in iblinks
      !   (1 = mym, 2 = myp, 3 = mzm, 4 = mzp)
      !4. The request array index is stored in iblinks for easy retrival of this data when needed
      if(jp1 < 1)then
        ipf_mymc = ipf_mymc + 1
        ipf_mym(ipf_mymc) = ipf_node(ipi, ip1, jp1, kp1)
        iblinks(0,1,n) = 1
        iblinks(0,2,n) = ipf_mymc
      elseif(jp1 > ly)then
        ipf_mypc = ipf_mypc + 1
        ipf_myp(ipf_mypc) = ipf_node(ipi, ip1, jp1, kp1)
        iblinks(0,1,n) = 2
        iblinks(0,2,n) = ipf_mypc
      elseif(kp1 < 1)then
        ipf_mzmc = ipf_mzmc + 1
        ipf_mzm(ipf_mzmc) = ipf_node(ipi, ip1, jp1, kp1)
        iblinks(0,1,n) = 3
        iblinks(0,2,n) = ipf_mzmc
      elseif(kp1 > lz)then
        ipf_mzpc = ipf_mzpc + 1
        ipf_mzp(ipf_mzpc) = ipf_node(ipi, ip1, jp1, kp1)
        iblinks(0,1,n) = 4
        iblinks(0,2,n) = ipf_mzpc
      endif

      if(jm1 < 1)then
        ipf_mymc = ipf_mymc + 1
        ipf_mym(ipf_mymc) = ipf_node(ip, im1, jm1, km1)
        iblinks(1,1,n) = 1
        iblinks(1,2,n) = ipf_mymc
        if(alpha > 0.5 .and. jm2 < 1 .and. im2f)then
          ipf_mymc = ipf_mymc + 1
          ipf_mym(ipf_mymc) = ipf_node(ip2, im2, jm2, km2)
          iblinks(2,1,n) = 1
          iblinks(2,2,n) = ipf_mymc
        endif
      elseif(jm1 > ly)then
        ipf_mypc = ipf_mypc + 1
        ipf_myp(ipf_mypc) = ipf_node(ip, im1, jm1, km1)
        iblinks(1,1,n) = 2
        iblinks(1,2,n) = ipf_mypc
        if(alpha > 0.5 .and. jm2 > ly .and. im2f)then
          ipf_mypc = ipf_mypc + 1
          ipf_myp(ipf_mypc) = ipf_node(ip2, im2, jm2, km2)
          iblinks(2,1,n) = 2
          iblinks(2,2,n) = ipf_mypc
        endif
      elseif(km1 < 1)then
        ipf_mzmc = ipf_mzmc + 1
        ipf_mzm(ipf_mzmc) = ipf_node(ip, im1, jm1, km1)
        iblinks(1,1,n) = 3
        iblinks(1,2,n) = ipf_mzmc
        if(alpha > 0.5 .and. im2f)then
          if(jm2 < 1)then
            ipf_mymc = ipf_mymc + 1
            ipf_mym(ipf_mymc) = ipf_node(ip2, im2, jm2, km2)
            iblinks(2,1,n) = 1
            iblinks(2,2,n) = ipf_mymc
          elseif(jm2 > ly)then
            ipf_mypc = ipf_mypc + 1
            ipf_myp(ipf_mypc) = ipf_node(ip2, im2, jm2, km2)
            iblinks(2,1,n) = 2
            iblinks(2,2,n) = ipf_mypc
          else
            ipf_mzmc = ipf_mzmc + 1
            ipf_mzm(ipf_mzmc) = ipf_node(ip2, im2, jm2, km2)
            iblinks(2,1,n) = 3
            iblinks(2,2,n) = ipf_mzmc
          endif
        endif
      elseif(km1 > lz)then
        ipf_mzpc = ipf_mzpc + 1
        ipf_mzp(ipf_mzpc) = ipf_node(ip, im1, jm1, km1)
        iblinks(1,1,n) = 4
        iblinks(1,2,n) = ipf_mzpc
        if(alpha > 0.5 .and. im2f)then
          if(jm2 < 1)then
            ipf_mymc = ipf_mymc + 1
            ipf_mym(ipf_mymc) = ipf_node(ip2, im2, jm2, km2)
            iblinks(2,1,n) = 1
            iblinks(2,2,n) = ipf_mymc
          elseif(jm2 > ly)then
            ipf_mypc = ipf_mypc + 1
            ipf_myp(ipf_mypc) = ipf_node(ip2, im2, jm2, km2)
            iblinks(2,1,n) = 2
            iblinks(2,2,n) = ipf_mypc
          else
            ipf_mzpc = ipf_mzpc + 1
            ipf_mzp(ipf_mzpc) = ipf_node(ip2, im2, jm2, km2)
            iblinks(2,1,n) = 4
            iblinks(2,2,n) = ipf_mzpc
          endif
        endif
      else
        if(alpha > 0.5 .and. im2f)then
          if(jm2 < 1)then
            ipf_mymc = ipf_mymc + 1
            ipf_mym(ipf_mymc) = ipf_node(ip2, im2, jm2, km2)
            iblinks(2,1,n) = 1
            iblinks(2,2,n) = ipf_mymc
          elseif(jm2 > ly)then
            ipf_mypc = ipf_mypc + 1
            ipf_myp(ipf_mypc) = ipf_node(ip2, im2, jm2, km2)
            iblinks(2,1,n) = 2
            iblinks(2,2,n) = ipf_mypc
          elseif(km2 < 1)then
            ipf_mzmc = ipf_mzmc + 1
            ipf_mzm(ipf_mzmc) = ipf_node(ip2, im2, jm2, km2)
            iblinks(2,1,n) = 3
            iblinks(2,2,n) = ipf_mzmc
          elseif(km2 > lz)then
            ipf_mzpc = ipf_mzpc + 1
            ipf_mzp(ipf_mzpc) = ipf_node(ip2, im2, jm2, km2)
            iblinks(2,1,n) = 4
            iblinks(2,2,n) = ipf_mzpc
           endif
        endif
      endif
116   continue

      end subroutine parse_MPI_links
!=============================================================================
!@subroutine beads_collision
!@desc Executes interpolated bounce back to represent the fluid bouncing off 
!      a solid particle boundary. Reference Lallemand and Luo (2003) and 
!      Bouzidia et al (2001) for information regard the interpolation scheme
!      used. Note that this is constructed to function post streaming.
!=============================================================================
      subroutine beads_collision
      use mpi 
      use var_inc
      implicit none

      integer id, i, j, k, ii, jj, kk,  ip, ipp, ix, iy, iz, ilen
      integer im1, im2, ip1, jm1, jm2, jp1, km1, km2, kp1, n
      integer ibm1, ibm2
      real alpha, xc, yc, zc, xx0, yy0, zz0
      real uwx, uwy, uwz, uwpro, ff1, ff2, ff3
      real w1, w2, w3, omg1, omg2, omg3 
      real c1, c2, c3, dff, dxmom, dymom, dzmom
      real xpnt, ypnt, zpnt

      real,dimension(lx,ly,lz):: f9print,alphaprint,ff1print
      real, dimension(3,npart):: fHIp0, torqp0
      real mymIpfRecv(ipf_mymc), mypIpfRecv(ipf_mypc), mzmIpfRecv(ipf_mzmc), mzpIpfRecv(ipf_mzpc)

      fHIp0 = 0.0
      torqp0 = 0.0

      !First get data needed from neighboring MPI tasks
      call exchng2direct(mymIpfRecv,mypIpfRecv,mzmIpfRecv,mzpIpfRecv)

      !Parse through all boundary links found in beads_links
      do n = 1,nlink

        i = xlink(n)
        j = ylink(n)
        k = zlink(n)
        ip = iplink(n)  

        id = mlink(n)
        alpha = alink(n)

        ipp = ipopp(ip)       

        xc = ypglb(1,id) 
        yc = ypglb(2,id) 
        zc = ypglb(3,id) 

        xpnt = real(i) - 0.5
        ypnt = real(j) - 0.5 + real(indy*ly)
        zpnt = real(k) - 0.5 + real(indz*lz)

        !Use the nearest particle center instead of the real center
        !if((xc - xpnt) > real(nxh)) xc = xc - real(nx) !Removed periodicity in x
        !if((xc - xpnt) < -real(nxh)) xc = xc + real(nx)

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

        ix = cix(ip)
        iy = ciy(ip)
        iz = ciz(ip)

        im1 = i - ix
        jm1 = j - iy
        km1 = k - iz
        im2 = i - 2*ix
        jm2 = j - 2*iy
        km2 = k - 2*iz
        ip1 = i + ix
        jp1 = j + iy
        kp1 = k + iz

        xx0 = xpnt + real(ix)*alpha - xc 
        yy0 = ypnt + real(iy)*alpha - yc 
        zz0 = zpnt + real(iz)*alpha - zc 

        uwx = w1 + omg2*zz0 - omg3*yy0
        uwy = w2 + omg3*xx0 - omg1*zz0
        uwz = w3 + omg1*yy0 - omg2*xx0

        uwpro = uwx*real(ix) + uwy*real(iy) + uwz*real(iz)

        !Get interpolation fluid node 1 inside the solid particle (fb)
        !Use iblinks(0,1,:) to determine which array to get data from. (1 = mym, 2 = myp, 3 = mzm, 4 = mzp)
        !Use iblinks(0,2,n) to get index of the distribution needed
        if(iblinks(0,1,n) == 0)then
          ff1 = f(ip,ip1,jp1,kp1,s)
        elseif(iblinks(0,1,n) == 1)then
          ff1 = mymIpfRecv(iblinks(0,2,n))
        elseif(iblinks(0,1,n) == 2)then
          ff1 = mypIpfRecv(iblinks(0,2,n))
        elseif(iblinks(0,1,n) == 3)then
          ff1 = mzmIpfRecv(iblinks(0,2,n))
        else
          ff1 = mzpIpfRecv(iblinks(0,2,n))
        endif

      !If the solid boundary's distance from the fluid node is > 0.5
       IF(alpha > 0.5)then
        !Mid-link bounce back handling
        if(im1 < 1 .or. im1 > lx)then
          ff2 = f(ip,i,j,k,s)
          goto 112
        endif
        !Get interpolation fluid node 2 (ff2)
        if(iblinks(1,1,n) == 0)then
          ff2 = f(ipp,im1,jm1,km1,s)
        elseif(iblinks(1,1,n) == 1)then
          ff2 = mymIpfRecv(iblinks(1,2,n))
        elseif(iblinks(1,1,n) == 2)then
          ff2 = mypIpfRecv(iblinks(1,2,n))
        elseif(iblinks(1,1,n) == 3)then
          ff2 = mzmIpfRecv(iblinks(1,2,n))
        else
          ff2 = mzpIpfRecv(iblinks(1,2,n))
        endif
        !Node is out of bonds so set to IBNODES_TRUE
        if(im2 < 0 .or. im2 > lx+1)then
          ff3 = IBNODES_TRUE
          goto 112
        endif

        !Get interpolation fluid node 3 (ff3)
        if(iblinks(2,1,n) == 0)then
          if(im2 == 0 .or. im2 == lx+1)then !If we're only one out use bounce back!
            ff3 = f(ip,im1,jm1,km1,s)
          else
            !Check Pre-Stream location for solid node conflicts
            if(ibnodes(im1,jm1,km1)>0)then
              ff3 = IBNODES_TRUE
            else
              ff3 = f(ipp,im2,jm2,km2,s)
            endif
          endif
        elseif(iblinks(2,1,n) == 1)then
          ff3 = mymIpfRecv(iblinks(2,2,n))
        elseif(iblinks(2,1,n) == 2)then
          ff3 = mypIpfRecv(iblinks(2,2,n))
        elseif(iblinks(2,1,n) == 3)then
          ff3 = mzmIpfRecv(iblinks(2,2,n))
        else
          ff3 = mzpIpfRecv(iblinks(2,2,n))
        endif

112     continue
        !If solid node conflict with ff3, use 2-point interpolation
        if(ff3 > IBNODES_TRUE - 1)then
            c1 = 0.5 / alpha
            c2 = 1.0 - c1
            f(ipp,i,j,k,s) = c1*ff1 + c2*ff2 - 6.0*wwp(ip)*c1*uwpro
        else !Use 3-point interpolation scheme
            c1 = 1.0 / alpha / (2.0*alpha + 1.0)
            c2 = (2.0*alpha - 1.0) / alpha
            c3 = 1.0 - c1 - c2
            f(ipp,i,j,k,s) = c1*ff1 + c2*ff2 + c3*ff3 - 6.0*wwp(ip)*c1*uwpro
        endif

       !If the solid boundary's distance from the fluid node is <= 0.5
       ELSE
        !Check Pre-Stream location for solid node conflicts
        if(ibnodes(im1,jm1,km1)>0)then
            ff2 = IBNODES_TRUE
            goto 113
        else
            ff2 = f(ip,i,j,k,s)
        endif
        !ff3 is out of bounds so set to IBNODES_TRUE
        if(im1 < 1 .or. im1 > lx)then
          ff3 = IBNODES_TRUE
          goto 113
        endif

        !Get interpolation fluid node 3 (ff3)
        if(iblinks(1,1,n) == 0)then
          !Check Pre-Stream location for solid node conflicts
          if(ibnodes(im2,jm2,km2)>0)then
            ff3 = IBNODES_TRUE
          else
            ff3 = f(ip,im1,jm1,km1,s)
          endif
        elseif(iblinks(1,1,n) == 1)then
          ff3 = mymIpfRecv(iblinks(1,2,n))
        elseif(iblinks(1,1,n) == 2)then
          ff3 = mypIpfRecv(iblinks(1,2,n))
        elseif(iblinks(1,1,n) == 3)then
          ff3 = mzmIpfRecv(iblinks(1,2,n))
        else
          ff3 = mzpIpfRecv(iblinks(1,2,n))
        endif
113     continue

        !If solid node conflict with ff2, use simple bounce back with momentum term
        if(ff2 > IBNODES_TRUE - 1)then
          f(ipp,i,j,k,s) = ff1 - 6.0*wwp(ip)*uwpro
        else if(ff3 > IBNODES_TRUE - 1)then !If solid node conflict with ff3, use 2-point interpolation
          f(ipp,i,j,k,s) = 2.0*alpha*(ff1 - ff2) + ff2 - 6.0*wwp(ip)*uwpro
        else !Use 3-point interpolation scheme
          c1 = alpha*(1.0 + 2.0*alpha)
          c2 = 1.0 - 4.0*alpha*alpha
          c3 = -alpha*(1.0 - 2.0*alpha)
          f(ipp,i,j,k,s) = c1*ff1 + c2*ff2 + c3*ff3 - 6.0*wwp(ip)*uwpro
        end if
       ENDIF

       !Just a back up check to ensure something didn't go wrong, can get removed...
        if(abs(f(ipp,i,j,k,s)) > 1e-2)then
          write(*,*)'=',myid,f(ipp,i,j,k,s),alpha,'='
          write(*,*)ip,i,j,k,iblinks(0,1,n),iblinks(1,1,n),iblinks(2,1,n)
          write(*,*),iblinks(0,2,n),iblinks(1,2,n),iblinks(2,2,n)
          write(*,*)im1,jm1,km1,im2,jm2,km2
          write(*,*)ff1,ff2,ff3,n
        endif

        !Compute force and torque acting on particles
        dff = ff1 + f(ipp,i,j,k,s)
        dxmom = dff*real(ix)
        dymom = dff*real(iy)
        dzmom = dff*real(iz)

        fHIp0(1,id) = fHIp0(1,id) + dxmom
        fHIp0(2,id) = fHIp0(2,id) + dymom
        fHIp0(3,id) = fHIp0(3,id) + dzmom

        torqp0(1,id) = torqp0(1,id) + dzmom*yy0 - dymom*zz0
        torqp0(2,id) = torqp0(2,id) + dxmom*zz0 - dzmom*xx0 
        torqp0(3,id) = torqp0(3,id) + dymom*xx0 - dxmom*yy0
 
      end do     

      !Sum up forcing and torque on all solid particles across all MPI tasks
      ilen = 3*npart
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(fHIp0, fHIp, ilen, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)      
      call MPI_ALLREDUCE(torqp0, torqp, ilen, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      end subroutine beads_collision
!=============================================================================
!@subroutine exchng2direct
!@desc Handles the exchange of distribution data between MPI tasks needed for
!      interpolation bounce back with the solid particles. Consult Nicholas
!      Geneva for further algorithm details.
!@param mymIpfRecv, mypIpfRecv, mzmIpfRecv, mzpIpfRecv = real array; stores
!       distributions requested from neighboring MPI tasks.
!=============================================================================
      subroutine exchng2direct(mymIpfRecv, mypIpfRecv, mzmIpfRecv, mzpIpfRecv)
      use mpi
      use var_inc
      implicit none

      type cornerNode
        integer index, mpidir
      endtype

      real mymIpfRecv(ipf_mymc), mypIpfRecv(ipf_mypc), mzmIpfRecv(ipf_mzmc), mzpIpfRecv(ipf_mzpc)
      real,allocatable,dimension(:):: mymIpfSend, mypIpfSend, mzmIpfSend, mzpIpfSend, ipfRecv
      integer ipf_mzmtc, ipf_mzptc, count, ymcount, ypcount, zmcount, zpcount
      integer status(MPI_STATUS_SIZE), status_array(MPI_STATUS_SIZE,10), req(10)
      integer ip, i, j, dir, xm1, ym1, zm1

      type(cornerNode), dimension(2*19*lx):: mzmt_index, mzpt_index
      type(ipf_node),allocatable,dimension(:):: ipfReq
      type(ipf_node), dimension(2*19*lx):: ipf_mzmt, ipf_mzpt

      ipf_mzmtc = 0; ipf_mzptc = 0;

      !Send request arrays to Y direction neighbors
      call MPI_ISEND(ipf_myp(1:ipf_mypc), ipf_mypc, MPI_IPF_NODE, myp, 99, MPI_COMM_WORLD, req(1), ierr)
      call MPI_ISEND(ipf_mym(1:ipf_mymc), ipf_mymc, MPI_IPF_NODE, mym, 99, MPI_COMM_WORLD, req(2), ierr)
      !Recieve requests from Y neighbors
      do i = 1, 2
        !Note, this process allows us to recieve a dynamic amount of data!
        call MPI_PROBE(MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, status, ierr)
        call MPI_GET_COUNT(status, MPI_IPF_NODE, count, ierr)
        allocate(ipfReq(count))
        call MPI_RECV(ipfReq, count, MPI_IPF_NODE, status(MPI_SOURCE), status(MPI_TAG), MPI_COMM_WORLD, status, ierr)
        !Allocate send buffers, and adjust cordinates to local grid
        if(status(MPI_SOURCE) == myp)then
          allocate(mypIpfSend(count))
          mypIpfSend = 0
          ypcount = count
          ipfReq(:)%y = ipfReq(:)%y + ly
          dir = 1
        else
          allocate(mymIpfSend(count))
          mymIpfSend = 0
          ymcount = count
          ipfReq(:)%y = ipfReq(:)%y - ly
          dir = -1
        endif
        ! Parse Recieved data
        do j = 1, count
          !If requested data is in the z neighbors (corner), add it to temporay array
          if(ipfReq(j)%z > lz)then 
            ipf_mzptc = ipf_mzptc + 1
            ipf_mzpt(ipf_mzptc) = ipfReq(j)
            mzpt_index(ipf_mzptc) = cornerNode(j,dir)
          elseif(ipfReq(j)%z < 1)then
            ipf_mzmtc = ipf_mzmtc + 1
            ipf_mzmt(ipf_mzmtc) = ipfReq(j)
            mzmt_index(ipf_mzmtc) = cornerNode(j,dir)
          else !If distribution is in local domain
            ip =  ipfReq(j)%ip
            !Get pre stream location
            xm1 = ipfReq(j)%x-cix(ipfReq(j)%ip)
            ym1 = ipfReq(j)%y-ciy(ipfReq(j)%ip)
            zm1 = ipfReq(j)%z-ciz(ipfReq(j)%ip)

            !Note that this assumes you have correctly handled wall boundary conditions in parseMPIlinks!
            if(status(MPI_SOURCE) == myp)then
              !Check Pre-Stream location for solid node conflicts (Ignores x wall)
              if(ibnodes(xm1,ym1,zm1) > 0)then
                mypIpfSend(j) = IBNODES_TRUE
              else
                mypIpfSend(j) = f(ip, ipfReq(j)%x, ipfReq(j)%y, ipfReq(j)%z,s)
              endif
            else
              !Check Pre-Stream location for solid node conflicts (Ignores x wall)
              if(ibnodes(xm1,ym1,zm1) > 0)then
                mymIpfSend(j) = IBNODES_TRUE
              else
                mymIpfSend(j) = f(ip, ipfReq(j)%x, ipfReq(j)%y, ipfReq(j)%z,s)
              endif
            endif
          endif         
        enddo
        deallocate(ipfReq)
      enddo

      ! If corner data was requested from Y neighbors add it to vertical requests
      ! Corner requests are always on the end of the array
      do j=1,ipf_mzptc
        ipf_mzp(ipf_mzpc+j) = ipf_mzpt(j)
      enddo
      do j=1,ipf_mzmtc
        ipf_mzm(ipf_mzmc+j) = ipf_mzmt(j)
      enddo

      !Send request arrays to Z direction neighbors
      call MPI_ISEND(ipf_mzp(1:ipf_mzpc+ipf_mzptc), ipf_mzpc+ipf_mzptc, MPI_IPF_NODE, mzp, 98, MPI_COMM_WORLD, req(3), ierr)
      call MPI_ISEND(ipf_mzm(1:ipf_mzmc+ipf_mzmtc), ipf_mzmc+ipf_mzmtc, MPI_IPF_NODE, mzm, 98, MPI_COMM_WORLD, req(4), ierr)
      !Recieve requests from Z neighbors
      do i = 1, 2
        call MPI_PROBE(MPI_ANY_SOURCE, 98, MPI_COMM_WORLD, status, ierr)
        call MPI_GET_COUNT(status, MPI_IPF_NODE, count, ierr)
        allocate(ipfReq(count))
        call MPI_RECV(ipfReq, count, MPI_IPF_NODE, status(MPI_SOURCE), status(MPI_TAG), MPI_COMM_WORLD, status, ierr)

        if(status(MPI_SOURCE) == mzp)then
          !Allocate send buffers, and adjust cordinates to local grid
          allocate(mzpIpfSend(count))
          zpcount = count
          ipfReq(:)%z = ipfReq(:)%z + lz

          do j=1, count
            ip =  ipfReq(j)%ip
            !Get pre stream location
            xm1 = ipfReq(j)%x-cix(ipfReq(j)%ip)
            ym1 = ipfReq(j)%y-ciy(ipfReq(j)%ip)
            zm1 = ipfReq(j)%z-ciz(ipfReq(j)%ip)

            if(ibnodes(xm1,ym1,zm1) > 0)then
              mzpIpfSend(j) = IBNODES_TRUE
            else
              mzpIpfSend(j) = f(ip, ipfReq(j)%x, ipfReq(j)%y, ipfReq(j)%z,s)
            endif
          enddo
          !Send distribution data back
          call MPI_ISEND(mzpIpfSend, zpcount, MPI_REAL8, mzp, 97, MPI_COMM_WORLD, req(5), ierr)
        else
          !Allocate send buffers, and adjust cordinates to local grid
          allocate(mzmIpfSend(count))
          zmcount = count
          ipfReq(:)%z = ipfReq(:)%z - lz

          do j=1, count
            ip =  ipfReq(j)%ip
            !Get pre stream location
            xm1 = ipfReq(j)%x-cix(ipfReq(j)%ip)
            ym1 = ipfReq(j)%y-ciy(ipfReq(j)%ip)
            zm1 = ipfReq(j)%z-ciz(ipfReq(j)%ip)

            if(ibnodes(xm1,ym1,zm1) > 0)then
              mzmIpfSend(j) = IBNODES_TRUE
            else
              mzmIpfSend(j) = f(ip, ipfReq(j)%x, ipfReq(j)%y, ipfReq(j)%z,s)
            endif
          enddo
          !Send distribution data back
          call MPI_ISEND(mzmIpfSend, zmcount, MPI_REAL8, mzm, 97, MPI_COMM_WORLD, req(6), ierr)
        endif        
        deallocate(ipfReq)
      enddo
      
      !Recieve z neighbor interpolation fluid node data
      do i = 1, 2
        call MPI_PROBE(MPI_ANY_SOURCE, 97, MPI_COMM_WORLD, status, ierr)
        call MPI_GET_COUNT(status, MPI_REAL8, count, ierr)
        allocate(ipfRecv(count))
        call MPI_RECV(ipfRecv, count, MPI_REAL8, status(MPI_SOURCE), status(MPI_TAG), MPI_COMM_WORLD, status, ierr)

        if(status(MPI_SOURCE) == mzp)then
          mzpIpfRecv = ipfRecv(1:ipf_mzpc)
          do j=1, ipf_mzptc !Move requested corner data into Y send buffs
            if(mzpt_index(j)%mpidir == 1)then
              mypIpfSend(mzpt_index(j)%index) = ipfRecv(ipf_mzpc+j)
            else
              mymIpfSend(mzpt_index(j)%index) = ipfRecv(ipf_mzpc+j)
            endif
          enddo
        else
          mzmIpfRecv = ipfRecv(1:ipf_mzmc)
          do j=1, ipf_mzmtc !Move requested corner data into Y send buffs
            if(mzmt_index(j)%mpidir == 1)then 
              mypIpfSend(mzmt_index(j)%index) = ipfRecv(ipf_mzmc+j)
            else
              mymIpfSend(mzmt_index(j)%index) = ipfRecv(ipf_mzmc+j)
            endif
          enddo
        endif        
        deallocate(ipfRecv)
      enddo

      !Send/receive Y interpolation fluid node data to
      call MPI_IRECV(mypIpfRecv,ipf_mypc,MPI_REAL8,myp,0,MPI_COMM_WORLD,req(7),ierr)
      call MPI_IRECV(mymIpfRecv,ipf_mymc,MPI_REAL8,mym,1,MPI_COMM_WORLD,req(8),ierr)

      call MPI_ISEND(mypIpfSend, ypcount, MPI_REAL8, myp, 1, MPI_COMM_WORLD, req(9), ierr)
      call MPI_ISEND(mymIpfSend, ymcount, MPI_REAL8, mym, 0, MPI_COMM_WORLD, req(10), ierr)

      call MPI_WAITALL(10,req,status_array,ierr)

      end subroutine exchng2direct
!=============================================================================
!@subroutine exchng2direct
!@desc Calculates the lubrication forces acting on the particle from either
!      other solid particles or walls
!=============================================================================
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

      !Particle-particle lubrication add-on
      !Parse through local solid particles
      do ii = 1,nps
        xc = yp(1,ii)
        yc = yp(2,ii)
        zc = yp(3,ii)
        id = ipglb(ii)
        
        if(id < npart)then
        !Check other particles with a greater index
        do jj = id+1,npart
          xb = ypglb(1,jj)
          yb = ypglb(2,jj)
          zb = ypglb(3,jj)

          !Use the nearest particle center instead of the real center
          !if((xb - xc) > real(nxh)) xb = xb - real(nx)
          !if((xb - xc) < -real(nxh)) xb = xb + real(nx)

          if((yb - yc) > real(nyh)) yb = yb - real(ny)
          if((yb - yc) < -real(nyh)) yb = yb + real(ny)

          if((zb - zc) > real(nzh)) zb = zb - real(nz)
          if((zb - zc) < -real(nzh)) zb = zb + real(nz)


          dxij = xc - xb
          dyij = yc - yb
          dzij = zc - zb
          dist = sqrt(dxij*dxij + dyij*dyij + dzij*dzij)
          hgap = dist - 2.0*rad

          !If solid particles are too close add some forcing
          !From Feng & Michaelides (2005) JCP 202, pp 20-51, eqn(28)
          if(dist < dist0)then
            dist1 = 0.0
            if(hgap < 0.0) dist1 = -hgap
            fp = (((dist0 - dist) / mingap)**2 / stf0 + dist1 / mingap / stf1)*fpc
            fpd = fp / dist 
          
            flubp0(1,id) = flubp0(1,id) + fpd*dxij          
            flubp0(2,id) = flubp0(2,id) + fpd*dyij          
            flubp0(3,id) = flubp0(3,id) + fpd*dzij          

            flubp0(1,jj) = flubp0(1,jj) - fpd*dxij          
            flubp0(2,jj) = flubp0(2,jj) - fpd*dyij          
            flubp0(3,jj) = flubp0(3,jj) - fpd*dzij          
          end if
          
        end do
        end if
      end do
      
      !Particle-wall lubrication add-on
      dist0 = mingap_w + rad
      !Parse through local solid particles
      do ii = 1,nps
        id = ipglb(ii)
        xd = yp(1,ii)
        !Bottom wall
        dist = xd
        hgap = dist - rad
        !If solid particle is too close to the the wall add some forcing
        if (dist < dist0) then 
          dist1 = 0.0
          if(hgap < 0.0)  dist1 = -hgap
          fp = (((dist0 - dist) / mingap_w)**2 / stf0_w + dist1 / mingap_w / stf1_w)*fpc
          flubp0(1,id) = flubp0(1,id) + fp
         end if

        !Top wall
        dist = real(nx) - xd
        hgap = dist - rad
        !If solid particle is too close to the the wall add some forcing
        if (dist < dist0) then
          dist1 = 0.0
          if(hgap < 0.0) dist1 = -hgap
          fp = (((dist0 - dist) / mingap_w)**2 / stf0_w + dist1 / mingap_w / stf1_w)*fpc
          flubp0(1,id) = flubp0(1,id) - fp
         end if
      end do

      !Sum up forcing on all solid particles across all MPI tasks
      ilen = 3*npart
      call MPI_ALLREDUCE(flubp0, flubp, ilen, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)      

      end subroutine beads_lubforce
!=============================================================================
!@subroutine beads_move
!@desc Uses forcing values to update solid particles position, angular 
!      displacment, velocity, and rotational velocity
!=============================================================================
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

      !Save as previous values for wp and omgp 
      wpp = wp
      omgpp = omgp
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      !Parse through local solid particles
      do ip = 1,nps
        id = ipglb(ip)
        !Combine fluid ounce back forces and lubrication forces
        forcep(:,ip) = fHIp(:,id) + flubp(:,id)
        
        !If gravity is desired in "-x" dir
        !forcep(1,ip) = forcep(1,ip) - volp*rhog

        !Add the pressure-gradient force for the particle in "+y" direction
        forcep(2,ip) = forcep(2,ip) + volp*force_in_y

        !Update particle translational velocity
        if(istep == irelease) forcepp(:,ip) = forcep(:,ip)
        dwdtn = 0.5*(forcep(:,ip) + forcepp(:,ip)) / amp
        if(istep == irelease) dwdt(:,ip) = dwdtn
        wp(:,id) = wpp(:,id) + 0.5*dt*(dwdt(:,ip) + dwdtn)

        !Update particle position
        yp(:,ip) = yp(:,ip) + 0.5*dt*(wp(:,id) + wpp(:,id))        
 
        !Update particle angular velocity
        if(istep == irelease) torqpp(:,ip) = torqp(:,id)
        domgdtn = 0.5*(torqp(:,id) + torqpp(:,ip)) / aip
        if(istep == irelease) domgdt(:,ip) = domgdtn
        omgp(:,id) = omgpp(:,id) + 0.5*dt*(domgdt(:,ip) + domgdtn)

        !Update particle angular displacement
        thetap(:,ip) = thetap(:,ip) + 0.5*dt*(omgp(:,id) + omgpp(:,id))

        !Update d(wp)/dt, d(omgp)/dt
        !dwdt(:,ip) = (wp(:,id) - wpp(:,id)) / dt  
        dwdt(:,ip) = dwdtn 

        !domgdt(:,ip) = (omgp(:,id) - omgpp(:,id)) / dt
        domgdt(:,ip) = domgdtn  

        !Save as previous values for central differencing
        forcepp(:,ip) = forcep(:,ip)
        torqpp(:,ip) = torqp(:,id)
 
        wp0(:,id) = wp(:,id)
        omgp0(:,id) = omgp(:,id)

      end do

      !Save as previous values for ibnodes and isnodes
      ibnodes0 = ibnodes

      !Collect velocity and rotational velocity data for all solid particles
      ilen = 3*npart
      call MPI_ALLREDUCE(wp0,wp,ilen,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(omgp0,omgp,ilen,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

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
!=============================================================================
!@subroutine beads_filling
!@desc Repopulates nodes that previously existed inside a solid particle.
!=============================================================================
      subroutine beads_filling
      use mpi
      use var_inc
      implicit none

      integer id, ix, iy, iz, ipop, ipmx, ii, nghb
      integer i, j, k, ip1, jp1, kp1, n
      integer ibp1,ib0p1

      real xc, yc, zc, xpnt, ypnt, zpnt
      real w1, w2, w3, omg1, omg2, omg3
      real aa, bb, cc, ddt, ddt0, ddt1  
      real xp1, yp1, zp1, xp2, yp2, zp2
      real xx0, yy0, zz0, prod0, prod
      real rho9, u9, v9, w9
      !Temporary array for sending filling request data      
      logical, dimension(nproc*8):: tempFlags

      real, dimension(0:npop-1):: f9, feq9 

      !Gather all fill flags for all MPI tasks
      call MPI_ALLGATHER(localReqData, 8, MPI_LOGICAL, tempFlags, 8, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
      !Gather dumps into 1d array, lets convert it to a 2d for better intuition
      do n=0, nproc-1
        fillMPIrequest(n,:) = tempFlags(n*8+1:n*8+8)
      enddo

      ! Send and recieve data for filling
      call exchngFill

      !Parse through fill nodes
      do n=1, nbfill
        id = idbfill(n)
        xc = ypglbp(1,id)*1.d0      
        yc = ypglbp(2,id)*1.d0    
        zc = ypglbp(3,id)*1.d0   

        i = xbfill(n)
        j = ybfill(n)
        k = zbfill(n)

        xpnt = dfloat(i) - 0.5d0 
        ypnt = dfloat(j) - 0.5d0 + dfloat(indy*ly)   
        zpnt = dfloat(k) - 0.5d0 + dfloat(indz*lz)

        !Use the nearest particle center instead of the real center
        !if((xc - xpnt) > dfloat(nxh)) xc = xc - dfloat(nx)
        !if((xc - xpnt) < -dfloat(nxh)) xc = xc + dfloat(nx)

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

        !Periodicity (Removed because of wall)
        !if(ip1 < 1) ip1 = ip1 + lx
        !if(ip1 > lx) ip1 = ip1 - lx

        ibp1 = ibnodes(ip1,jp1,kp1)
        ib0p1 = ibnodes0(ip1,jp1,kp1)

        IF(ibp1 < 0 .and. ib0p1 < 0)THEN
          rho9 = 0.0 
          u9 = 0.0
          v9 = 0.0
          w9 = 0.0
      
          if(jp1 > ly) then
          f9 = fillRecvYp(:,ip1,kp1)
          else if (jp1 < 1) then
          f9 = fillRecvYm(:,ip1,kp1)
          else
            if(kp1 > lz) then
            f9 = fillRecvZp(:,ip1,jp1)
            else if(kp1 < 1 ) then
            f9 = fillRecvZm(:,ip1,jp1)
            else
            f9 = f(:,ip1,jp1,kp1,s)
            end if
          end if
 
        do ipop = 0,npop-1
          rho9 = rho9 + f9(ipop)  
          u9 = u9 + real(cix(ipop))*f9(ipop)    
          v9 = v9 + real(ciy(ipop))*f9(ipop)   
          w9 = w9 + real(ciz(ipop))*f9(ipop)   
        end do
        call feqpnt(u9,v9,w9,rho9,feq9)
        f(:,i,j,k,s) = f9 - feq9
        ELSE
          f(:,i,j,k,s) = 0.0
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

          ibp1 = ibnodes(ip1,jp1,kp1)
          ib0p1 = ibnodes0(ip1,jp1,kp1)

          IF(ibp1 < 0 .and. ib0p1 < 0)THEN
            nghb = nghb + 1


            if(jp1 > ly) then
            f9 = fillRecvYp(:,ip1,kp1)
            else if (jp1 < 1) then
            f9 = fillRecvYm(:,ip1,kp1)
            else
              if(kp1 > lz) then
              f9 = fillRecvZp(:,ip1,jp1)
              else if(kp1 < 1 ) then
              f9 = fillRecvZm(:,ip1,jp1)
              else
              f9 = f(:,ip1,jp1,kp1,s)
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
        f(:,i,j,k,s) = f(:,i,j,k,s) + feq9
      end do

      end subroutine beads_filling
!=============================================================================
!@subroutine exchngFill
!@desc Exchanges information needed for repopulating fluid nodes in beads_filling
!      This follows a flag based algorithm. If a MPI needs some data from
!      a given neighbor it sets the respective flag to true in beads_links. 
!      These flags are distributed across all MPI task which allows any MPI 
!      task to determine which neighbors it will be sending or recieving data 
!      from. For further details contact Nicholas Geneva.
!=============================================================================
      subroutine exchngFill
      use mpi
      use var_inc
      implicit none

      integer ileny, ilenz, nreq
      logical utempinit, dtempinit
      real, dimension(0:npop-1,lx,0:lz+1):: tmpYp, tmpYm
      integer status_array(MPI_STATUS_SIZE,4), req(4)

      utempinit = .FALSE.
      dtempinit = .FALSE.
      ilenz = npop*lx*ly
      ileny = npop*lx*(lz+2)
      nreq = 0
      !Recieving Z+
      if(fillMPIrequest(myid,1) .or. fillMPIrequest(mym,6) .or. fillMPIrequest(myp,5))then
        nreq = nreq + 1
        utempinit = .TRUE.
        call MPI_IRECV(fillRecvZp,ilenz,MPI_REAL8,mzp,0,MPI_COMM_WORLD,req(nreq),ierr)
      endif
      !Recieving Z-
      if(fillMPIrequest(myid,2) .or. fillMPIrequest(mym,8) .or. fillMPIrequest(myp,7))then
        nreq = nreq + 1
        dtempinit = .TRUE.
        call MPI_IRECV(fillRecvZm,ilenz,MPI_REAL8,mzm,1,MPI_COMM_WORLD,req(nreq),ierr)
      endif
      !Sending Z+
      if(fillMPIrequest(mzp,2) .or. fillMPIrequest(mymzp,8) .or. fillMPIrequest(mypzp,7))then
        nreq = nreq + 1
        call MPI_ISEND(f(:,:,:,lz,s),ilenz,MPI_REAL8,mzp,1,MPI_COMM_WORLD,req(nreq),ierr)
      endif
      !Sending Z-
      if(fillMPIrequest(mzm,1) .or. fillMPIrequest(mymzm,6) .or. fillMPIrequest(mypzm,5))then
        nreq = nreq + 1
        call MPI_ISEND(f(:,:,:,1,s),ilenz,MPI_REAL8,mzm,0,MPI_COMM_WORLD,req(nreq),ierr)
      endif
      
      if(nreq > 0)then
        call MPI_WAITALL(nreq,req,status_array,ierr)
      endif
      nreq = 0

      !Receiving Y-
      if(fillMPIrequest(myid,3))then
        nreq = nreq + 1
        call MPI_IRECV(fillRecvYm,ileny,MPI_REAL8,mym,1,MPI_COMM_WORLD,req(nreq),ierr)
      endif
      !Recieving Y+
      if(fillMPIrequest(myid,4))then
        nreq = nreq + 1
        call MPI_IRECV(fillRecvYp,ileny,MPI_REAL8,myp,0,MPI_COMM_WORLD,req(nreq),ierr)
      endif
      !Sending Y-
      if(fillMPIrequest(mym,4))then
        nreq = nreq + 1
        !If we have corner data from the z neighbors that needs to be sent add it
        if(dtempinit) tmpYm(:,:,0) = fillRecvZm(:,:,1)
        tmpYm(:,:,1:lz) = f(:,:,1,:,s)
        if(utempinit) tmpYm(:,:,lz+1) = fillRecvZp(:,:,1)
        call MPI_ISEND(tmpYm,ileny,MPI_REAL8,mym,0,MPI_COMM_WORLD,req(nreq),ierr)
      endif
      !Sending Y+
      if(fillMPIrequest(myp,3))then
        nreq = nreq + 1
        !If we have corner data from the z neighbors that needs to be sent add it
        if(dtempinit) tmpYp(:,:,0) = fillRecvZm(:,:,ly)
        tmpYp(:,:,1:lz) = f(:,:,ly,:,s)
        if(utempinit) tmpYp(:,:,lz+1) = fillRecvZp(:,:,ly)
        call MPI_ISEND(tmpYp,ileny,MPI_REAL8,myp,1,MPI_COMM_WORLD,req(nreq),ierr)
      endif

      if(nreq > 0)then
        call MPI_WAITALL(nreq,req,status_array,ierr)
      endif

      end subroutine exchngFill
!=============================================================================
!@subroutine feqpnt
!@desc Calculates equilibirum distributions based of given input parameters
!@param u9,v9,w9 = real; macro-scopic velocities at the given point
!@param rho9 = real; density at the given point
!@param feq9 = real array; stores calculated equilibrium distribution
!=============================================================================
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