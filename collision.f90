!===========================================================================
! This subroutine is based on element-wise operations to increase efficiency.
! Tests show that this can reduce the running time compared
! with array operations, due to the elimination of unnecessary operations
! involving the "zero" elements of the sparce matrix [M], which transfers 
! velocity space into moment space.

      subroutine collision_MRT     
      use var_inc
      implicit none 

      real, dimension(0:npop-1) :: f9
      real rho9, ux9, uy9, uz9, ux9s, uy9s, uz9s
      real t1, tl1, tl2, tl3, tl4, tl5, tl6, tl7, tl8, tl9, tl10, tl11,&
           tl12, tl13, tl14, tl15, tl16, tl17, tl18, tl19, tl20, tl21
      real eqm1, eqm2, eqm3, eqm4, eqm5, eqm6, eqm7, eqm8, eqm9, eqm10,&
           eqm11, eqm12, eqm13, eqm14, eqm15 
      real sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10,&
           sum11
      real evlm1, evlm2, evlm3, evlm4, evlm5, evlm6, evlm7, evlm8,     &
           evlm9, evlm10, evlm11, evlm12, evlm13, evlm14, evlm15
      real eqmc1, eqmc2, eqmc3, eqmc4, eqmc5, eqmc6, eqmc7, eqmc8,     &
           eqmc9, eqmc10, eqmc11, eqmc12, eqmc13, eqmc14, eqmc15
      real suma, sumb, sumc, sumd, sume, sumf, sumg, sumh, sumi, sumk, &
           sump, sum67, sum89, sum1011
      integer ix, iy, iz 
      integer ip
      real fx9,fy9,fz9,G1,G2,G3
      real, dimension(0:npop-1) :: Fbar

      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx
      if(ibnodes(ix,iy,iz) < 0)then

        rho9 = rho(ix,iy,iz)
        ux9 = ux(ix,iy,iz)
        uy9 = uy(ix,iy,iz)
        uz9 = uz(ix,iy,iz)
        ux9s = ux9*ux9
        uy9s = uy9*uy9
        uz9s = uz9*uz9

!I'm merging the forcing term with the collision operator
        fx9 = force_realx(ix,iy,iz)
        fy9 = force_realy(ix,iy,iz)
        fz9 = force_realz(ix,iy,iz)
        G3 = ux9*fx9 + uy9*fy9 + uz9*fz9

        Fbar(0) = -G3

        do ip=1,6
          G1 = cix(ip)*fx9 + ciy(ip)*fy9 + ciz(ip)*fz9
          G2 = cix(ip)*ux9 + ciy(ip)*uy9 + ciz(ip)*uz9
          Fbar(ip) = ww1*(3.*G1 + 9.*G1*G2 - 3.*G3)
        enddo

        do ip=7,(npop-1)
          G1 = cix(ip)*fx9 + ciy(ip)*fy9 + ciz(ip)*fz9
          G2 = cix(ip)*ux9 + ciy(ip)*uy9 + ciz(ip)*uz9
          Fbar(ip) = ww2*(3.*G1 + 9.*G1*G2 - 3.*G3)
        enddo

        f9 = f(:,ix,iy,iz) + Fbar(:)/2.

        t1 = ux9s + uy9s + uz9s
        eqm1 = -11.0*rho9 + 19.0*t1
        eqm2 = omegepsl*rho9 + omegepslj*t1
        eqm3 = coef1*ux9
        eqm4 = coef1*uy9
        eqm5 = coef1*uz9
        eqm6 = 2.0*ux9s - uy9s - uz9s
        eqm7 = omegxx*eqm6
        eqm8 = uy9s - uz9s
        eqm9 = omegxx*eqm8
        eqm10 = ux9*uy9
        eqm11 = uy9*uz9
        eqm12 = ux9*uz9
        eqm13 = 0.0
        eqm14 = 0.0
        eqm15 = 0.0

        sum1 = f9(1) + f9(2) + f9(3) + f9(4) + f9(5) + f9(6)        
        sum2 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14) + f9(15) + f9(16) + f9(17) + f9(18)      
        sum3 = f9(7) - f9(8) + f9(9) - f9(10) + f9(11) - f9(12)        &
             + f9(13) - f9(14)
        sum4 = f9(7) + f9(8) - f9(9) - f9(10) + f9(15) - f9(16)        &
             + f9(17) - f9(18)
        sum5 = f9(11) + f9(12) - f9(13) - f9(14) + f9(15) + f9(16)     &
             - f9(17) - f9(18)
        sum6 = f9(1) + f9(2)
        sum7 = f9(3) + f9(4) + f9(5) + f9(6)
        sum8 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14)
        sum9 = f9(15) + f9(16) + f9(17) + f9(18)
        sum10 = f9(3) + f9(4) - f9(5) - f9(6)
        sum11 = f9(7) + f9(8) + f9(9) + f9(10) - f9(11) - f9(12)       &
              - f9(13) - f9(14)
    
        evlm1 = -30.0*f9(0) + coef2*sum1 + coef3*sum2
        evlm2 = 12.0*f9(0) + coef4*sum1 + sum2
        evlm3 = coef4*(f9(1) - f9(2)) + sum3
        evlm4 = coef4*(f9(3) - f9(4)) + sum4
        evlm5 = coef4*(f9(5) - f9(6)) + sum5
        evlm6 = coef5*sum6 - sum7 + sum8 - coef5*sum9
        evlm7 = coef4*sum6 + coef5*sum7 + sum8 - coef5*sum9
        evlm8 = sum10 + sum11
        evlm9 =-coef5*sum10 + sum11
        evlm10 = f9(7) - f9(8) - f9(9) + f9(10)
        evlm11 = f9(15) - f9(16) - f9(17) + f9(18)
        evlm12 = f9(11) - f9(12) - f9(13) + f9(14)
        evlm13 = f9(7) - f9(8) + f9(9) - f9(10) - f9(11) + f9(12)      &
               - f9(13) + f9(14)
        evlm14 =-f9(7) - f9(8) + f9(9) + f9(10) + f9(15) - f9(16)      &
               + f9(17) - f9(18)
        evlm15 = f9(11) + f9(12) - f9(13) - f9(14) - f9(15) - f9(16)   &
               + f9(17) + f9(18)
   
        eqmc1 = evlm1 - s1*(evlm1 - eqm1)
        eqmc2 = evlm2 - s2*(evlm2 - eqm2)
        eqmc3 = evlm3 - s4*(evlm3 - eqm3)
        eqmc4 = evlm4 - s4*(evlm4 - eqm4)
        eqmc5 = evlm5 - s4*(evlm5 - eqm5)
        eqmc6 = evlm6 - s9*(evlm6 - eqm6)
        eqmc7 = evlm7 - s10*(evlm7 - eqm7)
        eqmc8 = evlm8 - s9*(evlm8 - eqm8)
        eqmc9 = evlm9 - s10*(evlm9 - eqm9)
        eqmc10 = evlm10 - s13*(evlm10 - eqm10)
        eqmc11 = evlm11 - s13*(evlm11 - eqm11)
        eqmc12 = evlm12 - s13*(evlm12 - eqm12)
        eqmc13 = evlm13 - s16*(evlm13 - eqm13)
        eqmc14 = evlm14 - s16*(evlm14 - eqm14)
        eqmc15 = evlm15 - s16*(evlm15 - eqm15)


        tl1 = val1i*rho9
        tl2 = coef2*val2i*eqmc1
        tl3 = coef3*val2i*eqmc1
        tl4 = coef4*val3i*eqmc2
        tl5 = val3i*eqmc2
        tl6 = val4i*ux9
        tl7 = val5i*eqmc3
        tl8 = val4i*uy9
        tl9 = val5i*eqmc4
        tl10 = val4i*uz9
        tl11 = val5i*eqmc5
        tl12 = val6i*eqmc6
        tl13 = val7i*eqmc7
        tl14 = val8i*eqmc8
        tl15 = val9i*eqmc9
        tl16 = -coef4i*eqmc10
        tl17 = -coef4i*eqmc11
        tl18 = -coef4i*eqmc12
        tl19 = coef3i*eqmc13
        tl20 = coef3i*eqmc14
        tl21 = coef3i*eqmc15


        f9(0) = tl1 - 30.0*val2i*eqmc1 + val8*val3i*eqmc2

        suma = tl1 + tl2 + tl4
        sumb = tl1 + tl3 + tl5
        sumc = tl6 + coef4*tl7
        sumd = coef5*tl12 + coef4*tl13
        sume = tl8 + coef4*tl9
        sumf = -tl12 + coef5*tl13 + tl14 - coef5*tl15
        sumg = tl10 + coef4*tl11
        sumh = -tl12 + coef5*tl13 - tl14 + coef5*tl15
 
        sumi = tl12 + tl13 + tl14 + tl15
        sumk = tl12 + tl13 - tl14 - tl15
 
        sump = -coef5*tl12 - coef5*tl13

        sum67 = tl6 + tl7
        sum89 = tl8 + tl9
        sum1011 = tl10 + tl11

        f9(1) = suma + sumc + sumd
        f9(2) = suma - sumc + sumd
        f9(3) = suma + sume + sumf
        f9(4) = suma - sume + sumf
        f9(5) = suma + sumg + sumh
        f9(6) = suma - sumg + sumh

        f9(7) = sumb + sum67 + sum89 + sumi + tl16 + tl19 - tl20
        f9(8) = sumb - sum67 + sum89 + sumi - tl16 - tl19 - tl20
        f9(9) = sumb + sum67 - sum89 + sumi - tl16 + tl19 + tl20
        f9(10) = sumb - sum67 - sum89 + sumi + tl16 - tl19 + tl20
 
        f9(11) = sumb + sum67 + sum1011 + sumk + tl18 - tl19 + tl21
        f9(12) = sumb - sum67 + sum1011 + sumk - tl18 + tl19 + tl21
        f9(13) = sumb + sum67 - sum1011 + sumk - tl18 - tl19 - tl21
        f9(14) = sumb - sum67 - sum1011 + sumk + tl18 + tl19 - tl21
 
        f9(15) = sumb + sum89 + sum1011 + sump + tl17 + tl20 - tl21
        f9(16) = sumb - sum89 + sum1011 + sump - tl17 - tl20 - tl21
        f9(17) = sumb + sum89 - sum1011 + sump - tl17 + tl20 + tl21
        f9(18) = sumb - sum89 - sum1011 + sump + tl17 - tl20 + tl21

        f(:,ix,iy,iz) = f9 + 0.5*Fbar

      end if 
      end do
      end do
      end do

      end subroutine collision_MRT
!===================================================================

      subroutine rhoupdat 
      use var_inc
      implicit none 

!      rho = f(0,:,:,:)+f(1,:,:,:)+f(2,:,:,:)+f(3,:,:,:)+f(4,:,:,:)     &
!          +f(5,:,:,:)+f(6,:,:,:)+f(7,:,:,:)+f(8,:,:,:)+f(9,:,:,:)     &
!         +f(10,:,:,:)+f(11,:,:,:)+f(12,:,:,:)+f(13,:,:,:)+f(14,:,:,:)&
!        +f(15,:,:,:)+f(16,:,:,:)+f(17,:,:,:)+f(18,:,:,:)

      integer ip

      rho = f(0,:,:,:)
      do ip = 1,npop-1
        rho = rho + f(ip,:,:,:)
      end do

      end subroutine rhoupdat 
!===================================================================

      subroutine macrovar 
      use var_inc
      implicit none 

      integer ip, i, j, k, id             
      real xc, yc, zc, xpnt, ypnt, zpnt, xx0, yy0, zz0     
      real w1, w2, w3, omg1, omg2, omg3 


      integer  iz,iy,ix
      real  sum1,sum2,sum3,sum4,sum5,sum6,ux9,uy9,uz9
      real  rho9
      real, dimension(0:npop-1) :: f9

      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx
      f9 = f(:,ix,iy,iz)

      sum1 = f9(7) - f9(10)
      sum2 = f9(9) - f9(8)

      sum3 = f9(11) - f9(14)
      sum4 = f9(13) - f9(12)

      sum5 = f9(15) - f9(18)
      sum6 = f9(17) - f9(16)

      ux9 = f9(1) - f9(2) + sum1 + sum2 + sum3 + sum4
      uy9 = f9(3) - f9(4) + sum1 - sum2 + sum5 + sum6
      uz9 = f9(5) - f9(6) + sum3 - sum4 + sum5 - sum6

      rho9 = f9(0)+f9(1)+f9(2)+f9(3)+f9(4)+f9(5)+f9(6)&
           +f9(7)+f9(8)+f9(9)+f9(10)+f9(11)+f9(12)&
          + f9(13)+f9(14)+f9(15)+f9(16)+f9(17)+f9(18)

      ux(ix,iy,iz) = ux9 + force_realx(ix,iy,iz)/2.
      uy(ix,iy,iz) = uy9 + force_realy(ix,iy,iz)/2.
      uz(ix,iy,iz) = uz9 + force_realz(ix,iy,iz)/2.
      rho(ix,iy,iz) = rho9

      enddo
      enddo
      enddo

! for solid particle nodes
      if(ipart)then

      do k = 1,lz
      do j = 1,ly
      do i = 1,lx
      if(ibnodes(i,j,k) > 0)then
        id = isnodes(i,j,k) 

        xpnt = real(i) - 0.5
        ypnt = real(j) - 0.5 + real(indy*ly)
        zpnt = real(k) - 0.5 + real(indz*lz)

        xc = ypglb(1,id)
        yc = ypglb(2,id)
        zc = ypglb(3,id)

! use the nearest particle center instead of the real center
!       if((xc - xpnt) > real(nxh)) xc = xc - real(nx)
!       if((xc - xpnt) < -real(nxh)) xc = xc + real(nx)
!        if(abs(xc -xpnt) <= real(nxh)) xc = xc

        if((yc - ypnt) > real(nyh)) yc = yc - real(ny)
        if((yc - ypnt) < -real(nyh)) yc = yc + real(ny)
!        if(abs(yc -ypnt) <= real(nyh)) yc = yc

        if((zc - zpnt) > real(nzh)) zc = zc - real(nz)
        if((zc - zpnt) < -real(nzh)) zc = zc + real(nz)
!        if(abs(zc -zpnt) <= real(nzh)) zc = zc

        xx0 = xpnt - xc
        yy0 = ypnt - yc
        zz0 = zpnt - zc

        w1 = wp(1,id)
        w2 = wp(2,id)
        w3 = wp(3,id)

        omg1 = omgp(1,id)
        omg2 = omgp(2,id)
        omg3 = omgp(3,id)

        ux(i,j,k) = w1 + (omg2*zz0 - omg3*yy0)
        uy(i,j,k) = w2 + (omg3*xx0 - omg1*zz0)
        uz(i,j,k) = w3 + (omg1*yy0 - omg2*xx0)

        rho(i,j,k) = rhopart
      end if
      end do
      end do
      end do

      end if

     
      end subroutine macrovar
!===================================================================

!==================================================================
      SUBROUTINE FORCING
      use mpi
      use var_inc
      implicit none
      integer ixs,ihh,i,j,k,jj,kk
      real x9,y9,z9

      force_realx(:,:,:) = 0.0
      force_realy(:,:,:) = force_in_y*force_mag
      force_realz(:,:,:) = 0.0

      RETURN
      END SUBROUTINE FORCING
!=================================
      SUBROUTINE FORCINGP
      use mpi
      use var_inc
      implicit none
      integer ixs,ihh,i,j,k,jj,kk
      real x9,y9,z9,Amp0,beta9,gamma9,Tpd, phase9,Tpdp,ixs0

      Tpd = 2000.
      Tpdp = 1500.
      beta9 = 3.0
      gamma9=2.0
      ixs0 = 2
      Amp0 = 40.00*beta9/real(ny)*sin(pi2*real(istep)/Tpd)
!     phase9 = sin(pi2*real(istep)/Tpdp)
      phase9 = 0.25

      force_realx(:,:,:) = 0.0
      force_realy(:,:,:) = force_in_y
      force_realz(:,:,:) = 0.0
! add some perturbation
      ixs = ixs0
      ihh = lxh/2
      do k=1,lz
      kk = k + indz*lz
      z9 = pi2*(real(kk)-0.5)/real(nz)
      do j=1,ly
      jj = j + indy*ly
      y9 = pi2*(real(jj)-0.5)/real(ny)
      do i=1,ihh
      x9 = pi2*(real(i)-0.5)/real(ihh)
      force_realx(i+ixs,j,k) = force_in_y*0.5*Amp0*real(ihh)*(1.-cos(x9))  &
                                 *cos(beta9*y9)*cos(gamma9*z9)
      force_realy(i+ixs,j,k) = force_in_y*(1.0-Amp0*real(ny)/beta9    &
                       *sin(x9)*sin(beta9*y9)*cos(gamma9*z9))
      force_realz(i+ixs,j,k) = force_in_y*0.5*Amp0*real(nz)/gamma9*   &
                        sin(x9)*cos(beta9*y9)*sin(gamma9*z9)
      end do
      end do
      end do

      ihh = lxh/2
      ixs = nx - ixs0 - ihh
      do k=1,lz
      kk = k + indz*lz
      z9 = pi2*(real(kk)-0.5)/real(nz)
      do j=1,ly
      jj = j + indy*ly
      y9 = pi2*( (real(jj)-0.5)/real(ny) + phase9 )
      do i=1,ihh
      x9 = pi2*(real(i)-0.5)/real(ihh)
      force_realx(i+ixs,j,k) = -force_in_y*0.5*Amp0*real(ihh)*(1.-cos(x9))  &
                                 *cos(beta9*y9)*cos(gamma9*z9)
      force_realy(i+ixs,j,k) = force_in_y*(1.0+Amp0*real(ny)/beta9    &
                       *sin(x9)*sin(beta9*y9)*cos(gamma9*z9))
      force_realz(i+ixs,j,k) = -force_in_y*0.5*Amp0*real(nz)/gamma9*   &
                        sin(x9)*cos(beta9*y9)*sin(gamma9*z9)
      end do
      end do
      end do

      RETURN
      END SUBROUTINE FORCINGP

