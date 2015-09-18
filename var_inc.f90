!=============================================================================
!@module var_inc 
!@desc var_inc is a module that houses global variables that are used by
!      different subroutines in the program. The main purpose of this is to
!      allow easy access of certain variables and eleminate excessive
!      parameter passing when calling subroutines. Note that if a variable
!      does not need to be global, do not define it in var_inc and instead
!      define it locally in a subroutine.
!=============================================================================
      module var_inc
      implicit none
      
      !Custom datatypes
      type ipf_node !Interpolation fluid node
      sequence
        integer ip, x, y, z
      endtype

      !Fast Fourier Transfer Paramters (Not used currently)
      integer,parameter:: FFTW_FORWARD = -1, FFTW_BACKWARD = 1
      integer,parameter:: FFTW_REAL_TO_COMPLEX = -1,                   &
                          FFTW_COMPLEX_TO_REAL = 1 
      integer,parameter:: FFTW_ESTIMATE = 0, FFTW_MEASURE = 1
      integer,parameter:: FFTW_OUT_OF_PLACE = 0, FFTW_IN_PLACE = 8,    &
                          FFTW_USE_WISDOM = 16
      integer,parameter:: FFTW_THREADSAFE = 128 
      integer(kind = 8):: plan_RC, plan_CR, plan_F, plan_B  
      
      !Domain size paramters
      integer,parameter:: nx7=200,nx = nx7-1, ny = 2*nx7, nz = nx7
      integer,parameter:: lx = nx
      integer,parameter:: lxh = lx/2, lyh = ny/2
      integer,parameter:: nxh = nx7/2, nyh = ny/2, nzh = nz/2
      integer,parameter:: npop = 19

      !Diagnostic and data output paramters
      integer,parameter:: ndiag = 100, nstat = 100  , nspec=1000
      integer,parameter:: nflowout = 10000, npartout = 1000, ntime = 100
      integer,parameter:: nmovieout = 20000000, nsij = 100    

      !Other fluid/ particle related parameters
      integer,parameter:: iflowseed = 54321
      integer,parameter:: NTAB = 32
      real,parameter:: rho0 = 1.0, rhopart = 1.0
      integer,parameter:: npart = 270 
      real,parameter:: rad = 10.0, mingap = 2.0, mingap_w =2.0
      integer,parameter:: irelease = 10
      integer,parameter:: iprocrate = 2  
      real,parameter:: et0 = 2.354998E+01 

      !MPI, input/output, and runtime related variables
      integer ierr, myid, nproc
      integer nprocY, nprocZ
      integer istep, istep0, istep00, nsteps, istpload, imovie   
      integer lz, ly, lyext, lly, nek, MRTtype, mzp, mzm, istat, iseedf, iyf
      integer indy, indz, myp, mym, mypzp, mypzm, mymzp, mymzm

      integer iseedp, msize, nps, iyp, kpeak
      
      logical newrun, ipart, newinitflow, released

      !General LBM variables and constants
      real force_in_y,ustar,Rstar,ystar,force_mag
      real pi, pi2, anu, visc, tau, u0in 
      real vscale, escale, dscale, tscale 
      real omegepsl, omegepslj, omegxx 
      real s1, s2, s4, s9, s10, s13, s16
      real coef1, coef2, coef3, coef4, coef5, coef3i, coef4i
      real val1, val2, val3, val4, val5, val6, val7, val8, val9,       &
           val1i, val2i, val3i, val4i, val5i, val6i, val7i, val8i, val9i
      real ww0, ww1, ww2  
      real rhoerr, rhoerrmax, rhoepsl , rhomax, rhomin
      real volp, amp, aip, g_lbm, rhog, ws_normalized
      real stf0, stf1, stf0_w, stf1_w
      real time_start, time_end, time_diff, time_max,      &
                       time_lmt, time_buff, time_bond
  
      integer,dimension(0:npop-1):: cix, ciy, ciz, ipopp
      integer,dimension(NTAB):: ivp 
      real,dimension(npop-1):: wwp
      
      integer,dimension(NTAB):: ivf
      real,dimension(6,5,5)  :: a1r,a2r,a3r,b1r,b2r,b3r
      real,dimension(nx+2)   :: kxr
      real,dimension(ny)     :: kyr,kzr
      real,allocatable,dimension(:,:,:):: force_realx,force_realy,force_realz

      integer,allocatable,dimension(:,:,:):: ik2 
      integer,allocatable,dimension(:):: icouples 
      integer,allocatable,dimension(:):: ipglb, ovlpflag, ovlpflagt
      integer,allocatable,dimension(:,:,:):: ibnodes, ibnodes0
      integer,allocatable,dimension(:,:,:):: isnodes, isnodes0 

      integer s
      real,allocatable,dimension(:,:,:,:,:):: f

      !Particle Specific Arrays/ Variables
      integer nbfill, maxbfill
      integer nlink, maxlink
      integer ipf_mymc, ipf_mypc, ipf_mzmc, ipf_mzpc
      integer MPI_IPF_NODE, IBNODES_TRUE

      real,allocatable,dimension(:):: xlink, ylink, zlink
      real,allocatable,dimension(:):: iplink, alink, mlink
      integer,allocatable,dimension(:,:,:):: iblinks
      type(ipf_node), allocatable,dimension(:):: ipf_mym, ipf_myp, ipf_mzm, ipf_mzp

      integer,allocatable,dimension(:):: xbfill, ybfill, zbfill
      integer,allocatable,dimension(:):: idbfill
      logical,allocatable,dimension(:,:):: fillMPIrequest
      logical,allocatable,dimension(:):: localReqData
      real,allocatable,dimension(:,:,:):: fillRecvYm, fillRecvYp    
      real,allocatable,dimension(:,:,:):: fillRecvZm, fillRecvZp

      !Macroscopic variable arrays
      real,allocatable,dimension(:,:,:):: rho, rhop
      real,allocatable,dimension(:,:,:):: ux, uy, uz
      real,allocatable,dimension(:,:,:):: ox, oy, oz
      real,allocatable,dimension(:,:,:):: kx, ky, kz, k2
      real,allocatable,dimension(:,:):: yp, ypglb, ypglbp
      real,allocatable,dimension(:,:):: wp, wpp
      real,allocatable,dimension(:,:):: omgp, omgpp
      real,allocatable,dimension(:,:):: dwdt, domgdt
      real,allocatable,dimension(:,:):: forcep, forcepp 
      real,allocatable,dimension(:,:):: torqp,torqpp
      real,allocatable,dimension(:,:):: thetap 
      real,allocatable,dimension(:,:):: fHIp
      real,allocatable,dimension(:,:):: flubp

      ! note that to make use of FFTW library on bluefire, the complex numbers 
      ! must have a size of at least complex*16, or even higher. For current 
      ! real*4, double complex is equivalent to complex*16. 
      real,allocatable,dimension(:,:,:):: vx, vy, vz
      real,allocatable,dimension(:,:,:):: wx, wy, wz

      !File directory character arrays
      character(len=120):: dirgenr, dirinitflow
      character(len=120):: dirdiag, dirstat, dirprobe   
      character(len=120):: dircntdflow,dircntdflow0,dircntdpart,dircntdpart0
      character(len=120):: dirflowout, dirpartout    
      character(len=120):: dirmoviedata

      !Benchmarking variables
      real bnchstart
      real, allocatable,dimension(:):: collision_MRT_bnch, streaming_bnch, macrovar_bnch
      real, allocatable,dimension(:):: beads_collision_bnch, beads_lubforce_bnch, &
                beads_move_bnch, beads_redistribute_bnch, beads_links_bnch, beads_filling_bnch
      character(len=120):: dirbench, dirbenchmatlab, dirbenchbead, dirbenchflow, dirbenchstat

      end module var_inc
