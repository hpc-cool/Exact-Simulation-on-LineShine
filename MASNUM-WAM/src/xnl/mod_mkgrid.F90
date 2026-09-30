#define VALLOC(N,S) if(allocated(N)) deallocate(N) ; allocate (N S)
module m_mkgrid
  implicit none
  private
  public  m_tools_init,q_init,z_wnumb,q_makegrid !,q_ctrgrid
  public wws_type ,quads,quad_nloc,maq,mkq,klocus
  public naq,nkq,q_a,q_k,q_dk,q_k2,q_ka,q_delta
  public iq_err,iq_sym,mpi_id
  public qf_kn,qf_dn,q_depth
  public q_sig,q_sigr,q_kpow,q_tail,q_lamd,q_kfac,q_lkfac,qf_tail

  real,parameter::q_grav=9.81                    ! gravitational acceleration (Earth = 9.81 m/s^2)
  real,parameter::sqrtg = 3.1321     ! square root of grav
  real,parameter::gsq = q_grav*q_grav    ! square of grav
  ! mathematical constants
  real(8),parameter::pi = 3.14159265358979323846D0  ! circular constant, 3.1415...
  real(8),parameter::pi2 = 2.*pi   ! 2*pi
  real(8),parameter::pi4 = 0.25*pi   ! pi/2
  real(8),parameter::dera = pi/180.  ! conversion from degrees to radians
  real(8),parameter::rade = 180./pi  ! conversion from radians to degrees

  !============== General settings =================
  !============== method control =================
  integer::iq_type=2   !  method for computing the nonlinear interactions
  !                    depending on the value of iq_type a number of settings
  !                    for other processes or schematizations are set in Q_COMPU
  !  == 2, deep water computation with WAM depth scaling based on Herterich and Hasselmann (1980)
  !   iq_geom   = 1
  !   iq_dscale = 1
  !   iq_disp   = 1
  !   iq_cple   = 1

  integer::iq_grid  =3   !  type of spectral grid
  ! used by q_xnl4v4,q_t13v4,q_getlocus in mod_xnl4v5.f90
  ! used by m_tools_init in mod_mkgrid.f90
  !  == 1, sector & symmetric around zero
  !  == 2, sector & symmetric around zero & non-symmetric
  !  == 3, full circle & non-symmetric
  integer::iq_locus =3   !  Option for computation of locus
  ! set by m_tools_init in mod_mkgrid.f90
  ! used by q_polar2 in mod_mkgrid.f90
  !  ==1, explicit polar method with fixed k-step
  !  ==2, explicit polar method with adpative k-stepping
  !  ==3, explicit polar method with geometric k-spacing
  integer::iq_search=0 ! switch to determine search for a proper grid
  ! used by q_xnl4v4 in mod_xnl4v5.f90
  !  == 0, no search is carried out
  !  == 1, search nearest (relative) interaction grid
  integer::iq_sym   =1   ! switch to activate use of symmetry reduction
  !used by q_modify,in mod_mkgrid.f90
  !  == 0, no symmetries are used
  !  == 1, symmetry activated (default)
  integer::iq_syma  =1   ! switch to activate use of symmetry reduction
  !used by q_modify,in mod_mkgrid.f90
  !  == 0, no symmetries are used
  !  == 1, symmetry activated (default)
  integer::iq_compact=1  ! switch to compact data
  ! used by q_makegrid,in mod_mkgrid.f90
  !  == 0, do not compact
  !  == 1, compact data by elimiting zero contribution along locus
  integer::iq_mod =1  !  option to redistribute points on locus
  ! used by q_modify in mod_mkgrid.f90
  !  == 0, Points will be used as computed by tracing algortihm
  !  == 1, Equi-distant spacing on points along locus (NLOC1)

  character(len=60)::q_version ='GurboQuad  Version: 5.03 Build: 59 Date: 2003/09/15 [S]'
  integer::maq           ! number of theta elements in grid matrix
  integer::mkq           ! number of k-elements in grid matrix
  ! set absolute and relative accuracies
  real:: eps_q =1.e-3    ! absolute accuracy for check of Q
  real:: eps_qk=1.e-7    ! absolute accuracy for check of Q
  integer::id_facmax=2   ! Factor for determining range of depth search (Q_SEARCHGRID)
  integer::nlocus0=30    ! preferred number of points on locus
  integer::klocus        ! number of points on locus as stored in quadruplet database
  integer::nlocus        ! number of points on locus, equal to klocus
  integer::iamax         !  maximum difference in indices for sector grids
  integer::mpi_id
  real::q_depth=2000     ! local water depth in m
  real::q_sdepth=-1.     ! saved local water depth in m
  real::q_mindepth=0.1   ! minimum water depth, set in XNL_INIT, used in Q_CTRGRID
  real::q_maxdepth=2000  ! maximum water depth, set in XNL_INIT, used in Q_CTRGRID
  real::q_dstep     =0.1 ! step size for generating BQF files
  real::q_kfac           ! geometric factor between subsequent wave numbers
  real::q_lkfac          ! log geometric factor between subsequent wave numbers
  real::kqmin            ! lowest wave number
  real::kqmax            ! highest wave number
  real::wk_max           ! maximum weight for wave number interpolation, set in Q_INIT
  real::sk_max=50        ! maximum wave number in extended array
  real::qf_dmax=75       ! maximum directional difference between k1 and k3
  real::qf_amax          ! maximum directional difference between k1 and k3
  real::qf_frac=0.1      ! fraction of maximum action density to filter
  real::qf_tail=-5.      ! power of spectral tail of E(f), e.g. -4,, -4.5, -5
  real::qf_krat=6        ! 2.5 maximum ratio of the interacting wave numbers k1 and k3
  real::qk_tail          ! power of spectral tail of N(k), computed from qf_tail
  integer qf_dn          ! maximum diff between ia1,ia3
  integer qf_kn          ! maximum diff between ik1,ik3

  type wws_type 
    !integer nloc        !number of points on locus
    integer ik2,ia2      !lower wave number,direction  index of k2
    integer ik4,ia4      !lower wave number,direction  index of k2
    real wr_k2,wr_k4
    real w1k2,w2k2,w3k2,w4k2 !weights  of k2 
    real w1k4,w2k4,w3k4,w4k4 !weights  of k4 
    real zz !,sym        !base on depth, compound product of cple*ds*sym/jac
  end type wws_type 
  type(wws_type),target,allocatable :: quads(:,:,:)   ! 
  integer, allocatable :: quad_nloc(:,:)   ! number of points on locus
  !====== grid info ==============
  integer nkq     ! number of wave numbers of quad-grid
  integer naq     ! number of angles of quad-grad
  integer ndq     ! number of depth
  real(8)::q_deltad               !  directional spacing of angular grid in degrees
  real(8)::q_delta                !  directional spacing of angular grid in radians
  real(8),allocatable::q_a(:)     !  directions of quadruplet grid in radians
  real(8),allocatable::q_ad(:)    !  directions of quadruplet grid in radians
  real(8),allocatable::q_dk(:)    !  width of wave number bins [1/m]
  real(8),allocatable::q_k(:)     !  wave number grid [1/m]
  real(8),allocatable::q_k2(:)    !  k*k*deltha
  real(8),allocatable::q_ka(:)    !  dk*k*deltha
  real(8),allocatable::q_ar(:)    !  dk*k*deltha
  real(8),allocatable::q_km(:)    !  wave number grid [1/m]
  real(8),allocatable::q_kpow(:)  !  wave number to a certain power, used in filtering
  real(8),allocatable::q_lamd(:)  !  
  real(8),allocatable::q_tail(:) !  
  real   ,allocatable::q_ca(:),q_sa(:)!cos(q_a),sin(q_a) 
  !used by q_chkcons,q_t13v4,xnl_main
  real(8),allocatable::q_sig(:)   !  radian frequencies associated to wave number/depth
  real(8),allocatable::q_sigr(:)   !  radian frequencies associated to wave number/depth
  !used by q_dscale
  real(8),allocatable:: q_dsig(:) !  radian frequencies associated to wave number/depth
  real(8),allocatable::q_cg(:)    !  group velocity (m/s)
  !   result of subroutine Q_weight
  real, allocatable :: sym_loc(:) ! factor for symmetry between k3 and k4

  integer iq_warn   !  counts the number of warnings
  integer iq_err    !  counts the number of errors
  !------------------------------------------------------------------------------
  !  filtering coefficients
  !------------------------------------------------------------------------------
  real kmid        ! wave number at midpoint of locus along symmetry axis
  real kmidx       ! x-component of wave number at midpoint of locus along symmetry axis
  real kmidy       ! y-component of wave number at midpoint of locus along symmetry axis
  ! used by MAKEGRID and subs q_cmplocus,q_modify,
  real q           ! difference of radian frequencies, used in Resio-Tracy method
  real pmag        ! magnitude of P-vector
  real pang        ! angle related of P-vector, Pang = atan2(py,px), (radians)
  real px,py       ! components of difference k1-k3 wave number
  real k1x,k1y     ! components of k1 wave number
  real k2x,k2y     ! components of k2 wave number
  real k3x,k3y     ! components of k3 wave number
  real k4x,k4y     ! components of k4 wave number
  type gridvar
    !used makegrid only
    real wk_k2,wk_k4      ! position of k2 and k4 wave number ! w.r.t. discrete k-grid
    real wr_k2,wr_k4      ! k2 and k4 wave number 
    real wa_k2,wa_k4      ! position of k2 and k4 wave number ! w.r.t. discrete a-grid
    real k2m_mod,k4m_mod  ! k2 and k4 magnitude around locus
    real k2a_mod,k4a_mod  ! k2 and k4 angle around locus
  end type
  integer mlocus   ! maximum number of points on locus for defining arrays
  integer nlocus1  ! number of points on locus as computed in Q_CMPLOCUS
  type(gridvar),target,allocatable:: gv(:)
  real,allocatable:: x2_loc(:)    ! k2x coordinates around locus
  real,allocatable:: y2_loc(:)    ! k2y coordinates around locus
  real,allocatable:: x4_loc(:)    ! k4x coordinates around locus
  real,allocatable:: y4_loc(:)    ! k4y coordinates around locus
  real,allocatable:: x2_mod(:)    ! k2x coordinates along locus
  real,allocatable:: y2_mod(:)    ! k2y coordinates along locus
  real,allocatable:: x4_mod(:)    ! k4x coordinates along locus
  real,allocatable:: y4_mod(:)    ! k4y coordinates along locus
  real,allocatable:: s_loc (:)    ! coordinate along locus
  real,allocatable:: s_mod (:)    ! coordinate along locus
  real,allocatable:: ds_loc(:)    ! step size around locus
  real,allocatable:: ds_mod(:)    ! step size around locus
  real,allocatable:: sym_mod(:)   ! factor for symmetry between k3 and k4
  real,allocatable:: cple_loc(:)  ! coupling coefficient around locus
  real,allocatable:: cple_mod(:)  ! coupling coefficient around locus
  real,allocatable:: jac_loc (:)  ! jacobian term around locus
  real,allocatable:: jac_mod (:)  ! jacobian term around locus
contains
  subroutine m_mkgalloc
    mlocus = 1.3*nlocus0
    VALLOC(x2_loc ,(mlocus))
    VALLOC(y2_loc ,(mlocus))
    VALLOC(x4_loc ,(mlocus))
    VALLOC(y4_loc ,(mlocus))
    VALLOC(x2_mod ,(mlocus))
    VALLOC(y2_mod ,(mlocus))
    VALLOC(x4_mod ,(mlocus))
    VALLOC(y4_mod ,(mlocus))
    VALLOC(s_loc  ,(mlocus))
    VALLOC(s_mod  ,(mlocus))
    VALLOC(ds_loc ,(mlocus))
    VALLOC(ds_mod ,(mlocus))
    VALLOC(sym_mod,(mlocus))
    VALLOC(cple_loc,(mlocus))
    VALLOC(cple_mod,(mlocus))
    VALLOC(jac_loc ,(mlocus))
    VALLOC(jac_mod ,(mlocus))
    VALLOC(gv,(mlocus))
  end subroutine m_mkgalloc
  subroutine m_tools_init(ierror)
    !  1. Purpose:
    !     Initializing module for quadruplets
    !     and setting default settings
    !OUT: qk_tail<qf_tail  
    !OUT: sk_max 
    !OUT: wk_max ,q_sig, q_kpow,q_dsig,<=q_sig nkq 
    !OUT: q_deltad,q_delta,q_a ,1,oag2,iamax,<=naq
    !  2. Method
    !     Conversion of power of spectral tail from E(f) to N(k) using the following
    !     relations:
    !       E(f) ~ f^qf_tail
    !       N(k) ~ k^qk_tail
    !       qk_tail = qf_tail/2 -1
    !     See also Note 13 of G.Ph. van Vledder
    !  5. Called by
    !     XNL_INIT
    !  6. Subroutines used
    !     Z_STEPS
    ! 10. Source code
    integer,intent(out)::ierror
    integer idepth    ! index over water depths
    integer igrid     ! status of quadruplet grid file
    integer iaq,ikq    ! counters for loops over directions and wave numbers
    real ff            ! frequency
    real fqmin      ! lowest frequency in Hz
    real fqmax      ! highest frequency in Hz
    real(8) ::qkfach
    ierror=0
    call m_mkgalloc
    !--------------------------------------------------------------------------------
    ! convert power of E(f) f^qf_tail to power of N(k) k^qk_tail
    ! See Note 13 of G.Ph. van Vledder
    qk_tail = (qf_tail-2.)/2. ! power of spectral tail, of N(k)
    ! qk_tail = -7./2           !jiangxj:N->E
    wk_max = real(nkq+0.9999) ! set maximum wave number index
    ! compute frequency and wave number grid
    ! assume that frequencies are always geometrically spaced,
    ! in the case of deep water this also holds for the wave numbers
    !  Retrieve size of spectral grid from input
    kqmin = q_k(1); kqmax = q_k(nkq)
    q_deltad= 360./real(naq)            ! degrees
    q_delta = q_deltad*dera             ! directional step in radians
    qkfach= (q_k(nkq)/q_k(1))**dble(0.5/(nkq-1.))         ! geometric spacing factor of frequencies
    q_kfac= qkfach*qkfach
    q_lkfac=alog(q_kfac)
    qkfach=qkfach-1./qkfach
    qf_kn=alog(qf_krat)/q_lkfac
    ! =============== D I R E C T I O N S ===============================================
    ! the directions in the array ANGLE are running from 1 to NAQ
    ! for a sector definition the middle direction has index IAREF
    !  compute index IAREF of middle wave direction for sector grids
    qf_dn=qf_dmax/q_deltad
    call gxnldata_alloc
    q_km(1)=q_k(1)/qkfach
    do ikq = 1,nkq
      q_km(ikq+1)=q_k(ikq)*qkfach
      q_dk(ikq)=qkfach*q_k(ikq);
      !q_dk(ikq)=q_km(ikq+1)-qkfach*q_k(ikq-1);
      q_k2(ikq)=q_k(ikq)*q_k(ikq)*q_delta
      q_ka(ikq)=q_dk(ikq)*q_k(ikq)*q_delta
      q_ar(ikq)=q_dk(ikq)*q_delta
    enddo
    !call z_steps(q_k,  q_dk,  nkq)           ! step size of wave numbers
    !full circle & non-symmetric
    !ang=0;1-24;1-13;0,0
    iamax = maq -1+1;
    do iaq=1,naq
      q_ad(iaq) = q_deltad*(iaq-1.)
      q_a(iaq)  = q_ad(iaq)*dera
      q_ca(iaq) = cos(q_ad(iaq)*dera)
      q_sa(iaq) = sin(q_ad(iaq)*dera)
    end do
    !------------------------------------------------------------------------------
    !  Generate interaction grid and coefficients for each valid water depth
    !  Q_CTRGRID  controls grid generation
    !------------------------------------------------------------------------------
    !call q_ctrgrid(2,igrid)
    q_depth = q_maxdepth
  end subroutine m_tools_init
  real function calcg(wkk,wsk,d)
    real(8),intent(in):: wkk,wsk
    real,intent(in):: d
    real dk,cg
    dk=wkk*d;
    if(dk<40)then
      if(dk<0.001)then
        cg=sqrt(q_grav*d)
      else
        cg=0.5*wsk*(1.+2.*dk/sinh(2.*dk))/wkk
      endif
    else
      cg =0.5*wsk/wkk; 
    endif
    calcg=cg
  end function
  real function calwsk(wkk,d)
    real(8),intent(in):: wkk
    real   ,intent(in):: d
    real(8) dk,wsk
    dk=wkk*d;
    if(dk<40)then
      wsk=sqrt(q_grav*wkk*tanh(dk));
    else
      wsk=sqrt(q_grav*wkk)
    endif
    calwsk=wsk
  end function 
  subroutine q_init
    !  5. Called by
    !     m_tools_init,q_xnl4v4,
    !  6. Subroutines used
    !     Z_STEPS
    ! 10. Source code
    !------------------------------------------------------------------------------------------
    !     Local variables
    integer ikq    ! counters for loops over directions and wave numbers
    integer iuerr      ! error indicator for i/o
    real(8)::    dk,wkk,wsk,cg,dwk,dwkh,wsb,wsc
    !------------------------------------------------------------------------------
    !OUT: q_cg,q_kfac,q_k,q_dk,  <=q_sig,q_depth,q_grav
    if(q_depth==q_sdepth)return
    q_sdepth=q_depth
    dwkh=sqrt(q_kfac)
    wsb=calwsk(q_km(1),q_depth)
    do ikq = 1,nkq
      wkk=q_k(ikq); dk=wkk*q_depth
      ! pwk=pwkh*pwkh
      ! a,a*pwk,a*pwk^2,a*pwk^3...
      ! wk*(pwk-1)
      ! wk*(pwk-1/pwk)/2
      ! wk*(pwkh-1/pwkh)
      ![a/pwkh, a*pwkh],[a*pwkh,a*pwkh^3]
      !s(g*a),s(g*a)*pwkh,s(g*a)*pwkh^2,,s(g*a)*pwkh^3
      ! ws*(pwkh-1)       s(g*wk) s(pwk-1)  
      ! ws*(pwkh-1/pwkh)/2 
      ! ws*(pwkq-1/pwkq)   
      ! s(g*wk*pwkh*tanh(d*wk*pwkh))-s(g*wk/pwkh*tanh(d*wk/pwkh))
      ! ws*pwkhh*s(tanh(d*wk*pwkh))-1./pwkhh*s(tanh(d*wk/pwkh))
      ! 
      wsk=calwsk(q_k(ikq),q_depth)
      q_sig(ikq)=wsk;
      q_sigr(ikq)=1./wsk;
      q_cg(ikq)=calcg(wkk,wsk,q_depth)
      wsc=calwsk(q_km(ikq+1),q_depth)
      q_dsig(ikq)=wsc-wsb;wsb=wsc;
    enddo
    !call z_steps(q_sig,q_dsig,nkq)   !  compute step size of sigma's
    !0.99546220 - 0.99547695
    !do ikq = 1,nkq
    !  print*,'AB',ikq,q_cg(ikq),q_dsig(ikq)/q_cg(ikq),q_dk(ikq),q_dsig(ikq)/q_cg(ikq)/q_dk(ikq)
    !enddo
  end subroutine
  subroutine q_makegrid
    implicit none
    !  1. Purpose:
    !     Set-up grid for computation of loci
    !     Generate data file with basic loci for computation of
    !     nonlinear quadruplet interactions
    !  2. Method
    !  3. Parameter list:
    !     Name    I/O  Type  Description
    !  4. Error messages
    !  5. Called by:
    !     Q_CTRGRID
    !  6. Subroutines used
    !     Q_CPMLOCUS
    !     Q_MODIFY
    !     Q_WEIGHT
    !     Q_CHKRES
    !     Q_NEAREST
    !------------------------------------------------------------------------------
    !     Local variables
    integer iloc,jloc            ! counters
    integer iaq,ikq              ! counters
    integer iaq3,ikq3,nkq1  ! counters
    integer amem,kmem            ! index of angle and wave number in grid
    real(8) aa1,aa3,kk1,kk3         ! temporary wave number variables
    integer nzloc                ! counter for non-zero contributions along locus
    integer ik2,ia2              ! index of wave number k2
    integer ik4,ia4              ! index of wave number k4
    real wk,wa                   ! weights
    real w1k2,w2k2,w3k2,w4k2     ! interpolation weights
    real w1k4,w2k4,w3k4,w4k4     ! interpolation weights
    real qq,zz,m1zz,m2zz,m3zz,m4zz
    real szz,m1szz,m2szz,m3szz,m4szz
    type(wws_type),pointer:: pq
    type(gridvar),pointer:: pg
    !-------------------------------------------------------------------------------
    ! spectral interaction grid
    ! initializations
    quad_nloc = -1    ! number of points on all loci
    !  write header to grid file
    !  set range of do loops for computing interaction grid
    nkq1 = 1      ! use only first wave number for k1, since geometric scaling can be used
    !  compute components of reference wave number,
    !  for setting up interaction grid
    aa1   = q_a(1)
    m4zz=0; m4szz=0
    kk1   = q_k(1)
    k1x   = kk1*q_ca(1); k1y   = kk1*q_sa(1)
    m3zz=0;m3szz=0
    k3: do ikq3 = 1,qf_kn+1 !nkq   !
      kk3 = q_k(ikq3)
      m2zz=0;m2szz=0
      a3: do iaq3 = 1,qf_dn+1
        if(iaq3 == 1 .and. ikq3 == 1) cycle
        aa3 = q_a(iaq3)
        k3x = kk3*q_ca(iaq3); k3y = kk3*q_sa(iaq3)
        !   compute locus for a specified combination of k1 and k3
        !if(ikq3==2)call dbghere
        call q_cmplocus()
        if(iq_err/=0) return
        !     redistibute or filter data points along locus
        call q_modify !modi nlocus 
        if(iq_err > 0)return
        !     compute weights for interpolation in computational grid
        call q_weight
        if(iq_err > 0)return
        !    special storing mechanism for interactions per combination of k1 and k3
        !kmem  = (ikq3-1+1) - (1-2*nkq-2)*(1-1)/2;
        kmem = ikq3
        amem = iaq3        !index of direction, ensure that data stored in matrix start at index (1,1)
        !     Convert real indices to integer indexing and real weights
        !    3-----------4 ja2p         w1 = (1-wk)*(1-wa)
        !    |    .      |              w2 = wk*(1-wa)
        !    |. . + . . .| wa2   A      w3 = (1-wk)*wa
        !    |    .      |       |      w4 = wk*wa
        !    |    .      |       wa
        !    |    .      |       |
        !    1-----------2 ja2   V
        !   jk2  wk2  jk2p
        !    <-wk->
        !-------------------------------------------------------------------------------
        nzloc = 0;m1zz=0;m1szz=0
        !print'(a,4i3)','ZZ0:',1,1,iaq3,ikq3
        loc:  do iloc = 1,nlocus
          pg=>gv(iloc)
          ik2  = floor(pg%wk_k2);wk   = pg%wk_k2-ik2
          ia2  = floor(pg%wa_k2);wa   = pg%wa_k2-ia2
          ia2  = mod(ia2+naq,naq) !ZZM index from 0
          w1k2 = (1.-wk)*(1.-wa); w2k2 = wk*(1.-wa)
          w3k2 = (1.-wk)*    wa ; w4k2 = wk*    wa
          ik4  = floor(pg%wk_k4);wk   = pg%wk_k4-real(ik4)
          ia4  = floor(pg%wa_k4);wa   = pg%wa_k4-real(ia4)
          ia4  = mod(ia4+naq,naq) !ZZM index from 0
          w1k4 = (1.-wk)*(1.-wa); w2k4 = wk*(1.-wa)
          w3k4 = (1.-wk)*    wa ; w4k4 = wk*    wa
          !  compute combined tail factor and product of coupling coefficient, step size,
          !  symmetry factor, and tail factor divided by jacobian
          !tfac=sym_mod(iloc) 
          zz= cple_mod(iloc)/jac_mod(iloc)*ds_mod(iloc)
          szz=zz*sym_mod(iloc)
          if(m1zz<abs(zz))m1zz=zz;
          if(m1szz<abs(szz))m1szz=szz;
          !write(*,'(a,5i3,4e12.4," |",4E12.2,f4.1)') 'ZZA:',1,1,iaq3,ikq3,iloc,&
          !  pg%wk_k2,pg%wa_k2,pg%wk_k4,pg%wa_k4,&
          !  zz,cple_mod(iloc),ds_mod(iloc),jac_mod(iloc),sym_mod(iloc)
          !----------------------------------------------------------------------------------------
          !  compact data by elimating zero-contribution on locus
          !----------------------------------------------------------------------------------------
          if(iq_compact==1 .and. abs(zz)>1.e-15) then
            nzloc = nzloc + 1
            jloc  = nzloc
          else
            jloc = iloc
          end if
          !  shift data
          pq=>quads(jloc,kmem,amem)
          !pq%sym  = sym_mod(iloc)
          pq%wr_k2=pg%wr_k2
          pq%wr_k4=pg%wr_k4
          pq%zz   = szz
          pq%ik2  = ik2           ! lower wave number index of k2
          pq%ia2  = ia2           ! lower direction index of k2
          pq%ik4  = ik4           ! lower wave number index of k4
          pq%ia4  = ia4           ! lower direction index of k4
          pq%w1k2 = w1k2         ! weight 1 of k2
          pq%w2k2 = w2k2         ! weight 2 of k2
          pq%w3k2 = w3k2         ! weight 3 of k2
          pq%w4k2 = w4k2         ! weight 4 of k2
          pq%w1k4 = w1k4         ! weight 1 of k4
          pq%w2k4 = w2k4         ! weight 2 of k4
          pq%w3k4 = w3k4         ! weight 3 of k4
          pq%w4k4 = w4k4         ! weight 4 of k4
        end do loc
        if(iq_compact/=1) then
          nzloc = nlocus
        end if
        quad_nloc(kmem,amem) = nzloc               ! store number of points on locus
        !write(*,'(a,6i3,2i5,3e12.2)') 'ZZB:',1,1,iaq3,ikq3,nlocus,nzloc,kmem,amem,m1zz,m1szz,m1zz-m1szz
        if(m2zz<m1zz)m2zz=m1zz;
        if(m2szz<m1szz)m2szz=m1szz;
      end do a3
      !write(*,'(a,3i3,3e12.2)') 'ZZC:',1,1,iaq3,m2zz,m2szz,m2zz-m2szz
      if(m3zz<m2zz)m3zz=m2zz;
      if(m3szz<m2szz)m3szz=m2szz;
    end do k3
    !write(*,'(a,2i3,3e12.2)') 'ZZD:',1,1,m3zz,m3szz,m3zz-m3szz
    if(m4zz<m3zz)m4zz=m3zz;
    if(m4szz<m3szz)m4szz=m3szz;
    !write(*,'(a,i3,3e12.2)') 'ZZE:',1,m4zz,m4szz,m4zz-m4szz
    call q_fillangle
  end subroutine
  subroutine q_fillangle()
    integer kmem,amem,ams,nloc,adif
    type(wws_type),pointer:: t_q(:)
    type(wws_type),pointer:: r_q(:)
    !0-qf_dn(qf_dn+1)(iag2)  => -1 - -qf_dn (qf_dn)(iag2-1) 
    !0=>1 1=>2

    do kmem= 1,mkq
      !adif: 0,1,2,...na-1,na  ,-na ,-na+1,...  -2,-1
      !indx: 1,2,3,...na  ,na+1,na+2, na+3,...ma-1,ma
      do adif = -1,-qf_dn,-1
        !-1 =>maq -2=>maq-1
        amem=maq+1+adif
        !ams=maq-amem+1+1;
        ams=1-adif
        t_q=>quads(:,kmem,amem)
        r_q=>quads(:,kmem,ams)
        nloc = quad_nloc(kmem,ams)
        quad_nloc(kmem,amem)=nloc 
        t_q(1:nloc)%ik2  = r_q(1:nloc)%ik2
        t_q(1:nloc)%ik4  = r_q(1:nloc)%ik4
        t_q(1:nloc)%wr_k2= r_q(1:nloc)%wr_k2
        t_q(1:nloc)%wr_k4= r_q(1:nloc)%wr_k4
        t_q(1:nloc)%ia2  =-r_q(1:nloc)%ia2 
        t_q(1:nloc)%ia4  =-r_q(1:nloc)%ia4 

        t_q(1:nloc)%w1k2 = r_q(1:nloc)%w3k2
        t_q(1:nloc)%w2k2 = r_q(1:nloc)%w4k2
        t_q(1:nloc)%w3k2 = r_q(1:nloc)%w1k2
        t_q(1:nloc)%w4k2 = r_q(1:nloc)%w2k2

        t_q(1:nloc)%w1k4 = r_q(1:nloc)%w3k4
        t_q(1:nloc)%w2k4 = r_q(1:nloc)%w4k4
        t_q(1:nloc)%w3k4 = r_q(1:nloc)%w1k4
        t_q(1:nloc)%w4k4 = r_q(1:nloc)%w2k4
        t_q(1:nloc)%zz   = r_q(1:nloc)%zz
      enddo
    enddo
  end subroutine
  subroutine q_cmplocus()
    implicit none
    !  1. Purpose:
    !     Compute locus function used for the determination of the
    !     resonnance condition
    !  2. Method
    !     See ALKYON, 1999
    !  3. Parameter list:
    !Type   I/O          Name          Description
    !----------!----------------------------------------------------------------------------
    !  4. Error messages
    !  5. Called by:
    !    Q_MAKEGRID
    !-------------------------------------------------------------------------------
    !     Local variables
    real ka,kb       ! lower and higher wave number magnitude
    real kw          ! half width of locus
    real loclen      ! estimated length of locus
    real k1m          ! magnitude of wave number k1
    real k3m          ! magnitude of wave number k3
    real pcos,psin    ! cosine and sine of normalize angle of P
    real sang        ! angle of symmytry axis of locus, SANG = PANG +/ pi° (radians)
    real klen         ! total length of line locus for case w1=w3
    real kx_beg       ! x-component at start point
    real ky_beg       ! y-component at start point
    real kx_end       ! x-component at end point
    real ky_end       ! y-component at end point
    real dsp,dsm      ! distances in plus and minus direction
    real sum          ! total length of locus
    real eps          ! local accuracy for determination of roots
    real area1        ! area of locus as computed
    real area2        ! area of locus as based on LOCPOS and ellipse
    real ratio        ! maximum ratio between area1 and area2
    integer ierr      ! local error level
    integer iloc,jloc ! counters along locus
    integer ip1       ! index +1
    integer im1       ! index -1
    integer jj        ! counter
    real w1,w3        ! radian frequencies of wave numbers k1 and k3
    !------------------------------------------------------------------------------
    !  set initial values
    eps     = 10.*epsilon(1.)     ! set accuracy 10 times machine accuracy
    ! compute characteristics of configuration
    k1m  = sqrt(k1x**2 + k1y**2)
    k3m  = sqrt(k3x**2 + k3y**2)
    !DBG
    w1   = x_disper(k1m,q_depth)
    w3   = x_disper(k3m,q_depth)
    q    = w1-w3
    !qk   = z_wnumb(abs(q),q_depth,q_grav)
    !if(q<0)qk=-qk
    !  compute cosine and sine of direction of P-vector
    !  reverse direction for the case q<0
    !  first solution along locus: k2 = k3
    !  check for special case if q = 0
    if (abs(q) < eps_q) then
      call q_loc_w1w3(k1x,k1y,k3x,k3y,nlocus0,x2_loc,y2_loc,x4_loc,y4_loc,s_loc)
      nlocus1 = nlocus0
      ds_loc  = s_loc(2)-s_loc(1)
    else
      !------------------------------------------------------------------------------
      !  compute characteristics of locus, such as its position in
      !  wave number space
      !------------------------------------------------------------------------------
      px   = k1x - k3x; py   = k1y - k3y
      pmag = sqrt(px**2 + py**2)
      pang = atan2(py,px)
      if(q < 0) then
        sang = pang+pi
      else
        sang = pang
      end if
      call q_locpos(ka,kb,kw,loclen)
      if(iq_err/=0)return 
      !  compute position of start and end point for tracing the locus
      pcos = cos(sang); psin = sin(sang)
      kx_beg = ka*pcos; ky_beg = ka*psin
      kx_end = kb*pcos; ky_end = kb*psin
      !  compute position of locus using polar method see Van Vledder (2000)
      call q_polar2(ka,kb,kx_beg,ky_beg,kx_end,ky_end,loclen,ierr)
      ! check area of locus by a simple test  (added 12 June 2003)
      call z_polyarea(x2_loc,y2_loc,nlocus1,area1)
      area2 = pi*(kb-ka)*0.5*kw
      ratio = max(area1/area2,area2/area1)
      if(ratio>1.5 .and. k3m/k1m < 100.) then
        print*,"Error LOCUS:Severe problem in POLAR2, ratio > 1.5";
        stop
      end if
      !  compute position of k4 locus by a simple translation
      do iloc=1,nlocus1
        x4_loc(iloc) = x2_loc(iloc) + px
        y4_loc(iloc) = y2_loc(iloc) + py
      end do
    end if
    !----------------------------------------------------------------------------------
    ! compute characteristics around locus
    !----------------------------------------------------------------------------------
    s_loc(1) = 0.; sum      = 0
    do iloc=1,nlocus1
      !  compute step sizes
      if (abs(q) < eps_q) then
        !  for this case the sum of ds_loc is unequal to s_loc(nlocus1)
        sum = s_loc(nlocus1)
      else
        !   compute indices of previous and next index on locus
        ip1 = iloc+1; if (ip1 > nlocus1) ip1 = 1
        im1 = iloc-1; if (im1 < 1) im1 = nlocus1
        dsp = sqrt((x2_loc(iloc)-x2_loc(ip1))**2 + (y2_loc(iloc)-y2_loc(ip1))**2)
        dsm = sqrt((x2_loc(iloc)-x2_loc(im1))**2 + (y2_loc(iloc)-y2_loc(im1))**2)
        if(iloc < nlocus1) s_loc(iloc + 1) = s_loc(iloc) + dsp
        ds_loc(iloc) = 0.5*(dsp + dsm)
        sum = sum+ds_loc(iloc)
      end if
      k2x = x2_loc(iloc); k2y = y2_loc(iloc)
      k4x = x4_loc(iloc); k4y = y4_loc(iloc)
      !  compute gradient/Jacobian terms along locus
      jac_loc(iloc) = x_jacobian(k2x,k2y,k4x,k4y)
      !  compute coupling coefficients along locus
      cple_loc(iloc) = x_cple(k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,q_depth,q_grav)
    end do
  end subroutine
  !------------------------------------------------------------------------------
  subroutine q_modify
    implicit none
    !  1. Purpose:
    !     Modify points along the locus, such that they are evenly distributed
    !     Only when intented, i.e. when IQ_LOCUS==2
    !  2. Method
    !     Compute new spacing along locus
    !     Redistribute points and coefficient at new spacing using linear interpolation
    !     Output DIA configuration when also lumping active
    !     If no redistribution is needed, then copy relevant data
    !  3. Parameter list:
    !     Name    I/O  Type  Description
    !  4. Error messages
    !  5. Called by:
    !     Q_CMPLOCUS
    !  6. Subroutines used
    !     Q_SYMMETRY
    !     Z_INTP1
    !------------------------------------------------------------------------------
    !     Local parameters
    integer ierr         ! error indicators
    integer nold,nnew    ! old and new number of points on locus
    integer iold,inew    ! counter for loop along points
    integer iloc         ! counter for loop along locus
    real k2a,k2m         ! angle (deg) and wave number magnitude of wave number k2
    real k4a,k4m         ! angle (deg) and wave number magnitude of wave number k4
    real w2,w4           ! radian frequencies of wave numbers
    real dk13,dk14       ! difference wave number
    real dsnew,slen      ! new step size and length of locus
    real diold           ! 'real' old number of indices between succeeding lumped bins
    real dinew           ! 'real' new number of indices between succeeding lumped bins
    real, allocatable :: sold(:)     ! old coordinate along locus
    real, allocatable :: snew(:)     ! new coordinate along locus
    type(gridvar),pointer:: pg
    !------------------------------------------------------------------------------
    !  do not modify data when IQ_MOD==0
    !------------------------------------------------------------------------------
    sym_mod=1.
    if(iq_mod==0) then
      nlocus   = nlocus1
      x2_mod   = x2_loc
      y2_mod   = y2_loc
      x4_mod   = x4_loc
      y4_mod   = y4_loc
      s_mod    = s_loc
      ds_mod   = ds_loc
      jac_mod  = jac_loc
      cple_mod = cple_loc
      call q_symmetry(k1x,k1y,k3x,k3y,x4_mod,y4_mod,sym_mod,nlocus)
    else
      !------------------------------------------------------------------------------
      ! Modify spacing along locus
      !------------------------------------------------------------------------------
      nold = nlocus1
      ! close locus by adding one point, equal to first point
      ! only for normal locus
      if(abs(q)>eps_q) nold  = nold+1
      !------------------------------------------------------------------------------
      ! Determine new number of points along locus
      !------------------------------------------------------------------------------
      nnew  = nlocus0
      allocate (sold(nold),snew(nlocus0))
      !------------------------------------------------------------------------------
      !  Compute circumference of locus, distinguish 2 case, open or closed
      !------------------------------------------------------------------------------
      if(abs(q)<eps_q) then
        slen = s_loc(nold)
        sold = s_loc
        dsnew  = slen/real(nnew-1.)
        do inew=1,nnew
          snew(inew) = (inew-1)*dsnew
        end do
      else
        slen = 0
        do iold=1,nold-1               ! loop length minus one, since locus is closed
          sold(iold) = s_loc(iold)
          slen = slen + ds_loc(iold)
        end do
        !------------------------------------------------------------------------------
        !  close locus by copying first value in last value
        !------------------------------------------------------------------------------
        sold(nold)     = slen
        x2_loc(nold)   = x2_loc(1)
        y2_loc(nold)   = y2_loc(1)
        x4_loc(nold)   = x4_loc(1)
        y4_loc(nold)   = y4_loc(1)
        jac_loc(nold)  = jac_loc(1)
        cple_loc(nold) = cple_loc(1)
        dsnew  = slen/real(nnew)
        do inew=1,nnew
          snew(inew) = (inew-1.)*dsnew
        end do
      end if
      ! compute new spacing along loci and coordinates along locus
      ! Gauss-Legendre integration
      ds_mod = dsnew
      ! Compute characteristics of locus for special case q=0
      if(abs(q)<1.e-5) then
        call z_intp1(sold,x2_loc,snew,x2_mod,nold,nnew,ierr)
        call z_intp1(sold,y2_loc,snew,y2_mod,nold,nnew,ierr)
        call z_intp1(sold,x4_loc,snew,x4_mod,nold,nnew,ierr)
        call z_intp1(sold,y4_loc,snew,y4_mod,nold,nnew,ierr)
        call z_intp1(sold,s_loc,snew,s_mod,nold,nnew,ierr)
        call q_symmetry(k1x,k1y,k3x,k3y,x4_mod,y4_mod,sym_mod,nnew)
        ! ---  lumping along locus --------------------------------------------
        call z_intp1(sold,jac_loc,snew,jac_mod,nold,nnew,ierr)
        call z_intp1(sold,cple_loc,snew,cple_mod,nold,nnew,ierr)
        !------------------------------------------------------------------------------------------------
        !  compute characteristics for closed locus
        !------------------------------------------------------------------------------------------------
      else
        call z_intp1(sold,x2_loc,snew,x2_mod,nold,nnew,ierr)
        call z_intp1(sold,y2_loc,snew,y2_mod,nold,nnew,ierr)
        call z_intp1(sold,x4_loc,snew,x4_mod,nold,nnew,ierr)
        call z_intp1(sold,y4_loc,snew,y4_mod,nold,nnew,ierr)
        call z_intp1(sold,s_loc,snew,s_mod,nold,nnew,ierr)
        call q_symmetry(k1x,k1y,k3x,k3y,x4_mod,y4_mod,sym_mod,nnew)
        !  ----- Lumping along locus -----------------------------------
        call z_intp1(sold,jac_loc,snew,jac_mod,nold,nnew,ierr)
        call z_intp1(sold,cple_loc,snew,cple_mod,nold,nnew,ierr)
      end if
      nlocus = nnew
      deallocate(sold,snew)
    end if
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    !!  compute symmetry factor for reducing computational load
    !!
    !!call q_symmetry(k1x,k1y,k3x,k3y,x4_mod,y4_mod,sym,nnew)
    !!
    do iloc=1,nlocus
      pg=>gv(iloc)
      k2x = x2_mod(iloc)
      k2y = y2_mod(iloc)
      k4x = x4_mod(iloc)
      k4y = y4_mod(iloc)
      k2m = sqrt(k2x**2 + k2y**2)
      k4m = sqrt(k4x**2 + k4y**2)
      k2a = atan2(k2y,k2x)*rade
      k4a = atan2(k4y,k4x)*rade
      pg%k2m_mod = k2m
      pg%k4m_mod = k4m
      pg%k2a_mod = k2a
      pg%k4a_mod = k4a
    end do
  end subroutine
  !------------------------------------------------------------------------------
  subroutine q_locpos(ka,kb,kw,loclen)
    implicit none
    !  1. Purpose:
    !     Compute characteristics of locus used to optimize its acutal computation
    !  2. Method
    !  3. Parameter list:
    !Type    I/O          Name     Description
    !-----------------------------------------------------------------
    real, intent (out) :: ka     ! minimum k along symmetry axis
    real, intent (out) :: kb     ! maximum k along symmetry axis
    real, intent (out) :: kw     ! half width of locus at midpoint
    real, intent (out) :: loclen ! estimated length of locus
    !  4. Error messages
    !  5. Called by:
    !     Q_CMPLOCUS
    !  6. Subroutines used
    !     z_zero2   Root finding method
    !     x_locus1  Function of locus geometry, along symmetry axis
    !     x_locusk2  Function of locus geometry, perpendicular to symmetry axis
    !     x_flocus  Locus function
    !  7. Remarks
    !  8. Structure
    !  9. Switches
    !     /S  enable subroutine tracing
    !     /T  enable test output
    ! 10. Source code
    !------------------------------------------------------------------------------
    !     Local variables
    real kp           ! wave number at peak
    real kpx,kpy      ! wave number at peak maximum
    real zp           ! value of locus function at maximum
    real zz1,zz2      ! intermediate function values in interation process
    real kk1,kk2      ! start values for finding root of locus equation
    real kk1x,kk1y    ! wave number components at one side of root
    real kk2x,kk2y    ! wave number components at other side of root
    real beta1,beta2  ! parameters specifying cross component
    real betaw        ! parameter specifying iterated cross component
    real a1,a2,b1,b2  ! constants in polynomial approximation of elliptic function
    real aa,bb,mm1,mm2! semi-major exis of ellips and derived parameters
    real eps         ! local machine accuracy for determination of roots
    real bacc        ! accuracy for determination of beta
    real kacc        ! accuracy for determination of wave number roots
    real qs          ! (w1-w3)/sqrt(g)
    real qsq         ! gs^2
    ! Function declaration
    !real z_root2     ! root finding using Ridders method
    integer ierr     ! local error indicator, used in function Z-ZERO1
    integer iter     ! local iteration number
    integer maxiter  ! maximum number of iteration for determining starting points
    !  function declarations
    !!real x_flocus                 ! 2-d locus function
    !---------------------------------------------------------------------------------
    !  assign test options
    !  set initial values
    eps     = epsilon(1.)         ! determine machine accurcy
    maxiter = 20                  ! maximum number of iterations
    ! compute location of maximum, located at k_2 = P
    !q=|w1|-|w3|
    !(px,py)= k1-k3;pmag=|k1-k3|
    kpx  = -px; kpy  = -py
    kp   = pmag; !sqrt(px**2 + py**2)
    zp   = x_locus1(kp) !k=kp; |w1|-|w3|+|w(|k|)|-|w(|k|+|k1-k3|)|
    ! find location of points A and B on locus
    ! for deep water, explicit relations are available
    qs = q/sqrtg; qsq  = qs*qs
    ka = (0.5*(-qs+sqrt(2.0*pmag-qsq)))**2
    if(qs < 0) then
      kb = ((pmag+qsq)/(2.*qs))**2
    else
      kb = ((pmag-qsq)/(2.*qs))**2
      ka = -ka
    end if
    ! compute position of mid point
    kmid = 0.5*(ka+kb); 
    kmidx = kmid*cos(pang); kmidy = kmid*sin(pang)
    if(q < 0) then
      kmidx =-kmidx; kmidy =-kmidy
    end if
    ! compute width of locus near mid point of locus
    ! set starting values for determination of crossing point
    beta1 = 0. ; beta2 = 0.5;
    kk1x=kmidx; kk1y  = kmidy
    zz1=x_flocus(kk1x,kk1y);!k=kk1; |w1|-|w3|+|w(k)|-|w(k+k1-k3)|
    kk2x=kmidx-beta2*py;kk2y=kmidy+beta2*px
    zz2=x_flocus(kk2x,kk2y);!k=kk2; |w1|-|w3|+|w(k)|-|w(k+k1-k3)|
    do iter=0,maxiter
      if(zz1*zz2 <= 0 )exit
      kk2x=kmidx-beta2*py;kk2y=kmidy+beta2*px
      zz2 =x_flocus(kk2x,kk2y)
      beta2 = beta2*2
    end do
    ! call Ridders method to locate position of zero-crossing
    bacc = 10.*max(beta1,beta2)*eps
    ! k=kmid+l*rot((k1-k3),-90); |w1|-|w3|+|w(k)|-|w(k+k1-k3)|
    betaw = z_root2(x_locus2,beta1,beta2,bacc,ierr)
    !kwx = kmidx - betaw*py; kwy = kmidy + betaw*px
    kw  = betaw*pmag
    ! estimate circumference of locus, assuming it to be an ellips
    ! estimate axis, this seems to be a rather good estimate
    aa = 0.5*abs(ka-kb); bb = kw
    a1 = 0.4630151;  a2 = 0.1077812;
    b1 = 0.2452727;  b2 = 0.0412496;
    if (aa > bb) then
      mm1 = (bb/aa)**2; loclen = 4.*aa
    else
      mm1 = (aa/bb)**2; loclen = 4.*bb
    end if
    mm2=mm1*mm1
    loclen =loclen*((1.+a1*mm1+a2*mm2)-(b1*mm1+b2*mm2)*log(mm1))
  end subroutine
  !------------------------------------------------------------------------------
  real function x_flocus(kxx,kyy) ! |w1|-|w3|+|w(k)|-|w(k+k1-k3)|
    implicit none
    !  1. Purpose: 
    !     Compute locus function used for the determination of the
    !     resonance condition
    !  2. Method
    !     Explicit function evaluation
    !  3. Parameter list:
    !Type    I/O         Name       Description
    !-----------------------------------------------------------
    real, intent(in) ::  kxx      !  x-component of wave number
    real, intent(in) ::  kyy      !  y-component of wave number
    real w2,w4    ! radian frequencies of wave numbers k2 and k4
    !  5. Called by:
    !     Q_LOCPOS
    w2 = sqrtg * (kxx**2 + kyy**2)**(0.25)
    w4 = sqrtg * ((kxx+px)**2 + (kyy+py)**2)**(0.25)
    x_flocus = q + w2 - w4
    return
  end function
  real function x_locus1(k2) ! |w1|-|w3|+|w(|k|)|-|w(|k|+|k1-k3|)|
    implicit none
    !  1. Purpose:
    !     Compute locus function along symmetry axis
    !  2. Method
    !     See ALKYON, 1999
    !  3. Parameter list:
    !Type   I/O         name    Description
    !-------------------------------------------------------
    real, intent(in) :: k2   !  Magnitude of wave number k2
    !  5. Called by:
    !     Q_LOCPOS
    !  7. Remarks
    !     The routine assumes that w1 < w3 or q<0
    !     implying that the directions of k2 and P are opposite
    !  8. Structure
    !  9. Switches
    ! 10. Source code
    !------------------------------------------------------------------------------
    !     Local variables
    real k4       ! wave number magnitudes of k4
    real w2,w4    ! radian frequencies of wave numbers k2 and k4
    k4 = abs(-pmag+k2)
    w2 = sqrtg * sqrt(k2)
    w4 = sqrtg * sqrt(k4)
    x_locus1 = q + w2 - w4
    return
  end function
  !------------------------------------------------------------------------------
  real function x_locus2(lambda) ! k=kmid+l*rot((k1-k3),-90); |w1|-|w3|+|w(k)|-|w(k+k1-k3)|
    implicit none
    !  1. Purpose:
    !     Compute locus function perpendicluar to symmetry axis
    !  2. Method
    !     See ALKYON, 1999
    !  3. Parameter list:
    !     Name    I/O  Type  Description
    real, intent(in) ::  lambda
    !  4. Error messages
    !  5. Called by:
    !     Q_LOCPOS
    !  7. Remarks
    !     The routine assumes that w1 < w3 or q<0
    !     implying that the directions of k2 and P are opposite
    !  8. Structure
    !  9. Switches
    ! 10. Source code:
    !------------------------------------------------------------------------------
    !     local variables
    real kk2x,kk2y,kk2m  ! wave number components and magnitude for k2
    real kk4x,kk4y,kk4m  ! wave number components and magnitude for k4
    real w2,w4           ! radian frequencies of wave numbers k2 and k4
    real z               ! function value
    !px,py <= K1-K3
    kk2x = kmidx - lambda*py  
    kk2y = kmidy + lambda*px
    kk2m = sqrt(kk2x**2 + kk2y**2) ! |KMID-lm*(K1-K3)|
    kk4x = kk2x + px
    kk4y = kk2y + py
    kk4m = sqrt(kk4x**2 + kk4y**2)
    w2 = sqrtg * sqrt(kk2m)
    w4 = sqrtg * sqrt(kk4m)
    x_locus2 = q + w2 - w4
    return
  end function
  !------------------------------------------------------------------------------
  real function x_disper(k,d)
    implicit none
    !  1. Purpose:
    !     Compute radian frequency for a given wave number and water depth
    !  2. Method
    !     wave number is computed as:
    !     1) deep water
    !     2) finite depth linear dispersion relation
    !     3) finited depth non-linear dispersion relation (NOT YET implemented)
    !  3. Parameter list:
    ! Type    I/O          Name       Description
    !----------------------------------------------------------------------
    real, intent(in)   ::   k   !   wave number
    real, intent(in)   ::   d   !   water depth in m
    !  4. Error messages
    !  5. Called by:
    !     Q_CHKRES
    ! 10. Source code
    !------------------------------------------------------------------------------
    !     Local variables
    real kd      ! k*d
    kd = k * d
    x_disper = sqrt(q_grav*k)
  end function
  real function z_wsig(wkk,d,grav)
    implicit none
    !  1. Purpose:
    !     Compute wave number k for a given radian frequency and water depth
    !  2. Method
    !     finite depth linear dispersion relation, using a Pade approximation
    real, intent(in)  ::  wkk    ! radian frequency (rad)
    real, intent(in)  ::  d      ! water depth (m)
    real, intent(in)  ::  grav   ! graviational acceleration (m/s^2)
    !  6. Remarks
    !     The Pade approximation has been described in Hunt, 198.
    real dk,wsk
    if(d<=0 .or. wkk<= 0.) then
      wsk = -10.
    else
      dk=wkk*d;
      if(dk<40)then
        wsk=sqrt(grav*wkk*tanh(dk));
      else
        wsk=sqrt(grav*wkk)
      endif
    endif
    z_wsig=wsk
  end function
  real function z_wnumb(w,d,grav)
    implicit none
    !  1. Purpose:
    !     Compute wave number k for a given radian frequency and water depth
    !  2. Method
    !     finite depth linear dispersion relation, using a Pade approximation
    real, intent(in)  ::  w      ! radian frequency (rad)
    real, intent(in)  ::  d      ! water depth (m)
    real, intent(in)  ::  grav   ! graviational acceleration (m/s^2)
    !  6. Remarks
    !     The Pade approximation has been described in Hunt, 198.
    real x,xx,y,omega
    if(d<=0 .or. w<= 0.) then
      z_wnumb = -10.
    else
      !wk*tanh(d*wk)=sigma*sigma/g
      omega   = w**2/grav
      y       = omega*d
      xx      = y*(y+1./(1.+y*(0.66667+y*(0.35550+y*(0.16084+y*(0.06320+y* &
        (0.02174+y*(0.00654+y*(0.00171+y*(0.00039+y*0.00011))))))))))
      x       = sqrt(xx)
      z_wnumb = x/d
    end if
  end function
  !------------------------------------------------------------------------------
  real function x_cosk(k)
    !  1. Purpose:
    !     Compute cosine of points on locus for given wave number k
    !  2. Method
    !     Explicit polar method, see Van Vledder 2000, Monterey paper
    !     Optionally using a fixed k-step, geometric k-step or adaptive stepping
    real, intent(in) :: k  ! wave number along symmetry axis of locus
    !  5. Called by:
    !     Q_POLAR
    !  6. Subroutines used:
    !     Z_WNUMB   computation of wave number
    !  7. Remarks
    !     The variables q, pmag and q_depth are accessed from module m_xnldata
    !------------------------------------------------------------------------------
    !     Local variables
    real qq    ! constant in direct polar method qq=q/sqrt(g)
    real wk    ! intemediate radian frequency
    real kz    ! intermediate wave number
    !------------------------------------------------------------------------------
    !wk=wsk +wsk1-wsk3
    !kz=s2k(wk)
    !wk = x_disper(k,q_depth)+q !wsk+d_wsk 
    !kz = z_wnumb(wk,q_depth,q_grav) ! =>k 

    qq = q/sqrt(q_grav)
    x_cosk = ((qq+sqrt(k))**4 - k**2 - pmag**2)/(2.*k*pmag)
    x_cosk = min(1.,max(-1.,x_cosk))
  end function x_cosk
  !------------------------------------------------------------------------------
  !OUT: wa_k2,wa_k4 <= k2a_mod,k4a_mod
  !OUT: wk_k2, <=k2m_mod,q_k
  !OUT: wk_k4, <=k4m_mod,q_k
  subroutine q_weight
    !use m_mkgrid,only: wa_k2,wa_k4,k2a_mod,k4a_mod
    !use m_mkgrid,only: wk_k2,k2m_mod,q_k
    !use m_mkgrid,only: wk_k4,k4m_mod !,q_k
    !use m_mkgrid,only:nlocus,nkq,q_deltad,qk_tail,sk_max,wk_max
    implicit none
    !  1. Purpose:
    !     Compute interpolation weights of locus
    !  2. Method
    !     Compute position of wave number in wave number grid
    !     Usable for linear interpolation
    !  3. Parameter list:
    !     Name    I/O  Type  Description
    !  4. Error messages
    !  5. Called by:
    !     Q_MAKEGRID
    !  6. Subroutines used
    !  7. Remarks
    !     The tail factors wt_k2 and wt_k4 are valid for the decay of the action density spectrum
    !     N(kx,ky). With p (qk_tail) the power p of the tail of the N(k) spectrum, and q the power
    !     of the tail of the N(kx,ky) spectrum, we have q=p-1
    !     Since N(k) = k N(kx,ky) with k the Jacobian
    !     it follows that the tail functions are given by
    !     k^p = k k^q   =>   k^p = k^(q+1) => p=q+1  => q=p-1
    !  8. Structure
    !     Initialisations
    !     do for all points on locus
    !       compute directional index for k2 and k4
    !       if geometric scaling then
    !         compute wave number index directly
    !         convert log-scaling to linear scaling
    !       else
    !         search position of wave number in k-array
    !         if k < kmin then
    !           k-index = 1 and factor is 0
    !         elsif k < kmax then
    !           compute k-index for k2 and k4
    !         else
    !           compute tail factor
    !         end if
    !       end if
    !     end do
    !  9. Switches
    !     /T  enable test output
    ! 10. Source code:
    !------------------------------------------------------------------------------
    !     Local variables
    integer iloc      ! counter along locus loop
    integer jpos      ! index for interpolation and tracking of position in wave numebr array
    real k2a,k2m      ! angle (radians) and magnitude of wave number k2
    real k4a,k4m      ! angle (radians) and magnitude of wave number k2
    real dk           ! difference between two wave numbers
    real xtest        ! test value for checking computation of weights, by inversion test
    real ff,gg        ! variables in transformation of log-spacing to linear spacing
    type(gridvar),pointer:: pg
    ! initialisations
    do iloc=1,nlocus
      pg=>gv(iloc)
      k2a = pg%k2a_mod
      k4a = pg%k4a_mod
      ! compute directional weights
      pg%wa_k2 = (k2a)/q_deltad  !+1 !ZZM index from 0
      pg%wa_k4 = (k4a)/q_deltad  !+1 !ZZM index from 0

      k2m = pg%k2m_mod
      k4m = pg%k4m_mod
      !  compute position of k2 in wave number grid
      !  and compute weight function
      ! deep water is assumed and loci have geometric scaling
      pg%wk_k2 = alog(k2m/kqmin)/q_lkfac !+1 !ZZM index from 0
      pg%wk_k4 = alog(k4m/kqmin)/q_lkfac !+1 !ZZM index from 0
      if(pg%wk_k2<0)then
        pg%wr_k2 = k2m/kqmin
      else
        pg%wr_k2 = q_kfac**(-3.5*(pg%wk_k2-nkq))
      endif
      if(pg%wk_k4<0)then
        pg%wr_k4 = k4m/kqmin
      else
        pg%wr_k4 = q_kfac**(-3.5*(pg%wk_k4-nkq))
      endif
      !  Replace log-spacing by linear spacing
      ff = pg%wk_k2; gg = floor(ff)
      pg%wk_k2 = gg+(q_kfac**(ff-gg)-1.)/(q_kfac-1.)
      ff = pg%wk_k4; gg = floor(ff)
      pg%wk_k4 = gg+(q_kfac**(ff-gg)-1.)/(q_kfac-1.)
    end do
  end subroutine
  !------------------------------------------------------------------------------
  subroutine q_polar2(kmin,kmax,kx_beg,ky_beg,kx_end,ky_end,loclen,ierr)
    !use m_mkgrid,only:k_pol,a_pol,c_pol,y2_loc,x2_loc,
    !use m_mkgrid,only:x_cosk,dk0,mlocus,nlocus0,nlocus1,pang
    implicit none
    !  1. Purpose:
    !     Compute position of locus for given k1-k3 vector
    !  2. Method
    !     Explicit polar method, see Van Vledder 2000, Monterey paper
    !     Optionally using a fixed k-step, geometric k-step or adaptive stepping
    !  3. Parameters used:
    !Type    I/O        Name             Description
    !------------------------------------------------------------------------------
    real, intent(in) :: kmin           ! minimum wave number on locus
    real, intent(in) :: kmax           ! maximum wave number on locus
    real, intent(in) :: kx_beg         ! x-coordinate of begin point
    real, intent(in) :: ky_beg         ! y-coordinate of begin point
    real, intent(in) :: kx_end         ! x-coordinate of end point
    real, intent(in) :: ky_end         ! y-coordinate of end point
    real, intent(in) :: loclen         ! estimated length of locus
    integer, intent (out)  :: ierr     ! error condition
    !     Parameters with module
    !     nlocus0   Preferred number of points on locus
    !     q         w1-w3, difference of radian frequencies
    !     pmag      |k1-k3| (vector form)
    !     pdir      direction of difference vector k1-k3
    !  4. Error messages
    !  5. Called by:
    !     Q_CPMLOCUS
    !  6. Subroutines used:
    !     X_COSK
    !  7. Remarks
    !     The type of locus computation is controlled by the parameter IQ_LOCUS
    !     Set in Q_SETCFG
    !  8. Structure
    !  9. Switches
    !     /S  enable subroutine tracing
    !     /T  enable test output
    ! 10. Source code
    !------------------------------------------------------------------------------
    !     Local variabels
    integer ipol      ! counter
    integer jpol      ! counter
    integer ipass     ! counter for passes
    integer npol      ! number of points on locus
    integer npass     ! number of passes in iteration process
    integer mpol      ! maximum number of points on locus, related to MLOCUS
    real kold         ! temporary wave number
    real knew         ! temporary wave number
    real cosold       ! 'old' cosine of angle
    real cosnew       ! 'new' cosine of angle
    real dkpol        ! step in wave number
    real dkold        ! 'old' step in wave number
    real ang1         ! 'old' angle
    real ang2         ! 'new' angle
    real kratio       ! ratio between succesive k-values when IQ_LOCUS=3
    real arg          ! argument
    real dk           ! step in wave number
    real dsnew        ! new step size along locus
    real dsz          ! estimated step size along locus
    real kpol,cpol,apol
    !------------------------------------------------------------------------------
    ! initialisations
    !------------------------------------------------------------------------------
    ierr = 0                      ! set error code to zero
    npol = (nlocus0+1)/2+1        ! first estimate of number k-values along symmetry axis
    mpol = mlocus/2               ! set maximum number of points along locus axis
    select case(iq_locus)
    case(1)
      ! CASE = 1: Linear spacing of wave numbers along symmetry axis
      dk = (kmax-kmin)/real(npol-1)
      do ipol=1,npol
        jpol = 2*npol-ipol
        kpol = kmin + (ipol-1)*dk
        cpol = x_cosk(kpol)
        apol=  acos(cpol)
        x2_loc(ipol) = kpol*cos(pang + apol)
        y2_loc(ipol) = kpol*sin(pang + apol)
        x2_loc(jpol) = kpol*cos(pang - apol)
        y2_loc(jpol) = kpol*sin(pang - apol)
      end do
    case(2)
      !  Case = 2: Variable k-stepping along symmetry axis,
      !            such that step along locus is more or less constant
      !-------------------------------------------------------------------------------
      ! set first point on locus
      ! compute initial step size of polar wave number
      dsz  = loclen/real(nlocus0)          ! estimate of step size along locus
      npass= 3                             ! set number of passes in iteration
      dk   = (kmax - kmin)/real(npol)/2    ! estimate of step size of equidistant radii

      ipol = 1;jpol=mpol*2-ipol;

      kpol = kmin;ang1=pi;kold = kmin;

      x2_loc(ipol) = kpol*cos(pang + ang1)
      y2_loc(ipol) = kpol*sin(pang + ang1)

      do while (kpol < kmax .and.  ipol < mpol)
        do ipass=1,npass
          knew  = min(kmax,kpol+dk)
          dkold = knew - kpol
          cosnew = x_cosk(knew)
          ang2  =  acos(cosnew)
          arg   = kold**2 + knew**2 -2.*kold*knew*cos(ang2-ang1)
          dsnew = sqrt(abs(arg))
          if(dsnew>0) dk   = dk*dsz/dsnew
        end do
        kold  = knew; ang1=ang2
        !  assign new estimate and check value of IPOL
        ipol = ipol + 1;jpol=jpol-1
        kpol = kpol + dkold
        x2_loc(ipol) = kpol*cos(pang + ang1)
        y2_loc(ipol) = kpol*sin(pang + ang1)
        x2_loc(jpol) = kpol*cos(pang - ang1)
        y2_loc(jpol) = kpol*sin(pang - ang1)
        if (abs(dkold) < 0.0005*(kmax-kmin)) exit
      end do
      ! fill last bin with coordinates of end point
      if(kpol < kmax .and. ipol <  mpol) then
        ipol = ipol + 1;jpol=jpol-1
        kpol = kmax; ang1=acos(-1.)
        x2_loc(ipol) = kpol*cos(pang + ang1)
        y2_loc(ipol) = kpol*sin(pang + ang1)
        x2_loc(jpol) = kpol*cos(pang - ang1)
        y2_loc(jpol) = kpol*sin(pang - ang1)
      end if
      !  update the number of k-points on symmetry axis
      npol = ipol
      !------------------------------------------------------------------------------
      !  compute actual number of points on locus
      !  this will always be an even number
      !  mirror image the second half of the locus
      nlocus1 = 2*npol-2
      do ipol=npol+1,nlocus1 
        x2_loc(ipol) = x2_loc(jpol)
        y2_loc(ipol) = y2_loc(jpol)
        jpol = jpol+1
      end do
    case(3)
      !  Case 3: Geometric spacing of wave numbers along symmetry axis
      kratio = (kmax/kmin)**(1./(npol-1.))
      do ipol=1,npol
        jpol = 2*npol-ipol
        kpol = kmin*kratio**(ipol-1.)
        cpol = x_cosk(kpol)
        apol=  acos(cpol)
        x2_loc(ipol) = kpol*cos(pang + apol)
        y2_loc(ipol) = kpol*sin(pang + apol)
        x2_loc(jpol) = kpol*cos(pang - apol)
        y2_loc(jpol) = kpol*sin(pang - apol)
      end do
    end select
  end subroutine
  real function x_jacobian(x2,y2,x4,y4)
    !use m_mkgrid,only:q_depth,q_grav
    implicit none
    !  1. Purpose:
    !     Compute gradient/Jacobian term for a given point on the locus
    !  2. Method
    !     Explicit expressions for gradient term
    !     Using expression of Rasmussen (1998)
    !     J = |cg2-cg4|
    !  3. Parameter list:
    ! Type    I/O        Name    Description
    !--------------------------------------------------------------------
    real, intent(in) ::  x2   !  x-component of wave number k2
    real, intent(in) ::  y2   !  y-component of wave number k2
    real, intent(in) ::  x4   !  x-component of wave number k4
    real, intent(in) ::  y4   !  y-component of wave number k4
    !  4. Error messages
    !  5. Called by:
    !     Q_CMPLOCUS
    !  6. Subroutines used:
    !  7. Remarks
    !  8. Structure
    !  9. Switches
    ! 10. Source code:
    !------------------------------------------------------------------------------
    ! local variables
    real k2m,k4m    ! wave number magnitudes
    real k2md,k4md  ! k*d values
    real ang2,ang4  ! directions
    real cg2,cg4    ! group velocities
    real sig2,sig4  ! radian frequencies
    real cg2x,cg2y  ! components of group velocity cg2
    real cg4x,cg4y  ! components of group velocity cg4
    !------------------------------------------------------------------------------
    k2m = sqrt(x2**2 + y2**2)
    k4m = sqrt(x4**2 + y4**2)
    ang2 = atan2(x2,y2)
    ang4 = atan2(x4,y4)
    sig2 = sqrt(q_grav*k2m*tanh(k2m*q_depth))
    sig4 = sqrt(q_grav*k4m*tanh(k4m*q_depth))
    k2md = k2m*q_depth
    k4md = k4m*q_depth
    if(k2md > 20) then
      cg2 = 0.5*q_grav/sig2
    else
      cg2 = sig2/k2m*(0.5+k2md/sinh(2*k2md))
    end if
    if(k4md > 20) then
      cg4 = 0.5*q_grav/sig4
    else
      cg4 = sig4/k4m*(0.5+k4md/sinh(2*k4md))
    end if
    ! x_jacobian = sqrt(cg2**2+cg4**2-2*cg2*cg4*cos(ang2-ang4))
    !
    ! Modified 15/12/2011, based on problem reported by Ruslan Puscasu
    !
    cg2x = cg2*cos(ang2)
    cg2y = cg2*sin(ang2)
    cg4x = cg4*cos(ang4)
    cg4y = cg4*sin(ang4)
    x_jacobian = sqrt((cg2x-cg4x)**2 + (cg2y-cg4y)**2)
    return
  end function
  !------------------Not Use data -----------------------------------------------
  !-----------------------------------------------------------------------------!
  subroutine z_polyarea(xpol,ypol,npol,area)
    implicit none
    !  1. Purpose
    !     Computes area of a closed polygon
    !  2. Method
    !     The area of the polygon
    !  3. Parameter list
    !     Name    I/O  Type  Description
    integer, intent(in)  ::  npol       ! Number of points of polygon
    real, intent(in)     ::  xpol(npol) ! x-coodinates of polygon
    real, intent(in)     ::  ypol(npol) ! y-coordinates of polygon
    real, intent(out)    ::  area       ! area of polygon
    !  4. Subroutines used
    !  5. Error messages
    !  6. Remarks
    integer ipol,ipol1         ! counters
    real xmin,xmax,ymin,ymax   ! minima and maxima of polygon
    real xmean,ymean           ! mean values
    real xa,ya,xb,yb           ! temporary variables
    real darea                 ! piece of area
    !-------------------------------------------------------------------------------
    if(npol<=1) then
      area = 0.
      return
    end if
    ! compute minimum and maximum coordinates
    xmin = minval(xpol)
    xmax = maxval(xpol)
    ymin = minval(ypol)
    ymax = maxval(ypol)
    !  compute mean of range of x- and y-coordinates
    xmean = 0.5*(xmin + xmax)
    ymean = 0.5*(ymin + ymax)
    ! compute area and center of gravity
    ! do loop over all line pieces of polygon
    area = 0.
    do ipol=1,npol
      ipol1 = ipol + 1
      if(ipol==npol) ipol1 = 1
      xa = xpol(ipol)
      ya = ypol(ipol)
      xb = xpol(ipol1)
      yb = ypol(ipol1)
      darea = 0.5*((xa-xmean)*(yb-ymean) - (xb-xmean)*(ya-ymean))
      area  = area + darea
    end do
  end subroutine
  !-----------------------------------------------------------------
  subroutine q_loc_w1w3(k1x,k1y,k3x,k3y,npts,k2x,k2y,k4x,k4y,s)
    implicit none
    !  1. Purpose:
    !     Compute locus for the special case w1=w3
    !  2. Method
    !     For this case, the k2-locus consists of a straight line
    !  3. Parameter used:
    integer, intent(in) :: npts      ! Number of points
    real, intent(in)    :: k1x       ! x-component of wave number k1
    real, intent(in)    :: k1y       ! y-component of wave number k1
    real, intent(in)    :: k3x       ! x-component of wave number k3
    real, intent(in)    :: k3y       ! y-component of wave number k3
    real, intent(out)   :: k2x(npts) ! x-component of wave number k2
    real, intent(out)   :: k2y(npts) ! y-component of wave number k2
    real, intent(out)   :: k4x(npts) ! x-component of wave number k4
    real, intent(out)   :: k4y(npts) ! y-component of wave number k4
    real, intent(out)   :: s(npts)   ! distance along locus
    !  4. Error messages
    !  5. Caled by:
    !     Q_CMPLOCUS
    !  6. Subroutines used
    !  7. Remarks
    !     Routine based on modified version of routine SHLOCX of Resio and Tracy
    !     On 15/4/2002 a bug fixed in computation of THR when angle of k3 is larger than 90°
    !     In addition, the assumption that k1y=0 and thus dir1=0 is removed
    !     In bug fix of 20/8/2002 this restriction is removed.
    !  8. Structure
    !     Compute angle of symmetry axis
    !     Compute distance between 2 lines of solution
    !     compute wave numbers along locus
    !     rotate angles
    !  9. Switches
    ! 10. Source code
    !------------------------------------------------------------------------------
    !     Local variables
    integer ipt   ! counter of points along locus
    real dirs     ! angle of symmetry axis
    real dir1     ! direction of wave number k1
    real dir3     ! direction of wave number k3
    real dk0      ! step size along locus
    real xk0      ! x-component
    real yk0      ! y-component
    real w2       ! radian frequency
    real xx2      ! values along k2-locus
    real k1m      ! magnitude of wave number k1
    real cosd,sind,xx,xy,yx,yy,dxk0,xoff
    integer ipt2
    !------------------------------------------------------------------------------
    !    dirs is the angle of rotation from the x-axis to the "bisecting" angle
    dir1 = atan2(k1y,k1x)
    dir3 = atan2(k3y,k3x)
    dirs = 0.5*(pi-abs(pi-abs(dir3-dir1))) !fixme
    k1m  = sqrt(k1x**2 + k1y**2)
    !  k1x is the total length of the wavenumber vector
    !  xk0 is the length of this vector in the rotated coordinate system
    xk0 = k1m * cos(dirs)
    yk0 = k1m * sin(dirs)
    ! Specify step size for solution of singular case
    dk0 = kqmax/real(npts-1.)       !this value depends on actual grid
    !dk0 = 3.01/real(npts-1.)       !  this is test value
    !  modify rotation angle
    dirs = dirs + dir1
    !  generate sequence of parallel lines rotate lines over modified angle DIRS
    cosd=cos(dirs); sind=sin(dirs);
    yy=yk0*sind; yx=yk0*cosd 
    dxk0=dk0*xk0;
    ipt2=(-npts/2)*2 ; !fixme for compat ,should be ipt2=-npts
    do ipt=1,npts
      s(ipt)   = (ipt-1)*dxk0
      ipt2     = ipt2+2  ; xx2      = ipt2*dxk0
      xx       = xx2*cosd; xy       = xx2*sind ;
      k2x(ipt) = xx - yy ; k2y(ipt) = xy + yx 
      k4x(ipt) = xx + yy ; k4y(ipt) = xy - yx 
      !|k2|==|k4|
    end do
  end subroutine
  !------------------------------------------------------------------------------
  real function x_cple(k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,depth,grav)
    implicit none
    !  1. Purpose:
    !     Compute coupling coefficient between a quadruplet of
    !     interacting wave numbers
    !  2. Method
    !  3. Parameter list:
    !     Name    I/O  Type  Description
    !Type    I/O              Name      Description
    !-----------------------------------------------------------------------------
    real, intent(in) ::       k1x     !  x-component of wave number k1
    real, intent(in) ::       k1y     !  y-component of wave number k1
    real, intent(in) ::       k2x     !  x-component of wave number k2
    real, intent(in) ::       k2y     !  y-component of wave number k2
    real, intent(in) ::       k3x     !  x-component of wave number k3
    real, intent(in) ::       k3y     !  y-component of wave number k3
    real, intent(in) ::       k4x     !  x-component of wave number k4
    real, intent(in) ::       k4y     !  y-component of wave number k4
    real, intent(in) ::       depth   !  Water depth in meters
    real, intent(in) ::       grav    !  Gravitational acceleration
    !  5. Called by:
    !     Q_CMPLOCUS
    !  6. Subroutines used
    !     X_WEBB
    !     X_HH
    !------------------------------------------------------------------------------
    !  1) Deep water coupling coefficient of Webb
    x_cple = xc_webb(k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,grav)
  end function
  !------------------------------------------------------------------------------
  real function xc_webb(k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,grav)
    !  1. Purpose:
    !     Compute deep water coupling coefficient for
    !     non-linear quadruplet interactions
    !  2. Method
    !     Webb (1978) and modified and corrected by Dungey and Hui (1979)
    !  3. Parameter list:
    ! Type    I/O        Name    Description
    real, intent(in) ::  k1x   !  x-component of wave number k1
    real, intent(in) ::  k1y   !  y-component of wave number k1
    real, intent(in) ::  k2x   !  x-component of wave number k2
    real, intent(in) ::  k2y   !  y-component of wave number k2
    real, intent(in) ::  k3x   !  x-component of wave number k3
    real, intent(in) ::  k3y   !  y-component of wave number k3
    real, intent(in) ::  k4x   !  x-component of wave number k4
    real, intent(in) ::  k4y   !  y-component of wave number k4
    real, intent(in) ::  grav  !  gravitational acceleration m/s^2
    !  4. Error messages
    !  5. Called by:
    !     X_CPLE
    ! 10. Source code:
    ! local variables
    double precision wsqp12         ! derived variable
    double precision wsqm13         ! derived variable
    double precision wsq13          ! derived variable
    double precision wsqm14         ! derived variable
    double precision wsq14          ! derived variable
    double precision wsq12          ! derived variable
    real z,z12,z13,z14              ! derived variables
    real dwebb                      ! final coefficient
    real p1,p2,p3,p4,p5,p6,p7,p8,p9 ! partial summations
    real w1,w2,w3,w4                ! radian frequencies
    real k1,k2,k3,k4                ! wave number magnitudes
    real dot12                      ! k1*k2
    real dot13                      ! k1*k3
    real dot14                      ! k1*k4
    real dot23                      ! k2*k3
    real dot24                      ! k2*k4
    real dot34                      ! k3*k4
    real eps                        ! internal accuracy
    ! initialisations
    eps = 1.0e-30
    k1 = sqrt(k1x*k1x + k1y*k1y);w1 = sqrt(k1)
    k2 = sqrt(k2x*k2x + k2y*k2y);w2 = sqrt(k2)
    k3 = sqrt(k3x*k3x + k3y*k3y);w3 = sqrt(k3)
    k4 = sqrt(k4x*k4x + k4y*k4y);w4 = sqrt(k4)

    dot12 = k1x*k2x + k1y*k2y; dot23 = k2x*k3x + k2y*k3y
    dot13 = k1x*k3x + k1y*k3y; dot24 = k2x*k4x + k2y*k4y
    dot14 = k1x*k4x + k1y*k4y; dot34 = k3x*k4x + k3y*k4y
    
    wsqp12= sqrt((k1x+k2x)*(k1x+k2x)+(k1y+k2y)*(k1y+k2y))
    wsqm13= sqrt((k1x-k3x)*(k1x-k3x)+(k1y-k3y)*(k1y-k3y))
    wsqm14= sqrt((k1x-k4x)*(k1x-k4x)+(k1y-k4y)*(k1y-k4y))
    wsq12 = (w1+w2)*(w1+w2)
    wsq13 = (w1-w3)*(w1-w3)
    wsq14 = (w1-w4)*(w1-w4)
    z12   = wsqp12-wsq12
    z13   = wsqm13-wsq13
    z14   = wsqm14-wsq14
    p1    = 2.  *wsq12*(k1*k2-dot12)*(k3*k4-dot34)/(z12+eps)
    p2    = 2.  *wsq13*(k1*k3+dot13)*(k2*k4+dot24)/(z13+eps)
    p3    = 2.  *wsq14*(k1*k4+dot14)*(k2*k3+dot23)/(z14+eps)
    p4    = 0.5 *(dot12*dot34 + dot13*dot24 + dot14*dot23)
    p6    =-0.25*wsq12**2 *(dot12+dot34) 
    p5    = 0.25*wsq13**2 *(dot13+dot24) 
    p7    = 0.25*wsq14**2 *(dot14+dot23) 
    p8    = 2.5*k1*k2*k3*k4
    p9    = wsq12*wsq13*wsq14* (k1 + k2 + k3 + k4)
    dwebb  = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9
    xc_webb = grav**2*pi*0.25*dwebb*dwebb/(w1*w2*w3*w4+eps)
    return
  end  function
  !-----------------------------------------------------------------------------!
  subroutine z_intp1(x1,y1,x2,y2,n1,n2,ierr)                                    !
    implicit none
    !  1. Purpose
    !     Interpolate function values
    !  2. Method
    !     Linear interpolation

    !     If a requested point falls outside the input domain, then
    !     the nearest point is used (viz. begin or end point of x1/y1 array
    !     If the input array has only one point. A constant value is assumed
    !  3. Parameter list
    !     Name    I/O  Type  Description
    integer, intent(in) ::  n1   !   number of data points in x1-y1 arrays
    integer, intent(in) ::  n2   !   number of data points in x2-y2 arrays
    real, intent(in) ::  x1(n1)  !   x-values of input data
    real, intent(in) ::  y1(n1)  !   y-values of input data
    real, intent(in) ::  x2(n2)  !   x-values of output data
    real, intent(out) :: y2(n2)  !   y-values of output data
    integer, intent(out) :: ierr !   Error indicator
    !  4. Subroutines used
    !  5. Error messages
    !     ierr = 0    No errors detected
    !          = 1    x1-data not monotonic increasing
    !          = 10   x2-data not monotonic increasing
    !          = 11   x1- and x2 data not monotonic increasing
    !          = 2    x1-data not monotonic decreasing
    !          = 20   x1-data not monotonic decreasing
    !          = 22   x1- and x2 data not monotonic decreasing
    !          = 2    No variation in x1-data
    !          = 3    No variation in x2-data is allowed if n2=1
    !  6. Remarks
    !     It is assumed that the x1- and x2-data are either
    !     monotonic increasing or decreasing
    !     If a requested x2-value falls outside the range of x1-values
    !     it is assumed that the corresponding y2-value is equal to
    !     the nearest boundary value of the y1-values
    !     Example: x1 = [0 1 2 3]
    !              y1 = [1 2 1 0]
    !              x2 = -1,  y2 = 1
    !              x2 =  5,  y2 = 0
    !------------------------------------------------------------------------------
    integer i1,i2        ! counters
    real ds            ! step size
    real fac           ! factor in linear interpolation
    real s1,s2         ! search values
    real xmin1,xmax1   ! minimum and maximum of x1-data
    real xmin2,xmax2   ! minimum and maximum of x2-data
    real, parameter :: eps=1.e-20
    !------------------------------------------------------------------------------
    !   initialisation
    ierr = 0
    !  check number of points of input array
    if(n1==1) then
      y2 = y1(1)
      return
    end if
    !  check minimum and maximum data values
    xmin1 = minval(x1)
    xmax1 = maxval(x1)
    xmin2 = minval(x2)
    xmax2 = maxval(x2)
    if (abs(xmin1-xmax1) < eps .or. abs(x1(1)-x1(n1)) < eps) then
      ierr = 2
      return
    end if
    if ((abs(xmin2-xmax2) < eps .or. abs(x2(1)-x2(n2)) < eps) .and. n2 > 1) then
      ierr = 3
      return
    end if
    ! check input data for monotonicity
    if(x1(1) < x1(n1)) then             ! data increasing
      do i1=1,n1-1
        if(x1(i1) > x1(i1+1)) then
          ierr=1
          write(*,*) 'z_intp1: i1 x1(i1) x1(i1+1):',i1,x1(i1),x1(i1+1)
          return
        end if
      end do
      do i2=1,n2-1
        if(x2(i2) > x2(i2+1)) then
          ierr=ierr+10
          write(*,*) 'z_intp1: i2 x2(i2) x2(i2+1):',i2,x2(i2),x2(i2+1)
          return
        end if
      end do
    else                                 ! data decreasing
      do i1=1,n1-1
        if(x1(i1) < x1(i1+1)) then
          ierr=2
          write(*,*) 'z_intp1: i1 x1(i1) x1(i1+1):',i1,x1(i1),x1(i1+1)
          return
        end if
      end do
      do i2=1,n2-1
        if(x2(i2) < x2(i2+1)) then
          ierr=ierr + 20
          write(*,*) 'z_intp1: i2 x2(i2) x2(i2+1):',i2,x2(i2),x2(i2+1)
          return
        end if
      end do
    end if
    !------------------------------------------------------------------------------
    ! initialize
    !------------------------------------------------------------------------------
    if(ierr==0) then
      i1 = 1
      s1 = x1(i1)
      do i2 = 1,n2
        s2 = x2(i2)
        do while (s1 <= s2 .and. i1 < n1)
          i1 = i1 + 1
          s1 = x1(i1)
        end do
        !  special point
        !  choose lowest s1-value if x2(:) < x1(1)
        if(i1 ==1) then
          y2(i2) = y1(i1)
        else
          ds = s2 - x1(i1-1)
          fac = ds/(x1(i1)-x1(i1-1))
          y2(i2) = y1(i1-1) + fac*(y1(i1)-y1(i1-1))
        end if
        ! special case at end: choose s2(n2) > s1(n1), choose last value of y1(1)
        if(i2==n2 .and. s2>s1) y2(n2) = y1(n1)
      end do
    end if
  end subroutine
  !------------------------------------------------------------------------------
  ! set symfac to 1|0 
  subroutine q_symmetry(k1x,k1y,k3x,k3y,k4x,k4y,symfac,nloc)
    implicit none
    !  1. Purpose:
    !     Compute symmetry factor to reduce integration
    !  2. Method
    !     Compute distance between k1 and k3, and between k4 and k1
    !  3. Parameter list:
    ! Type   i/o             Name           Description
    !----------------------------------------------------------------------------------
    integer, intent(in)   :: nloc         ! number of points in array with wave number
    real, intent(in)      :: k1x          ! x-component  of wave number k1
    real, intent(in)      :: k1y          ! y-component  of wave number k1
    real, intent(in)      :: k3x          ! x-component  of wave number k3
    real, intent(in)      :: k3y          ! y-component  of wave number k3
    real, intent(in)      :: k4x(nloc)    ! x-components of wave number k4
    real, intent(in)      :: k4y(nloc)    ! y-components of wave number k4
    real, intent(out)     :: symfac(nloc) ! symmetry factor
    !----------------------------------------------------------------------------------
    !  5. Called by:
    !     Q_MODIFY
    integer iloc      ! counter
    real dk13         ! distance between k1 and k3
    real dk14         ! distance between k1 and k4
    ! evaluate criterion |k3-k1| < |k4-k1|
    ! if false then symfac=0
    if(iq_syma==0) return
    dk13 = (k1x-k3x)**2 + (k1y-k3y)**2
    do iloc=1,nloc
      dk14 = (k1x-k4x(iloc))**2 + (k1y-k4y(iloc))**2
      if (dk13 >= dk14) symfac(iloc) = 0.
    end do
  end subroutine
  !-----------------------------------------------------------------------------!
  real function z_root2(func,x1,x2,xacc,ierr)
    implicit none
    !  1. Purpose
    !     Find zero crossing point of function FUNC between the
    !     initial values on either side of zero crossing
    !  2. Method
    !     Ridders method of root finding
    !     adapted from routine zridddr
    !     Numerical Recipes
    !     The art if scientific computing, second edition, 1992
    !     W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery
    !  3. Parameter list
    !     Name    I/O  Type  Description
    !     func     i    r    real function
    !     x1       i    r    initial x-value on left/right side of zero-crossing
    !     x2       i    r    initial x-value on right/left side of zero-crossing
    !     xacc     i    r    accuracy, used as |x1(i)-x2(i)|< xacc
    !     ierr     o    i    Error indicator
    !  4. Subroutines used
    !     Func      user supplied real function
    !  5. Error messages
    !     ierr = 0   No errors occured during iteration process
    !            1   Iteration halted in dead end, this combination may NEVER occur
    !            2   Maximum number of iterations exceeded
    !            3   Solution jumped outside interval
    !  6. Remarks
    !     It is assumed that the x1- and x2-coordinate lie
    !     on different sides of the actual zero crossing
    real func                         ! external function
    real, intent (in) :: x1           ! x-value at one side of interval
    real, intent (in) :: x2           ! x-value at other side of interval
    real, intent (in) :: xacc         ! requested accuracy
    integer, intent (out) :: ierr     ! error indicator
    real unused                       ! default value
    real zriddr                       ! intermediate function value
    real xx1,xx2,xx                   ! local boundaries during iteration
    integer maxit                     ! maximum number of iteration
    logical lopen                     ! check if a file is opened

    parameter (maxit = 20)
    external func
    integer iter      ! counter for number of iterations
    real fh           ! function value FUNC(xh)
    real fl           ! function value FUNC(xl)
    real fm           ! function value FUNC(xm)
    real fnew         ! function value FUNC(xnew)
    real s            ! temp. function value, used for inverse quadratic interpolation
    real xh           ! upper (high) boundary of interval
    real xl           ! lower boundary of interval
    real xm           ! middle point of interval
    real xnew         ! new estimate according to Ridders method
    ierr   = 0        ! set error level
    unused =-1.11e30  ! set start value
    xx1 = x1          ! copy boundaries of interval to local variables
    xx2 = x2
    ! check boundaries on requirement x2 > x1
    if(xx1 > xx2) then
      xx  = xx1
      xx1 = xx2
      xx2 = xx
    end if
    fl = func(xx1)
    fh = func(xx2)
    if(fl*fh <= 0. )then
      xl = xx1
      xh = xx2
      zriddr = unused
      do iter=1,maxit
        xm = 0.5*(xl+xh)
        fm = func(xm)
        s = sqrt(fm**2-fl*fh)
        if(s == 0.) exit
        xnew = xm+(xm-xl)*(sign(1.,fl-fh)*fm/s)
        if (abs(xnew-zriddr) <= xacc) exit
        zriddr = xnew
        fnew = func(zriddr)
        if (fnew == 0.)exit
        if(sign(fm,fnew) /= fm) then
          xl = xm; fl = fm; xh = zriddr; fh = fnew
        elseif(sign(fl,fnew) /= fl) then
          xh = zriddr; fh = fnew
        elseif(sign(fh,fnew) /= fh) then
          xl = zriddr; fl = fnew
        else
          ierr = 1;exit
        endif
        if(abs(xh-xl) <= xacc)exit
      end do
      if(iter>maxit)ierr = 2
      z_root2 = zriddr
    else if (fl == 0.) then
      z_root2 = xx1
    else if (fh == 0.) then
      z_root2 = xx2
    else
      ierr = 3
    endif
    return
  end function
  real function tanz(x)
    real x
    if (x.gt.20.) x=25.
    tanz=tanh(x)
    return
  end  function
  real function cosz(x)
    real x
    if (x.gt.20.) x=25.
    cosz=cosh(x)
    return
  end  function
  !=================================================
  subroutine setxnltype(iq_type_,iq_geom_,iq_grid_,iq_locus_,iq_search_,iq_sym_,iq_compact_,iq_mod_)
    integer::iq_type_,iq_geom_,iq_grid_,iq_locus_,iq_search_,iq_sym_,iq_compact_,iq_mod_
    if(iq_type_/=2)then !  deep water computation with WAM depth scaling based on Herterich and Hasselmann (1980)
      print*,"Error wrtvv:Problem in setxntype"
      stop
    endif
    if(iq_grid_   >0)iq_grid   =iq_grid_
    if(iq_locus_  >0)iq_locus  =iq_locus_
    if(iq_search_ >0)iq_search =iq_search_
    if(iq_sym_    >0)iq_sym    =iq_sym_
    if(iq_compact_>0)iq_compact=iq_compact_
    if(iq_mod_    >0)iq_mod    =iq_mod_
  end subroutine
  subroutine dbghere
    print*,'dbghere'
  end subroutine 
  subroutine gxnldata_alloc
    mkq = 1+qf_kn       !nkq qf_kn=alog(qf_krat)/alog(q_kfac),qf_krat=2.5,6,q_kfac=1.21
    maq = 1+qf_dn+qf_dn !naq qf_dn=qf_dmax/q_deltad,naq=360/q_deltad,qf_dmax=60,75 <180
    klocus = nlocus0
    VALLOC(q_km  ,(nkq+1))
    VALLOC(q_ca  ,(naq))
    VALLOC(q_sa  ,(naq))
    VALLOC(q_ad  ,(naq))
    VALLOC(q_dk  ,(nkq))
    VALLOC(q_k2  ,(nkq))
    VALLOC(q_ka  ,(nkq))
    VALLOC(q_ar  ,(nkq))
    VALLOC(q_cg  ,(nkq))
    VALLOC(q_sig ,(nkq))
    VALLOC(q_sigr,(nkq))
    VALLOC(q_dsig,(nkq))
    VALLOC(q_kpow,(nkq))
    VALLOC(q_lamd,(nkq))
    VALLOC(q_tail,(nkq))
    VALLOC(quad_nloc,(mkq,maq))
    VALLOC(quads,(klocus,mkq,maq))
  end subroutine gxnldata_alloc
  !------------------------------------------------------------------------------
  subroutine q_chkcons(xnl,sum_e,sum_a,sum_mx,sum_my)
    implicit none
    !  1. Purpose:
    !     Check conservation laws of non-linear transfer
    !  2. Method
    !     The following conservation laws should be fulfilled:
    !     Wave Energy        SUME=0
    !     Wave Action        SUMA=0
    !     Momentum vector    SUMMX,SUMMY=0
    real, intent(in)  :: xnl(nkq,naq) ! transfer rate
    real, intent(out) :: sum_e        ! sum of wave energy
    real, intent(out) :: sum_a        ! sum of wave action
    real, intent(out) :: sum_mx       ! sum of momentum in x-direction
    real, intent(out) :: sum_my       ! sum of momentum in y-direction
    real aa    ! action density
    real ee    ! energy density
    real kk    ! wave number
    real momx  ! momentum in x-direction
    real momy  ! momentum in y-direction
    real qq,qqk! bin size
    integer ia ! counter over directions
    integer ik ! counter over wave numbers
    sum_a  = 0.; sum_e  = 0.; sum_mx = 0.; sum_my = 0.
    do ik=1,nkq
      qq = q_ar(ik)
      qqk= q_ka(ik)
      kk = q_k(ik)
      do ia = 1,naq
        aa = xnl(ik,ia)
        ee = aa*q_sig(ik)
        momx = aa*q_ca(ia)
        momy = aa*q_sa(ia)
        sum_a  = sum_a  + aa*qq
        sum_e  = sum_e  + ee*qq
        sum_mx = sum_mx + momx*qqk
        sum_my = sum_my + momy*qqk
      end do
    end do
  end subroutine
end module m_mkgrid
