#define VALLOC(N,S) if(allocated(N)) deallocate(N) ; allocate (N S)
#define GETLOCN
module m_xnldata
  use m_mkgrid
  implicit none
  public xnl_init,xnlcal,q_dscale
  private
  ! set by q_getlocus ,used by q_t13v4
  type(wws_type),pointer :: t_q(:)! transformed quads   
  !real, allocatable :: aspec(:,:)     !  Action density on wave number grid A(sigma,theta)
  real, allocatable :: nspec(:,:) !  Action density on wave number grid N(kx,ky)
  real q_scale                    ! additional scale factor resulting from SEARCH for neasrest grid
contains
  !----------------------------------------------------------------------------------
  subroutine xnl_init(k,dirr,nk,ndir,pftail,iquad,iproc,ierror)
    !  1. Purpose:
    !     Initialize coefficients, integration space, file i/o for computation
    !     nonlinear quadruplet wave-wave interaction
    !  2. Method
    !     Set version number
    !     Set unit unit numbers
    !     Open quad related files
    !     Optionally reset configuration by a back door option
    !     Compute integration spaces for given water depths
    !  3. Parameter list:
    !Type   I/O              Name          Description
    !------------------------------------------------------------------------------
    integer, intent(in)  ::  nk            ! Number of sigma values
    integer, intent(in)  ::  ndir          ! Number of directions
    real(8), intent(in)  ::  k(nk)         ! Radian frequencies
    real(8), intent(in)  ::  dirr(ndir)    ! Directions (degrees)
    real(8), intent(in)  ::  pftail        ! power of spectral tail, e.g. -4 or -5
    integer, intent(in)  ::  iquad         ! Type of method for computing nonlinear interactions
    integer, intent(in)  ::  iproc         ! Processor number, controls output file for MPI
    integer, intent(out) ::  ierror        ! Error indicator. If no errors are detected IERR=0
    !  5. Called by:
    !     host program, e.g. SWANQUAD4
    !  6. Subroutines used:
    !     Q_INIT
    real(8) qkfach
    integer ikq
    ! Initialisations
    ierror     = 0                 ! set error condition
    !  set values of physical quantities
    !  and store them in quad data area
    mpi_id=iproc
    if(pftail/=0)qf_tail=pftail ! Power of parametric spectral tail
    nkq   = nk                  ! number of frequencies/wave numbers
    naq   = ndir                ! number of directions
    !VALLOC(nspec ,(nkq,naq))
    !VALLOC(aspec ,(nkq,naq))
#ifndef GETLOCN
    VALLOC(t_q ,(klocus))
#endif
    VALLOC(nspec ,(nkq,naq))
    VALLOC(q_a   ,(naq))
    VALLOC(q_k   ,(nkq))
    q_a=dirr; q_k=k; 
    call m_tools_init(ierror)
    call q_init
    if(ierror>0)return 
    do ikq = 1,nkq
      !q_lamd(ikq) = (q_k(1)/q_k(ikq))**7.5             ! used in filtering
      q_tail(ikq) =q_kfac**(-3.5*ikq) 
      q_kpow(ikq) =q_kfac**((ikq-1)) 
      q_lamd(ikq) =q_kfac**((ikq-1)*7.5) 
    enddo
    call q_makegrid
    call xnlsetgrd( naq,nkq,mkq,maq,klocus,qf_dn,qf_kn, &
         quads,quad_nloc,q_k2,q_ka,q_sig,q_sigr,q_kpow,q_tail,q_lamd);
  end subroutine
  !------------------------------------------------------------------------------
  subroutine xnlcal(ee,ks,q_dfac,xnl,diag, ierror)
    !use m_gxnldata ,only:q_k,q_dk,q_kpow,q_a,qf_krat,qf_frac,qf_dmax,q_depth
    !use m_gxnldata ,only:naq,nkq,iq_err,iq_search,q_delta
    !use m_xnldata,only:q_searchgrid,q_t13v4
    !use m_xnldata,only:nspec,q_scale
    !use m_mkgrid,only:
    !use m_tools,only:q_ctrgrid,q_error,q_init
    !use m_tools,only:q_cg
    implicit none
    !  1. Purpose:
    !     Compute nonlinear transfer for a given action density spectrum
    !     on a given wave number and direction grid
    !  2. Method
    !     Compute nonlinear transfer in a surface gravity wave spectrum
    !     due to resonant four wave-wave interactions
    !     Methods: Webb/Resio/Tracy/VanVledder
    !  3. Parameter list:
    ! Type    I/O          Name               Description
    !------------------------------------------------------------------------------
    real,   intent(in)  :: ee(nkq,naq)    ! Action density spectrum as a function of (sigma,theta)
    integer,intent(in)  :: ks
    real,   intent(in)  :: q_dfac
    real,   intent(out) :: xnl(nkq,naq)   ! nonlinear quadruplet interaction computed with
    !                                         a certain exact method (k,theta)
    real,   intent(out) :: diag(nkq,naq)  ! Diagonal term for WAM based implicit integration scheme
    integer, intent(out) :: ierror          ! error indicator
    !  4. Error messages
    !  5. Called by:
    !     XNL_MAIN
    !  7. Remarks
    !     The external action density spectrum is given as N(sig,dir)
    !     The internal action density spectrum is given as N(kx,ky)
    !     These 2 spectra are conected via the Jacobian transformation
    !                cg
    !     N(kx,ky) = -- N(sig,theta)
    !                 k
    !  8. Structure
    !  9. Switches
    ! 10. Source code
    !------------------------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------------------------------
    real diagk1      ! diagonal term for k1
    real diagk3      ! diagonal term for k3
    integer ia1,ia3,la3   ! limits for directional loops
    integer ik1,ik3,lk3   ! counters for wave number loop
    integer igrid         ! status of grid file
    real t13              ! value of sub-integral
    real qn1,qn3          ! action densities in k1 and k3
    !======================= for q_t13v4 =======================
    type(wws_type),pointer:: tq
    integer iloc          ! counter along locus
    integer nloc           ! indicator if correct locus is found
    integer ja2,ja2p      ! direction indices for interpolation of k2
    integer jk2,jk2p      ! wave number indices for interpolation of k2
    integer ja4,ja4p      ! direction indices for interpolation of k4
    integer jk4,jk4p      ! wave number indices for interpolation of k4
    real qn2,qn4          ! action densities at wave numbers k1, k2, k3 and k4
    real nprod            ! wave number product
    real zz,t2,t4            ! tail factors for k2 and k4
    real qd1,qd3          ! contribution to diagonal term
    real qn13p            ! product of N1 and N3
    real qn13d,z_lambda   ! difference of N1 and N3
    integer ibeta
    !------------------------------------------------------------------------------
    ! initialisations
    !------------------------------------------------------------------------------
    ierror = 0              ! error status
    xnl = 0.;diag=0.
    !call q_init
    !call q_ctrgrid(2,igrid)
    !if(iq_err /= 0) then
    !  ierror = 1; return
    !end if
    do ik1 = 1,nkq
      !aspec(ikq,iaq) = ee(ikq,iaq)*q_k(ikq)/q_sig(ikq)/q_cg(ikq)  
      nspec(ik1,:) = ee(ik1,:)*q_sigr(ik1)
      !aspec(ikq,iaq) = nspec(ikq,iaq)*q_k(ikq)/q_cg(ikq)  
      !nspec(ikq,iaq) = aspec(ikq,iaq)/q_k(ikq)*q_cg(ikq)  
    enddo
    !------------------------------------------------------------------------------
    !  set overall scale factor resulting from optional SEARCH for nearest grid
    !------------------------------------------------------------------------------
    q_scale = 1.; 
    !--------------------------------------------------------------------------------------
    do lk3 = 0,qf_kn                   ! 
      do la3 = -qf_dn,qf_dn            ! loop over all possible wave directions
        call q_getlocus(lk3,la3,nloc)
        if(lk3==0 .and. la3 ==0)cycle  ! skip routine if k1=k3
        do ia1 = 1,naq                 ! loop over selected part of grid, set in q_init
          ia3=mod(ia1+la3+naq-1,naq)+1
          ibeta=ia1; if(lk3<0) ibeta=ia3;
          do ik1 = 1,nkq
            ik3=ik1+lk3
            if(ik3>nkq)cycle 
            ! assume ik1<ik3
            !z_lambda=q_kpow(min(ik1,ik3)) ! q_kfac**((ik1 - 1)*7.5)
            z_lambda=q_lamd(ik1) ! q_kfac**((ik1 - 1)*7.5)
            ! 11 33 ; 13 31 ; 31 13 ; 33 11
            ! perform integration along locus
            !call q_t13v4(nspec,ik1,ia1,ik3,ia3,t13,diagk1,diagk3)   !jiangxj:diag
            t13    = 0.; diagk1 = 0.; diagk3 = 0.
            if(nloc/=0 ) then
              qn1 = nspec(ik1,ia1); qn3 = nspec(ik3,ia3)
              qn13p = qn1*qn3; qn13d = qn3-qn1
              !---------------------------------------------------------------------------------------
              !    3-----------4 ja2p         w1 = (1-wk)*(1-wa)
              !    |    .      |              w2 = wk*(1-wa)
              !    |. . + . . .| wa2   A      w3 = (1-wk)*wa
              !    |    .      |       |      w4 = wk*wa
              !    |    .      |       wa
              !    |    .      |       |
              !    1-----------2 ja2   V
              !   jk2  wk2  jk2p
              !    <-wk->
              !-----------------------------------------------------------------------------------
              !  loop over the locus
              do iloc=1,nloc
                tq=>t_q(iloc)
                zz=tq%zz;
                ja2 = tq%ia2+ibeta;ja2p = mod(ja2+naq,naq)+1;ja2  = mod(ja2-1+naq,naq)+1; 
                ja4 = tq%ia4+ibeta;ja4p = mod(ja4+naq,naq)+1;ja4  = mod(ja4-1+naq,naq)+1; 
                jk2 = tq%ik2+ik1;
                if(jk2>=nkq)then
                  !t2 = tq%wr_k2* (q_kfac)**(-3.5*ik1)
                  t2 = tq%wr_k2* q_tail(ik1)
                  qn2 =((tq%w1k2+tq%w2k2)*nspec(nkq,ja2)+(tq%w3k2+tq%w4k2)*nspec(nkq,ja2p))*t2
                else if(jk2>0)then
                  jk2p= jk2+1
                  qn2 = tq%w1k2*nspec(jk2,ja2)+tq%w2k2*nspec(jk2p,ja2)+tq%w3k2*nspec(jk2,ja2p)+tq%w4k2*nspec(jk2p,ja2p)
                else
                  !t2 = exp((tq%wr_k2+ik1-1)*q_lkfac)
                  t2  = tq%wr_k2*q_kpow(ik1)
                  qn2 =((tq%w1k2+tq%w2k2)*nspec(1,ja2)+(tq%w3k2+tq%w4k2)*nspec(1,ja2p))*t2
                endif
                jk4 = tq%ik4+ik1;
                if(jk4>=nkq)then
                  !t4 = q_kfac**(-3.5*(tq%wr_k4+ik1-nkq))
                  !t4 = tq%wr_k4* q_kfac**(-3.5*ik1)
                  t4  = tq%wr_k4* q_tail(ik1)
                  qn4 =((tq%w1k4 +tq%w2k4)*nspec(nkq,ja4)+(tq%w3k4+tq%w4k4)*nspec(nkq,ja4p))*t4
                else if(jk4>0)then
                  jk4p= jk4+1
                  qn4 = tq%w1k4*nspec(jk4,ja4)+tq%w2k4*nspec(jk4p,ja4)+tq%w3k4*nspec(jk4,ja4p)+ tq%w4k4*nspec(jk4p,ja4p)
                else
                  !t4 = exp((tq%wr_k4+ik1-1)*q_lkfac)
                  t4  = tq%wr_k4*q_kpow(ik1)
                  qn4 = ((tq%w1k4+tq%w2k4)*nspec(1,ja4) + (tq%w3k4 + tq%w4k4)*nspec(1,ja4p))*t4
                endif
                nprod = qn13p*(qn4-qn2) + qn2*qn4*qn13d !*ws
                t13   = t13 + zz*z_lambda*nprod
                qd1   = qn3*(qn4-qn2) - qn2*qn4
                qd3   = qn1*(qn4-qn2) + qn2*qn4
                diagk1= diagk1 + qd1*zz
                diagk3= diagk3 + qd3*zz
              end do
            end if
            !end q_t13v4
            !  check contribution T13 as computed with triplet method
            !  take care of additional scale factor aring from search of nearest grid
            !t13    = t13*q_scale
            !diagk1 = diagk1*q_scale
            !diagk3 = diagk3*q_scale
            if(iq_sym==1) then
              t13 = 2.*t13
              diagk1 = 2.*diagk1
              diagk3 = 2.*diagk3
            end if
            xnl (ik1,ia1) = xnl (ik1,ia1) + t13*q_k2(ik3)          !user mannual (2.33)
            xnl (ik3,ia3) = xnl (ik3,ia3) - t13*q_k2(ik1)          
            diag(ik1,ia1) = diag(ik1,ia1) + diagk1*q_ka(ik3)      
            diag(ik3,ia3) = diag(ik3,ia3) - diagk3*q_ka(ik1)      
          end do
        end do
      end do
    end do
    !call q_dscale(nspec,depth,g,q_dfac)
    do ik1 = 1,nkq
      !jacb = ccg(k,iax,icx)*sigma(k)/wk(k)
      xnl(ik1,:) = xnl(ik1,:) *q_dfac*q_sig(ik1)
    enddo
  end subroutine
  !------------------------------------------------------------------------------
#ifdef GETLOCN
  subroutine q_getlocus(kdif,adif,nloc_)
    integer, intent(in)  ::  kdif    !  k-index of wave number k1
    integer, intent(in)  ::  adif    !  k-index of wave number k1
    integer, intent(out) :: nloc_   !  indicator if reference locus exists in database
    integer kmem                 ! index for storing 2-d matrix in 1-d array
    integer amem                 ! index for storing direction of reference wave number k3
    kmem = kdif + 1;!assume kdif>=0
    if(adif>=0)then ! 0-iag2-1  maq=iag2+iag2-1
      amem=adif+1   ! 1-iag2 
    else
      amem=maq+1+adif ! -1=>maq,-2=>maq-1
    endif
    amem = mod(adif+naq,naq)+1    ! difference index
    t_q=>quads(:,kmem,amem);
  end subroutine
#else
  subroutine q_getlocus(kdif,adif,nloc_)
    implicit none
    !----------------------------------------------------------------------------------
    !  1. Purpose:
    !     Retrieve locus from basic locus as stored in the database
    !  2. Method
    !     In the case of geometric scaling, k-scaling is used using scale laws
    !     described by Tracy
    !     Directional transformation using linear transformations, shifting and mirror
    !     imaging.
    !  3. Parameter list:
    !Type      I/O           name       Description
    !------------------------------------------------------------------------------
    integer, intent(in)  ::  kdif    !  k-index of wave number k1
    integer, intent(in)  ::  adif    !  k-index of wave number k1
    !integer, intent(in)  ::  ia1    !  theta-index of wave number k1
    !integer, intent(in)  ::  ia3    !  theta-index of wave number k3
    integer, intent(out) :: nloc_   !  indicator if reference locus exists in database
    !  5. Called by
    !     Q_T13V4
    !------------------------------------------------------------------------------
    !     Local variables
    integer iadif             ! difference in angular and k-index
    integer nloc             ! number of points on locus
    integer ierr
    integer kmem                 ! index for storing 2-d matrix in 1-d array
    integer amem                 ! index for storing direction of reference wave number k3
    type(wws_type),pointer:: r_q(:)   
    integer nhalf
    integer irev,krev
    !------------------------------------------------------------------------------
    ! initialisations
    !------------------------------------------------------------------------------
    !kdif=ik3-ik1
    kmem = abs(kdif) + 1;
    !------------------------------------------------------------------------------
    !  circle grid, modify ranges and transformation variables
    !------------------------------------------------------------------------------
    !adif=ia3-ia1
    iadif = abs(adif)    ! difference index
    !0- 5  naq-5 - naq
    ! 1-naq     it1-a it1+a
    nhalf = naq/2
    ! +2 23 1
    ! -2 1 23
    if (iadif > nhalf) then
      if(adif>0)then
        adif=adif-naq
      else
        adif=adif+naq
      endif
      iadif=abs(adif)
    end if
    !A 1 3 1 3 2a 4a
    !B 1 3 3 1 2b 4b 
    !C 3 1 1 3 4b 2b
    !D 3 1 3 1 4a 2a
    ! nA = (n1*n3*(n4a-n2a) + n2a*n4a*(n3-n1))
    ! nB = (n1*n3*(n4b-n2b) + n2b*n4b*(n3-n1))
    ! nC =-(n1*n3*(n4b-n2b) + n2b*n4b*(n3-n1))
    ! nD =-(n1*n3*(n4a-n2a) + n2a*n4a*(n3-n1))
    ! NC=-NB;ND=-NA
    ! KK=k*k*delth
    ! 
    ! q=q_dfac**((min(ik1,ik3)-1)*7.5)
    ! x(1,1)+=   NA*q*KK3; ! x(3,3)-=   NA*q*KK1;
    ! x(1,3)+=   NB*q*KK3; ! x(3,1)-=   NB*q*KK1;

    ! x(3,1)+=  -NB*q*KK1; ! x(1,3)-=  -NB*q*KK3;
    ! x(3,3)+=  -NA*q*KK1; ! x(1,1)-=  -NA*q*KK3;

    ! x(1,1)+=  NA*q1*KK3 +NA*q1*KK3
    ! x(3,3)-=  NA*q1*KK1 +NA*q1*KK1
    ! x(1,3)+=  NB*q1*KK3 +NB*q1*KK3
    ! x(3,1)-=  NB*q1*KK1 +NB*q1*KK1

    ! x(1,1)+=2*NA*q1*KK3; ! x(3,3)-=2*NA*q1*KK1;
    ! x(1,3)+=2*NB*q1*KK3; ! x(3,1)-=2*NB*q1*KK1;

    !AO  a1 a3  da  t1 t3 dt iad;itm  
    ! +2 10 12;  2; 10 12; 2; 2 ;10 
    ! -2 10  8; -2; 10  8;-2; 2 ; 8 
    ! +2 23 1 ;-22; -1, 1; 2; 2 ; -1
    ! -2 1 23 ; 22;  1,-1;-2; 2 ; -1
    !(naq/2 - abs(naq/2-abs(it1-it3)))  ! compute shortest difference in indices
    !iadif = (naq - abs(naq-2*abs(it1-it3)))/2   ! compute shortest difference in indices
    !extend to all angle
    if(adif>=0)then
      amem = adif + 1   ! compute index of reference wave number k3 in interaction grid
    else
      amem = adif+naq + 1  ! compute index of reference wave number k3 in interaction grid
    endif
    !------------------------------------------------------------------------------
    ! retrieve info from reference locus in
    ! get actual number of valid points along locus (NLOCUSZ)
    ! depending on value of switch IQ_COMPACT
    !------------------------------------------------------------------------------
    nloc    = quad_nloc(kmem,amem)
    !  short-cut when number of NON-ZERO points on locus is ZERO [27/8/2003]
    if(nloc==0) return
    nloc_=nloc
    !  compute combined scale factor
    r_q=>quads(:,kmem,amem);
    !------------------------------------------------------------------------------
    ! select case to transform reference locus
    ! Transform of weigths reduces to an addition or subtraction
    ! because of log-spacing of wave numbers in the case of deep water
    if(kdif>=0 ) then      
      krev=0; irev = 0;   
      t_q(1:nloc)  = r_q(1:nloc) 
      t_q(1:nloc)  = r_q(1:nloc) 
    else 
      krev=1;  irev = 1; 
      t_q(1:nloc)%ia2  =-r_q(1:nloc)%ia4 
      t_q(1:nloc)%ia4  =-r_q(1:nloc)%ia2 
      t_q(1:nloc)%ik2  = r_q(1:nloc)%ik4
      t_q(1:nloc)%ik4  = r_q(1:nloc)%ik2
      t_q(1:nloc)%w1k2 = r_q(1:nloc)%w3k2
      t_q(1:nloc)%w2k2 = r_q(1:nloc)%w4k2
      t_q(1:nloc)%w3k2 = r_q(1:nloc)%w1k2
      t_q(1:nloc)%w4k2 = r_q(1:nloc)%w2k2
      t_q(1:nloc)%w1k4 = r_q(1:nloc)%w3k4
      t_q(1:nloc)%w2k4 = r_q(1:nloc)%w4k4
      t_q(1:nloc)%w3k4 = r_q(1:nloc)%w1k4
      t_q(1:nloc)%w4k4 = r_q(1:nloc)%w2k4
      t_q(1:nloc)%zz   = r_q(1:nloc)%zz
    endif
    !t_q(1:nloc)%sym = r_q(1:nloc)%sym
    !lambda = q_kfac**(ik1 - 1)
    !c_lambda = lambda**6
    !z_lambda=lambda*c_lambda*sqrt(lambda)

    !z_lambda = q_kfac**((ik1 - 1)*7.5)
    !t_q(1:nloc)%zz  = r_q(1:nloc)%zz*z_lambda
    !fixme

  end subroutine
#endif
  !================== NO DATA======================
  subroutine q_dscale(nspec,depth,grav,q_dfac)
    implicit none
    !  1. Purpose:
    !     Compute scaling factor for nonlinear transfer in finite depth
    !  2. Method
    !     Compute mean wave number km
    !     Compute scale factor based on parameterized function of (km*d)
    !     according to Herterich and Hasselmann
    !     and parameterisation from WAM model
    !  3. Interface parameter list:
    ! Type          I/O     Name           Description
    !-------------------------------------------------------------------------
    real,   intent(in)  :: nspec(nkq,naq)    ! Action density spectrum as a function of (sigma,theta)
    real, intent(in)     :: depth         ! Depth (m)
    real, intent(in)     :: grav          ! Gravitational acceleration
    real, intent(out)    :: q_dfac        ! scale factor
    !  5. Called by:
    !     XNL_MAIN
    !------------------------------------------------------------------------------
    real w          ! radian frequency
    real kk         ! local wave number
    real sqkk       ! square root of local wave number
    real dnn        ! summation quantity
    real kms        ! mean wave number
    real kd         ! depth*mean wave number product
    real sum0       ! summation variable for total energy
    real sumk       ! summation variable for wave number
    integer isig    ! counter over sigma loop
    integer iang    ! counter over direction loop
    !------------------------------------------------------------------------------
    sum0 = 0.
    sumk = 0.
    !  compute sums for total energy andwave number
    do isig = 1,nkq
      kk   = q_k(isig) 
      sqkk = sqrt(kk)
      do iang=1,naq
        !nwam aspec(k,:) = e(k,:,iax,icx)*wk(k)/sigma(k)/ccg(k,iax,icx)  
        !nwam nspec(ikq,:) = aspec(ikq,:)/q_k(ikq)*cg(ikq)
        !O dnn  = aspec(isig,iang)*q_dsig(isig)*q_delta      
        !T aspec(ikq,iaq) = nspec(ikq,iaq)*q_k(ikq)/q_cg(ikq)  
        !T dnn  = nspec(ikq,iaq)*q_k(ikq)/q_cg(ikq)*q_dsig(isig)*q_delta      
        dnn  = nspec(isig,iang)*q_dk(isig)*q_delta   
        sum0 = sum0 + dnn   !sum(e*dsig*dthet)
        sumk = sumk + 1./sqkk*dnn !sum(e*dsig*dthet/sqrt(k))
      end do
    end do
    !  compute mean wave number and scale factor based
    !  on the WAM approximation
    if(sum0 > 0) then
      kms = (sum0/sumk)**2   !ark
      kd = max(0.5,0.75*kms*depth)
      q_dfac = 1+5.5/kd*(1.-5./6.*kd)*exp(-5./4.*kd)
    else
      kms = 0.
      kd  = 0.
      q_dfac = 1.
    end if
  end subroutine
end module
subroutine xnl_inite(k,dirr,nk,ndir,pftail,iquad,iproc,ierror)
  use m_xnldata
  !  1. Purpose:
  !     Initialize coefficients, integration space, file i/o for computation
  !     nonlinear quadruplet wave-wave interaction
  !  2. Method
  !     Set version number
  !     Set unit unit numbers
  !     Open quad related files
  !     Optionally reset configuration by a back door option
  !     Compute integration spaces for given water depths
  !  3. Parameter list:
  !Type   I/O              Name          Description
  !------------------------------------------------------------------------------
  integer, intent(in)  ::  nk            ! Number of sigma values
  integer, intent(in)  ::  ndir          ! Number of directions
  real(8), intent(in)  ::  k(nk)         ! Radian frequencies
  real(8), intent(in)  ::  dirr(ndir)    ! Directions (degrees)
  real(8), intent(in)  ::  pftail        ! power of spectral tail, e.g. -4 or -5
  integer, intent(in)  ::  iquad         ! Type of method for computing nonlinear interactions
  integer, intent(in)  ::  iproc         ! Processor number, controls output file for MPI
  integer, intent(out) ::  ierror        ! Error indicator. If no errors are detected IERR=0
  call xnl_init(k,dirr,nk,ndir,pftail,iquad,iproc,ierror)
end subroutine
