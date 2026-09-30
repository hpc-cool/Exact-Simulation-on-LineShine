module wamxnl_mod
  !================================================================
  ! -- module for exact calculation of wave-wave non-linear energy transfer 
  !    with the WRT method, codes based on ww3 ver. 3.14
  ! ******
  !    module used:
  !                               ------m_constants
  !                    (ww3)      |
  !     wamxnl_mod ---mod_xnl4v5--------m_fileio
  !                               |
  !                               ------serv_xnl4v5
  ! ******
  !    main subroutine structrue:
  !     
  !      wamxnl ----- xnl_init
  !               |
  !               --- xnl_main-----q_xnl4v4-----q_t12v4
  ! ******
  !    usage: used in implsch instead of snl codes
  !              wamxnl(ia,ic,se,dse)
  ! ******
  !    *NOTE* : 1.ftail = 4, according to k**-7/2 in MASNUM-WAM
  !             2.iqgrid = 3, for full circle and non-symmetric
  !             3.iquad = 2, defined in namelist of *.inp files:
  !                      1-deep water/2-deep water & WAM scaling/3-finite depth
  !                      when 1 or 2 is chosen, depth will be forced to 1000m
  !             4.ndepth = 1, only 1 can be chosen becuse only 1 level of depth 
  !                      could be found in MASNUM-WAM
  !             5.in subr. xnl_init, sigma is in rad/s,
  !                                  dird is in degree which means 360 is a circle
  !             6.in subr. xnl_main, sigma is in rad/s,
  !                                  angle(dird) is in rad, means 2*pi is a circle
  !                                  aspec is the action density spectrum, as a 
  !                                  function of (sig,theta) in witch sig is in rad/s,
  !                                  theta is in rad.
  !             7.xnl: calculated by subr. q_xnl4v4: 
  !                                  xnl = dN(sig,theta)/dt  
  !                                  se = dE(kx,ky)/dt
  !                             ***  se(k,:) = xnl(k,:)*ccg(k,ia,ic)*sigma(k)/wk(k)
  !
  !               diag: calculated by subr. q_xnl4v4:
  !                                  dN       1    dE     1
  !                         SS(N) = ----- = -----x---- = ----xSS(E)
  !                                  dt      sig   dt    sig
  !                                 dSS(N(kx,ky))   d[1/sig x SS(E)]
  !                    ***  diag = ---------------=----------------
  !                                       dN         d[1/sig x E]
  !                              = dSS(E(kx,ky))/dE = dse  ***
  !                         that is diag == dse          
  ! ******                                         --------------jiangxj----------------
  !                                          -----ver. 1.0---2014.6.9-------------
  !                         and so does diag(sig,theta) = diag(kx,ky) 
  !================================================================ 
  implicit none
  public wamxnl,init

  private
  real(8), parameter :: wkmax     = 0.6894
  real(8), parameter :: wkmin     = 0.0071
  real(8), parameter :: pwk       = 1.21
  real(8),parameter::pi=3.14159265358979D0
  real(8),parameter::zpi=pi*2,pid=1./pi
  real,parameter:: g=9.81
  integer ::kl=25,jl=24,ndep=1,kld=30,kldp1,jlp1
  integer ::ixs=1,iys=1,ixl=1,iyl=1;
  real,pointer::wk(:),thet(:),d(:,:)
  real,pointer::dwk(:),wkh(:),grolim(:),cosths(:),sinths(:)
  real,pointer::ws(:,:,:),wf(:,:,:),ccg(:,:,:),dwf(:,:,:)
  real,pointer::se(:,:),dse(:,:),e(:,:,:,:)
contains
  subroutine wamxnl(iax,icx)  !,xnl,diag)
    use m_xnldata, only: xnl_init,xnl_main,q_dscale
    !use wamvar_mod, only: kl,jl,wf,wk,ccg,thet,e,d,zpi,g,se,dse
    implicit none 
    integer,intent(in)        ::    iax,icx
    real                      ::    nspec(kl,jl)
    real                      ::    xnl(kl,jl),diag(kl,jl)  !jiangxj:se & dse in masnum-wam implsch
    real,      parameter      ::    ftail = -5.     !~k**-7/2 ((-7/2)+1)*2
    integer,   parameter      ::    iqgrid = 3,  &  !full circle and non-symmetric 
      iquad = 2,   &  !1-deep water/2-deep water & WAM scaling/3-finite depth
      ndepth = 1      !jiangxj:IQUAD = 3 can not be used because ndepth always .eq. 1 in MASNUM-WAM, so shallow water is not considered.
    integer   :: ierr   !error infomation from 'xnl_init'
    real      :: jacb,q_dfac
    real      :: depth(1)
    integer   :: i,j,k,ks
    ! --------------------------------------------------------------------------
    !sigma(:) = wf(1:kl,iax,icx)*zpi   !in rad/s
    !sigma=sqrt(g*wk*tanh(d*wk))
    !wk*tanh(d*wk)=sigma*sigma/g
    !sig(i+1)-dig(i)=sqrt(g*wk(i+1)*tanh(d*wk(i+1)))-sqrt(g*wk(i)*tanh(d*wk(i)))
    !DS=0.5/sqrt(g*wk*tanh(d*wk))*(g*tanh(d*wk)+g*wk*d*(cosh(d*wk)^-2)
    xnl=0.;diag=0.
    ks=kl;
    call xnl_init ( wk, thet, kl, jl, ftail,  iquad, 1, ierr )

    !EK=ES*CG     ES=EK/CG
    !WW3
    !  ES*CG/SIG  AW 
    !  ES/SIG     A 
    !  ES/WK/SIG*CG
    !NWAM
    !  EK/WK        ES*CG/WK
    !  EK/SIG/CG    ES/SIG
    !  EK/WK/SIG    ES/WK/SIG*CG
    !SWAN
    !  ES/SIG
    !  ES/SIG
    !  ES/WK/SIG*CG
    !  convert input action density spectrum from A(sigma,theta) -> N(kx,ky)
    !do k = 1,kl
    !  !aspec(ikq,iaq) = ee(ikq,iaq)*q_k(ikq)/q_sig(ikq)/q_cg(ikq)  
    !  nspec(k,:) = e(k,:,iax,icx)/ws(k,iax,icx)
    !  !aspec(ikq,iaq) = nspec(ikq,iaq)*q_k(ikq)/q_cg(ikq)  
    !  !nspec(ikq,iaq) = aspec(ikq,iaq)/q_k(ikq)*q_cg(ikq)  
    !enddo
    call q_dscale(nspec,depth(1),g,q_dfac)
    call xnlcal ( e(:,:,iax,icx),ks,q_dfac,xnl, diag, ierr )
    !xnl = xnl*q_dfac
    !do k = 1,kl
    !  !jacb = ccg(k,iax,icx)*sigma(k)/wk(k)
    !  jacb=ws(k,iax,icx)
    !  se(k,:) = xnl(k,:) *jacb               !jiangxj:dN/dt,N(sig,theta),se~E(kx,ky)
    !  dse(k,:) = diag(k,:)
    !enddo
  end subroutine wamxnl
  subroutine init(dep)
    real dep
    jlp1=jl+1;kldp1=kld+1
    allocate(wk(kldp1),thet(jlp1),d(ixl,iyl),cosths(jl),sinths(jl))
    allocate(dwk(kldp1),wkh(kldp1),grolim(kldp1))
    allocate(ws(kldp1,ixl,iyl),wf(kldp1,ixl,iyl),ccg(kldp1,ixl,iyl),dwf(kldp1,ixl,iyl))
    allocate(se(kl,jl),dse(kl,jl),e(kl,jl,ixl,iyl))
    d=dep;
    call setwave
    call setspec(10.,0.)
  end subroutine init

  subroutine setwave
    implicit none
    integer :: j, k,ia,ic
    real(8) :: wh, di, wkk, dk, tanhdk, wfk, wsk
    real(8) :: cgro,deltth
    deltth=zpi/float(jl)
    do j=1,jlp1
      thet(j)=(j-1)*deltth
      cosths(j)=cos(thet(j))
      sinths(j)=sin(thet(j))
    enddo
    wh=sqrt((1./pwk)**7)
    wkh(1)=1.
    do k=1,kldp1
      wk(k)=wkmin*(pwk**(k-1)) !discretion of wave number
      if (k.le.kld) dwk(k)=(pwk-1.)*(wk(k)**2)*deltth/2
      if (k.ge.2) wkh(k)=wkh(k-1)*wh
    enddo

    do ic=iys,iyl
      do ia=ixs,ixl
        di=d(ia,ic)
        do k=1,kldp1
          wkk=wk(k)
          dk=di*wkk
          tanhdk=1.
          if (dk.lt.4.) tanhdk=tanh(dk)
          wsk=sqrt(g*wkk*tanhdk)
          wfk=wsk/zpi
          ws(k,ia,ic)=wsk
          wf(k,ia,ic)=wfk
          if (dk.gt.4.) then
            ccg(k,ia,ic)=0.5*wsk/wkk
          else
            if (dk.lt.0.14) then
              ccg(k,ia,ic)=sqrt(g*di)
            else
              ccg(k,ia,ic)=0.5*wsk*(1.+2.*dk/sinh(2.*dk))/wkk
            endif
          endif
        enddo
      enddo
    enddo
    do ic=iys,iyl
      do ia=ixs,ixl
        do k=1,kld
          dwf(k,ia,ic)=(wf(k+1,ia,ic)-wf(k,ia,ic))*deltth/2
        enddo
      enddo
    enddo
  end subroutine setwave
  subroutine setspec(vx,vy)
    implicit none
    real ,intent(in):: vx, vy
    real, parameter :: gama=3.3, sq3=3.0**0.5
    real::ww,wwr,wwr2, xj, xj0, arlfa, wsj, wkj,wsjr
    real::theta0, sinth, costh, wk0r,  ws0, wl, sigmar, alpha
    integer :: j, k,ia,ic
    do ia=1,ixl
      do ic=1,iyl
        ww=vx**2+vy**2
        if (ww.le.0.) ww=0.9
        wwr=1./ww; wwr2=wwr*wwr
        xj0=200.*1000.; xj=g*xj0*wwr2
        arlfa=(0.076*(xj**(-0.4)))*pid*wwr2
        wsj=22.*(xj**(-0.33))*g*wwr
        wsjr=1./wsj; wkj=wsj**2/g
        do k=1,kl
          wk0r=1./wk(k); ws0=ws(k,ia,ic)
          if (ws0.le.wsj) then
            sigmar=1./0.07
          else
            sigmar=1./0.09
          endif
          alpha=arlfa*wk0r**4* exp(-1.25*(wkj*wk0r)**2)     &
            *gama**(exp(-0.5*((1.-ws0*wsjr)*sigmar)**2))
          do j=1,jl
            wl=vx*cosths(j)+vy*sinths(j)
            e(k,j,ia,ic)=alpha*(wl*wl)
          enddo
        enddo
      enddo
    enddo
  end subroutine setspec
end module
program xnltest
  use wamxnl_mod
  call init(1000.)
  call wamxnl(1,1)
end program
