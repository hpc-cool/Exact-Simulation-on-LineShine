#include "wavedeff.h"
#define USE_SMOOTH8
module propagat_mod
  use varcommon_mod
  use partition_mod
IMPLICIT NONE
  real   (8),allocatable::wkd(:)
  real   (8),allocatable::dxym(:,:)
  real (8),pointer::pvkdp(:,:,:)
  type  spc_interg
    real(4) pab(4)
    integer*2 kpos(4)
  end type spc_interg
  type  geo_interg
    real(4) pabp(4)
    integer*4 iquad
  end type geo_interg
  type(spc_interg),pointer:: ps_vs(:,:)
  type(geo_interg),pointer:: pg_vs(:,:)
  integer mkj,iact
!NSMOOTH =4 | 8
#define NSMOOTH 4

 integer*1 LID(1)
 integer kjt,ixx,jyy,ixx1,jyy1
#define ICCG   1
#define IDSIDD 2
#define IMDPSP 2
! ICCG   = 1, cg                  | used implsch_mod & propagat_mod.F90
! IDSIDD = 2, dsidd               | used propagat_mod
contains
  SUBROUTINE InitPropagata !{
    integer k
    mkj=kl*jnthet
    ALLOCATE(wkd(kld+1)) ;
    DO  k=1,kld !{
      if(k>1)wkd(k-1)=(wk(k)-wk(k-1))
    end do !}
  END SUBROUTINE InitPropagata !{

#ifndef C_PROPAGAT
  
  SUBROUTINE propagats_spec1(ed,es,iac) !{
    REALD,intent(in ):: es(kl*jnthet)
    REALD,intent(out):: ed(kl*jnthet)
    integer,intent(in):: iac
    integer kj
    type(spc_interg),pointer:: ppv
#define PE(L) ppv%pab(L)*es(ppv%kpos(L))
    DO  kj=1,jnthet*kl !{
      ppv=>ps_vs(kj,iac)
      ed(kj)=PE(1)+PE(2)+PE(3)+PE(4)
    end do
#undef PE
  END SUBROUTINE 
  SUBROUTINE propagats_spec(ed,es,iacb,iace) !{
    REALD,intent(in ):: es(kl*jnthet,0:nwpc)
    REALD,intent(out):: ed(kl*jnthet,0:nwpc)
    integer,intent(in ):: iacb,iace
    integer iac
    do iac=iacb,iace
      call propagats_spec1(ed(1,iac),es(1,iac),iac)
    end do
  END SUBROUTINE 
  SUBROUTINE propagats_geo(ed,es,iacb,iace) !{
    REALD,intent(in ):: es(kl*jnthet,0:nwpc)
    REALD,intent(out):: ed(kl*jnthet,0:nwpc)
    integer,intent(in ):: iacb,iace
    integer kj,iac,ppos(4),iquad,iquad_pre
    type(geo_interg),pointer:: ppv
    do iac=iacb,iace
      iquad_pre=-1
      ppos(1)=iac
      DO  kj=1,jnthet*kl !{
        ppv=>pg_vs(kj,iac)
        iquad=ppv%iquad*3+1;
        ppos(2:4)=ipos12(iquad:iquad+2,iac);
        iquad_pre=iquad
#define PE(L) ppv%pabp(L)*es(kj,ppos(L))
        ed(kj,iac)=PE(1)+PE(2)+PE(3)+PE(4)
#undef PE
      end do
    end do
  END SUBROUTINE !}

#endif
  SUBROUTINE GetIKJ(k,j,k_d,th_d,kpos,pab) !{
    integer k,j
    real(8) k_d,th_d
    real(8) wks,ths_
    real(4) pab(4)
    integer kpos(4)
    integer iwk,iwk1,jth,jth1
    real(8) p,q,r,ths
    integer ik,kj
    kj=(j-1)*kl+k
    if(abs(th_d)<1e-10)th_d=0
    ths=thet(j)+th_d;
    if(ths>=zpi)then
      ths=ths-zpi
    else if(ths<0)then
      ths=ths+zpi
    endif
    jth=idint(ths/deltth)+1;
    if(ths<thet(jth))then
      jth=jth-1
    endif
    jth1=jth+1;
    if(jth>=jnthet)then
      jth=jnthet
      jth1=1
    endif
    r=(ths-thet(jth))*ddeltth
    if(r>1)r=1
    if(abs(k_d)<1e-10)k_d=0
    wks=wk(k)+k_d
    if(k_d>=0)then
      if(wks>=wk(kld))then
        iwk=kl;iwk1=kl;
        p=wkh(kld-kl+1)
        q=0
      else
        iwk=k+1;    
        do while(wks>=wk(iwk))
          iwk=iwk+1;
        enddo
        iwk=iwk-1
        if(iwk<kl)then
          iwk1=iwk+1;
          q=(wks-wk(iwk))/wkd(iwk)
          p=1-q
        else if(iwk<kld)then
          q=(wks-wk(iwk))/wkd(iwk)
          p=(1-q)*wkh(iwk-kl+1)+q*wkh(iwk-kl+1+1)
          iwk=kl;iwk1=kl;
          q=0
        end if
      endif
    else
      if(wks<wk(1))then
        iwk=1;iwk1=1;
        if(wks>0)then
          q=wks/wkmin
          p=0;
        else
          p=0;q=0
        endif
      else
        iwk=k-1
        
        do while(iwk>0 .and. wks<wk(iwk))
          iwk=iwk-1;
        enddo
        q=(wks-wk(iwk))/wkd(iwk)
        iwk1=iwk+1;p=1-q
      endif
    endif
    pab(1)=p*(1-r);  pab(2)=q*(1-r)
    pab(3)=p*r    ;  pab(4)=q*r
    kpos(1)=(jth -1)*kl+iwk;  kpos(2)=(jth -1)*kl+iwk1;
    kpos(3)=(jth1-1)*kl+iwk;  kpos(4)=(jth1-1)*kl+iwk1;
    if(minval(pab)<-1e-10)then
      write(6,*)"K",k_d,th_d,q,r,iwk,iwk1,jth,jth1,wkd(iwk),wks,thet(jth),ths,zpi
    endif

  end SUBROUTINE GetIKJ !}
  SUBROUTINE GetIPOS(x_d,y_d,iac,iquad,pabp) !{
    real(8) x_d,y_d
    real(4) pabp(4)
    integer iac,iquad,ia,ic
    real(8) q,r
      !   4---3---2
      !   | 1 | 0 |
      !   5---+---1
      !   | 3 | 2 |
      !   6---7---8
    if(abs(x_d)<1e-5)x_d=0
    if(abs(y_d)<1e-5)y_d=0
    if(y_d>=0)then !{
      if(x_d>=0)then !{
        iquad=0
        q= x_d/dxm(iac)
        r= y_d/dym(iac)
      else     !}{
        iquad=1
        q=-x_d/dxm(iac)
        r= y_d/dym(iac)
      end if
    else
      if(x_d>=0)then !{
        iquad=2
        q= x_d/dxm(iac)
        r=-y_d/dym(iac)
      else
        iquad=3
        q=-x_d/dxm(iac)
        r=-y_d/dym(iac)
      end if !}
    end if !}
    pabp(1)=(1-q)*(1-r)
    pabp(2)=(q)*(1-r)
    pabp(3)=(1-q)*(r)
    pabp(4)=(q)*(r)
    if(minval(pabp)<-1e-10)then
      call id2pos(iac,ia,ic)
      write(6,*)"P",x_d,y_d,iquad,r,q,dxym(5,iac),dxym(6,iac),ia,ic
      call flush(6)
      stop
    endif
  end SUBROUTINE GetIPOS !}
  SUBROUTINE propagats_pre(iacb,iace)
    integer iacb,iace
    real(8) d0,dddx0,dddy0
    real(8) wk0,dk0,cg,cgx,cgy
    real(8) dsidd,ssrwk,ssrth,ssrwkd,ssrthd,ssrwkbc,ssrthbc
    real(4) pab(4),pabp(4)
    real(8) rt0,tRsd_tanLat
    integer iac,ia,ic,iquad,kpos(4),ppos(3)
    real(8)::x_d,y_d,k_d,th_d,rsec
    real(8) ttd
    integer j,k,kj,isec
    integer nnt(10),nt
#if LOGSCURR==1
    real(8) duxdx0,duxdy0,duydx0,duydy0,ux,uy,dudsl,dudnl
    ux=0;uy=0;duxdx0=0;duxdy0=0;duydx0=0;duydy0=0
#endif
    nnt=0; pab=0;pabp=0
    !ttd=Difftimer(2)
    !call zlock
    !write(*,'(i8,f8.3,a,5i10.6)')iwalltime(),ttd,' ABB',iacb,iace,iacb,iacb-iacb,iace-iacb+1
    !call flush(6)
    !call zunlock
    do iac=iacb,iace
      !call zlock
      !print'(a,3i8.6,a)','AB',iacb,iace,iac,'B';
      !call flush(6)
      !call zunlock

      nt=0;
      iact=iac
      tRsd_tanLat=Rsd_tanLat(iac)
      d0   =dep(iac)
      dddx0=dddx(iac)
      dddy0=dddy(iac)
#if LOGSCURR==1
      if(logscurr/=0)then
        ux    =ucur(1,iac)
        uy    =ucur(2,iac)
        duxdx0=ucur(3,iac)
        duxdy0=ucur(4,iac)
        duydx0=ucur(5,iac)
        duydy0=ucur(6,iac)
      end if
#endif

      DO  j=1,jnthet !{
        ssrwkd=-( dddx0*cosths(j)+dddy0*sinths(j))
        ssrthd=-(-dddx0*sinths(j)+dddy0*cosths(j))
        rt0=cosths(j)*tRsd_tanLat !全球模式，大圆 系数  ! Big Circle
        ssrthbc=0 ;       ssrwkbc=0

#if LOGSCURR==1
          dudsl=duxdx0*cosths2(j)+(duxdy0+duydx0)*sincosths(j)+duydy0*sinths2(j)
          dudnl=(duydy0-duxdx0)*sincosths(j)+duxdy0*cosths2(j)-duydx0*sinths2(j)
          !波流 大圆  ! Big Circle
          ! !1:ux;2:uy;3:uxx;4:uxy;5:uyx;6:uyy
          ssrwkbc=(ux*sinths(j)-uy*cosths(j))*rt0
          ssrthbc=(ux*cosths(j)+uy*sinths(j))*rt0  ! + uy*tRsd_tanLat
#endif
        DO  k=1,kl !{
          wk0=wk(k)
          cg=pvkdp(k,ICCG,iac)
          ssrth=0 +cg*rt0 !大圆  ! Big Circle
          ssrwk=0;
          cgx=cg*cosths(j);cgy=cg*sinths(j)
#if LOGSCURR==1
            !******  1.  "the calculation of wave engery-current spreading"
            cgx=cgx+ux;cgy=cgy+uy
#endif
          !======================================================================c
          !******  2.  "the effect of refraction caused by topography and current"
            dsidd=pvkdp(k,IDSIDD,iac)
            ssrwk=ssrwk+ssrwkd*dsidd
            ssrth=ssrth+ssrthd*dsidd/wk0
          !if(iact==iacp(1))then
          !! write(3163,'(2i5,20e26.13)')k,j,ssrwk,dddx0,dddy0,dsidd
          !  !write(6,'(6i5," A",8e13.7)')k,j,kpos,wk0+k_d,pab
          !endif

#if LOGSCURR==1
          ssrwk=ssrwk-(dudsl*wk0)
          ssrth=ssrth-(dudnl) !*wk0/wk0
          ssrwk=ssrwk+ssrwkbc
          ssrth=ssrth+ssrthbc
#endif

          x_d=-deltts*cgx;y_d=-deltts*cgy
          k_d=-deltts*ssrwk;th_d=-deltts*ssrth

          call id2pos(iac,ia,ic)


          IF(abs(x_d)>dxm(iac).or.abs(y_d)>dym(iac).or.abs(th_d)>pi/2) THEN !{
            call id2pos(iac,ia,ic)
            DBGO(0,*) 'deltts too big, be careful !!! program have to stop !!!'
            DBGO(0,*) 'decrease deltts in nwamctl.ini please '
            DBGO(0,*) 'Note:86400( seconds of 1day ) MUST Divide Exactly deltts  '
            !86400=128*27*25  ;3600 =16*9*25 1440=32*9*5
            ! 1 2 4 8 16 32
            ! 3 6 12 24 48 96
            ! 9 18 36 72 144 288
            ! 5 10 20 40 80 160
            ! 15 30 60 120 240 480
            ! 45 90 180 360 720 1440
            ! 1 2 3 4 5 6 10 12 15 20 30 60 min
            !   8   9 16 18 24 32 36 40 45 48 72 80 90 96 120 144 160 180 240 288 360 480 720 1440 min
            ! 180 160 90 80 60 45 40 36 32 30
            ! 1 2 3 4 5 6 8 9 12 15 18,20,24,25,27,30,36,40,
            ! 32
            DBGO(0,*) 'at mpi_id:',mpi_id,iac,kj
            DBGO(0,*) 'ia=',ia,'ic=',ic,'k=',k,'j=',j,'iac=',iac,'x=',sxcord(iac),'y=',sycord(iac)
            DBGO(0,*) 'abs(x_d):abs(',x_d,')must <',dxm(iac)
            DBGO(0,*) 'abs(y_d):abs(',y_d,')must <',dym(iac)
            DBGO(0,*) 'abs(th_d):abs(',th_d,')Must <',pi/2
            rsec=abs(dxm(iac)/cgx);
            rsec=min(rsec,abs(dym(iac)/cgy));
            rsec=min(rsec,abs(pi/ssrth));
            isec=rsec;
            do while(mod(86400,isec)/=0)
              isec=isec-1
            enddo
            DBGO(0,*) 'try deltts from ',deltts,' to ',isec,'Sec',rsec,'sec'
            call wav_abort('deltts too big')
            stop
          end if !}
          call id2pos(iac,ia,ic)
          call GetIPOS(x_d,y_d,iac,iquad,pabp)
          call GetIKJ(k,j,k_d,th_d,kpos,pab)
          kj=(j-1)*kl+k
          if(mod(k,4)/=1)THEN
            if(kpos(1)-ps_vs(kj-1,iac)%kpos(1)/=1)THEN
              nt=nt+1;nnt(1)=nnt(1)+1
            endif
          endif
          pg_vs(kj,iac)%iquad=iquad ;
          ps_vs(kj,iac)%kpos =kpos  ! not -1 for fortran
          ps_vs(kj,iac)%pab  =pab   ;
          pg_vs(kj,iac)%pabp =pabp  ;
          if(minval(pab)<-1e-10.or.minval(pabp)<-1e-10)THEN
            write(6,'(4i5," A",16e11.3)')mpi_id,k,j,iac,k_d,th_d,thet(j)+th_d,x_d,y_d,pab,pabp
            call flush(6)
            stop
          endif
        end do
      end do
      if(nt/=0)THEN
        nnt(2)=nnt(2)+1
      endif
      !if(mod(iac-iacb,1000)==0)then
      !  call zlock
      !  print'(a,5i8.6,a)','AB',iacb,iace,iac,iac-iacb,iace-iacb,'E';
      !  call flush(6)
      !  call zunlock
      !endif
    end do
    !ttd=Difftimer(2)
    !call zlock
    !write(*,'(i8,f8.3,a,5i10.6)')iwalltime(),ttd,' ABE',iacb,iace,iac,iac-iacb,iace-iacb+1
    !call flush(6)
    !call zunlock
  end SUBROUTINE propagats_Pre !}
  SUBROUTINE InitSmooth
#ifdef USE_SMOOTH8
    smooth_p0=1./16./8.
#else
    smooth_p0=1./16./4.
#endif
  end SUBROUTINE InitSmooth
  SUBROUTINE smooth_e(ed,es,nwpb,nwpe) !{
    REALD ed(kl,jnthet,0:nwpc)
    REALD es(kl,jnthet,0:nwpa)
    integer nwpb,nwpe
    real(8) et,ev
    integer iac,j,k,iac1,iac3,iac5,iac7
    real(8) sp1,sp3,sp5,sp7

#ifdef USE_SMOOTH8
    integer iac2,iac4,iac6,iac8
    real(8) sp2,sp4,sp6,sp8
#endif
      !   4---3---2
      !   |   |   |
      !   5---+---1
      !   |   |   |
      !   6---7---8
    !e_c_n=e_c_n*(1-nn*p0) +sum(e_o_t)*p0 +(e_c_n-e_c_o)*nn*p0
    !e_c_n=e_c_n +(sum(e_o_t)-e_c_o*nn)*p0
    !e_c_n=e_c_n +(sum(e_o_t-e_c_o))*p0
    do iac=nwpb,nwpe
      iac1=ipos8(iac)%i8(1);if(iac1>0)then;sp1=smooth_p0;else;sp1=0.;endif
      iac3=ipos8(iac)%i8(3);if(iac3>0)then;sp3=smooth_p0;else;sp3=0.;endif
      iac5=ipos8(iac)%i8(5);if(iac5>0)then;sp5=smooth_p0;else;sp5=0.;endif
      iac7=ipos8(iac)%i8(7);if(iac7>0)then;sp7=smooth_p0;else;sp7=0.;endif
#ifdef USE_SMOOTH8
      iac2=ipos8(iac)%i8(2);if(iac2>0)then;sp2=smooth_p0;else;sp2=0.;endif
      iac4=ipos8(iac)%i8(4);if(iac4>0)then;sp4=smooth_p0;else;sp4=0.;endif
      iac6=ipos8(iac)%i8(6);if(iac6>0)then;sp6=smooth_p0;else;sp6=0.;endif
      iac8=ipos8(iac)%i8(8);if(iac8>0)then;sp8=smooth_p0;else;sp8=0.;endif
#endif

      DO  j=1,jnthet !{
        DO  k=1,kl !{
          et=0.;ev=es(k,j,iac)
          et=et+(es(k,j,iac1)-ev)*sp1
          et=et+(es(k,j,iac3)-ev)*sp3
          et=et+(es(k,j,iac5)-ev)*sp5
          et=et+(es(k,j,iac7)-ev)*sp7
#ifdef USE_SMOOTH8
          et=et+(es(k,j,iac2)-ev)*sp2
          et=et+(es(k,j,iac4)-ev)*sp4
          et=et+(es(k,j,iac6)-ev)*sp6
          et=et+(es(k,j,iac8)-ev)*sp8
#endif
          ed(k,j,iac)=ev+et
        enddo
      enddo
    enddo
  end SUBROUTINE !}
end module propagat_mod
SUBROUTINE InitPropagat(iacb,iace) !{
  use varcommon_mod
  !use partition_mod
  use propagat_mod
  IMPLICIT NONE
  integer iacb,iace
  integer k,iac,ia,ic,iad,icd,iaa,ica,i
  real(8) dxmt,dymt
  real(8) wkk,dk,wsk,cg,tanhdk
   
  do iac=iacb,iace
    IF(nsp(iac)/=1) CYCLE
    call id2pos(iac,ia,ic)
    if(ia<=0.or.ic<=0)then
        write(6,*)mpi_id,'InitPropagat Error:',iac,nsp(iac),ia,ic
        stop
        call wav_abort('InitPropagat Error:1')
    endif
    !iad=adjix(ia-1);icd=adjiy(ic-1)
    !iaa=adjix(ia+1);ica=adjiy(ic+1)
    !if(icd<rectc(2).or.icd>rectc(4))icd=ic; !ZZZEERR
    !if(ica<rectc(2).or.ica>rectc(4))ica=ic; !ZZZEERR
    !   4---3---2
    !   |   |   |
    !   5---+---1
    !   |   |   |
    !   6---7---8
    !i  2  3  4
    !132534178576

    dymt=deg2m*grdszy
    dxmt=deg2m*grdszx
    dxym(1,iac)=  dxm(iac);  dxym(2,iac)=dym(iac)
    dxym(3,iac)=  dxm(iac);  dxym(4,iac)=dym(iac)
    dxym(5,iac)=  dxm(iac);  dxym(6,iac)=dym(iac)
    dxym(7,iac)=  dxm(iac);  dxym(8,iac)=dym(iac)
    DO  k=1,kl !{
      wkk=wk(k)
      dk=dep(iac)*wkk
      !avoid float overflow
      !  sinh(80)   =2.87e+34
      !  cosh(40)**2=1.38e+34,
      !1-tanh(40)   =3.61e-35
      if(dk>40)then  
        wsk=sqrt(g*wkk)
        pvkdp(k,ICCG  ,iac)=0.5*wsk/wkk
        pvkdp(k,IDSIDD,iac)=0
      else
        wsk=sqrt(g*wkk*tanh(dk))
        pvkdp(k,ICCG  ,iac)=0.5*wsk*(1.+2.*dk/sinh(2.*dk))/wkk
        pvkdp(k,IDSIDD,iac)=0.5*g/cosh(dk)*wkk**2/wsk/cosh(dk)
      endif
    enddo
  end do
  call propagats_Pre(iacb,iace)
  if(iacb<10) call InitSmooth
end SUBROUTINE InitPropagat !}
