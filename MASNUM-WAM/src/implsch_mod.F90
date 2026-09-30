#include "wavedeff.h"
#define USEINTERACTNEW
!源函数
Module implsch_mod
use varcommon_mod
use partition_mod
IMPLICIT NONE
    real(8),allocatable:: wpsk(:,:),cwks17(:)
    real(8),allocatable:: grolim(:)
    integer,allocatable::wmmask(:)
    !setted by nlweight & used by implsch
    integer ,allocatable::jps(:,:,:),ikps(:,:)

    real(8) wps(16)
    real(8):: ads=1.,abo=1.,sbo
    real(8)::cksp=14.0,cksa=4.5
    real(8) alv(5)

  type  implsch_pack_type
    real(8) WS    (CKL) !sqrt(g*wkk*tanhdk)  |
    real(8) IWSDK (CKL) !(1/wsk)*dwk         | used mean2
    real(8) SSBO  (CKL) !ssbo                |
!   real(8) CCG   (CKL) !cg                  |
!#if LOGSCURR==1
!   real(8) CCGD  (CKL) !cg*wkk/wsk          | used sscu when use current
!#endif
    real(8) depth,Rsd_tanLat
  end type implsch_pack_type
  type(implsch_pack_type) ,pointer::ip(:)!计算优化 深度相关
  integer mkj
#ifdef C_IMPLSCH
  type  implsch_gpar_type
    ! CKL must be divisible by 4
    real(8) wkal(4,12);
    integer(8) kjps(2,12,CJNTHET);
    integer(1) kjs(2,2,12,CJNTHET);
    real(8) cwks17(CKL),grolim(CKL);
    real(8) wk(CKL),wkh(CKL),dwk(CKL),wkdk(CKL),wkibdk(CKL);
#if SIGMALVL==2
    double expdep(CKL*MBVDEP);
#endif

    real(8)    cosths(CJNTHET),sinths(CJNTHET);
    real(8)    lmdpd,lmdmd,lmdpd20,lmdmd20,Dm2,D1;
    real(4)    beta10,beta11;
    real(4)    ads,brkd1,brkd2,deltts;
    real(4)    spdeltts;
    integer(4) ntsplit;
    real(8) sds,awk,enh;
    real(4) zpi;
    !here Must aligned 8 Byte
  end type implsch_gpar_type
  type(implsch_gpar_type) gpar
#endif

#ifdef USEINTERACTNEW
  real(8) lmdpd,lmdmd,lmdpd20,lmdmd20,wkal(4,12)
  integer::kjps(2,12,CJNTHET);
  integer(1) kjs(2,2,12,CJNTHET);
#endif
CONTAINS
#ifndef C_IMPLSCH
SUBROUTINE mean2(et,pvkdt,sds,awk) !{
  REALD et(kl,jnthet)
  type(implsch_pack_type) pvkdt
  real(8) sds,awk
  real(8) ae,asi,ark,ekspm
  integer k,j
  real(8) ekjs,ekjst
    ae=0.;   asi=0.;      awk=0.;   ark=0.
    !click (4*4*kl/4+4)*jnthet+=1392
    ! op 4*2*kl*jnthet 2688
    DO  k=1,kl-1 !{
      ekjs=0;
      DO  j=1,jnthet !{
        ekjs =ekjs +et(k,j)
      end do !}
      ae =ae +ekjs*dwk(k)
      awk=awk+ekjs*wkdk(k)
      asi=asi+ekjs*pvkdt%IWSDK(k)
      ark=ark+ekjs*wkibdk(k)
    end do !}
    ekjst=0
    DO  j=1,jnthet !{
      ekjst =ekjst +et(kl,j)
    end do !}
    DO  k=kl,kld !{
      ekjs=ekjst*wkh(k-kl+1)
      ae =ae +ekjs*dwk(k)
      awk=awk+ekjs*wkdk(k)
      asi=asi+ekjs*pvkdt%IWSDK(k)
      ark=ark+ekjs*wkibdk(k)
    end do !}
    if(ae>1.d-100)then
      awk=awk/ae;     asi=ae/asi;     ark=(ae/ark)**2
  !     sds=2.36e-5*asi*ark**3*ae**2/alpm2
  !         2.36e-5/alpm2=2.587605807
      ekspm=ae*ark*ark/0.0030162
      sds=ads*brkd1*asi/ark*sqrt(ekspm)*exp(-brkd2*0.64/ekspm)
    else
      sds=0
      awk=(wk(kl)+wk(1))/2
    endif
  return
end SUBROUTINE mean2 !}
SUBROUTINE mean3(et,depth) !{
  REALD et(kl,jnthet)
  real(8) depth
  real(8) ae,awk
  integer k,j
  real(8) ekjs,ekjst
  real(8) hs,hb,hbb,chbh
    ae=0.;    awk=0.;
    DO  k=1,kl-1 !{
      ekjs=0;
      DO  j=1,jnthet !{
        ekjs =ekjs +et(k,j)
      end do !}
      ae =ae +ekjs*dwk(k)
      awk=awk+ekjs*wkdk(k)
    end do !}
    ekjst=0
    DO  j=1,jnthet !{
      ekjst =ekjst +et(kl,j)
    end do !}
    DO  k=kl,kld !{
      ekjs=ekjst*wkh(k-kl+1)
      ae =ae +ekjs*dwk(k)
      awk=awk+ekjs*wkdk(k)
    end do !}
    if(ae>1e-30)then
      awk=awk/ae
      hs=4.*sqrt(ae)
      hb=zpi/awk*0.142*tanh(depth*awk)  !hb=zpi/awk*0.12*tanh(dep(ia,ic)*awk)/1.6726
      hbb =depth*(0.78125/1.6726)       !hbb=0.78125*dep(ia,ic)/1.5864792
      IF(hb>hbb) hb=hbb
      IF(hs>hb) THEN !{
        chbh=(hb/hs)**2
        ae=ae*chbh
        et=chbh*et
      end if !}
    end if
  return
end SUBROUTINE mean3 !}
integer function serkj(k,j)
  integer k,j
  serkj=(j-1)*kl+k
  if(k>kl)serkj=mkj+1
end function serkj
SUBROUTINE checkse(se,se1,mv,loc)
  REALD  se(kl,jnthet)
  REALD  se1(kl,jnthet)
  REALD  se2(kl,jnthet)
  real(8) mv
  integer loc(2)
  se2=abs(se1-se)
  mv=maxval(se2)
  loc=maxloc(se2)
end SUBROUTINE checkse
#ifdef USEINTERACTNEW
SUBROUTINE wwinteractionnew(ee,se,enh)
  REALD ee(kl*jnthet)
  REALD se(kl*jnthet)
  real(8) enh
  integer k,j,kj,jt,jtp
  real(8) eij,eij2,set,e0,e1,e2
  REALD eet(-(klp-kl-1):klp*jnthet)
  real(8) ecwks17(kl)
  do k=1,kl
    ecwks17(k)=enh*cwks17(k)
  enddo
  DO  k=-(klp-kl-1),0
    eet(k)=0
  enddo
  DO  j=1,jnthet
    jtp=(j-1)*klp;jt=(j-1)*kl
    DO  k=1,kl
      eet(jtp+k)=ee(jt+k)
    enddo
    DO  k=kl+1,klp
      eet(jtp+k)=0;
    enddo
  enddo
#define Calc_e1(id) e1=( wkal(1,id+1)*eet(kjps(1,id+1,j)+k)+wkal(2,id+1)*eet(kjps(1,id+1,j)+k+1)+wkal(3,id+1)*eet(kjps(2,id+1,j)+k  )+wkal(4,id+1)*eet(kjps(2,id+1,j)+k+1))
#define Calc_e2(id) e2=( wkal(1,id+2)*eet(kjps(1,id+2,j)+k)+wkal(2,id+2)*eet(kjps(1,id+2,j)+k+1)+wkal(3,id+2)*eet(kjps(2,id+2,j)+k  )+wkal(4,id+2)*eet(kjps(2,id+2,j)+k+1))
#define Calc(id,k1,k2,k3,k4,coff0,ect)         \
  jt=0                                        ;\
  DO  j=1,jnthet                              ;\
    DO  k=k1,k2-1                             ;\
      kj=k+jt;e0=(coff0)*ee(kj)               ;\
      Calc_e1(id) ;e2=0                       ;\
      se(kj)=se(kj)+(ect)*(e0*e1+e0*e2+e1*e2) ;\
    enddo                                     ;\
    DO  k=k2,k3-1                             ;\
      kj=k+jt;e0=(coff0)*ee(kj)               ;\
      Calc_e1(id) ;Calc_e2(id)                ;\
      se(kj)=se(kj)+(ect)*(e0*e1+e0*e2+e1*e2) ;\
    enddo                                     ;\
    DO  k=k3,k4                               ;\
      kj=k+jt;e0=(coff0)*ee(kj)               ;\
      Calc_e2(id) ;e1=0                       ;\
      se(kj)=se(kj)+(ect)*(e0*e1+e0*e2+e1*e2) ;\
    enddo                                     ;\
    jt=jt+kl                                  ;\
  END DO
  ! e1=lmdp     =1.25  e0=       1      e2=lmdm     =0.75  2,3  0    -4,-3
  ! e0=          1     e1=1/lmdp=0.8    e2=lmdm/lmdp=0.6   0   -3,-2 -6,-5
  ! e1=lmdp/lmdm=1.67  e2=1/lmdm=1.33   e0=          1     5,6  3,4  0
    Calc(0*2,1,3,kl-1  ,kl  ,   1.,  (-2.)*e0  ) !* resonate 1 koff: 2, 3 : -4,-3
    Calc(1*2,1,3,kl-1  ,kl  ,   1.,  (-2.)*e0  ) !* resonate 2 koff: 2, 3 : -4,-3
    Calc(2*2,2,5,kl    ,0   ,lmdpd,lmdpd20*e1  ) !* resonate 3 koff:-3,-2 : -6,-5
    Calc(3*2,2,5,kl    ,0   ,lmdpd,lmdpd20*e1  ) !* resonate 4 koff:-3,-2 : -6,-5
    Calc(4*2,1,1,kl-4  ,kl-3,lmdmd,lmdmd20*e2  ) !* resonate 5 koff: 5, 6 :  3, 4
    Calc(5*2,1,1,kl-4  ,kl-3,lmdmd,lmdmd20*e2  ) !* resonate 6 koff: 5, 6 :  3, 4
  jt=0
  DO  j=1,jnthet !{
    DO  k=1,kl !{
      kj=k+jt;
      se(kj)=se(kj)*enh*cwks17(k)
    END DO !}
    jt=jt+kl;
  END DO !}
end SUBROUTINE wwinteractionnew

#else 
!USEINTERACTNEW
!#define SERKJ
#ifdef SERKJ
SUBROUTINE wwinteraction(ee,se,dse,enh,ks)
  REALD ee(mkj)
  REALD se(mkj+4),dse(mkj+4)
  real(8) enh
  integer ks

  integer k,j,mr,kj,iid(8)
  integer im,im1,kp,kp1,j11,j12,j21,j22
  real(8) eij,eij2
  real(8) sap,sam,zua,ad,adp,adm,delad,deladp,deladm
  real(8) enht
  real(8) dv_t,dv_t1,dv_t2
      !*      1.2 angular loop
    DO  j=1,jnthet !{
        DO  mr=1,2 !{
          j11=jps(1,mr,j);        j12=jps(2,mr,j)
          j21=jps(3,mr,j);        j22=jps(4,mr,j)
          DO  k=1,ks !{
            kj=serkj(k,j)
            eij=ee(kj)
            IF (eij< 1.e-20) CYCLE
            kp =ikps(1,k);      kp1=ikps(2,k)
            im =ikps(3,k);      im1=ikps(4,k)
              eij2=eij*eij
              zua=2.*eij*alv(3) !### alv(3)=1/((1.+alamd)*(1.-alamd))

            !****************************************************************
            !*      nonlinear transfer
            iid(1)=kp +j11
            iid(2)=kp +j12
            iid(3)=kp1+j11
            iid(4)=kp1+j12
            iid(5)=im +j21
            iid(6)=im +j22
            iid(7)=im1+j21
            iid(8)=im1+j22
            if(kp<kl)then
                sap =wps(1)*ee(iid(1)) &
                    +wps(2)*ee(iid(2)) &
                    +wps(3)*ee(iid(3)) &
                    +wps(4)*ee(iid(4))
            else
                sap =(wpsk(1,k)*ee(serkj(kl,j11)) + wpsk(2,k)*ee(serkj(kl ,j12)))
            end if

            if(wmmask(k)/=0)then
              sam =(wps(5)*ee(iid(5)) &
                   +wps(6)*ee(iid(6))) &
                  +(wps(7)*ee(iid(7)) &
                   +wps(8)*ee(iid(8)))
            else
              sam=0
            endif
            dv_t1=sap*alv(1)+sam*alv(2)      !### alv(1)=1/(1.+alamd) alv(2)=1/(1.-alamd)
            dv_t2=-2.*sap*sam*alv(3)     !### alv(3)=alv(1)*alv(2)
            zua=2.*eij*alv(3)
            enht=enh*cwks17(k)
            ad    =(eij2*dv_t1   +dv_t2*eij)*enht
            delad =(eij *dv_t1*2.+dv_t2    )*enht
             se(kj  )= se(kj  )-2.0*ad
            dse(kj  )=dse(kj  )-2.0*delad
            adp   =ad*alv(4)  ! ## alv(4)=1/(1.+alamd)**3
            deladp=(eij2*alv(1)-zua*sam)*enht*alv(4)  !### alv(1)=1/(1.+alamd) alv(4)=1/(1.+alamd)**3
            if(kp<=kl)then
               se(iid(1))= se(iid(1))+   adp*wps( 1)
               se(iid(2))= se(iid(2))+   adp*wps( 2)
              dse(iid(1))=dse(iid(1))+deladp*wps( 9)
              dse(iid(2))=dse(iid(2))+deladp*wps(10)
              if(kp1<=kl)then
                 se(iid(3))= se(iid(3))+   adp*wps( 3)
                 se(iid(4))= se(iid(4))+   adp*wps( 4)
                dse(iid(3))=dse(iid(3))+deladp*wps(11)
                dse(iid(4))=dse(iid(4))+deladp*wps(12)
              endif
            endif
             if(wmmask(k)/=0)then
                adm=ad*alv(5) ! ## alv(5)=1/(1.-alamd)**3
                deladm=(eij2*alv(2)-zua*sap)*enht*alv(5)  !### alv(2)=1/(1.-alamd) alv(5)=1/(1.-alamd)**3
                 se(iid(5))= se(iid(5))+   adm*wps( 5)
                 se(iid(6))= se(iid(6))+   adm*wps( 6)
                 se(iid(7))= se(iid(7))+   adm*wps( 7)
                 se(iid(8))= se(iid(8))+   adm*wps( 8)
                dse(iid(5))=dse(iid(5))+deladm*wps(13)
                dse(iid(6))=dse(iid(6))+deladm*wps(14)
                dse(iid(7))=dse(iid(7))+deladm*wps(15)
                dse(iid(8))=dse(iid(8))+deladm*wps(16)
             end if
          end do !}
      end do !}
    end do !}
end SUBROUTINE wwinteraction
#else 
!SERKJ
SUBROUTINE wwinteraction(ee,se,dse,enh,ks)
  REALD ee(kl,jnthet)
  REALD se(kl,jnthet),dse(kl,jnthet)
  real(8) enh
  integer ks
  integer k,j,mr,kj
  integer im,im1,kp,kp1,j11,j12,j21,j22
  real(8) eij,eij2
  real(8) sap,sam,zua,ad,adp,adm,delad,deladp,deladm
  real(8) enht
  real(8) dv_t,dv_t1,dv_t2

    DO  k=1,ks !{
      kp =ikps(1,k);      kp1=ikps(2,k)
      im =ikps(3,k);      im1=ikps(4,k)
      !*      1.2 angular loop
      DO  j=1,jnthet !{
        eij=ee(k,j)
        IF (eij< 1.e-20) CYCLE
        eij2=eij*eij
        zua=2.*eij*alv(3) !### alv(3)=1/((1.+alamd)*(1.-alamd))
          DO  mr=1,2 !{
            j11=jps(1,mr,j);        j12=jps(2,mr,j)
            j21=jps(3,mr,j);        j22=jps(4,mr,j)
            !****************************************************************
            !*      nonlinear transfer
            if(kp<kl)then
                sap =wps(1)*ee(kp ,j11) + wps(2)*ee(kp ,j12) &
                    +wps(3)*ee(kp1,j11) + wps(4)*ee(kp1,j12)
            else
                sap =(wpsk(1,k)*ee(kl ,j11) + wpsk(2,k)*ee(kl ,j12))
            end if

            if(wmmask(k)/=0)then
              sam =(wps(5)*ee(im ,j21) &
                   +wps(6)*ee(im ,j22)) &
                  +(wps(7)*ee(im1,j21) &
                   +wps(8)*ee(im1,j22))
            else
              sam=0
            endif
            dv_t1=sap*alv(1)+sam*alv(2)      !### alv(1)=1/(1.+alamd) alv(2)=1/(1.-alamd)
            dv_t2=-2.*sap*sam*alv(3)     !### alv(3)=alv(1)*alv(2)
            zua=2.*eij*alv(3)
            enht=enh*cwks17(k)
            ad    =(eij2*dv_t1   +dv_t2*eij)*enht
            delad =(eij *dv_t1*2.+dv_t2    )*enht
            se (k  ,j  )= se(k  ,j  )-2.0*ad
            dse(k  ,j  )=dse(k  ,j  )-2.0*delad
            adp   =ad*alv(4)  ! ## alv(4)=1/(1.+alamd)**3
            deladp=(eij2*alv(1)-zua*sam)*enht*alv(4)  !### alv(1)=1/(1.+alamd) alv(4)=1/(1.+alamd)**3
            if(kp<=kl)then
              if(kp>kl.or.im>kl.or.im1>kl)then
                DBGO(0,*) 'Error1',kp,kp1,im,im1
                stop
              endif
               se(kp ,j11)= se(kp ,j11)+   adp*wps( 1)
               se(kp ,j12)= se(kp ,j12)+   adp*wps( 2)
              dse(kp ,j11)=dse(kp ,j11)+deladp*wps( 9)
              dse(kp ,j12)=dse(kp ,j12)+deladp*wps(10)
              if(kp1<=kl)then
                se (kp1,j11)= se(kp1,j11)+   adp*wps( 3)
                se (kp1,j12)= se(kp1,j12)+   adp*wps( 4)
                dse(kp1,j11)=dse(kp1,j11)+deladp*wps(11)
                dse(kp1,j12)=dse(kp1,j12)+deladp*wps(12)
              endif
            endif
             if(wmmask(k)/=0)then
                adm=ad*alv(5) ! ## alv(5)=1/(1.-alamd)**3
                deladm=(eij2*alv(2)-zua*sap)*enht*alv(5)  !### alv(2)=1/(1.-alamd) alv(5)=1/(1.-alamd)**3
                 se(im ,j21)= se(im ,j21)+   adm*wps( 5)
                 se(im ,j22)= se(im ,j22)+   adm*wps( 6)
                 se(im1,j21)= se(im1,j21)+   adm*wps( 7)
                 se(im1,j22)= se(im1,j22)+   adm*wps( 8)
                dse(im ,j21)=dse(im ,j21)+deladm*wps(13)
                dse(im ,j22)=dse(im ,j22)+deladm*wps(14)
                dse(im1,j21)=dse(im1,j21)+deladm*wps(15)
                dse(im1,j22)=dse(im1,j22)+deladm*wps(16)
             end if
          end do !}
      end do !}
    end do !}
end SUBROUTINE wwinteraction

#endif 
!SERKJ


#endif 
!USEINTERACTNEW
SUBROUTINE implsch_1b(ed,es,pvkdt,wvs,ucurc) !{
  ! NOTE :ed es maybe equal
  REALD ed(mkj)
  REALD es(mkj)
  real(4) wvs(4),ucurc(6)
  type(implsch_pack_type) pvkdt
  !Only Here
  REALD se(mkj+4),dse(mkj+4)
  integer k,j,ks
  integer im,im1,kp,kp1,j11,j12,j21,j22
  real(8) eij,eij2
  real(8) sap,sam,zua,ad,adp,adm,delad,deladp,deladm
  real(8) ww,scd,bett,betta,beta,wstarm
  real(8) sds,awk,ssds(kl+1)
  real(8) gadiag,deltee,dset,set
  real(8) enh,enht
  real(8) dv_t,dv_t1,dv_t2
  integer iid(8),kj
#if LOGSCURR==1
  real(8) sscu,sscut(3)
#endif
  REALD mv
  integer loc(2)
    call mean2(es,pvkdt,sds,awk)
    ww=wvs(3)**2
    IF (ww<0.05) ww=0.05
    dv_t1=cksp*gc2/ww
    dv_t=cksa*awk
    if(dv_t<dv_t1)dv_t=dv_t1
    ks=IndWk(dv_t*dwkmin)+1
    if(ks>kl)ks=kl
    se=0;   dse=0

    !*************************************
    !      call snonlin(e)
    !*************************************
    dv_t=0.75*pvkdt%depth*awk
    IF (dv_t<0.5) dv_t=0.5
    enh=1.+(5.5/dv_t)*(1.-0.833*dv_t)*exp(-1.25*dv_t)
    k=0;j=0
#ifdef USEINTERACTNEW
    call wwinteractionnew(es,se,enh)
#else
    call wwinteraction(es,se,dse,enh,ks)
#endif

    !Wind Input

    scd=sqrt((0.80+0.065*wvs(3))*0.001)
    wstarm=wvs(3)*scd
    bett =beta10*(1.+beta1*wvs(4))
    !breaken
    ssds=0
    DO  k=1,ks !{
      ssds(k)=-sds*wk(k)
    end do  !}
#if LOGSCURR==1
    if(logscurr/=0)then !{
            ! bcosths2(j)=1+cosths(j)*cosths(j)
            ! bsinths2(j)=1+sinths(j)*sinths(j)
            ! sincosths(j)=cosths(j)*sinths(j)
      ! !1:ux;2:uy;3:uxx;4:uxy;5:uyx;6:uyy
      sscut(1)=-acu*((bcosths2 (j))* ucurc(3) + (sincosths(j))*(ucurc(4)+ucurc(5)) + (bsinths2 (j))* ucurc(6) )
      sscut(2)=-acu*(-0.5)*(ucurc(3)+ucurc(6))  !!!??ZZ?
      sscut(3)=0 !+ucurc(2)*pvkdt%Rsd_tanLat  ! Big Circle
    end if !}
#endif

    DO  j=1,jnthet !{
      !sinput(e)
      betta=bett*28.*(wvs(1)*cosths(j)+wvs(2)*sinths(j))*scd
      DO  k=1,ks !{
        ! call sinput(e)
        beta=betta*wk(k)-bett*pvkdt%WS(k) ! 4 3
        if(beta<0)beta=0;
        !*************************************
        !      call source terms end
        !beta wind input
        !ssds:sdissip
        !ISSBO:sbottom
        !sscu:scurrent
        !cg*cosths(j)*Rsd_tanLat spherical coords
        ! ZZP cg=pvkdt%CCG(k)  ! Big Circle
        dset=beta+ssds(k)+pvkdt%SSBO(k)  !4 3 ! +pvkdt%CCG(k)*cosths(j)*pvkdt%Rsd_tanLat

#if LOGSCURR==1
        if(logscurr/=0)then !{
          dset=dset+sscut(1)*pvkdt%CCGD(k)+sscut(2) +sscut(3)
        end if !}
#endif
        kj=(j-1)*kl+k
        dv_t=es(kj)
        set= se(kj)+dset*dv_t !3 2
#ifndef USEINTERACTNEW
        dset=dse(kj)+dset
        gadiag=1.-deltts*0.5*dset
        gadiag=max(gadiag,1.)
        set=set/gadiag
#endif

#ifndef NODELTEE
        deltee=wstarm*grolim(k)
        if(set<-deltee)then !{
          set=-deltee;
        else if(set>deltee)then !}{
          set=deltee
        end if !}
#endif
        dv_t=dv_t+set*deltts ! 3 2
        if(dv_t<0)dv_t=0;
        ed(kj)=dv_t
      end do !}
      kj=(j-1)*kl+ks
      dv_t=es(kj)
      DO  k=1,kl-ks !{
        ed(kj+k)=dv_t*wkh(k+1)
      end do !}
    end do !}
    call mean3(ed,pvkdt%depth)
    return
end SUBROUTINE implsch_1b !}

SUBROUTINE implschs(ed,es,iacb,iace) !{
  REALD ed(kl,jnthet,0:*)
  REALD es(kl,jnthet,0:*)
  integer iacb,iace
  !Only Here
  integer iac
  !real(4) ::ucur1(6,10)  !1:ux;2:uy;3:uxx;4:uxy;5:uyx;6:uyy
  do iac=iacb,iace
    IF(nsp(iac)/=1) CYCLE
#if LOGSCURR==1
    if(logscurr/=0)then
      call implsch_1b(ed(:,:,iac),es(:,:,iac),ip(iac),wxy(:,iac),ucur(:,iac) )
    else
      call implsch_1b(ed(:,:,iac),es(:,:,iac),ip(iac),wxy(:,iac),ucur)
    end if
#else
    call implsch_1b(ed(:,:,iac),es(:,:,iac),ip(iac),wxy(:,iac),ucur)
#endif
  end do !}
  return
end SUBROUTINE implschs !}
#endif 
!C_IMPLSCH

#ifdef USEINTERACTNEW
SUBROUTINE nlweight !{
  real(8) lmdp,lmdm,alamd,cont,delphi1,delphi2,cl1,cl2,cl3,con,ch,wkk,cong
  real(8) scales(6),wkl(2,6),angls(6),wal(2,6)
  real(8) coffs(6)
  real(8) acl1,acl2,acl3
  integer j,k,j1,j2,k1,k2
  !real(8) wklp,wklm,wkp,wklap,wklap1,wklam,wklam1,wkm,wkm1
  !*      1. setting initial value
  !
  real(8) cgro
  integer kpst(6),l,i
  !kpst:  -6  -4  -3   2   3   5
  !kpst(ms):2 -4   2  -4  -3  -6 -3 -6 5 3 5 3

  !  integer ::ns(12)=[1,5,4,2,4,6,1,3,5,6,2,3]
  integer ms(12),ns(12)
  ! e1=lmdp     =1.25  e0=       1      e2=lmdm     =0.75  2,3  0    -4,-3
  ! e0=          1     e1=1/lmdp=0.8    e2=lmdm/lmdp=0.6   0   -3,-2 -6,-5
  ! e1=lmdp/lmdm=1.67  e2=1/lmdm=1.33   e0=          1     5,6  3,4  0
  data ms/4,2,4,2,3,1,3,1,6,5,6,5/
  data ns/1,5,4,2,4,6,1,3,5,6,2,3/
  !when pwk=1.21
  !data koff/2,-4,2,-4,-3,-6,-3,-6,3,5,3,5/
  !data koff/2,-4,2,-4,-3,-6,-3,-6,5,3,5,3/
  if(.not.allocated(cwks17))then
    ALLOCATE(cwks17(kl),grolim(kl))
  endif

  alamd=0.25
  con=7.862532087
  cong=con*sqrt(g)

  delphi1=11.48
  delphi2=33.56
  cgro=0.0000091*0.025 !*deltts

  DO  k=1,kl
    cwks17(k)=cong*wk(k)**8.5
    grolim(k)=cgro/wk(k)**4
  end do
  angls(1)=delphi1/(deltth*360/zpi)
  angls(2)=delphi2/(deltth*360/zpi)
  angls(3)=(delphi1+delphi2)/(deltth*360/zpi)
  angls(4)=-angls(1)
  angls(5)=-angls(2)
  angls(6)=-angls(3)

  lmdp=1+alamd;lmdpd=1/lmdp;lmdpd20=lmdpd**20
  lmdm=1-alamd;lmdmd=1/lmdm;lmdmd20=lmdmd**20

  scales(1)= (lmdm/lmdp)**2       !0.60^2=0.36   <1    1.21^-5.36   -6,-5
  scales(2)= (lmdm     )**2       !0.75^2=0.5625 >1    1.21^-3.02   -4,-3
  scales(3)= (1   /lmdp)**2       !0.80^2=0.64   >1    1.21^-2.34   -3,-2
  scales(4)= (lmdp     )**2       !1.25^2=1.5625 <1    1.21^ 2.34    2, 3
  scales(5)= (1   /lmdm)**2       !1.33^2=1.78   >1    1.21^ 3.02    3, 4
  scales(6)= (lmdp/lmdm)**2       !1.67^2=2.78   >1    1.21^ 5.36    5, 6

  coffs(1)=lmdmd;coffs(2)=lmdmd;coffs(3)=1;
  coffs(4)=lmdpd;coffs(5)=lmdpd;coffs(6)=1;

  do i=1,6
    wal(2,i) =angls(i)-floor(angls(i))
    wal(1,i) =1-wal(2,i)
  enddo
  do i=1,6
    kpst(i)= IndWk(scales(i)) !(-6:5)
    wkl(2,i)=(scales(i)-pwk**kpst(i))/((pwk-1)*pwk**kpst(i)) !!ZZ
    wkl(1,i)=1-wkl(2,i)
  enddo
  do l=1,12
    wkal(1,l)=coffs(ms(l))*wkl(2,ms(l))*wal(2,ns(l))
    wkal(2,l)=coffs(ms(l))*wkl(1,ms(l))*wal(2,ns(l))
    wkal(3,l)=coffs(ms(l))*wkl(2,ms(l))*wal(1,ns(l))
    wkal(4,l)=coffs(ms(l))*wkl(1,ms(l))*wal(1,ns(l))
  enddo
  do j=1,jnthet
    do l=1,12
        call GetJs( angls(ns(l)),j,j1,j2,jnthet)
        kjps(1,l,j)=(j1-1)*klp+kpst(ms(l))-1   !(-6:5)
        kjps(2,l,j)=(j2-1)*klp+kpst(ms(l))-1
        kjs(1,1,l,j)=kpst(ms(l))-1   !(-6:5)
        kjs(2,1,l,j)=j1-1
        kjs(1,2,l,j)=kpst(ms(l))-1
        kjs(2,2,l,j)=j2-1
    enddo
  enddo
#if 0
  if(1==0)then
    print'("kpst:",$)'
    do i=1,6
      print"(i4,$)",kpst(i);
    enddo
    print*
    do j=1,jnthet
      print'("kjs (1)",i3,":" ,$)',j
      do l=1,12
        print'(i4,$)',kjs(1,1,l,j)
        print'(i4,$)',kjs(1,2,l,j)
      enddo
      print*
      print'("kjs (2)",i3,":" ,$)',j
      do l=1,12
        print'(i4,$)',kjs(2,1,l,j)
        print'(i4,$)',kjs(2,2,l,j)
      enddo
      print*
      print'("kjsc(0)",i3,":" ,$)',j
      do l=1,12
        print'(i4,$)',kjs(1,1,l,j)+kjs(2,1,l,j)*klp
        print'(i4,$)',kjs(1,2,l,j)+kjs(2,2,l,j)*klp
      enddo
      print*
      print'("kjps(0)",i3,":" ,$)',j
      do l=1,12
        print'(i4,$)',kjps(1,l,j)
        print'(i4,$)',kjps(2,l,j)
      enddo
      print*
    enddo
    print*,klp
  endif
#endif
  !! R(kd)*Cl*sqrt(g)*k**8.5
  !  c0(e0*e0*(c1* e1+c2*e2)-2*c3*e0*e1*e2 )  :c3=c1*c2
  !     -2,lmdpd,lmdmd,lmdpd*lmdmd
  !     -2,lmdpd,lmdmd,lmdpd*lmdmd
  !  lmdpd20,lmdpd,lmdmd,lmdpd*lmdmd
  !  lmdpd20,lmdpd,lmdmd,lmdpd*lmdmd
  !  lmdmd20,lmdpd,lmdmd,lmdpd*lmdmd
  !  lmdmd20,lmdpd,lmdmd,lmdpd*lmdmd

  !  e(k),e(lp2*k1p),e(lm2*k2m)          00,41,25
  !  e(k),e(lp2*k1m),e(lm2*k2p)          00,44,22
  !  e(1/lp2*k1m),e(k),e(lm2/lp2*k3m)    34,00,16
  !  e(1/lp2*k1p),e(k),e(lm2/lp2*k3p)    31,00,13
  !  e(1/lm2*k2m),e(lp2/lm2*k3m),e(k)    55,66,00
  !  e(1/lm2*k2p),e(lp2/lm2*k3p),e(k)    52,63,00
  !
  !
  ! 12*4*jn*kl=48*kjn points
  ! 12*4*2*(jn+kl)=>! 12*4*4*(jn+kl)
  ! 12*4=48 weights
  ! 12*
  ! e(mn:=:L) !:L=1:12 m=1,6,n=1,6
  !   24w,96*(kl+jn)
  !    = wkl(1,m)*wal(1,n)*e(kps(1,m)+k+jps(1,n,j))
  !     +wkl(2,m)*wal(1,n)*e(kps(2,m)+k+jps(1,n,j))
  !     +wkl(2,m)*wal(2,n)*e(kps(2,m)+k+jps(2,n,j))
  !     +wkl(1,m)*wal(2,n)*e(kps(1,m)+k+jps(2,n,j))
  !   48w,192*(kl+jn)
  !    = wkal(1,L)*e(kps4(1,m)+k+jps(1,n,j))
  !     +wkal(2,L)*e(kps4(2,m)+k+jps(2,n,j))
  !     +wkal(3,L)*e(kps4(3,m)+k+jps(3,n,j))
  !     +wkal(4,L)*e(kps4(4,m)+k+jps(4,n,j))
  !   48w,192*kl*jn
  !    = wkal(1,L)*e(kjps4(1,L,j)+k)
  !     +wkal(2,L)*e(kjps4(2,L,j)+k)
  !     +wkal(3,L)*e(kjps4(3,L,j)+k)
  !     +wkal(4,L)*e(kjps4(4,L,j)+k)
  ! (w4 i4)*2=16  k1  =9: 1 1 3
  ! (w0 e4 o4 )*2

  !*      indices and weights f frquency grid
  return
end SUBROUTINE nlweight !}

#else 
!USEINTERACTNEW
SUBROUTINE nlweight !{
  real(8) alamd,cont,delphi,delphi2,cl1,cl2,con,ch,wkk
  real(8) acl1,acl2,cl12,cl21,cl22
  integer j,j1,j2,ikn,k
  real(8) wklp,wklm,wkp,wklap,wklap1,wklam,wklam1,wkm,wkm1,cong
  !*      1. setting initial value
  !
  real(8) cgro
  ALLOCATE(jps(4,2,jnthet),ikps(4,kl),wmmask(kl),wpsk(2,kl))
  if(.not.allocated(cwks17))then
    ALLOCATE(cwks17(kl),grolim(kl))
  endif
  alamd=0.25
  con=7.862532087
  cong=con*sqrt(g)

  delphi=-11.48
  delphi2=33.56
  cgro=0.0000091*0.025 !*deltts

#ifndef NODELTEE
  DO  k=1,kl !{
    grolim(k)=cgro/wk(k)**4
  end do !}
#endif
  !*      2.      computation of weights f angular grid
  cl1=delphi /(deltth*360/zpi)
  cl2=delphi2/(deltth*360/zpi)
  DO  j=1,jnthet !{
    call GetJs( cl1,j,jps(1,1,j),jps(2,1,j),jnthet)
    call GetJs( cl2,j,jps(3,1,j),jps(4,1,j),jnthet)
    call GetJs(-cl1,j,jps(1,2,j),jps(2,2,j),jnthet)
    call GetJs(-cl2,j,jps(3,2,j),jps(4,2,j),jnthet)
  end do !}
#ifdef SERKJ
  DO  j=1,jnthet !{
    jps(1,1,j)=(jps(1,1,j)-1)*kl
    jps(2,1,j)=(jps(2,1,j)-1)*kl
    jps(3,1,j)=(jps(3,1,j)-1)*kl
    jps(4,1,j)=(jps(4,1,j)-1)*kl
    jps(1,2,j)=(jps(1,2,j)-1)*kl
    jps(2,2,j)=(jps(2,2,j)-1)*kl
    jps(3,2,j)=(jps(3,2,j)-1)*kl
    jps(4,2,j)=(jps(4,2,j)-1)*kl
  end do !}
#endif
  acl1=abs(cl1-floor(cl1))
  acl2=abs(cl2-floor(cl2))
  alv(1)=1/(1.+alamd)
  alv(2)=1/(1.-alamd)
  alv(3)=alv(1)*alv(2)
  alv(4)=1/((1/alv(1))**3)
  alv(5)=1/((1/alv(2))**3)

  !*      indices and weights f frquency grid
  !高频插值系数
  ikn=IndWk(1/alv(1)**2)
  !important ZZ
  !(1/alv(1))**2-(1/alv(1)**2)=0.22204460492503131E-015 ,alv(1)=0.8

  wklap=((1/alv(1))**2-pwk**ikn)/((pwk-1)*pwk**ikn) !!ZZ
  wps(1)=(1.-wklap)*(1.-acl1)
  wps(2)=(1.-wklap)*acl1
  wps(3)=wklap *(1.-acl1)
  wps(4)=wklap *acl1
  wps( 9)=wps(1)**2
  wps(10)=wps(2)**2
  wps(11)=wps(3)**2
  wps(12)=wps(4)**2

  !低频插值系数
  ikn=IndWk((1/alv(2))**2)-1
  wklam=((1/alv(2))**2-pwk**ikn)/((pwk-1)*pwk**ikn)
  wps(5)=(1-wklam)*(1.-acl2)
  wps(6)=(1-wklam)*acl2
  wps(7)=   wklam *(1.-acl2)
  wps(8)=   wklam *acl2
  wps(13)=wps(5)**2
  wps(14)=wps(6)**2
  wps(15)=wps(7)**2
  wps(16)=wps(8)**2

  DO   k=1,kl !{
    wkk=wk(k)
    cwks17(k)=cong*sqrt(wkk**17)
    wklp=wkk/alv(1)**2
    wklm=wkk/alv(2)**2
    ikn=IndWk(1/alv(1)**2)
    ikn=k+ikn

    ikps(1,k)=ikn
    wkp=wk(ikn)
    ikps(2,k)=ikn+1

    IF(ikn>kl)then
      ikps(1,k)=kl+1
      ikps(2,k)=kl+1
    end if
    IF(ikn>=kl)THEN !{
      wpsk(1,k)=wps(1)*wkh(ikn-kl+1)+wps(3)*wkh(ikn-kl+1+1)
      wpsk(2,k)=wps(2)*wkh(ikn-kl+1)+wps(4)*wkh(ikn-kl+1+1)
    end if !}
    !低频插值系数
    IF(wklm<wk(1))THEN !{
      !超低频插值系数等于0
      wmmask(k)=0
      ikps(3,k)=1 ;     ikps(4,k)=1
    ELSE !} else {
      wmmask(k)=1
      ikn=IndWk((1/alv(2))**2)
      ikn=k+ikn-1
      IF(ikn<1) ikn=1
      ikps(3,k)=ikn
      ikps(4,k)=ikn+1
    end if !}
  end do !}
  return
end SUBROUTINE nlweight !}
#endif 
!USEINTERACTNEW
SUBROUTINE GetJs(ch,j,j1,j2,jnthet) !{
  real(8) ch
  integer j,j1,j2,jnthet
  j1=mod(j+floor(ch)+jnthet-1,jnthet)+1
  j2=mod((j1+1+jnthet-1),jnthet)+1
  return
end SUBROUTINE GetJs !}s

SUBROUTINE InitImplscha !{
  sbo=0.038*2./g
#ifdef C_IMPLSCH
  if(kl/=CKL)then
    DBGO0(0,*) 'CKl=',CKL,' not equal kl=',kl
    stop
  endif
#endif
  mkj=kl*jnthet
  call nlweight
#if SIGMALVL==2
  do ih=1,ndep
    dept=-vdep(ih)
    if(dept>0)then
      print*,"SIGMALVL==2 ,only can use fixed bv lvl,but vdep(",id,")=", vdep(ih),"<0,means use SIGMALVL"
    endif
    do k=1,kld
      gpar%expdep(k+CKL*(ih-1))=exp(2*dept*wk(k))*dwk(k)
    enddo
  end do
#endif
#ifdef C_IMPLSCH
  gpar%wk     =wk(1:kl)
  gpar%wkh    =wkh(1:kl)
  gpar%dwk    =dwk
  gpar%wkdk   =wkdk
  gpar%wkibdk =wkibdk
  gpar%cwks17 =cwks17
#ifndef NODELTEE
  gpar%grolim =grolim
#endif
  gpar%cosths =cosths(1:jnthet)
  gpar%sinths =sinths(1:jnthet)
  gpar%wkal   =wkal
  gpar%kjps   =kjps
  gpar%kjs    =kjs
  gpar%lmdpd  =lmdpd
  gpar%lmdmd  =lmdmd
  gpar%lmdpd20=lmdpd20
  gpar%lmdmd20=lmdmd20
  gpar%Dm2    =-2
  gpar%D1     =1
  gpar%deltts =deltts
  gpar%ntsplit=16
  gpar%spdeltts=deltts/gpar%ntsplit
  gpar%ads    =ads
  gpar%brkd1  =brkd1
  gpar%brkd2  =brkd2
  gpar%beta10 =beta10
  !gpar%beta1  =beta1
  gpar%beta11 =beta10*beta1
#endif
#ifdef C_IMPLSCH
  call C_InitImplsch(gpar,ip)
#endif
END SUBROUTINE InitImplscha !{

SUBROUTINE CalPRC(pvkdt,idep) !{
  type(implsch_pack_type),intent(out)::pvkdt
  integer ,intent(in)::idep
  integer k
  real(8) wkk,dk,tanhdk,wfk,wsk,dwkk,cg,ssbo
  DO  k=1,kld !{
    wkk=wk(k)
    dwkk=dwk(k)
    dk=idep*wkk
    if(dk<40)then
      wsk=sqrt(g*wkk*tanh(dk))
      cg=0.5*wsk*(1.+2.*dk/sinh(2.*dk))/wkk
      ssbo=-abo*sbo*wk(k)/sinh(2.*dk)
    else
      wsk=sqrt(g*wkk)
      ssbo=0
      cg=0.5*wsk/wkk
    endif

#ifdef OCHECK
    tanhdk=1.
    if(dk<4)tanhdk=tanh(dk)
    wsk=sqrt(g*wkk*tanh(dk))
    if (dk.gt.4.) then
      cg=0.5*wsk/wkk
    else
      if (dk.lt.0.14) then
        cg=sqrt(g*idep)
      else
        cg=0.5*wsk*(1.+2.*dk/sinh(2.*dk))/wkk
      endif
    endif
    ssbo=0
    if(dk<30)ssbo=-abo*sbo*wk(k)/sinh(2.*dk)
#endif
    pvkdt%SSBO(k)=ssbo
    pvkdt%WS(k)=wsk
    pvkdt%IWSDK(k)=(1/wsk)*dwkk
    !pvkdt%CCG(k)=cg
#if LOGSCURR==1
    pvkdt%CCGD(k)=cg*wkk/wsk
# endif
  end do !}
end SUBROUTINE CalPRC !}

end Module implsch_mod
SUBROUTINE InitImplschs(iacb,iace) !{
  use varcommon_mod
  !use partition_mod
  use implsch_mod
  IMPLICIT NONE
  integer iacb,iace
  integer iac,ierr,idep
  DO  iac=iacb,iace !{
    ip(iac)%Rsd_tanLat=Rsd_tanLat(iac)
    ip(iac)%depth=dep(iac)
    idep=dep(iac)
    call CalPRC(ip(iac),idep)
  end do !}
end SUBROUTINE InitImplschs !}
