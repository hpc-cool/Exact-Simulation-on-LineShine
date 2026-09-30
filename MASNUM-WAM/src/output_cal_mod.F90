#include "wavedeff.h"
Module output_cal_mod
  use varcommon_mod
  IMPLICIT NONE
  !=======================
  ! private
  ! public::InitOutPutCal,mean1,mean1_1,ACCUMEA,mixture
  real(8),pointer::pebdep(:,:,:) !计算优化 深度相关
  real(8),pointer::pvkdo(:,:,:) !计算优化 深度相关

  real(8)::tztp=1.2,tztz=1.099314
#define  IWFDK   1
#define  IIWF2DK 2
#define  IWS2DK  3
#define  IWF     4
#define  BVIWS2  1
integer ::iact,OutPutCalInited=0;
CONTAINS
  SUBROUTINE DeinitOutPutCal
    if(OutPutCalInited/=0)then
    OutPutCalInited=0;
    endif
  end SUBROUTINE DeinitOutPutCal

  !-----------------------------------------------
  !计算波高周期波向等
  !-----------------------------------------------
#ifndef C_IMPLSCH
  SUBROUTINE mean1t(et,iac,apet,tpft,aett,h1_3t)
    REALD,intent(in ) :: et(kl,jnthet)
    integer,intent(in)  :: iac
    real(8),intent(out) :: apet,tpft,aett,h1_3t
    iact=iac
    call mean1_1(et,pvkdo(:,:,iac),apet,tpft,aett,h1_3t)
  end SUBROUTINE mean1t
  SUBROUTINE mean1_1(et,pvkdot,apet,tpft,aett,h1_3t)
    REALD,intent(in ) :: et(kl,jnthet)
    real(8)::pvkdot(kld,IMDPSO)
    real(8),intent(out) :: apet,tpft,aett,h1_3t
    INTEGER k,j
    real(8) aess,awfss,asiss,apess,aets,aetc
    real(8) wkdkk
    real(8) ekj,ekjs,ekjst,aetst,aetct
    !aet  :波向
    !tpf  :谱峰周期
    !h1_3 :波高
    !ape  :跨零中期
    real(8) wfk ,wfk1,wsk ,wsk1,wkk ,wkk1
    aess=0.;asiss=0.;awfss=0.;apess=0.; aets=0.0; aetc=0.0
    DO  k=1,kl !{
      ekjs=0;aetst=0.;aetct =0.;
      wfk =pvkdot(k  ,IWF)
      wfk1=pvkdot(k+1,IWF)
      wsk =zpi*wfk
      wsk1=zpi*wfk1
      wkk =wk(k)
      wkk1=wk(k+1)
      DO  j=1,jnthet !{
        ekj =et(k,j)
        ekjs =ekjs +ekj
        aetst=aetst+ekj*sinths(j)
        aetct=aetct+ekj*cosths(j)
      end do !}
      wkdkk=wkdk(k)
      aess =aess +ekjs  *dwk(k)
      apess=apess+ekjs  *pvkdot(k,IWS2DK )
      asiss=asiss+ekjs  *pvkdot(k,IIWF2DK)
      awfss=awfss+ekjs  *pvkdot(k,IWFDK )
      aets =aets +aetst *wkdkk
      aetc =aetc +aetct *wkdkk
    end do !}

    ekjst=0;aetst=0.;aetct =0.
    DO  j=1,jnthet !{
      ekj =et(kl,j)
      ekjst =ekjst +ekj
      aetst=aetst+ekj*sinths(j)
      aetct=aetct+ekj*cosths(j)
    end do !}
    DO  k=kl+1,kld !{
      ekjs=ekjst*wkh(k-kl+1)
      wkdkk=wkdk(k)*wkh(k-kl+1)
      aess =aess +ekjs  * dwk(k)
      apess=apess+ekjs  * pvkdot(k,IWS2DK )
      asiss=asiss+ekjs  * pvkdot(k,IIWF2DK)
      awfss=awfss+ekjs  * pvkdot(k,IWFDK )
      aets =aets +aetst * wkdkk
      aetc =aetc +aetct * wkdkk
    end do !}
    !^^^^^^^^^^^^^^^^^^
    !ape: tz; tpf: tp; aet: th; h1_3:hs
    aett=atan2d(aets,aetc)
    IF (aett<0.)aett=360.+aett
    h1_3t=4.*sqrt(aess)

    if(aess>1e-30)then
      apet=tztz*zpi*sqrt(aess/apess)
      tpft=asiss*awfss/aess**2
    else
      apet=0
      tpft=0
    endif
  end SUBROUTINE mean1_1
  SUBROUTINE mean1(ea,iacb,iace,ape,tpf,aet,h1_3)
    REALD,intent(in) :: ea(kl,jnthet,0:*)
    INTEGER,intent(in) :: iacb,iace
    real(4)::aet(0:*),tpf(0:*),h1_3(0:*),ape(0:*)
    !aet  :波向
    !tpf  :谱峰周期
    !h1_3 :波高
    !ape  :跨零中期
    real(8) ::apet,tpft,aett,h1_3t
    INTEGER iac
    do iac=iacb,iace
      iact=iac
      call mean1_1(ea(:,:,iac),pvkdo(:,:,iac),apet,tpft,aett,h1_3t)
      ape(iac)=apet; tpf(iac)=tpft;
      aet(iac)=aett;h1_3(iac)=h1_3t
    end do !}
    RETURN
  END SUBROUTINE mean1 !}
  subroutine mixture_1(et,pvkdot,pebdept,bvo)
    REALD,intent(in) ::et(kl,jnthet)
    real(8)::pvkdot(kld,IMDPSO),pebdept(kld,BVIMDPS+ndep)
    real(4),intent(out)::bvo(ndep)
    real(8) ebdep,ekj,ekjs
    real(8) wkk,wsk2
    real(8) bv1(ndep),bv2(ndep),bv3(ndep),bvt
    INTEGER kh,k,j
    bv1=0;bv2=0;bv3=0

    do k=1,kl-1 !{
      ekjs=0;
      do j=1,jnthet !{
        ekjs=ekjs+et(k,j)
      end do !}
      wkk=wk(k)
      wsk2=pvkdot(k,BVIWS2)

      do kh=1,ndep!{
        ebdep=pebdept(k,BVIMDPS+kh)
        bvt=ekjs*ebdep  ;bv1(kh)=bv1(kh)+bvt
        bvt=bvt*wsk2    ;bv2(kh)=bv2(kh)+bvt
        bvt=bvt*wkk     ;bv3(kh)=bv3(kh)+bvt
      end do !}
    end do !}

    ekjs=0;
    do j=1,jnthet !{
      ekjs=ekjs+et(kl,j)
    end do !}
    do k=kl,kld !{
      wsk2=pvkdot(k,BVIWS2)
      wkk=wk(k)
      do kh=1,ndep!{
        ebdep=pebdept(k,BVIMDPS+kh)*wkh(k-kl+1) !ZZZZZZZZZZZZZZZZZ ERROR here
        bvt=ekjs*ebdep  ;bv1(kh)=bv1(kh)+bvt
        bvt=bvt*wsk2    ;bv2(kh)=bv2(kh)+bvt
        bvt=bvt*wkk     ;bv3(kh)=bv3(kh)+bvt
      end do !}
    end do !}

    do kh=1,ndep!{
      if(bv2(kh)>1e-30)then
        bvo(kh)=bv1(kh)/sqrt(bv2(kh))*bv3(kh)
      else
        bvo(kh)=0
      endif
    end do !}
  end subroutine mixture_1
  subroutine mixture(ea,iacb,iace,bv) !{
    REALD ,intent(in)::ea(kl,jnthet,0:*)
    INTEGER ,intent(in)::iacb,iace
    real(4)::bv(0:nwpc,ndep)
    INTEGER iac
    real (4) bvt(ndep)
    do iac=iacb,iace
      call mixture_1(ea(:,:,iac),pvkdo(:,:,iac),pebdep(:,:,iac),bvt)
      bv(iac,:)=bvt
    end do !}
    RETURN
  END subroutine mixture !}

  subroutine AverageEA(ea,nea,iacb,iace)
    REALD:: ea(kl,jnthet,0:nwpc)
    INTEGER nea,iacb,iace
    INTEGER iac
    do iac=iacb,iace
      ea(:,:,iac)=ea(:,:,iac)/nea
    enddo
  end subroutine AverageEA
  subroutine ZeroEA(ea,iacb,iace)
    REALD :: ea(kl,jnthet,0:nwpc)
    INTEGER iacb,iace
    INTEGER iac
    do iac=iacb,iace
      ea(:,:,iac)=0.
    enddo
  end subroutine ZeroEA

  subroutine addea(ea,ee,iacb,iace)
    REALD:: ea(kl,jnthet,0:nwpc)
    REALD:: ee(kl,jnthet,0:nwpc)
    INTEGER iacb,iace
    INTEGER iac
    do iac=iacb,iace
      ea(:,:,iac)=ea(:,:,iac)+ee(:,:,iac)
    enddo
  end subroutine addea

#endif

END Module output_cal_mod
SUBROUTINE InitOutPutCal(iacb,iace,dep)
  use varcommon_mod
  use output_cal_mod
  IMPLICIT NONE
  integer iacb,iace
  real(4) dep(*)
  integer k,iac,ih,ierr,ko
  type (AVHIST),pointer::ph
  real(8)dept,depv
  real(8) wkk,dk,wfk,wsk,dwkk
  do iac=iacb,iace
    depv=dep(iac)
    DO  k=1,kld !{
      wkk=wk(k)
      dwkk=dwk(k)
      dk=depv*wkk
      if(dk<40)then
        wsk=sqrt(g*wkk*tanh(dk))
      else
        wsk=sqrt(g*wkk)
      endif
      wfk=wsk/zpi
      pvkdo  (k,IWF    ,iac)=wfk
      pvkdo  (k,IWFDK  ,iac)=wfk*dwkk
      pvkdo  (k,IIWF2DK,iac)=1./(wfk**2)*dwkk
      pvkdo  (k,IWS2DK ,iac)=wsk**2*dwkk
      pebdep(k,BVIWS2 ,iac)=wsk**2
    enddo
    do ih=1,ndep
      if(vdep(ih)<0)then
        dept=depv*vdep(ih)
      else
        dept=-vdep(ih)
      end if
      do k=1,kld
        pebdep(k,BVIMDPS+ih,iac)=exp(2*dept*wk(k))*dwk(k)
      enddo
    end do
  end do
#ifdef C_IMPLSCH
  call initc_output_cal(nwpc,ndep,pvkdo,pebdep);
#endif
  OutPutCalInited=1;
end SUBROUTINE InitOutPutCal
#ifndef MTHREAD
subroutine accumea(ind,ee,nwpb,nwpe)
  use varcommon_mod
  IMPLICIT NONE
  REALD ,intent(in):: ee(kl,jnthet,0:nwpc)
  integer,intent(in)::nwpb,nwpe,ind
  integer ko,iac
  do ko=1,Navhists
    if(avhists(ko)%mean/=0)then
      do iac=nwpb,nwpe
        avhists(ko)%ea(:,:,iac)=avhists(ko)%ea(:,:,iac)+ee(:,:,iac);
      enddo
      if(nwpb==1)avhists(ko)%nea=avhists(ko)%nea+1
    endif
  enddo
end subroutine accumea
subroutine CheckOutPut(igrp,ind,RFB,ee,iacb,iace,iacb2,iace2,nout)
  use varcommon_mod
  use output_cal_mod
  IMPLICIT NONE
  integer,intent(in):: iacb,iace ,nout
  integer,intent(in):: igrp,ind,RFB,iacb2,iace2
  REALD ,target,intent(in):: ee(kl,jnthet,0:nwpc)
  integer ko,ahist,newfile,ntpf
  integer noutc
  real(8) ::timet
  REALD ,pointer::ea(:,:,:)
  type (AVHIST),pointer::ph
  noutc=nout
  do ko=1,Navhists
    if(avhists(ko)%hist_ect/=0)then
      ph=>avhists(ko)
      if(ph%mean/=0)then
        ea=>ph%ea
        call AverageEA(ph%ea,ph%nea,iacb,iace)
      else
        ea=>ee
      endif
      IF(ph%wrwa>0)  call mean1(ea,iacb,iace,ph%ape,ph%tpf,ph%aet,ph%h1_3) 
      IF(ph%wrbv>0)  call mixture(ea,iacb,iace,ph%bv)
      if(ph%mean/=0) call ZeroEA(ea,iacb,iace) !clear
      call HistOutPut(ph)
      ph%nea=0; ph%hist_ect=0
      noutc=noutc-1
      if(noutc<=0)exit
    endif
  enddo
end subroutine CheckOutPut
#endif
