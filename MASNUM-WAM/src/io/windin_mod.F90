#include "wavedeff.h"
Module windin_mod
  use varcommon_mod
  use wav_cal_mod
  use DataIO_mod
#ifdef COUPLE_CCSM
  use msg_mod
#endif
  IMPLICIT NONE
  real(8)::etime0,etime1,rdtime
  integer ::constwind=0
  !风数据 wind Data used by windinput
  real(4) ,POINTER::wind0x(:),wind0y(:),windx(:),windy(:),wvt(:)
  real(4) ,POINTER::dwindx(:),dwindy(:)
  integer ::firstcall_InitDataIn=1;
CONTAINS
  SUBROUTINE InitDataIn !{
    integer ierr
    if(.not.associated(wxy))then
      ALLOCATE(wxy(4,0:nwpc),STAT=ierr)
    endif
    ALLOCATE(wxyo(4,0:nwpc),STAT=ierr)
    if(.not.associated(wind0x))then
      ALLOCATE(wind0x(0:nwpc),wind0y(0:nwpc),windx(0:nwpc),windy(0:nwpc),wvt(0:nwpc),STAT=ierr)
      ALLOCATE(dwindx(0:nwpc),dwindy(0:nwpc))
    endif
    if(.not.associated(ucur))then
#if LOGSCURR==1
      if(logscurr/=0)then
        ALLOCATE(ucur(6,0:nwpc),STAT=ierr)
      else
        ALLOCATE(ucur(6,1),STAT=ierr) 
      end if
#else
      ALLOCATE(ucur(6,1),STAT=ierr) 
#endif
    endif
    wxy=0;
    wind0x=0;wind0y=0;windx=0;windy=0;wvt=0;
    etime0=0;etime1=0;
    ucur=0;
  end SUBROUTINE InitDataIn
!======================================
! mode:
!   0: Init MPI msg etc
!   1: Init before startup
!   2: Get First Data for cold start
!   3: Get First Data for warm start
!   4: Get Data For for timeStep
!======================================
integer function IDataIO(mode) !{
  integer mode
  real(8) wt,cc,cosd,sind
  real(8) sh
  integer idate,i,j ,key,iac
  real*8 etimet,etimep,dtm
  integer acwv;
  IDataIO=0
  IF(mode == 4)then 
    runstate%timecur=runstate%timecur+delltday 
    wxyo=wxy;
    if(constwind>0)then 
      IDataIO=1
      return ;
    endif
    IF(runstate%edaycur+runstate%timecur>etime1+1e-9)then
        wind0x=windx ;wind0y=windy;
        etime0=etime1
        if(IRDataIO(4,runstate%edaycur+runstate%timecur,windx,windy,etime1,etimep)>0)then
          rdtime=1./(etime1-etime0);
          dwindx=windx-wind0x
          dwindy=windy-wind0y
        else
          etime1=etime0; rdtime=0; dwindx=0.; dwindy=0.
        endif
    end if
    dtm=(runstate%edaycur+runstate%timecur-etime0)*rdtime
    if(dtm>1.e-6)then
      wxy(1,:)=wind0x+dtm*dwindx;
      wxy(2,:)=wind0y+dtm*dwindy;
    else
      wxy(1,:)=wind0x
      wxy(2,:)=wind0y
    endif
    wvt=sqrt(wxy(1,:)**2+wxy(2,:)**2)
    wxy(4,:)=wvt-wxy(3,:)
    wxy(3,:)=wvt
  else IF(mode == 0)then 
    constwind=0;
    ! here only First para mode is used
    key=IRDataIO(0,0.D0,windx,windy,etime1,etimep) 
  else IF(mode == 1)then
    call InitDataIn
    acwv=(abs(constwindx)+abs(constwindy))*100
    if(acwv>=10)then
      constwind=2
    else if(acwv>1)then
      constwind=1
    else 
      constwind=0 
    endif
    if(constwind<=1)key=IRDataIO(1,runstate%edaycur+runstate%timecur,windx,windy,etime1,etimep) 
  else IF(mode == -1)then 
    if(constwind<=1)key=IRDataIO(-1,runstate%edaycur+runstate%timecur,windx,windy,etime1,etimep) 
  else 
    if(constwind<=1)then
      if(IRDataIO(2,runstate%edaycur+runstate%timecur,windx,windy,etime1,etimep)<=0)return
      etime0=etimep;
      if(IRDataIO(2,etimep-1e-6,wind0x,wind0y,etime0,etimep)<=0)then
        etime0=etime1; wind0x=windx;wind0y=windy
      endif
      if(etime1>etime0+1e-4)then
        dwindx=windx-wind0x
        dwindy=windy-wind0y
        rdtime=1./(etime1-etime0);
        dtm=(runstate%edaycur+runstate%timecur-etime0)*rdtime
        wxy(1,:)=wind0x+dtm*dwindx;
        wxy(2,:)=wind0y+dtm*dwindy;
      else
        rdtime=0.;dwindx=0; dwindy=0;
        wxy(1,:)=wind0x;
        wxy(2,:)=wind0y;
      endif
      wvt=sqrt(wxy(1,:)**2+wxy(2,:)**2)
      wxy(4,:)=wvt-wxy(3,:)
      wxy(3,:)=wvt
    else
        wxy(1,:)=constwindx;
        wxy(2,:)=constwindy;
    endif
    IDataIO=1 
  end if
  IDataIO=1
END function IDataIO

End Module windin_mod
