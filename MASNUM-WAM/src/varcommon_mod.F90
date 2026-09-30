#include "wavedeff.h"
#ifndef KLPAT
#define KLPAT 8
#endif
#ifndef MAXAVHIST
#define MAXAVHIST 16
#endif
Module varcommon_mod
#ifndef NO_MPI
    use wav_mpi_mod
#endif
    use debughlp_mod
    use netcdf_mod
    use wav_cal_mod
    IMPLICIT NONE
public
    integer(4) :: kl=CKL,klp=CKL+KLPAT,kld=CKL,jnthet=CJNTHET 
    real(8) :: pwk=1.21,wkmin=0.0071
    integer(4) :: gridtype=1 ! 0:均匀经纬度网格;1:不均匀经纬度网格；2正交网格
    !均匀经纬度网格参数
    integer(4) :: Gcircle=1,TESTDEP=1 ! 全球模式,环绕'
    integer(4) :: gixl,giyl ! 计算区域网格数
    integer(4) :: nwpa ! all Water Points Array number
    integer(4) :: nwpc ! Inner Water Points Array number
    integer(4) :: nwps ! Need Send Water Points Array number
    !integer ixd(6),iyd(6)  !分区网格定义
    integer(4) :: rectc(4)
    integer ::NCorePClu=38,MCorePClu=37,NCluPNode=16,NManageCore=0
    integer ::NGrpPProc=16 ,NThPGrp=38 ,NProcPNode=1,ManageCoreId=-1
    integer ::OffClu=0, SkipClu=1, OffCore=0, SkipCore=1
    !16*4
    real(8) :: grdszx=2,grdszy=2 !网格大小
    real(8) :: gxlon0=-1.,gylat0=14.,gxlon1=122.,gylat1=27. !(gxlon0-gxlon1,gylat0-gylat1):计算区域
    ! Constant
    !g=9.78049(1 + 0.0052884*Sin(y(ic)**2 - 0.0000059 *Sin(2*y(ic)) ** 2) - 0.00000286
    !g=9.80665 when y(ic)=45
    !real(8):: g=9.80665,gc2=0.877**2*g,pi=3.1415926535897932384626433832795d0,zpi=2*pi
    real(8) :: g=9.81,pi=3.1415926535897932384626433832795d0
    real(8) :: gc2,zpi
    !地球参数
    !deg2m 111194.4
    real(8) :: Rs=6367451.637,deg2m
    !real(8):: Rs=6371300,deg2m=111200.4
    real(8) :: alog10pwk,dlog10pwk,dwkmin,wkmax
    real(8) :: deltth ,ddeltth !角度划分
    real(8) :: windfield=50    !风区km
    real(8) ::beta0=1.,beta1=0.,acu=1.,brkd1=0.000132,brkd2=2.61,beta10
    real(8) :: deltts=1800.,delltday   ! deltts:时间步长(秒)      

    real(8) :: RunDays=10
    real   (8) smooth_p0
    real(8) :: timeInit
    integer :: dayinit
    !16*4 + 29*8
    integer(4) :: StartDate,StartTime
    integer(4) :: logscurr =0!THE SWITCH FOR WAVE-CURRENT INTERACTION
    integer(4) :: logsdiss =0!THE SWITCH FOR Sds,1-Yuan''s,0-Komen''s
    !地形
    !DEPTHMOD: 0  use Input DEPTH
    !          1  Inter Input DEPTH
    !      10 flat Depth
    integer(4) :: DEPTHMOD=0
    integer(4) :: Logau=0,Logbv=2,LogWind=2,Loghs=2,Logtp=2,Logtz=2,Logth=2
    integer(4) :: ndep=14,sigmalvl=0, itt__
    !16*4 + 29*8 +14*4
    real(4)::vdep(160)=(/  0., 10., 20., 30., 50., 75., 100.,125.,150.,200., &
                                   250.,300.,400.,500.,(0.,itt__=15,160) /)
    !16*4 + 29*8 + 14*4 + 160*4
    ! Test Para
    real(4) :: plo(16)=1000,pla(16)=1000
    !16*4 + 29*8 + 14*4 + 160*4+2*4
    integer(4) :: nmonitor,gnmonitor
    integer(4) :: iap(16),icp(16),iacp(16)
    integer(4) :: iCheckRestart=0
    integer(4) :: outputNEP=4,OutType=0
    integer(4) :: noOutPut=0
    real   (4) :: constwindx=0,constwindy=0   
    integer(4) :: info_dbug
    integer :: ncpl = 4,NSCPL,Delttscpl
    type TRunState
      integer(4) :: state,stop_now,stop_eod,rest_now,rest_eod,rest_eot,hist_now,hist_eod,hist_eot;
      integer(4) :: nTimeStep
      integer ::NextOutStep=-1;
      integer ::ONextOutStep=-1;
      integer(4) :: cdatecur,ctimecur,edaycur,ndays;
      real(8)    :: timecur;
      integer(4) :: iPreCalT,mTimeStep;
    end type TRunState
    type(TRunState),target :: runstate_t        ! avoid Error TempVar when runstate Not Inited
    type(TRunState),pointer:: runstate

  real(4),allocatable::xcord(:,:),ycord(:,:)       ! 网格坐标(度)
  real(8),allocatable::sxcord(:),sycord(:) !nwpa
  integer(1),allocatable:: nsp(:) ! nwpc   !格点类型
  real(4),pointer:: depl(:)     ! nwpa   !水深
  real(4),pointer:: dep(:)     ! nwpa   !水深
  !谱空间定义参数
  real(8) ,allocatable::wk(:),thet(:) !谱空间坐标

  integer,allocatable::nspg(:,:)
  real(4) ,allocatable::depg(:,:)
  real(8),allocatable::dxm(:),dym(:) !网格大小(米)
  real(8),allocatable::dddx(:),dddy(:) !地形梯度
  real(4),pointer::Rsd_tanLat(:) ! 计算优化 tand(y)/rs
  !real(8),allocatable::Cosths_Rsd_tanLat(:,:) ! 计算优化 (cosd(thet)*Rsd_tanLat
  integer,parameter:: imaxwfd=1000 !最大影响深度

  real(8) ,allocatable::wkh(:),dwk(:),wkdk(:),wkibdk(:),sinths(:),cosths(:),bcosths2(:),bsinths2(:),cosths2(:),sinths2(:),sincosths(:)
  !风数据 wind Data used by Model
  real(4) ,pointer::wxy(:,:)
  real(4) ,pointer::wxyo(:,:)
  !流数据 Current Data
  real(4) ,pointer::ucur(:,:)  !1:ux;2:uy;3:uxx;4:uxy;5:uyx;6:uyy

  REALD ,pointer::eet(:,:,:),eec(:,:,:)
  character*(256) :: startType='startup',caseid='test'
  character*(256) :: pathinit='init/',pathrest='rest/',pathhist='hist/',pathLoc='/tmp/'
  character*(256) :: pathwind='wind/',pathcur='cur/'
  character*(256) :: topofn='init/topography.dat'
  character*(256) :: ufnform='%Pwind%4.4Y%2.2M.nc'
  character*(256) :: vfnform='%Pwind%4.4Y%2.2M.nc'
  character*(16 ) :: uwind_name='windu',vwind_name='windv'
  !character*(256):: ufnform='%Puwnd.sig995.%4.4Y.nc'
  !character*(256):: vfnform='%Pvwnd.sig995.%4.4Y.nc'
  !character*(16 ):: uwind_name='uwnd',vwind_name='vwnd'
  real*8          :: WTIME_PATCH=0
  integer         :: topofType=1
  character(8)    :: REST_OPTION="NMONTH"
  INTEGER         :: REST_N=1,REST_TYPE
  real(8)::RestTime=0
  type AVHIST
      INTEGER ::NN=1,itype=1000,NTPF=1,mean=1,hist_ect=0
      real(8)::nextTime=0
      INTEGER ::oldID=-100000
      INTEGER ::nea,itpf=1000000 !,HISTIDS(6)=-1
      INTEGER ::istime=-1
      INTEGER ::wrwa=0,ncidwa=-1
      INTEGER ::wrbv=0,ncidbv=-1
      INTEGER ::Logbv=2,LogWind=2,Loghs=2,Logtp=2,Logtz=2,Logth=2,Logau=0
      integer*8::waoff=-1,bvoff=-1
      character*256 desc,ftimestr,wafn,bvfn
      character(8)::OPTION='NONE',TAG='ah'
      REALD ,pointer::ea(:,:,:)
      real(4),pointer::aet(:),tpf(:),h1_3(:),ape(:),bv(:,:)
  end type AVHIST

  INTEGER,parameter:: IMAXAVHIST = MAXAVHIST 
  INTEGER Navhists
  type (AVHIST),allocatable,target::avhists(:)
  INTEGER,private::OutputAllPara=1
  INTEGER,private::ctlOpened=0,ctlfid=8

  type propinf_type
    !real(4) dep
    integer i8(8)
  end type propinf_type
  type (propinf_type),pointer::ipos8(:)
  integer(4),pointer::ipos12(:,:)
  real(8),allocatable::nfpcs(:)
  !for output,moved from out_cal_mod.F90
  !real(4),pointer,public::aet(:),tpf(:),h1_3(:),ape(:),bv(:,:)
!===================================================
!----                      THE         END                     -------C
CONTAINS
  SUBROUTINE SkipLines(ifd) !{
      integer ifd
      character line*100
      read(ifd,'(a)')line
      do while(line(1:1) /='@')
          read(ifd,'(a)')line
      end do
  end SUBROUTINE SkipLines !}

  SUBROUTINE InitPara !{
    real(8) ::PreCalTime=0
    character*(8)::AVHIST_OPTIONS(IMAXAVHIST)='NNNN',AVHIST_OTAGS(IMAXAVHIST)='ah'
    INTEGER ::AVHIST_NS(IMAXAVHIST)=24,AVHIST_NTPFS(IMAXAVHIST)=-1,AVHIST_TYPES(IMAXAVHIST)
    INTEGER ::LogauS(IMAXAVHIST)=0,LogbvS(IMAXAVHIST)=0,LogWindS(IMAXAVHIST)=0,LoghsS(IMAXAVHIST)=0,LogtpS(IMAXAVHIST)=0,LogtzS(IMAXAVHIST)=0,LogthS(IMAXAVHIST)=0
    INTEGER ::MinvS(IMAXAVHIST)
    NAMELIST /ParaCtl/OutputAllPara
    NAMELIST /runDef/startType,caseid,pathinit,pathrest,pathhist,pathLoc,PATHCUR,topofn,topofType,DEPTHMOD,DbgLvl,StartDate,StartTime,PreCalTime,RunDays,outputNEP,OutType
    NAMELIST /WindForm/PATHWIND,ufnform,uwind_name,vfnform,vwind_name,WTIME_PATCH
    NAMELIST /GPARA/deltts,gxlon0,gylat0,gxlon1,gylat1,grdszx,grdszy,Gcircle,windfield,logscurr,logsdiss
    NAMELIST /BVDepDef/ndep,vdep
    NAMELIST /OUTDEF/REST_OPTION,REST_N,AVHIST_OTAGS,AVHIST_OPTIONS,AVHIST_NS,AVHIST_NTPFS,LogbvS,LogWindS,LoghsS,LogtpS,LogtzS,LogthS,LogauS
    NAMELIST /implschPara/beta0,beta1,acu,brkd1,brkd2
    NAMELIST /PKPARA/pwk,wkmin,KL,kld,jnthet
    NAMELIST /OUTPUTPnt/plo,pla
    NAMELIST /testOpt/iCheckRestart,noOutPut,constwindx,constwindy
    NAMELIST /CPLDEF/NCPL
    NAMELIST /THEADSEF/NCorePClu,MCorePClu,NCluPNode,NManageCore,&
      NGrpPProc ,NThPGrp ,NProcPNode,ManageCoreId,&
      OffClu, SkipClu, OffCore, SkipCore

    INTEGER ::i,j,k,minv,minp,minpt(10)

    integer ios
#ifndef NO_MPI
    type(mpipacket) pk
#endif
    runstate=>runstate_t
    gc2=0.877**2*g;    zpi=2*pi;deg2m=Rs*pi/180.
    if(mpi_id==0)then
      call OpenCtl;read(ctlfid,ParaCtl,IOSTAT=ios)
      if(ios/=0)then
        call OpenCtl(1)
        write(ctlfid,*)'# =====ParaCtl:是否输出默认控制参数========'
        write(ctlfid,*)'# OutputAllPara:not Equal 0,OutPut All Control para'
        write(ctlfid,*)'# =================================='
        write(ctlfid,ParaCtl)
      end if
      call OpenCtl;      read(ctlfid,runDef,IOSTAT=ios)
      if(OutputAllPara/=0.and.ios/=0)then
        call OpenCtl(1)
        write(ctlfid,*)'# Rundef:运行参数'
        write(ctlfid,*)'# startType:startup:Cold Start 0;'
        write(ctlfid,*)'#     hybrid:restart Use Data Only;'
        write(ctlfid,*)'#     branch:restart Change Case;'
        write(ctlfid,*)'#     continue: restart;'
        write(ctlfid,*)'# caseid   :caseName'
        write(ctlfid,*)'# pathinit :Init File Directory'
        write(ctlfid,*)'# pathrest :restart File Directory'
        write(ctlfid,*)'# pathhist :History File Directory'
        write(ctlfid,*)'# pathloc  :History File Local Directory'
        write(ctlfid,*)'# pathwind :Wind Data Directory'
        write(ctlfid,*)'# pathcur  :Cur  Data Directory'
        write(ctlfid,*)'# topofn   :topography File'
        write(ctlfid,*)'# topofType:topography File Type 1:text,3:NetCdf'
        write(ctlfid,*)'# DEPTHMOD : 0'
        write(ctlfid,*)'# DBGLVL   :Debug level'
        write(ctlfid,*)'# RUNDAYS  :RUNDAYS'
        write(ctlfid,*)'# StartDate:StartDate'
        write(ctlfid,*)'# StartTime:StartTime'
        write(ctlfid,*)'# PreCalTime:Define PreCal Time in hour'
        write(ctlfid,runDef)
      end if
      call OpenCtl;      read(ctlfid,WindForm,IOSTAT=ios)
      if(OutputAllPara/=0.and.ios/=0)then
        call OpenCtl(1)
        write(ctlfid,WindForm)
      end if
      call OpenCtl;      read(ctlfid,GPARA,IOSTAT=ios)
      if(OutputAllPara/=0.and.ios/=0)then
        call OpenCtl(1)
        write(ctlfid,*)'#!!!!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!!!!!!!!!!'
        write(ctlfid,*)'# ===========GPARA: 全局参数========='
        write(ctlfid,*)'# (gxlon0-gxlon1,gylat0-gylat1):计算区域'
        write(ctlfid,*)'# deltts:时间步长(秒钟),grdszx,grdszy:x、y网格大小'
        write(ctlfid,*)'# 1deg:1h/3,0.5deg:1h/5,1/4deg:1h/10(1h/9不行)'
        write(ctlfid,*)'# Gcircle:全球模式,环绕'
        write(ctlfid,*)'# logscurr:THE SWITCH FOR WAVE-CURRENT INTERACTION'
        write(ctlfid,*)'# logsdiss:THE SWITCH FOR Sds,1-Yuan''s,0-Komen''s'
        write(ctlfid,*)'# =================================='
        write(ctlfid,GPARA)
      end if
      call OpenCtl;      read(ctlfid,BVDepDef,IOSTAT=ios)
      if(OutputAllPara/=0.and.ios/=0)then
        call OpenCtl(1)
        write(ctlfid,*)'# BVDepDef:bv layer define '
        write(ctlfid,*)'# ndep :BV N layers'
        write(ctlfid,*)'# vdep :bv layers '
        write(ctlfid,BVDepDef)
      end if
      call OpenCtl;      read(ctlfid,OUTdef,IOSTAT=ios)
      if(OutputAllPara/=0.and.ios/=0)then
        call OpenCtl(1)
        write(ctlfid,*)'# OUTdef:输出参数'
        WRITE(ctlfid,*)'# AVHIST_OPTIONS,AVHIST_NS:'
        WRITE(ctlfid,*)'# REST_OPTION,REST_N:'
        WRITE(ctlfid,*)'#   "AUTO":use CPL Control;"NONE","NEVER":Not Output;'
        WRITE(ctlfid,*)'#   "NSTEPS":N Steps;"NHOURS: N Hours;"NDAYS":N Days;"NMONTHS":N months;"NYEARS":N Years'
        WRITE(ctlfid,*)'#   "daily": 1day;"monthly":1 month,"yearly":1 year'
        WRITE(ctlfid,*)'#   if(N<0):Mean OutPut not Average,transient field'
        WRITE(ctlfid,*)'# AVHIST_NTPFS :每文件输出次数,>0: n Times;0: 1 day;-1: 1 month;-2: 1 Season;-3:half year:<=-4: -n-4 year '
        WRITE(ctlfid,*)'# AVHIST_OTAGS: OutName Tag,default:"ahi"'
        WRITE(ctlfid,*)'# AVHIST_TAGS,AVHIST_OPTION,AVHIST_NS,AVHIST_NTPFS:is array surport Multi(max 16) Output'
        write(ctlfid,OUTdef)
      end if
      !ZZZ
      call OpenCtl;      read(ctlfid,implschPara,IOSTAT=ios)
      if(OutputAllPara/=0.and.ios/=0)then
        call OpenCtl(1)
        write(ctlfid,*)'# =====implschPara:源函数系数========'
        write(ctlfid,*)'# beta0,beta1:WIND INPUT COEFFICIENT'
        write(ctlfid,*)'# brkd1,brkd2:WIND INPUT COEFFICIENT'
        write(ctlfid,*)'# acu:Current COEFFICIENT'
        write(ctlfid,*)'# =================================='
        write(ctlfid,implschPara)
      end if
#if 0
      call OpenCtl;      read(ctlfid,PKPARA,IOSTAT=ios)
      if(OutputAllPara/=0.and.ios/=0)then
        call OpenCtl(1)
        write(ctlfid,*)'# PKPARA:波束划分参数'
        write(ctlfid,*)'# KL波束划分数 ,kld,总波束划分数,jnthet 方向划分数'
        write(ctlfid,*)'# wk(k)=wkmin*(pwk**(k-1))'
        write(ctlfid,PKPARA)
      end if
#endif
      call OpenCtl;      read(ctlfid,OUTPUTPnt,IOSTAT=ios)
      if(OutputAllPara/=0.and.ios/=0)then
        call OpenCtl(1)
        WRITE(ctlfid,*)'# OUTPUTPnt:监测/工作点坐标'
        WRITE(ctlfid,OUTPUTPnt)
      end if
      call OpenCtl;      read(ctlfid,testOpt,IOSTAT=ios)
      if(OutputAllPara/=0.and.ios/=0)then
        call OpenCtl(1)
        write(ctlfid,testOpt)
      end if
      call OpenCtl;      read(ctlfid,CPLDEF,IOSTAT=ios)
      if(OutputAllPara/=0.and.ios/=0)then
        call OpenCtl(1)
        write(ctlfid,*)'# CPLDEF:CCSM 耦合参数'
        write(ctlfid,*)'# NCPL :每天耦合次数'
        write(ctlfid,CPLDEF)
      end if
      call OpenCtl;      read(ctlfid,THEADSEF,IOSTAT=ios)
      if(OutputAllPara/=0.and.ios/=0)then
        call OpenCtl(1)
        write(ctlfid,*)'# 系统参数'
        write(ctlfid,*)'# NCorePClu  :每簇(NUMA)核数(38)'
        write(ctlfid,*)'# MCorePClu  ;每簇(NUMA)可用核数(38),假设不可用核在最后'
        write(ctlfid,*)'# NCluPNode  :每节点簇数(16)'
        write(ctlfid,*)'# NManageCore;管理核核数，编号>NCorePClu*NCluPNode'
        write(ctlfid,*)'# 软件设置参数'
        write(ctlfid,*)'# NThPGrp :每组启动线程数(38),2<=NProcPNode<=MCorePClu'
        write(ctlfid,*)'# NGrpPProc:每进程组数(16),<=NCluPNode ,支持 4,8,16'
        write(ctlfid,*)'# NProcPNode:每节点进程数(1),<=NCluPNode/NGrpPProc'
        write(ctlfid,*)'# ManageCoreId:管理核(-1)'
        write(ctlfid,*)'# 测试用参数'
        write(ctlfid,*)'# OffClu   :第一组使用的簇(0)'
        write(ctlfid,*)'# SkipClu  :几个簇启动1个组间(1)'
        write(ctlfid,*)'# OffCore  :每个组使用的第一个核(0)'
        write(ctlfid,*)'# SkipCore :几个核启动一个线程(1)'
        write(ctlfid,THEADSEF)
      end if
      close(ctlfid)
      call UPPERCASE(startType);
      call UPPERCASE( REST_OPTION);
      call ADJ_OPTION(REST_OPTION, REST_N,REST_TYPE);
      do i=1,IMAXAVHIST
        call UPPERCASE(AVHIST_OPTIONS(i));
      enddo
      !write(6,*)AVHIST_OPTIONS,AVHIST_NS
    endif
#ifndef NO_MPI
    call initmpipacket(pk,10240)
    if(mpi_id==0)then
      call PackVars(0)
      call SetMpiPacketDSize(pk)
      call wav_mpi_bcast(pk);
    else
      call wav_mpi_bcast(pk);
      call PackVars(1)
    endif
    call InitMpiPacket(pk,-1)
#endif

    deltts=(24*3600)/int((24*3600)/deltts)
    NSCPL=(86400)/NCPL/deltts
    Delttscpl=86400/NCPL
    delltday=deltts/86400.d0
    !if(NCPL*NSCPL*deltts/=86400)then
    !    print*,'Couple Time Not divisible Model time step'
    !    call wav_abort('Couple Time Not divisible Model time step')
    !endif
    alog10pwk=log10(pwk)
    dlog10pwk=1/alog10pwk
    dwkmin=1/wkmin
    klp=kl+KLPAT;! for c code Pache
    deltth=zpi/float(jnthet)
    ddeltth=1/deltth
    beta10=beta0*0.25*1.25*0.001
    gixl=(gxlon1-gxlon0)/grdszx+0.5
    if(Gcircle==0)gixl=gixl+1
    giyl=(gylat1-gylat0)/grdszy+1.5
    call SetBaseYear(1,1950)
    !DBGO(0,*)StartDate
    call wav_cal_date2eday(StartDate, dayinit)
    timeInit=ctime2day(StartTime)

    runstate%iPreCalT=3600*PreCalTime/deltts
    runstate%mTimeStep=RunDays*24*3600/deltts+runstate%iPreCalT
    runstate%nTimeStep=-1
    runstate%edaycur=dayinit
    runstate%timecur=timeinit
    runstate%ndays=0
    if(Gcircle==0)gixl=gixl+1
#if LOGSCURR==0
    if(logscurr/=0)then
      DBGO(0,*)  'LOGSCURR Not Defined,but logscurr used'
      call wav_abort('LOGSCURR Not Defined,but logscurr used')
    end if
#endif
 !! write(6,ParaCtl,IOSTAT=ios)
 !! write(6,Rundef,IOSTAT=ios)
 !! write(6,GPARA,IOSTAT=ios)
 !! write(6,BVDepDef,IOSTAT=ios)
 !! write(6,OUTPUTPnt,IOSTAT=ios)
 !! write(6,implschPara,IOSTAT=ios)
 !! write(6,PKPARA,IOSTAT=ios)
 !! write(6,CPLDEF,IOSTAT=ios)
    !InitInOutOpt
    Navhists=IMAXAVHIST
    Navhists=0
    MinvS=10000000
    do i=1,IMAXAVHIST
      call ADJ_OPTION(AVHIST_OPTIONS(i),AVHIST_NS(i),AVHIST_TYPES(i))
      if(AVHIST_TYPES(i)>=0)then
        Navhists=Navhists+1
        if(Navhists/=i)then
         AVHIST_OTAGS  (Navhists)=AVHIST_OTAGS  (i)
         AVHIST_OPTIONS(Navhists)=AVHIST_OPTIONS(i)
         AVHIST_NS     (Navhists)=AVHIST_NS     (i)
         AVHIST_NTPFS  (Navhists)=AVHIST_NTPFS  (i)
         LogbvS        (Navhists)=LogbvS        (i)
         LogWindS      (Navhists)=LogWindS      (i)
         LoghsS        (Navhists)=LoghsS        (i)
         LogtpS        (Navhists)=LogtpS        (i)
         LogtzS        (Navhists)=LogtzS        (i)
         LogthS        (Navhists)=LogthS        (i)
         LogauS        (Navhists)=LogauS        (i)
        endif
        MinvS(Navhists)=AVHIST_TYPES(Navhists)*10000+abs(AVHIST_NS(Navhists))
      endif
    enddo
    Logau=0;Logbv=0;LogWind=0;Loghs=0;Logtp=0;Logtz=0;Logth=0
    if(Navhists>0)then
      allocate(avhists(Navhists)) ;
      do i=1,Navhists
        minpt(1:1)=minloc(MinvS)
        minp=minpt(1)
        MinvS(minp)=10000000
        nullify(avhists(i)%ea)
        nullify(avhists(i)%aet )
        nullify(avhists(i)%tpf )
        nullify(avhists(i)%h1_3)
        nullify(avhists(i)%ape )
        nullify(avhists(i)%bv  )
        avhists(i)%TAG   =AVHIST_OTAGS  (minp)
        avhists(i)%OPTION=AVHIST_OPTIONS(minp)
        avhists(i)%NN    =abs(AVHIST_NS (minp))
        avhists(i)%itype =AVHIST_TYPES  (minp)
        avhists(i)%NTPF  =AVHIST_NTPFS  (minp)
        if(AVHIST_NS(minp)<0)then
          avhists(i)%mean=0
        else
          avhists(i)%mean=0
        endif
        avhists(i)%Logbv  =LogbvS  (minp)
        avhists(i)%LogWind=LogWindS(minp)
        avhists(i)%Loghs  =LoghsS  (minp)
        avhists(i)%Logtp  =LogtpS  (minp)
        avhists(i)%Logtz  =LogtzS  (minp)
        avhists(i)%Logth  =LogthS  (minp)
        avhists(i)%Logau  =LogauS  (minp)
        if(avhists(i)%Logbv  /=0)Logbv  =Logbv  +1
        if(avhists(i)%LogWind/=0)LogWind=LogWind+1
        if(avhists(i)%Loghs  /=0)Loghs  =Loghs  +1
        if(avhists(i)%Logtp  /=0)Logtp  =Logtp  +1
        if(avhists(i)%Logtz  /=0)Logtz  =Logtz  +1
        if(avhists(i)%Logth  /=0)Logth  =Logth  +1
        if(avhists(i)%Logau  /=0)Logau  =Logau  +1 
        ! write(6,*)i,avhists(i)%OPTION,avhists(i)%NN,minp,minv
      enddo
      !DBGO0(0,*)avhists
    endif
#ifndef NO_MPI
    CONTAINS
!#define S(v) #v
!#define PKV(v) print*,__LINE__,shape(v);call wav_mpi_pack(pk,v,iunpack );print*,__LINE__,shape(v)
#define PKV(v) call wav_mpi_pack(pk,v,iunpack );
    SUBROUTINE PackVars(iunpack)
      integer iunpack
      PKV(gxlon0    );PKV(gylat0    )
      PKV(gxlon1    );PKV(gylat1    )
      PKV(grdszx    );PKV(grdszy    )
      PKV(Gcircle   );PKV(windfield );
      PKV(DEPTHMOD  )
      PKV(deltts    );
      PKV(LogWind   )
      PKV(Logau     );PKV(Logbv     )
      PKV(Loghs     );PKV(Logtp     )
      PKV(Logtz     );PKV(Logth     )

      PKV(ndep      );PKV(vdep      )

      PKV(plo       );PKV(pla       )

      PKV(logscurr  );PKV(logsdiss  )
      PKV(startType );PKV(caseid    )
      PKV(pathinit  );PKV(pathrest  )
      PKV(pathhist  );PKV(pathLoc   );
      PKV(pathcur   );
      PKV(topofn    );PKV(DbgLvl    );
      PKV(RunDays   );PKV(StartDate );
      PKV(StartTime );PKV(PreCalTime);
      PKV(outputNEP );PKV(OutType   );
      PKV(PATHWIND  );PKV(WTIME_PATCH);
      PKV(uwind_name);PKV(vwind_name);

      PKV(ufnform   );PKV(vfnform   );
      PKV(noOutPut  );PKV(iCheckRestart);
      PKV(constwindx);PKV(constwindy);

      PKV(pwk       );PKV(wkmin      )
      PKV(KL        );PKV(kld        )
      PKV(jnthet    );

      PKV(beta0     );PKV(beta1      )
      PKV(acu       );PKV(brkd1      )
      PKV(brkd2     )

      PKV(REST_OPTION );PKV(REST_N        );
      PKV(AVHIST_OTAGS);PKV(AVHIST_OPTIONS);
      PKV(AVHIST_NS   );PKV(AVHIST_NTPFS  );
      PKV(LogbvS    );  PKV(LogWindS);
      PKV(LoghsS    );  PKV(LogtpS  );
      PKV(LogtzS    );  PKV(LogthS  );
      PKV(LogauS    );
      PKV(NCPL      );
      PKV(NCorePClu );PKV(MCorePClu   );
      PKV(NCluPNode );PKV(NManageCore );
      PKV(NGrpPProc );PKV(NThPGrp     );
      PKV(NProcPNode);PKV(ManageCoreId);
      PKV(OffClu    );PKV(SkipClu     );
      PKV(OffCore   );PKV(SkipCore    );
    end SUBROUTINE PackVars
#endif
  end SUBROUTINE InitPara !}

  SUBROUTINE ADJ_OPTION( opt,N,itype)
    character*(*) opt
    integer N,itype
    itype=0
    if(opt=='AUTO')then
      N=0;itype=0
    else if(opt=='NONE')then
      opt='NONE';N=1000;itype=-1
    else if(opt=='NEVER')then
      opt='NONE';N=1000;itype=-1
    else if(opt=='DAILY')then
      opt='NDAYS'; if(N<0)then; N=-1;else; N=1;endif
    else if(opt=='MONTHLY')then
      opt='NMONTHS';if(N<0)then;N=-1;else; N=1;endif
    else if(opt=='YEARLY')then
      opt='NYEARS'; if(N<0)then;N=-1;else; N=1;endif
    endif
    if(opt=='NSTEPS')then
      itype=1
    else if(opt=='NHOURS')then
      itype=2
    else if(opt=='NDAYS')then
      itype=3
    else if(opt=='NMONTHS')then
      itype=4
    else if(opt=='NYEARS')then
      itype=5
    else
      opt='NONE';N=0;itype=-1
    endif
  end SUBROUTINE ADJ_OPTION

!===============================================================================================
  SUBROUTINE setwave !{
    integer j,k
    real(8) wh
    integer ierr
    ALLOCATE(wk(kld+1),wkh(kld+1),dwk(kld),wkdk(kld),wkibdk(kld),STAT=ierr)
    ALLOCATE(thet(jnthet+1),sinths(jnthet+1),cosths(jnthet+1),cosths2(jnthet+1),sinths2(jnthet+1),STAT=ierr)
    ALLOCATE(bcosths2(jnthet+1),bsinths2(jnthet+1),sincosths(jnthet+1),STAT=ierr)
    DO  j=1,jnthet+1 !{
        thet(j)=(j-1)*deltth
        sinths(j)=sin(thet(j))
        cosths(j)=cos(thet(j))

        bcosths2(j)=1+cosths(j)*cosths(j)
        bsinths2(j)=1+sinths(j)*sinths(j)
        cosths2(j)=cosths(j)*cosths(j)
        sinths2(j)=sinths(j)*sinths(j)
        sincosths(j)=cosths(j)*sinths(j)
    end do !}
    wkh(1)=1.;    wh=sqrt((1./pwk)**7)
    !(wk(k+1)-wk(k-1))=wkmin*(pwk**(k)-pwk**(k-2)=(pwk-1/pwk)*wk(k)
    !(wk(k+1)-wk(k))=wkmin*(pwk**(k)-pwk**(k-1)=(pwk-1)*wk(k)
    DO  k=1,kld !{
        wk(k)=wkmin*(pwk**(k-1)) !discretion of wave number 
        wkh(k)=wh**(k-1)
         if(k==1)then
           dwk(k)=(pwk-1)*(wk(k)**2) *deltth/2
         else if(k==kld)then
           dwk(k)=(1-1/pwk)*(wk(k)*wk(k)) *deltth/2
         else
           dwk(k)=(pwk-1/pwk)*(wk(k)**2) *deltth/2
         endif
         wkdk(k)=wk(k)*dwk(k)
         wkibdk(k)=1/sqrt(wk(k))*dwk(k)
    end do !}
    wk(kld+1)=wkmin*(pwk**(kld)) !discretion of wave number
    wkh(kld+1)=wh**(kld)
    return
  end SUBROUTINE setwave !}
  integer function IndWk(w)
    real(8) w
    IndWk=floor(log10(w)*dlog10pwk)
  end function IndWk
  integer function IndWkd(w)
    real(8) w
    IndWkd=log10(w*dwkmin)*dlog10pwk
  end function IndWkd

  subroutine OpenCtl(endf) !{
    integer,optional::endf
    integer::endf_=0
    integer ios
    character*(8) pos
    if(ctlOpened/=0)close(ctlfid)
    if(present(endf))endf_=endf
    if(endf_/=0)then
      pos="append";
    else
      pos="rewind"
    endif
    open(ctlfid,file='nwamctl.ini',status='old',position=pos,DELIM='QUOTE',IOSTAT=ios)
    if(ios/=0)then
      open(ctlfid,file='nwamctl.ini',status='new',DELIM='QUOTE')
    end if
    ctlOpened=1
  end subroutine OpenCtl

  SUBROUTINE UPPERCASE(s)
    character*(*) s
    integer len,i,c
    len=len_trim(s)
    do i=1,len
      c=ichar(s(i:i))
      if(c>=97.and.c<=122)s(i:i)=char(c-32)
    end do
  end SUBROUTINE UPPERCASE
  SUBROUTINE LOWERCASE(s)
    character*(*) s
    integer len,i,c
    len=len_trim(s)
    do i=1,len
      c=ichar(s(i:i))
      if(c>=65.and.c<=90)s(i:i)=char(c+32)
    end do
  end SUBROUTINE LOWERCASE
  subroutine NextTime(time,opt,N,mean)
    real(8) time
    character*(*) opt
    integer N,mean
    integer y,m,d,ed
    if(time<1e-4)then
      time=runstate%edaycur+runstate%timecur
    endif
    if(opt=='AUTO')then
      time=2.d300
    else if(opt=='NONE')then
      time=1.d300
    else if(opt=='NSTEPS')then
      if(time<1e-4)then
        time=int((runstate%edaycur+runstate%timecur)/(N*deltts/86400.)+0.5)*(N*deltts/86400.)
      else
        time=time+N*deltts/86400.
      endif
    else if(opt=='NHOURS')then
      time=int(time*24+N)/24.
    else if(opt=='NDAYS')then
      if(mean/=0)then
        time=int(time)+N
      else
        time=int(time-0.5)+N+0.5
      endif
    else if(opt=='NMONTHS')then
      call wav_cal_eday2ymd(int(time),y,m,d)
      m=m+N;
      if(m>12)then
        y=y+(m-1)/12;
        m=mod(m-1,12)+1
      endif
      call wav_cal_ymd2eday(y,m,d,ed)
      time=ed+0.5
    else if(opt=='NYEARS')then
      call wav_cal_eday2ymd(int(time),y,m,d)
      y=y+N;
      call wav_cal_ymd2eday(y,1,1,ed)
      time=ed
    endif
  end subroutine NextTime
  SUBROUTINE ResetNspg(nspt)
    INTEGER nspt(:,:)
    integer i1,j1
      DO  j1=1,giyl !{
      DO  i1=1,gixl !{        
        IF(depg(i1,j1)>1.e-5)then
          nspt(i1,j1)=1
        else
          nspt(i1,j1)=0
        endif
      end do !}
      end do !}
  
      if(Gcircle==0)then
        DO  j1=1,giyl !{
          if(nspt(   1,j1)/=0)nspt(   1,j1)=2
          if(nspt(gixl,j1)/=0)nspt(gixl,j1)=2
        end do !}
      endif
      DO  i1=1,gixl !{
        if(nspt(i1,   1)/=0)nspt(i1,   1)=0 !ZZZ
        if(nspt(i1,giyl)/=0)nspt(i1,giyl)=0 !ZZZ
      end do !}     
      if(gylat1>=78)then
        DO  i1=11,gixl !{
          nspt(i1,giyl)=0
        end do !}
      endif  
  end SUBROUTINE ResetNspg 
subroutine InterGrd(sv,dv,lon,nlon,lat,nlat,iv_)
integer nlon,nlat
real(4),intent(in)::sv(:,:),lat(nlat),lon(nlon)
real(4),intent(out)::dv(:,:)
real(4),target,optional::iv_(:,:)
real(4),pointer::iv(:,:)
integer i,j,ij,it,it1,jt,jt1
real x,y,p
  it1=1
  if(present(iv_))then
    iv=>iv_
  else
    allocate(iv(gixl,nlat))   
  endif
  do i=1,gixl
    x=xcord(i,1)
    do ij=it1,nlon
      if(x<lon(ij))exit
    enddo
    it1=ij;it=it1-1
    if(it<=0)then
      it=it+nlon
      p=(x+360-lon(it))/(lon(it1)+360-lon(it))
    else
      p=(x-lon(it))/(lon(it1)-lon(it))
    endif
    do j=1,nlat
        iv(i,j)=sv(it,j)*(1-p)+sv(it1,j)*p    
        ! write(6,'(4i4,4f8.2)')i,j,it,it1,sv(it,j),sv(it1,j),p,iv(i,j)
    enddo
  enddo
  if(lat(1)<lat(nlat))then
    jt1=1
    do j=1,giyl
      y=ycord(1,j)
      do ij=jt1,nlat
        if(y<lat(ij))exit
      enddo
      jt1=ij;jt=jt1-1
      if(jt<=0)then
        jt=jt1;p=0
      else if(jt1>nlat)then
        jt=nlat;jt1=nlat;p=0
      else
        p=(y-lat(jt))/(lat(jt1)-lat(jt))
      endif
      do i=1,gixl
          dv(i,j)=iv(i,jt)*(1-p)+iv(i,jt1)*p      
      enddo
    enddo
  else
    jt1=nlat
    do j=1,giyl
      y=ycord(1,j)
      do ij=jt1,1,-1
        if(y<lat(ij))exit
      enddo
      jt1=ij;jt=jt1+1
      if(jt>nlat)then
        jt=jt1;p=0
      else if(jt1<1)then
        jt=1;jt1=1;p=0
      else
        p=(y-lat(jt))/(lat(jt1)-lat(jt))
      endif
      if(jt==jt1)then
        dv(:,j)=iv(:,jt)
      else
        do i=1,gixl
            dv(i,j)=iv(i,jt)*(1-p)+iv(i,jt1)*p      
        enddo
      endif
    enddo
  endif
  if(.not. present(iv_))deallocate(iv)
end subroutine InterGrd
  subroutine checkee(ee,n,fn,ln)
    integer n,ln
    character *(*) fn
    REALD ee(CKL*CJNTHET*n)
    integer k,j,kj,iac,ie
    call ccheckee(ee,n*CKL*CJNTHET,ie);
    if(ie/=0)then
      iac=ie/(CKL*CJNTHET)
      kj=mod(ie,CKL*CJNTHET);
      print'(a,a,2i4,f0.5,3i4)',"EE ERR",fn,ln,mpi_id,ee(ie),iac,kj
      call wav_abort("EE Check Error");
      stop
    endif
  end subroutine 
#ifndef NO_MPI
  subroutine gathercheck(line)
    integer iv,i,line,ierr
    integer,allocatable:: ivs(:)
    allocate(ivs(0:mpi_npe-1))
    iv=mpi_id
    call wav_mpi_gather(iv,ivs)
    if(mpi_id==0)then
      do i=0,mpi_npe-1
        if(ivs(i)/=i)print*,"gather error",i,ivs(i);
      enddo
    endif
    deallocate(ivs)
  end subroutine gathercheck
#endif
end Module varcommon_mod

