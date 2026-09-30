#include "wavedeff.h"
!#define ALLLOG
!# define DBGINF   !call wav_mpi_barrier("AA");write(*,'(a,i6,i4.4,";",18i8)')__FILE__,__LINE__,mpi_id

#define IMDPSP  2
#define IMDPSO  4
#define BVIMDPS 1
Module waveinit_mod
#ifndef NO_MPI
  use wav_mpi_mod
#endif
  use varcommon_mod
  use implsch_mod
  use propagat_mod
  use boundary_mod
  use windin_mod
  use output_mod
  use partition_mod
  use restart_mod
  use platform_init_mod
IMPLICIT NONE
private
  integer::hist,rest
  public::InitWaveMdl,EndWaveMdl,DecIPrecalT,SetModelTime
  integer ::ndstep,oldday=-1
  character(152):: headf=" pid    date   time Days                  dtime    adtime  alltime exchange  propspec propgeo smooth implschs    other   accuma    Input   output rest"
  character(60):: headp=" date   time   Days       dtime      adtime     alltime" !exchange  propspec propgeo smooth implschs    other   accuma    Input   output rest"
  real(8)  dfpc
  type FC_PARA
    integer(4) kl_,jnthet_,kld_,NBVDEP_,sigmalvl_,logscurr_
    integer(4) gixl,giyl;
    integer(4) nwps,nwpc,nwpa;
    integer LXB,LYB,LXN,LYN
    integer(4) mpi_id,mpi_comm_wav,mpi_npe;
    real   (4) :: constwindx,constwindy
    integer(4) ::mTimeStep
  end type FC_PARA
contains
#ifdef ALLLOG
#define RINFO(TAG) write(*,'(i8," ",a,f8.3," ",2i8)')iwalltime(),TAG,Difftimer(2),mpi_id,mpi_npe;  call flush(6);
#else
#define RINFO(TAG) if(mpi_id==0)then;write(*,'(i8," ",a,f8.3," ",i8)')iwalltime(),TAG,Difftimer(2),mpi_npe;  call flush(6);endif
#endif
#define DBGA call wav_mpi_barrier;write(*,'("DBGA",i4,":",i4.3)')__LINE__,mpi_id;  call flush(6);
#define DBGB(tag) write(*,'("DBGA",i4,a,":",4i4.3)')__LINE__,tag,mpi_id,mpi_npe;  call flush(6);
  SUBROUTINE InitWaveMdl
    integer key,ir,out,ierr
    real*8 tt
    !ir=ieee_flags('clear', 'exception', 'all',out )
    call Resettimer(-1)
    call starttimer(1)
    call starttimer(2)
    if(IDataIO(0)<0)then   ! IGetWind(0) Init MPI msg etc,must before partition
      !DBGO(0,*)"Init MPI msg  Error"
      stop
    endif
    call InitMpi(0)      ! Mpi Init
    call zf_setpid(mpi_id) ! Set Mpi_id for zfile.c
    RINFO('InitMpi End')
#ifndef NO_MPI
    !call wav_mpi_gather(nodeid,NodeIds)
    !RINFO('Gather nodeid End')
    call gathercheck(__LINE__)
    RINFO('gathercheck A End')
#endif
    !call testtimes
    !call zftest;         ! call zftest;         ! stop
    Call InitPara      !varcommon_mod Init
    RINFO('InitPara End')
    call InitAfterMpi
    DBGLVL=20
    ! call wav_chStdOut
    RINFO('InitAfterMpi End')
    if(mpi_id>0) then
      if(DBGLvl<5)DBGLvl=0
    endif
    call setwave
    RINFO('setwave End')
    !call nlweight;    stop
    !RINFO('ALLOCATE depg,nspg begin')
    ALLOCATE(depg(gixl,giyl),nspg(gixl,giyl));
    !RINFO('ALLOCATE depg,nspg End')
    call ReadTopog
    RINFO('ReadTopog End')
#ifndef NO_MPI
    call gathercheck(__LINE__)
    RINFO('gathercheck B End')
#endif
    !partition !!!!!
    call partition
    RINFO('partition End')
    !call ReadTopog;    RINFO('ReadTopog End')
    !call InitMonitor
    !RINFO('InitMonitor End')
    if(mpi_id==0)then
      DBGO(0,'(a,i7  ,a,i7  )') 'KL     =',kl     ,' DBGLVL  =',DBGLVL
      DBGO(0,'(a,f7.2,a,f7.2)') 'GXLON0 =',GXLON0 ,' GYLAT0  =',GYLAT0
      DBGO(0,'(a,f7.2,a,f7.2)') 'GXLON1 =',GXLON1 ,' GYLAT1  =',GYLAT1
      DBGO(0,'(a,f7.4,a,f7.4)') 'GRDSZX =',GRDSZX ,' GRDSZY  =',GRDSZY
      DBGO(0,'(a,i7  ,a,i7  )') 'GCIRCLE=',GCIRCLE,' DEPTHMOD=',DEPTHMOD
      DBGO(0,'(a,f7.2       )') 'DELTTS =',DELTTS
      DBGO(0,'(a,2i5,f12.2,4i10 )')'ix/yl,nwpx',gixl,giyl,(gnwpc+0.)/mpi_npe,nwpa,nwpc,nwps,gnwpc
      DBGO(4,'(10i9)')nwpcs
      call flush(6);
    endif
    ! sum case Needed allocate Mem use C Program for mem align
    ! if Not ,InitAfterPartition  is null
    !print*,__FILE__,__LINE__,mpi_id
    !call wav_mpi_barrier();
    RINFO('InitAfterPartition start')
    call InitAfterPartition
    RINFO('InitAfterPartition')
    call set_neighbor
    RINFO('set_neighbor')
    call InitOutPut(0)   !Init OutPut
    RINFO('InitOutPut End')
    key=IDataIO(1)            !must Init before InitImplsch

    call InitImplscha
#ifdef USEXNL
    call xnl_inite( wk, thet, kl, jnthet, -5.,  2, 1, ierr )
#endif
    call InitPropagata
    call InitSetspeca
    call InitCheckOutPut(0)
#ifndef MTHREAD
    call InitSetspec(1,nwpc)
    RINFO('InitSetspec End')
    call InitImplschs(1,nwpc)
    RINFO('InitImplsch End')
    call InitPropagat(1,nwpc)
    RINFO('InitPropagat End')
    call InitOutPutCal(1,nwpc,dep);
    RINFO('InitOutPutCal End')
#endif
    !call wav_mpi_barrier();
    ! For every case ,sum Distinct Init here
    call Init_Distinct
    RINFO('Init_DistinctEnd')
    call endtimer(1)
    ! Init run
    !if(iCheckRestart>0)then
    !  call InitCheckOutPut(1)
    !  runstate%edaycur=dayinit
    !  runstate%timecur=timeinit
    !  runstate%ndays=0
    !  call setmodeltime
    !  call TestRestart(iCheckRestart,eec)
    !  RINFO('TestRestart End')
    !  stop
    !endif
    RINFO( "wave model Init End")
    call NextTime(resttime,rest_option,rest_N,0)
    ndstep=3600/DELTTS
    if(ndstep<=1)ndstep=1
    DBGO0(4,'(a)')trim(startType)
    RINFO( "wave model AA")
    call Init_Startup
    RINFO( "wave model AB")
  end SUBROUTINE InitWaveMdl !}

  SUBROUTINE DecIPrecalT
    integer nn
    if( mod(runstate%iPreCalT,ndstep)==0)then
      !DBGO0(0,'("p",$)')
    endif
    nn=3600*4/DELTTS;if(nn<1)nn=1
    runstate%iPreCalT=runstate%iPreCalT-1
    if( mod(runstate%iPreCalT,nn)==0)then
       DBGO0(1,'(a,f10.2,a,i8,a,i8,a)')'runstate%iPreCalT use',gettimer(1),"s",runstate%iPreCalT,' leave ',int(runstate%iPreCalT*DELTTS/3600),'hours'
       !call Monitor(eec,-1,0)

    endif
  end SUBROUTINE DecIPrecalT

  SUBROUTINE EndWaveMdl
    integer ir
    ir=IDataIO(-1)
    call InitOutPut(-1)  ! Close  OutPut
    DBGO0(0,*) "end wave model A",runstate%nTimeStep
    call End_Distinct
    call InitMpi(-1)      ! Mpi Deinit
    call flush6
  end SUBROUTINE EndWaveMdl

  subroutine SetModelTime(RFB)
    integer RFB
    real*8 dt,adt,tt,nfpca,ta
    integer  eday
    if(runstate%timecur>1-delltday*0.01)then
      runstate%timecur=runstate%timecur-1
      runstate%edaycur=runstate%edaycur+1
    endif
    call day2cdatetime(runstate%edaycur+runstate%timecur,runstate%cdatecur,runstate%ctimecur)
    if(noOutPut==0 .and.runstate%iPreCalT<=0)then  !
      if(runstate%NextOutStep<runstate%nTimeStep)then
        call CheckNeedOutPut(runstate%hist_eot)
      endif
      call CheckNeedRestart(runstate%rest_eot)
    endif
#define RT(l) gettimer(l)
    if(oldday/=runstate%edaycur)then
      runstate%ndays=runstate%ndays+1
      oldday=runstate%edaycur
      dt =Difftimer(1); call endstarttimer(1,1);
      if(runstate%ndays>1)then
        adt=RT(1)/(runstate%ndays-1)
      else
        adt=dt;
      endif
      ta=RT(2)
      if(runstate%ndays>RunDays)runstate%stop_now=1
      if(mpi_id==0)then
        if(runstate%ndays==0.or.mod(runstate%ndays,50)==1)then
          !write(1,*)headp
          write(*,*)headp
        endif
        !write(1,'(i8.8,i7.6,i5,3f12.5,":","T:",i8,7f10.5)')runstate%cdatecur,runstate%ctimecur,runstate%ndays-1 ,dt,adt, ta,iwalltime()  
        write(*,'(i8.8,i7.6,i5,3f12.5,":","T:",2i8)')runstate%cdatecur,runstate%ctimecur,runstate%ndays-1 ,dt,adt, ta,runstate% nTimeStep,iwalltime()  
        call flush6
#ifdef ALLLOG
      else
        if(runstate%ndays==0.or.mod(runstate%ndays,50)==1)then
          write(*,*)headp
        endif
        write(*,'(i8.8,i7.6,i5,3f12.5,":","T:",i8,7f10.5)')runstate%cdatecur,runstate%ctimecur,runstate%ndays-1 ,dt,adt, ta,iwalltime()  
        call flush6
#endif
      endif
    endif
    if(1==1)then
      if( mpi_id==0)then
        dt =Difftimer(1);
        write(*,'(i8.8,i7.6,i5,2f12.5,":","S:",3i8)')runstate%cdatecur,runstate%ctimecur,runstate%ndays-1 ,dt,ta,runstate%nTimeStep, iwalltime(),RFB  
#ifdef ALLLOG
      else
        dt =Difftimer(1);
        write(*,'(i8.8,i7.6,i5,2f12.5,":","S:",2i8)')runstate%cdatecur,runstate%ctimecur,runstate%ndays-1 ,dt,ta,runstate%nTimeStep,iwalltime()  
#endif
      endif
    endif
  END subroutine SetModelTime
  SUBROUTINE Init_Startup
    integer ir
    if(startType=='STARTUP')then
      ir=IDataIO(2)          !Get First Data for cold start
      !call wav_mpi_barrier("wind");
#ifndef MTHREAD
      call setspec(0,eec,1,nwpc,1)     ! cpp Not surport setspec
#endif
      call InitCheckOutPut(1)
      RINFO( "InitCheckOutPut End")
    else
      call read_restart(eec)
      runstate%iPreCalT=0
      ir=IDataIO(3)        !Get First Data for warm start
    end if 
  end SUBROUTINE Init_Startup
  SUBROUTINE InitAfterPartition
#ifdef C_CALCULATE 
    type(FC_PARA) gd
    integer sigmlvl
    integer res;
    !call c_checkstructs ;
    nullify(eec,eet,ip,ps_vs,pg_vs,ipos12,ipos8);
    sigmlvl=0
    !if(vdep[1]<0)  sigmlvl=1
    gd%kl_        =kl
    gd%jnthet_    =jnthet
    gd%kld_       =kld
    gd%NBVDEP_    =ndep
    gd%sigmalvl_  =sigmlvl
    gd%logscurr_  =logscurr
    gd%mpi_id     =mpi_id
    gd%mpi_npe    =mpi_npe
    gd%mpi_comm_wav=mpi_comm_wav
    gd%gixl       =gixl
    gd%giyl       =giyl
    gd%nwps       =nwps
    gd%nwpc       =nwpc
    gd%nwpa       =nwpa
    gd%lxb       =lxb
    gd%lyb       =lyb
    gd%lxn       =lxn
    gd%lyn       =lyn
    gd%constwindx =constwindx
    gd%constwindy =constwindy
    gd%mTimeStep = runstate%mTimeStep
    !call C_CheckcGrid(gd)
    call C_SetGPar(gd,iepos,ieind,nsp,dep)
#ifdef USE_C_ALLOC
    call C_Allocate_Vee(res)! Set Gvar same time
    if(res<0)then
      print*,'C_Allocate_Vee allocate mem error,',mpi_id
      call InitMpi(-1)      ! Mpi Deinit
      stop
    endif
#endif
#endif
    if(.not.associated(eet))then
      ALLOCATE(eet(kl,jnthet,0:nwpa),eec(kl,jnthet,0:nwpa));
      ALLOCATE(pvkdo(kld,IMDPSO,0:nwpc), pebdep(kld,BVIMDPS+ndep,0:nwpc));
      ALLOCATE(ip(0:nwpc))
      ALLOCATE(ps_vs(kl*jnthet,0:nwpc))
      ALLOCATE(pg_vs(kl*jnthet,0:nwpc))
      ALLOCATE(pvkdp(kl,IMDPSP,0:nwpc))
      ALLOCATE(pvws(kl,0:nwpc))
    endif
    ALLOCATE(dxym(8,0:nwpc))
    allocate(nfpcs(mpi_npe*2))
    nfpcs=0
    RINFO('InitAfterPartition A')
    eec=0;
    eet=0;
    pvkdo=0; pebdep=0;
    RINFO('InitAfterPartition B')
    pvws=0;
    pvkdp=0; ipos12=0;
    RINFO('InitAfterPartition C')
  end SUBROUTINE InitAfterPartition
end Module waveinit_mod
! called by C code for oprate Fortran variable
#ifdef USE_C_ALLOC
subroutine F_SetEVar(cec,cet)
  use varcommon_mod
IMPLICIT NONE
  REALD,target:: cec (kl,jnthet,0:nwpa)
  REALD,target:: cet (kl,jnthet,0:nwpa)
  eec=>cec;   
  eet=>cet;
    eec=0;
    eet=0;
end subroutine F_SetEVar
subroutine F_SetImpVar(cip,cpvws ,cwxy);
  use varcommon_mod
  use implsch_mod
  use boundary_mod
IMPLICIT NONE
  type(implsch_pack_type),target::cip(0:nwpc)
  real(8),target::cpvws(kl,0:nwpc)
  real(4),target::cwxy(4,0:nwpc)
  ip=>cip;
  pvws   =>cpvws
  wxy    =>cwxy;
  wxy=0.;
end subroutine F_SetImpVar
subroutine F_SetOutput(cpvkdo,cpebdep)
  use varcommon_mod
  use output_cal_mod
  IMPLICIT NONE
  real(8),target::cpvkdo(kld,IMDPSO,0:nwpc)
  real(8),target::cpebdep(kld,BVIMDPS+ndep,0:nwpc)
  pvkdo  =>cpvkdo
  pebdep =>cpebdep

end subroutine F_SetOutput
subroutine F_SetProp(psis,pgis,cpvkdp,cipos8,cipos12)
  use varcommon_mod
  use propagat_mod
  IMPLICIT NONE
  type(spc_interg),target::psis(kl*jnthet,0:nwpc)
  type(geo_interg),target::pgis(kl*jnthet,0:nwpc)
  real(8),target::cpvkdp(kl,IMDPSP,0:nwpc)
  integer(4),target::cipos12(12,0:nwpc)
  type (propinf_type),target::cipos8(0:nwpc)

  ps_vs=>psis;
  pg_vs=>pgis;
  pvkdp  =>cpvkdp
  ipos12=>cipos12
  ipos8=>cipos8
  !ps_vs=0;pg_vs=0; ipos8=0;
end subroutine F_SetProp
#endif
