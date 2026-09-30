#include "wavedeff.h"
#ifdef MTHREAD
# include "mthreadf.h"
#else
# define TPT() 
#endif
#ifdef DBGINF
# undef DBGINF
#endif
#ifdef DEBUG
# define CKEE(et) call checkee(et,nwpa+1,__FILE__,__LINE__)
# ifdef NO_MPI
#  define DBGINF   write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
# else
#  define DBGINF   call wav_mpi_barrier("AA");write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
# endif
#else
# define DBGINF
# define CKEE(e)
#endif
# undef DBGINF
# define DBGINF   ! write(*,'(a,3i6,";",18i8)')"MMM",__LINE__,mpi_id,RFB,RFG
!# define DBGINF0   !if(mpi_id==0)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
Module wavemdl_mod
  use varcommon_mod
  use implsch_mod
  use propagat_mod
  use boundary_mod
  use windin_mod
  use output_mod
  use output_cal_mod
  use output_nc_mod
  use partition_mod
  use restart_mod
  use waveinit_mod
  use platform_init_mod
  IMPLICIT NONE
  private
  public::RunWaveMdl
  integer ::curexgid=0
  real(8),external::dclktime
contains
  ! Main program
#ifdef ALLLOG
#define RINFO(TAG) write(*,'(i8," ",a,f8.3," ",2i8)')iwalltime(),TAG,Difftimer(2),mpi_id,mpi_npe;  call flush(6);
#else
#define RINFO(TAG) if(mpi_id==0)then;write(*,'(i8," ",a,f8.3," ",i8)')iwalltime(),TAG,Difftimer(2),mpi_npe;  call flush(6);endif
#endif
  SUBROUTINE exg_boundary
    curexgid=2-(curexgid-1) ! 1=>2 2=>1
#ifndef NO_MPI
    call exchange_boundary_start(eec,curexgid);
    call exchange_boundary_end()
#endif
  END SUBROUTINE 

#if defined(MMTHREAD)
  SUBROUTINE RunWaveMdl !{ !MMTHREAD
    integer key,ith,ko,iid,srcid,i
    real(8),allocatable:: irtcfs(:)
    real(8) rttt,t0,td
    integer RFB,RFG
    integer n,dnwp,cnwp,nnwp,exg_state,nwpb,ncheck
    call InitWaveMdl   ! wavemdl Init
    call Resettimer(-1)
    RINFO('Main A')
    RFB=1; 
    MSS(RFB);MWS(RFB)
    RFB=3
    NMWS(RFB,20);
    MSS(RFB);
    curexgid=1;
    RINFO('Main B')
    ! Main Loop
    BT(1)
    !wait threads ready
    ! W 1   25 L  r4 25
    ! S 1/2 25 L  15 25
    !TW 1   25 L r15 25
    !TS 1   25 L   4 25
    RFB=9; MSS(RFB);MWS(RFB);
    RINFO('Main C')
#ifndef NO_MPI
    if(mpi_npe>1)then
      RFG=1; call c_checkboundary(RFG) 
      call c_sendboundary(RFG)
    endif
#endif
    ! set modeltime
    ! RFB+1,for same RFB
    ! checkboundary data,with sync 
    ! setflag ,grp last loop check flag end 
    ! exchange ,must after checkboundary 
    ! setboundary
    ! RFB++,MSS(RFB),
    ! RFB++
    ! getnew wind,must after grp copy wind end
    ! check hist,rest,stop
    !if(mpi_id==0)
    RFB=RFB+1; MWS(RFB);MSS(RFB);
    RINFO('Main loop begin')
    t0=dclktime()

    DO   !{
      BT(2);
      if(RFB>MAXRFB )then
        RFB=10; MWSR(RFB);MSS(RFB);
      endif
      runstate%nTimeStep=runstate%nTimeStep+1
      !set Model time ,not
      call SetModelTime(RFB)
      RFB=RFB+1;MSS(RFB); !S1 for same flag,grp Setwind OK to sub,mainthread read wind at init & loop end
#ifndef NO_MPI
      ! wait data ok
      if(mpi_npe>1)then
        BT(14);RFG=RFG+1;call c_checkboundary(RFG) ;!last step
        ET(14)
      else
        MWS(RFB); ! no_mpi 
      endif
#else
      MWS(RFB); ! no_mpi 
#endif
      !set hist stop rest flags,must after check data
      call setflag(runstate%NextOutStep,runstate%NextOutStep,runstate%hist_eot,runstate%stop_now,runstate%rest_eot)
#ifndef NO_MPI
      if(mpi_npe>1) then
        BT( 4); call exg_boundary; ET(4) ! exchange data in MPI
      endif
#endif
      if(mpi_npe>1)then
        BT(15);call c_sendboundary(RFG) ;ET(15); !//last step//fixme:RFG
      endif
      RFB=RFB+1; MSS(RFB); !S2 !set boundary ok
      if(runstate%iPreCalT>0)then
        call DecIPrecalT;cycle
      endif
      !get next wind ,double buf,not need wait
      BT( 12);key=IDataIO(4); ET(12)
      RFB=RFB+1; !S3 grp wait subthread compute end
      if(runstate%nTimeStep==runstate%NextOutStep)then ! check output have wait
        BT(13);call CheckOutPut(-1,-1,RFB,eec,1,0,1,0,runstate%hist_eot); ET(13)
      endif
      if(runstate%rest_eot/=0)then
        !BT(14);call CheckRestart; ET(14)
        runstate%rest_eot=0
        call setflag(runstate%NextOutStep,runstate%NextOutStep,runstate%hist_eot,runstate%stop_now,runstate%rest_eot)
      endif
      ET(2)
#ifdef DTIMES
      if(runstate%nTimeStep==0) PT();
#endif
      if(runstate%stop_now/=0)exit
      !RFB=RFB+1; MWS(RFB);MSS(RFB); !S2 !set boundary ok
    end do !}
    td=dclktime()
    if(mpi_id==0)print*,"Main wait end ",RFB,td-t0,runstate%mTimeStep
    RFB=RFB+1000;MSS(RFB);
    RFB=RFB-10; MWS(RFB);
    PT();
    call EndWaveMdl
  end SUBROUTINE RunWaveMdl !}
#elif defined(MTHREAD)
  SUBROUTINE RunWaveMdl !{ !MTHEAD 
    integer key,ith,ko,iid,srcid,i
    real(8),allocatable:: irtcfs(:)
    real(8) rttt,t0,td
    integer RFB,RFG
    integer n,dnwp,cnwp,nnwp,exg_state,nwpb,ncheck
    call InitWaveMdl   ! wavemdl Init
    RINFO('Main A0')
    call Resettimer(-1)
    RINFO('Main A')
    RFB=1; 
    MSS(RFB);MWS(RFB)
    RFB=3
    NMWS(RFB,20);
    MSS(RFB);
    curexgid=1;
    RINFO('Main B')
    ! Main Loop
    BT(1)
    RFB=9; MSS(RFB);MWS(RFB);
    RFB=RFB+1;
    RINFO('Main C')
    ! set modeltime
    ! RFB+1,for same RFB
    ! checkboundary data,with sync 
    ! setflag ,grp last loop check flag end 
    ! exchange ,must after checkboundary 
    ! setboundary
    ! RFB++,MSS(RFB),
    ! RFB++
    ! getnew wind,must after grp copy wind end
    ! check hist,rest,stop
    !if(mpi_id==0)
    RINFO('Main loop begin')
    t0=dclktime()
    call testmpi(0)
    MSS(RFB);
    DO   !{
      BT(2);
      if(RFB>MAXRFB )then
        RFB=10; MWSR(RFB);MSS(RFB);
      endif
      runstate%nTimeStep=runstate%nTimeStep+1
      !set Model time ,not
      call SetModelTime(RFB)
      RFB=RFB+1;MSS(RFB); !S1 for same flag,grp Setwind OK to sub,mainthread read wind at init & loop end
      MWS(RFB); ! 1 wait sub OK
      !set hist stop rest flags,must after check data
      call setflag(runstate%NextOutStep,runstate%hist_eot,runstate%stop_now,runstate%rest_eot)
#ifndef NO_MPI
      if(mpi_npe>1) then
        ! exchange data in MPI
        BT( 4); call exg_boundary; ET(4)
      endif
#endif
      RFB=RFB+1; MSS(RFB); !S2 !set boundary ok
      if(runstate%iPreCalT>0)then
        call DecIPrecalT;cycle
      endif
      !call ppp(eet,eec);
      !RFB=RFB+1;MWS(RFB); call outhss(runstate%nTimeStep+1000,eec); MSS(RFB);!O0
      RFB=RFB+1;MWS(RFB); MSS(RFB);!O1
      !RFB=RFB+1;MWS(RFB); call outhss(runstate%nTimeStep+2000,eec); MSS(RFB);!O1
      !RFB=RFB+1;MWS(RFB); call outhss(runstate%nTimeStep+3000,eec); MSS(RFB);!O2
      !RFB=RFB+1;MWS(RFB); call outhss(runstate%nTimeStep+4000,eet); MSS(RFB);!O3
      !RFB=RFB+1;MWS(RFB); call outhss(runstate%nTimeStep+5000,eec); MSS(RFB);!O4
      !RFB=RFB+1;MWS(RFB); call outhss(runstate%nTimeStep+6000,eec); MSS(RFB);!O5
      !stop
      !get next wind ,double buf,not need wait
      BT( 12);key=IDataIO(4); ET(12)
      !MWS(RFB); ! 3 hist have sync
      RFB=RFB+1; !S3 grp wait subthread compute end
      if(runstate%nTimeStep==runstate%NextOutStep)then ! check output have wait
        BT(13);call CheckOutPut(-1,-1,RFB,eec,1,0,1,0,runstate%hist_eot); ET(13)
      endif
      if(runstate%rest_eot/=0)then
        !BT(14);call CheckRestart; ET(14)
        runstate%rest_eot=0
        call setflag(runstate%NextOutStep,runstate%hist_eot,runstate%stop_now,runstate%rest_eot)
      endif
      ET(2)
#ifdef DTIMES
      if(runstate%nTimeStep==0) PT();
#endif
      if(runstate%stop_now/=0)exit
    end do !}
    td=dclktime()
    if(mpi_id==0)print*,"Main wait end ",RFB,td-t0,runstate%mTimeStep
    RFB=RFB+1000;MSS(RFB);
    RFB=RFB-10; MWS(RFB);
    PT();
    call EndWaveMdl
  end SUBROUTINE RunWaveMdl !}
#else
  SUBROUTINE RunWaveMdl !{ !NO MTHREAD
    integer key,ith,ko,iid,srcid,i
    real(8),allocatable:: irtcfs(:)
    real(8) rttt,t0,td
    integer RFB,RFG
    integer n,dnwp,cnwp,nnwp,exg_state,nwpb,ncheck
    call InitWaveMdl   ! wavemdl Init
#ifdef USEXNL
    call init_xnl
#else
    call init_dia
#endif
    RINFO('Main A0')
    call Resettimer(-1)
    RINFO('Main A')
    curexgid=1;
    ! Main Loop
    BT(1)
    RFB=10;
    call testmpi(0)
    call setspec(0,eec,1,nwpc,1)   ;CKEE(eec) !set ee !      set water boundary
    !print*,"AAAAAAAAAAA",runstate%nTimeStep
    DO   !{
      BT(2);
      runstate%nTimeStep=runstate%nTimeStep+1
      !set Model time ,not
      call SetModelTime(RFB)
      !set wind ok
      !z call setflag(runstate%NextOutStep,runstate%hist_eot,runstate%stop_now,runstate%rest_eot)
      ! exchange data in MPI
#ifndef NO_MPI
      BT( 4); call exg_boundary; ET(4)
#endif
      CKEE(eec)
      !CALZ NWPC*(199*NK*NJ+6*NJ+12*NK +28 +0.1*NK*NJ div 8 sqrt 3 exp 2 tanh 1)
      !CALZ NWPC*(76938 sqrt 3 exp 2 tanh 1)
      !CALA 36*NK*NJ
      !call prtee(1,eec(1,1,1),eet(1,1,1));
      CKEE(eec)
      !call ppp(eet,eec);
      !call outhss(runstate%nTimeStep+1000,eec);!O0
      BT( 5);call propagats_spec(eet,eec,1,nwpc);CKEE(eet)
      !call outhss(runstate%nTimeStep+2000,eec);!O0
      EBT(6);call propagats_geo (eec,eet,1,nwpc);CKEE(eec)
      !call outhss(runstate%nTimeStep+3000,eec);!O0
      EBT(7);call smooth_e(eet,eec,1,nwpc);CKEE(eet)
      !call outhss(runstate%nTimeStep+4000,eec);!O0
      EBT(8);call implschs(eec,eet,1,nwpc);CKEE(eec)
      !call outhss(runstate%nTimeStep+5000,eec);!O0
      EBT(9);call setspec(0,eec,1,nwpc,2)   ;CKEE(eec) !set ee !      set water boundary
      !call outhss(runstate%nTimeStep+6000,eec);!O5 
      !stop
      EBT(10);!call Monitor(eec,-1,0);
      if(runstate%iPreCalT>0)then
        call DecIPrecalT;cycle
      endif
      EBT(11);call accumea(-1,eec,1,nwpc) ;
      !get next wind ,
      EBT( 12);key=IDataIO(4);
      ET(12)
      if(runstate%nTimeStep==runstate%NextOutStep)then ! check output have wait
        BT(13);call CheckOutPut(-1,-1,RFB,eec,1,nwpc,1,0,runstate%hist_eot); ET(13)
      endif
      if(runstate%rest_eot/=0)then
        !BT(14);call CheckRestart; ET(14)
      endif
      ET(2)
      if(runstate%stop_now/=0)exit
    end do !}
    call EndWaveMdl
  end SUBROUTINE RunWaveMdl !}
#endif
end Module wavemdl_mod
!===============================
! called by C code for oprate Fortran variable
#ifdef MTHREAD
subroutine accumea(ind,ee,nwpb,nwpe)
  use varcommon_mod
  IMPLICIT NONE
  REALD ,intent(in):: ee(kl,jnthet,0:nwpc)
  integer,intent(in)::nwpb,nwpe,ind
  integer ko,iac
  type (AVHIST),pointer::ph
  do ko=1,Navhists
    ph=>avhists(ko)
    if(ph%mean/=0)then
      if(ind>0)then
        call addea(ph%ea(:,:,iac),ee,nwpb,nwpe);
      else
        ph%nea=ph%nea+1
      endif
    endif
  enddo
end  subroutine accumea
#define CKINF(tag,ind,rfb)  !call ckinf(tag,ind,rfb);
subroutine CheckOutPut(igrp,ind,RFBt,ee,iacb,iace,iacb2,iace2,nout)
  use varcommon_mod
  IMPLICIT NONE
  integer,intent(in):: iacb,iace ,nout
  integer,intent(in):: igrp,ind,iacb2,iace2
  integer,intent(inout):: RFBt
  REALD ,target,intent(in):: ee(kl,jnthet,0:nwpc)
  integer ko,ahist,newfile,ntpf
  integer noutc,RFB
  real(8) ::timet
  REALD ,pointer::ea(:,:,:)
  type (AVHIST),pointer::ph
  timet=runstate%edaycur+runstate%timecur+deltts/10/86400.
  noutc=nout
  RFB=RFBt;
  do ko=1,Navhists
    ph=>avhists(ko)
    !if(ph%nextTime<=timet)then
      !if(ph%hist_ect/=0)then
      RFB=RFB+1
      if(ind>=0)then
#ifdef MMTHREAD
        if(ind==0)then
          GWS(RFB) ; 
          SSM(RFB)
        else
          SSG(RFB);
        endif
#else
        SSM(RFB)
#endif
        if(ph%mean/=0)then
          ea=>ph%ea
          call AverageEA(ea,ph%nea,iacb,iace)
          call AverageEA(ea,ph%nea,iacb2,iace2)
        else
          ea=>ee
        endif
        IF(ph%wrwa>0)THEN
          call mean1(ea,iacb ,iace ,ph%ape,ph%tpf,ph%aet,ph%h1_3) 
          call mean1(ea,iacb2,iace2,ph%ape,ph%tpf,ph%aet,ph%h1_3)
        endif
        IF(ph%wrbv>0)THEN
          call mixture(ea,iacb ,iace ,ph%bv)
          call mixture(ea,iacb2,iace2,ph%bv)
        endif
        RFB=RFB+1
#ifdef MMTHREAD
        if(ind==0)then
          GWS(RFB) ;!l2 
          CKINF(125*10+ko,ind,rfb);
          SSM(RFB);!l3
          SWM(RFB);!l6
          GSS(RFB);!l7
        else
          CKINF(124*10+ko,ind,rfb);
          SSG(RFB);!l1
          SWG(RFB);!l8
        endif
#else
        SSM(RFB);SWM(RFB);
#endif
        if(ph%mean/=0)then
          call ZeroEA(ea,iacb,iace) !clear
          call ZeroEA(ea,iacb2,iace2) !clear
        endif
      else
        MSS(RFB); MWS(RFB);
        RFB=RFB+1;
        runstate%hist_eot=0 ! reset it,
        call setflag(runstate%NextOutStep, runstate%hist_eot,runstate%stop_now,runstate%rest_eot)
        MWS(RFB);!l4
        MSS(RFB);!l5 
        CKINF(127*10+ko,ind,rfb);
        if(mpi_id==0)then;
          write(*,'(i8," ",a,f8.3," ",i8)')iwalltime(),"Output",Difftimer(2),runstate%nTimeStep;  
          call flush(6);
        endif
        call HistOutPut(ph)
        if(mpi_id==0)then;
          write(*,'(i8," ",a,f8.3," ",i8)')iwalltime(),"Outputend",Difftimer(2),runstate%nTimeStep;  
          call flush(6);
        endif
        ph%nea=0;
        ph%hist_ect=0
      endif
      noutc=noutc-1
      if(noutc<=0)exit
    !endif
  enddo
  RFBt=RFB;
end subroutine 
#endif
