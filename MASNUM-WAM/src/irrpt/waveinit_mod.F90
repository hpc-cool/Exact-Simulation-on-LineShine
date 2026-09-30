#include "wavedeff.h"
!#define ALLLOG
!# define DBGINF   !call wav_mpi_barrier("AA");write(*,'(a,i6,i4.4,";",18i8)')__FILE__,__LINE__,mpi_id

Module waveinit_mod
#ifndef NO_MPI
  use wav_mpi_mod
#endif
  use varcommon_mod
  use partition_mod
  use platform_init_mod
IMPLICIT NONE
private
  integer::hist,rest
  public::InitWaveMdl,EndWaveMdl,DecIPrecalT,SetModelTime
  integer ::ndstep,oldday=-1
  real(8)  dfpc
contains
#define RINFO(TAG) if(mpi_id==0)write(*,'(i8," ",a,f8.3," ",i8)')iwalltime(),TAG,Difftimer(2),mpi_npe;  call flush(6);
#define DBGA call wav_mpi_barrier;write(*,'("DBGA",i4,":",i4.3)')__LINE__,mpi_id;  call flush(6);
#define DBGB write(*,'("DBGA",i4,":",4i4.3)')__LINE__,mpi_id,mpi_npe;  call flush(6);
  SUBROUTINE InitWaveMdl
    integer key,ir,out
    real*8 tt
    !ir=ieee_flags('clear', 'exception', 'all',out )
    call Resettimer
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
    RINFO('gathercheck End')
#endif
    !call testtimes
    !call zftest;         ! call zftest;         ! stop
    Call InitPara      !varcommon_mod Init
    call InitAfterMpi
    DBGLVL=20
    ! call wav_chStdOut
    if(mpi_id==0) then
      RINFO('InitPara End')
    else
      if(DBGLvl<5)DBGLvl=0
    endif
    call setwave
    RINFO('setwave End')
    !call nlweight;    stop
    RINFO('ALLOCATE depg,nspg begin')
    ALLOCATE(depg(gixl,giyl),nspg(gixl,giyl));
    RINFO('ALLOCATE depg,nspg End')
    call ReadTopog
    RINFO('ReadTopog End')
#ifndef NO_MPI
    call gathercheck(__LINE__)
    RINFO('gathercheck End')
#endif
    !partition !!!!!
    call partition
    RINFO('partition End')
    !call ReadTopog;    RINFO('ReadTopog End')
    !call InitMonitor
    !RINFO('InitMonitor End')
    if(mpi_id==0)then
      DBGO(0,*) 'KL       = ',kl
      DBGO(0,*) 'DBGLVL   = ',DBGLVL
      DBGO(0,*) 'DELTTS   = ',DELTTS
      DBGO(0,*) 'GXLON0   = ',GXLON0
      DBGO(0,*) 'GYLAT0   = ',GYLAT0
      DBGO(0,*) 'GXLON1   = ',GXLON1
      DBGO(0,*) 'GYLAT1   = ',GYLAT1
      DBGO(0,*) 'GRDSZX   = ',GRDSZX
      DBGO(0,*) 'GRDSZY   = ',GRDSZY
      DBGO(0,*) 'GCIRCLE  = ',GCIRCLE
      DBGO(0,*) 'DEPTHMOD = ',DEPTHMOD
      DBGO(0,*) 'LOGSCURR = ',LOGSCURR
      DBGO(0,*) 'LOGSDISS = ',LOGSDISS
      DBGO(0,*) 'nwpa = ',nwpa,nwpc,nwps,gnwpc,(gnwpc+0.)/mpi_npe,gixl,giyl
      DBGO(4,'(10i9)')nwpcs
      call flush(6);
    endif
    ! sum case Needed allocate Mem use C Program for mem align
    ! if Not ,InitAfterPartition  is null
    !print*,__FILE__,__LINE__,mpi_id
    !call wav_mpi_barrier();
    RINFO('InitAfterPartition start')
    call flush6

    call InitAfterPartition
    RINFO('InitAfterPartition')
    call set_neighbor
    RINFO('set_neighbor')
  end SUBROUTINE InitWaveMdl !}

  SUBROUTINE DecIPrecalT
  end SUBROUTINE DecIPrecalT

  SUBROUTINE EndWaveMdl
    integer ir
    call End_Distinct
    call InitOutPut(-1)  ! Close  OutPut
    ir=IDataIO(-1)
    call InitMpi(-1)      ! Mpi Deinit
    DBGO0(0,*) "end wave model A",runstate%nTimeStep
    call flush6
  end SUBROUTINE EndWaveMdl

  subroutine SetModelTime
  END subroutine SetModelTime
end Module waveinit_mod

! called by C code for oprate Fortran variable

