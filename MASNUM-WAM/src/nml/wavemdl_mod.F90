#include "wavedeff.h"
#ifdef DBGINF
#undef DBGINF
#endif
#ifdef DEBUG
# define CKEE(et) call checkee(et,nwpa+1,__FILE__,__LINE__)
#ifdef NO_MPI
# define DBGINF   write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
#else
# define DBGINF   call wav_mpi_barrier("AA");write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
#endif
#else
# define DBGINF
# define CKEE(e) !call checkee(e,nwpa+1,__FILE__,__LINE__)
#endif
!# define DBGINF0   !if(mpi_id==0)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
Module wavemdl_mod
  use varcommon_mod
  use implsch_mod
  use propagat_mod
  use boundary_mod
  use windin_mod
  use output_mod
  use output_cal_mod
  use partition_mod
  use restart_mod
  use waveinit_mod
  use platform_init_mod
  IMPLICIT NONE
  private
  public::RunWaveMdl
  integer ::curexgid=0
contains
  ! Main program
  SUBROUTINE RunWaveMdl !{
    integer key,ith,ko,iid,srcid,i
    real(8),allocatable:: irtcfs(:)
    real(8) rttt
    integer n,dnwp,cnwp,nnwp,exg_state,nwpb,ncheck
    !program Init
    call InitWaveMdl   ! wavemdl Init
    call Resettimer
    curexgid=1;
    ! Main Loop
    BT(1)
    DO   !{
      BT(2);
      runstate%nTimeStep=runstate%nTimeStep+1
      !set Model time ,not
      call SetModelTime
      curexgid=2-(curexgid-1) ! 1=>2 2=>1
      ! exchange data in MPI
#ifndef NO_MPI
      BT( 4);
      call exchange_boundary_start(eec,curexgid);
      call exchange_boundary_end()
      ET(4)
#endif
      CKEE(eec)
      !CALZ NWPC*(199*NK*NJ+6*NJ+12*NK +28 +0.1*NK*NJ div 8 sqrt 3 exp 2 tanh 1)
      !CALZ NWPC*(76938 sqrt 3 exp 2 tanh 1)
      !CALA 36*NK*NJ
      !call prtee(1,eec(1,1,1),eet(1,1,1));
      CKEE(eec)
      BT( 5);call propagats_spec(eec,eet,1,nwpc);CKEE(eet)
      EBT(6);call propagats_geo (eec,eet,1,nwpc);CKEE(eec)
      EBT(7);call smooth_e(eet,eec,1,nwpc);CKEE(eet)
      !CALA pm NWPC*(6*NJ+12*NK+163*NK*NJ+28+0.1*NK*NJ div 8 sqrt 3 exp 2 tanh 1)
      !call prtee(2,eet(1,1,1),eec(1,1,1));
      EBT(8);call implschs(eec,eet,1,nwpc);CKEE(eec)
      EBT(9);call setspec(0,eec,1,nwpc,2)   ;CKEE(eec) !set ee !      set water boundary
      EBT(10);!call Monitor(eec,-1,0);
      if(runstate%iPreCalT>0)then
        call DecIPrecalT;cycle
      endif
      !for output ,
      EBT(11);call ACCUMEA(eec,1,nwpc) ;
      !get next wind ,
      EBT( 12);key=IDataIO(4);
      ET(12)
      if(runstate%hist_eot>0)then
        !BT(13);call CheckOutPut(eec,1,nwpc,runstate%hist_eot); ET(13)
      endif
      if(runstate%rest_eot/=0)then
        !BT(14);call CheckRestart; ET(14)
      endif
      ET(2)
      if(runstate%stop_now/=0)exit
    end do !}
    call EndWaveMdl
  end SUBROUTINE RunWaveMdl !}
end Module wavemdl_mod
