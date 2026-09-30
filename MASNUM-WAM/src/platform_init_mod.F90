#include "wavedeff.h"
#ifdef MTHREAD
#include "mthreadf.h"
#endif
module platform_init_mod
#ifndef NO_MPI
  use wav_mpi_mod
#endif
  use varcommon_mod
  use partition_mod
  use windin_mod
  use output_cal_mod
  use output_mod
  use restart_mod
  use propagat_mod
  use boundary_mod
  use implsch_mod
  !use,intrinsic::iso_c_binding
  !use wavemdl_mod
  IMPLICIT NONE

  public::Init_Distinct,InitAfterMpi
  real(8),external::c_getfpcd
  ! 必须与std.h mthread.h 中相同
  type CPUINFO
    integer NCorePClu   ;! 每簇(NUMA)核数(38)
    integer MCorePClu   ;! 每簇(NUMA)可用核数(38),假设不可用核在最后
    integer NCluPNode   ;! 每节点簇数(16)
    integer NManageCore ;! 管理核核数，编号>NCorePClu*NCluPNode    
    integer NThPGrp     ;! 每组启动线程数(38)
    integer NGrpPProc   ;! 每进程组数(16)
    integer NProcPNode  ;! 每节点进程数(38)
    integer ManageCoreId;! 管理核(-1)
    integer OffClu      ;! 第一组使用的簇(0)
    integer SkipClu     ;! 几个簇启动1个组间(1)
    integer OffCore     ;! 每个组使用的第一个核(0)
    integer SkipCore    ;! 几个核启动一个线程(1)
  end type CPUINFO
contains
#ifdef MMTHREAD
#define RMODE 2
#elif defined(MTHREAD)
#define RMODE 1
#else
#define RMODE 0
#endif
  SUBROUTINE   InitAfterMpi
    integer ::ierr=0;
  type(CPUINFO) ::cpuinf
    if(NProcPNode <0)NProcPNode = NThPGrp*NGrpPProc;   
    if(NProcPNode >MCorePClu*NCluPNode)then
      print*,'NProcPNode ',NProcPNode ,"> MCorePClu*NCluPNode",MCorePClu,NCluPNode
      stop
    endif
#define VD(v) cpuinf%v=v
    VD(NCorePClu   );
    VD(MCorePClu   );
    VD(NCluPNode   );
    VD(NManageCore );
    VD(NThPGrp     );
    VD(NGrpPProc   );
    VD(NProcPNode  );
    VD(ManageCoreId);
    VD(OffClu      );
    VD(SkipClu     );
    VD(OffCore     );
    VD(SkipCore    );
#undef VD

#ifndef NO_MPI
    call init_cpu(RMODE,mpi_comm_wav,mpi_id,mpi_npe,cpuinf,ierr)
#else
    call init_cpu(RMODE,-1,mpi_id,mpi_npe,cpuinf,ierr)
#endif
#ifdef MTHREAD
    call InitThreads(mpi_id,cpuinf,ierr)
#endif
#ifndef NO_MPI
    if(ierr<0) call wav_mpi_abort("Init_cpu")
#else
    if(ierr<0)stop
#endif
    ManageCoreId=cpuinf%ManageCoreId
  end SUBROUTINE   InitAfterMpi
  SUBROUTINE End_Distinct
#ifdef MTHREAD
    call endthreads
#endif
  end SUBROUTINE End_Distinct
  real(8) function GetNfpc(nt)
    integer i,nt
  end function GetNfpc

  SUBROUTINE Init_Distinct
#ifndef USE_C_ALLOC
    call c_setpointers(eec,eet,ip,ps_vs,pg_vs,ipos12,ipos8,wxy,pvkdo,pebdep)
#endif
    call c_Init_Distinct
    !print*,mpi_id,"Start Threads=============";call flush(6);
#ifdef MTHREAD
    call startthreads
#endif
  end SUBROUTINE Init_Distinct
end module platform_init_mod
