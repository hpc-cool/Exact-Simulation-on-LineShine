#include "wavedeff.h"
#ifdef DEBUG
!#  define DBGINF
#  define CKEE(e)
#else
!#  define DBGINF   write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
#  define CKEE(e)
#endif
!# define DBGINF0   !if(mpi_id==0)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
Module WaveMdl_mod
  use varcommon_mod
  use partition_mod
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
    call InitWaveMdl   ! wavemdl Init
    DBGO0(4,'(a)')trim(startType)
    call exchange_boundary_start(eec,curexgid)
    call exchange_boundary_end()       
    call EndWaveMdl
  end SUBROUTINE RunWaveMdl !}

end Module WaveMdl_mod
