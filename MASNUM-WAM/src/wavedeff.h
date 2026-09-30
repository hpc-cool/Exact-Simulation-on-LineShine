#include "wavedef.h"
!  fortran code
#  define REALD real(8)
#  define DBGO0(lvl,format)  if(mpi_id==0.and.lvl<=DBGLvl)write(*,format)
#  define DBGO(lvl,format)  if(lvl<=DBGLvl)write(*,format)
! #define NODBGINFO
#  ifdef NODBGINFO
#   define OUT_ALLOCSIZE(v)   !
#   define DBGINFN(n) !
#   define DBGINF0  !
#   define DBGINF   !
#   define DBGINFA  !
#   define SDBGINF  !
#   define FFO      !
#   define FFO0     !
#   define DBG(lvl) !
#   define DBARRIER
#  else
#   define OUT_ALLOCSIZE(v) !write(6,'(a,a8,a,i8,a,18i6)')"ALLOCATE:size(","v",")=",size(v),__FILE__,__LINE__,mpi_id;
#   define DBGINFN(n)  if(mpi_id==n)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id,runstate%nTimeStep
#   define DBGINF0  if(mpi_id==0.and.DbgLvl>8)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id,runstate%nTimeStep
#   define DBGINFA0 if(mpi_id==0.and.DbgLvl>8)write(*,*)__FILE__,__LINE__,mpi_id,runstate%nTimeStep,";"
#   define DBGINF   if(DbgLvl>8)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id,runstate%nTimeStep
#   define DBGINFA  if(DbgLvl>8)write(*,*)__FILE__,__LINE__,mpi_id,runstate%nTimeStep,";"
#   define SDBGINF  if(DbgLvl>8)write(*,'(a,3i6,";",a,18i8)')__FILE__,__LINE__,mpi_id,runstate%nTimeStep
#   define FFO0     if(mpi_id==0)call flush6
#   define FFO      call flush6
#   define DBG(lvl) if(lvl<=DBGLvl)
#   define DBARRIER call wav_mpi_barrier
#  endif
#  ifndef MTHREAD
#   define TBT(l)  
#   define TET(m)  
#   define TEBT(l) 
#   define TPT()   
#  endif
#   define BT(l)  TBT(l) ;call starttimer( l);
#   define ET(m)  TET(m) ;call endtimer( m);
#   define EBT(l) TEBT(l);call endstarttimer(l-1, l);
#   define PT()   TPT()

