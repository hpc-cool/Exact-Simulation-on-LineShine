Module debughlp_mod
#ifndef NO_MPI
  use wav_mpi_mod
#endif
  IMPLICIT NONE
public
  integer,external::iwalltime
  real(8),external::Difftimer,gettimer
  integer ::OutDType,NTYPE
  integer::dbglvl=10
#ifdef NO_MPI
   integer,public:: mpi_comm_wav=-1,mpi_id=0,mpi_npe=1
#endif
#define MAXTIMER 30
  real*8:: rtcb(0:MAXTIMER+1),rtcl(0:MAXTIMER+1)
!#define HAVE_CPU_TIME
#define HAVE_SYSTEM_CLOCK
contains
#  define DBGINF   write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__
!==========================================================================
  subroutine wav_stdio()
#ifdef USE_SHR_MODS
      !call shr_msg_chdir   ('wav') ! changes cwd
      call shr_msg_stdio('wav')
#endif
  end subroutine wav_stdio
  SUBROUTINE wav_chStdOut
    character*256 pid
    call flush6
    if(mpi_id>=0)then
      close(6)
      if(mpi_id==0)then
        write(pid,'("wav.log.",i6.6)')mpi_npe
        open(unit=6,file=trim(pid),status='unknown')
        write(6,*)mpi_id,'Redirect Out to ',trim(pid)
      else
        if(DBGLvl<=5)then
          DBGLvl=0
        else
          ! write(pid,'("logs/log",i6.6,"/wav.log.",i6.6)')mpi_npe,mpi_id
          write(pid,'("wav.log.",i6.6)')mpi_id
          open(unit=6,file=trim(pid),status='unknown')
          write(6,*)mpi_id,'Redirect Out to ',trim(pid)
        endif
      endif
      call flush6
    endif
  END SUBROUTINE wav_chStdOut

  subroutine wav_abort(msg)
    character*(*) msg
#ifndef NO_MPI
    call wav_mpi_abort(msg,-1)
#else
    write(*,*)msg
    stop
#endif
  end subroutine wav_abort

  subroutine WaitPrev
    integer it
#ifndef NO_MPI
    if(mpi_id>0)CALL wav_mpi_RECV(it,mpi_id-1,100)
#endif
  end subroutine WaitPrev
  subroutine NEXTCPU
    integer it
#ifndef NO_MPI
    if(mpi_id<mpi_npe-1)call wav_mpi_send(it,mpi_id+1,100)
#endif
  end subroutine NEXTCPU
  subroutine Zbarrier
    integer it
    integer,allocatable,save::vi(:)
#ifndef NO_MPI
  if(.not.allocated(vi))allocate(vi(mpi_npe))
  call wav_mpi_gather(it,vi,0)
  call wav_mpi_scatter(vi,it,0)
!  deallocate(vi)
#endif
  end subroutine Zbarrier
  subroutine flush6
    call flush (6)
    !call flush (0)
  end subroutine flush6
  subroutine flusht(ift)
  integer ift
    call flush (ift)
  end subroutine flusht
end Module debughlp_mod
