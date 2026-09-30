!#################################################################################################
!-------------------------------------------------------------------------------------------------

  module irrp_smpi_mod

#ifndef SERIAL_TEST
#ifdef GF13
  use mpi
#endif
#endif

  public

!-------------------------------------------------------------------------------------------------

#ifndef SERIAL_TEST
#ifndef GF13
  include "mpif.h"
#endif
  integer,allocatable::vi(:)
#else
  integer, parameter :: MPI_COMM_WORLD = 1140850688
#endif

  integer(8),private::irtcb,clkrate
  type mpipacket
    integer             :: bsize = 0
    integer             :: pos = 0
    integer             :: dsize = 0
    integer(4), pointer :: buf(:)
  end type mpipacket

!-------------------------------------------------------------------------------------------------

  contains

!-------------------------------------------------------------------------------------------------
  SUBROUTINE starttimer
    irtcb=wav_irtc(clkrate)
  end SUBROUTINE starttimer
  real*4 function Difftimer()
    Difftimer=(wav_irtc(clkrate)-irtcb)/float(clkrate)
  end function Difftimer
  integer*8 function wav_irtc(rate)
   integer(8), optional :: rate
   !----- local -----
   integer      :: count
   integer      :: count_rate
   integer      :: count_max
   integer,save :: last_count = -1
   integer(8),save :: count_offset = 0

   call system_clock(count=count,count_rate=count_rate, count_max=count_max)
   if ( present(rate) ) rate = count_rate
   wav_irtc  = count
   if ( last_count /= -1 )then
     if ( count < last_count ) count_offset = count_offset + count_max
   end if
   wav_irtc = wav_irtc  + count_offset
   last_count = count
  end function wav_irtc
		integer function iwalltime() 
      integer vals(8)
      call date_and_time(VALUES=vals)
      iwalltime=vals(5)*10000+vals(6)*100+vals(7)
		end function iwalltime

#ifndef SERIAL_TEST

  subroutine irrp_abort(file, line, msg)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    character(len=*), intent(in), optional :: msg
    integer :: it
    print*,"aa"//"bb"//file
    if(present(msg))then
      write(6,'(a,i4,2a )')'IRRP Error ',line,' ',msg
    else
      write(6,'(a,i4)')'IRRP Error ',line
    endif
    call MPI_FINALIZE(it); 
    stop
  end subroutine irrp_abort

#else

  subroutine irrp_abort(file, line, msg)
    character(len=*),           intent(in) :: file
    integer,                    intent(in) :: line
    character(len=*), optional, intent(in) :: msg
    integer :: it
    if(present(msg))then
      write(6,'(a,i,2a )')'IRRP Error ',line,' ',msg
    else
      write(6,'(a,i)')'IRRP Error ',line
    endif
    stop
  end subroutine irrp_abort

#endif

#ifndef SERIAL_TEST

!-------------------------------------------------------------------------------------------------
! Set dsize and pos equal to zero as initial of mpipacket
! If lsize is given, allocate memory for the buffer.

  subroutine InitMpiPacket(pk, lsize)
    type(mpipacket),     intent(inout) :: pk
    integer*4, optional, intent(in)    :: lsize
    if(present(lsize))then
      if(lsize<=0 .or. pk%bsize<lsize)then
        if(pk%bsize > 0)then
          deallocate(pk%buf) ! If buf is not empty, deallocate it first.
          pk%bsize = 0       ! And set bsize equal to zero.
        endif
        if(lsize>0)then     ! If lsize is given and greater than 0,
          pk%bsize = lsize  ! set it as bsize and use it to allocate buf.
          allocate(pk%buf((pk%bsize+4)/4))
        endif
      endif
    endif
    pk%dsize = 0; pk%pos = 0   ! Initialize dsize and pos for the begining of packet.
  end subroutine InitMpiPacket

! Get dsize before using it (i.e. send, isend, bcast, unpack, recv, irecv, ...)
  subroutine SetMpiPacketDSize(pk, dsize)
     type(mpipacket),     intent(inout) :: pk
     integer*4, optional, intent(in)    :: dsize
     if(present(dsize))then
       pk%dsize = dsize
     else
       pk%dsize = pk%pos
     endif
     pk%pos = 0
  end subroutine SetMpiPacketDSize

  subroutine bcast_packet(pk, root, pid, mpicomm)
    type(mpipacket), intent(inout) :: pk
    integer, intent(in) :: root, pid, mpicomm
    integer :: lsize, ierr
    !if(pid == root)then
    call SetMpiPacketDSize(pk)
    lsize = pk%dsize;
    !endif
    call MPI_BCAST(lsize, 1, MPI_INTEGER4, root, mpicomm, ierr)
    if(pid /= root)then
      call InitMpiPacket(pk, lsize+100)
    endif
    call MPI_BCAST(pk%buf, lsize, MPI_PACKED, root, mpicomm, ierr)
    !if(pid /= root)then
    call SetMpiPacketDSize(pk, lsize)
    !endif
  end subroutine bcast_packet

#endif

!-------------------------------------------------------------------------------------------------

#ifndef SERIAL_TEST

  subroutine next(flg, pid, npe, mpicomm)
    integer, intent(in) :: flg, pid, npe, mpicomm
    integer ierr,v
    call flush(6)
    if(flg==2)then
        if(pid==npe-1)call MPI_SEND(v,1,MPI_INTEGER,0,1000,mpicomm,ierr)
    else
      if(pid<npe-1)then
        call MPI_SEND(v,1,MPI_INTEGER,pid+1,1000,mpicomm,ierr)
      else if(flg/=0)then
        call MPI_SEND(v,1,MPI_INTEGER,0,1000,mpicomm,ierr)
      endif
    endif
  end subroutine next

  subroutine prev(flg, pid, npe, mpicomm)
    integer, intent(in) :: flg, pid, npe, mpicomm
    integer ierr,v
    integer :: status(MPI_STATUS_SIZE)  ! mpi status info
    if(flg==2)then
      if(pid==0)call MPI_RECV(v,1,MPI_INTEGER,npe-1,1000,mpicomm,status,ierr)
    else
      if(pid>0)then
        call MPI_RECV(v,1,MPI_INTEGER,pid-1,1000,mpicomm,status,ierr)
      else if(flg/=0)then
        call MPI_RECV(v,1,MPI_INTEGER,npe-1,1000,mpicomm,status,ierr)
      endif
    endif
  end subroutine prev

  subroutine prevnext(flg, pid, npe, mpicomm)
    integer, intent(in) :: flg, pid, npe, mpicomm
    call prev(flg, pid, npe, mpicomm)
    call Next(flg, pid, npe, mpicomm)
  end subroutine prevnext

#endif

!-------------------------------------------------------------------------------------------------

#ifndef SERIAL_TEST

  subroutine barr(pid, npe, mpicomm)
    integer, intent(in) :: pid, npe, mpicomm
    integer it,ierr
		allocate(vi(npe))
   	call mpi_Gather(it,1,mpi_integer,vi,1,mpi_integer,0,mpicomm,ierr)
		deallocate(vi)
  end subroutine barr

#endif

!-------------------------------------------------------------------------------------------------

  end module irrp_smpi_mod

!-------------------------------------------------------------------------------------------------
!#################################################################################################
