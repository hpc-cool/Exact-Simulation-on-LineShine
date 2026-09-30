#define MPI_GATHERV MYMPI_GATHERV
#define MPI_SCATTERV MYMPI_SCATTERV
Module wav_mpi_mod
! -------------------------------------------------------------------------------
! PURPOSE: general layer on MPI functions
! -------------------------------------------------------------------------------
   implicit none
   private
   integer::DBGLVL=10
! mpi library include file
   integer,public:: mpi_comm_wav=-1,mpi_id=0,mpi_npe=1
#ifndef NO_MPI
   integer,public:: NextProc=mpi_proc_null,PrevProc=mpi_proc_null
   public::mpi_any_source,mpi_real8,mpi_real4,MPI_STATUS_SIZE,MPI_INTEGER,MPI_INTEGER2
!  PUBLIC: Public interfaces
   character*(8),public:: mpi_ids
   character*(MPI_MAX_PROCESSOR_NAME),public:: NodeName
   public :: wav_mpi_chkerr
   public :: wav_mpi_bcast ,wav_mpi_gather,wav_mpi_gatherv
   public :: wav_mpi_scatter,wav_mpi_scatterv
   public :: wav_mpi_send  ,wav_mpi_recv
!   public :: wav_mpi_isend ,wav_mpi_irecv
   public :: wav_mpi_sum   ,wav_mpi_min  ,wav_mpi_max
   public :: wav_mpi_commsize
   public :: wav_mpi_commrank
   public :: wav_init_mpi
   public :: wav_mpi_init
   public :: wav_mpi_initialized
   public :: wav_mpi_abort
   public :: wav_mpi_barrier
   public :: wav_mpi_finalize
   public :: InitMpiPacket,wav_mpi_pack,wav_mpi_unpack,SetMpiPacketDSize
   public :: SerRun,SerRB,SerRE,SerRunr,SerRBr,SerREr,SerRunf,wav_mpi_check
   type,public:: mpipacket
    integer bsize,pos,dsize
    integer*4,allocatable::Buf(:)
   end type mpipacket
   interface wav_mpi_gather ; module procedure &
    wav_mpi_gatheri0, wav_mpi_gatheri1, &
    wav_mpi_gatherr0, wav_mpi_gatherr1, &
    wav_mpi_gatherd0, wav_mpi_gatherd1
   end interface
   interface wav_mpi_scatter ; module procedure &
    wav_mpi_scatteri0 , wav_mpi_scatteri1, &
    wav_mpi_scatteri80, wav_mpi_scatteri81, &
    wav_mpi_scatterr0 , wav_mpi_scatterr1, &
    wav_mpi_scatterd0 , wav_mpi_scatterd1
   end interface
   interface wav_mpi_gatherv ; module procedure &
    wav_mpi_gathervi1, wav_mpi_gathervi2, &
    wav_mpi_gathervr1, wav_mpi_gathervr2, &
    wav_mpi_gathervd1
   end interface
   interface wav_mpi_scatterv ; module procedure &
    wav_mpi_scattervi1, wav_mpi_scattervi2, &
    wav_mpi_scattervr1, wav_mpi_scattervr2, &
    wav_mpi_scattervd1
   end interface
   interface wav_mpi_send ; module procedure &
     wav_mpi_sendi0, wav_mpi_sendi1, &
     wav_mpi_sendr0, wav_mpi_sendr1
   end interface
   interface wav_mpi_recv ; module procedure &
     wav_mpi_recvi0, wav_mpi_recvi1, &
     wav_mpi_recvr0, wav_mpi_recvr1
   end interface
   interface wav_mpi_pack  ; module procedure &
     wav_mpi_packi0, wav_mpi_packi1, &
     wav_mpi_packr0, wav_mpi_packr1, &
     wav_mpi_packd0, wav_mpi_packd1, &
     wav_mpi_packc0, wav_mpi_packc1
   end interface
   interface wav_mpi_unpack; module procedure &
     wav_mpi_unpacki0, wav_mpi_unpacki1, &
     wav_mpi_unpackr0, wav_mpi_unpackr1, &
     wav_mpi_unpackd0, wav_mpi_unpackd1, &
     wav_mpi_unpackc0, wav_mpi_unpackc1
   end interface
   interface wav_mpi_bcast ; module procedure &
     wav_mpi_bcasts0 , wav_mpi_bcasts1, &
     wav_mpi_bcastl0 , wav_mpi_bcastb2, &
     wav_mpi_bcastw0 , wav_mpi_bcastw1, wav_mpi_bcastw2, &
     wav_mpi_bcasti0 , wav_mpi_bcasti1, wav_mpi_bcasti2, &
     wav_mpi_bcasti80, wav_mpi_bcasti81, &
     wav_mpi_bcastr0 , wav_mpi_bcastr1, wav_mpi_bcastr2 , &
     wav_mpi_bcastp
   end interface
   interface wav_mpi_sum ; module procedure &
     wav_mpi_sumi0, wav_mpi_sumi1, &
     wav_mpi_sumr0, wav_mpi_sumr1, wav_mpi_sumr2, wav_mpi_sumr3
   end interface
   interface wav_mpi_min ; module procedure &
     wav_mpi_mini0, wav_mpi_mini1, &
     wav_mpi_minr0, wav_mpi_minr1
   end interface
   interface wav_mpi_max ; module procedure &
     wav_mpi_maxi0, wav_mpi_maxi1, &
     wav_mpi_maxr0, wav_mpi_maxr1
   end interface
CONTAINS
! ===============================================================================
    subroutine wav_init_mpi(mode )
    integer mode
        integer ierr,nl
        logical Inited
! int COUPLE_CCSM   Init Mpi In msg_pass('connect') & msg_pass('disconnect'
				DBGINF
        call wav_mpi_initialized(Inited)
				DBGINF
        if(mode<0)then
            ! ============  Shut down Parallel environment
            if(Inited)then
              call wav_mpi_finalize('Deinit')
          endif
        else
            ! ===========  Initialise Parallel environment
				DBGINF
            if(.not. Inited)then
              call wav_mpi_init()
          endif
				DBGINF
            if(mpi_comm_wav==-1)mpi_comm_wav=MPI_COMM_WORLD
				DBGINF
            call wav_mpi_commrank(mpi_id)
				DBGINF
            call wav_mpi_commsize(mpi_npe)
				DBGINF
            call MPI_Get_processor_name(NodeName,nl,ierr);
				DBGINF
            if(mpi_id<10)write(6,*)"MPI",mpi_id,"NodeName=",trim(NodeName)
				DBGINF
            write(mpi_ids,'(".",i4.4)')mpi_id
            call flush(6)
        end if
    end subroutine wav_init_mpi
SUBROUTINE wav_mpi_chkerr(rcode,sname,string)
   ! ----- arguments ---
   integer, intent(in) :: rcode  ! input MPI error code
   character(*),         intent(in) :: sname  ! message
   character(*),optional,intent(in) :: string ! message
   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_chkerr) '
   integer             :: len

! -------------------------------------------------------------------------------
! PURPOSE: layer on MPI error checking
! -------------------------------------------------------------------------------

  if (rcode /= MPI_SUCCESS) then
      if (present(string)) then
        write(6,*) trim(subName),":",trim(sname),rcode,trim(string)
    else
        write(6,*) trim(subName),":",trim(sname),rcode
    endif
      call wav_mpi_abort(subName)
    endif

END SUBROUTINE wav_mpi_chkerr

! ===============================================================================

SUBROUTINE wav_mpi_sendi0(lvec,pid,tag,string)
   ! ----- arguments ---
   integer, intent(in) :: lvec     ! send value
   integer, intent(in) :: pid      ! pid to send to
   integer, intent(in) :: tag      ! tag
   character(*),optional,intent(in) :: string   ! message

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_sendi0) '
   integer             :: lsize
   integer             :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: Send a vector of reals
! -------------------------------------------------------------------------------

   lsize = 1

   call MPI_SEND(lvec,lsize,MPI_INTEGER,pid,tag,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_sendi0

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_sendi1(lvec,pid,tag,string)
   ! ----- arguments ---
   integer, intent(in) :: lvec(:)  ! in/out local values
   integer, intent(in) :: pid      ! pid to send to
   integer, intent(in) :: tag      ! tag
   character(*),optional,intent(in) :: string   ! message

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_sendi1) '
   integer             :: lsize
   integer             :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: Send a vector of reals
! -------------------------------------------------------------------------------

   lsize = size(lvec)

   call MPI_SEND(lvec,lsize,MPI_INTEGER,pid,tag,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)


END SUBROUTINE wav_mpi_sendi1

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_sendr0(lvec,pid,tag,string)
   ! ----- arguments ---
   real(8),    intent(in) :: lvec     ! in/out local values
   integer, intent(in) :: pid      ! pid to send to
   integer, intent(in) :: tag      ! tag
   character(*),optional,intent(in) :: string   ! message

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_sendr0) '
   integer             :: lsize
   integer             :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: Send a vector of reals
! -------------------------------------------------------------------------------

   lsize = 1

   call MPI_SEND(lvec,lsize,MPI_REAL8,pid,tag,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_sendr0

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_sendr1(lvec,pid,tag,string)
   ! ----- arguments ---
   real(8),    intent(in) :: lvec(:)  ! in/out local values
   integer, intent(in) :: pid      ! pid to send to
   integer, intent(in) :: tag      ! tag
   character(*),optional,intent(in) :: string   ! message

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_sendr1) '
   integer             :: lsize
   integer             :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: Send a vector of reals
! -------------------------------------------------------------------------------

   lsize = size(lvec)

   call MPI_SEND(lvec,lsize,MPI_REAL8,pid,tag,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_sendr1

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_recvi0(lvec,pid,tag,string,srcid)
   ! ----- arguments ---
   integer, intent(out):: lvec     ! in/out local values
   integer, intent(in) :: pid      ! pid to recv from
   integer, intent(in) :: tag      ! tag
   character(*),optional,intent(in) :: string   ! message
   integer,optional::srcid
   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_recvi0) '
   integer             :: lsize
   integer             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer             :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
! -------------------------------------------------------------------------------

   lsize = 1

   call MPI_RECV(lvec,lsize,MPI_INTEGER,pid,tag,mpi_comm_wav,status,ierr)
   if (present(srcid)) srcid=status(MPI_SOURCE)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_recvi0

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_recvi1(lvec,pid,tag,string,srcid)
   ! ----- arguments ---
   integer, intent(out):: lvec(:)  ! in/out local values
   integer, intent(in) :: pid      ! pid to recv from
   integer, intent(in) :: tag      ! tag
   character(*),optional,intent(in) :: string   ! message
   integer,optional::srcid

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_recvi1) '
   integer             :: lsize
   integer             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer             :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
! -------------------------------------------------------------------------------

   lsize = size(lvec)

   call MPI_RECV(lvec,lsize,MPI_INTEGER,pid,tag,mpi_comm_wav,status,ierr)
   if (present(srcid)) srcid=status(MPI_SOURCE)
   call wav_mpi_chkerr(ierr,subName,string)
    lvec=0

END SUBROUTINE wav_mpi_recvi1

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_recvr0(lvec,pid,tag,string,srcid)
   ! ----- arguments ---
   real(8),    intent(out):: lvec     ! in/out local values
   integer, intent(in) :: pid      ! pid to recv from
   integer, intent(in) :: tag      ! tag
   character(*),optional,intent(in) :: string   ! message
   integer,optional::srcid

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_recvr0) '
   integer             :: lsize
   integer             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer             :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
! -------------------------------------------------------------------------------

   lsize = 1
   call MPI_RECV(lvec,lsize,MPI_REAL8,pid,tag,mpi_comm_wav,status,ierr)
   if (present(srcid)) srcid=status(MPI_SOURCE)
   call wav_mpi_chkerr(ierr,subName,string)
    lvec=0

END SUBROUTINE wav_mpi_recvr0

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_recvr1(lvec,pid,tag,string,srcid)

   ! ----- arguments ---
   real(8),    intent(out):: lvec(:)  ! in/out local values
   integer, intent(in) :: pid      ! pid to recv from
   integer, intent(in) :: tag      ! tag
   character(*),optional,intent(in) :: string   ! message
   integer,optional::srcid

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_recvr1) '
   integer             :: lsize
   integer             :: status(MPI_STATUS_SIZE)  ! mpi status info
   integer             :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: Recv a vector of reals
! -------------------------------------------------------------------------------

   lsize = size(lvec)

   call MPI_RECV(lvec,lsize,MPI_REAL8,pid,tag,mpi_comm_wav,status,ierr)
   if (present(srcid)) srcid=status(MPI_SOURCE)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_recvr1
! ===============================================================================
SUBROUTINE wav_mpi_bcasts0(vec,root_,string)
   ! ----- arguments ---
   character*(*), intent(inout):: vec      ! vector of 1
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_bcasts) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_
   lsize = len(vec)
   call MPI_BCAST(vec,lsize,MPI_CHARACTER,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_bcasts0
! ===============================================================================
SUBROUTINE wav_mpi_bcasts1(vec,root_,string)
   ! ----- arguments ---
   character*(*), intent(inout):: vec(:)      ! vector of 1
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message
   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_bcasts) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_
! -------------------------------------------------------------------------------
! PURPOSE: Broadcast an integer
! -------------------------------------------------------------------------------
   lsize = len(vec(1))*size(vec)
   call MPI_BCAST(vec,lsize,MPI_CHARACTER,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_bcasts1

! ===============================================================================
SUBROUTINE wav_mpi_bcastl0(vec,root_,string)
   ! ----- arguments ---
   logical, intent(inout):: vec      ! vector of 1
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_bcastl0) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Broadcast an integer
! -------------------------------------------------------------------------------

   lsize = 1

   call MPI_BCAST(vec,lsize,MPI_LOGICAL,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_bcastl0

! ===============================================================================

SUBROUTINE wav_mpi_bcastw0(vec,root_,string)
   integer(2), intent(inout):: vec      ! vector of 1
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message
   character(*),parameter             :: subName = '(wav_mpi_bcasti0) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_
   lsize = 1
   call MPI_BCAST(vec,lsize,MPI_INTEGER2,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_bcastw0

! ===============================================================================
SUBROUTINE wav_mpi_bcastw1(vec,root_,string)
   integer(2), intent(inout):: vec(:)   ! vector
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message
   character(*),parameter             :: subName = '(wav_mpi_bcasti1) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_
   lsize = size(vec)
   call MPI_BCAST(vec,lsize,MPI_INTEGER2,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_bcastw1
SUBROUTINE wav_mpi_bcastw2(vec,root_,string)
   integer(2), intent(inout):: vec(:,:)   ! vector
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message
   character(*),parameter             :: subName = '(wav_mpi_bcasti1) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_
   lsize = size(vec)
   call MPI_BCAST(vec,lsize,MPI_INTEGER2,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_bcastw2

SUBROUTINE wav_mpi_bcasti0(vec,root_,string)
   ! ----- arguments ---
   integer, intent(inout):: vec      ! vector of 1
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_bcasti0) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Broadcast an integer
! -------------------------------------------------------------------------------

   lsize = 1

   call MPI_BCAST(vec,lsize,MPI_INTEGER,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_bcasti0

! ===============================================================================
SUBROUTINE wav_mpi_bcasti1(vec,root_,string)
   ! ----- arguments ---
   integer, intent(inout):: vec(:)   ! vector
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_bcasti1) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_
! -------------------------------------------------------------------------------
! PURPOSE: Broadcast a vector of integers
! -------------------------------------------------------------------------------

   lsize = size(vec)

   call MPI_BCAST(vec,lsize,MPI_INTEGER,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_bcasti1
SUBROUTINE wav_mpi_bcasti2(vec,root_,string)
   ! ----- arguments ---
   integer, intent(inout):: vec(:,:)   ! vector
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_bcasti1) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Broadcast a vector of integers
! -------------------------------------------------------------------------------

   lsize = size(vec)
   call MPI_BCAST(vec,lsize,MPI_INTEGER,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_bcasti2
! ===============================================================================
SUBROUTINE wav_mpi_bcasti80(vec,root_,string)
   integer(8), intent(inout):: vec      ! vector of 1
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message
   character(*),parameter             :: subName = '(wav_mpi_bcasti0) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_
   lsize = 1
   call MPI_BCAST(vec,lsize,MPI_INTEGER8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_bcasti80
SUBROUTINE wav_mpi_bcasti81(vec,root_,string)
   integer(8), intent(inout):: vec(:)   ! vector
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message
   character(*),parameter             :: subName = '(wav_mpi_bcasti1) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_
   lsize = size(vec)
   call MPI_BCAST(vec,lsize,MPI_INTEGER8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_bcasti81
! ===============================================================================

SUBROUTINE wav_mpi_bcastr0(vec,root_,string)
   ! ----- arguments ---
   real(8),    intent(inout):: vec      ! vector of 1
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_bcastr0) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Broadcast a real
! -------------------------------------------------------------------------------
   lsize = 1

   call MPI_BCAST(vec,lsize,MPI_REAL8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_bcastr0

! ===============================================================================
! ===============================================================================


! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_bcastr1(vec,root_,string)
   ! ----- arguments ---
   real(8),    intent(inout):: vec(:)   ! vector
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_bcastr1) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Broadcast a vector of reals
! -------------------------------------------------------------------------------

   lsize = size(vec)

   call MPI_BCAST(vec,lsize,MPI_REAL8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_bcastr1
SUBROUTINE wav_mpi_bcastb2(vec,root_,string)
   ! ----- arguments ---
   integer*1, intent(inout):: vec(:,:)   ! vector
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_bcasti1) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Broadcast a vector of integers
! -------------------------------------------------------------------------------

   lsize = size(vec)

   call MPI_BCAST(vec,lsize,MPI_BYTE,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_bcastb2
SUBROUTINE wav_mpi_bcastr2(vec,root_,string)
   ! ----- arguments ---
   real(8),    intent(inout):: vec(:,:)   ! vector
   integer,optional::root_
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_bcastr2) '
   integer  :: ierr,lsize,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Broadcast a vector of reals
! -------------------------------------------------------------------------------
   lsize = size(vec)
   call MPI_BCAST(vec,lsize,MPI_REAL8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_bcastr2
SUBROUTINE wav_mpi_bcastp(pk,root_,string)
   ! ----- arguments ---
   type(mpipacket),intent(inout):: pk   ! vector
   integer,optional::root_
   character(*),optional,intent(in):: string   ! message
   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_bcastr2) '
   integer               :: ierr
   integer               :: root
   root=0;if(present(root_))root=root_
! -------------------------------------------------------------------------------
! PURPOSE: Broadcast a vector of reals
! -------------------------------------------------------------------------------
   call MPI_BCAST(pk%dsize,1,MPI_INTEGER,root,mpi_comm_wav,ierr)
   call MPI_BCAST(pk%Buf,pk%dsize,MPI_PACKED,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_bcastp

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_sumi0(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   integer, intent(in) :: lvec     ! in/out local values
   integer, intent(out):: gvec     ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_sumi0) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_sumi0

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_sumi1(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   integer, intent(in) :: lvec(:)  ! in/out local values
   integer, intent(out):: gvec(:)  ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_sumi1) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_sumi1

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_sumr0(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   real(8),    intent(in) :: lvec     ! in/out local values
   real(8),    intent(out):: gvec     ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_sumr0) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_sumr0

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_sumr1(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   real(8),    intent(in) :: lvec(:)  ! in/out local values
   real(8),    intent(out):: gvec(:)  ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_sumr1) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_sumr1

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_sumr2(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   real(8),    intent(in) :: lvec(:,:)! in/out local values
   real(8),    intent(out):: gvec(:,:)! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_sumr2) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_sumr2

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_sumr3(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   real(8),    intent(in) :: lvec(:,:,:) ! in/out local values
   real(8),    intent(out):: gvec(:,:,:) ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_sumr3) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_
! -------------------------------------------------------------------------------
! PURPOSE: Finds sum of a distributed vector of values, assume local sum
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_SUM
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_sumr3

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_mini0(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   integer, intent(in) :: lvec     ! in/out local values
   integer, intent(out):: gvec     ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_mini0) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_mini0

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_mini1(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   integer, intent(in) :: lvec(:)  ! in/out local values
   integer, intent(out):: gvec(:)  ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_mini1) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_mini1

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_minr0(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   real(8),    intent(in) :: lvec     ! in/out local values
   real(8),    intent(out):: gvec     ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_minr0) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_minr0

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_minr1(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   real(8),    intent(in) :: lvec(:)  ! in/out local values
   real(8),    intent(out):: gvec(:)  ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_minr1) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Finds min of a distributed vector of values, assume local min
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_MIN
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_minr1

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_maxi0(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   integer, intent(in) :: lvec     ! in/out local values
   integer, intent(out):: gvec     ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_maxi0) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif
END SUBROUTINE wav_mpi_maxi0

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_maxi1(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   integer, intent(in) :: lvec(:)  ! in/out local values
   integer, intent(out):: gvec(:)  ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_maxi1) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_
! -------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_INTEGER,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_maxi1

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_maxr0(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   real(8),    intent(in) :: lvec     ! in/out local values
   real(8),    intent(out):: gvec     ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_maxr0) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = 1
   gsize = 1

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_maxr0

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_maxr1(lvec,gvec,root_,string,all)
   ! ----- arguments ---
   real(8),    intent(in) :: lvec(:)  ! in/out local values
   real(8),    intent(out):: gvec(:)  ! in/out global values
   integer,optional::root_
   character(*),optional,intent(in) :: string   ! message
   logical,     optional,intent(in) :: all      ! allreduce if true

   ! ----- local ---
   character(*),parameter           :: subName = '(wav_mpi_maxr1) '
   logical                          :: lall
   integer             :: reduce_type  ! mpi reduction type
   integer :: lsize,gsize,ierr,root
   root=0;if(present(root_))root=root_

! -------------------------------------------------------------------------------
! PURPOSE: Finds max of a distributed vector of values, assume local max
!          already computed
! -------------------------------------------------------------------------------

   reduce_type = MPI_MAX
   if (present(all)) then
     lall = all
   else
     lall = .false.
   endif

   lsize = size(lvec)
   gsize = size(gvec)

   if (lsize /= gsize) then
     call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string))
   endif

   if (lall) then
     call MPI_ALLREDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string)
   else
     call MPI_REDUCE(lvec,gvec,gsize,MPI_REAL8,reduce_type,root,mpi_comm_wav,ierr)
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string)
   endif

END SUBROUTINE wav_mpi_maxr1

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_commsize(size,string)
   ! ----- arguments ---
   integer,intent(out)                :: size
   character(*),optional,intent(in)   :: string   ! message
   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_commsize) '
   integer               :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: MPI commsize
! -------------------------------------------------------------------------------

   call MPI_COMM_SIZE(mpi_comm_wav,size,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_commsize

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_commrank(rank,string)
   ! ----- arguments ---
   integer,intent(out)                :: rank
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_commrank) '
   integer               :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: MPI commrank
! -------------------------------------------------------------------------------

   call MPI_COMM_RANK(mpi_comm_wav,rank,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_commrank

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_initialized(flag,string)
   ! ----- arguments ---
   logical,intent(out)                :: flag
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_initialized) '
   integer               :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: MPI initialized
! -------------------------------------------------------------------------------

				DBGINF
   call MPI_INITIALIZED(flag,ierr)
				DBGINF
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_initialized

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_abort(string,rcode)
   ! ----- arguments ---
   character(*),optional,intent(in)   :: string   ! message
   integer,optional,intent(in)        :: rcode    ! optional code

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_abort) '
   integer               :: ierr
   integer rcode_
! -------------------------------------------------------------------------------
! PURPOSE: MPI abort
! -------------------------------------------------------------------------------
    rcode_=-1
   if(present(string))then
     if(present(rcode))then
       write(6,*) subName,":",trim(string),rcode
       rcode_=rcode
     else
       write(6,*) subName,":",trim(string)
     endif
   else
     if(present(rcode))then
       write(6,*) subName," Err Code:",rcode
       rcode_=rcode
     else
       write(6,*) subName
     endif
   endif
   call MPI_ABORT(MPI_COMM_WORLD,rcode_,ierr)
   stop

END SUBROUTINE wav_mpi_abort

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_barrier(string)
   ! ----- arguments ---
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_barrier) '
   integer               :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: MPI wav_mpi_barrier
! -------------------------------------------------------------------------------

   call MPI_BARRIER(mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_barrier

! ===============================================================================
! ===============================================================================

SUBROUTINE wav_mpi_init(string)
   ! ----- arguments ---
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_init) '
   integer               :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: MPI init
! -------------------------------------------------------------------------------

   call MPI_INIT(ierr)
   call wav_mpi_chkerr(ierr,subName,string)
    if(mpi_comm_wav==-1)mpi_comm_wav=MPI_COMM_WORLD
    call wav_mpi_commrank(mpi_id)
    call wav_mpi_commsize(mpi_npe)
    NextProc=mpi_id+1;if(NextProc>=mpi_npe)NextProc=mpi_proc_null
    PrevProc=mpi_id-1;if(PrevProc< 0      )PrevProc=mpi_proc_null
END SUBROUTINE wav_mpi_init
! ===============================================================================

SUBROUTINE wav_mpi_finalize(string)
   ! ----- arguments ---
   character(*),optional,intent(in)   :: string   ! message

   ! ----- local ---
   character(*),parameter             :: subName = '(wav_mpi_finalize) '
   integer               :: ierr

! -------------------------------------------------------------------------------
! PURPOSE: MPI finalize
! -------------------------------------------------------------------------------

   call MPI_FINALIZE(ierr)
   call wav_mpi_chkerr(ierr,subName,string)

END SUBROUTINE wav_mpi_finalize

subroutine wav_mpi_gatheri0(sbuf,rbuf,root_,string)
    integer sbuf,rbuf(:)
    integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_gatheri0) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
   call mpi_Gather(sbuf,1,mpi_integer,rbuf,1,mpi_integer,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_gatheri0


subroutine wav_mpi_gatheri1(sbuf,rbuf,root_,string)
    integer sbuf(:),rbuf(:)
    integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    integer ::count
    character(*),parameter             :: subName = '(wav_mpi_gatheri1) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
    count=size(sbuf)
    call mpi_Gather(sbuf,count,mpi_integer,rbuf,count,mpi_integer,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_gatheri1
subroutine wav_mpi_gatherr0(sbuf,rbuf,root_,string)
    real*4 sbuf,rbuf(:)
    integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_gatherr0) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
   call mpi_Gather(sbuf,1,MPI_REAL4,rbuf,1,MPI_REAL4,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_gatherr0


subroutine wav_mpi_gatherr1(sbuf,rbuf,root_,string)
    real*4 sbuf(:),rbuf(:)
      integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    integer ::count
    character(*),parameter             :: subName = '(wav_mpi_gatherr1) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
    count=size(sbuf)
    call mpi_Gather(sbuf,count,MPI_REAL4,rbuf,count,MPI_REAL4,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_gatherr1
subroutine wav_mpi_gatherd0(sbuf,rbuf,root_,string)
    real*8 sbuf,rbuf(:)
    integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_gatherr0) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
   call mpi_Gather(sbuf,1,MPI_REAL8,rbuf,1,MPI_REAL8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_gatherd0


subroutine wav_mpi_gatherd1(sbuf,rbuf,root_,string)
    real*8 sbuf(:),rbuf(:)
    integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    integer ::count
    character(*),parameter             :: subName = '(wav_mpi_gatherr1) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
    count=size(sbuf)
    call mpi_Gather(sbuf,count,MPI_REAL8,rbuf,count,MPI_REAL8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_gatherd1
subroutine wav_mpi_scatteri0(sbuf,rbuf,root_,string)
    integer sbuf(:),rbuf
    integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_scatteri0) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
   call mpi_scatter(sbuf,1,mpi_integer,rbuf,1,mpi_integer,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scatteri0
subroutine wav_mpi_scatteri1(sbuf,rbuf,root_,string)
    integer sbuf(:),rbuf(:)
    integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    integer ::count
    character(*),parameter             :: subName = '(wav_mpi_scatteri1) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
    count=size(rbuf)
    call mpi_scatter(sbuf,count,mpi_integer,rbuf,count,mpi_integer,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scatteri1
subroutine wav_mpi_scatteri80(sbuf,rbuf,root_,string)
    integer*8 sbuf(:),rbuf
    integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_scatteri0) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
   call mpi_scatter(sbuf,1,mpi_integer8,rbuf,1,mpi_integer8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scatteri80
subroutine wav_mpi_scatteri81(sbuf,rbuf,root_,string)
    integer*8 sbuf(:),rbuf(:)
    integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    integer ::count
    character(*),parameter             :: subName = '(wav_mpi_scatteri1) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
    count=size(rbuf)
    call mpi_scatter(sbuf,count,mpi_integer8,rbuf,count,mpi_integer8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scatteri81
subroutine wav_mpi_scatterr0(sbuf,rbuf,root_,string)
    real*4 sbuf(0),rbuf
    integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_scatterr0) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
   call mpi_scatter(sbuf,1,MPI_REAL4,rbuf,1,MPI_REAL4,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scatterr0


subroutine wav_mpi_scatterr1(sbuf,rbuf,root_,string)
    real*4 sbuf(:),rbuf(:)
      integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    integer ::count
    character(*),parameter             :: subName = '(wav_mpi_scatterr1) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
    count=size(rbuf)
    call mpi_scatter(sbuf,count,MPI_REAL4,rbuf,count,MPI_REAL4,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scatterr1
subroutine wav_mpi_scatterd0(sbuf,rbuf,root_,string)
    real*8 sbuf(:),rbuf
    integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_scatterr0) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
   call mpi_scatter(sbuf,1,MPI_REAL8,rbuf,1,MPI_REAL8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scatterd0


subroutine wav_mpi_scatterd1(sbuf,rbuf,root_,string)
    real*8 sbuf(:),rbuf(:)
    integer,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    integer ::count
    character(*),parameter             :: subName = '(wav_mpi_scatterr1) '
    integer :: ierr, root
   root=0;if(present(root_))root=root_
    count=size(rbuf)
    call mpi_scatter(sbuf,count,MPI_REAL8,rbuf,count,MPI_REAL8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scatterd1

subroutine wav_mpi_gathervi1(sbuf,scount,rbuf,rcounts,displs,root_,string)
    integer sbuf(:),rbuf(:)
    integer ::scount,rcounts(:),displs(:)
    integer ,optional::root_
  character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_gatheri1) '
    integer :: ierr,root
    root=0;if(present(root_))root=root_
   call MPI_GATHERV(sbuf,scount,mpi_integer,rbuf,rcounts,displs,mpi_integer,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_gathervi1
subroutine wav_mpi_gathervi2(sbuf,scount,rbuf,rcounts,displs,root_,string)
    integer sbuf(:,:),rbuf(:,:)
    integer ::scount,rcounts(:),displs(:)
    integer ,optional::root_
    character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_gatheri2) '
    integer :: ierr,root
    root=0;if(present(root_))root=root_
   call MPI_GATHERV(sbuf,scount,mpi_integer,rbuf,rcounts,displs,mpi_integer,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_gathervi2
subroutine wav_mpi_gathervr1(sbuf,scount,rbuf,rcounts,displs,root_,string)
    real*4 sbuf(:),rbuf(:)
    integer ::scount,rcounts(:),displs(:)
    integer ,optional::root_
    character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_gatherr1) '
    integer :: ierr,root
    root=0;if(present(root_))root=root_
   call MPI_GATHERV(sbuf,scount,MPI_REAL4,rbuf,rcounts,displs,MPI_REAL4,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_gathervr1
subroutine wav_mpi_gathervr2(sbuf,scount,rbuf,rcounts,displs,root_,string)
    real*4 sbuf(:,:),rbuf(:,:)
    integer ::scount,rcounts(:),displs(:)
    integer ,optional::root_
    character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_gatherr2) '
    integer :: ierr,root
    root=0;if(present(root_))root=root_
   call MPI_GATHERV(sbuf,scount,MPI_REAL4,rbuf,rcounts,displs,MPI_REAL4,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_gathervr2
subroutine wav_mpi_gathervd1(sbuf,scount,rbuf,rcounts,displs,root_,string)
    real*8 sbuf(:),rbuf(:)
    integer ::scount,rcounts(:),displs(:)
    integer ,optional::root_
    character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_gatherr1) '
    integer :: ierr,root
    root=0;if(present(root_))root=root_
   call MPI_GATHERV(sbuf,scount,MPI_REAL8,rbuf,rcounts,displs,MPI_REAL8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_gathervd1

subroutine wav_mpi_scattervi1(sbuf,scounts,displs,rbuf,rcount,root_,string)
    integer sbuf(:),rbuf(:)
    integer ::scounts(:),displs(:),rcount
    integer ,optional::root_
    character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_scatteri1) '
    integer :: ierr,root
    root=0;if(present(root_))root=root_
    call MPI_SCATTERV(sbuf,scounts,displs,mpi_integer,rbuf,rcount,mpi_integer,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scattervi1
subroutine wav_mpi_scattervi2(sbuf,scounts,displs,rbuf,rcount,root_,string)
    integer sbuf(:,:),rbuf(:,:)
    integer ::scounts(:),displs(:),rcount
    integer ,optional::root_
    character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_scatteri2) '
    integer :: ierr,root
    root=0;if(present(root_))root=root_
    call MPI_SCATTERV(sbuf,scounts,displs,mpi_integer,rbuf,rcount,mpi_integer,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scattervi2
subroutine wav_mpi_scattervr1(sbuf,scounts,displs,rbuf,rcount,root_,string)
    real*4 sbuf(:),rbuf(:)
    integer ::scounts(:),displs(:),rcount
    integer ,optional::root_
    character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_scatterr1) '
    integer :: ierr,root
    root=0;if(present(root_))root=root_
    ! DBGINFA,'wav_mpi_scattervr1',shape(sbuf),shape(rbuf),rcount,scounts,sbuf(rcount)
    call MPI_SCATTERV(sbuf,scounts,displs,MPI_REAL4,rbuf,rcount,MPI_REAL4,root,mpi_comm_wav,ierr)
    ! SDBGINF,'wav_mpi_scattervr1 End'
    ! DBGINFA,rbuf(rcount)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scattervr1
subroutine wav_mpi_scattervr2(sbuf,scounts,displs,rbuf,rcount,root_,string)
    real*4 sbuf(:,:),rbuf(:,:)
    integer ::scounts(:),displs(:),rcount
    integer ,optional::root_
    character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_scatterr2) '
    integer :: ierr,root
    root=0;if(present(root_))root=root_
    call MPI_SCATTERV(sbuf,scounts,displs,MPI_REAL4,rbuf,rcount,MPI_REAL4,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scattervr2
subroutine wav_mpi_scattervd1(sbuf,scounts,displs,rbuf,rcount,root_,string)
    real*8 sbuf(:),rbuf(:)
    integer ::scounts(:),displs(:),rcount
    integer ,optional::root_
    character(*),optional,intent(in)   :: string   ! message
    character(*),parameter             :: subName = '(wav_mpi_scatterr1) '
    integer :: ierr,root
    root=0;if(present(root_))root=root_
    call MPI_SCATTERV(sbuf,scounts,displs,MPI_REAL8,rbuf,rcount,MPI_REAL8,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
end subroutine wav_mpi_scattervd1

! ===============================================================================
SUBROUTINE InitMpiPacket(pk,lsize)
  type( mpipacket) pk
  integer*4,optional::lsize
  if(present(lsize))then
    if(allocated(pk%buf))then
      deallocate(pk%buf)
    endif
    if(lsize>0)then
      pk%bsize=lsize
      allocate(pk%buf((pk%bsize+4)/4) )
    endif
  endif
  pk%dsize=0
  pk%pos=0
end SUBROUTINE InitMpiPacket
SUBROUTINE SetMpiPacketDSize(pk,dsize)
   type( mpipacket) pk
   integer*4,optional::dsize
   if(present(dsize))then
    pk%dsize=dsize
   else
    pk%dsize=pk%pos
   endif
end SUBROUTINE SetMpiPacketDSize
SUBROUTINE wav_mpi_packi0(pk,iobuf,iunpack)
   type( mpipacket) pk
   integer(4),intent(inout) :: iobuf
   integer,optional,intent(in)::iunpack
   integer(4):: iblen,ierr
   iblen=1
   if(present(iunpack).and.iunpack/=0)then
    call MPI_UNPACK(pk%buf,pk%dsize,pk%pos,iobuf,iblen,MPI_INTEGER,mpi_comm_wav,ierr)
   else
     call MPI_PACK(iobuf,iblen,MPI_INTEGER,pk%buf,pk%bsize,pk%pos,mpi_comm_wav,ierr)
   endif
END SUBROUTINE wav_mpi_packi0
SUBROUTINE wav_mpi_packi1(pk,iobuf,iunpack)
   type( mpipacket) pk
   integer(4), intent(inout) :: iobuf(:)
   integer,optional,intent(in)::iunpack
   integer(4) :: iblen,ierr
   iblen=size(iobuf)
   if(present(iunpack).and.iunpack/=0)then
     call MPI_UNPACK(pk%buf,pk%dsize,pk%pos,iobuf,iblen,MPI_INTEGER,mpi_comm_wav,ierr)
   else
     call MPI_PACK(iobuf,iblen,MPI_INTEGER,pk%buf,pk%bsize,pk%pos,mpi_comm_wav,ierr)
   endif
END SUBROUTINE wav_mpi_packi1
SUBROUTINE wav_mpi_packr0(pk,iobuf,iunpack)
   type( mpipacket) pk
   real(4), intent(inout) :: iobuf
   integer,optional,intent(in)::iunpack
   integer(4) :: iblen,ierr
   iblen=1
   if(present(iunpack).and.iunpack/=0)then
     call MPI_UNPACK(pk%buf,pk%dsize,pk%pos,iobuf,iblen,MPI_REAL4,mpi_comm_wav,ierr)
   else
     call MPI_PACK(iobuf,iblen,MPI_REAL4,pk%buf,pk%bsize,pk%pos,mpi_comm_wav,ierr)
   endif
END SUBROUTINE wav_mpi_packr0
SUBROUTINE wav_mpi_packr1(pk,iobuf,iunpack)
   type( mpipacket) pk
   real(4), intent(inout) :: iobuf(:)
   integer,optional,intent(in)::iunpack
   integer(4) :: iblen,ierr
   iblen=size(iobuf)
   if(present(iunpack).and.iunpack/=0)then
     call MPI_UNPACK(pk%buf,pk%dsize,pk%pos,iobuf,iblen,MPI_REAL4,mpi_comm_wav,ierr)
   else
     call MPI_PACK(iobuf,iblen,MPI_REAL4,pk%buf,pk%bsize,pk%pos,mpi_comm_wav,ierr)
   endif
END SUBROUTINE wav_mpi_packr1
SUBROUTINE wav_mpi_packd0(pk,iobuf,iunpack)
   type( mpipacket) pk
   real(8), intent(inout) :: iobuf
   integer,optional,intent(in)::iunpack
   integer(4) :: iblen,ierr
   iblen=1
   if(present(iunpack).and.iunpack/=0)then
     call MPI_UNPACK(pk%buf,pk%dsize,pk%pos,iobuf,iblen,MPI_REAL8,mpi_comm_wav,ierr)
   else
     call MPI_PACK(iobuf,iblen,MPI_REAL8,pk%buf,pk%bsize,pk%pos,mpi_comm_wav,ierr)
   endif
END SUBROUTINE wav_mpi_packd0
SUBROUTINE wav_mpi_packd1(pk,iobuf,iunpack)
   type( mpipacket) pk
   real(8), intent(inout) :: iobuf(:)
   integer,optional,intent(in)::iunpack
   integer(4) :: iblen,ierr
   iblen=size(iobuf)
   if(present(iunpack).and.iunpack/=0)then
     call MPI_UNPACK(pk%buf,pk%dsize,pk%pos,iobuf,iblen,MPI_REAL8,mpi_comm_wav,ierr)
   else
     call MPI_PACK(iobuf,iblen,MPI_REAL8,pk%buf,pk%bsize,pk%pos,mpi_comm_wav,ierr)
   endif
END SUBROUTINE wav_mpi_packd1
SUBROUTINE wav_mpi_packc0(pk,iobuf,iunpack)
   type( mpipacket) pk
   character*(*), intent(inout) :: iobuf
   integer,optional,intent(in)::iunpack
   integer(4) :: iblen,ierr
   if(present(iunpack).and.iunpack/=0)then
      iobuf=''
     call MPI_UNPACK(pk%buf,pk%bsize,pk%pos,iblen,1,MPI_INTEGER,mpi_comm_wav,ierr)
     call MPI_UNPACK(pk%buf,pk%bsize,pk%pos,iobuf,iblen,MPI_CHARACTER,mpi_comm_wav,ierr)
   else
     iblen=len_trim(iobuf)
     call MPI_PACK(iblen,1,MPI_INTEGER,pk%buf,pk%bsize,pk%pos,mpi_comm_wav,ierr)
     call MPI_PACK(iobuf,iblen,MPI_CHARACTER,pk%buf,pk%bsize,pk%pos,mpi_comm_wav,ierr)
   endif
END SUBROUTINE wav_mpi_packc0
SUBROUTINE wav_mpi_packc1(pk,iobuf,iunpack)
   type( mpipacket) pk
   character*(*), intent(inout) :: iobuf(:)
   integer,optional,intent(in)::iunpack
   integer(4) :: iblen,i,nn,ierr
   nn=size(iobuf)
   if(present(iunpack).and.iunpack/=0)then
     do i=1,nn
      iobuf(i)=''
       call MPI_UNPACK(pk%buf,pk%bsize,pk%pos,iblen,1,MPI_INTEGER,mpi_comm_wav,ierr)
       call MPI_UNPACK(pk%buf,pk%bsize,pk%pos,iobuf(i),iblen,MPI_CHARACTER,mpi_comm_wav,ierr)
     enddo
   else
     do i=1,nn
       iblen=len_trim(iobuf(i))
       call MPI_PACK(iblen,1,MPI_INTEGER,pk%buf,pk%bsize,pk%pos,mpi_comm_wav,ierr)
       call MPI_PACK(iobuf(i),iblen,MPI_CHARACTER,pk%buf,pk%bsize,pk%pos,mpi_comm_wav,ierr)
     enddo
   endif
END SUBROUTINE wav_mpi_packc1

! ==========================================================
SUBROUTINE wav_mpi_unpacki0(pk,obuf)
   type( mpipacket) pk
   integer(4), intent(out)::obuf
   call wav_mpi_pack(pk,obuf,1)
END SUBROUTINE wav_mpi_unpacki0
SUBROUTINE wav_mpi_unpacki1(pk,obuf)
   type( mpipacket) pk
   integer(4), intent(out)::obuf(:)
   call wav_mpi_pack(pk,obuf,1)
END SUBROUTINE wav_mpi_unpacki1
SUBROUTINE wav_mpi_unpackr0(pk,obuf)
   type( mpipacket) pk
   real(4), intent(out)::obuf
   call wav_mpi_pack(pk,obuf,1)
END SUBROUTINE wav_mpi_unpackr0
SUBROUTINE wav_mpi_unpackr1(pk,obuf)
   type( mpipacket) pk
   real(4), intent(out)::obuf(:)
   call wav_mpi_pack(pk,obuf,1)
END SUBROUTINE wav_mpi_unpackr1
SUBROUTINE wav_mpi_unpackd0(pk,obuf)
   type( mpipacket) pk
   real(8), intent(out)::obuf
   call wav_mpi_pack(pk,obuf,1)
END SUBROUTINE wav_mpi_unpackd0
SUBROUTINE wav_mpi_unpackd1(pk,obuf)
   type( mpipacket) pk
   real(8), intent(out)::obuf(:)
   call wav_mpi_pack(pk,obuf,1)
END SUBROUTINE wav_mpi_unpackd1
SUBROUTINE wav_mpi_unpackc0(pk,obuf)
   type( mpipacket) pk
   character*(*), intent(out) :: obuf
   call wav_mpi_pack(pk,obuf,1)
END SUBROUTINE wav_mpi_unpackc0
SUBROUTINE wav_mpi_unpackc1(pk,obuf)
   type( mpipacket) pk
   character*(*), intent(out) :: obuf(:)
   call wav_mpi_pack(pk,obuf,1)
END SUBROUTINE wav_mpi_unpackc1

subroutine wav_mpi_check
  integer iv
  if (mpi_id>0)then
    call wav_mpi_recv(iv,PrevProc,1010)
  else
    write(6,*)'wav_mpi_check begin'; call flush(6)
  endif
  write(6,*)'wav_mpi_check ',mpi_id,NextProc,trim(NodeName);  call flush(6)
  if(mpi_id<mpi_npe-1)then
    call wav_mpi_send(iv,NextProc,1010)
  else
    call wav_mpi_send(iv,0,1010)
  endif
  if(mpi_id==0)then
    call wav_mpi_recv(iv,mpi_npe-1,1010)
    write(6,*)'wav_mpi_check End'; call flush(6)
  endif
end SUBROUTINE wav_mpi_check

subroutine flush6_
  !call flush (6)
  !call flush (0)
end subroutine flush6_
SUBROUTINE SerRunf
   integer:: is,ir,ierr,status(MPI_STATUS_SIZE)  ! mpi status info
  is=0;
  call mpi_Sendrecv(is,1,MPI_INTEGER,PrevProc,1011,ir,1,MPI_INTEGER,NextProc,1011,mpi_comm_wav,status,ierr)
  call mpi_Sendrecv(is,1,MPI_INTEGER,NextProc,1010,ir,1,MPI_INTEGER,PrevProc,1010,mpi_comm_wav,status,ierr)
End SUBROUTINE SerRunf

SUBROUTINE SerRun
   integer:: is,ir,ierr,status(MPI_STATUS_SIZE)  ! mpi status info
  is=0;
  call mpi_Sendrecv(is,1,MPI_INTEGER,NextProc,1010,ir,1,MPI_INTEGER,PrevProc,1010,mpi_comm_wav,status,ierr)
End SUBROUTINE SerRun

SUBROUTINE SerRB
  integer ::iv
  call wav_mpi_recv(iv,PrevProc,1010)
End SUBROUTINE SerRB
SUBROUTINE SerRE
  integer ::iv=0
  call wav_mpi_send(iv,NextProc,1010)
End SUBROUTINE SerRE
SUBROUTINE SerRunr
   integer:: is=0,ir,ierr,status(MPI_STATUS_SIZE)  ! mpi status info
  call mpi_Sendrecv(is,1,MPI_INTEGER,PrevProc,1011,ir,1,MPI_INTEGER,NextProc,1011,mpi_comm_wav,status,ierr)
  call flush6_
End SUBROUTINE SerRunr
SUBROUTINE SerRBr
  integer ::iv=0
  call wav_mpi_recv(iv,NextProc,1011)
End SUBROUTINE SerRBr
SUBROUTINE SerREr
  integer ::iv
  call wav_mpi_send(iv,PrevProc,1011)
End SUBROUTINE SerREr
#endif
END MODULE wav_mpi_mod
