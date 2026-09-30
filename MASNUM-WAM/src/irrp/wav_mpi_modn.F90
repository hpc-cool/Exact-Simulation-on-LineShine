!#ifdef NO_MPI
!#include "wav_mpi_mod_dump.F90"
!#else
!#define MPI_GATHERV MYMPI_GATHERV
!#define MPI_SCATTERV MYMPI_SCATTERV
Module wav_mpi_mod
! -------------------------------------------------------------------------------
! PURPOSE: general layer on MPI functions
! -------------------------------------------------------------------------------
#ifndef NO_MPI
#ifdef GF13
  use mpi
#endif
#endif
   implicit none
   private
   integer::DBGLVL=10
   integer,public:: mpi_comm_wav=-1,mpi_id=0,mpi_npe=1
! mpi library include file
#ifndef NO_MPI
#ifndef GF13
  include "mpif.h"
#endif
   integer,public:: NextProc=mpi_proc_null,PrevProc=mpi_proc_null
   public::mpi_any_source,mpi_real8,mpi_real4,MPI_STATUS_SIZE,MPI_INTEGER,MPI_INTEGER2
!  PUBLIC: Public interfaces
   character*(8),public:: mpi_ids
   character*(MPI_MAX_PROCESSOR_NAME),public:: NodeName
   !character*(16),public:: NodeName
   !integer(4),public::nodeid
   !integer(4),allocatable,public:: NodeIds(:)
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
   public :: InitMpiPacket,wav_mpi_pack,SetMpiPacketDSize
   public :: SerRun,SerRB,SerRE,SerRunr,SerRBr,SerREr,SerRunf,wav_mpi_check,SerRRun,zbarr
   type,public:: mpipacket
    integer bsize,pos,dsize,id
    integer*4,allocatable::Buf(:)
   end type mpipacket
   type MPI_Common
    integer igrp
    integer id,npe,comm,pcomm,pcommid
    character*(32) cname;
    integer NextProc,PrevProc
   end type MPI_Common
#define WMIP(N) interface N; module procedure  \
    N##i, N##l, N##r, N##d ;\
   end interface
#define WMIPC(N) interface N; module procedure  \
    N##w,N##i, N##l, N##r, N##d ,N##c0,N##c1 ;\
   end interface
#define WMIPP(N) interface N; module procedure  \
    N##w,N##i, N##l, N##r, N##d ,N##p ,N##c;\
   end interface

   WMIP(wav_mpi_send )
   WMIP(wav_mpi_recv )
   WMIP(wav_mpi_gather )
   WMIP(wav_mpi_scatter )
   WMIP(wav_mpi_gatherv )
   WMIP(wav_mpi_scatterv )
   WMIP(wav_mpi_reduce)
   WMIP(wav_mpi_sum)
   WMIP(wav_mpi_min)
   WMIP(wav_mpi_max)
   WMIPC(wav_mpi_pack)
   WMIPP(wav_mpi_bcast)

CONTAINS
SUBROUTINE wav_mpi_split()
  type(MPI_Common) zcomm
  integer ierr
  if(zcomm%igrp>0)then
    call MPI_COMM_SPLIT(zcomm%pcomm,zcomm%igrp,0,zcomm%comm,ierr);
    call MPI_COMM_RANK (zcomm%comm ,zcomm%id  ,             ierr)
    call MPI_COMM_SIZE (zcomm%comm ,zcomm%npe ,             ierr)
  endif
end SUBROUTINE wav_mpi_split
! ===============================================================================
  subroutine wav_init_mpi(mode )
    integer mode
    integer ierr,nl
    logical Inited
! int COUPLE_CCSM   Init Mpi In msg_pass('connect') & msg_pass('disconnect'
    call wav_mpi_initialized(Inited)
    if(mode<0)then
      ! ============  Shut down Parallel environment
      if(Inited)then
        call wav_mpi_finalize('Deinit')
      endif
    else
      ! ===========  Initialise Parallel environment
      if(.not. Inited)then
        call wav_mpi_init()
      endif
      if(mpi_comm_wav==-1)mpi_comm_wav=MPI_COMM_WORLD
      call wav_mpi_commrank(mpi_id)
      call wav_mpi_commsize(mpi_npe)
      call MPI_Get_processor_name(NodeName,nl,ierr);
      !if(mpi_id<10)write(6,*)"MPI",mpi_id,"NodeName=",trim(NodeName)
      !read(NodeName(11:12),'(i)')nodeid
      !allocate(NodeIds(0:mpi_npe-1))
      !call mpi_Gather(nodeid,1,mpi_integer,NodeIds,1,mpi_integer,0,mpi_comm_wav,ierr)
      write(mpi_ids,'(".",i4.4)')mpi_id
      call flush(6)
    end if
  end subroutine wav_init_mpi
  SUBROUTINE wav_mpi_chkerr(rcode,sname,string)
  integer, intent(in) :: rcode  ! input MPI error code
  character(*),         intent(in) :: sname  ! message
  character(*),optional,intent(in) :: string ! message
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
! PURPOSE: Send a vector of reals
integer function strsize(lvec,ni)
  character(*) ::lvec(ni)
  integer ni
  strsize=ni* len(lvec(1))
end function
#define SEND(TAG,VT,MVT)  \
SUBROUTINE wav_mpi_send##TAG(lvec,pid,tag,string) ;\
   VT, intent(in) :: lvec(..) ;\
   integer, intent(in) :: pid      ;\
   integer, intent(in) :: tag      ;\
   character(*),optional,intent(in) :: string;\
   character(*),parameter           :: subName = "(wav_mpi_send" // #TAG // ")";\
   integer             :: lsize ,i,ierr;\
   lsize = 1; if(rank(lvec)>0) lsize = size(lvec);\
   call MPI_SEND(lvec,lsize,MVT,pid,tag,mpi_comm_wav,ierr);\
   call wav_mpi_chkerr(ierr,subName,string);\
END SUBROUTINE
SEND(i,integer(4),MPI_INTEGER)
SEND(l,integer(8),MPI_INTEGER8)
SEND(r,real(4),MPI_REAL4)
SEND(d,real(8),MPI_REAL8)
! PURPOSE: Recv a vector of reals
#define RECV(TAG,VT,MVT)  \
SUBROUTINE wav_mpi_recv##TAG(lvec,pid,tag,string,srcid) ;\
   VT, intent(in) :: lvec(..) ;\
   integer, intent(in) :: pid      ;\
   integer, intent(in) :: tag      ;\
   character(*),optional,intent(in) :: string;\
   integer,optional::srcid ;\
   character(*),parameter           :: subName = "(wav_mpi_recv"// #TAG // ")";\
   integer        :: lsize ,ierr, status(MPI_STATUS_SIZE)  ;\
   lsize = 1; if(rank(lvec)>0) lsize = size(lvec);\
   call MPI_RECV(lvec,lsize,MVT,pid,tag,mpi_comm_wav,status,ierr);\
   if (present(srcid)) srcid=status(MPI_SOURCE);\
   call wav_mpi_chkerr(ierr,subName,string);\
END SUBROUTINE
RECV(i,integer(4),MPI_INTEGER)
RECV(l,integer(8),MPI_INTEGER8)
RECV(r,real(4),MPI_REAL4)
RECV(d,real(8),MPI_REAL8)
#define BCAST(TAG,VT,MVT)  \
SUBROUTINE wav_mpi_bcast##TAG(lvec,root_,string) ;\
   VT, intent(in) :: lvec(..) ;\
   integer,optional::root_ ;\
   character(*),optional,intent(in) :: string;\
   character(*),parameter:: subName = "(wav_mpi_bcast"// #TAG // ")";\
   integer        :: lsize ,ierr,root;\
   root=0;if(present(root_))root=root_;\
   lsize = 1; if(rank(lvec)>0) lsize = size(lvec);\
   call MPI_BCAST(lvec,lsize,MVT,root,mpi_comm_wav,ierr);\
   call wav_mpi_chkerr(ierr,subName,string);\
END SUBROUTINE
! ===============================================================================
!BCAST(s,character(*),MPI_CHARACTER)
!BCAST(g,logical,MPI_LOGICAL)
!BCAST(b,integer(1),MPI_BYTE)
!BCAST(w,integer(2),MPI_INTEGER2)
BCAST(i,integer(4),MPI_INTEGER)
BCAST(l,integer(8),MPI_INTEGER8)
BCAST(r,real(4),MPI_REAL4)
BCAST(d,real(8),MPI_REAL8)
SUBROUTINE wav_mpi_bcastc(vec,root_,string)
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
END SUBROUTINE
SUBROUTINE wav_mpi_bcastp(pk,root_,string)
   type(mpipacket),intent(inout):: pk   ! vector
   integer,optional::root_
   character(*),optional,intent(in):: string   ! message
   character(*),parameter             :: subName = '(wav_mpi_bcastr2) '
   integer               :: ierr
   integer               :: root
   root=0;if(present(root_))root=root_
   call MPI_BCAST(pk%dsize,1,MPI_INTEGER,root,mpi_comm_wav,ierr)
   call MPI_BCAST(pk%Buf,pk%dsize,MPI_PACKED,root,mpi_comm_wav,ierr)
   call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_bcastp
#define REDUCE(TAG,VT,MVT)  \
SUBROUTINE wav_mpi_reduce##TAG(lvec,gvec,op,all,root_,string) ;\
   VT, intent(in) :: lvec(..)  ;\
   VT, intent(out):: gvec(..)  ;\
   integer::op                 ;\
   integer,optional::root_     ;\
   character(*),optional,intent(in) :: string   ;\
   logical,     optional,intent(in) :: all      ;\
   character(*),parameter           :: subName = "(wav_mpi_sum"// #TAG // ") ";\
   logical :: lall;\
   integer :: reduce_type,lsize,gsize,ierr,root;\
   root=0;if(present(root_))root=root_;\
   reduce_type=op;\
   lall = .false.; if (present(all)) lall = all;\
   lsize = 1;gsize=1; if(rank(lvec)>0) then; lsize = size(lvec); gsize = size(gvec);endif;\
   if (lsize /= gsize) call wav_mpi_abort(subName//" lsize,gsize incompatable "//trim(string)) ;\
   if (lall) then ;\
     call MPI_ALLREDUCE(lvec,gvec,gsize,MVT,reduce_type,mpi_comm_wav,ierr) ;\
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_ALLREDUCE",string) ;\
   else ;\
     call MPI_REDUCE(lvec,gvec,gsize,MVT,reduce_type,root,mpi_comm_wav,ierr) ;\
     call wav_mpi_chkerr(ierr,trim(subName)//":"//" MPI_REDUCE",string) ;\
   endif ;\
END SUBROUTINE
REDUCE(i,integer(4),MPI_INTEGER)
REDUCE(l,integer(8),MPI_INTEGER8)
REDUCE(r,real(4),MPI_REAL4)
REDUCE(d,real(8),MPI_REAL8)
#define REDUCEO(TAG,OPN,OP,VT,MVT)  \
SUBROUTINE wav_mpi_##OPN##TAG(lvec,gvec,root_,string,all) ;\
   VT, intent(in) :: lvec(..)  ;\
   VT, intent(out):: gvec(..)  ;\
   integer::op                 ;\
   integer,optional::root_     ;\
   character(*),optional,intent(in) :: string   ;\
   logical,     optional,intent(in) :: all      ;\
   call wav_mpi_reduce##TAG(lvec,gvec,OP,all,root_,string);\
END SUBROUTINE
REDUCEO(i,sum,MPI_SUM,integer(4),MPI_INTEGER)
REDUCEO(l,sum,MPI_SUM,integer(8),MPI_INTEGER8)
REDUCEO(r,sum,MPI_SUM,real(4),MPI_REAL4)
REDUCEO(d,sum,MPI_SUM,real(8),MPI_REAL8)
REDUCEO(i,min,MPI_MIN,integer(4),MPI_INTEGER)
REDUCEO(l,min,MPI_MIN,integer(8),MPI_INTEGER8)
REDUCEO(r,min,MPI_MIN,real(4),MPI_REAL4)
REDUCEO(d,min,MPI_MIN,real(8),MPI_REAL8)
REDUCEO(i,max,MPI_MAX,integer(4),MPI_INTEGER)
REDUCEO(l,max,MPI_MAX,integer(8),MPI_INTEGER8)
REDUCEO(r,max,MPI_MAX,real(4),MPI_REAL4)
REDUCEO(d,max,MPI_MAX,real(8),MPI_REAL8)
! ===============================================================================
SUBROUTINE wav_mpi_commsize(size,string)
  integer,intent(out)                :: size
  character(*),optional,intent(in)   :: string   ! message
  character(*),parameter             :: subName = '(wav_mpi_commsize) '
  integer               :: ierr
  ! PURPOSE: MPI commsize
  call MPI_COMM_SIZE(mpi_comm_wav,size,ierr)
  call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_commsize
SUBROUTINE wav_mpi_commrank(rank,string)
  integer,intent(out)                :: rank
  character(*),optional,intent(in)   :: string   ! message
  character(*),parameter             :: subName = '(wav_mpi_commrank) '
  integer               :: ierr
  ! PURPOSE: MPI commrank
  call MPI_COMM_RANK(mpi_comm_wav,rank,ierr)
  call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_commrank
! ===============================================================================
SUBROUTINE wav_mpi_initialized(flag,string)
  logical,intent(out)                :: flag
  character(*),optional,intent(in)   :: string   ! message
  character(*),parameter             :: subName = '(wav_mpi_initialized) '
  integer               :: ierr
  ! PURPOSE: MPI initialized
  call MPI_INITIALIZED(flag,ierr)
  call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_initialized
! ===============================================================================
SUBROUTINE wav_mpi_abort(string,rcode)
  character(*),optional,intent(in)   :: string   ! message
  integer,optional,intent(in)        :: rcode    ! optional code
  character(*),parameter             :: subName = '(wav_mpi_abort) '
  integer               :: ierr
  integer rcode_
  ! PURPOSE: MPI abort
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
  #ifndef NO_MPI
  call MPI_ABORT(MPI_COMM_WORLD,rcode_,ierr)
  #endif
  stop
END SUBROUTINE wav_mpi_abort
! ===============================================================================
SUBROUTINE wav_mpi_barrier(string)
  character(*),optional,intent(in)   :: string   ! message
  character(*),parameter             :: subName = '(wav_mpi_barrier) '
  integer               :: ierr
  ! PURPOSE: MPI wav_mpi_barrier
  call MPI_BARRIER(mpi_comm_wav,ierr)
  call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_barrier
! ===============================================================================
SUBROUTINE wav_mpi_init(string)
  character(*),optional,intent(in)   :: string   ! message
  character(*),parameter             :: subName = '(wav_mpi_init) '
  integer               :: ierr
  ! PURPOSE: MPI init
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
  character(*),optional,intent(in)   :: string   ! message
  character(*),parameter             :: subName = '(wav_mpi_finalize) '
  integer               :: ierr
  ! PURPOSE: MPI finalize
  call MPI_FINALIZE(ierr)
  call wav_mpi_chkerr(ierr,subName,string)
END SUBROUTINE wav_mpi_finalize
#define GATHER(TAG,VT,MVT,IL)                    \
subroutine wav_mpi_gather##TAG(sbuf,rbuf,root_,string) ;\
  VT ::rbuf(:),sbuf(..);\
  integer,optional::root_;\
  character(*),optional,intent(in)   :: string  ;\
  character(*),parameter             :: subName = "(wav_mpi_gather"// #TAG // ") ";\
  integer :: lsize,ierr, root;\
  root=0;if(present(root_))root=root_;\
  lsize = 1; if(rank(sbuf)>0) lsize = size(sbuf);\
  call mpi_Gather(sbuf,lsize,MVT,rbuf,lsize,MVT,root,mpi_comm_wav,ierr);\
  call wav_mpi_chkerr(ierr,subName,string);\
end subroutine
GATHER(i,integer(4),MPI_INTEGER,4)
GATHER(l,integer(8),MPI_INTEGER8,8)
GATHER(r,real(4),MPI_REAL4,4)
GATHER(d,real(8),MPI_REAL8,8)
#define SCATTER(TAG,VT,MVT,IL)                    \
subroutine wav_mpi_scatter##TAG(sbuf,rbuf,root_,string);\
  VT:: sbuf(:),rbuf(..) ;\
  integer,optional::root_;\
  character(*),optional,intent(in)   :: string  ;\
  character(*),parameter             :: subName = "(wav_mpi_scatter"// #TAG // ") ";\
  integer :: lsize,ierr, root;\
  root=0;if(present(root_))root=root_;\
  lsize = 1; if(rank(sbuf)>0) lsize = size(sbuf);\
  call mpi_scatter(sbuf,lsize,MVT,rbuf,lsize,MVT,root,mpi_comm_wav,ierr);\
  call wav_mpi_chkerr(ierr,subName,string);\
end subroutine
SCATTER(i,integer(4),MPI_INTEGER,4)
SCATTER(l,integer(8),MPI_INTEGER8,8)
SCATTER(r,real(4),MPI_REAL4,4)
SCATTER(d,real(8),MPI_REAL8,8)

#define GATHERV(TAG,VT,MVT,IL)                    \
subroutine wav_mpi_gatherv##TAG(sbuf,scount,rbuf,rcounts,displs,root_,string) ;\
  VT:: sbuf(..),rbuf(..) ;\
  integer ::scount,rcounts(:),displs(:);\
  integer,optional::root_;\
  character(*),optional,intent(in)   :: string  ;\
  character(*),parameter             :: subName = "(wav_mpi_gatherv"// #TAG // ") ";\
  integer :: lsize,ierr, root;\
  root=0;if(present(root_))root=root_;\
  lsize = 1; if(rank(sbuf)>0) lsize = size(sbuf);\
  call MPI_GATHERV(sbuf,scount,MVT,rbuf,rcounts,displs,MVT,root,mpi_comm_wav,ierr);\
  call mpi_Gather(sbuf,lsize,MVT,rbuf,lsize,MVT,root,mpi_comm_wav,ierr);\
  call wav_mpi_chkerr(ierr,subName,string);\
end subroutine

GATHERV(i,integer(4),MPI_INTEGER,4)
GATHERV(l,integer(8),MPI_INTEGER8,8)
GATHERV(r,real(4),MPI_REAL4,4)
GATHERV(d,real(8),MPI_REAL8,8)

#define SCATTERV(TAG,VT,MVT,IL)                    \
subroutine wav_mpi_scatterv##TAG(sbuf,scounts,displs,rbuf,rcount,root_,string) ;\
  VT:: sbuf(..),rbuf(..) ;\
  integer ::scounts(:),displs(:),rcount;\
  integer,optional::root_;\
  character(*),optional,intent(in)   :: string  ;\
  character(*),parameter             :: subName = "(wav_mpi_scatterv"// #TAG // ") ";\
  integer :: lsize,ierr, root;\
  root=0;if(present(root_))root=root_;\
  lsize = 1; if(rank(sbuf)>0) lsize = size(sbuf);\
  call MPI_SCATTERV(sbuf,scounts,displs,MVT,rbuf,rcount,MVT,root,mpi_comm_wav,ierr);\
  call mpi_scatter(sbuf,lsize,MVT,rbuf,lsize,MVT,root,mpi_comm_wav,ierr);\
  call wav_mpi_chkerr(ierr,subName,string);\
end subroutine
SCATTERV(i,integer(4),MPI_INTEGER,4)
SCATTERV(l,integer(8),MPI_INTEGER8,8)
SCATTERV(r,real(4),MPI_REAL4,4)
SCATTERV(d,real(8),MPI_REAL8,8)
! ===============================================================================
SUBROUTINE InitMpiPacket(pk,lsize)
  type( mpipacket) pk
  integer*4,optional::lsize
  pk%id=12340+mpi_id
  if(present(lsize))then
    if(allocated(pk%buf))then
      deallocate(pk%buf)
    endif
    if(lsize>0)then
      pk%bsize=((lsize+8)/4)*4
      allocate(pk%buf(pk%bsize/4+4) )
    endif
  endif
  pk%dsize=0
  pk%pos=0
end SUBROUTINE InitMpiPacket
SUBROUTINE VMpiPacket(pk,lsize)
  type( mpipacket) pk
  integer*4,optional::lsize
  integer*4,allocatable::Buf(:)
  integer mm;
  !print*,"VV",pk%pos,lsize,pk%bsize
  if(pk%pos+lsize>pk%bsize)then
    mm=pk%bsize/4
    !print*,'===C',mpi_id,pk%id,pk%id,mm
    allocate(Buf(mm+4));
    Buf=pk%Buf;
    deallocate(pk%buf);
    pk%bsize=pk%bsize+(lsize/4)*4+10240;
    !print*,'===B',mpi_id,pk%id,pk%bsize/4
    allocate(pk%buf(pk%bsize/4+4) )
    pk%Buf(1:mm)=Buf;
  endif
  !print*,"VE",pk%pos,lsize,pk%bsize
end SUBROUTINE VMpiPacket
SUBROUTINE SetMpiPacketDSize(pk,dsize)
  type( mpipacket) pk
  integer*4,optional::dsize
  if(present(dsize))then
    pk%dsize=dsize
  else
    pk%dsize=pk%pos
  endif
end SUBROUTINE SetMpiPacketDSize
#define PACK(TAG,VT,MVT,IL)                    \
SUBROUTINE wav_mpi_pack##TAG(pk,iobuf,iunpack);\
  type( mpipacket) pk                         ;\
  VT,intent(inout) :: iobuf(..)               ;\
  integer,optional,intent(in)::iunpack        ;\
  integer(4):: iblen,ierr                     ;\
    iblen= 1; if(rank(iobuf)>0) iblen= size(iobuf);\
    if(present(iunpack).and.iunpack/=0)then     ;\
      call MPI_UNPACK(pk%buf,pk%dsize,pk%pos,iobuf,iblen,MVT,mpi_comm_wav,ierr);\
    else;\
      call VMpiPacket(pk,IL*iblen);\
      call MPI_PACK(iobuf,iblen,MVT,pk%buf,pk%bsize,pk%pos,mpi_comm_wav,ierr);\
    endif;\
END SUBROUTINE
PACK(w,integer(2),MPI_INTEGER2,2)
PACK(i,integer(4),MPI_INTEGER,4)
PACK(l,integer(8),MPI_INTEGER8,8)
PACK(r,real(4),MPI_REAL4,4)
PACK(d,real(8),MPI_REAL8,8)
SUBROUTINE wav_mpi_packc0(pk,iobuf,iunpack)
   type( mpipacket) pk
   character*(*), intent(inout) :: iobuf
   integer,optional,intent(in)::iunpack
   integer(4) :: iblen,ierr
   integer:: iunpack_=0
   if(present(iunpack))iunpack_=iunpack
   if(iunpack_/=0)then
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
   integer:: iunpack_=0
   nn=size(iobuf)
   if(present(iunpack))iunpack_=iunpack
   if(iunpack_/=0)then
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

subroutine flush6
  call flush (6)
  call flush (0)
end subroutine flush6
SUBROUTINE SerRunf
  integer:: is,ir,ierr,status(MPI_STATUS_SIZE)  ! mpi status info
  is=0;
  call mpi_Sendrecv(is,1,MPI_INTEGER,PrevProc,1011,ir,1,MPI_INTEGER,NextProc,1011,mpi_comm_wav,status,ierr)
  call mpi_Sendrecv(is,1,MPI_INTEGER,NextProc,1010,ir,1,MPI_INTEGER,PrevProc,1010,mpi_comm_wav,status,ierr)
End SUBROUTINE SerRunf

SUBROUTINE SerRun
  integer:: iv,is,ir,ierr,status(MPI_STATUS_SIZE)  ! mpi status info
  is=0;
  call flush (6)
  call wav_mpi_recv(iv,PrevProc,1010)
  call wav_mpi_send(iv,NextProc,1010)
End SUBROUTINE SerRun

SUBROUTINE zbarr
  integer:: iv,is,ir,ierr,status(MPI_STATUS_SIZE)  ! mpi status info
  is=0;
  call flush (6)
  call wav_mpi_recv(iv,NextProc,1010)
  call wav_mpi_send(iv,PrevProc,1010)
  call wav_mpi_recv(iv,PrevProc,1010)
  call wav_mpi_send(iv,NextProc,1010)
End SUBROUTINE zbarr

SUBROUTINE SerRRun
  integer:: iv,is,ir,ierr,status(MPI_STATUS_SIZE)  ! mpi status info
  is=0;
  call wav_mpi_recv(iv,NextProc,1010)
  call wav_mpi_send(iv,PrevProc,1010)
End SUBROUTINE SerRRun

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
  call flush6
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
!#endif
