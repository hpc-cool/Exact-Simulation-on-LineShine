#ifndef NO_MPI
!#define MPI_GATHERV MYMPI_GATHERV
!#define MPI_SCATTERV MYMPI_SCATTERV
Module zmph_mod
  ! -------------------------------------------------------------------------------
  ! PURPOSE: general layer on MPI functions
  ! -------------------------------------------------------------------------------
#ifdef GF13
  use mpi
#endif
  implicit none
  type MPI_Common
    integer id,npe,comm,pcomm,pcommid
    integer NextProc,PrevProc
    integer igrp
    character*(32) cname;
  end type MPI_Common

  type(MPI_Common),public,pointer ::mpicomms(:)
  character*(MPI_MAX_PROCESSOR_NAME),public:: NodeName

  private
  integer ::ncom=0,mcom=0

  integer::DBGLVL=10
#ifndef GF13
  include 'mpif.h'
#endif
CONTAINS
  ! ===============================================================================
  subroutine resizecom(mcom_)
    integer mcom_
    type(MPI_Common),pointer ::mpicomms_t
    integer i;
    if(mcom_<mcom)return
    if(associated(mpicomms))then
      mpicomms_t=>mpicomms
    endif
    mcom=mcom_;
    allocate(mpicomms(0:mcom))
    if(associated(mpicomms))then
      if(ncom>0)then
        mpicomms(0:ncom-1)=mpicomms_t(0:ncom-1)
      endif
      do i=ncom,mcom-1
        mpicomms(i)%cname='';
        mpicomms(i)%pcommid=0;
        mpicomms(i)%pcomm=MPI_COMM_WORLD;
        mpicomms(i)%comm=MPI_COMM_WORLD;
        mpicomms(i)%id=0
        mpicomms(i)%npe=0
        NextProc=mpi_proc_null
        PrevProc=mpi_proc_null
      enddo
      deallocate(mpicomms_t);
    endif
  end subroutine resizecom
  subroutine zmph_init(mode )
    integer mode
    integer ierr,nl,cid
    logical Inited
    ! int COUPLE_CCSM   Init Mpi In msg_pass('connect') & msg_pass('disconnect'
    if(mode<0)then
      ! ============  Shut down Parallel environment
      call MPI_INITIALIZED(Inited,ierr)
      if(Inited)then
        call wav_mpi_finalize('Deinit')
      endif
    else
      ! ===========  Initialise Parallel environment
      if(ncom>=1)return
      call MPI_INITIALIZED(Inited,ierr)
      if(.not. Inited)then
        call MPI_INIT(ierr)
      endif
      call MPI_Get_processor_name(NodeName,nl,ierr);
      call resizecom(1)
      mpicomms(0)%comm=MPI_COMM_WORLD
      call MPI_COMM_RANK(mpicomms(0)%comm,mpicomms(0)%id,ierr)
      call MPI_COMM_SIZE(mpicomms(0)%comm,mpicomms(0)%npe,ierr)
      cid=zmph_components_file_id(
    end if
  end subroutine zmph_init
  subroutine zmph_Global()
    end
    integer function FindByComm(comm)
      integer comm
      integer i
      FindByComm=0
      do i=0,ncom-1
        if(mpicomms(i)%comm==comm)then
          FindByComm=i
          exit
        endif
      enddo
    end function FindByComm

    integer function FindByName(name,pcomm)
      integer comm
      character*(*)name
      integer i
      FindByName=-1
      do i=0,ncom-1
        if(mpicomms(i)%pcomm==pcomm.and. mpicomms(i)%name==name )then
          FindByName=i
          exit
        endif
      enddo
    end function FindByName
    integer function FindByGrp(name,pcomm)
      integer comm,igrp
      integer i
      FindByGrp=-1
      do i=0,ncom-1
        if(mpicomms(i)%pcomm==pcomm.and. mpicomms(i)%igrp==igrp )then
          FindByGrp=i
          exit
        endif
      enddo
    end function FindByGrp
    integer  function zmph_components_id(name,pid_)
      character*(*) name
      integer,optional::pid_
      character*(32),pointer::names(:), nnames(:)
      integer ,pointer::igrps(:)
      integer iGrp,ngrp,i,j
      integer pid,cid
      call zmph_init(1);
      pid=0
      if(present(pid_))pid=pid_
      if(pid==-1)pid=ncom-1; ! error
    else if(pid<0)return pid; ! error
    cid=FindByName(name,mpicomms(pid)%comm)
    if(cid>=0)return cid;
    ncom=ncom+1;  cid=ncom-1;
    mpicomms(cid)%name=name
    allocate(names(mpicomms(pid)%npe))
    allocate(nnames(mpicomms(pid)%npe))
    allocate(iGrps(mpicomms(pid)%npe))
    call mpi_Gather(mpicomms(cid)%name,32,MPI_CHARACTER,names,32,MPI_CHARACTER,0,comm)
    if(mpicomms(pid)%id==0)then
      ngrp=0;
      do i=1,mpicomms(pid)%npe
        iGrp=0
        do j=1,ngrp
          if(names[i]==nnames[j])then
            iGrp=j;exit
          endif
        enddo
        if(iGrp==0)then
          ngrp=ngrp+1
          nnames[ngrp]=names[i]
          iGrp=ngrp
        endif
        iGrps(i)=iGrp
      enddo
    endif
    call MPI_SCATTER(iGrps,1,mpi_integer,iGrp,1,mpi_integer,0,comm)
    mpicomms(cid)%igrp=igrp
    mpicomms(cid)%pcomm=mpicomms(pid)%comm
    call MPI_COMM_SPLIT(mpicomms(cid)%pcomm,mpicomms(cid)%igrp,0,mpicomms(cid)%comm,ierr);
    call MPI_COMM_RANK (mpicomms(cid)%comm ,mpicomms(cid)%id,ierr)
    call MPI_COMM_SIZE (mpicomms(cid)%comm ,mpicomms(cid)%npe,ierr)
    deallocate(iGrps,names,nnames)
  end function zmph_components_id
  integer function zmph_components(name,pcomm_)
    character*(*) name
    integer,optional::pcomm_
    integer pcomm,pid,cid
    call zmph_init(1);
    pcomm=MPI_COMM_WORLD
    if(present(pcomm_))pcomm=pcomm_
    pid=FindByComm(pcomm)
    if(pid<0)return pid; ! error
    cid=zmph_components_id(name,pid)
    return mpicomms(cid)%comm
  end function zmph_components
  integer function zmph_components_file_id(filename){
    character*(*) filename
    integer ios,i,cid
    character*(32) tag
    call zmph_init(1);
    open(10,file=filename,status='old',IOSTAT=ios)
    if(ios/=0)then
      zmph_components_file_id=ncom-1
      rturn;
    endif
    read(10,*)tag
    if(tag=='ZMPH_SPLIT_DEF')then
      do i=1,100
        read(10,*)tag
        if(tag=='ZMPH_SPLIT_DEF_END')exit
        cid=zmph_components_id(tag,cid)
      enddo
    endif
    close(10);
    zmph_components_file_id=mpicomms(cid)%comm
  end function zmph_components_file_id
  integer function zmph_components_file(filename){
    character*(*) filename
    integer ios,i,cid
    cid=zmph_components_file_id(filename)
    zmph_components_file=mpicomms(cid)%comm
  end function zmph_components_file
END MODULE zmph_mod
#endif
