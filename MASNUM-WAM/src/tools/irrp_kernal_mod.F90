!#define MPI_GATHER   MYMPI_GATHER
!#define MPI_GATHERV  MYMPI_GATHERV
!#define MPI_SCATTER  MYMPI_SCATTER
!#define MPI_SCATTERV  MYMPI_SCATTERV

#define DBGINF   write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__
#define DBGB0(msg)    if(gsi%pid==0)then; write(6,'(i4,f8.3,"s ",i7," ",a)') gsi%pid,Difftimer(),iwalltime(),msg;call flush(6);endif
!#define DBGP if(gsi%pid==509)print*,'AAAAA ',gsi%pid,__LINE__

!#################################################################################################
!-------------------------------------------------------------------------------------------------
!*************************************************************************************************
!-------------------------------------------------------------------------------------------------
!                                                     Copyright (C) 2013 Xunqiang Yin & Wei Zhao
!                                                     MODULE NAME : irrp_kernal_mod
!                                                     PRESENT VERSION : 2013-04-30
!
! --- USAGE : To convert time among different types: datestr, datevec, datenum.
! --- DEPEND: irrp_smpi_mod, irrp_pemask_mod
!
! --- NOTE for describing of subroutine / function :
!  A. The parameters bracketed with [], means optional parameter.
!  B. The describe for the parameters of subroutine / function, started with:
!   * It means input prameter;
!   # It means output prameter;
!   @ It means input and output prameter(it will be changed inside).
!
!-------------------------------------------------------------------------------------------------
! ***                                 INTERFACE DESCRIBE                                       ***
!-------------------------------------------------------------------------------------------------
!
!  1. sub. irrp_set_pemask : Do irregular partition using MASK and return the results in PEMASK.
!
!     irrp_part_init(partmode, npe, pid, mpi_comm, mask, halosize, cycle_flag, scycle)
!
!    * integer :: mask_   = Mask for horizonal space/points in which 0 means land/useless points.
!    # integer :: pemask_ = Arrary to save partition reaults each point with the value equal to
!                           the ID of the PE that this point is belonged. This ID of PEs is
!                           started from 1 which is different with the way for MPI.
!    * integer :: im_     = The first dimension size of mask_ and pemask_.
!    * integer :: jm_     = The second dimension size of mask_ and pemask_.
!    * integer :: npe_    = The total number of PEs used for this partition.
!
!-------------------------------------------------------------------------------------------------
!*************************************************************************************************
!-------------------------------------------------------------------------------------------------

  module irrp_kernal_mod

!-------------------------------------------------------------------------------------------------

  use irrp_smpi_mod
  use irrp_split_mod
  !use debughlp_mod
  implicit none

!-------------------------------------------------------------------------------------------------

  public :: irrp_part_init, irrp_part_final,gs2ij
  public :: irrp_SetPartMatrix, irrp_SetPartSerial,irrp_getrects

  public :: pi_pos_type, pi_rsinfo_type, nb_8pnts_def_type
  public :: globle_info_type,pi_partinfo_type,pi_rsbuf_type,pi_rplist_type
  public :: gsi, cpart

  private
  !integer ierr
  !call MPI_BARRIER(gsi%mpi_comm,ierr)

!-------------------------------------------------------------------------------------------------

!  type pi_pos_type                          ! Self-defined type to record 2 dimensional index.
!                                            ! This will be fast than 2 dimensional arrary.
!    integer :: i                            ! Index of the first dimension
!    integer :: j                            ! Index of the second dimension
!  end type pi_pos_type

  type nb_8pnts_def_type
    integer :: r                            ! Right neighbor point.
    integer :: ur                           ! Up-right neighbor point.
    integer :: u                            ! Up neighbor point.
    integer :: ul                           ! Up-left neighbor point.
    integer :: l                            ! Left neighbor point.
    integer :: dl                           ! Down-left neighbor point.
    integer :: d                            ! Down neighbor point.
    integer :: dr                           ! Down-right neighbor point.
  end type nb_8pnts_def_type

  type pi_rsbuf_type
    integer          :: nrb                 ! Size of receive buffer.
    integer          :: nsb                 ! Size of send buffer.
    integer, pointer :: rbuf(:)             ! (nrb), buffer for receiving.
    integer, pointer :: sbuf(:)             ! (nsb), buffer for sending.
  end type pi_rsbuf_type

  type pi_rsinfo_type
    integer :: id                           ! ID of the destination PE.
    integer :: sn                           ! Number of sending points.
    integer :: rn                           ! Number of receiving points.
    integer, pointer :: sspnts(:)           ! List of sending   points for serial part.
    integer, pointer :: srpnts(:)           ! List of receiving points for serial part.
                                            ! on part End :spnts==sspnts rpnts==srpnts
    integer, pointer :: mspnts(:)           ! List of sending   points for matrix part.
    integer, pointer :: mrpnts(:)           ! List of receiving points for matrix part
    integer, pointer :: spnts(:)            ! List of sending   points for current part.
    integer, pointer :: rpnts(:)            ! List of receiving points for current part
  end type pi_rsinfo_type

!-------------------------------------------------------------------------------------------------

  type globle_info_type                     ! Self-defined type to record partition information.
    integer          :: mpi_comm   = MPI_COMM_WORLD
    integer          :: parttype   = -1     ! Default is for serial, 1 for matrix.
                                            ! model part mem shape,
                                            !       0: serial without dump points;
                                            !       1: serial with dump points,for nopack sending;
                                            !       2: rectangle shape
    integer          :: balance    = 0      ! Type of balance. Default is 0. 1: absolute balance.
    integer          :: npe        = 0      ! Numbers of PEs.
    integer          :: pid        = 0      ! Current PE id.
    integer          :: halosize   = 1      ! size of halo (outer boundary)
    integer          :: cycle_flag = 0      ! Cycle type for 4 boundaries of grid matrix.
    integer          :: im, jm              ! Size of grid matrix, first & second dimension.
    integer          :: ijm                 ! Maximum computing points.
    integer          :: sumnp               ! Sum of all points.
    real(8)          :: avenp               ! Averaged points number for each PE.
    integer          :: gnpc                ! Total computing points. It may include dump points.
    integer, pointer :: npcs(:)             ! (npe)
    integer, pointer :: mask(:, :)          ! (im, jm)
    integer, pointer :: rects(:, :)         ! (4, npe)
                                            !-----------------------------------------------------
                                            ! Temp Array no needed when partion end.
    integer                    :: mnprs     ! Buf size for partion exchange Data
                                            !-----------------------------------------------------
  end type globle_info_type

!-------------------------------------------------------------------------------------------------
  type pi_rplist_type
    integer gis,i,j,pe
  end type pi_rplist_type

  type pi_partinfo_type
    integer :: snpc                         ! All computing points os all partions.
    integer :: snp                          ! All points of all patitions.
    integer :: nps                          ! All sending pnts of current partition.
    integer :: npr                          ! All receiving points of current partion.
    integer :: npc                          ! All conputing points of current partition.
    integer :: np                           ! All points for current partion.
    integer :: iblv = 1                     ! The start of sequentializing: from 1 or 0.
    integer :: nx, ny                       ! Size of matrix: first and second dimension.
    integer :: nxy                          ! nxy = nx * ny.
    !---------------------------------------------------------------------------------------------
    ! pnttype:
    ! serial ,nodump :(0:np),0:space/land,1:open bound             3:calc,4:recv1, 4+:recv1+
    ! serial ,dump   :(0:np),0:space/land,1:open bound 2:dump pnts,3:calc,4:recv1, 4+:recv1+
    ! matrix : (0:im*jm)
    !---------------------------------------------------------------------------------------------
    integer, pointer :: spnttype(:)         ! Point type for serial part
    integer, pointer :: mpnttype(:)         ! Point type for matrix part
    integer, pointer :: pnttype(:)          ! Point type for current part
    integer, pointer :: calpnts(:)          ! (snp), ! for matrix part
                                            ! calpnts SN;
                                            !      on parttype=0,noneeded,calpnts(1:snpc)=1:snpc
                                            !      othertype Need to Set for gather/scatter
                                            !   sort:nps,npi=snpc-nps,npr1,npr2 ...
                                            !  needed for gather/scatter
    integer          :: recti(4)            ! rectangle of the partion only for computing points.
    integer          :: recto(4)            ! rectangle of the partion for all points.
    integer          :: ib ,ie ,jb ,je
    integer          :: ibi,iei,jbi,jei
    !ZNOUSE integer, pointer :: ij2s(:, :)
    integer, pointer :: ij2s(:, :)! (1-halosize:im+halosize, 1-halosize:jm+halosize)
    type(nb_8pnts_def_type), pointer :: nb8(:)
                                            !(0:np), to record neighbor points for each points.
                                            !-----------------------------------------------------
                                            ! neibhor part info needed for exchange
    integer          :: inited_rsv = 0      ! Check flag for initializing of following arraries.
    integer          :: nids                ! Number of PEs for sending or receiving.
    integer, pointer :: ppos(:)             ! (0:np)
    integer, pointer :: requests(:)         ! (0:nids)
    integer, pointer :: requestr(:)         ! (0:nids)
    integer, pointer :: statuss(:, :)       ! (MPI_STATUS_SIZE,0:nids)
    integer, pointer :: statusr(:, :)       ! (MPI_STATUS_SIZE,0:nids)
    type(pi_rsinfo_type), pointer :: rsinfo(:)     ! (0:nids)
    type(pi_rsbuf_type) , pointer :: rsbufs(:)     ! (0:nids), buffer for recv and send.
                                            !-----------------------------------------------------
    type(pi_rplist_type), pointer :: rplist(:) ! (nxy)
  end type pi_partinfo_type

!-------------------------------------------------------------------------------------------------

  type(Globle_Info_type) :: gsi
#ifdef SERIAL_TEST
  type(pi_partinfo_type), allocatable :: parts(:)
#else
  type(pi_partinfo_type) :: cpart
#endif

  integer, allocatable :: tempmask(:,:)

!-------------------------------------------------------------------------------------------------

  contains
  subroutine pmask(tag,mask)
    character*(*) tag
    integer, intent(in),target :: mask(:, :)
    character*256 fn
    integer ix,jy
    write(fn,'(a,i2.2,".txt")')trim(tag),gsi%pid
    open(unit=1006,file=trim(fn),status='unknown')
    do jy=gsi%jm-1,1,-1
      write(1006,'(1000i2)')(mask(ix,jy),ix=1,gsi%im)
    enddo
    close(1006)
  end subroutine 
  ! getrects
  ! Added by zhaowei for get rect & rects
  subroutine irrp_getrects(recti,recto,rects)
    integer,intent(out), optional:: recti(4)
    integer,intent(out), optional:: recto(4)
    integer,intent(out), pointer, optional:: rects(:,:)
    if(present(recti))recti = cpart%recti
    if(present(recto))recto = cpart%recto
    if(present(rects))then
      allocate(rects(4,gsi%npe))          ! allocate rects
      rects=gsi%rects;
    endif
  end subroutine irrp_getrects
  !-----------------------------------------------------------------------------------------------
  !SetPartMatrix
  !return :
  !   recto:
  !   ptype: pointer Type
  subroutine irrp_SetPartMatrix(recto, ptype)
    integer,          intent(out), optional:: recto(4)
    integer, pointer, intent(out), optional:: ptype(:, :)
    type(pi_rsinfo_type), pointer :: rsi
    integer, pointer :: pnts(:)
    integer n,i,j,i1,j1,ij,ni
    i1 = cpart%recto(1); j1 = cpart%recto(2)
    if(gsi%parttype /= 1)then
      gsi%parttype = 1
      cpart%nx  = cpart%recto(3) - cpart%recto(1) + 1
      cpart%ny  = cpart%recto(4) - cpart%recto(2) + 1
      cpart%nxy = cpart%nx * cpart%ny
      cpart%np  = cpart%nxy
      cpart%npc = cpart%nxy
      cpart%iblv=1
      if(associated(cpart%calpnts))deallocate(cpart%calpnts)
      if(associated(cpart%mpnttype))deallocate(cpart%mpnttype)
      allocate(cpart%calpnts(cpart%snp))
      allocate(cpart%mpnttype(cpart%nxy))
      cpart%mpnttype = 0
      cpart%calpnts = 0
      do n = 1, cpart%snp
        ij = cpart%ppos(n)
        call gs2ij(ij,i,j)
        ij = i - i1 + 1 + (j - j1) * cpart%nx
        cpart%calpnts(n) = ij
        cpart%mpnttype(ij) = 1                      !?????
        if(n > cpart%snpc)cpart%mpnttype(ij) = 3    !?????
        if(n > cpart%nps)cpart%mpnttype(ij) = 2     !?????
      enddo
      do ni=1,cpart%nids
        rsi=>cpart%rsinfo(ni)
        if(associated(rsi%mrpnts))deallocate(rsi%mrpnts)
        if(associated(rsi%mspnts))deallocate(rsi%mspnts)
        allocate(rsi%mrpnts(rsi%rn))
        allocate(rsi%mspnts(rsi%sn))
        do n = 1, rsi%rn
          ij = cpart%ppos(rsi%srpnts(n))
          call gs2ij(ij,i,j)
          rsi%mrpnts(n)=i-i1+1+(j-j1)*cpart%nx
        enddo
        do n=1,rsi%sn
          ij=cpart%ppos(rsi%sspnts(n))
          call gs2ij(ij,i,j)
          rsi%mspnts(n)=i-i1+1+(j-j1)*cpart%nx
        enddo
        rsi%rpnts => rsi%mrpnts
        rsi%spnts => rsi%mspnts
      enddo
    endif
    if(present(recto))recto = cpart%recto
    if(present(ptype))then
      allocate(ptype(cpart%nx, cpart%ny))
      do j = 1, cpart%ny
      do i = 1, cpart%nx
        ij = i  + (j-1) * cpart%nx
        ptype(i, j) = cpart%mpnttype(ij)
      enddo
      enddo
    endif
  end subroutine irrp_SetPartMatrix

  !-------------------------------------------------------------------------------------------------
  !SetPartSerial
  !input:
  ! iblv_:The start of sequentializing: from 1 or 0.
  !return :
  !   gnpc:number of Calpnts in all part
  !   npc:number of Calpnts in cur part
  !    np:number of pnts(Include recv pnts) in cur part
  ! plist:attribs of pnts( np)
  ! nb8  : neibhor of pnts( npc)
  subroutine irrp_SetPartSerial(gnpc, npc, np, plist, nb8, iblv_,nwps_)
    integer                ,intent(out), optional :: npc, np, gnpc,nwps_
    integer                ,intent(in),  optional :: iblv_
    type(pi_pos_type)      ,intent(out), optional,  pointer:: plist(:)
    type(nb_8pnts_def_type),intent(out), optional,  pointer:: nb8(:)
    type(pi_rsinfo_type), pointer :: rsi
    integer, pointer :: pnts(:)
    integer :: ni, n,inb
    integer,allocatable :: itmp(:) ! yinxq 2015-4-4 16:39:49
    !DBGB0('irrp_SetPartSerial 1')
    if(gsi%parttype /= 0)then
      gsi%parttype = 0
      cpart%np = cpart%snp ;
      cpart%npc = cpart%snpc
      if(associated(cpart%calpnts))deallocate(cpart%calpnts)
      do ni = 1, cpart%nids
        rsi => cpart%rsinfo(ni)
        if(associated(rsi%mrpnts))deallocate(rsi%mrpnts)
        if(associated(rsi%mspnts))deallocate(rsi%mspnts)
        rsi%rpnts => rsi%srpnts
        rsi%spnts => rsi%sspnts
          !    :: spnttype(:) ! for matrix part
          !    :: mpnttype(:) ! for matrix part
          !    :: pnttype(:)  ! for matrix part
      enddo
    endif
    !DBGB0('irrp_SetPartSerial 2')
    if(present(npc )) npc  = cpart%snpc
    if(present(np  )) np   = cpart%snp
    if(present(gnpc)) gnpc = gsi%gnpc
    if(present(plist))then
      if(associated(plist))deallocate(plist)
      allocate(plist(0:cpart%np));

      do n = 1, cpart%np
        call gs2ij(cpart%ppos(n),plist(n)%i,plist(n)%j)
      enddo

    endif
    if(present(nb8))then
      if(associated(nb8))deallocate(nb8)
      allocate(nb8(0:cpart%np)); nb8 = cpart%nb8
    endif
    cpart%iblv=1
    if(present(iblv_))then
      if(iblv_==0)then
        cpart%iblv=0
      else
        cpart%iblv=1
      endif
    endif
    if(present(nwps_))then
      !     nwps_=0
      !     do inb=1,cpart%nids
      !       nwps_=nwps_+cpart%rsinfo(inb)%sn
      !     enddo
      allocate(itmp(cpart%nxy));itmp=0 ! yinxq 2015-4-4 16:45:00
      do inb=1,cpart%nids
        do n=1,cpart%rsinfo(inb)%sn
          itmp(cpart%rsinfo(inb)%sspnts(n))=1
        enddo
      enddo
      nwps_=sum(itmp);deallocate(itmp)  ! yinxq 2015-4-4 16:44:54
    endif
    DBGB0('irrp_SetPartSerial End')
  end subroutine irrp_SetPartSerial

  !-------------------------------------------------------------------------------------------------
  ! init Part:
  ! input:
  !  partmode: 0:Normal part;other part method,Input pemask
  !  npe: number of PE
  !  pid: cur pid
  !  mpi_comm:mpi comm
  !    mask_:Mask for horizonal space/points in which 0 means land/useless points.????
  ! halosize:
  ! cycle_flag: Open boundary type.
  !   0 - None cycle boundary.
  !   1 - Cycled boundary in x-direction for (normal/displaced) polar grid.
  !   2 - Cycled boundary in x-direction and y-direction for triple polar grid.
  !   3 - Cycled boundaries are given manually.
  ! scycle:
  !
  !-------------------------------------------------------------------------------------------------

  subroutine irrp_part_init(partmode, npe, pid, mpi_comm, mask, halosize, cycle_flag, scycle)
    integer, intent(in)           :: partmode, npe, pid, mpi_comm
    integer, intent(inout),target :: mask(:, :)
    integer, intent(in), optional :: halosize, cycle_flag
    integer, intent(in), optional :: scycle(:, :)
                                            ! (2,(im+jm+2)*2)
    integer :: i, j, ij, nn, m
    !integer, allocatable :: pemask(:, :)    ! (im, jm), same shape as mask.
    call starttimer

    gsi%npe      = npe                      ! Number of PEs.
    gsi%pid      = pid                      ! Current ID of PE, start from 0.
    gsi%mpi_comm = mpi_comm                 ! Common world of MPI.
    gsi%im       = size(mask, 1)            ! Size of grid matrix: first dimension.
    gsi%jm       = size(mask, 2)            ! Size of grid matrix: second dimension.
    gsi%ijm      = gsi%im*gsi%jm            ! Maximum computing points.
    gsi%halosize = 1                        ! Halo size, default is 1.
    if(present(halosize))gsi%halosize = halosize
                                            ! Record halo size given from outside.
    !-------------------------------------------------------------------------------------------
    ! cycle flag is given from outside.
    !   0 - None cycle boundary.
    !   1 - Cycled boundary in x-direction for (normal/displaced) polar grid.
    !   2 - Cycled boundary in x-direction and y-direction for triple polar grid.
    !   3 - Cycled boundaries is given manually.
    !-------------------------------------------------------------------------------------------
    if(present(cycle_flag))then             ! Record cycle flag: default is 0.
      gsi%cycle_flag = cycle_flag
    elseif(present(scycle))then
      gsi%cycle_flag = 3                    ! Determined by scycle:
                                            ! Cycled boundaries are given manually.
    else
      gsi%cycle_flag = 0                    ! Default setting: None cycle boundary.
    endif
    if(gsi%cycle_flag /= 0 .and. gsi%halosize /= 1)then
                                            ! Check the matchment of cycle flag and halo size.
      write(6,*) "The case of halosize /=1 not surport yet when cycle_flag /=0. "
      call irrp_abort(__FILE__, __LINE__)
    endif

    DBGB0('irrp_part_init begin')
    allocate(gsi%rects(4,gsi%npe))          ! allocate rects
    allocate(gsi%npcs(gsi%npe))             ! allocate npcs
    !DBGB0('irrp_part_init 1')
    nn = 0
    do j = 1, gsi%jm
      do i = 1, gsi%im                        ! sequentializing points from matrix form.
        if(mask(i, j) > 0)nn = nn + 1;
      enddo
    enddo
    !if(pid==0) print*,"gnpc=",nn,"========="
    !DBGB0('irrp_part_init 2 ')
    gsi%sumnp = nn                          ! Set number of sequentialized points.
    gsi%avenp = gsi%sumnp / float(gsi%npe)  ! Compute averaged number of points for each PE.
    gsi%mnprs = (gsi%avenp + 1 + 4) + 4     ! ??
    DBGB0('init_part begin')
#ifdef SERIAL_TEST
    allocate(parts(gsi%npe))                ! allocate parts.
    do m = 1, gsi%npe
      call init_part(parts(m))              ! initialize parts for all partion by serial.
    enddo
#else
    call init_part(cpart)                   ! initialize cpart.
#endif
    DBGB0('init_part end')
    !call pmask("mki",mask);
    if(partmode == 0)then                   ! Do partion using irrp_pemask_mod
    !      allocate(tempmask(gsi%im,gsi%jm)); tempmask = mask
      DBGB0('irrp_set_pemask begin')
      !allocate(pemask(gsi%im, gsi%jm))      ! allocate for pemask.
      call irrp_set_pemask(mask,  gsi%im, gsi%jm, gsi%npe,gsi%pid,gsi%mpi_comm)
      gsi%balance = 1
      DBGB0('irrp_set_pemask end')
      !mask = pemask
      !deallocate(pemask)                    ! release pemask.
    else
      gsi%balance = 0
    endif
    gsi%mask=>mask
    !call pmask("mko",mask);
    !call set_cycle_bound(scycle)            ! Deal with cycle boundaries.
#ifdef SERIAL_TEST
    DBGB0('SetPartInf')
    call SetPartInf                     ! Set rects and npcs.Set information for each partition.
    DBGB0('set_part_ij2s')
    do m = 1, gsi%npe
      call set_part_ij2s(m,parts(m),scycle)
    enddo
    DBGB0('get_nbp_recv')
    do m = 1, gsi%npe
      call get_nbp_recv(m, parts(m))        ! Set recv points for each partition.
    enddo
    DBGB0('SetPnts')
    do m = 1, gsi%npe
      call SetPnts(m, parts(m))             ! Set all points for each partition
    enddo
#else
    DBGB0('SetPartInf')
    call SetPartInf(gsi%pid+1, cpart)        ! Set information for current partition.
    DBGB0('set_part_ij2s')
    call set_part_ij2s(gsi%pid+1,cpart,scycle)
    DBGB0('get_nbp_recv')
    call get_nbp_recv(gsi%pid+1, cpart)     ! Set recv points for current partition.
    DBGB0('SetPnts')
    call SetPnts(gsi%pid+1, cpart)          ! Set all points for current partition.
    DBGB0('SetPnts End ')
#endif
    gsi%gnpc = sum(gsi%npcs)                ! Number of all computer points for all PEs.
    DBGB0('irrp_SetPartSerial')
    call irrp_SetPartSerial                 ! Prepare point list in serial case as default.
                                            ! It can be changed into matrix case later.
    DBGB0('End irrp_part_init')

  end subroutine irrp_part_init


!-------------------------------------------------------------------------------------------------
  ! Init part struct
  subroutine init_part(part)
    type(pi_partinfo_type), intent(inout) :: part
    nullify(part%ppos    )
    nullify(part%nb8     )
    nullify(part%rsinfo  )
    nullify(part%requests)
    nullify(part%requestr)
    nullify(part%statuss )
    nullify(part%statusr )
    nullify(part%rsbufs  )
    part%nids       = 0
    part%npr        = 0
    part%nps        = 0
    part%snp        = 0
    part%snpc       = 0
    part%inited_rsv = 0
  end subroutine init_part

!-------------------------------------------------------------------------------------------------
  ! End Part
  ! mode: 0:part end ,leave vars for use
  !       1:free all var,not use ????
  subroutine irrp_part_final(mode)
    integer, optional, intent(in) :: mode
#ifdef SERIAL_TEST
    integer :: m
    if(allocated(parts))then
      do m = 1, gsi%npe
        call release_part(parts(m), mode)
      enddo
      if(present(mode))then
        if(mode > 0)deallocate(parts)
      endif
    endif
#else
    call release_part(cpart, mode)          ! Finalize cpart.
#endif
    if(present(mode))then
      if(mode == 0)return
      if(associated(gsi%rects  ))deallocate(gsi%rects  )
      if(associated(gsi%npcs   ))deallocate(gsi%npcs   )
    endif
  end subroutine irrp_part_final

  !-----------------------------------------------------------------------------------------------
  ! sub. release_part: release buffers and point lists in part.
  !
  !    release_part(part, mode)
  !
  !   @ type(pi_partinfo_type) :: part = partinfo need to be released.
  !   # integer                :: mode = method of release. If not given, only release buffers.
  !
  !-----------------------------------------------------------------------------------------------

  subroutine release_part(part, mode)
    type(pi_partinfo_type), intent(inout)         :: part
    integer,                intent(in),  optional :: mode
    integer :: i
    if(associated(part%rsbufs))then
      do i = 1, part%nids
        if(associated(part%rsbufs(i)%rbuf))deallocate(part%rsbufs(i)%rbuf)
        if(associated(part%rsbufs(i)%sbuf))deallocate(part%rsbufs(i)%sbuf)
      enddo
      deallocate(part%rsbufs)
    endif
    if(present(mode))then                   ! If mode is not given or equal to 0,
      if(mode == 0)return                   ! only release those rs buffers.
      if(associated(part%rsinfo))then
        do i = 1, part%nids
          if(associated(part%rsinfo(i)%mrpnts))deallocate(part%rsinfo(i)%mrpnts)
          if(associated(part%rsinfo(i)%mspnts))deallocate(part%rsinfo(i)%mspnts)
          if(associated(part%rsinfo(i)%srpnts))deallocate(part%rsinfo(i)%srpnts)
          if(associated(part%rsinfo(i)%sspnts))deallocate(part%rsinfo(i)%sspnts)
        enddo
      endif
      if(associated(part%rsinfo  ))deallocate(part%rsinfo  )
      if(associated(part%ppos    ))deallocate(part%ppos    )
      if(associated(part%spnttype))deallocate(part%spnttype)
      if(associated(part%mpnttype))deallocate(part%mpnttype)
      if(associated(part%nb8     ))deallocate(part%nb8     )
      if(associated(part%requests))deallocate(part%requests)
      if(associated(part%requestr))deallocate(part%requestr)
      if(associated(part%statuss ))deallocate(part%statuss )
      if(associated(part%statusr ))deallocate(part%statusr )
    endif
  end subroutine release_part

  !-----------------------------------------------------------------------------------------------
  ! sub. SetPartInf: Prepare rects and npcs for all PEs.
  !-----------------------------------------------------------------------------------------------

#ifdef SERIAL_TEST
  subroutine SetPartInf
    integer :: n, m, i, j
    type(pi_partinfo_type),pointer :: cpart
    gsi%rects = -1                          ! Initiale value for rects.
    gsi%npcs = 0                            ! Initiale value for npcs.
    do j = 1, gsi%jm
    do i = 1, gsi%im                        ! sequentializing points from matrix form.
      m = gsi%mask(i,j)                    ! Get ID for this point.
      if(m<=0)cycle
      if(gsi%rects(1, m) < 0)then
         gsi%rects(1, m) = i                ! First guess of minimum i.
         gsi%rects(2, m) = j                ! First guess of minimum j.
         gsi%rects(3, m) = i                ! First guess of maximum i.
         gsi%rects(4, m) = j                ! First guess of maximum j.
      else                                  ! Adjust for minimum and maximum i/j.
        if(gsi%rects(1, m) > i)gsi%rects(1, m) = i
        if(gsi%rects(2, m) > j)gsi%rects(2, m) = j
        if(gsi%rects(3, m) < i)gsi%rects(3, m) = i
        if(gsi%rects(4, m) < j)gsi%rects(4, m) = j
      endif
      gsi%npcs(m) = gsi%npcs(m) + 1         ! Accumulation for npcs.
      enddo
    enddo
    do m=1,gsi%npe
      part=>parts(m)
      nullify(part%ppos, part%nb8, part%rsinfo)
      nullify(part%requests, part%requestr)
      nullify(part%statuss, part%statusr, part%rsbufs)
      part%snp      = 0
      part%nps      = 0
      part%npr      = 0
      part%nids     = 0
      part%snpc     = gsi%npcs(m)
      part%recti    = gsi%rects(:, m )
      part%rect     = gsi%rects(:,m)
      part%snpc     = gsi%npcs(m)
      part%recto(1) = part%recti(1) - gsi%halosize
      part%recto(2) = part%recti(2) - gsi%halosize
      part%recto(3) = part%recti(3) + gsi%halosize
      part%recto(4) = part%recti(4) + gsi%halosize
      part%nx       = part%recto(3) - part%recto(1) + 1
      part%ny       = part%recto(4) - part%recto(2) + 1
      part%nxy      = part%nx       * part%ny
    enddo

    !---------------------------------------------------------------------------------------------
    ! Check balance and points omitted.
    ! Only the first PE do the following checks.
    !---------------------------------------------------------------------------------------------
    if(gsi%pid == 0)then
      if(gsi%balance == 1)then               ! This check is only for irregular partition.
        do n = 1, gsi%npe                   ! Check absolute balance for each PE.
          if(abs(gsi%npcs(n) - gsi%avenp) > 1)then
            if(gsi%pid == 0)then
              write(6, *)'Error A Not Absolute balance ', n, gsi%npcs(n), gsi%avenp, gsi%pid
            endif
            !call irrp_abort(__FILE__, __LINE__)
          endif
        enddo
      endif
    endif
  end subroutine SetPartInf
#else
  subroutine SetPartInf(m,part)
    integer, intent(in) :: m
    type(pi_partinfo_type), intent(inout) :: part
    integer :: n, mt, i, j,ierr
    nullify(part%ppos, part%nb8, part%rsinfo)
    nullify(part%requests, part%requestr)
    nullify(part%statuss, part%statusr, part%rsbufs)
    part%snp      = 0
    part%nps      = 0
    part%npr      = 0
    part%nids     = 0
#if 0
    call countrect(gsi%mask,m,gsi%im,gsi%jm,part%recti,part%snpc)
#else
    part%snpc=0
    part%recti=-1
    !print*,gsi%pid,i1,i2,gsi%im
    do j = 1, gsi%jm
      do i = 1, gsi%im                        ! sequentializing points from matrix form.

        if(m == gsi%mask(i,j))then                    ! Get ID for this point.
          if(part%recti(1) < 0)then
             part%recti(1) = i                ! First guess of minimum i.
             part%recti(2) = j                ! First guess of minimum j.
             part%recti(3) = i                ! First guess of maximum i.
             part%recti(4) = j                ! First guess of maximum j.
          else                                  ! Adjust for minimum and maximum i/j.
            if(part%recti(1) > i)part%recti(1) = i
            if(part%recti(2) > j)part%recti(2) = j
            if(part%recti(3) < i)part%recti(3) = i
            if(part%recti(4) < j)part%recti(4) = j
          endif
          part%snpc = part%snpc + 1         ! Accumulation for npcs.
        endif
      enddo
    enddo
#endif
    part%recto(1) = part%recti(1) - gsi%halosize
    part%recto(2) = part%recti(2) - gsi%halosize
    part%recto(3) = part%recti(3) + gsi%halosize
    part%recto(4) = part%recti(4) + gsi%halosize
    part%nx  = part%recto(3) - part%recto(1) + 1
    part%ny  = part%recto(4) - part%recto(2) + 1
    part%nxy = part%nx*part%ny
    DBGB0('gather npc')
    !print'("AAAAAAAAA", 2i3.2,4i8)',gsi%pid,gsi%npe,part%snpc,gsi%sumnp
    call MPI_GATHER(part%snpc,1,mpi_integer,gsi%npcs ,1,mpi_integer,0,gsi%mpi_comm,ierr)
    DBGB0('gather rect')
    call MPI_GATHER(part%recti,4,mpi_integer,gsi%rects,4,mpi_integer,0,gsi%mpi_comm,ierr)
    DBGB0('BCAST npcs')
    call MPI_BCAST(gsi%npcs,gsi%npe,MPI_INTEGER4,0,gsi%mpi_comm,ierr)
    DBGB0('BCAST rects')
    call MPI_BCAST(gsi%rects,4*gsi%npe,MPI_INTEGER4,0,gsi%mpi_comm,ierr)

    !---------------------------------------------------------------------------------------------
    ! Check balance and points omitted.
    ! Only the first PE do the following checks.
    !---------------------------------------------------------------------------------------------
    if(gsi%pid == 0)then
      if(gsi%balance == 1)then               ! This check is only for irregular partition.
        do n = 1, gsi%npe                   ! Check absolute balance for each PE.
          if(abs(gsi%npcs(n) - gsi%avenp) > 1)then
            write(6, *)'Error B Not Absolute balance ', n,gsi%npe, gsi%npcs(n), gsi%avenp, "nwpcs",gsi%npcs,"sum:",sum(gsi%npcs),gsi%sumnp
            !call irrp_abort(__FILE__, __LINE__)
          endif
        enddo
      endif
    endif
  end subroutine SetPartInf
#endif

!-------------------------------------------------------------------------------------------------
  subroutine set_part_ij2s(m,part,scycle)
    integer, intent(in) :: m
    type(pi_partinfo_type), intent(inout) :: part
    integer, intent(in),optional:: scycle(:, :) ! [4,N]
    integer :: i, j, i0,j0,i1,j1,nn,NS3,n
    !   need to adjust for halosize.
    NS3=size(scycle, 2)
    part%ib=part%recto(1);part%ibi=part%ib
    part%jb=part%recto(2);part%jbi=part%jb
    part%ie=part%recto(3);part%iei=part%ie
    part%je=part%recto(4);part%jei=part%je

    if(part%ibi<1     )part%ibi=1
    if(part%iei>gsi%im)part%iei=gsi%im
    if(part%jbi<1     )part%jbi=1
    if(part%jei>gsi%jm)part%jei=gsi%jm
    allocate(part%ij2s( part%ib:part%ie,part%jb:part%je) )
    part%ij2s = 0;
    nn=0;
    do i = part%ibi, part%iei                        ! sequentializing points from matrix form.
      do j = part%jbi, part%jei
        if(gsi%mask(i,j)>0)part%ij2s(i, j)=ij2gs(i,j)
      enddo
    enddo
    if(gsi%cycle_flag == 1)then
      if(part%jb<1     )part%jb=1
      if(part%je>gsi%jm)part%jb=gsi%jm
      if(part%ib==0)then
        do j = part%jbi, part%jei
          if(gsi%mask(gsi%im,j)>0) part%ij2s(0, j)=ij2gs( gsi%im,j)
        enddo
      endif
      if(part%ie==gsi%im+1)then
        do j = part%jbi, part%jei
          if(gsi%mask(1,j)>0)part%ij2s(gsi%im+1, j)=ij2gs(1,j)
        enddo
      endif
    elseif(gsi%cycle_flag == 2)then
      if(part%ib<1     )part%ib=1
      if(part%ie>gsi%im)part%ib=gsi%im
      if(part%jb<1     )part%jb=1
      if(part%je>gsi%jm)part%jb=gsi%jm
      if(part%ib==0)then
        do j = part%jbi, part%jei
          if(gsi%mask(gsi%im,j)>0) part%ij2s(0, j)=ij2gs( gsi%im,j)
        enddo
      endif
      if(part%ie==gsi%im+1)then
        do j = part%jbi, part%jei
          if(gsi%mask(1,j)>0)part%ij2s(gsi%im+1, j)=ij2gs(1,j)
        enddo
      endif

      if(part%jb==0)then
        do i = part%ibi, part%iei
          if(gsi%mask(i,gsi%jm)>0)part%ij2s(i,0)=ij2gs(i,gsi%jm)
        enddo
      endif
      if(part%je==gsi%jm+1)then
        do i = part%ibi, part%iei
          if(gsi%mask(i,1)>0)part%ij2s(i, gsi%jm+1)=ij2gs(i,1)
        enddo
      endif

    elseif(gsi%cycle_flag == 3)then
      if(.not. present(scycle))then
        write(6,*) "The scycle is not given when cycle_flag == 3. "
        call irrp_abort(__FILE__, __LINE__)
      endif
      do n=1,NS3
        i1=scycle(1,n);j1=scycle(2,n)
        if(i1>=part%ib.and.i1<=part%ie.and.j1>=part%jb.and.j1<=part%je)then
          i0=scycle(2,n);j0 =scycle(3,n)
          if(gsi%mask(i0,j0)>0)part%ij2s(i1, j1) = ij2gs(i0,j0)
        endif
      enddo
    elseif(gsi%cycle_flag == 0)then
    else
      write(6,*) "Wrong cycle_flag, the given value is ", gsi%cycle_flag
      call irrp_abort(__FILE__, __LINE__)
    endif
  end subroutine set_part_ij2s

  !-----------------------------------------------------------------------------------------------
  ! sub. set_recv_rplist: prepare the list of the receiving points for part.
  !
  !    set_recv_rplist(m, part, idx)
  !
  !   * type(pi_partinfo_type) :: part = the information of the given partition.
  !   * integer                :: m    = ID of current PE, this number is begin from 1.
  !   # integer                :: idx  = the number of receiving points.
  !-----------------------------------------------------------------------------------------------

  subroutine set_recv_rplist(m, part, idx)
    integer, intent(in)                :: m
    integer, intent(out)               :: idx
    type(pi_partinfo_type), intent(inout) :: part
    integer :: i, j, is, i1, j1, ih
    real :: dm, dmm
    real, allocatable :: mt(:, :)
    real, allocatable :: mvs(:)
    allocate(mt(part%recto(1):part%recto(3), part%recto(2):part%recto(4)))
    allocate(mvs(gsi%halosize))
    mt = 0
    do j = part%recti(2), part%recti(4)
    do i = part%recti(1), part%recti(3)
      is = part%ij2s(i, j)
      if(is /= 0 .and. gsi%mask(i,j) == m)then
        dm = 1                              ! The added offset for surrounding points is 1 for
                                            ! the nearest outer points.
        mt(i, j) = mt(i, j) - 1e30          ! The point belong to m-th PE is set as nearly infinite.
        do ih = gsi%halosize, 1, -1         ! For each circle of surrounding points.
          do j1 = j-ih, j+ih                ! For the receiving points neighbored with inner points,
            do i1 = i-ih, i+ih              ! the added value is greater than mvs(1). For the second
              mt(i1, j1) = mt(i1, j1) + dm  ! nearest ones, the added value is greater than mvs(2).
            enddo                           ! And so on ......
          enddo
          dmm = 2 * (gsi%halosize + 1 - ih) + 1
          dmm = dmm * dmm + 1
          dm = dm * dmm                     ! Shift the offset for next circle of points.
          mvs(ih) = dm                      ! Record upper limits for each circle of outer points.
          dm = dm * (1 + 1. / dmm)          ! Add a very small value for the next circle.
        enddo
      endif
    enddo
    enddo
    idx = 0;
    !if(gsi%pid==6)print*,'6=>1ppp',part%recto,part%ij2s(361, 57),gs2i(part%ij2s(361, 57)),gs2j(part%ij2s(361, 57))
    !if(gsi%pid==8 )print*,'8=>?',part%recto,part%ij2s(0, 2434),mt(0, 2434),gsi%mask(10800, 2434),gsi%mask(1, 2434)
    do j = part%recto(2), part%recto(4)
    do i = part%recto(1), part%recto(3)
      is = part%ij2s(i, j)                   ! part%ij2s(i, j) == 0 means land point.
      if(is ==0)cycle
      call gs2ij(is,i1,j1)
      if(gsi%mask(i1,j1) == m)cycle
      if(gsi%mask(i1,j1) == 0)cycle            ! Z must not land
      if(is /= 0 .and. mt(i, j) > 0)then    ! The mt(i, j) > 0 means receiving point.
        idx = idx + 1;
        part%rplist(idx)%gis = is ! Record this receiving point.
        part%rplist(idx)%i   = gs2i(is)
        part%rplist(idx)%j   = gs2j(is)
        part%rplist(idx)%pe  = gsi%mask(gs2i(is),gs2j(is))
        !if(gsi%pid==509.and.part%rplist(idx)%pe==9)print*,'509=>8',i,j,gs2i(is),gs2j(is),part%recto
        !if(gsi%pid==8.and.part%rplist(idx)%pe==510)print*,'8=>509',i,j,gs2i(is),gs2j(is),part%recto

        ! Using the following method can find out which circle
        ! is belonged of this receiving point.
        !     do ih = 1, gsi%halosize
        !       if(mt(i, j) < mvs(ih))then
        !         ! ih buf
        !       endif
        !     enddo
      endif
    enddo
    enddo
    call sort_rplist(part,idx)                   ! yinxq: 2015-3-3 16:30:42
    deallocate(mt, mvs)
  end subroutine set_recv_rplist

!-------------------------------------------------------------------------------------------------

  subroutine psort_rplist(part,idx,seqs,nna)
    integer,intent(in) :: idx,nna
    type(pi_partinfo_type), intent(in) :: part
    integer,intent(in):: seqs(:)
    integer :: i,ii,j
    integer :: idxn, pidc
    integer, allocatable :: nbps(:), nnbps(:), rpind(:)
    character*256 fn
    type(pi_rplist_type), pointer :: rplist(:) ! (nxy)
    type(pi_rplist_type), pointer ::rp
    integer :: ix,jy
    integer :: pe,pc,pu,pd
    rplist=>part%rplist
    pc=gsi%pid+1
    pu=pc+1;pd=pc-1;
    if(part%recto(2)==1)pd=-1;
    if(part%recto(4)==gsi%jm)pu=-1
    allocate(nbps(gsi%npe), nnbps(gsi%npe), rpind(idx))
    idxn = 0; nbps = 0; nnbps = 0
    do ii = 1, idx
      pidc = rplist(ii)%pe
      do j = idxn, 1, -1
        if(pidc==nbps(j))exit
      enddo

      if(j < 1)then
        idxn = idxn + 1; j = idxn;
        nbps(idxn)=pidc;
      endif
      nnbps(j) = nnbps(j) + 1;
      rpind(ii) = j
    enddo
    write(fn,'("ssorti_",i2.2,".txt")')gsi%pid
    open(unit=1006,file=trim(fn),status='unknown')
    write(*,'("RECT:",4i4,":",2i4,i8,":",4i3,": ",16i5)') ,part%recto,gsi%im,gsi%jm,nna,gsi%pid,pc,pu,pd,nbps(1:gsi%npe)
    write(1006,'("RECT:",4i4,":",2i4,i8,":",4i3,": ",16i5)') ,part%recto,gsi%im,gsi%jm,nna,gsi%pid,pc,pu,pd,nbps(1:gsi%npe)
    do jy=gsi%jm-1,1,-1
      write(1006,'(1000i2)')(gsi%mask(ix,jy),ix=1,gsi%im)
    enddo
    do ii=1,idx
      rp=>rplist(ii)
      write(1006,'(i8,2i9,3i6)')ii,rp%gis,seqs(ii),rp%i,rp%j,rp%pe
    enddo
    close(1006)
  end subroutine 
  subroutine sort_rplist(part,idx)
    integer,intent(in) :: idx
    type(pi_partinfo_type), intent(inout) :: part
    integer :: i,ii,exchange
    integer,allocatable :: seqs(:)
    integer :: lf,nna,pe,pc,pu,pd,j,ixm,ixb,ixe
    type(pi_rplist_type), pointer :: rplist(:) ! (nxy)
    type(pi_rplist_type) :: rp ! (nxy)
    allocate(seqs(idx))
    rplist=>part%rplist
    ! Set a scale value for sorting.
    ! If i is same, place greater j first; if j is same, place smaller i first.
    ! Use gsi%jm*2-j as the j-index, and gsi%im*2 as a scale for j to distinct with i.
    !
    nna=(gsi%im+gsi%jm)*16
    nna=100000
    pc=gsi%pid+1
    pu=pc+1;pd=pc-1;
    if(part%recto(2)==1)pd=-1;
    if(part%recto(4)==gsi%jm)pu=-1
    ixb=part%recto(1)
    ixe=part%recto(3)
    ixm=(ixb+ixe)/2
    do ii=1,idx
      ! yinxq 2015-4-4 16:52:08 !ZZ
      ! zhaowei 2025-3-21
      i=rplist(ii)%i
      j=rplist(ii)%j
      pe=rplist(ii)%pe
      !assume max 8 holo
      ! nna 是足够大的数，使得四条边不会混淆

      ! 下边界排第一,从左到右，如果i相等，则从上到下
      ! down boundart 1st ,left to right,if i is same,up to down
      ! sample: 1 holo
      !+++++++++++++++++++++++++++++++++++++++
      !01234567890123456789+++++++++++++++++++
      !                   01234567890123456789
      ! sample: 2 holo
      !+++++++++++++++++++++++++++++++++++++++
      !OOOOOOOOOOOOOOOOOOOO+++++++++++++++++++
      !OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      !                  OOOOOOOOOOOOOOOOOOOOO

      if(pe==pd)then ;lf=nna*0+i*8-j
      ! down boundart last ,left to right,if i is same,up to down
      else if(pe==pu)then ;lf=nna*3+i*8-j
      ! right boundart 2nd ,down to up,if j is same,right to left
      ! fixme if width <=3 will be error
      else if(i<ixb.or.(i<=ixe.and.i>=ixm))then;lf=nna*1+j*8-i
      ! left boundart 3rd ,down to up,if j is same,right to left
      else               ;lf=nna*2+j*8-i
      endif
      seqs(ii)=lf
    enddo
    exchange=1
    do while(exchange/=0)
      exchange=0
      do ii=2,idx
        if(seqs(ii-1)>seqs(ii))then
          i=seqs(ii-1);   seqs(ii-1)=seqs(ii);    seqs(ii)=i
          rp=rplist(ii-1);rplist(ii-1)=rplist(ii);rplist(ii)=rp
          exchange=1
        endif
      enddo
    enddo
    !call psort_rplist(part,idx,seqs,nna)
    deallocate(seqs)
  end subroutine sort_rplist

!-------------------------------------------------------------------------------------------------

  subroutine get_nbp_recv(m, part)
    integer, intent(in)                   :: m
    type(pi_partinfo_type), intent(inout) :: part
    integer, allocatable :: nbps(:), nnbps(:), rpind(:)
    integer :: idx, idxn, ii, pidc, j, is

    allocate(part%rplist(part%nxy))
    call set_recv_rplist(m, part, idx)
    part%npr = idx
    part%snp  = part%snpc + part%npr

    allocate(nbps(gsi%npe), nnbps(gsi%npe), rpind(idx))
    idxn = 0; nbps = 0; nnbps = 0
    do ii = 1, idx
      pidc = part%rplist(ii)%pe
      do j = idxn, 1, -1
        if(pidc==nbps(j))exit
      enddo

      if(j < 1)then
        idxn = idxn + 1; j = idxn;
        nbps(idxn)=pidc;
      endif
      nnbps(j) = nnbps(j) + 1;
      rpind(ii) = j
    enddo
    !print*,gsi%pid,idxn,nbps
    part%nids = idxn
    allocate(part%rsinfo(part%nids))
    do ii = 1, part%nids
      part%rsinfo(ii)%id = nbps(ii) - 1
      part%rsinfo(ii)%rn = nnbps(ii)
      allocate(part%rsinfo(ii)%srpnts(part%rsinfo(ii)%rn))
    enddo
    nnbps=0
    do ii=1,idx
      j= rpind(ii); nnbps(j) = nnbps(j) + 1
      part%rsinfo(j)%srpnts(nnbps(j)) = part%rplist(ii)%gis

    enddo
    allocate(part%rsbufs(part%nids))
    do j=1,part%nids
      allocate(part%rsbufs(j)%sbuf(gsi%mnprs), part%rsbufs(j)%rbuf(gsi%mnprs))
    enddo
    deallocate(nbps, rpind,part%rplist)
  end subroutine get_nbp_recv

#ifndef SERIAL_TEST

!-------------------------------------------------------------------------------------------------

  subroutine exchange(m, part)
    integer, intent(in)                   :: m
    type(pi_partinfo_type), intent(inout) :: part
    integer i, ierr
    integer, allocatable :: requests(:), requestr(:), statusr(:, :)
    allocate(requestr(part%nids), requests(part%nids), statusr(MPI_STATUS_SIZE, part%nids))
    do i = 1, part%nids
      !print*,gsi%pid,"=>",part%rsinfo(i)%id,part%rsbufs(i)%sbuf(2)
      !if(gsi%pid==509.or.part%rsinfo(i)%id==509)print*,gsi%pid,"=>",part%rsinfo(i)%id,part%rsbufs(i)%sbuf(2)
      call MPI_ISEND(part%rsbufs(i)%sbuf, part%rsbufs(i)%sbuf(2), MPI_INTEGER, &
                     part%rsinfo(i)%id, 10000, gsi%mpi_comm, requests(i), ierr)
    enddo
    DBGB0('irrp exchange 1')
    do i = 1, part%nids
      !if(gsi%pid==509.or.part%rsinfo(i)%id==509)print*,gsi%pid,"<=",part%rsinfo(i)%id,gsi%mnprs
      call MPI_IRECV(part%rsbufs(i)%rbuf, gsi%mnprs, MPI_INTEGER, part%rsinfo(i)%id, 10000, &
                     gsi%mpi_comm, requestr(i), ierr)
    enddo
    DBGB0('irrp exchange 2')

    call MPI_WAITALL(part%nids, requestr, statusr, ierr)
    call MPI_WAITALL(part%nids, requests, statusr, ierr)
    DBGB0('irrp exchange 3')
    deallocate(requestr,requests,statusr)
  end subroutine exchange

!-------------------------------------------------------------------------------------------------

  subroutine Getnbrecv(m, part)
    integer, intent(in)                   :: m
    type(pi_partinfo_type), intent(inout) :: part
    integer :: i,k
    do i = 1, part%nids
      part%rsbufs(i)%sbuf = 0
      part%rsbufs(i)%sbuf(1) = m
      part%rsbufs(i)%sbuf(2) = part%rsinfo(i)%rn + 3
      part%rsbufs(i)%sbuf(3) = part%rsinfo(i)%rn
      do k = 1, part%rsinfo(i)%rn
        part%rsbufs(i)%sbuf(3+k) = part%rsinfo(i)%srpnts(k)
      enddo
    enddo

    call exchange(m, part)

  end subroutine Getnbrecv

#else

!-------------------------------------------------------------------------------------------------
  !! AAAAAAAAAAAAAAAAAAAAA
  subroutine Getnbrecv(m, part)
    integer, intent(in)                   :: m
    type(pi_partinfo_type), intent(inout) :: part
    integer :: i, j, k, m1, rs
    do i = 1, part%nids
      m1 = part%rsinfo(i)%id + 1
      part%rsbufs(i)%rbuf = 0
      part%rsbufs(i)%rbuf(1) = m1
      do j = 1, parts(m1)%nids
        if(parts(m1)%rsinfo(i)%id+1 == m)then
          part%rsbufs(i)%rbuf(1) = m1
          part%rsbufs(i)%rbuf(2) = parts(m1)%rsinfo(i)%rn + 3
          part%rsbufs(i)%rbuf(3) = parts(m1)%rsinfo(i)%rn
          do k = 1, parts(m1)%rsinfo(j)%rn
            part%rsbufs(i)%rbuf(3+k) = parts(m1)%rsinfo(j)%srpnts(k)
          enddo
          exit
        endif
      enddo
    enddo
  end subroutine Getnbrecv
#endif
!-------------------------------------------------------------------------------------------------
  subroutine SetPnts(m, part)
    integer, intent(in)                   :: m
    type(pi_partinfo_type), intent(inout) :: part
    integer, allocatable :: iplist(:)     ! (npr)
    integer, allocatable :: aplist(:)     ! (npr)
    integer, allocatable :: g2l(:,:)        ! (npr)
    type(nb_8pnts_def_type), pointer :: nb8
    integer, pointer :: ij2s(:, :)
    integer :: i, ii, i1, i2, m1, m2, j, j1, npidxm, k, is, ls
    integer :: iptype, maxiptype, ips
    !DBGB0('SetPnts1')
    call Getnbrecv(m, part)
    !DBGB0('SetPnts2')
    do i=1,part%nids
      part%rsinfo(i)%sn    = part%rsbufs(i)%rbuf(3)
      allocate(part%rsinfo(i)%sspnts(part%rsinfo(i)%sn))
      part%rsinfo(i)%sspnts = part%rsbufs(i)%rbuf(4:3+part%rsinfo(i)%sn)
    enddo
    !print*,'AAAY',gsi%pid,part%rsinfo(2)%sspnts(229),gs2i(31680),gs2j(31680),part%recti,gsi%mask(gs2i(31680),gs2j(31680))

    !DBGB0('SetPnts3')
    part%snp = part%snpc + part%npr
    allocate(aplist(part%snp))
    allocate(iplist(part%snpc))
    allocate(g2l(part%ib:part%ie,part%jb:part%je))
    npidxm=0;iplist=0;g2l=0
    !DBGB0('SetPnts4')
    do j = part%recti(2), part%recti(4)
    do i = part%recti(1), part%recti(3)
      is=part%ij2s(i,j)
      if(gsi%mask(i,j) == m)then
        npidxm = npidxm + 1;
        iplist(npidxm) = is;
        g2l(i,j) = npidxm
      endif
    enddo
    enddo
    !DBGB0('SetPnts5')
    ! Sort Send
    npidxm=0
    if(allocated(tempmask))then
      maxiptype = maxval(tempmask)
      do iptype = 1, maxiptype !  yinxq 2014-3-23 20:04:32
        ! Sort Send pnts
        do i = 1, part%nids
          do j = 1, part%rsinfo(i)%sn
            is = part%rsinfo(i)%sspnts(j);
            ls = g2l(i,j)
            if(iplist(ls)/=0)then
              ips = iplist(ls);  call gs2ij(ips,i1,j1)
              if(iptype == tempmask(i1,j1))then
                npidxm = npidxm + 1
                aplist(npidxm) = iplist(ls)
                iplist(ls) = 0
                ! part%rsinfo(i)%sspnts(j)=npidxm !??
              endif
            endif
          enddo
        enddo
        part%nps = npidxm
        ! For inner pnts
        do ls = 1, part%snpc
          if(iplist(ls) /= 0)then
            ips = iplist(ls); call gs2ij(ips,i1,j1)
            if(iptype == tempmask(i1,j1))then
              npidxm = npidxm + 1
              aplist(npidxm) = iplist(ls)
              iplist(ls) = 0
            endif
          endif
        enddo
      enddo
      deallocate(tempmask)
    else ! if(allocated(tempmask))then
      ! Sort Send pnts
      do i = 1, part%nids
        do j = 1, part%rsinfo(i)%sn
          is = part%rsinfo(i)%sspnts(j);
          ls = g2l(gs2ip(is,part),gs2j(is))
          if(iplist(ls)/=0)then
            npidxm = npidxm + 1
            aplist(npidxm) = iplist(ls)
            iplist(ls) = 0
            ! part%rsinfo(i)%sspnts(j)=npidxm !??
          endif
        enddo
      enddo
      !DBGB0('SetPnts6')
      part%nps = npidxm

      ! For inner pnts
      do ls = 1, part%snpc
        if(iplist(ls) /= 0)then
          npidxm = npidxm + 1
          aplist(npidxm) = iplist(ls)
          iplist(ls) = 0
        endif
      enddo
    endif ! if(allocated(tempmask))then
    ! Sort recv pnts
    !DBGB0('SetPnts7')
    do i = 1, part%nids
      !part%rsinfo(i)%re = npidxm + 1
      do j = 1, part%rsinfo(i)%rn
        npidxm = npidxm + 1
        aplist(npidxm) = part%rsinfo(i)%srpnts(j)
      enddo
      !part%rsinfo(i)%re = npidxm
    enddo
    !DBGB0('SetPnts8')
    allocate(part%nb8(0:part%snp))
    allocate(part%ppos(0:part%snp))
    part%ppos(0)=0
    part%nb8(0)=nb_8pnts_def_type(0,0,0,0,0,0,0,0)
    g2l = 0
    ij2s => part%ij2s
    !ZNOUSE allocate(part%ij2s(gsi%im, gsi%jm))
    !DBGB0('SetPnts9')

    do ls = 1, part%snp
      is = aplist(ls);
      call gs2ij(is,i,j)
      g2l(gs2ip(is,part),gs2j(is)) = ls;
      part%ppos(ls) = aplist(ls)
    enddo
    !DBGB0('SetPnts10')
    do ls = 1, part%snpc
      is = aplist(ls)     ; nb8 => part%nb8(ls)
      call gs2ijp(is,i1,j1,part)
      nb8%ul = g2l(i1-1, j1+1); nb8%u = g2l(i1, j1+1); nb8%ur = g2l(i1+1, j1+1)
      nb8% l = g2l(i1-1, j1  );                        nb8% r = g2l(i1+1, j1  )
      nb8%dl = g2l(i1-1, j1-1); nb8%d = g2l(i1, j1-1); nb8%dr = g2l(i1+1, j1-1)
    enddo
    !DBGB0('SetPnts11')
    do i = 1, part%nids
      do j = 1, part%rsinfo(i)%sn
        is = part%rsinfo(i)%sspnts(j)
        part%rsinfo(i)%sspnts(j) = g2l(gs2ip(is,part),gs2j(is))
      enddo
      do j=1,part%rsinfo(i)%rn
        is = part%rsinfo(i)%srpnts(j)
        part%rsinfo(i)%srpnts(j) = g2l(gs2ip(is,part),gs2j(is))
      enddo
    enddo
    !DBGB0('SetPnts12')
    do i = 1, part%nids
      nullify(part%rsinfo(i)%mrpnts,part%rsinfo(i)%mspnts)
    enddo
    !DBGB0('SetPnts13')
    deallocate(iplist, aplist, g2l)
    allocate(part%spnttype(part%snp))
    !DBGB0('SetPnts14')
    part%spnttype=0
    !do i=1,part%snpc
    !enddo
    !DBGB0('SetPnts15')
  end subroutine SetPnts
  subroutine gs2ij(igs,i,j)
    integer igs,i,j
    i=mod(igs-1,gsi%im)+1
    j=(igs-1)/gsi%im+1
  end subroutine gs2ij
  integer function gs2i(igs)
  integer igs
    gs2i=mod(igs-1,gsi%im)+1
  end function gs2i
  integer function gs2ip(igs,part)
  integer igs,i
    type(pi_partinfo_type), intent(inout) :: part
    i=mod(igs-1,gsi%im)+1
    if(i<part%ib)i=i+gsi%im
    if(i>part%ie)i=i-gsi%im
    gs2ip=i
  end function gs2ip
  subroutine gs2ijp(igs,i,j,part)
    integer igs,i,j
    type(pi_partinfo_type), intent(inout) :: part
    i=mod(igs-1,gsi%im)+1
    j=(igs-1)/gsi%im+1
    if(i<part%ib)i=i+gsi%im
    if(i>part%ie)i=i-gsi%im
  end subroutine gs2ijp
  integer function gs2j(igs)
  integer igs
    gs2j=(igs-1)/gsi%im+1
  end function gs2j
  integer function ij2gs(i,j)
  integer i,j
    ij2gs=i-1+(j-1)*gsi%im+1
  end function ij2gs

!-------------------------------------------------------------------------------------------------
  end module irrp_kernal_mod

!-------------------------------------------------------------------------------------------------
!#################################################################################################
