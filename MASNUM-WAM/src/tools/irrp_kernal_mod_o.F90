#define DBGB0(msg)		! call barr(gsi%pid,gsi%npe,gsi%mpi_comm);if(gsi%pid==0)  write(6,'(f8.3,"s ",a)') Difftimer(),msg;call flush(6)
#define MPI_GATHER   MYMPI_GATHER
#define MPI_GATHERV  MYMPI_GATHERV
#define MPI_SCATTER  MYMPI_SCATTER
#define MPI_SCATTERV  MYMPI_SCATTERV

#define DBGINF   write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__

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
  implicit none

!-------------------------------------------------------------------------------------------------

  public :: irrp_part_init, irrp_part_final
  public :: irrp_SetPartMatrix, irrp_SetPartSerial,irrp_getrects

  public :: pi_pos_type, pi_rsinfo_type, nb_8pnts_def_type
  public :: globle_info_type,pi_partinfo_type,pi_rsbuf_type
  public :: gsi, cpart

  private

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
    integer, pointer :: rects(:, :)         ! (4, npe)
                                            !-----------------------------------------------------
                                            ! Temp Array no needed when partion end.
    integer                    :: mnprs     ! Buf size for partion exchange Data
    integer,           pointer :: spemask(:)! (0:sumnp)
    integer,           pointer :: ij2s(:, :)! (1-halosize:im+halosize, 1-halosize:jm+halosize)
    type(pi_pos_type), pointer :: s2ij(:)   ! (0:sumnp)
    integer,           pointer :: rplist(:) ! (sumnp)
                                            !-----------------------------------------------------
  end type globle_info_type

!-------------------------------------------------------------------------------------------------

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
    integer, pointer :: ij2s(:, :)
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
        i = mod(ij-1 , gsi%im) + 1
        j =    (ij-1)/ gsi%im  + 1
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
          i = mod(ij-1  , gsi%im) + 1
          j =    (ij-1) / gsi%im  + 1
          rsi%mrpnts(n)=i-i1+1+(j-j1)*cpart%nx
        enddo
        do n=1,rsi%sn
          ij=cpart%ppos(rsi%sspnts(n))
          i=mod(ij-1 ,gsi%im)+1
          j=   (ij-1)/gsi%im +1
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
    DBGB0('irrp_SetPartSerial 1')
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
    DBGB0('irrp_SetPartSerial 2')
    if(present(npc )) npc  = cpart%snpc
    if(present(np  )) np   = cpart%snp
    if(present(gnpc)) gnpc = gsi%gnpc
    if(present(plist))then
      if(associated(plist))deallocate(plist)
      allocate(plist(0:cpart%np));
      do n = 1, cpart%np
        plist(n)%i = mod(cpart%ppos(n)-1 ,gsi%im)+1
        plist(n)%j =    (cpart%ppos(n)-1)/gsi%im +1
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
!	    nwps_=0
!	    do inb=1,cpart%nids
!	      nwps_=nwps_+cpart%rsinfo(inb)%sn
!	    enddo
	    allocate(itmp(gsi%im*gsi%jm));itmp=0 ! yinxq 2015-4-4 16:45:00
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
    integer, intent(inout)        :: mask(:, :)
    integer, intent(in), optional :: halosize, cycle_flag
    integer, intent(in), optional :: scycle(:, :)
                                            ! (2,(im+jm+2)*2)
    integer :: i, j, ij, nn, m
    integer, allocatable :: pemask(:, :)    ! (im, jm), same shape as mask.
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
    allocate(gsi%s2ij(0:gsi%im*gsi%jm))           ! allocate temporary arraries.
    allocate(gsi%ij2s( (1-gsi%halosize-1):(gsi%im+gsi%halosize+1), &
                       (1-gsi%halosize-1):(gsi%jm+gsi%halosize+1)) )
    gsi%ij2s = 0; nn = 0
    do i = 1, gsi%im                        ! sequentializing points from matrix form.
      do j = 1, gsi%jm
        if(mask(i, j) > 0)then
          nn = nn + 1; gsi%s2ij(nn)%i = i; gsi%s2ij(nn)%j = j; gsi%ij2s(i, j) = nn
        endif
      enddo
    enddo
    call set_cycle_bound(scycle)            ! Deal with cycle boundaries.
    gsi%sumnp = nn                          ! Set number of sequentialized points.
    gsi%avenp = gsi%sumnp / float(gsi%npe)  ! Compute averaged number of points for each PE.
    gsi%mnprs = (gsi%avenp + 1 + 4) + 4     ! ??
#ifdef SERIAL_TEST
    allocate(parts(gsi%npe))                ! allocate parts.
    do m = 1, gsi%npe
      call init_part(parts(m))              ! initialize parts for all partion by serial.
    enddo
#else
    call init_part(cpart)                   ! initialize cpart.
#endif
		DBGB0('irrp_set_pemask begin')
    if(partmode == 0)then                   ! Do partion using irrp_pemask_mod
      !allocate(tempmask(gsi%im,gsi%jm)); tempmask = mask
      do i = 1, gsi%im                        ! sequentializing points from matrix form.
        do j = 1, gsi%jm
          if(mask(i, j) > 0)mask(i, j) = 1  ! yinxq: 2014-3-23 19:43:49
        enddo
      enddo
      allocate(pemask(gsi%im, gsi%jm))      ! allocate for pemask.
      call irrp_set_pemask(mask, pemask,gsi%s2ij, gsi%im, gsi%jm, gsi%npe)
      gsi%balance = 1
      mask = pemask
      deallocate(pemask)                    ! release pemask.
    else
      gsi%balance = 0
    endif

    allocate(gsi%spemask(0:gsi%sumnp))      ! allocate for spemask (sequentialized pemask).
    do nn=1,gsi%sumnp                     ! Re-arrange pemask into sequentialized arrary.
      gsi%spemask(nn) = mask(gsi%s2ij(nn)%i, gsi%s2ij(nn)%j)
    enddo

		DBGB0('irrp_set_pemask end')
    allocate(gsi%rects(4,gsi%npe))          ! allocate rects
    allocate(gsi%npcs(gsi%npe))             ! allocate npcs
		DBGB0('set_rects_npcs')
    call set_rects_npcs                     ! Set rects and npcs.
    allocate(gsi%rplist(gsi%sumnp))
#ifdef SERIAL_TEST
    do m = 1, gsi%npe
      call SetPartInf(m, parts(m))          ! Set information for each partition.
    enddo
    do m = 1, gsi%npe
      call get_nbp_recv(m, parts(m))        ! Set recv points for each partition.
    enddo
    do m = 1, gsi%npe
      call SetPnts(m, parts(m))             ! Set all points for each partition
    enddo
#else
		DBGB0('SetPartInf')
    call SetPartInf(gsi%pid+1, cpart)        ! Set information for current partition.
		DBGB0('get_nbp_recv')
    call get_nbp_recv(gsi%pid+1, cpart)     ! Set recv points for current partition.
		DBGB0('SetPnts')
    call SetPnts(gsi%pid+1, cpart)          ! Set all points for current partition.
		DBGB0('SetPnts End ')
#endif
    if(associated(gsi%spemask))deallocate(gsi%spemask)
    if(associated(gsi%s2ij   ))deallocate(gsi%s2ij   )
    if(associated(gsi%rplist ))deallocate(gsi%rplist )
    gsi%gnpc = sum(gsi%npcs)                ! Number of all computer points for all PEs.
		DBGB0('irrp_SetPartSerial')
    call irrp_SetPartSerial                 ! Prepare point list in serial case as default.
                                            ! It can be changed into matrix case later.
		DBGB0('End irrp_part_init')
                                            
  end subroutine irrp_part_init

  subroutine set_cycle_bound(scycle)
    integer, intent(in), optional :: scycle(:, :) ! [2,(im+jm+2)*2]
    integer :: i, j, ij
    !   need to adjust for halosize.
    if(gsi%cycle_flag == 1)then
      do j = 1, gsi%jm
        gsi%ij2s(0, j) = gsi%ij2s(gsi%im, j); gsi%ij2s(gsi%im+1, j) = gsi%ij2s(1, j)
      enddo
    elseif(gsi%cycle_flag == 2)then
      do i = 1, gsi%im
        gsi%ij2s(i, j) = gsi%ij2s(gsi%im-i+1, j)
      enddo
      do j = 0, gsi%jm+1
        gsi%ij2s(0, j) = gsi%ij2s(gsi%im, j); gsi%ij2s(gsi%im+1, j) = gsi%ij2s(1, j)
      enddo
    elseif(gsi%cycle_flag == 3)then
      if(.not. present(scycle))then
        write(6,*) "The scycle is not given when cycle_flag == 3. "
        call irrp_abort(__FILE__, __LINE__)
      endif
      j = 0; ij = 0
      do i = 0, gsi%im+1
        ij = ij + 1; gsi%ij2s(i, j) = gsi%ij2s(scycle(ij, 1), scycle(ij, 2))
      enddo
      i = gsi%im + 1
      do j = 1, gsi%jm+1
        ij = ij + 1; gsi%ij2s(i, j) = gsi%ij2s(scycle(ij, 1), scycle(ij, 2))
      enddo
      j = gsi%jm + 1
      do i = gsi%im, 0, -1
        ij = ij + 1; gsi%ij2s(i, j) = gsi%ij2s(scycle(ij, 1), scycle(ij, 2))
      enddo
      i = 0
      do j = gsi%jm, 1, -1
        ij = ij + 1; gsi%ij2s(i, j) = gsi%ij2s(scycle(ij, 1), scycle(ij, 2))
      enddo
    elseif(gsi%cycle_flag == 0)then
      return
    else
      write(6,*) "Wrong cycle_flag, the given value is ", gsi%cycle_flag
      call irrp_abort(__FILE__, __LINE__)
    endif
  end subroutine set_cycle_bound

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
    if(associated(gsi%spemask))deallocate(gsi%spemask)
    if(associated(gsi%s2ij   ))deallocate(gsi%s2ij   )
    if(associated(gsi%rplist ))deallocate(gsi%rplist )
    if(associated(gsi%ij2s   ))deallocate(gsi%ij2s   )
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
      if(associated(part%ij2s    ))deallocate(part%ij2s    )
      if(associated(part%requests))deallocate(part%requests)
      if(associated(part%requestr))deallocate(part%requestr)
      if(associated(part%statuss ))deallocate(part%statuss )
      if(associated(part%statusr ))deallocate(part%statusr )
    endif
  end subroutine release_part

  !-----------------------------------------------------------------------------------------------
  ! sub. set_rects_npcs: Prepare rects and npcs for all PEs.
  !-----------------------------------------------------------------------------------------------

  subroutine set_rects_npcs
    integer :: n, m, i, j
    gsi%rects = -1                          ! Initiale value for rects.
    gsi%npcs = 0                            ! Initiale value for npcs.
    do n = 1, gsi%sumnp
      m = gsi%spemask(n)                    ! Get ID for this point.
      i = gsi%s2ij(n)%i; j = gsi%s2ij(n)%j  ! Get i and j for this point.
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
    !---------------------------------------------------------------------------------------------
    ! Check balance and points omitted.
    ! Only the first PE do the following checks.
    !---------------------------------------------------------------------------------------------
    if(gsi%pid == 0)then
      if(gsi%balance == 1)then               ! This check is only for irregular partition.
        do n = 1, gsi%npe                   ! Check absolute balance for each PE.
          if(abs(gsi%npcs(n) - gsi%avenp) > 1)then
            if(gsi%pid == 0)then
              write(6, *)'Error Not Absolute balance ', n, gsi%npcs(n), gsi%avenp, gsi%pid
            endif
            call irrp_abort(__FILE__, __LINE__)
          endif
        enddo
      endif
      do n = 1, gsi%sumnp                   ! Check the omitted points.
        if(gsi%spemask(n) <= 0)then
          if(gsi%pid == 0)write(6, '(a,i9,2i5,i6,i9,f10.1)')'some Points Not Processed ', n, &
                          gsi%s2ij(n), gsi%spemask(n), gsi%sumnp, gsi%avenp
          call irrp_abort(__FILE__, __LINE__)
        endif
      enddo
    endif
  end subroutine set_rects_npcs

!-------------------------------------------------------------------------------------------------

  subroutine SetPartInf(m, part)
    integer, intent(in) :: m
    type(pi_partinfo_type), intent(inout) :: part
    nullify(part%ppos, part%nb8, part%rsinfo)
    nullify(part%requests, part%requestr)
    nullify(part%statuss, part%statusr, part%rsbufs)
    part%snpc     = gsi%npcs(m)
    part%snp      = 0
    part%nps      = 0
    part%npr      = 0
    part%nids     = 0
    part%recti    = gsi%rects(:, gsi%pid + 1)
    part%recto(1) = part%recti(1) - gsi%halosize
    part%recto(2) = part%recti(2) - gsi%halosize
    part%recto(3) = part%recti(3) + gsi%halosize
    part%recto(4) = part%recti(4) + gsi%halosize
  end subroutine SetPartInf

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
    type(pi_partinfo_type), intent(in) :: part
    integer :: i, j, is, i1, j1, ih
    real :: dm, dmm
    real, allocatable :: mt(:, :)
    real, allocatable :: mvs(:)
    allocate(mt(part%recto(1):part%recto(3), part%recto(2):part%recto(4)))
    allocate(mvs(gsi%halosize))
    mt = 0
    do j = part%recti(2), part%recti(4)
    do i = part%recti(1), part%recti(3)
      is = gsi%ij2s(i, j)
      if(is /= 0 .and. gsi%spemask(is) == m)then
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
    idx = 0; gsi%rplist = 0
    do j = part%recto(2), part%recto(4)
    do i = part%recto(1), part%recto(3)
      is = gsi%ij2s(i, j)                   ! gsi%ij2s(i, j) == 0 means land point.
      if(gsi%spemask(is) == m)cycle
      if(is /= 0 .and. mt(i, j) > 0)then    ! The mt(i, j) > 0 means receiving point.
        idx = idx + 1; gsi%rplist(idx) = is ! Record this receiving point.
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
    call sort_rplist(idx)                   ! yinxq: 2015-3-3 16:30:42
    deallocate(mt, mvs)
  end subroutine set_recv_rplist

!-------------------------------------------------------------------------------------------------

  subroutine sort_rplist(idx)
    integer,intent(in) :: idx
    integer :: i,ii,exchange
    integer,allocatable :: seqs(:)
    allocate(seqs(idx))
    ! Set a scale value for sorting. 
    ! If i is same, place greater j first; if j is same, place smaller i first.
    ! Use gsi%jm*2-j as the j-index, and gsi%im*2 as a scale for j to distinct with i.
    ! 
    ! is=gsi%rplist(ii);id=gsi%spemask(is);ix=gsi%s2ij(is)%i;iy=gsi%s2ij(is)%j
    ! seqs(ii)=(gsi%jm*2-iy) * gsi%im*2 + ix
    !
    do ii=1,idx
!      seqs(ii)=(gsi%jm*2-gsi%s2ij(gsi%rplist(ii))%j)*gsi%im*2+gsi%s2ij(gsi%rplist(ii))%i
! yinxq 2015-4-4 16:52:08
      seqs(ii)=(gsi%im*2-gsi%s2ij(gsi%rplist(ii))%i)*gsi%jm*2 + gsi%s2ij(gsi%rplist(ii))%j
    enddo
    exchange=1
    do while(exchange/=0)
      exchange=0
      do ii=2,idx
        if(seqs(ii-1)>seqs(ii))then
          i=seqs(ii-1);            seqs(ii-1)=seqs(ii);            seqs(ii)=i
          i=gsi%rplist(ii-1);gsi%rplist(ii-1)=gsi%rplist(ii);gsi%rplist(ii)=i
          exchange=1
        endif
      enddo
    enddo
    deallocate(seqs)
  end subroutine sort_rplist

!-------------------------------------------------------------------------------------------------

  subroutine get_nbp_recv(m, part)
    integer, intent(in)                   :: m
    type(pi_partinfo_type), intent(inout) :: part
    integer, allocatable :: nbps(:), nnbps(:), rpind(:)
    integer :: idx, idxn, ii, pidc, j, is
    call set_recv_rplist(m, part, idx)
    part%npr = idx
    part%snp  = part%snpc + part%npr
    allocate(nbps(gsi%npe), nnbps(gsi%npe), rpind(idx))
    idxn = 0; nbps = 0; nnbps = 0
    do ii = 1, idx
      pidc = gsi%spemask(gsi%rplist(ii))
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
      part%rsinfo(j)%srpnts(nnbps(j)) = gsi%rplist(ii)
    enddo
    allocate(part%rsbufs(part%nids))
    do j=1,part%nids
      allocate(part%rsbufs(j)%sbuf(gsi%mnprs), part%rsbufs(j)%rbuf(gsi%mnprs))
    enddo
    deallocate(nbps, rpind)
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
      call MPI_ISEND(part%rsbufs(i)%sbuf, part%rsbufs(i)%sbuf(2), MPI_INTEGER, &
                     part%rsinfo(i)%id, 10000, gsi%mpi_comm, requests(i), ierr)
    enddo
		DBGB0('irrp exchange 1')
    do i = 1, part%nids
      call MPI_IRECV(part%rsbufs(i)%rbuf, gsi%mnprs, MPI_INTEGER, part%rsinfo(i)%id, 10000, &
                     gsi%mpi_comm, requestr(i), ierr)
    enddo
		DBGB0('irrp exchange 2')
    call MPI_WAITALL(part%nids, requestr, statusr, ierr)
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
    integer, allocatable :: g2l(:)        ! (npr)
    type(nb_8pnts_def_type), pointer :: nb8
    integer, pointer :: ij2s(:, :)
    integer :: i, ii, i1, i2, m1, m2, j, j1, npidxm, k, is, ls
    integer :: iptype, maxiptype, ips
		DBGB0('SetPnts1')
    call Getnbrecv(m, part)
		DBGB0('SetPnts2')
    do i=1,part%nids
      part%rsinfo(i)%sn    = part%rsbufs(i)%rbuf(3)
      allocate(part%rsinfo(i)%sspnts(part%rsinfo(i)%sn))
      part%rsinfo(i)%sspnts = part%rsbufs(i)%rbuf(4:3+part%rsinfo(i)%sn)
    enddo
		DBGB0('SetPnts3')
    part%snp = part%snpc + part%npr
    allocate(aplist(part%snp))
    allocate(iplist(part%snpc))
    allocate(g2l(0:gsi%sumnp))
    npidxm=0;iplist=0;g2l=0
		DBGB0('SetPnts4')
		do j = part%recti(2), part%recti(4)
    do i = part%recti(1), part%recti(3)
      is=gsi%ij2s(i,j)
      if(is /= 0 .and. gsi%spemask(is) == m)then
        npidxm = npidxm + 1;iplist(npidxm) = is;g2l(is) = npidxm
      endif
    enddo
    enddo
		DBGB0('SetPnts5')
    ! Sort Send
    npidxm=0
    if(allocated(tempmask))then
      maxiptype = maxval(tempmask)
      do iptype = 1, maxiptype !  yinxq 2014-3-23 20:04:32
        ! Sort Send pnts
        do i = 1, part%nids
          do j = 1, part%rsinfo(i)%sn
            is = part%rsinfo(i)%sspnts(j);
            ls = g2l(is)
            if(iplist(ls)/=0)then
              ips = iplist(ls); i1 = gsi%s2ij(ips)%i ; j1 = gsi%s2ij(ips)%j
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
            ips = iplist(ls); i1 = gsi%s2ij(ips)%i ; j1 = gsi%s2ij(ips)%j
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
	        ls = g2l(is)
	        if(iplist(ls)/=0)then
	          npidxm = npidxm + 1
	          aplist(npidxm) = iplist(ls)
	          iplist(ls) = 0
	          ! part%rsinfo(i)%sspnts(j)=npidxm !??
	        endif
	      enddo
	    enddo
			DBGB0('SetPnts6')
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
		DBGB0('SetPnts7')
    do i = 1, part%nids
      !part%rsinfo(i)%re = npidxm + 1
      do j = 1, part%rsinfo(i)%rn
        npidxm = npidxm + 1
        aplist(npidxm) = part%rsinfo(i)%srpnts(j)
      enddo
      !part%rsinfo(i)%re = npidxm
    enddo
		DBGB0('SetPnts8')
    allocate(part%nb8(0:part%snp))
    allocate(part%ppos(0:part%snp))
    part%ppos(0)=0
    part%nb8(0)=nb_8pnts_def_type(0,0,0,0,0,0,0,0)
    g2l = 0
    ij2s => gsi%ij2s
    allocate(part%ij2s(gsi%im, gsi%jm))
		DBGB0('SetPnts3')
    do ls = 1, part%snp
      is = aplist(ls); g2l(is) = ls;
      i=gsi%s2ij(is)%i;j=gsi%s2ij(is)%j;
      part%ppos(ls) = i+(j-1)*gsi%im
      part%ij2s(i,j)= ls
    enddo
		DBGB0('SetPnts9')
    do ls = 1, part%snp
      is = aplist(ls)     ; nb8 => part%nb8(ls)
      i1 = gsi%s2ij(is)%i ; j1 = gsi%s2ij(is)%j
      nb8%ul = g2l(ij2s(i1-1, j1+1)); nb8%u = g2l(ij2s(i1, j1+1)); nb8%ur = g2l(ij2s(i1+1, j1+1))
      nb8% l = g2l(ij2s(i1-1, j1  ));                              nb8% r = g2l(ij2s(i1+1, j1  ))
      nb8%dl = g2l(ij2s(i1-1, j1-1)); nb8%d = g2l(ij2s(i1, j1-1)); nb8%dr = g2l(ij2s(i1+1, j1-1))
    enddo
		DBGB0('SetPnts10')
    do i = 1, part%nids
      do j = 1, part%rsinfo(i)%sn
        is = part%rsinfo(i)%sspnts(j)
        part%rsinfo(i)%sspnts(j) = g2l(is)
      enddo
      do j=1,part%rsinfo(i)%rn
        is = part%rsinfo(i)%srpnts(j)
        part%rsinfo(i)%srpnts(j) = g2l(is)
      enddo
    enddo
		DBGB0('SetPnts11')
    do i = 1, part%nids
      nullify(part%rsinfo(i)%mrpnts,part%rsinfo(i)%mspnts)
    enddo
		DBGB0('SetPnts12')
    deallocate(iplist, aplist, g2l)
    allocate(part%spnttype(part%snp))
		DBGB0('SetPnts13')
    part%spnttype=0
    do i=1,part%snpc
    enddo
		DBGB0('SetPnts3')
  end subroutine SetPnts


!-------------------------------------------------------------------------------------------------

  end module irrp_kernal_mod

!-------------------------------------------------------------------------------------------------
!#################################################################################################
