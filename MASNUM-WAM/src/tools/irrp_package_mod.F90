!#################################################################################################
!-------------------------------------------------------------------------------------------------
#define DBGINF   print*,__FILE__,__LINE__,gsi%pid
#ifdef ALLLOG
#define DBGB0(msg) write(6,'(2i5,f8.3,"s ",a)') __LINE__,gsi%pid,Difftimer(),msg;call flush(6)
#else
#define DBGB0(msg) if(gsi%pid==0)  write(6,'(2i5,f8.3,"s ",a)') __LINE__,gsi%pid,Difftimer(),msg;call flush(6)
#endif
#define DBGPV(tag)    
!print'(a,i8)',tag,gsi%pid;call flush(6)                           
#define MPI_GATHER   MYMPI_GATHER
#define MPI_GATHERV  MYMPI_GATHERV
#define MPI_SCATTER  MYMPI_SCATTER
#define MPI_SCATTERV  MYMPI_SCATTERV

  module irrp_package_mod

!-------------------------------------------------------------------------------------------------

  use irrp_smpi_mod
  use irrp_kernal_mod

  implicit none

!-------------------------------------------------------------------------------------------------
  public :: pi_pos_type                     ! Self-defined type to record 2 dimensional index.
  public :: nb_8pnts_def_type               ! Self-defined type for 8 neighbored pnts.
                                            ! pi_pos_type & nb_8pnts_def_type are from
                                            !       irrp_kernal_mod
	public :: irrp_output_pposg
  public :: irrp_init                       ! Initialize the irrp package.
  public :: irrp_init_data
  public :: irrp_exginf
  public :: irrp_SetPartMatrix              ! Set the pointers arranging method as Matrix way.
  public :: irrp_SetPartSerial              ! Set the pointers arranging method as Serial way.
  public :: irrp_getrects
  public :: irrp_final                      ! Finalize/Deinit the irrp package.

  public :: irrp_gather                     ! Gathering info/data from all PEs to root PE.
  public :: irrp_scatter                    ! Scatering info/data from root PE to all PEs.

  public :: irrp_exg_init                   ! Initialize exchage functions
  public :: irrp_exg_setvar                 ! Set/append var into a special exchange group.
  public :: irrp_exg_action                 ! Act the exchange by group.
  public :: irrp_exg_final                  ! Finalize the exchange functions.
  public :: irrp_exg_start                  ! start the exchange by group,equal irrp_exg_action(1).
  public :: irrp_exg_check                  ! test mpi irecv  is end ,
                                            ! in some realizations of MPI standard ,big data packet
                                            ! completed at wait or test
  public :: irrp_exg_end                    ! end the exchange by group,equal irrp_exg_action(6).

  public :: irrp_scatter_force_init         ! Initialize the function of forcing scatter
  public :: irrp_scatter_force              ! Act the forcing scattering.
  public :: irrp_scatter_force_final        ! Finalize the scatter_force Function.

  public :: irrp_scatter_ext_init           ! Initialize scatter with extended pnts.
  public :: irrp_scatter_ext                ! Scattering the info with extended pnts.
  public :: irrp_scatter_ext_final          ! Finalize scatter_ext Function

  private

!-------------------------------------------------------------------------------------------------

#ifndef typ

!-------------------------------------------------------------------------------------------------

  type globle_buf_type                      ! --- Variables related with global buffer.
    integer             :: gnpc             ! Number of global calculation points.
    integer,    pointer :: gdispls(:)       ! Displacements of each partition (npe).
                                            ! Being used together with gsi%npcs in functions of
                                            ! MPI_GATHERV & MPI_SCATTERV (npe).
    integer,    pointer :: pposg(:)         ! Relationship between global storage and partition
                                            ! storage (gnpc).
                                            ! The follows are 4 kinds of ponter variables.
                                            ! gsv** variables continuely stored in global. (gnpc)
                                            ! lv**  local variables in current partition. (snpc)
    integer(2), pointer :: gsvi2(:), lvi2(:)! In type of integer(kind=2).
    integer(4), pointer :: gsvi4(:), lvi4(:)! In type of integer(kind=4).
    real(4)   , pointer :: gsvr4(:), lvr4(:)! In type of real(kind=4).
    real(8)   , pointer :: gsvr8(:), lvr8(:)! In type of real(kind=8).
  end type globle_buf_type

!-------------------------------------------------------------------------------------------------

  type vars_def_type                        ! Variables used for exchange.
    character(len=16)   :: vname            ! Name of variable.
    integer             :: vartype = 0      ! Type of variable: 1 = r4, 2 = r8, 3 = i2, 4 = i4
    integer             :: km = 0           ! Number of Vertical layers.
                                            ! This should be in front, i.e. (km, im, jm)
                                            ! If it is not the case, it can be done one by one.
                                            ! The follows are 4 kinds of ponter variables.
                                            ! v**: the variable need to be exchanged.
                                            ! They need to be setted by irrp_exg_setvar, and use
                                            ! nullify to clean up.
    integer(2), pointer :: vi2(:, :)        ! Type of integer(kind=2), dims of (km,snpc).
    integer(4), pointer :: vi4(:, :)        ! Type of integer(kind=4), dims of (km,snpc).
    real(4)   , pointer :: vr4(:, :)        ! Type of real(kind=4), dims of (km,snpc).
    real(8)   , pointer :: vr8(:, :)        ! Type of real(kind=8), dims of (km,snpc).
  end type vars_def_type

!-------------------------------------------------------------------------------------------------

  type exchange_group_type
    integer                      :: nv = 0  ! Number of variables.
    integer                      :: mv = 0  ! Max-length of exchange buffer.
    integer                      :: state = 0
                                            ! 0:for not initting . 0 means not.
                                            ! 1:inited
                                            ! 2:sending
                                            ! 3:recving & sending
                                            ! 4:check end
                                            ! 5:wait end
                                            ! 6:exchange end
    type(vars_def_type), pointer :: vardef(:)
                                            ! Variables defining. (nv)
    type(mpipacket), pointer     :: sbuf(:) ! Buffer for sending info/data. (mv)
    type(mpipacket), pointer     :: rbuf(:) ! Buffer for receiving info/data. (mv)
  end type exchange_group_type

!-------------------------------------------------------------------------------------------------

  type force_info_type                      ! needed for scatter_ext, scatter_force
    integer            :: nx, ny, nxy       ! Number of input dim-size of forcing data.
    integer            :: nfp   = 0         ! Number points of input forcing data
    integer            :: infp  = 0         ! input pnts for cpart
    integer            :: infpg = 0         ! input pnts for all part infpg = sum(infp)
    integer,   pointer :: fpos(:)           ! ????(infp)
    integer,   pointer :: fposg(:)          ! Start position of each part (infpg), only master id
    integer,   pointer :: infps(:)          ! Input points of each part (npe).
    integer,   pointer :: fdispls(:)        ! Displacements of each PEs (npe).
                                            ! The follows are 4 kinds of ponter variables.
                                            ! gsv** variables continuely stored in global. (gnpc)
                                            ! lv**  local variables in current partition. (snpc)
    integer(2),pointer :: gsvi2(:), lvi2(:) ! In type of integer(kind=2).
    integer(4),pointer :: gsvi4(:), lvi4(:) ! In type of integer(kind=4).
    real(4)   ,pointer :: gsvr4(:), lvr4(:) ! In type of real(kind=4).
    real(8)   ,pointer :: gsvr8(:), lvr8(:) ! In type of real(kind=8).
    integer   ,pointer :: idx(:, :)         ! Indeces of each related points (4 * nfp) for interp.
    real(8)   ,pointer ::   w(:, :)         ! Weighting of each related points (4 * nfp)
  end type force_info_type

!-------------------------------------------------------------------------------------------------

  ! !!!!!!!!!!!!!!!!?????????????
  type obc_info_type
    integer            :: nob, gnob
    integer,   pointer :: nobs(:)           ! (npe)
    integer,   pointer :: obdispls(:)       ! (npe)
    integer(2),pointer :: gsvi2(:)          ! (gnob)
    integer(4),pointer :: gsvi4(:)
    real(4)   ,pointer :: gsvr4(:)
    real(8)   ,pointer :: gsvr8(:)
  end type obc_info_type

#endif
!-------------------------------------------------------------------------------------------------

  type(globle_buf_type) :: gsd

  integer :: maxgroups = 0
  type(exchange_group_type), pointer :: exgroups(:)

  integer :: maxforces = 0
  type(force_info_type), pointer :: forces(:)
  type(force_info_type), target  :: scatterext

!-------------------------------------------------------------------------------------------------

  interface irrp_exg_setvar; module procedure                               &
     eg_setvar_ks_i2,  eg_setvar_ks_i4,  eg_setvar_ks_r4,  eg_setvar_ks_r8, &
      eg_setvar_s_i2,   eg_setvar_s_i4,   eg_setvar_s_r4,   eg_setvar_s_r8, &
      eg_setvar_m_i2,   eg_setvar_m_i4,   eg_setvar_m_r4,   eg_setvar_m_r8, &
    eg_setvar_kls_i2, eg_setvar_kls_i4, eg_setvar_kls_r4, eg_setvar_kls_r8, &
    eg_setvar_klm_i2, eg_setvar_klm_i4, eg_setvar_klm_r4, eg_setvar_klm_r8, &
     eg_setvar_km_i2,  eg_setvar_km_i4,  eg_setvar_km_r4,  eg_setvar_km_r8
  end interface

  interface irrp_gather          ; module procedure             &
    gather_klss_i2, gather_klss_i4, gather_klss_r4, gather_klss_r8, &
    gather_klms_i2, gather_klms_i4, gather_klms_r4, gather_klms_r8, &
    gather_klmm_i2, gather_klmm_i4, gather_klmm_r4, gather_klmm_r8, &
     gather_kss_i2,  gather_kss_i4,  gather_kss_r4,  gather_kss_r8, &
     gather_kms_i2,  gather_kms_i4,  gather_kms_r4,  gather_kms_r8, &
     gather_kmm_i2,  gather_kmm_i4,  gather_kmm_r4,  gather_kmm_r8, &
      gather_ss_i2,   gather_ss_i4,   gather_ss_r4,   gather_ss_r8, &
      gather_ms_i2,   gather_ms_i4,   gather_ms_r4,   gather_ms_r8, &
      gather_mm_i2,   gather_mm_i4,   gather_mm_r4,   gather_mm_r8
  end interface

  interface irrp_scatter         ; module procedure                 &
    scatter_klss_i2, scatter_klss_i4, scatter_klss_r4, scatter_klss_r8, &
    scatter_klms_i2, scatter_klms_i4, scatter_klms_r4, scatter_klms_r8, &
    scatter_klmm_i2, scatter_klmm_i4, scatter_klmm_r4, scatter_klmm_r8, &
     scatter_kss_i2,  scatter_kss_i4,  scatter_kss_r4,  scatter_kss_r8, &
     scatter_kms_i2,  scatter_kms_i4,  scatter_kms_r4,  scatter_kms_r8, &
     scatter_kmm_i2,  scatter_kmm_i4,  scatter_kmm_r4,  scatter_kmm_r8, &
      scatter_ss_i2,   scatter_ss_i4,   scatter_ss_r4,   scatter_ss_r8, &
      scatter_ms_i2,   scatter_ms_i4,   scatter_ms_r4,   scatter_ms_r8, &
      scatter_mm_i2,   scatter_mm_i4,   scatter_mm_r4,   scatter_mm_r8

  end interface irrp_scatter

  interface irrp_scatter_ext     ; module procedure                             &
    scatter_ext_klss_i2, scatter_ext_klss_i4, scatter_ext_klss_r4, scatter_ext_klss_r8, &
    scatter_ext_klms_i2, scatter_ext_klms_i4, scatter_ext_klms_r4, scatter_ext_klms_r8, &
    scatter_ext_klmm_i2, scatter_ext_klmm_i4, scatter_ext_klmm_r4, scatter_ext_klmm_r8, &
     scatter_ext_kss_i2,  scatter_ext_kss_i4,  scatter_ext_kss_r4,  scatter_ext_kss_r8 , &
     scatter_ext_kms_i2,  scatter_ext_kms_i4,  scatter_ext_kms_r4,  scatter_ext_kms_r8, &
     scatter_ext_kmm_i2,  scatter_ext_kmm_i4,  scatter_ext_kmm_r4,  scatter_ext_kmm_r8, &
      scatter_ext_ss_i2,   scatter_ext_ss_i4,   scatter_ext_ss_r4,   scatter_ext_ss_r8, &
      scatter_ext_ms_i2,   scatter_ext_ms_i4,   scatter_ext_ms_r4,   scatter_ext_ms_r8, &
      scatter_ext_mm_i2,   scatter_ext_mm_i4,   scatter_ext_mm_r4,   scatter_ext_mm_r8
  end interface

  interface irrp_scatter_force   ; module procedure      &
    scatter_force_klss_r4, scatter_force_klss_r8, &
    scatter_force_klms_r4, scatter_force_klms_r8, &
    scatter_force_klmm_r4, scatter_force_klmm_r8, &
     scatter_force_kss_r4,  scatter_force_kss_r8, &
     scatter_force_kms_r4,  scatter_force_kms_r8, &
     scatter_force_kmm_r4,  scatter_force_kmm_r8, &
      scatter_force_ss_r4,   scatter_force_ss_r8, &
      scatter_force_ms_r4,   scatter_force_ms_r8, &
      scatter_force_mm_r4,   scatter_force_mm_r8
  end interface

!-------------------------------------------------------------------------------------------------

  contains

!-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  ! init Part:
  ! input:
  !  partmode: ??????????????!!!!!!!!!!!!!!!
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

  subroutine irrp_init(partmode, npe, pid, mpi_comm, mask, halosize, cycle_flag, scycle)
    integer, intent(in) :: partmode
    integer, intent(in) :: npe
    integer, intent(in) :: pid
    integer, intent(in) :: mpi_comm
    integer, intent(inout) :: mask(:, :) ! (im, jm)
    integer, intent(in), optional :: halosize
    integer, intent(in), optional :: cycle_flag
    integer, intent(in), optional :: scycle(:, :) ! (2,(im+jm+2)*2)
    call irrp_part_init(partmode,   &
                        npe,        &
                        pid,        &
                        mpi_comm,   &
                        mask,       &
                        halosize,   &
                        cycle_flag, &
                        scycle      )
!    DBGB0('Initgsd')
!    call initgsd
!    DBGB0('Initgsd end')
  end subroutine irrp_init
	subroutine irrp_init_data
    call initgsd
	end subroutine irrp_init_data
!-------------------------------------------------------------------------------------------------

  subroutine irrp_final(mode)
    integer, optional :: mode
    call irrp_part_final(mode)
    if(present(mode))then
      if(mode>0)then
        call irrp_exg_final_all
        call irrp_scatter_force_final_all
        call gsd_final
      endif
    endif
  end subroutine irrp_final

!-------------------------------------------------------------------------------------------------

  subroutine gsd_final
    if(gsi%npe==0)return
    gsi%npe = 0
    if(associated(gsd%gdispls))deallocate(gsd%gdispls)
    if(associated(gsd%  pposg))deallocate(gsd% pposg)
    if(associated(gsd%  gsvi2))deallocate(gsd% gsvi2)
    if(associated(gsd%   lvi2))deallocate(gsd%  lvi2)
    if(associated(gsd%  gsvi4))deallocate(gsd% gsvi4)
    if(associated(gsd%   lvi4))deallocate(gsd%  lvi4)
    if(associated(gsd%  gsvr4))deallocate(gsd% gsvr4)
    if(associated(gsd%   lvr4))deallocate(gsd%  lvr4)
    if(associated(gsd%  gsvr8))deallocate(gsd% gsvr8)
    if(associated(gsd%   lvr8))deallocate(gsd%  lvr8)
  end subroutine gsd_final

!-------------------------------------------------------------------------------------------------
  ! Init gsd Data
  subroutine initgsd
    integer :: i,j,l, ierr
    if(associated(gsd%gdispls))return

    allocate(gsd%gdispls(gsi%npe))
    gsd%gdispls(1) = 0
    DBGB0('Initgsd 1')
    do i = 1, gsi%npe-1
      gsd%gdispls(i+1) = gsd%gdispls(i) + gsi%npcs(i)
    enddo
    DBGB0('Initgsd 2')

     allocate(gsd%pposg(0:gsi%gnpc)); gsd%pposg = 0
    DBGB0('Initgsd 3')
    call MPI_GATHERV(cpart%ppos(1), cpart%snpc, MPI_INTEGER4, gsd%pposg(1), &
                     gsi%npcs, gsd%gdispls, MPI_INTEGER4, 0, gsi%mpi_comm, ierr)
    DBGB0('Initgsd 4')
    !do l=1,gsi%gnpc
    !  call gs2ij(gsd%pposg(l),i,j)
    !  gsd%pposg(l)=i-1+(j-1)*gsi%im+1
    !enddo
    call MPI_BCAST(gsd%pposg, gsi%gnpc+1, MPI_INTEGER4, 0, gsi%mpi_comm, ierr)
    DBGB0('Initgsd 5')
    nullify(gsd%gsvi2, gsd%lvi2)
    nullify(gsd%gsvi4, gsd%lvi4)
    nullify(gsd%gsvr4, gsd%lvr4)
    nullify(gsd%gsvr8, gsd%lvr8)
  end subroutine initgsd
	subroutine irrp_output_pposg
		integer i
    if(gsi%pid==0)then
      open(11,file="pposg.dat",form='unformatted')
      WRITE(11)gsi%gnpc
      WRITE(11)gsd%pposg
      close(11)
    endif
	end subroutine irrp_output_pposg
!-------------------------------------------------------------------------------------------------

#ifndef NOEXCHANGE__
  ! Init Exchange
  ! group_ind:index of Exchange group,<=0 ,Find Index
  ! nvar_: Number Vars in the group,optional ,default 4
  subroutine irrp_exg_init(group_ind, nvar_) ! default nvar = 16
    ! Init exchange_group,can
    integer,intent(inout) ::group_ind
    integer,intent(in),optional ::nvar_
    integer :: mm, i, nvar
    type(exchange_group_type), pointer :: exgt(:)
    nvar=4
    if(present(nvar_))nvar=nvar_
    if(nvar<=0)nvar=4
    ! check group_ind
    if(group_ind<=0)then
      group_ind=maxgroups+1
      do i=1,maxgroups
        if(exgroups(i)%mv==0)then
          group_ind=i
          exit
        endif
      enddo
    endif
    if(group_ind>maxgroups)then
      mm=maxgroups;maxgroups=group_ind+4
      if(mm/=0)then
        exgt=>exgroups
        allocate(exgroups(maxgroups))
        exgroups(1:mm)=exgt
        deallocate(exgt)
      else
        allocate(exgroups(maxgroups))
      endif
      do i=mm+1,maxgroups
        call Init_exchange_group(exgroups(i),0)
      enddo
    endif
    ! Clear the group
    call irrp_exg_final(group_ind)
    ! Init the group
    call Init_exchange_group(exgroups(group_ind),nvar)
  end subroutine irrp_exg_init
  ! Clear all exchange group
  subroutine irrp_exg_final_all
    integer group_ind
    do group_ind=1,maxgroups
      call irrp_exg_final(group_ind)
    enddo
  end subroutine irrp_exg_final_all
  ! check var id of exchange group
  subroutine VerifyVarDef(exg,ivar)
    type(exchange_group_type):: exg
    integer ,intent(inout)::ivar
    type(vars_def_type), pointer :: vardef(:) ! nv
    integer mm,i
    if(ivar==0)then
      if(associated(exg%vardef))deallocate(exg%vardef)
      nullify(exg%vardef); exg%mv=0;exg%nv=0
      return
    endif
    vardef=>exg%vardef
    if(ivar<0)then
      ivar=exg%mv+1
      do i=1,exg%mv
        if(vardef(i)%vartype<=0)then
          ivar=i
          exit
        endif
      enddo
    endif
    if(exg%mv<ivar)then
      mm=exg%mv;if(mm<0)mm=0
      exg%mv=ivar
      if(mm>0)then
        vardef=>exg%vardef
        allocate(exg%vardef(exg%mv))
        exg%vardef(1:mm)=vardef
        deallocate(vardef)
      else
        allocate(exg%vardef(exg%mv))
      endif
      vardef=>exg%vardef
      do i=mm+1,exg%mv
        vardef(i)%vartype=0
        vardef(i)%vname='FREE'
        nullify(vardef(i)%vi2)
        nullify(vardef(i)%vi4)
        nullify(vardef(i)%vr4)
        nullify(vardef(i)%vr8)
      enddo
    endif
  end subroutine VerifyVarDef

  subroutine Init_exchange_group(exg,nv)
    type(exchange_group_type):: exg
    integer,intent(in)::nv
    integer :: i, iv
    exg%mv=0;exg%nv=0;exg%state=0
    nullify(exg%vardef)
    nullify(exg%sbuf)
    nullify(exg%rbuf)
    iv=nv;if(iv<0)iv=0
    call VerifyVarDef(exg,iv)
  end subroutine Init_exchange_group
  ! Clear exchange group
  subroutine irrp_exg_final(group_ind) ! default nvar = 16
    integer,intent(in):: group_ind
    type(exchange_group_type), pointer :: exg
    integer :: iv
    exg=>exgroups(group_ind)
    if(associated(exg%sbuf  ))deallocate(exg%sbuf)
    if(associated(exg%rbuf  ))deallocate(exg%rbuf)
    iv = 0
    call VerifyVarDef(exg, iv)
    call Init_exchange_group(exg,0)
  end subroutine irrp_exg_final
  ! Set var info into exchange group
  ! group_ind:exchange group index,<=0 auto find it
  ! ivar index of var in the group
  ! vname:var name
  ! ivt:var type
  ! km :km
  subroutine eg_setvar(group_ind, ivar, vname,ivt, km)
    integer,intent(inout):: group_ind, ivar
    integer,intent(in):: ivt
    integer,intent(in),optional:: km
    character*(*),intent(in):: vname
    type(vars_def_type), pointer :: vardef(:) ! nv
    type(exchange_group_type), pointer :: exg
    integer :: i,mm,mv,needInit
    needInit=0
    if(maxgroups==0 .or. (group_ind<=0.or. group_ind>maxgroups) )then
      call irrp_exg_init(group_ind,ivar)
    endif
    exg=>exgroups(group_ind)
    call VerifyVarDef(exg,ivar)
    vardef=>exg%vardef
    if(vardef(ivar)%vartype==0)exg%nv=exg%nv+1
    vardef(ivar)%vartype=ivt
    vardef(ivar)%vname=vname
    vardef(ivar)%km=km
    exg%state=0
  end subroutine eg_setvar

! eg_setvar_xx_nn
! Set var info into exchange group,4 var type,6 var shape
! all other shape reshape to (km,np)
! group_ind:exchange group index,<=0 auto find it
! ivar index of var in the group
! vname:var name
! var:var
! km :km
# define TMPL_eg_setvar_ks(NA,VT,MVT)                                          \
  subroutine eg_setvar_ks_##NA (group_ind, ivar, vname, var, km)              ;\
    integer,intent(inout):: group_ind, ivar                                   ;\
    integer,intent(in):: km                                                   ;\
    character*(*),intent(in):: vname                                          ;\
    VT,target :: var(km,cpart%iblv:cpart%np)                                  ;\
    call eg_setvar(group_ind, ivar, vname, MVT, km)                           ;\
    exgroups(group_ind)%vardef(ivar)%v##NA=>var                               ;\
  end subroutine eg_setvar_ks_##NA
# define _ ,
# define TMPL_eg_setvar(TAG,NA,VT,MVT,VSHAPE,EXARGS,EXARGS1)                   \
  subroutine eg_setvar_##TAG##NA (group_ind, ivar, vname, var EXARGS)      ;\
    integer,intent(inout):: group_ind, ivar EXARGS                            ;\
    character*(*),intent(in):: vname                                          ;\
    VT,TARGET :: var VSHAPE                                                   ;\
    call eg_setvar_ks_##NA(group_ind, ivar, vname, var EXARGS1)               ;\
  end subroutine eg_setvar_##TAG##NA

#define LVSA_S   (      cpart%iblv:cpart%np)
#define LVSA_M   (        cpart%nx,cpart%ny)
#define LVSA_KS  (km,   cpart%iblv:cpart%np)
#define LVSA_KM  (km,     cpart%nx,cpart%ny)
#define LVSA_KLS (km,lm,cpart%iblv:cpart%np)
#define LVSA_KLM (km,lm,  cpart%nx,cpart%ny)
  TMPL_eg_setvar_ks(r4,real   (4),MPI_REAL4   )
  TMPL_eg_setvar_ks(r8,real   (8),MPI_REAL8   )
  TMPL_eg_setvar_ks(i2,integer(2),MPI_INTEGER2)
  TMPL_eg_setvar_ks(i4,integer(4),MPI_INTEGER4)
  TMPL_eg_setvar(  s_, r4, real   (4), MPI_REAL4   ,LVSA_S,     ,_ 1)
  TMPL_eg_setvar(  s_, r8, real   (8), MPI_REAL8   ,LVSA_S,     ,_ 1)
  TMPL_eg_setvar(  s_, i2, integer(2), MPI_INTEGER2,LVSA_S,     ,_ 1)
  TMPL_eg_setvar(  s_, i4, integer(4), MPI_INTEGER4,LVSA_S,     ,_ 1)

  TMPL_eg_setvar(  m_, r4, real   (4), MPI_REAL4   ,LVSA_M,     ,_ 1)
  TMPL_eg_setvar(  m_, r8, real   (8), MPI_REAL8   ,LVSA_M,     ,_ 1)
  TMPL_eg_setvar(  m_, i2, integer(2), MPI_INTEGER2,LVSA_M,     ,_ 1)
  TMPL_eg_setvar(  m_, i4, integer(4), MPI_INTEGER4,LVSA_M,     ,_ 1)

  TMPL_eg_setvar( km_, r4, real   (4), MPI_REAL4   ,LVSA_KM, _ km,_ km)
  TMPL_eg_setvar( km_, r8, real   (8), MPI_REAL8   ,LVSA_KM, _ km,_ km)
  TMPL_eg_setvar( km_, i2, integer(2), MPI_INTEGER2,LVSA_KM, _ km,_ km)
  TMPL_eg_setvar( km_, i4, integer(4), MPI_INTEGER4,LVSA_KM, _ km,_ km)

  TMPL_eg_setvar(kls_, r4, real   (4), MPI_REAL4   ,LVSA_KLS,_ km _ lm,_ km*lm)
  TMPL_eg_setvar(kls_, r8, real   (8), MPI_REAL8   ,LVSA_KLS,_ km _ lm,_ km*lm)
  TMPL_eg_setvar(kls_, i2, integer(2), MPI_INTEGER2,LVSA_KLS,_ km _ lm,_ km*lm)
  TMPL_eg_setvar(kls_, i4, integer(4), MPI_INTEGER4,LVSA_KLS,_ km _ lm,_ km*lm)

  TMPL_eg_setvar(klm_, r4, real   (4), MPI_REAL4   ,LVSA_KLM,_ km _ lm,_ km*lm)
  TMPL_eg_setvar(klm_, r8, real   (8), MPI_REAL8   ,LVSA_KLM,_ km _ lm,_ km*lm)
  TMPL_eg_setvar(klm_, i2, integer(2), MPI_INTEGER2,LVSA_KLM,_ km _ lm,_ km*lm)
  TMPL_eg_setvar(klm_, i4, integer(4), MPI_INTEGER4,LVSA_KLM,_ km _ lm,_ km*lm)

  subroutine irrp_exg_enddef(group_ind)
    integer,intent(in)::group_ind
    type(pi_rsinfo_type), pointer :: rsinfo(:),rsi ! nids
    type(exchange_group_type), pointer :: exg
    type(vars_def_type), pointer :: var
    integer inb,lsize,iv
    rsinfo=>cpart%rsinfo
    exg=>exgroups(group_ind)
    if(.not. associated(exg%sbuf)) allocate(exg%sbuf(cpart%nids))
    if(.not. associated(exg%rbuf)) allocate(exg%rbuf(cpart%nids))
    do inb=1,cpart%nids
      rsi=>rsinfo(inb)
      lsize=0
      do iv=1,exg%mv
        var=>exg%vardef(iv)
        if(var%vartype==MPI_REAL4)then
          lsize=lsize+4*var%km
        else if(var%vartype==MPI_REAL8)then
          lsize=lsize+8*var%km
        else if(var%vartype==MPI_INTEGER2)then
          lsize=lsize+2*var%km
        else if(var%vartype==MPI_INTEGER4)then
          lsize=lsize+4*var%km
        endif
      enddo
      call InitMpiPacket(exg%sbuf(inb),lsize*rsi%sn)
      call InitMpiPacket(exg%rbuf(inb),lsize*rsi%rn)
    enddo
    exg%state=1
  end subroutine irrp_exg_enddef

  subroutine exchange_group_action_prepare(group_ind)
    integer,intent(in):: group_ind
    if(group_ind<0.or. group_ind>maxgroups )then
      write(6,*)'irrp_exg_action:group_ind is invalid',group_ind,maxgroups
      call irrp_abort(__FILE__, __LINE__)
    endif
    if(exgroups(group_ind)%state==0)then
      call irrp_exg_enddef(group_ind)
    endif
    if(cpart%Inited_RSV==0)then
      cpart%Inited_RSV=1
      allocate(cpart%requests(cpart%nids),cpart%statuss(MPI_STATUS_SIZE,cpart%nids));
!ZZ use 32 CPU at cpu=20   DBGT,cpart%nids,MPI_STATUS_SIZE
      allocate(cpart%requestr(cpart%nids),cpart%statusr(MPI_STATUS_SIZE,cpart%nids));
    endif
    cpart%statusr=0;cpart%statuss=0;cpart%requestr=0;cpart%requests=0;
  end subroutine exchange_group_action_prepare
  subroutine irrp_exginf(nst)
    integer nst
    integer inb,i,iac,ia,ic
    type(pi_rsinfo_type), pointer :: rsinfo(:), rsi ! nids
    rsinfo=>cpart%rsinfo
    do inb=1,cpart%nids
      rsi=>rsinfo(inb)
      write(6,'(a,i5,a,5i7)')'========',gsi%pid,'>',rsi%id,rsi%sn*nst
      !do i=1,rsi%sn
      !  iac=rsi%spnts(i)
      !  ia = mod(cpart%ppos(iac)-1 ,gsi%im)+1
      !  ic =    (cpart%ppos(iac)-1)/gsi%im +1
      !  write(50+gsi%pid,'(i5,a,5i5)')gsi%pid,'>',rsi%id,iac,ia,ic
      !enddo
      !do i=1,rsi%rn
      !  iac=rsi%rpnts(i)
      !  ia = mod(cpart%ppos(iac)-1 ,gsi%im)+1
      !  ic =    (cpart%ppos(iac)-1)/gsi%im +1
      !  write(50+gsi%pid,'(i5,a,5i5)')gsi%pid,'<',rsi%id,iac,ia,ic
      !enddo
    enddo
  end subroutine irrp_exginf
  subroutine exchange_group_action_send(group_ind)
    integer,intent(in):: group_ind
    integer inb,dstid,istart,isize,direct,tag
    integer :: ierr,esize
    type(pi_rsinfo_type), pointer :: rsinfo(:), rsi ! nids
    type(exchange_group_type), pointer :: exg
    type(vars_def_type), pointer :: var
    integer ::iv,i,km,lsize
    type(mpipacket), pointer :: pk
    integer, pointer :: pnts(:)
    exg=>exgroups(group_ind)
    exg%state=2
    rsinfo=>cpart%rsinfo

    do inb=1,cpart%nids
      rsi=>rsinfo(inb)
      dstid =rsi%id
      tag=1000;
      pk=>exg%sbuf(inb)
      pnts=>rsi%spnts
      call InitMpiPacket(pk)
      do iv=1,exg%mv
        var=>exg%vardef(iv)
        km=var%km
        if(var%vartype==MPI_REAL4)then
          do i=1,rsi%sn
            call MPI_PACK(var%vr4(1,pnts(i)),km,MPI_REAL4,pk%buf,pk%bsize,pk%pos,gsi%mpi_comm,ierr)
          enddo
        else if(var%vartype==MPI_REAL8)then
          do i=1,rsi%sn
            call MPI_PACK(var%vr8(1,pnts(i)),km,MPI_REAL8,pk%buf,pk%bsize,pk%pos,gsi%mpi_comm,ierr)
          enddo
        else if(var%vartype==MPI_INTEGER2)then
          do i=1,rsi%sn
            call MPI_PACK(var%vi2(1,pnts(i)),km,MPI_INTEGER2,pk%buf,pk%bsize,pk%pos,gsi%mpi_comm,ierr)
          enddo
        else if(var%vartype==MPI_INTEGER4)then
          do i=1,rsi%sn
            call MPI_PACK(var%vi4(1,pnts(i)),km,MPI_INTEGER4,pk%buf,pk%bsize,pk%pos,gsi%mpi_comm,ierr)
          enddo
        endif
      enddo
      call SetMpiPacketDSize(pk)
      ! pack Data
      !print'(a,i6,a,i6,i8)','Isend',gsi%pid,"=>",dstid,pk%dsize;call flush(6)
      call MPI_ISEND(pk%buf,pk%dsize,MPI_PACKED,dstid,tag,gsi%mpi_comm,cpart%requests(inb),ierr)
    end do
  end subroutine exchange_group_action_send

  subroutine exchange_group_action_recvs(group_ind)
    integer,intent(in):: group_ind
    integer :: inb,dstid,tag,ierr
    type(mpipacket), pointer ::  pk
    type(exchange_group_type), pointer :: exg
    integer     :: status(MPI_STATUS_SIZE)  ! mpi status info

    exg=>exgroups(group_ind)
    exg%state=3
    do inb=1,cpart%nids
      dstid =cpart%rsinfo(inb)%id
      tag=1000;
      pk=>exgroups(group_ind)%rbuf(inb)
      call InitMpiPacket(pk)
      !print'(a,i6,a,i6,i8)','recv ',gsi%pid,"<=",dstid,pk%bsize;call flush(6)
      call MPI_RECV(pk%buf,pk%bsize,MPI_PACKED,dstid,tag,gsi%mpi_comm,status,ierr)
    end do
    !call MPI_WAITALL(cpart%nids,cpart%requestr,cpart%statusr,ierr)
    call exchange_group_action_recv_unpack(group_ind)
    call MPI_WAITALL(cpart%nids,cpart%requests,cpart%statuss,ierr)
  end subroutine exchange_group_action_recvs
  subroutine exchange_group_action_recv(group_ind)
    integer,intent(in):: group_ind
    integer :: inb,dstid,tag,ierr
    type(mpipacket), pointer ::  pk
    type(exchange_group_type), pointer :: exg
    exg=>exgroups(group_ind)
    exg%state=3
    do inb=1,cpart%nids
      dstid =cpart%rsinfo(inb)%id
      tag=1000;
      pk=>exgroups(group_ind)%rbuf(inb)
      call InitMpiPacket(pk)
      !print'(a,i6,a,i6,i8)','irecv ',gsi%pid,"<=",dstid,pk%bsize;call flush(6)
      call MPI_IRECV(pk%buf,pk%bsize,MPI_PACKED,dstid,tag,gsi%mpi_comm,cpart%requestr(inb),ierr)
    end do
  end subroutine exchange_group_action_recv
  subroutine exchange_group_action_recv_unpack(group_ind)
    integer,intent(in):: group_ind
    integer :: inb,dstid,istart,isize,direct,tag
    integer :: ierr,esize
    type(pi_rsinfo_type), pointer :: rsi ! nids
    type(exchange_group_type), pointer :: exg
    type(vars_def_type), pointer :: var
    integer ::iv,i,km,lsize
    type(mpipacket), pointer ::  pk
    integer, pointer :: pnts(:)
    exg=>exgroups(group_ind)
    exg%state=5
    do inb=1,cpart%nids
      rsi=>cpart%rsinfo(inb)
      pk=>exg%rbuf(inb)
      !call MPI_GET_COUNT(cpart%statusr(1,inb),MPI_PACKED,lsize,ierr)
      call SetMpiPacketDSize(pk,pk%bsize)
      pnts=>rsi%rpnts
      ! unpack
      do iv=1,exg%mv
        var=>exg%vardef(iv)
        km=var%km
        if(var%vartype==MPI_REAL4)then
          do i=1,rsi%rn
            call MPI_UNPACK(pk%buf,pk%dsize,pk%pos,var%vr4(1,pnts(i)),km   ,MPI_REAL4  ,gsi%mpi_comm,ierr)
          enddo
        else if(var%vartype==MPI_REAL8)then
          do i=1,rsi%rn
            call MPI_UNPACK(pk%buf,pk%dsize,pk%pos,var%vr8(1,pnts(i)),km,MPI_REAL8   ,gsi%mpi_comm,ierr)
          enddo
        else if(var%vartype==MPI_INTEGER2)then
          do i=1,rsi%rn
            call MPI_UNPACK(pk%buf,pk%dsize,pk%pos,var%vi2(1,pnts(i)),km,MPI_INTEGER2,gsi%mpi_comm,ierr)
          enddo
        else if(var%vartype==MPI_INTEGER4)then
          do i=1,rsi%rn
            call MPI_UNPACK(pk%buf,pk%dsize,pk%pos,var%vi4(1,pnts(i)),km,MPI_INTEGER4,gsi%mpi_comm,ierr)
          enddo
        endif
      enddo
    enddo
  end subroutine exchange_group_action_recv_unpack
  subroutine irrp_exg_start(group_ind)
    integer,intent(in):: group_ind
    call exchange_group_action_prepare(group_ind)
    call exchange_group_action_send(group_ind)
    call exchange_group_action_recv(group_ind)
  end subroutine irrp_exg_start
  subroutine irrp_exg_check(group_ind,state)
    integer,intent(in):: group_ind
    integer,intent(out)::state
    logical flag
    integer ierr
    state=1
    if(exgroups(group_ind)%state==3)then
      call MPI_TESTALL(cpart%nids,cpart%requestr,flag,cpart%statusr,ierr)
      if(flag)then
        state=1
        exgroups(group_ind)%state=4
      else
        state=0
      endif
    endif
  end subroutine irrp_exg_check
  subroutine irrp_exg_end(group_ind)
    integer,intent(in):: group_ind
    integer ierr
    call MPI_WAITALL(cpart%nids,cpart%requestr,cpart%statusr,ierr)
    call MPI_WAITALL(cpart%nids,cpart%requests,cpart%statuss,ierr)
    call exchange_group_action_recv_unpack(group_ind)
  end subroutine irrp_exg_end
  subroutine irrp_exg_action(group_ind,mode_)
    integer,intent(in):: group_ind
    integer,intent(in),optional::mode_
    integer mode,ierr,n
    mode=7;
    if(present(mode_))mode=mode_
    if(mode==7)then
      call exchange_group_action_prepare(group_ind)
      call exchange_group_action_send(group_ind)
      call exchange_group_action_recvs(group_ind)
      call exchange_group_action_recv_unpack(group_ind)
      return
    endif
    if(iand(mode,1)/=0)then
      !if(gsi%pid==0)write(6,*)'sendrecv',mode
      call exchange_group_action_prepare(group_ind)
      call exchange_group_action_send(group_ind)
      call exchange_group_action_recv (group_ind)
    endif
    if(iand(mode,2)/=0)then
      !if(gsi%pid==0)write(6,*)'wait',mode
      do n=1,cpart%nids
       call MPI_WAIT(cpart%requestr(n),cpart%statusr(:,n),ierr)
       call MPI_WAIT(cpart%requests(n),cpart%statuss(:,n),ierr)
      enddo
      !call MPI_WAITALL(cpart%nids,cpart%requestr,cpart%statusr,ierr)
      !call MPI_WAITALL(cpart%nids,cpart%requests,cpart%statuss,ierr)
    endif
    if(iand(mode,4)/=0)then
      !if(gsi%pid==0)write(6,*)'unpack',mode
      call exchange_group_action_recv_unpack(group_ind)
    endif
  end subroutine irrp_exg_action

#endif

!-------------------------------------------------------------------------------------------------

#ifndef NOGATHER

  integer function Verifyroot(root_)
    integer ,optional::root_
    integer i,ierr
    Verifyroot=0;
    if(present(root_))then
      Verifyroot=root_
      if(Verifyroot<0 .or. Verifyroot>=gsi%npe)Verifyroot=0
    endif
  end function Verifyroot

  subroutine vg_lfvar(root,vt,flg)
    type(force_info_type)::fi
    integer,intent(in)::root,vt
    logical,intent(in)::flg
    logical isr
    isr = (gsi%pid==root)
    if(vt==MPI_REAL4)then
      if(.not.associated(gsd%gsvr4).and.isr)allocate(gsd%gsvr4(gsi%gnpc))
      if(.not.associated(gsd% lvr4).and.flg)allocate(gsd% lvr4(cpart%snpc ))
    elseif(vt==MPI_REAL8)then
      if(.not.associated(gsd%gsvr8).and.isr)allocate(gsd%gsvr8(gsi%gnpc))
      if(.not.associated(gsd% lvr8).and.flg)allocate(gsd% lvr8(cpart%snpc ))
    elseif(vt==MPI_INTEGER2)then
      if(.not.associated(gsd%gsvi2).and.isr)allocate(gsd%gsvi2(gsi%gnpc))
      if(.not.associated(gsd% lvi2).and.flg)allocate(gsd% lvi2(cpart%snpc ))
    elseif(vt==MPI_INTEGER4)then
      if(.not.associated(gsd%gsvi4).and.isr)allocate(gsd%gsvi4(gsi%gnpc))
      if(.not.associated(gsd% lvi4).and.flg)allocate(gsd% lvi4(cpart%snpc ))
    endif
  end subroutine vg_lfvar

#define TMPL_GATHER__kss(NA,VT,MVT)                                            \
  subroutine gather_kss_##NA  (lvar, gvar, km, k,root_)                       ;\
    integer,intent(in) ::km,k                                                 ;\
    VT,intent(in)      ::lvar(km,cpart%iblv:cpart%npc)                        ;\
    VT,intent(out)     ::gvar(gsi%ijm)                                        ;\
    integer ,optional  ::root_                                                ;\
    integer :: ierr,root,i                                                    ;\
    root=Verifyroot(root_)                                                    ;\
    call vg_lfvar(root,MVT,km/=1.or.gsi%parttype/=0)                          ;\
    if(km==1.and. gsi%parttype==0)then                                        ;\
      call MPI_GATHERV(lvar(1,1)  ,cpart%snpc,MVT,gsd%gsv##NA,                 \
        gsi%npcs,gsd%gdispls,MVT,root,gsi%mpi_comm,ierr)                      ;\
    else                                                                      ;\
      if(gsi%parttype==0)then                                                 ;\
        gsd%lv##NA=lvar(k,1:)                                                 ;\
      else                                                                    ;\
        do i=1,cpart%snpc                                                     ;\
          gsd%lv##NA(i)=lvar(k,cpart%calpnts(i))                              ;\
        enddo                                                                 ;\
      endif                                                                   ;\
      call MPI_GATHERV(gsd%lv##NA,cpart%snpc,MVT,gsd%gsv##NA,                  \
            gsi%npcs,gsd%gdispls,MVT,root,gsi%mpi_comm,ierr)                  ;\
    endif                                                                     ;\
    if(gsi%pid==root)then                                                     ;\
      do i=1,gsi%gnpc                                                         ;\
        gvar(gsd%pposg(i))=gsd%gsv##NA(i)                                     ;\
      enddo                                                                   ;\
    endif                                                                     ;\
  end subroutine gather_kss_##NA

#define TMPL_SCATTER_kss(NA,VT,MVT)                                            \
  subroutine scatter_kss_##NA(gvar, lvar, km, k,root_)                        ;\
    integer,intent(in) ::k,km                                                 ;\
    VT,intent(in)      ::gvar(gsi%ijm)                                        ;\
    VT,intent(out)     ::lvar(km,cpart%iblv:cpart%npc)                        ;\
    integer ,optional  ::root_                                                ;\
    integer :: ierr,root,i                                                    ;\
    root=Verifyroot(root_)                                                    ;\
    call vg_lfvar(root,MVT,km/=1.or.gsi%parttype/=0)                          ;\
    if(gsi%pid==root)then                                                     ;\
      do i=1,gsi%gnpc                                                         ;\
        gsd%gsv##NA(i)=gvar(gsd%pposg(i))                                     ;\
      enddo                                                                   ;\
    endif                                                                     ;\
    if(km==1.and. gsi%parttype==0)then                                        ;\
      DBGPV("SCATTERV A");\
      call MPI_SCATTERV(gsd%gsv##NA,gsi%npcs,gsd%gdispls,MVT,                  \
                 lvar(1,1),cpart%snpc,MVT,root,gsi%mpi_comm,ierr)             ;\
      DBGPV("SCATTERV AE");\
    else                                                                      ;\
      DBGPV("SCATTERV B");\
      call MPI_SCATTERV(gsd%gsv##NA,gsi%npcs,gsd%gdispls,MVT,                  \
             gsd%lv##NA,cpart%snpc,MVT,root,gsi%mpi_comm,ierr)                ;\
      DBGPV("SCATTERV BE");\
      if(gsi%parttype==0)then                                                 ;\
        lvar(k,1:)=gsd%lv##NA                                                 ;\
      else                                                                    ;\
        do i=1,cpart%snpc                                                     ;\
          lvar(k,cpart%calpnts(i))=gsd%lv##NA(i)                              ;\
        enddo                                                                 ;\
      endif                                                                   ;\
    endif                                                                     ;\
  end subroutine scatter_kss_##NA

#define TMPL_GATHER_(TAG,NA,VT,MVT,EXARGS,EXARGS1,LVSHAPE,GVSHAPE)             \
  subroutine gather_##TAG##NA  (lvar, gvar EXARGS ,root_)                     ;\
    VT,intent( in)::lvar LVSHAPE                                              ;\
    VT,intent(out)::gvar GVSHAPE                                              ;\
    integer ::root_ EXARGS                                                    ;\
    call gather_kss_##NA(lvar, gvar EXARGS1,root_)                            ;\
  end subroutine gather_##TAG##NA

#define TMPL_SCATTER(TAG,NA,VT,MVT,EXARGS,EXARGS1,LVSHAPE,GVSHAPE)             \
  subroutine scatter_##TAG##NA(gvar, lvar EXARGS ,root_)                      ;\
    VT,intent( in)::gvar GVSHAPE                                              ;\
    VT,intent(out)::lvar LVSHAPE                                              ;\
    integer ::root_ EXARGS                                                    ;\
    call scatter_kss_##NA(gvar, lvar EXARGS1,root_)                           ;\
  end subroutine scatter_##TAG##NA
TMPL_GATHER__kss(r4,real   (4),MPI_REAL4   )
TMPL_GATHER__kss(r8,real   (8),MPI_REAL8   )
TMPL_GATHER__kss(i2,integer(2),MPI_INTEGER2)
TMPL_GATHER__kss(i4,integer(4),MPI_INTEGER4)
TMPL_SCATTER_kss(r4,real   (4),MPI_REAL4   )
TMPL_SCATTER_kss(r8,real   (8),MPI_REAL8   )
TMPL_SCATTER_kss(i2,integer(2),MPI_INTEGER2)
TMPL_SCATTER_kss(i4,integer(4),MPI_INTEGER4)

#define GVS_S   (gsi%ijm      )
#define GVS_M   (gsi%im,gsi%jm)
#define LVSC_S   (   cpart%iblv:cpart%npc)
#define LVSC_M   (   cpart%nx  ,cpart%ny )
#define LVSC_KS  (km,cpart%iblv:cpart%npc)
#define LVSC_KM  (km,cpart%nx  ,cpart%ny )
#define LVSC_KLS (km,lm,cpart%iblv:cpart%npc)
#define LVSC_KLM (km,lm,cpart%nx  ,cpart%ny )
#define EXARG_KK  _ km _ k
#define EXARG_kk  _ km _ k
#define EXARG_KL  _ km _ k _ lm _ l
#define EXARG_kl  _ km*lm _ (l-1)*km+k
TMPL_GATHER_(  ss_,r4,real   (4),MPI_REAL4   ,        ,_  1 _ 1,LVSC_S,GVS_S)
TMPL_GATHER_(  ss_,r8,real   (8),MPI_REAL8   ,        ,_  1 _ 1,LVSC_S,GVS_S)
TMPL_GATHER_(  ss_,i2,integer(2),MPI_INTEGER2,        ,_  1 _ 1,LVSC_S,GVS_S)
TMPL_GATHER_(  ss_,i4,integer(4),MPI_INTEGER4,        ,_  1 _ 1,LVSC_S,GVS_S)

TMPL_SCATTER(  ss_,r4,real   (4),MPI_REAL4   ,        ,_  1 _ 1,LVSC_S,GVS_S)
TMPL_SCATTER(  ss_,r8,real   (8),MPI_REAL8   ,        ,_  1 _ 1,LVSC_S,GVS_S)
TMPL_SCATTER(  ss_,i2,integer(2),MPI_INTEGER2,        ,_  1 _ 1,LVSC_S,GVS_S)
TMPL_SCATTER(  ss_,i4,integer(4),MPI_INTEGER4,        ,_  1 _ 1,LVSC_S,GVS_S)

TMPL_GATHER_(  ms_,r4,real   (4),MPI_REAL4   ,        ,_  1 _ 1,LVSC_S,GVS_M)
TMPL_GATHER_(  ms_,r8,real   (8),MPI_REAL8   ,        ,_  1 _ 1,LVSC_S,GVS_M)
TMPL_GATHER_(  ms_,i2,integer(2),MPI_INTEGER2,        ,_  1 _ 1,LVSC_S,GVS_M)
TMPL_GATHER_(  ms_,i4,integer(4),MPI_INTEGER4,        ,_  1 _ 1,LVSC_S,GVS_M)

TMPL_SCATTER(  ms_,r4,real   (4),MPI_REAL4   ,        ,_  1 _ 1,LVSC_S,GVS_M)
TMPL_SCATTER(  ms_,r8,real   (8),MPI_REAL8   ,        ,_  1 _ 1,LVSC_S,GVS_M)
TMPL_SCATTER(  ms_,i2,integer(2),MPI_INTEGER2,        ,_  1 _ 1,LVSC_S,GVS_M)
TMPL_SCATTER(  ms_,i4,integer(4),MPI_INTEGER4,        ,_  1 _ 1,LVSC_S,GVS_M)

TMPL_GATHER_(  mm_,r4,real   (4),MPI_REAL4   ,        ,_  1 _ 1,LVSC_M,GVS_M)
TMPL_GATHER_(  mm_,r8,real   (8),MPI_REAL8   ,        ,_  1 _ 1,LVSC_M,GVS_M)
TMPL_GATHER_(  mm_,i2,integer(2),MPI_INTEGER2,        ,_  1 _ 1,LVSC_M,GVS_M)
TMPL_GATHER_(  mm_,i4,integer(4),MPI_INTEGER4,        ,_  1 _ 1,LVSC_M,GVS_M)

TMPL_SCATTER(  mm_,r4,real   (4),MPI_REAL4   ,        ,_  1 _ 1,LVSC_M,GVS_M)
TMPL_SCATTER(  mm_,r8,real   (8),MPI_REAL8   ,        ,_  1 _ 1,LVSC_M,GVS_M)
TMPL_SCATTER(  mm_,i2,integer(2),MPI_INTEGER2,        ,_  1 _ 1,LVSC_M,GVS_M)
TMPL_SCATTER(  mm_,i4,integer(4),MPI_INTEGER4,        ,_  1 _ 1,LVSC_M,GVS_M)
!=== kss is real function====
TMPL_GATHER_( kms_,r4,real   (4),MPI_REAL4   ,EXARG_KK,EXARG_kk,LVSC_KS,GVS_M)
TMPL_GATHER_( kms_,r8,real   (8),MPI_REAL8   ,EXARG_KK,EXARG_kk,LVSC_KS,GVS_M)
TMPL_GATHER_( kms_,i2,integer(2),MPI_INTEGER2,EXARG_KK,EXARG_kk,LVSC_KS,GVS_M)
TMPL_GATHER_( kms_,i4,integer(4),MPI_INTEGER4,EXARG_KK,EXARG_kk,LVSC_KS,GVS_M)

TMPL_SCATTER( kms_,r4,real   (4),MPI_REAL4   ,EXARG_KK,EXARG_kk,LVSC_KS,GVS_M)
TMPL_SCATTER( kms_,r8,real   (8),MPI_REAL8   ,EXARG_KK,EXARG_kk,LVSC_KS,GVS_M)
TMPL_SCATTER( kms_,i2,integer(2),MPI_INTEGER2,EXARG_KK,EXARG_kk,LVSC_KS,GVS_M)
TMPL_SCATTER( kms_,i4,integer(4),MPI_INTEGER4,EXARG_KK,EXARG_kk,LVSC_KS,GVS_M)

TMPL_GATHER_( kmm_,r4,real   (4),MPI_REAL4   ,EXARG_KK,EXARG_kk,LVSC_KM,GVS_M)
TMPL_GATHER_( kmm_,r8,real   (8),MPI_REAL8   ,EXARG_KK,EXARG_kk,LVSC_KM,GVS_M)
TMPL_GATHER_( kmm_,i2,integer(2),MPI_INTEGER2,EXARG_KK,EXARG_kk,LVSC_KM,GVS_M)
TMPL_GATHER_( kmm_,i4,integer(4),MPI_INTEGER4,EXARG_KK,EXARG_kk,LVSC_KM,GVS_M)

TMPL_SCATTER( kmm_,r4,real   (4),MPI_REAL4   ,EXARG_KK,EXARG_kk,LVSC_KM,GVS_M)
TMPL_SCATTER( kmm_,r8,real   (8),MPI_REAL8   ,EXARG_KK,EXARG_kk,LVSC_KM,GVS_M)
TMPL_SCATTER( kmm_,i2,integer(2),MPI_INTEGER2,EXARG_KK,EXARG_kk,LVSC_KM,GVS_M)
TMPL_SCATTER( kmm_,i4,integer(4),MPI_INTEGER4,EXARG_KK,EXARG_kk,LVSC_KM,GVS_M)

TMPL_GATHER_(klss_,r4,real   (4),MPI_REAL4   ,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_S)
TMPL_GATHER_(klss_,r8,real   (8),MPI_REAL8   ,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_S)
TMPL_GATHER_(klss_,i2,integer(2),MPI_INTEGER2,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_S)
TMPL_GATHER_(klss_,i4,integer(4),MPI_INTEGER4,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_S)

TMPL_SCATTER(klss_,r4,real   (4),MPI_REAL4   ,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_S)
TMPL_SCATTER(klss_,r8,real   (8),MPI_REAL8   ,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_S)
TMPL_SCATTER(klss_,i2,integer(2),MPI_INTEGER2,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_S)
TMPL_SCATTER(klss_,i4,integer(4),MPI_INTEGER4,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_S)

TMPL_GATHER_(klms_,r4,real   (4),MPI_REAL4   ,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_M)
TMPL_GATHER_(klms_,r8,real   (8),MPI_REAL8   ,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_M)
TMPL_GATHER_(klms_,i2,integer(2),MPI_INTEGER2,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_M)
TMPL_GATHER_(klms_,i4,integer(4),MPI_INTEGER4,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_M)

TMPL_SCATTER(klms_,r4,real   (4),MPI_REAL4   ,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_M)
TMPL_SCATTER(klms_,r8,real   (8),MPI_REAL8   ,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_M)
TMPL_SCATTER(klms_,i2,integer(2),MPI_INTEGER2,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_M)
TMPL_SCATTER(klms_,i4,integer(4),MPI_INTEGER4,EXARG_KL,EXARG_kl,LVSC_KLS,GVS_M)

TMPL_GATHER_(klmm_,r4,real   (4),MPI_REAL4   ,EXARG_KL,EXARG_kl,LVSC_KLM,GVS_M)
TMPL_GATHER_(klmm_,r8,real   (8),MPI_REAL8   ,EXARG_KL,EXARG_kl,LVSC_KLM,GVS_M)
TMPL_GATHER_(klmm_,i2,integer(2),MPI_INTEGER2,EXARG_KL,EXARG_kl,LVSC_KLM,GVS_M)
TMPL_GATHER_(klmm_,i4,integer(4),MPI_INTEGER4,EXARG_KL,EXARG_kl,LVSC_KLM,GVS_M)

TMPL_SCATTER(klmm_,r4,real   (4),MPI_REAL4   ,EXARG_KL,EXARG_kl,LVSC_KLM,GVS_M)
TMPL_SCATTER(klmm_,r8,real   (8),MPI_REAL8   ,EXARG_KL,EXARG_kl,LVSC_KLM,GVS_M)
TMPL_SCATTER(klmm_,i2,integer(2),MPI_INTEGER2,EXARG_KL,EXARG_kl,LVSC_KLM,GVS_M)
TMPL_SCATTER(klmm_,i4,integer(4),MPI_INTEGER4,EXARG_KL,EXARG_kl,LVSC_KLM,GVS_M)

#undef GVS_S
#undef GVS_M
#undef LVSC_S
#undef LVSC_M
#undef LVSC_KS
#undef LVSC_KM
#undef LVSC_KLS
#undef LVSC_KLM
#undef EXARG_KK
#undef EXARG_kk
#undef EXARG_KL
#undef EXARG_kl


#endif

!-------------------------------------------------------------------------------------------------

#ifndef NOSCATTER_EXT

  subroutine irrp_scatter_ext_init
    integer:: forceid
    type(force_info_type),pointer::fi
    fi=>scatterext
    fi%infp=cpart%snp
    fi%nfp=cpart%snpc
    allocate(fi%fpos(0:fi%infp))
    fi%fpos=cpart%ppos
    deallocate(fi%fpos)
    allocate(fi%fpos(0:fi%infp))
    fi%fpos=cpart%ppos
    call force_init(fi)
  end subroutine irrp_scatter_ext_init

  subroutine vf_lfvar(fi,vt,flg)
    type(force_info_type)::fi
    integer,intent(in)::vt
    logical,intent(in)::flg
    if(vt==MPI_REAL4)then
      if(.not.associated(fi%gsvr4)        ) allocate(fi%gsvr4(fi%infpg))
      if(.not.associated(fi% lvr4).and.flg)allocate(fi% lvr4(fi%infp ))
    elseif(vt==MPI_REAL8)then
      if(.not.associated(fi%gsvr8)        )allocate(fi%gsvr8(fi%infpg))
      if(.not.associated(fi% lvr8).and.flg)allocate(fi% lvr8(fi%infp ))
    elseif(vt==MPI_INTEGER2)then
      if(.not.associated(fi%gsvi2)        )allocate(fi%gsvi2(fi%infpg))
      if(.not.associated(fi% lvi2).and.flg)allocate(fi% lvi2(fi%infp ))
    elseif(vt==MPI_INTEGER4)then
      if(.not.associated(fi%gsvi4)        )allocate(fi%gsvi4(fi%infpg))
      if(.not.associated(fi% lvi4).and.flg)allocate(fi% lvi4(fi%infp ))
    endif
  end subroutine vf_lfvar

#define TMPL_SCATTER_EXT_kss(NA,VT,MVT)                                        \
  subroutine scatter_ext_kss_##NA(gvar, lvar, km, k, root_)                    ;\
    integer,intent(in) ::k, km                                                ;\
    VT,intent(in)      ::gvar(gsi%ijm)                                        ;\
    VT,intent(out)     ::lvar(km, cpart%iblv:cpart%np)                        ;\
    integer, intent(in), optional :: root_                                    ;\
    type(force_info_type), pointer :: fi                                      ;\
    integer :: idx, i, j, n, st, ierr, root                                   ;\
    root = Verifyroot(root_)                                                  ;\
    fi => scatterext;                                                         ;\
    if(fi%infp==0)call irrp_scatter_ext_init                                  ;\
    call vf_lfvar(fi, MVT, km/=1.or. gsi%parttype/=0 )                        ;\
    if(gsi%pid==root)then                                                     ;\
      do n = 1, fi%infpg                                                      ;\
        fi%gsv##NA(n) = gvar(fi%fposg(n))                                     ;\
      enddo                                                                   ;\
    endif                                                                     ;\
    if(km==1 .and. gsi%parttype==0)then                                       ;\
      DBGPV("SCATTERV C");\
      call MPI_SCATTERV( fi%gsv##NA,fi%infps,fi%fdispls,MVT,                   \
                         lvar(1,1),fi%infp,MVT,root,gsi%mpi_comm,ierr)        ;\
      DBGPV("SCATTERV CE");\
    else                                                                      ;\
      DBGPV("SCATTERV D");\
      call MPI_SCATTERV( fi%gsv##NA,fi%infps,fi%fdispls,MVT,                   \
                    fi%lv##NA,fi%infp,MVT,root,gsi%mpi_comm,ierr)             ;\
      DBGPV("SCATTERV DE");\
      if(gsi%parttype==0)then                                                 ;\
        lvar(k,1:)=fi%lv##NA                                                  ;\
      else                                                                    ;\
        do i=1,cpart%snp                                                      ;\
          lvar(k,cpart%calpnts(i))=fi%lv##NA(i)                               ;\
        enddo                                                                 ;\
      endif                                                                   ;\
    endif                                                                     ;\
  end subroutine scatter_ext_kss_##NA

#define TMPL_SCATTER_EXT(TAG,NA,VT,MVT,EXARGS,EXARGS1,LVSHAPE,GVSHAPE)         \
  subroutine scatter_ext_##TAG##NA(gvar, lvar EXARGS ,root_)                ;\
    VT,intent(in)      ::gvar GVSHAPE                                         ;\
    VT,intent(out)     ::lvar LVSHAPE                                         ;\
    integer :: root_ EXARGS                                                   ;\
    call scatter_ext_kss_##NA(gvar,lvar EXARGS1,root_)                         ;\
  end subroutine scatter_ext_##TAG##NA

  TMPL_SCATTER_EXT_kss(r4,real   (4),MPI_REAL4   )
  TMPL_SCATTER_EXT_kss(r8,real   (8),MPI_REAL8   )
  TMPL_SCATTER_EXT_kss(i2,integer(2),MPI_INTEGER2)
  TMPL_SCATTER_EXT_kss(i4,integer(4),MPI_INTEGER4)
#define GVS_S    (gsi%ijm      )
#define GVS_M    (gsi%im,gsi%jm)
#define LVSS_S   (      cpart%iblv:cpart%snp)
#define LVSS_M   (        cpart%nx,cpart%ny )
#define LVSS_KS  (km,   cpart%iblv:cpart%snp)
#define LVSS_KM  (km,     cpart%nx,cpart%ny )
#define LVSS_KLS (km,lm,cpart%iblv:cpart%snp)
#define LVSS_KLM (km,lm,  cpart%nx,cpart%ny )
#define EXARG_KK  _ km _ k
#define EXARG_kk  _ km _ k
#define EXARG_KL  _ km _ k _ lm _ l
#define EXARG_kl  _ km*lm _ (l-1)*km+k

TMPL_SCATTER_EXT(  ss_,r4,real   (4),MPI_REAL4   ,        ,_  1 _ 1,LVSS_S,GVS_S)
TMPL_SCATTER_EXT(  ss_,r8,real   (8),MPI_REAL8   ,        ,_  1 _ 1,LVSS_S,GVS_S)
TMPL_SCATTER_EXT(  ss_,i2,integer(2),MPI_INTEGER2,        ,_  1 _ 1,LVSS_S,GVS_S)
TMPL_SCATTER_EXT(  ss_,i4,integer(4),MPI_INTEGER4,        ,_  1 _ 1,LVSS_S,GVS_S)

TMPL_SCATTER_EXT(  ms_,r4,real   (4),MPI_REAL4   ,        ,_  1 _ 1,LVSS_S,GVS_M)
TMPL_SCATTER_EXT(  ms_,r8,real   (8),MPI_REAL8   ,        ,_  1 _ 1,LVSS_S,GVS_M)
TMPL_SCATTER_EXT(  ms_,i2,integer(2),MPI_INTEGER2,        ,_  1 _ 1,LVSS_S,GVS_M)
TMPL_SCATTER_EXT(  ms_,i4,integer(4),MPI_INTEGER4,        ,_  1 _ 1,LVSS_S,GVS_M)

TMPL_SCATTER_EXT(  mm_,r4,real   (4),MPI_REAL4   ,        ,_  1 _ 1,LVSS_M,GVS_M)
TMPL_SCATTER_EXT(  mm_,r8,real   (8),MPI_REAL8   ,        ,_  1 _ 1,LVSS_M,GVS_M)
TMPL_SCATTER_EXT(  mm_,i2,integer(2),MPI_INTEGER2,        ,_  1 _ 1,LVSS_M,GVS_M)
TMPL_SCATTER_EXT(  mm_,i4,integer(4),MPI_INTEGER4,        ,_  1 _ 1,LVSS_M,GVS_M)
!=== kss is real function====

TMPL_SCATTER_EXT( kms_,r4,real   (4),MPI_REAL4   ,EXARG_KK,EXARG_KK,LVSS_KS,GVS_M)
TMPL_SCATTER_EXT( kms_,r8,real   (8),MPI_REAL8   ,EXARG_KK,EXARG_KK,LVSS_KS,GVS_M)
TMPL_SCATTER_EXT( kms_,i2,integer(2),MPI_INTEGER2,EXARG_KK,EXARG_KK,LVSS_KS,GVS_M)
TMPL_SCATTER_EXT( kms_,i4,integer(4),MPI_INTEGER4,EXARG_KK,EXARG_KK,LVSS_KS,GVS_M)

TMPL_SCATTER_EXT( kmm_,r4,real   (4),MPI_REAL4   ,EXARG_KK,EXARG_KK,LVSS_KM,GVS_M)
TMPL_SCATTER_EXT( kmm_,r8,real   (8),MPI_REAL8   ,EXARG_KK,EXARG_KK,LVSS_KM,GVS_M)
TMPL_SCATTER_EXT( kmm_,i2,integer(2),MPI_INTEGER2,EXARG_KK,EXARG_KK,LVSS_KM,GVS_M)
TMPL_SCATTER_EXT( kmm_,i4,integer(4),MPI_INTEGER4,EXARG_KK,EXARG_KK,LVSS_KM,GVS_M)

TMPL_SCATTER_EXT(klss_,r4,real   (4),MPI_REAL4   ,EXARG_KL,EXARG_kl,LVSS_KLS,GVS_S)
TMPL_SCATTER_EXT(klss_,r8,real   (8),MPI_REAL8   ,EXARG_KL,EXARG_kl,LVSS_KLS,GVS_S)
TMPL_SCATTER_EXT(klss_,i2,integer(2),MPI_INTEGER2,EXARG_KL,EXARG_kl,LVSS_KLS,GVS_S)
TMPL_SCATTER_EXT(klss_,i4,integer(4),MPI_INTEGER4,EXARG_KL,EXARG_kl,LVSS_KLS,GVS_S)

TMPL_SCATTER_EXT(klms_,r4,real   (4),MPI_REAL4   ,EXARG_KL,EXARG_kl,LVSS_KLS,GVS_M)
TMPL_SCATTER_EXT(klms_,r8,real   (8),MPI_REAL8   ,EXARG_KL,EXARG_kl,LVSS_KLS,GVS_M)
TMPL_SCATTER_EXT(klms_,i2,integer(2),MPI_INTEGER2,EXARG_KL,EXARG_kl,LVSS_KLS,GVS_M)
TMPL_SCATTER_EXT(klms_,i4,integer(4),MPI_INTEGER4,EXARG_KL,EXARG_kl,LVSS_KLS,GVS_M)

TMPL_SCATTER_EXT(klmm_,r4,real   (4),MPI_REAL4   ,EXARG_KL,EXARG_kl,LVSS_KLM,GVS_M)
TMPL_SCATTER_EXT(klmm_,r8,real   (8),MPI_REAL8   ,EXARG_KL,EXARG_kl,LVSS_KLM,GVS_M)
TMPL_SCATTER_EXT(klmm_,i2,integer(2),MPI_INTEGER2,EXARG_KL,EXARG_kl,LVSS_KLM,GVS_M)
TMPL_SCATTER_EXT(klmm_,i4,integer(4),MPI_INTEGER4,EXARG_KL,EXARG_kl,LVSS_KLM,GVS_M)
#undef GVS_S
#undef GVS_M
#undef LVSS_S
#undef LVSS_M
#undef LVSS_KS
#undef LVSS_KM
#undef LVSS_KLS
#undef LVSS_KLM
#undef EXARG_KK
#undef EXARG_kk
#undef EXARG_KL
#undef EXARG_kl
  subroutine irrp_scatter_ext_final
    call scatter_force_final(scatterext)
  end subroutine irrp_scatter_ext_final

#endif

!-------------------------------------------------------------------------------------------------

#ifndef NOSCATTER_FORCE

  integer function sf_init(forceid)
    integer,intent(inout):: forceid
    type(force_info_type), pointer:: fs(:)
    integer mm,i
    if(forceid<1)then
      if(maxforces>0)then
        do i=1,maxforces
          if(forces(i)%nfp==0)then
            forceid=i
            exit
          endif
        enddo
      endif
      if(forceid<1)forceid=maxforces+1
    endif
    if(forceid>maxforces)then
      mm=maxforces;
      maxforces=forceid+4
      if(mm/=0)then
        fs=>forces
        allocate(forces(maxforces))
        forces(1:mm)=fs
        deallocate(fs)
      else
        allocate(forces(maxforces))
      endif
    endif
    sf_init=forceid
  end function sf_init

  subroutine force_init(fi)
    type(force_info_type),pointer::fi
    integer k,ierr
    allocate(fi%infps(gsi%npe),fi%fdispls(gsi%npe))
    fi%infps=-1
    call MPI_GATHER(fi%infp,1,mpi_integer4,fi%infps,1,mpi_integer4,0,gsi%mpi_comm,ierr)
    fi%infpg=sum(fi%infps)
    call MPI_BCAST(fi%infpg,1,MPI_INTEGER4,0,gsi%mpi_comm,ierr)
    call MPI_BCAST(fi%infps,gsi%npe,MPI_INTEGER4,0,gsi%mpi_comm,ierr)
    fi%fdispls(1)=0
    do k=1,gsi%npe-1
      fi%fdispls(k+1)=fi%fdispls(k)+fi%infps(k)
    enddo
    allocate(fi%fposg(0:fi%infpg))
    call MPI_GATHERV(fi%fpos(1),fi%infp,MPI_INTEGER,fi%fposg(1),fi%infps,fi%fdispls, &
                     MPI_INTEGER,0,gsi%mpi_comm,ierr)
    call MPI_BCAST(fi%fposg,(fi%infpg+1),MPI_INTEGER4,0,gsi%mpi_comm,ierr)
  end subroutine force_init

  subroutine scatter_force_final(fi)
    type(force_info_type)::fi
    fi%nfp=0;fi%infp=0;fi%infpg=0
    if(associated(fi%   fpos))deallocate(fi%   fpos)
    if(associated(fi%  fposg))deallocate(fi%  fposg)
    if(associated(fi%  infps))deallocate(fi%  infps)
    if(associated(fi%fdispls))deallocate(fi%fdispls)
    if(associated(fi%  gsvi2))deallocate(fi%  gsvi2)
    if(associated(fi%   lvi2))deallocate(fi%   lvi2)
    if(associated(fi%  gsvi4))deallocate(fi%  gsvi4)
    if(associated(fi%   lvi4))deallocate(fi%   lvi4)
    if(associated(fi%  gsvr4))deallocate(fi%  gsvr4)
    if(associated(fi%   lvr4))deallocate(fi%   lvr4)
    if(associated(fi%  gsvr8))deallocate(fi%  gsvr8)
    if(associated(fi%   lvr8))deallocate(fi%   lvr8)
    if(associated(fi%    idx))deallocate(fi%    idx)
    if(associated(fi%      w))deallocate(fi%      w)
  end subroutine scatter_force_final

  subroutine irrp_scatter_force_final_all
    integer:: forceid
    call irrp_scatter_ext_final
    do forceid=1,maxforces
      call scatter_force_final(forces(forceid))
    enddo
  end subroutine irrp_scatter_force_final_all

  subroutine irrp_scatter_force_final(forceid)
    integer,intent( in):: forceid
    if(forceid<0 .or. forceid>maxforces)return
    call scatter_force_final(forces(forceid))
  end subroutine irrp_scatter_force_final

  subroutine findi1i2(fp,n,x,i1,i2,p,fcycle)
    integer,intent(in):: n
 real   (8) ,intent(in):: x,fp(n),fcycle
    integer,intent(out):: i1,i2
 real   (8) ,intent(out):: p
    integer :: i
    if(fp(2)>=fp(1))then
      i1 = n
      do i=1,n
        if(x<fp(i))then
          i1 = i - 1; exit
        endif
      enddo
    else
      i1 = n
      do i=1,n
        if(x>fp(i))then
          i1 = i - 1; exit
        endif
      enddo
    endif
    i2=i1+1
    if(i1==0 .or. i1==n)then
      if(fcycle>1e-5)then
        i1 = n; i2 = 1; p=(x-fp(i1))/(fp(i2)+fcycle-fp(i1))
      else
        if(i1 == 0)i1 = 1
        i2 = i1; p = 0
      endif
    else
      p=(x-fp(i1))/(fp(i2)-fp(i1))
    endif
  end subroutine findi1i2

  subroutine InterInit1(fi,fix, fiy, fox, foy,fm,nix, niy,np,xcycle)
    integer,intent(in):: nix, niy, np
 real   (8),intent(in):: xcycle,fix(nix), fiy(niy),fox(cpart%iblv:cpart%np), foy(cpart%iblv:cpart%np)
    integer,intent(out):: fm(nix*niy)
    type(force_info_type)::fi
    integer i,j,n,i1,i2,j1,j2,m,ij
 real   (8) x,y,p,q
    fm=0;fi%idx=0;fi%w=0
    do n=1,np
      if(gsi%parttype==0)then
        x=fox(n);y=foy(n)
      else
        ij=cpart%calpnts(n)
        x=fox(ij);y=foy(ij)
      endif
      call findi1i2(fix,nix,x,i1,i2,p,xcycle)
      call findi1i2(fiy,niy,y,j1,j2,q,0.d0)
      ij=i1+(j1-1)*nix;fm(ij)=1; fi%idx(1,n)=ij; fi%w(1,n)=(1-p)*(1-q)
      ij=i2+(j1-1)*nix;fm(ij)=1; fi%idx(2,n)=ij; fi%w(2,n)=   p *(1-q)
      ij=i1+(j2-1)*nix;fm(ij)=1; fi%idx(3,n)=ij; fi%w(3,n)=(1-p)*   q
      ij=i2+(j2-1)*nix;fm(ij)=1; fi%idx(4,n)=ij; fi%w(4,n)=   p *   q
    enddo
  end subroutine InterInit1

  subroutine InterInit2(fi,fix, fiy, fox, foy,fm,nix, niy,np,xcycle)
    integer,intent(in):: nix, niy, np
 real   (8),intent(in):: xcycle, fix(nix, niy), fiy(niy, niy)
 real   (8),intent(in):: fox(cpart%iblv:cpart%np), foy(cpart%iblv:cpart%np)
    integer,intent(out):: fm(nix*niy)
    type(force_info_type)::fi
    write(6,*)'Orthogonal grid Force interpolation,Not surport Now'
    fm=0
  end subroutine InterInit2

   ! interface for r8
  subroutine irrp_scatter_force_init(forceid, fix_, fiy_, fox, foy, nix_, niy_, nnx_, nny_, npflg, xcycle, root_)
    integer, intent(inout) :: forceid       ! Index of forcing variable.
    integer, intent(in)    :: nix_, niy_    ! Size of data matrix.
    integer, intent(in)    :: npflg         ! flag for pnts need to prepare forcing.
                                            ! 0 for computer pnts only, else for all pnts.
    integer, intent(in)    :: nnx_, nny_    ! size of input coordinate data of forcing.
 real   (8), intent(in)    :: fix_(nnx_)    ! x-coordinate of input forcing. For curvlinear grid, nnx=nix*niy
 real   (8), intent(in)    :: fiy_(nny_)    ! y-coordinate of input forcing. For curvlinear grid, nnx=nix*niy
 real   (8), intent(in)    :: xcycle        ! The value of cycled in x direction, i.e. 360.
 real   (8), intent(in)    :: fox(cpart%iblv:cpart%np)  ! The x-coordinate of model.
 real   (8), intent(in)    :: foy(cpart%iblv:cpart%np)  ! The y-coordinate of model.
    integer, intent(in), optional :: root_  ! Root PE to scatter data.

    type(force_info_type), pointer :: fi
    integer, allocatable :: fm(:), ij2s(:)
    integer :: idx, i, j, n, k, st, ierr, ij, root
    integer :: nix, niy, nnx, nny
 real   (8), allocatable :: fix(:), fiy(:)
    type(mpipacket) :: pk
    root = Verifyroot(root_)
    forceid = sf_init(forceid)
    fi => forces(forceid)
    if(npflg==0)then
      fi%nfp = cpart%snpc
    else
      fi%nfp = cpart%snp
    endif
    if(associated(fi%idx))deallocate(fi%idx,fi%w)
    allocate(fi%idx(4,fi%nfp),fi%w(4,fi%nfp))
    if(gsi%pid==root)then
      call initmpipacket(pk,(4+nnx_+nny_+100)*8)  ! Use a larger number for safe.
      call MPI_PACK(nix_,   1,  MPI_INTEGER4, pk%buf, pk%bsize, pk%pos, gsi%mpi_comm, ierr)
      call MPI_PACK(niy_,   1,  MPI_INTEGER4, pk%buf, pk%bsize, pk%pos, gsi%mpi_comm, ierr)
      call MPI_PACK(nnx_,   1,  MPI_INTEGER4, pk%buf, pk%bsize, pk%pos, gsi%mpi_comm, ierr)
      call MPI_PACK(nny_,   1,  MPI_INTEGER4, pk%buf, pk%bsize, pk%pos, gsi%mpi_comm, ierr)
      call MPI_PACK(fix_, nnx_, MPI_REAL8   , pk%buf, pk%bsize, pk%pos, gsi%mpi_comm, ierr)
      call MPI_PACK(fiy_, nny_, MPI_REAL8   , pk%buf, pk%bsize, pk%pos, gsi%mpi_comm, ierr)
    endif
    call bcast_packet(pk, root, gsi%pid, gsi%mpi_comm)
    call MPI_UNPACK(pk%buf, pk%dsize, pk%pos, nix,   1, MPI_INTEGER4, gsi%mpi_comm, ierr)
    call MPI_UNPACK(pk%buf, pk%dsize, pk%pos, niy,   1, MPI_INTEGER4, gsi%mpi_comm, ierr)
    call MPI_UNPACK(pk%buf, pk%dsize, pk%pos, nnx,   1, MPI_INTEGER4, gsi%mpi_comm, ierr)
    call MPI_UNPACK(pk%buf, pk%dsize, pk%pos, nny,   1, MPI_INTEGER4, gsi%mpi_comm, ierr)
    allocate(fix(nnx), fiy(nny), ij2s(nix*niy), fm(nix*niy)); fm = 0
    call MPI_UNPACK(pk%buf, pk%dsize, pk%pos, fix, nnx, MPI_REAL8   , gsi%mpi_comm, ierr)
    call MPI_UNPACK(pk%buf, pk%dsize, pk%pos, fiy, nny, MPI_REAL8   , gsi%mpi_comm, ierr)
    fi%nx = nix; fi%ny = niy; fi%nxy = nix * niy
    if(nnx==nix)then
      call InterInit1(fi, fix, fiy, fox, foy, fm, nix, niy, fi%nfp, xcycle)
    else
      call InterInit2(fi, fix, fiy, fox, foy, fm, nix, niy, fi%nfp, xcycle)
    endif
    fi%infp = sum(fm)
    allocate(fi%fpos(0:fi%infp))
    idx = 0
    do ij = 1, nix*niy
        if(fm(ij) /= 0)then
          idx = idx + 1; ij2s(ij) = idx
          fi%fpos(idx) = ij
        endif
    enddo
    do n = 1, fi%nfp
      do k = 1, 4
        st = fi%idx(k, n) ; ! fi%idx(k,n) Set At InterInitx
        fi%idx(k, n) = ij2s(st)
      enddo
    enddo
    deallocate(fm, ij2s, fix, fiy)
    call force_init(fi)
  end subroutine irrp_scatter_force_init

  function GetForceInfo(forceid)
  type(force_info_type),pointer ::GetForceInfo
    !type(force_info_type),pointer::fi
    integer ,intent(in)::forceid
    if(forceid<1.or.forceid>maxforces)then
      write(6,*)'the forces Not Inited',forceid,maxforces
      call irrp_abort(__FILE__, __LINE__)
    endif
    !fi=>forces(forceid)
    GetForceInfo=>forces(forceid)
    if(GetForceInfo%nfp<=0)then
      write(6,*)'the force Not Inited',forceid,maxforces
      call irrp_abort(__FILE__, __LINE__)
    endif
  end function GetForceInfo

#define TMPL_SCATTER_FORCE_kssfi(NA,VT,MVT)                                    \
  subroutine sforce_kssfi_##NA(fi,gvar, lvar,km,k,root_)                      ;\
    type(force_info_type) ::fi                                                ;\
    integer, intent(in) ::k,km                                                ;\
    VT     , intent( in)::gvar(fi%nxy)                                        ;\
    VT     , intent(out)::lvar(km,cpart%iblv:cpart%np)                        ;\
    integer, intent(in), optional :: root_                                    ;\
    VT     , pointer    :: lv(:)                                              ;\
    real(8), pointer    ::  w(:,:)                                            ;\
    integer, pointer    ::idx(:,:)                                            ;\
    integer i,n,st,ierr,ij,root                                               ;\
    root=Verifyroot(root_)                                                    ;\
    call vf_lfvar(fi, MVT, .True. )                                           ;\
    if(gsi%pid==root)then                                                     ;\
      do n=1,fi%infpg                                                         ;\
        fi%gsv##NA(n)=gvar(fi%fposg(n))                                       ;\
      enddo                                                                   ;\
    endif                                                                     ;\
    lv=>fi%lv##NA                                                             ;\
    DBGPV("SCATTERV D");\
    call MPI_SCATTERV(fi%gsv##NA,fi%infps,fi%fdispls,MVT,                      \
                    lv,fi%infp,MVT,root,gsi%mpi_comm,ierr)                    ;\
    DBGPV("SCATTERV DE");\
    idx=>fi%idx;w=>fi%w                                                       ;\
    if(gsi%parttype==0)then                                                   ;\
        do n=1,fi%nfp                                                         ;\
          lvar(k, n)=lv(idx(1,n))*w(1,n)+lv(idx(2,n))*w(2,n)                   \
                  +lv(idx(3,n))*w(3,n)+lv(idx(4,n))*w(4,n)                    ;\
        enddo                                                                 ;\
      else                                                                    ;\
        do n=1,fi%nfp                                                         ;\
          ij=cpart%calpnts(n)                                                 ;\
          lvar(k,ij)=lv(idx(1,n))*w(1,n)+lv(idx(2,n))*w(2,n)                   \
                    +lv(idx(3,n))*w(3,n)+lv(idx(4,n))*w(4,n)                  ;\
        enddo                                                                 ;\
      endif                                                                   ;\
  end subroutine sforce_kssfi_##NA
  TMPL_SCATTER_FORCE_kssfi(r4,real   (4),MPI_REAL4   )
  TMPL_SCATTER_FORCE_kssfi(r8,real   (8),MPI_REAL8   )
  TMPL_SCATTER_FORCE_kssfi(i2,integer(2),MPI_INTEGER2)
  TMPL_SCATTER_FORCE_kssfi(i4,integer(4),MPI_INTEGER4)

#define TMPL_SCATTER_FORCE(TAG,NA,VT,MVT,EXARGS,EXARGS1,LVSHAPE,GVSHAPE)       \
  subroutine scatter_force_##TAG##NA(fid,gvar, lvar EXARGS,root_)          ;\
    integer,intent( in):: fid EXARGS                                          ;\
    VT     , intent( in)::gvar GVSHAPE                                        ;\
    VT     , intent(out)::lvar LVSHAPE                                        ;\
    integer, intent(in), optional :: root_                                    ;\
    call sforce_kssfi_##NA(forces(fid),gvar, lvar EXARGS1,root_)              ;\
  end subroutine scatter_force_##TAG##NA
#define FVS_S    (:  )
#define FVS_M    (:,:)
#define LVSS_S   (      cpart%iblv:cpart%snp)
#define LVSS_M   (        cpart%nx,cpart%ny )
#define LVSS_KS  (km,   cpart%iblv:cpart%snp)
#define LVSS_KM  (km,     cpart%nx,cpart%ny )
#define LVSS_KLS (km,lm,cpart%iblv:cpart%snp)
#define LVSS_KLM (km,lm,  cpart%nx,cpart%ny )
#define EXARG_KK  _ km _ k
#define EXARG_kk  _ km _ k
#define EXARG_KL  _ km _ k _ lm _ l
#define EXARG_kl  _ km*lm _ (l-1)*km+k
TMPL_SCATTER_FORCE(  ss_,r4,real   (4),MPI_REAL4   ,        ,_  1 _ 1,LVSS_S  ,FVS_S)
TMPL_SCATTER_FORCE(  ss_,r8,real   (8),MPI_REAL8   ,        ,_  1 _ 1,LVSS_S  ,FVS_S)
TMPL_SCATTER_FORCE(  ss_,i2,integer(2),MPI_INTEGER2,        ,_  1 _ 1,LVSS_S  ,FVS_S)
TMPL_SCATTER_FORCE(  ss_,i4,integer(4),MPI_INTEGER4,        ,_  1 _ 1,LVSS_S  ,FVS_S)

TMPL_SCATTER_FORCE(  ms_,r4,real   (4),MPI_REAL4   ,        ,_  1 _ 1,LVSS_S  ,FVS_M)
TMPL_SCATTER_FORCE(  ms_,r8,real   (8),MPI_REAL8   ,        ,_  1 _ 1,LVSS_S  ,FVS_M)
TMPL_SCATTER_FORCE(  ms_,i2,integer(2),MPI_INTEGER2,        ,_  1 _ 1,LVSS_S  ,FVS_M)
TMPL_SCATTER_FORCE(  ms_,i4,integer(4),MPI_INTEGER4,        ,_  1 _ 1,LVSS_S  ,FVS_M)

TMPL_SCATTER_FORCE(  mm_,r4,real   (4),MPI_REAL4   ,        ,_  1 _ 1,LVSS_M  ,FVS_M)
TMPL_SCATTER_FORCE(  mm_,r8,real   (8),MPI_REAL8   ,        ,_  1 _ 1,LVSS_M  ,FVS_M)
TMPL_SCATTER_FORCE(  mm_,i2,integer(2),MPI_INTEGER2,        ,_  1 _ 1,LVSS_M  ,FVS_M)
TMPL_SCATTER_FORCE(  mm_,i4,integer(4),MPI_INTEGER4,        ,_  1 _ 1,LVSS_M  ,FVS_M)

TMPL_SCATTER_FORCE( kss_,r4,real   (4),MPI_REAL4   ,EXARG_KK,EXARG_kk,LVSS_KS ,FVS_S)
TMPL_SCATTER_FORCE( kss_,r8,real   (8),MPI_REAL8   ,EXARG_KK,EXARG_kk,LVSS_KS ,FVS_S)
TMPL_SCATTER_FORCE( kss_,i2,integer(2),MPI_INTEGER2,EXARG_KK,EXARG_kk,LVSS_KS ,FVS_S)
TMPL_SCATTER_FORCE( kss_,i4,integer(4),MPI_INTEGER4,EXARG_KK,EXARG_kk,LVSS_KS ,FVS_S)

TMPL_SCATTER_FORCE( kms_,r4,real   (4),MPI_REAL4   ,EXARG_KK,EXARG_kk,LVSS_KS ,FVS_M)
TMPL_SCATTER_FORCE( kms_,r8,real   (8),MPI_REAL8   ,EXARG_KK,EXARG_kk,LVSS_KS ,FVS_M)
TMPL_SCATTER_FORCE( kms_,i2,integer(2),MPI_INTEGER2,EXARG_KK,EXARG_kk,LVSS_KS ,FVS_M)
TMPL_SCATTER_FORCE( kms_,i4,integer(4),MPI_INTEGER4,EXARG_KK,EXARG_kk,LVSS_KS ,FVS_M)

TMPL_SCATTER_FORCE( kmm_,r4,real   (4),MPI_REAL4   ,EXARG_KK,EXARG_kk,LVSS_KM ,FVS_M)
TMPL_SCATTER_FORCE( kmm_,r8,real   (8),MPI_REAL8   ,EXARG_KK,EXARG_kk,LVSS_KM ,FVS_M)
TMPL_SCATTER_FORCE( kmm_,i2,integer(2),MPI_INTEGER2,EXARG_KK,EXARG_kk,LVSS_KM ,FVS_M)
TMPL_SCATTER_FORCE( kmm_,i4,integer(4),MPI_INTEGER4,EXARG_KK,EXARG_kk,LVSS_KM ,FVS_M)

TMPL_SCATTER_FORCE(klss_,r4,real   (4),MPI_REAL4   ,EXARG_KL,EXARG_kl,LVSS_KLS,FVS_S)
TMPL_SCATTER_FORCE(klss_,r8,real   (8),MPI_REAL8   ,EXARG_KL,EXARG_kl,LVSS_KLS,FVS_S)
TMPL_SCATTER_FORCE(klss_,i2,integer(2),MPI_INTEGER2,EXARG_KL,EXARG_kl,LVSS_KLS,FVS_S)
TMPL_SCATTER_FORCE(klss_,i4,integer(4),MPI_INTEGER4,EXARG_KL,EXARG_kl,LVSS_KLS,FVS_S)

TMPL_SCATTER_FORCE(klms_,r4,real   (4),MPI_REAL4   ,EXARG_KL,EXARG_kl,LVSS_KLS,FVS_M)
TMPL_SCATTER_FORCE(klms_,r8,real   (8),MPI_REAL8   ,EXARG_KL,EXARG_kl,LVSS_KLS,FVS_M)
TMPL_SCATTER_FORCE(klms_,i2,integer(2),MPI_INTEGER2,EXARG_KL,EXARG_kl,LVSS_KLS,FVS_M)
TMPL_SCATTER_FORCE(klms_,i4,integer(4),MPI_INTEGER4,EXARG_KL,EXARG_kl,LVSS_KLS,FVS_M)

TMPL_SCATTER_FORCE(klmm_,r4,real   (4),MPI_REAL4   ,EXARG_KL,EXARG_kl,LVSS_KLM,FVS_M)
TMPL_SCATTER_FORCE(klmm_,r8,real   (8),MPI_REAL8   ,EXARG_KL,EXARG_kl,LVSS_KLM,FVS_M)
TMPL_SCATTER_FORCE(klmm_,i2,integer(2),MPI_INTEGER2,EXARG_KL,EXARG_kl,LVSS_KLM,FVS_M)
TMPL_SCATTER_FORCE(klmm_,i4,integer(4),MPI_INTEGER4,EXARG_KL,EXARG_kl,LVSS_KLM,FVS_M)
#undef FVS_S
#undef FVS_M
#undef LVSS_S
#undef LVSS_M
#undef LVSS_KS
#undef LVSS_KM
#undef LVSS_KLS
#undef LVSS_KLM
#endif
!-------------------------------------------------------------------------------------------------

  end module irrp_package_mod

!-------------------------------------------------------------------------------------------------
!#################################################################################################
