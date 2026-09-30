!#################################################################################################
!-------------------------------------------------------------------------------------------------
#define DBGINF   print*,__FILE__,__LINE__,gsi%pid
#define DBGO0(fmt) if(gsi%pid>=0)  write(6,fmt)
#define DBGB0(msg) !if(gsi%pid>=0)  write(6,'(2i5,f8.3,"s ",a)') __LINE__,gsi%pid,Difftimer(),msg;call flush(6)
module irrp_package_mod
  !-------------------------------------------------------------------------------------------------
  !use irrp_smpi_mod
  !use irrp_kernal_mod
  implicit none
  integer,parameter::kl=32,jnthet=12
  type pi_pos_type                          ! Self-defined type to record 2 dimensional index.
    integer :: i , j
  end type pi_pos_type
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
  public :: irrp_setpartserial              ! Set the pointers arranging method as Serial way.
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
#define TMPL_PUB(NA) public:: NA##_i2,NA##_i4,NA##_r4,NA##_r8;

  !TMPL_PUB(irrp_exg_setvar)
  !TMPL_PUB(irrp_gather)
  !TMPL_PUB(irrp_scatter)
  !TMPL_PUB(irrp_scatter_ext)
  !TMPL_PUB(irrp_scatter_force)

  integer,public,parameter::VT_R4=1
  integer,public,parameter::VT_R8=2
  integer,public,parameter::VT_I4=3
  integer,public,parameter::VT_I2=4
  private
!-------------------------------------------------------------------------------------------------
#define TMPL_INTERF(NA) \
  interface NA; module procedure \         
    NA##_i2,   NA##_i4,   NA##_r4,   NA##_r8   ; \
  end interface
  ! ks s m kls klm km
  TMPL_INTERF(irrp_exg_setvar)
  !ss ms mm kss kms kmm klss klms klmm
  TMPL_INTERF(irrp_gather)
  TMPL_INTERF(irrp_scatter)
  TMPL_INTERF(irrp_scatter_ext)
  TMPL_INTERF(irrp_scatter_force)
!-------------------------------------------------------------------------------------------------
#define OPTARG 
!,optional 
  interface
    subroutine irrp_exg_setvart(group_ind, ivar,var,ivt, km, vname)
      integer,intent(inout):: group_ind, ivar
      integer,intent(in):: ivt,km
      character*(*),intent(in):: vname
      class(*),dimension(*),intent(in) ::var
    end subroutine 
    subroutine irrp_gathert  (lvar, gvar,ivt, km, k,root_)
      class(*),dimension(*),intent(in) ::lvar
      class(*),dimension(*)::gvar
      integer,intent(in) ::ivt,km,k
      integer,intent(in)::root_
    end subroutine 
    subroutine irrp_scattert(gvar, lvar,ivt, km, k,root_) 
      integer, intent(in) ::ivt,k,km                      
      class(*),dimension(*),intent(in) ::gvar
      class(*),dimension(*)::lvar
      integer ,intent(out)  ::root_                       
    end subroutine 
    subroutine irrp_scatter_extt(gvar, lvar,ivt, km, k, root_)
      integer,intent(in) ::k, km                         
      class(*), dimension(*),intent(in) ::gvar
      class(*), dimension(*)::lvar
      integer, intent(in) :: root_           
    end subroutine 
    subroutine irrp_scatter_forcet(fid,gvar, lvar,ivt,km,k,root_)
      integer, intent(in) ::fid,ivt,k,km
      class(*),dimension(*),intent(in) ::gvar
      class(*),dimension(*)::lvar
      integer, intent(in) :: root_
    end subroutine 

    subroutine irrp_output_pposg()
    end subroutine 
    subroutine irrp_init_data()
    end subroutine 

    subroutine irrp_exginf(nst)
      integer nst
    end subroutine 
    subroutine irrp_final(mode)
      integer :: mode
    end subroutine 
    subroutine irrp_exg_init(group_ind, nvar_) ! default nvar = 16
      ! Init exchange_group,can
      integer,intent(inout) ::group_ind
      integer,intent(in) ::nvar_
    end subroutine 
    subroutine irrp_exg_action(group_ind,mode_)
      integer,intent(in)::group_ind
      integer,intent(in)::mode_
    end subroutine 
    subroutine irrp_exg_final(group_ind) ! default nvar = 16
      integer,intent(in)::group_ind
    end subroutine 
    subroutine irrp_exg_start(group_ind)
      integer,intent(in)::group_ind
    end subroutine 
    subroutine irrp_exg_check(group_ind,state)
      integer,intent(in)::group_ind
      integer,intent(out)::state
    end subroutine 
    subroutine irrp_exg_end(group_ind)
      integer,intent(in):: group_ind
    end subroutine 
    subroutine irrp_scatter_force_init(forceid, fix_, fiy_, fox, foy, nix_, niy_, nnx_, nny_, npflg, xcycle, root_)
      integer, intent(inout) :: forceid       ! Index of forcing variable.
      integer, intent(in)    :: nix_, niy_    ! Size of data matrix.
      integer, intent(in)    :: npflg         ! flag for pnts need to prepare forcing.
      ! 0 for computer pnts only, else for all pnts.
      integer, intent(in)    :: nnx_, nny_    ! size of input coordinate data of forcing.
      real(8), intent(in)    :: fix_(*)    ! x-coordinate of input forcing. For curvlinear grid, nnx=nix*niy
      real(8), intent(in)    :: fiy_(*)    ! y-coordinate of input forcing. For curvlinear grid, nnx=nix*niy
      real(8), intent(in)    :: xcycle        ! The value of cycled in x direction, i.e. 360.
      real(8), intent(in)    :: fox(*)  ! The x-coordinate of model.
      real(8), intent(in)    :: foy(*)  ! The y-coordinate of model.
      integer, intent(in)    :: root_  ! Root PE to scatter data.
    end subroutine 
    subroutine irrp_scatter_force_final(forceid)
      integer ::forceid
    end subroutine 
    subroutine irrp_scatter_ext_init()
    end subroutine 
    subroutine irrp_scatter_ext_final()
    end subroutine 
    subroutine irrp_initt(partmode, npe, pid, mpi_comm, mask,im,jm, halosize, cycle_flag, scycle,NS3)
      integer, intent(in   ) :: partmode
      integer, intent(in   ) :: npe
      integer, intent(in   ) :: pid,im,jm
      integer, intent(in   ) :: mpi_comm
      integer, intent(inout) :: mask(*) ! (im, jm)
      integer, intent(in   ) :: halosize
      integer, intent(in   ) :: cycle_flag
      integer, intent(in   ) :: scycle(*) ! (2,(im+jm+2)*2)
      integer, intent(in   ) :: NS3
    end subroutine 
    subroutine irrp_SetPartMatrixt(recto, ptype)
      integer, intent(out),optional:: recto(*)
      integer, intent(out),optional:: ptype(*)
    end subroutine 
    subroutine irrp_setpartserialt(gnpc, npc, np, plist, nb8, iblv_,nwps_)
      integer ,intent(out):: npc, np, gnpc,nwps_
      integer ,intent( in):: iblv_
      class(*) :: plist(*)
      class(*) :: nb8(*)
    end subroutine 
    subroutine irrp_getrectst(recti,recto,rects)
      integer,intent(out):: recti(4)
      integer,intent(out):: recto(4)
      integer,intent(out):: rects(*)
    end subroutine 
  end interface
contains
  subroutine irrp_init(partmode, npe, pid, mpi_comm, mask, halosize, cycle_flag, scycle)
    integer, intent(in) :: partmode
    integer, intent(in) :: npe
    integer, intent(in) :: pid
    integer, intent(in) :: mpi_comm
    integer, intent(inout) :: mask(:, :) ! (im, jm)
    integer, intent(in) ,optional :: halosize
    integer, intent(in) ,optional :: cycle_flag
    integer, intent(in) ,optional :: scycle(:, :) ! (2,(im+jm+2)*2)

    integer :: halosize_=1,NS3=0
    integer :: cycle_flag_=0,im=0,jm=0
    im       = size(mask, 1)            ! Size of grid matrix: first dimension.
    jm       = size(mask, 2)            ! Size of grid matrix: second dimension.
    if(present(halosize))halosize_ = halosize
    if(present(cycle_flag))then             ! Record cycle flag: default is 0.
      cycle_flag_ = cycle_flag
    elseif(present(scycle))then
      NS3=size(scycle, 2)
      cycle_flag_ = 3                    ! Determined by scycle:
                                         ! Cycled boundaries are given manually.
    else
      cycle_flag_ = 0                    ! Default setting: None cycle boundary.
    endif
    call irrp_initt(partmode, npe, pid, mpi_comm, mask,im,jm, halosize_, cycle_flag_, scycle,NS3)
  end subroutine 
  subroutine irrp_SetPartMatrix(recto, ptype)
    integer, intent(out),optional:: recto(4)
    integer, intent(out),optional,pointer:: ptype(:, :)
    integer ::nx,ny
    if(present(ptype))then
      call getnxny(nx,ny)
      allocate(ptype(nx, ny))
    endif
    call irrp_SetPartMatrixt(recto, ptype)
  end subroutine 
    subroutine irrp_setpartserial(gnpc, npc, np, plist, nb8, iblv_,nwps_)
      integer ,intent(out),optional :: npc, np, gnpc,nwps_
      integer ,intent( in),optional :: iblv_
      type(pi_pos_type)       ,intent(out),optional,  pointer:: plist(:)
      type(nb_8pnts_def_type) ,intent(out),optional,  pointer:: nb8(:)
      integer,external::getnp;
      integer np_
      np_=getnp();
      if(present(plist))then
        if(associated(plist))deallocate(plist)
        allocate(plist(0:np_));
      endif
      if(present(nb8))then
        if(associated(nb8))deallocate(nb8)
        allocate(nb8(0:np_)); 
      endif
    end subroutine 
    subroutine irrp_getrects(recti,recto,rects)
      integer,intent(out), optional:: recti(4)
      integer,intent(out), optional:: recto(4)
      integer,intent(out), optional,pointer:: rects(:,:)
      integer,external::getnpe;
      integer::npe
      if(present(rects))then
        npe=getnpe();
        allocate(rects(4,npe))       
      endif
      call irrp_getrectst(recti,recto,rects)
    end subroutine 
  ! exg_setvar
  ! Set var info into exchange group,4 var type,6 var shape
  ! all other shape reshape to (km,np)
  ! group_ind:exchange group index,<=0 auto find it
  ! ivar index of var in the group
  ! vname:var name
  ! var:var
  ! km :km
# define _ ,
#define EXARG_kk  _ km _ k
# define TMPL_exg_setvar(NA,VT)                            \
  subroutine irrp_exg_setvar_##NA(group_ind,ivar,vname,var,km);\
    integer:: group_ind, ivar                                 ;\
    integer,optional:: km                                     ;\
    character*(*),intent(in):: vname                          ;\
    VT:: var(*)                                               ;\
    call irrp_exg_setvart(group_ind, ivar, var,VT_##NA ,km, vname) ;\
  end subroutine 
  TMPL_exg_setvar(r4,real   (4)   ) 
  TMPL_exg_setvar(r8,real   (8)   )
  TMPL_exg_setvar(i2,integer(2)   )
  TMPL_exg_setvar(i4,integer(4)   )
#define TMPL_GATHER_(NA,VT)             \
  subroutine irrp_gather_##NA(lvar, gvar ,root_,km,k) ;\
    VT,intent( in)::lvar(*)                              ;\
    VT,intent(out)::gvar(*)                              ;\
    integer,optional::root_,km,k                         ;\
    call irrp_gathert(lvar, gvar,VT_##NA ,km,k,root_)   ;\
  end subroutine 
#define TMPL_SCATTER(NA,VT)             \
  subroutine irrp_scatter_##NA(gvar, lvar ,root_,km,k)  ;\
    VT,intent( in)::gvar(*)                              ;\
    VT,intent(out)::lvar(*)                              ;\
    integer,optional::root_,km,k                         ;\
    call irrp_scattert(gvar, lvar,VT_##NA ,km,k,root_)  ;\
  end subroutine 

  TMPL_GATHER_(r4,real   (4)) 
  TMPL_GATHER_(r8,real   (8))
  TMPL_GATHER_(i2,integer(2))
  TMPL_GATHER_(i4,integer(4))
  TMPL_SCATTER(r4,real   (4))
  TMPL_SCATTER(r8,real   (8))
  TMPL_SCATTER(i2,integer(2))
  TMPL_SCATTER(i4,integer(4))

#define TMPL_SCATTER_EXT(NA,VT)           \
  subroutine irrp_scatter_ext_##NA(gvar, lvar ,root_,km,k);\
    VT,intent(in)      ::gvar(*)                           ;\
    VT,intent(out)     ::lvar(*)                           ;\
    integer,optional::root_,km,k                           ;\
    call irrp_scatter_extt(gvar,lvar,VT_##NA,km,k ,root_) ;\
  end subroutine 
  TMPL_SCATTER_EXT(r4,real   (4))
  TMPL_SCATTER_EXT(r8,real   (8))
  TMPL_SCATTER_EXT(i2,integer(2))
  TMPL_SCATTER_EXT(i4,integer(4))
#define TMPL_SCATTER_FORCE(NA,VT)           \
  subroutine irrp_scatter_force_##NA(fid,gvar, lvar ,root_,km,k) ;\
    integer,intent( in):: fid                               ;\
    VT, intent( in)::gvar(*)                                ;\
    VT, intent(out)::lvar(*)                                ;\
    integer,optional::root_,km,k                            ;\
    call irrp_scatter_forcet(fid,gvar, lvar,VT_##NA,km,k ,root_);\
  end subroutine 
  TMPL_SCATTER_FORCE(r4,real   (4))
  TMPL_SCATTER_FORCE(r8,real   (8))
  TMPL_SCATTER_FORCE(i2,integer(2))
  TMPL_SCATTER_FORCE(i4,integer(4))
end module irrp_package_mod
!#################################################################################################
