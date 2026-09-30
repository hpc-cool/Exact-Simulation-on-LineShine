#include "wavedeff.h"
#define RINFO(TAG) if(mpi_id==0)write(*,'(i8," ",a,f8.3," ",i8)')iwalltime(),TAG,Difftimer(2),mpi_npe;  call flush(6);
Module partition_mod
#ifndef NO_MPI
  use wav_mpi_mod
  use irrp_kernal_mod
  use irrp_package_mod
#endif
  use debughlp_mod
  use varcommon_mod
  IMPLICIT NONE
  integer ::gnwpc
  ! Globle Work Point Pos (2,gnwpc)
  integer ,pointer:: ieposc(:,:)
  ! Area Rects (1:npe]
  integer,pointer::rects(:,:)
  integer ::rect(4),recto(4)
  integer LXB,LYB,LXN,LYN
  integer,allocatable:: ieind(:,:) ! gixl giyl
  integer,allocatable:: iepos(:,:) ! (2,nwpa)
#ifndef NO_MPI
  type(pi_pos_type), pointer:: plist(:)
#else
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
#endif
  type(nb_8pnts_def_type), pointer :: nb8(:)
  INTEGER ,pointer::nwpcs(:)
  public  :: nb8,partition,DeinitPart,exchange_boundary,GatherVar2d,scattervar2d
  public  :: scatter_wind_init
  INTEGER ::exginited(2)=0,partition_firstcall=1
  INTEGER ::windforceid=-1
  INTEGER ::curexg=1

CONTAINS
#ifndef NO_MPI
  subroutine partition
    character(256) fn
    INTEGER iac,n,i,j
    if(partition_firstcall/=0)then
      partition_firstcall    =0
      nullify(ieposc ,nwpcs  ,plist  ,nb8    ,rects  )
    endif
    exginited=0
    !call SerRB;    print*,"partition  begin",mpi_id,trim(nodename);call flush(6);call SerRE
    DBGO0(1,*)'irrp_init begin'
#define DPART
#ifdef DPART
    if(mpi_id==0)then
      open(11,file="pemaskn.dat",status="unknown")
      WRITE(11,*)gixl,giyl
    endif
#endif
    call irrp_init(0, mpi_npe, mpi_id, mpi_comm_wav, nspg, 1, Gcircle)    !here nspg  is pemask
    DBGO0(1,*)'irrp_init end'
#ifdef DPART
    if(mpi_id==0)then
      WRITE(11,*)"nspg"
      DO  j=1,giyl !{
      WRITE(11,'(10000i2)')nspg(:,j)
      end do !}
      WRITE(11,*)"depg"
      DO  j=1,giyl !{
      WRITE(11,'(10000f7.1)')depg(:,j)
      end do !}
    endif
#endif
    call ResetNspg(nspg) !reset nspg for SetIDPOS
    call irrp_SetPartSerial(gnwpc, nwpc, nwpa, plist, nb8,0,nwps)
    call irrp_getrects(rect,recto,rects)
    call irrp_final(0)
    call SetLocalArea(recto,nwpa)

    allocate(ieposc(2,0:nwpc));
    allocate(nwpcs(mpi_npe))
    !call SerRB;    print*,"gather nwpc begin",mpi_id,trim(nodename);call flush(6);call SerRE
    RINFO('gather nwpc begin')
    call wav_mpi_gather(nwpc,nwpcs)
    RINFO('gather nwpc end')
    call wav_mpi_bcast(nwpcs)
    RINFO('bcast nwpc end')
    gnwpc=sum(nwpcs)
    do iac=1,nwpc
      ieposc(1,iac)=plist(iac)%i
      ieposc(2,iac)=plist(iac)%j
      call SetIDPOS(0,iac,plist(iac)%i,plist(iac)%j,0,nspg(plist(iac)%i,plist(iac)%j))
    enddo
    do iac=nwpc+1,nwpa
      call SetIDPOS(1,iac,plist(iac)%i,plist(iac)%j,0,nspg(plist(iac)%i,plist(iac)%j))
    enddo
    call settopog
    RINFO('settopog End')
    if(allocated(depg)) deallocate(depg)
    if(allocated(nspg)) deallocate(nspg) ! nspg is no used
    call irrp_init_data
    RINFO('irrp_init_data End')

    !DBGINF;call zbarr
    DBGO0(1,*)'partition End'
    !call InitMpi(-1)      ;    call flush(6)   ;    stop
    !call irrp_exginf;    stop
    !   4---3---2
    !   |   |   |
    !   5---+---1
    !   |   |   |
    !   6---7---8
    !call irrp_output_pposg
    call wav_mpi_abort("TTT");

  end subroutine partition
  subroutine DeinitPart
    if(associated(ieposc ))deallocate(ieposc )
    if(associated(nwpcs  ))deallocate(nwpcs  )
    if(associated(plist  ))deallocate(plist  )
    if(associated(nb8    ))deallocate(nb8    )
    if(associated(rects  ))deallocate(rects  )
    call irrp_final(1)
    exginited=0
    partition_firstcall=1
  end subroutine DeinitPart
  subroutine testgather
    integer it,its(16)
    call wav_mpi_gather(it,its)
  end subroutine testgather
  subroutine GatherVar2d(gvar,lvar,root_)
    real(4),intent(in)::lvar(0:nwpc)
    real(4),intent(out)::gvar(gixl,giyl)
    integer,optional::root_
    call irrp_gather(lvar,gvar,root_)
  end subroutine GatherVar2d
  subroutine ScatterVar2d(gvar,lvar,root_)
    real(4),intent(out)::lvar(0:nwpc)
    real(4),intent(in)::gvar(gixl,giyl)
    integer,optional::root_
    call irrp_scatter(gvar,lvar,root_)
  end subroutine ScatterVar2d
  subroutine exchange_boundary_start(ee,exg)
    integer exg
    REALD ,intent(inout)::ee(kl,jnthet,0:*)
    curexg=exg
    if(exginited(curexg)==0)then
      call SetExgVars1(ee)
    endif
    call irrp_exg_start(curexg)
  end subroutine exchange_boundary_start
  subroutine exchange_boundary_check(exg_state)
    integer ,intent(out):: exg_state
    if(exg_state==0)call irrp_exg_check(curexg,exg_state)
  end subroutine exchange_boundary_check
  subroutine exchange_boundary_end()
    call irrp_exg_end(curexg)
  end subroutine exchange_boundary_end
  subroutine exchange_boundary(ee,exg)
    integer exg
    REALD ,intent(inout)::ee(kl,jnthet,0:*)
    curexg=exg
    if(exginited(curexg)==0)then
      call SetExgVars1(ee)
    endif
    call irrp_exg_action(curexg,7)
  end subroutine exchange_boundary
  subroutine SetExgVars1(ee)
    REALD ,intent(inout):: ee(kl,jnthet,0:*)
    integer vid
    vid=-1
    call irrp_exg_init(curexg,1)
    if(curexg==1)then
      call irrp_exg_setvar(curexg,vid,'ee1',ee,kl*jnthet)
    else
      call irrp_exg_setvar(curexg,vid,'ee2',ee,kl*jnthet)
    endif
    exginited(curexg)=1;
  end subroutine SetExgVars1
  subroutine scatter_wind_init(fix_, fiy_, nix_, niy_, nnx_, nny_,  xcycle)
    integer, intent(in)    :: nix_, niy_    ! Size of data matrix.
    integer, intent(in)    :: nnx_, nny_    ! size of input coordinate data of forcing.
    real(8), intent(in)    :: fix_(nnx_)    ! x-coordinate of input forcing. For curvlinear grid, nnx=nix*niy
    real(8), intent(in)    :: fiy_(nny_)    ! y-coordinate of input forcing. For curvlinear grid, nnx=nix*niy
    real(8), intent(in)    :: xcycle        ! The value of cycled in x direction, i.e. 360.
    call irrp_scatter_force_init(windforceid, fix_, fiy_, sxcord, sycord, nix_, niy_, nnx_, nny_, 0, xcycle,0)
  end subroutine scatter_wind_init
  subroutine scatter_wind(windi, windx_ )
    real*4 , intent( in)::windi (:,:)
    real*4 , intent(out)::windx_ (:)
    call irrp_scatter_force(windforceid,windi,windx_,0)
  end subroutine scatter_wind
  subroutine InitMpi(mode)
    integer mode
    call wav_init_mpi(mode)      ! Mpi Deinit
  end subroutine InitMpi

  SUBROUTINE Barrier
    call wav_mpi_barrier
  end SUBROUTINE Barrier


#else
  !ZZZ here
  subroutine partition
    INTEGER iac,n,ia,ic,ial,iar,icu,icd
    if(partition_firstcall/=0)then
      partition_firstcall    =0
      nullify(ieposc ,nwpcs  ,nb8    ,rects  )
    endif
    nwpc=0
    call ResetNspg(nspg) !reset nspg for SetIDPOS
    do ic=1,giyl
      do ia=1,gixl
        if(nspg(ia,ic)/=0)nwpc=nwpc+1
      enddo
    enddo
    gnwpc=nwpc;nwpa=nwpc;nwps=0
    rect=[1,1,gixl,giyl]
    recto=rect
    call SetLocalArea(recto,nwpa)
    allocate(ieposc(2,0:nwpc),nb8(0:nwpc),rects(4,1));
    rects(:,1)=rect
    iac=0
    do ic=1,giyl
      do ia=1,gixl
        if(nspg(ia,ic)/=0)then
          iac=iac+1
          ieposc(1,iac)=ia
          ieposc(2,iac)=ic
          call SetIDPOS(1,iac,ia,ic,0,nspg(ia,ic))
        endif
      enddo
    enddo
    call settopog
    do iac=1,nwpc
      call id2pos(iac,ia,ic)
      ial=adjix(ia-1);icd=adjiy(ic-1)
      iar=adjix(ia+1);icu=adjiy(ic+1)
      if(Gcircle/=0)then
        if(ial<rectc(1))ial=rectc(3); !ZZZEERR
        if(iar>rectc(3))iar=rectc(1); !ZZZEERR
      else
        if(ial<rectc(1))ial=ia; !ZZZEERR
        if(iar>rectc(3))iar=ia; !ZZZEERR
      endif
      if(icd<rectc(2).or.icd>rectc(4))icd=ic; !ZZZEERR
      if(icu<rectc(2).or.icu>rectc(4))icu=ic; !ZZZEERR
      nb8(iac)%r =pos2id(iar,ic )
      nb8(iac)%ur=pos2id(iar,icu)
      nb8(iac)%u =pos2id(ia ,icu)
      nb8(iac)%ul=pos2id(ial,icu)
      nb8(iac)%l =pos2id(ial,ic )
      nb8(iac)%dl=pos2id(ial,icd)
      nb8(iac)%d =pos2id(ia ,icd)
      nb8(iac)%dr=pos2id(iar,icd)
    enddo
    DBGO0(1,*)'partition End'
    !call irrp_exginf;    stop
  end subroutine partition
  subroutine DeinitPart
    if(associated(ieposc ))deallocate(ieposc )
    if(associated(nwpcs  ))deallocate(nwpcs  )
    if(associated(nb8    ))deallocate(nb8    )
    if(associated(rects  ))deallocate(rects  )
    partition_firstcall=1
  end subroutine DeinitPart
  subroutine GatherVar2d(gvar,lvar,root_)
    real(4),intent(in)::lvar(0:nwpc)
    real(4),intent(out)::gvar(gixl,giyl)
    integer,optional::root_
    integer iac,ia,ic
    do iac=1,nwpc
      ia=ieposc(1,iac)
      ic=ieposc(2,iac)
      gvar(ia,ic)=lvar(iac)
    enddo
  end subroutine GatherVar2d
  subroutine ScatterVar2d(gvar,lvar,root_)
    real(4),intent(out)::lvar(0:nwpc)
    real(4),intent(in)::gvar(gixl,giyl)
    integer,optional::root_
    integer iac,ia,ic
    do iac=1,nwpc
      ia=ieposc(1,iac)
      ic=ieposc(2,iac)
      lvar(iac)=gvar(ia,ic)
    enddo
  end subroutine ScatterVar2d
  subroutine exchange_boundary(ee,mode)
    REALD ee(kl,jnthet,0:*)
    integer mode
  end subroutine exchange_boundary
  subroutine scatter_wind_init(fix_, fiy_, nix_, niy_, nnx_, nny_,  xcycle)
    integer, intent(in)    :: nix_, niy_    ! Size of data matrix.
    integer, intent(in)    :: nnx_, nny_    ! size of input coordinate data of forcing.
    real(8), intent(in)    :: fix_(nnx_)    ! x-coordinate of input forcing. For curvlinear grid, nnx=nix*niy
    real(8), intent(in)    :: fiy_(nny_)    ! y-coordinate of input forcing. For curvlinear grid, nnx=nix*niy
    real(8), intent(in)    :: xcycle        ! The value of cycled in x direction, i.e. 360.
    write(6,*)'Not Surport '
    stop
  end subroutine scatter_wind_init
  subroutine scatter_wind(rwind, swind )
    real*4 , intent( in)::rwind (:,:)
    real*4 , intent(out)::swind (:)
    integer iac,ia,ic
    do iac=1,nwpc
      ia=ieposc(1,iac)
      ic=ieposc(2,iac)
      swind(iac)=rwind(ia,ic)
    enddo
  end subroutine scatter_wind
  SUBROUTINE Barrier
  end SUBROUTINE Barrier
  subroutine InitMpi(mode)
    integer mode
    mpi_id=0
    mpi_npe=1
  end subroutine InitMpi
#endif
  subroutine set_neighbor
    INTEGER iac
    if(.not. associated(ipos8))then
      ALLOCATE(ipos8(0:nwpc))
    endif

    if(.not. associated(ipos12))then
      ALLOCATE(ipos12(12,0:nwpc))
    endif
    ipos12=0;
    do iac=1,nwpc
      !ipos8(iac)%dep=dep(iac)
      !ipos8(iac)%i8( 0)=    iac;
      ipos8(iac)%i8( 1)=nb8(iac)%r;
      ipos8(iac)%i8( 2)=nb8(iac)%ur;
      ipos8(iac)%i8( 3)=nb8(iac)%u;
      ipos8(iac)%i8( 4)=nb8(iac)%ul;
      ipos8(iac)%i8( 5)=nb8(iac)%l;
      ipos8(iac)%i8( 6)=nb8(iac)%dl;
      ipos8(iac)%i8( 7)=nb8(iac)%d;
      ipos8(iac)%i8( 8)=nb8(iac)%dr;
      !i  2  3  4
      !132534178576
      ipos12( 1,iac)=ipos8(iac)%i8(1);! nb8(iac)%r;
      ipos12( 2,iac)=ipos8(iac)%i8(3);! nb8(iac)%u;
      ipos12( 3,iac)=ipos8(iac)%i8(2);! nb8(iac)%ur;
      ipos12( 4,iac)=ipos8(iac)%i8(5);! nb8(iac)%l;
      ipos12( 5,iac)=ipos8(iac)%i8(3);! nb8(iac)%u;
      ipos12( 6,iac)=ipos8(iac)%i8(4);! nb8(iac)%ul;
      ipos12( 7,iac)=ipos8(iac)%i8(1);! nb8(iac)%r;
      ipos12( 8,iac)=ipos8(iac)%i8(7);! nb8(iac)%d;
      ipos12( 9,iac)=ipos8(iac)%i8(8);! nb8(iac)%dr;
      ipos12(10,iac)=ipos8(iac)%i8(5);! nb8(iac)%l;
      ipos12(11,iac)=ipos8(iac)%i8(7);! nb8(iac)%d;
      ipos12(12,iac)=ipos8(iac)%i8(6);! nb8(iac)%dl;
    enddo
  end subroutine set_neighbor
  integer function AdjIx(ix)
    integer ix
    AdjIx=ix
    if(Gcircle/=0)then
      if(AdjIx>gixl)AdjIx=AdjIx-gixl
      if(AdjIx<1   )AdjIx=AdjIx+gixl
    else
      !if(AdjIx>gixl.or.AdjIx<1   )then
      !  DBGO(0,*)mpi_id,'AdjIx Error',ix,gixl
      !endif
      if(AdjIx>gixl)AdjIx=gixl
      if(AdjIx<1   )AdjIx=1
    end if
  end function AdjIx
  integer function AdjIy(iy)
    integer iy
    AdjIy=iy
    !if(AdjIy>giyl.or.AdjIy<1   )then
    !  DBGO(0,*)mpi_id,'AdjIy Error',iy,giyl
    !endif
    if(AdjIy>giyl)AdjIy=giyl
    if(AdjIy<1   )AdjIy=1
  end function AdjIy
  logical function intersect(rectt)
    integer rectt(4)
    intersect=.false.;
    if(rectc(1)>rectt(3))return
    if(rectc(3)<rectt(1))return
    if(rectc(2)>rectt(4))return
    if(rectc(4)<rectt(2))return
    intersect=.true.
  end function intersect
  logical function innerArea(ix,iy)
    integer ix,iy
    innerarea=.false.
    if(rectc(1)>ix.or.rectc(3)<ix)return
    if(rectc(2)>iy.or.rectc(4)<iy)return
    innerArea=.true.
  end function innerArea
  subroutine SetLocalArea(rect,Nwp_)
    integer rect(4),Nwp_
    rectc=rect
    LXB=rect(1);LYB=rect(2);LXN=rect(3)-rect(1)+1;LYN=rect(4)-rect(2)+1
    allocate(ieind(LXN,LYN),iepos(2,0:nwp_),nsp(0:nwp_));
    !write(*,*)__FILE__,__LINE__,rect,LXB,LXN
    ieind=0;iepos=0;nsp=0
  end subroutine SetLocalArea
  subroutine SetIDPOS(flag,ind,ix_,iy,eqid,insp)
    integer ind,ix_,iy,eqid,insp,flag
    integer ix,it
    ix=ix_;
    if(ind<=0.and.ind>nwpa)then
      SDBGINF,'IDPOS Error ',mpi_id,ind,ix,iy,eqid,insp
    endif
    if(ix<LXB)then
      ix=ix+gixl
    endif
    if(ix>LXB+LXN)then
      ix=ix-gixl
    endif
    if(ix<LXB)then
      write(6,*) 'Error SetIDPOS ',mpi_id,ind,nwps,nwpc,nwpa,"ix=",ix,ix_,' gixl=',gixl,LXB,LXB+LXN
    endif
    if(flag/=0 .and. LXN>gixl)then
      if    (ix==1   )then;ix=gixl+1
      elseif(ix==gixl)then;ix=0
      endif
    endif
    ! if(mpi_id==0)then
    !   write(6,*)'VV',ind,ix,iy,LXB,LXB+LXN
    ! endif
    it=ieind(ix-LXB+1,iy-LYB+1)
    !if(ix_==360.and.iy==9)then
    ! print*,ix,ix_,iy,ind,insp,eqid,flag,gixl,LXN
    !endif
    !if(ix_==360.and.iy==9)then
    ! print*,ix,ix_,iy,ind,insp,eqid,flag,gixl,LXN
    !endif
    if(it/=0)then !Fixme
      !write(*,'(a,20i8)')'SetIDPOS' ,mpi_id,ix_,iy,ind,nwpa,it,nwpc,iepos(:,it),eqid
    endif
    if(eqid==0)then
      ieind(ix-LXB+1,iy-LYB+1)=ind
      nsp(ind)=insp
    endif
    iepos(1,ind)=ix
    iepos(2,ind)=iy
  end subroutine SetIDPOS
  integer function pos2id(ix,iy)
    integer ix,iy
    if(innerArea(ix,iy))then
      pos2id=ieind(ix-LXB+1,iy-LYB+1)
    else
      pos2id=0
    endif
  end function pos2id
  SUBROUTINE id2pos(id,ix,iy)
    integer id,ix,iy
    ix=iepos(1,id)
    iy=iepos(2,id)
  end SUBROUTINE id2pos
  SUBROUTINE SetTopog
    character fmtf*32
    real(8) dymt,dxmt
    real(8) dx,dy,dxyl
    integer ia,ic,iac,i,j,i1,j1
    integer ierr
    ALLOCATE(xcord(gixl,1),ycord(1,giyl),STAT=ierr) ;
    ALLOCATE(sxcord(0:nwpa),sycord(0:nwpa),STAT=ierr)
    ALLOCATE(Rsd_tanLat(0:nwpc),dxm(0:nwpc),dym(0:nwpc),dddx(0:nwpc),dddy(0:nwpc),STAT=ierr)
    if(.not. associated(dep))then
      ALLOCATE(dep(0:nwpc),STAT=ierr)
    endif
    xcord=0;ycord=0;dep=0;dxm=0;dddx=0;dddy=0
    ! extract Data
    DO  ic=1,giyl !{
      DO  ia=1,gixl !{
        IF(depg(ia,ic)>0 .and. depg(ia,ic)<5.) depg(ia,ic)=5.
      end do !}
    end do !}
    do iac=1,nwpc
      if(nsp(iac)<=0)cycle
      call id2pos(iac,ia,ic)
      i1=adjix(ia);j1=AdjIy(ic)
      dep(iac)=depg(i1,j1)
    end do !}
    dymt=deg2m*grdszy
    dxmt=deg2m*grdszx
    DO  ia=1,gixl !{
      xcord(ia,1)=gxlon0+grdszx*(ia-1)
    end do !}
    DO  ic=1,giyl !{
      ycord(1,ic)=gylat0+grdszy*(ic-1)
    end do !}
    do iac=1,nwpc
      if(nsp(iac)<=0)cycle
      call id2pos(iac,ia,ic)
      if(ic<=0 .or. ic>giyl)then
        write(6,*)'id2pos Error',iac,ia,ic
      endif
      sxcord(iac)=xcord(ia,1)
      sycord(iac)=ycord(1,ic)
      dym(iac)=dymt
      Rsd_tanLat(iac)=tand(ycord(1,ic))/Rs
      dx=dxmt*(cosd(ycord(1,ic)))
      dxm(iac)=dx
    end do !}
    do iac=1,nwpc
      if(nsp(iac)<=0)cycle
      call id2pos(iac,ia,ic)

      if(ia==1)then
        dddx(iac)=(depg(ia+1,ic)-depg(gixl,ic))/(2*dxm(iac))
      else if(ia==gixl)then
        dddx(iac)=(depg(1,ic)-depg(ia-1,ic))/(2*dxm(iac))
      else
        dddx(iac)=(depg(ia+1,ic)-depg(ia-1,ic))/(2*dxm(iac))
      endif
      if(ic==1)then
        dddy(iac)=(depg(ia,ic+1)-depg(ia,1))/(dym(iac))
      else if(ic==giyl)then
        dddy(iac)=(depg(ia,giyl)-depg(ia,ic-1))/(dym(iac))
      else
        dddy(iac)=(depg(ia,ic+1)-depg(ia,ic-1))/(2*dym(iac))
      endif
      !if(iac>36800.and.iac<36810)print*,iac,dddx(iac),dddy(iac)
    end do !}
#ifdef FOR_YYZ
    if(mpi_id==0)call OutPutdep4yyz
#endif
  end SUBROUTINE settopog !}

end Module partition_mod
