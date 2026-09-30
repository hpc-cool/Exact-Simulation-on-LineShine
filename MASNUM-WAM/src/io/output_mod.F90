#include "wavedeff.h"
#define DBGINFZ   if(mpi_id==2)write(*,'(a,3i6,";",18i8)')__FILE__,__LINE__,mpi_id
Module output_mod
#ifndef NO_MPI
  use wav_mpi_mod
#endif
  use netcdf_mod
  use varcommon_mod
  use output_cal_mod
  use output_loc_mod
  use output_nc_mod
  use partition_mod
  IMPLICIT NONE
  !=======================
  private
  public::InitCheckOutPut,OpenHist,OutPutHist,CheckNeedOutPut,HistOutPutr
  public::InitOutPut,DeinitOutPut,InitMonitor,Monitor,ReadTopog,settopog
  !public::ACCUMEA
  public::zftest,OutType
  INTEGER ::OutPutInited  =0
  real(8) PntAb(4)
  real(8) phs,ptz,ptp,pth,pwv,pwd,pcu,pcd
  real(4),allocatable::work2d(:,:)
  INTEGER ::ncidw=-1,irecw=-1,ilvlw=-1
  character*16 ::ncvnamew
  INTEGER ::ncid
  INTEGER ::InitedLoge=0,InitedLog=0,InitedLoget=0
CONTAINS

  SUBROUTINE DeinitOutPut
    integer ko
    do ko=1,Navhists
      call closeOutpput(avhists(ko))
    enddo
  end SUBROUTINE DeinitOutPut
  SUBROUTINE InitOutPut(mode)!{
    INTEGER mode
    INTEGER ios,ih,iac,k
    real(8)dept
    integer ko
    if(OutType==0)then
      call InitOutPutNc(mode)
    else if(OutType==1)then
      call InitOutPutLoc(mode)
    else if(OutType==2)then
      call InitOutPutLoc(mode)
    endif
    IF(mode==-1)THEN
      call DeinitOutPut
      OutPutInited=0
      RETURN
    end if
    IF(OutPutInited/=0)RETURN
    do ko=1,Navhists
      avhists(ko)%wrwa=avhists(ko)%LogWind+avhists(ko)%Loghs+avhists(ko)%Logtz+avhists(ko)%Logtp+avhists(ko)%Logth
      avhists(ko)%wrbv=avhists(ko)%logbv
    enddo
    NTYPE=4
    OutPutInited=1
  END SUBROUTINE InitOutPut !}
  SUBROUTINE closeOutpput(avh)
    type(AVHIST) avh
    if(OutType==0)then
      call closeOutpputNc(avh)
    else if(OutType==1)then
      call closeOutpputLoc(avh)
    else if(OutType==2)then
      call closeOutpputLoc(avh)
    endif
  end SUBROUTINE closeOutpput

  subroutine OpenHist(avh,newfile)
    type(AVHIST) avh
    INTEGER newfile
    if(OutType==0)then
      call OpenHistNc(avh,newfile)
    else if(OutType==1)then
      call OpenHistLoc(avh,newfile)
    else if(OutType==2)then
      call OpenHistLoc(avh,newfile)
    endif
  end subroutine OpenHist
  subroutine OutPutHist(timet,avh)  !{
    real(8) timet !_,timet    
    type(AVHIST) avh
    !timet=idint(timet_*10000+0.5)/10000.
    if(OutType==0)then
      call OutPutHistnc(timet,avh)
    else if(OutType==1)then
      call OutPutHistLoc(timet,avh)
    else if(OutType==2)then
      call OutPutHistLoc(timet,avh)
    endif
  end subroutine OutPutHist
  ! =======================================================================
  integer function GussFileType(fn)
    character*(*) fn
    integer*1 buf(16)
    integer ifi,i,it
    integer(8) off
    DBGO0(0,*)'GussFileType '//trim(fn);call flush(6)
    ifi=zf_open(3,fn)
    GussFileType=0
    if(ifi<=0)return ;
    off=0;
    buf=0
    off=zf_read(ifi,off,buf,1,16)
    call zf_close(ifi)
    DBGO0(0,'(16i4)')buf
    GussFileType=3
    if(buf(1)==67.and.buf(2)==68.and.buf(3)==70.and.buf(4)<7)return ;!netcdf
    do i=1,16
      it=buf(i)
      if(it<32.and.it/=9.and.it/=10.and.it/=13)exit
    enddo
    GussFileType=2;if(i<16)return; !binary
    GussFileType=1;! Text
  end function GussFileType

#define DBGA call wav_mpi_barrier;write(*,'("DBGA",i4,":",i4.3)')__LINE__,mpi_id;  call flush(6);    
#define DBGB write(*,'("DBGA",i4,":",i4.3)')__LINE__,mpi_id;  call flush(6);    
SUBROUTINE ptop(tag,nx,ny,dept)
  character(256) fn
  integer nx,ny,i1,j1
  real(4)::dept(nx,ny)
  character(*) tag
  integer ,allocatable::zsg(:,:)
  write(fn,'(a,i2.2,".dat")')trim(tag),mpi_id
  open(11,file=fn,status='unknown')
  WRITE(11,*)nx,ny
  ALLOCATE(zsg(nx,ny));
  DO  j1=1,ny !{
    DO  i1=1,nx !{
      zsg(i1,j1)=1
      if(dept(i1,j1)<0.01)zsg(i1,j1)=0;
    enddo
  enddo
  DO  j1=ny,1,-1 !{
    WRITE(11,'(10000i1)')zsg(:,j1)
  end do !}
  DO  j1=ny,1,-1 !{
    WRITE(11,'(10000f7.1)')dept(:,j1)
  end do !}
  close(11)
  deallocate(zsg);
end SUBROUTINE 
  SUBROUTINE ReadTopog !{
    character fmtf*10
    real(4) ,allocatable::dept(:,:)
    integer ixltr,ixlt,iylt,ixo,iyo,i,j,i1,j1,ierr,ix,iy
    integer Gcirclet,stat,ftype
    real(4),allocatable:: xts(:),yts(:)
    real(4) xt0,xt1,yt0,yt1,x,y
    real(8) grdszxt,grdszyt,depv
    real(8) rgrdszxt,rgrdszyt

    integer ncid,dlatid,dlonid,vlatid,vlonid
    integer bctopomod
#ifndef NO_MPI
    type(mpipacket) pk
#endif
  ! SDBGINF,'gixl,giyl=',gixl,giyl
    
    if(DEPTHMOD>=10)then
      !nspg=1
      depg=DEPTHMOD;
      bctopomod=1;
    else
      if(mpi_id==0)then
        !write(*,*)'Read Topo File ',trim(topofn)," DEPTHMOD=",DEPTHMOD,"TOPOFTYPE=",TOPOFTYPE
        !DBGINF,DEPTHMOD,TOPOFTYPE
        ftype=TOPOFTYPE
        !ftype=GussFileType(topofn)
        ! ftype=3;
        if(ftype==1)then
          open(unit=1,file=topofn,status='old',IOSTAT=stat)
          if(stat/=0)then
            write(6,*)'Open Topo File error',trim(topofn)
            call wav_abort('Open Topo File error')
          endif
          call SkipLines(1)
          ! Read Grid Define
          read(1,*)xt0,xt1,yt0,yt1,grdszxt,grdszyt,Gcirclet
          rgrdszxt=1./grdszxt; rgrdszyt=1./grdszyt
          if(Gcircle==0 .and.(xt0 > gxlon0 .or. xt1<gxlon1))then
            DBGO(0,*)'Water Depth Data Area Not Match Setting',xt0 ,gxlon0 ,xt1,gxlon1
            call wav_abort('Water Depth Data Area Not Match Setting')
          end if
          if(yt0>gylat0 .or. yt1 < gylat1)then
            DBGO(0,*)'Water Depth Data Area Not Match Setting',yt0,gylat0 ,yt1 ,gylat1
            call wav_abort('Water Depth Data Area Not Match Setting')
          end if
          ixlt=(xt1-xt0)*rgrdszxt  ;iylt=(yt1-yt0)*rgrdszyt+1
          if(Gcirclet==0)ixlt=ixlt+1          
          ALLOCATE(xts(ixlt),yts(iylt));          
          ALLOCATE(dept(ixlt,iylt));dept=0;
          
          call SkipLines(1)
          read(1,*)((dept(i,j),i=1,ixltr),j=iylt,1,-1)
          close(1)
          do i=1,ixlt
            xts(i)=(i-1)*grdszxt+xt0
          enddo
          do i=1,iylt
            yts(i)=(i-1)*grdszyt+yt0
          enddo
        else if(ftype==3)then
          Gcirclet=0
          !DBGO(0,*)'open_nc '//trim(topofn)
          call open_nc(ncid,topofn,'read')
          ixlt=get_dimension_len(ncid,'lon')
          iylt=get_dimension_len(ncid,'lat')
          ALLOCATE(dept(ixlt,iylt));
          ALLOCATE(xts(ixlt),yts(iylt));
          call readnc(ncid,'lon',xts)
          call readnc(ncid,'lat',yts)
          ! SDBGINF,'xts';DBGO(8,'(20f8.3)')xts
          ! SDBGINF,'yts';DBGO(8,'(20f8.3)')yts
          call readnc(ncid,'etop5',dept)
          call close_nc(ncid)
          !write(101,*)"etop5",ixlt,iylt
          !write(101,*)xts
          !write(101,*)yts
          !write(101,*)dept
          dept=-dept
          do i1=2,ixlt
            if(xts(i1)<xts(i1-1))xts(i1)=xts(i1)+360
          enddo

          !if(mpi_id==0) call ptop("top1",ixlt,iylt,dept)
          !write(*,*)'Read Topo File END'
          xt0=xts(1);xt1=xts(ixlt)
          yt0=yts(1);yt1=yts(iylt)
          grdszxt=(xt1-xt0)/(ixlt-1)
          grdszyt=(yt1-yt0)/(iylt-1)
          rgrdszxt=1./grdszxt; rgrdszyt=1./grdszyt
          ixltr=ixlt
          if(abs(xt1-xt0-360)<1e-10)then
            Gcirclet=1;
            ixlt=ixlt-1
            xt0=xts(1);xt1=xts(ixlt)
            yt0=yts(1);yt1=yts(iylt)
          else if(abs(xt1-xt0+grdszxt-360)<1e-10)then
            Gcirclet=1;
          endif
          deallocate(xts,yts)
        else
          DBGO(0,*)'Topo File Type Error',ftype
          call wav_abort('Topo File Type Error')
        endif
        ! extract Data
        gnwpc=0
        DO  j1=1,iylt !{
        DO  i1=1,ixlt !{
          if(dept(i1,j1)>1e-5)then
            gnwpc=gnwpc+1
          endif
        end do !}
        end do !}
        !DBGO0(0,*)'Read Top InfoA',xt0,xt1,yt0,yt1,ixlt,iylt,gnwpc,shape(dept)
        !DBGO0(0,*)'Topotinfo',DEPTHMOD,grdszxt ,grdszx,grdszyt ,grdszy,yt0,yt1
        if(DEPTHMOD==0)then
          if(abs(idint(grdszx*rgrdszxt+0.5)*grdszxt-grdszx)>1e-5 .or.idint(grdszy*rgrdszyt+0.5)*grdszyt-grdszy>1e-5 )then
            DBGO(0,*)'Water Depth Data GridSize Not Match Setting',grdszxt,grdszyt, grdszx, grdszy
            call wav_abort('Water Depth Data GridSize Not Match Setting')
          end if
        endif
        bctopomod=0;
        if(grdszxt*grdszyt<grdszx*grdszy)bctopomod=1;
      endif
      !DBGO0(0,*)abs(idint(grdszx/grdszxt+0.5)*grdszxt-grdszx),idint(grdszy/grdszyt+0.5)*grdszyt-grdszy
        
#ifndef NO_MPI
      call initmpipacket(pk,10240)
      if(mpi_id==0)then
        call PackVars(0)
        call SetMpiPacketDSize(pk)
      endif
      
      call wav_mpi_bcast(pk);
      if(mpi_id/=0)then
        call PackVars(1)
      endif
      call InitMpiPacket(pk,-1)
#endif
      if(bctopomod==1)then
          if(mpi_id/=0)then
            ALLOCATE(dept(ixltr,iylt));
          endif
#ifndef NO_MPI
          if(mpi_id==0.or.mpi_id==mpi_npe-1)print*,mpi_id,'bcast dept begin'
          call wav_mpi_bcast(dept,  0,"AAAAAA"  )
          if(mpi_id==0.or.mpi_id==mpi_npe-1)            print*,mpi_id,'bcast dept end'
          !call ptop("topB",ixltr,iylt,dept)
#endif
      endif
      !if(mpi_id==0.or.mpi_id==mpi_npe-1)print*,mpi_id,'InterDep begin'
      if(mpi_id==0 .or. bctopomod==1)then
        if(DEPTHMOD==0)then
          DO  j1=1,giyl !{
            y=(j1-1)*grdszy+gylat0
            iy=idint((y-yt0)*rgrdszyt+1.5)
            !DBGO0(0,*)yt0,y,iy,yts(iy),(y-yt0)/grdszyt+1.5
            DO  i1=1,gixl !{
              x=(i1-1)*grdszx+gxlon0
              ix=idint((x-xt0)*rgrdszxt+1.5)
              !if(abs(x-xts(ix))>1e-5.or.abs(y-yts(iy))>1e-5)write(6,*)x,y,xts(ix),yts(iy)
              depv=dept(ix,iy)
              IF(depv>1.e-5)then
                if(depv<10.) then
                  depg(i1,j1)=10
                else if(depv>200.) then
                  depg(i1,j1)=200
                else
                  depg(i1,j1)=depv
                endif
              else
                depg(i1,j1)=0
              endif
            end do !}
          end do !}
          !if(mpi_id==0) call ptop("top2.dat",gixl,giyl,depg)
        else if(DEPTHMOD==1)then
          !call InterGrd(dept,depg,xts,ixlt,yts,iylt)
          !print*,"AAAAAAAAAAA",gxlon0,gylat0,grdszx,grdszy,xt0,yt0,1/rgrdszxt,1/rgrdszyt,rgrdszxt,rgrdszyt
          DO  j1=1,giyl !{
          DO  i1=1,gixl !{
            depv=InterDep(i1,j1)
            IF(depv>1.e-5)then
              if(depv<10.) then
                depg(i1,j1)=10
              else if(depv>200.) then
                depg(i1,j1)=200
              else
                depg(i1,j1)=depv
              endif
            else
              depg(i1,j1)=0
            endif
          end do !}
          end do !}
        endif
        !if(mpi_id==0) call ptop("topC",gixl,giyl,depg)
  
        DEALLOCATE(dept)
      endif
    endif
    !if(mpi_id==0.or.mpi_id==mpi_npe-1)print*,mpi_id,'InterDep end'
#ifndef NO_MPI
    if(bctopomod==0)then
      !if(mpi_id==0.or.mpi_id==mpi_npe-1)print*,mpi_id,'bcast depg begin'
      call wav_mpi_bcast(depg)
      !call ptop("topA",gixl,giyl,depg)
      !if(mpi_id==0.or.mpi_id==mpi_npe-1)print*,mpi_id,'bcast depg end'
    endif
#endif
    !if(mpi_id==0.or.mpi_id==mpi_npe-1)print*,mpi_id,'ResetNspg begin'
    call ResetNspg(nspg)
    !if(mpi_id==0.or.mpi_id==mpi_npe-1)print*,mpi_id,'ResetNspg end'
    

    if(1==0 .and. mpi_id==0)then
      open(11,file="nsg.dat")
      WRITE(11,*)gixl,giyl
        DO  j=giyl,1,-1
          WRITE(11,'(10000i1)')(nspg(i,j),i=1,gixl)
        end do !}
      close(11)
    endif
    CONTAINS
    integer function InterDep(i,j)
      integer::i,j, i1,j1,i11,j11
      real*8 x,y,a,b
      real pxt0,pyt0
      pxt0=xt0
      pyt0=yt0
      y=(j-1)*grdszy+gylat0
      j1=(y-yt0)*rgrdszyt+1 ;
      b=(y-yt0)*rgrdszyt-(j1-1)

      x=(i-1)*grdszx+gxlon0
      if(x<xt0)then
        x=x+360.
      else if(x>xt1)then
        x=x-360.
      endif 
      i1=(x-xt0)*rgrdszxt+1 ;
      a=(x-xt0)*rgrdszxt-(i1-1)
      InterDep=dept(i1  ,j1  )*(1-a)*(1-b)+ &
               dept(i1+1,j1  )*(  a)*(1-b)+ &
               dept(i1  ,j1+1)*(1-a)*(  b)+ &
               dept(i1+1,j1+1)*(  a)*(  b)
    end function InterDep
#ifndef NO_MPI
    SUBROUTINE PackVars(iunpack)
      integer iunpack
      call wav_mpi_pack(pk,bctopomod,iunpack );
      call wav_mpi_pack(pk,DEPTHMOD ,iunpack );
      call wav_mpi_pack(pk,grdszxt  ,iunpack );
      call wav_mpi_pack(pk,grdszyt  ,iunpack );
      call wav_mpi_pack(pk,rgrdszxt  ,iunpack );
      call wav_mpi_pack(pk,rgrdszyt  ,iunpack );
      call wav_mpi_pack(pk,xt0      ,iunpack );
      call wav_mpi_pack(pk,xt1      ,iunpack );
      call wav_mpi_pack(pk,yt0      ,iunpack );
      call wav_mpi_pack(pk,yt1      ,iunpack );
      call wav_mpi_pack(pk,ixltr    ,iunpack );
      call wav_mpi_pack(pk,ixlt     ,iunpack );
      call wav_mpi_pack(pk,iylt     ,iunpack );
    end SUBROUTINE PackVars
#endif    
  end SUBROUTINE ReadTopog
#ifdef FOR_YYZ
  SUBROUTINE OutPutdep4yyz
    integer ncid,latID,lonID,vlonID,vlatID,zyyzid,vzyyzid,timeID,vtimeID
    integer dims(8)
    integer nlon,nlat,nzyyz,ic
    real(4),allocatable::lon(:),lat(:),zyyz(:),depyyz(:,:),nspyyz(:,:),rslat(:)
    nlon=gixl;nlat=giyl;nzyyz=ndep
    !if(Gcircle/=0)nlon=nlon+1
    allocate(lon(nlon),lat(nlat),rslat(nlat),zyyz(nzyyz),depyyz(nlon,nlat),nspyyz(nlon,nlat))
    lat=ycord(1,:);
    zyyz=-vdep(1:ndep);
    lon(1:gixl)=xcord(:,1);
    depyyz(1:gixl,:)=depg
    nspyyz(1:gixl,:)=nspg
    if(Gcircle/=0)then
      !lon(nlon)=lon(1)
      !depyyz(nlon,:)=depyyz(1,:)
      !nspyyz(nlon,:)=nspyyz(1,:)
    endif
    do ic=1,nlat
      rslat(ic)=rs*cosd(lat(ic))
    enddo
    call open_nc(ncid,trim(pathinit)//"wamyyz.nc",'create')
    call dimension_define(NcID,'lon'   , nlon,lonID, 'lon', nf_float,vlonID)
    call dimension_define(NcID,'lat'   , nlat,latID, 'lat', nf_float,vlatID)
    call dimension_define(NcID,'zyyz'  , nzyyz,zyyzid, 'zyyz', nf_float,vzyyzid)
    call variable_define (ncid,'rslat' ,nf_float,[latID])
    call variable_define (ncid,'depyyz',nf_float,[lonID,latID])
    call variable_define (ncid,'nspyyz',nf_int,[lonID,latID])
    call end_define(ncid)
    call writenc(ncid,'lon'   ,lon   )
    call writenc(ncid,'lat'   ,lat   )
    call writenc(ncid,'zyyz'  ,zyyz  )
    call writenc(ncid,'rslat' ,rslat )
    call writenc(ncid,'depyyz',depyyz)
    call writenc(ncid,'nspyyz',nspyyz)
    call close_nc(ncid)
  end SUBROUTINE OutPutdep4yyz
#endif  

  subroutine Monitor(ea,iacp_,iid)
    REALD,intent(in) :: ea(kl,jnthet,0:nwpc)
    integer iacp_,iid
    real(8) ::wxt,wyt,apet,tpft,aett,h1_3t ! aet  :波向   tpf  :谱峰周期    h1_3 :波高    ape  :跨零中期
    integer i,iit,iac,ia,ic
    if(iacp_>0)then
      if(iid<0)then
        WRITE(116,*)ea(:,:,iac)
      endif
      iac=iacp_
      call mean1t(ea(1,1,iac),iac,apet,tpft,aett,h1_3t)
      wxt=wxy(1,iac);wyt=wxy(2,iac);
      write(116,'(i6,3i6,I5,i7,7f12.6)')mpi_id, ia,ic,iac,runstate%ndays,runstate%ctimecur,h1_3t,apet,tpft,aett,wxt,wyt
      call flush(116)
    else
      do i=1,nmonitor
        ia =iap (nmonitor)
        ic =icp (nmonitor)
        iac=iacp(nmonitor)
        call mean1t(ea(1,1,iac),iac,apet,tpft,aett,h1_3t)
        wxt=wxy(1,iac);wyt=wxy(2,iac);
      if(iid<0)then
        WRITE(116,'(32f8.3)')ea(:,:,iac)
      endif
        !write(6,'(i6,3i6,I10,I5,i7.6,7f12.6)')mpi_id, ia,ic,iac,runstate%nTimeStep,runstate%ndays,runstate%ctimecur,h1_3t,apet,tpft,aett,wxt,wyt
        write(116,'(i6,3i6,I10,I5,i7.6,7f12.6)')mpi_id, ia,ic,iac,runstate%nTimeStep,runstate%ndays,runstate%ctimecur,h1_3t,apet,tpft,aett,wxt,wyt
        call flush(116)
      enddo
    endif
  end subroutine Monitor
  subroutine InitMonitor
    integer i,ia,ic,iac
    nmonitor=0
    iacp=0
    do i=1,size(pla)
      if(pla(i)>360.)exit
      gnmonitor=i
      ia=(plo(i)-gxlon0)/GRDSZX+1
      ic=(pla(i)-gylat0)/GRDSZY+1
      if(Gcircle/=0)then
        if(ia<0)ia=ia+gixl
        if(ia>gixl)ia=ia-gixl
      endif
      iac=pos2id(ia,ic)
      !write(*,'(a,2f8.3,12i6)')'working point site:',plo(i),pla(i),mpi_id,ia,ic,nspg(ia,ic),iac,rectc
      
      if(iac/=0)then
        !write(*,'(a,i5,a,8i6)')'working point site: AA',iac,'mpi_id:',mpi_id,nsp(iac),nwpc,nwpa
        if(iac>nwpc.or.iac<0)then
          DBGO(0,*)'pos2id Error', ia,ic,iac
          iac=0
        else
          if(iac/=0)then
            if(nsp(iac)<=0)iac=0
          endif
        endif
        if(iac/=0)then
          nmonitor=nmonitor+1
          iap(nmonitor)=ia;
          icp(nmonitor)=ic
          iacp(nmonitor)=iac
          write(*,'(a,i8,a,i6)')'working point site: ',iac,' of mpi_id:',mpi_id
          write(*,'(a,2f9.3,3i8,a,i6)')'working point site: ',plo(i),pla(i),ia,ic,iac,' @ mpi_id:',mpi_id
          open(116,file="monitor.dat",status="unknown")
          write(6,*)mpi_id,'InitMonitor Open  monitor.dat OK'
          write(116,'(a,2f9.3,3i5,a,i6)')'working point site: ',plo(i),pla(i),ia,ic,iac,'@ mpi_id:',mpi_id
          call flush(116)
        endif
      endif
    enddo
    return
  end subroutine InitMonitor
  subroutine InitCheckOutPut(mode)
    integer ko
    integer mode
    type(AVHIST),pointer:: avh

    !DBGO0(0,*)' InitCheckOutPut ',Navhists,mode
    if(mode==0)then
      do ko=1,Navhists
        avh=>avhists(ko)
        if(avh%mean/=0)then
          if(.not. associated(avhists(ko)%ea) )then
            ALLOCATE(avhists(ko)%ea(kl,jnthet,0:nwpc));
          endif
          avhists(ko)%ea=0;
        endif
#define VALLOC(ln,vn,S)\
        IF(avhists(ko)%ln/=0)then;\
          if(.not.associated(avhists(ko)%vn))ALLOCATE(avhists(ko)%vn S );\
        endif
        VALLOC(loghs,h1_3,(0:nwpc))
        VALLOC(logtz,ape ,(0:nwpc))
        VALLOC(logtp,tpf ,(0:nwpc))
        VALLOC(logth,aet ,(0:nwpc))
        VALLOC(logbv,bv  ,(0:nwpc,ndep))
        avhists(ko)%nea=0;
        avhists(ko)%nextTime=0;
        avhists(ko)%itpf=0
        ! DBGO0(0,*)"InitCheckOutPut ",ko,runstate%timecur,avhists(ko)%nextTime
      enddo
    else
      do ko=1,Navhists
        call NextTime(avhists(ko)%nextTime,avhists(ko)%OPTION,avhists(ko)%NN,avhists(ko)%mean)
        ! DBGO0(0,*)"InitCheckOutPut ",ko,timecur,avhists(ko)%nextTime
      enddo
    endif
  end subroutine InitCheckOutPut

  subroutine CheckNeedOutPut(flag) 
    integer flag
    integer ko,nt
    real(8) ::timet
    type (AVHIST),pointer::ph
    flag=0
    timet=runstate%edaycur+runstate%timecur+deltts/10/86400.
    if(runstate%NextOutStep<=runstate%nTimeStep)runstate%NextOutStep=runstate%nTimeStep+1000000;
    do ko=1,Navhists      
      ph=>avhists(ko)
      if(ph%itype==0)then
        if(runstate%hist_now/=0)flag=flag+1
        ph%hist_ect=1
      else if(ph%nextTime<=timet)then
        ph%hist_ect=1
        flag=flag+1
      endif     
      if(ph%hist_ect/=0)then
        call NextTime(ph%nextTime,ph%OPTION,ph%NN,ph%mean)
      endif
      nt=runstate%nTimeStep+((ph%nextTime-timet)*86400-deltts*0.1)/deltts+1;
      if(nt<runstate%NextOutStep)runstate%NextOutStep=nt;
    enddo   
  end subroutine CheckNeedOutPut
  
  subroutine HistOutPutr(avh)
    type(AVHIST) avh
    integer ::newfile,ntpf
    newfile=0
    avh%nextTime=runstate%edaycur+runstate%timecur;
    ! AVHIST_NTPFS :每文件输出次数,>0: n Times;0: 1 day;-1: 1 month;-2: 1 Season;-3:half year:<=-4: -n-4 year '
    ntpf=avh%NTPF
    if(ntpf>0)then
      if(avh%itpf>=ntpf)then
        write(avh%ftimestr,'(i8.8,"-",i6.6)')runstate%cdatecur,runstate%ctimecur
        newfile=1
      endif
    else if(ntpf==0)then ! day
        if(runstate%cdatecur /=avh%oldid)then
          write(avh%ftimestr,'(i8.8)')runstate%cdatecur
          avh%oldid=runstate%cdatecur
          newfile=1
        endif
    else if(ntpf>-4)then ! month
        if(runstate%cdatecur/100 /=avh%oldid)then
          write(avh%ftimestr,'(i6.6)')runstate%cdatecur/100
          avh%oldid=runstate%cdatecur/100
          newfile=1
        endif
    else if(ntpf<=-4)then ! years
        if(runstate%cdatecur/10000-avh%oldid>-ntpf-3)then
          write(avh%ftimestr,'(i4.4)')runstate%cdatecur/10000
          avh%oldid=runstate%cdatecur/10000
          newfile=1
        endif
    endif
    !DBGO0(5,'("O",$)')
    !DBGO0(1,'(" OutPutHist  ",4i4,2i9," ",a)')newfile,ntpf,avh%itpf,runstate%cdatecur/10000-avh%oldid,runstate%cdatecur,avh%oldid,trim(avh%ftimestr)
    call OpenHist(avh,newfile)
    if(newfile/=0)then
      avh%itpf=0
    endif
    avh%itpf=avh%itpf+1
    if(avh%mean/=0)then
      call OutPutHist(runstate%edaycur+runstate%timecur,avh) ! OutPut
    else
      call OutPutHist(runstate%edaycur+runstate%timecur,avh) ! OutPut
    endif
    call NextTime(avh%nextTime,avh%OPTION,avh%NN,avh%mean)
    call endtimer(10);call starttimer(10)  
  end subroutine HistOutPutr
END Module output_mod
subroutine HistOutPut(avh)
  use varcommon_mod
  use output_mod
  type(AVHIST) avh
  call HistOutPutr(avh)
end subroutine HistOutPut
