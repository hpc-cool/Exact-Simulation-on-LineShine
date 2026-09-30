#include "wavedeff.h"
Module output_nc_mod 
  use, intrinsic :: iso_c_binding 
#ifndef NO_MPI
  use wav_mpi_mod
#endif
  use netcdf_mod
  use varcommon_mod
  use partition_mod
  IMPLICIT NONE
  !=======================
  private
  public::InitOutPutNc,OpenHistNc,OutPutHistNc,closeoutpputNc ,Write_ncV2d 
  public::writegrd,writeee,outhss,testmpi
  INTEGER ::OutPutInited  =0
  real(4),allocatable::work2d(:,:)
  character*(*),parameter::ONAMEWindx='windx'
  character*(*),parameter::ONAMEwindy='windy'
  character*(*),parameter::ONAMEhs='HS'
  character*(*),parameter::ONAMEtz='TZ'
  character*(*),parameter::ONAMEtp='TP'
  character*(*),parameter::ONAMEth='TH'
  character*(*),parameter::ONAMEbv='BV'
  character*(*),parameter::ONAMEau11='AU11'
  character*(*),parameter::ONAMEau12='AU12'
  character*(*),parameter::ONAMEau22='AU22'
  character*(*),parameter::ONAMEau33='AU33'
  integer(1) cstrbuf(16)
  INTEGER DT(8)
  real(8) PntAb(4)
  real(8) phs,ptz,ptp,pth,pwv,pwd,pcu,pcd
  INTEGER ::DeltaNP,IOCP=0,NeedNext=0,MaxIOCP=0
  INTEGER ::ncidw=-1,irecw=-1,ilvlw=-1,ncvidw=-1
  character*16 ::ncvnamew
  INTEGER ::ncid,ncvids(30)
CONTAINS
  SUBROUTINE DeinitOutPut
    integer ko
    do ko=1,Navhists
      call closeOutpputNc(avhists(ko))
    enddo
  end SUBROUTINE DeinitOutPut
  SUBROUTINE InitOutPutNc(mode)!{
    INTEGER mode
    INTEGER ios,ih,iac,k
    real(8)dept
    integer ko
    IF(mode==-1)THEN
      call DeinitOutPut
      OutPutInited=0
      RETURN
    end if
    IF(OutPutInited/=0)RETURN
    IOCP=0
    if(outputNEP<=1)then
      DeltaNP=0;
      MaxIOCP=0
    else 
      DeltaNP=mpi_npe/outputNEP
      if(DeltaNP<1)DeltaNP=1
      MaxIOCP=outputNEP*DeltaNP
    endif
    ALLOCATE(work2d(gixl,giyl));
    work2d=nf_fill_real
    OutPutInited=1
    do ko=1,Navhists
      avhists(ko)%ftimestr='tt'
      avhists(ko)%ncidwa=-1
      avhists(ko)%ncidbv=-1
    enddo
  END SUBROUTINE InitOutPutNc !}
  SUBROUTINE closeOutpputNc(avh)
    type(AVHIST) avh
    call WriteOut
    if(avh%ncidwa>=0)call close_nc(avh%ncidwa)
    if(avh%ncidbv>=0)call close_nc(avh%ncidbv)
    avh%ncidwa=-1;avh%ncidbv=-1
  end SUBROUTINE closeOutpputNc
  SUBROUTINE NextIOP
    if(NeedNext==0)return 
    IF(IOCP>=0)then
      IOCP=IOCP+DeltaNP
      if(IOCP>=MaxIOCP)IOCP=0
    else
      IOCP=0
    endif
  end SUBROUTINE NextIOP
  SUBROUTINE WriteOut
    INTEGER nstart(4)
    if(ncidw<0)return
    if(irecw<=0.or.irecw>100)then
      ! return
    endif
    nstart=1;nstart(3)=ilvlw;
    if(ilvlw>=0)then
      call writenc(ncIDw,ncvnamew,work2d,irecw,nstart)
    else
      call writenc(ncIDw,ncvnamew,work2d,irecw)
    endif
    ncidw=-1
  end SUBROUTINE WriteOut
  SUBROUTINE SetDbuf(ncid,vname,irec,ilvl)
    INTEGER ncid,irec,ilvl
    character*(*)vname
    NeedNext=1
    if(mpi_id/=IOCP)return
    if(ncidw>=0)call WriteOut
    ncidw=ncid;  ncvnamew=vname
    irecw=irec;  ilvlw=ilvl
    !DBGO(0,'(i8,a,i6,a,6i6)')nTimeStep,'SetDbuf  ',irecw,ncvnamew,ilvlw,mpi_id
  end SUBROUTINE SetDbuf

  SUBROUTINE GetHistFn(fn,tag,vname,ftimestr)
    character*(*) fn,tag,vname,ftimestr
    call flush(6);call SerRun
    write(fn,'(a,a,".wav.",a,".",a,".",a,".nc")')trim(pathhist),trim(caseid),trim(tag),trim(vname),trim(ftimestr)
    call flush(6);call SerRun
  end SUBROUTINE GetHistFn
  logical function IsNewF(fn,newfile)
    character*(*) fn
    integer::newfile
    logical::fexist
    fexist=.false.
    if(newfile==0)then
      inquire(file=trim(fn), exist = fexist)
    endif
    IsNewF=.not.fexist
  end function IsNewF
  SUBROUTINE HistFileDef(fn,ncid,udep,dims)
    character*(*)fn
    integer dims(:),udep
    integer status
    INTEGER latID,lonID,depID,timeID,vlatID,vlonID,vdepID,vtimeID,ncid
    call open_nc(ncid,trim(fn),'create')
    call dimension_define(NcID, 'lon', gixl,lonID) !, 'lon', nf_float,vlonID)
    call dimension_define(NcID, 'lat', giyl,latID) !, 'lat', nf_float,vlatID)  
    if(udep/=0)then
      call dimension_define(NcID, 'dep', ndep,depID, 'dep', nf_float,vdepID)
    endif
    call dimension_define(NcID, 'time',NCUNLIM,timeID, 'time', nf_double,vtimeID)
    Dims(1)=lonID;  Dims(2)=latID;           
    if(gridtype<2)then
      call variable_define(ncid,'lon',nf_double,Dims(1:1))  
      call variable_define(ncid,'lat',nf_double,Dims(2:2))  
    else
      call variable_define(ncid,'lon',nf_double,Dims(1:2))  
      call variable_define(ncid,'lat',nf_double,Dims(1:2))        
    endif
    status=set_attribute(NcID,"units"    ,"degrees east","lon")
    status=set_attribute(NcID,"long_name", "coordiante longitude","lon")
    status=set_attribute(NcID,"units"    ,"degrees north","lat")  
    status=set_attribute(NcID,"long_name","coordiante latitude" ,"lat")
    status=set_attribute(NcID,"long_name","time","time")
    status=set_attribute(NcID,"units"    ,"days since 1950-01-01 00:00:00","time")      
    Dims(1)=lonID;  Dims(2)=latID;           
    if (udep==0) then
      Dims(3)=timeID
    else
      Dims(3)=depID;Dims(4)=timeID
    endif  
  end SUBROUTINE HistFileDef
  subroutine writexy(ncid)
    INTEGER ncid
    real(8),allocatable:: x1(:),y1(:)
    if(gridtype<2)then
      ALLOCATE(x1(gixl),y1(giyl))
      x1=xcord(:,1);
      y1=ycord(1,:);
      call writenc(NcID,'lon',x1)
      call writenc(NcID,'lat',y1)
      deallocate(x1,y1)
    else
      call writenc(NcID,'lon',xcord)
      call writenc(NcID,'lat',ycord)    
    endif
  end subroutine writexy
  subroutine OpenHistNc(avh,newfile)
    type(AVHIST) avh
    INTEGER newfile
    INTEGER Dims(4),mode
    INTEGER latID,lonID,depID,timeID,vlatID,vlonID,vdepID,vtimeID,ncid
    INTEGER ios
    logical alive
    INTEGER vid
    ! SDBGINF,' OpenHist '//avh%tag//trim(ftimestr),newfile,avh%LogWind
    !DBGINF0;    call flush(6);call SerRun
    if(newfile/=0)then
      call closeOutpputNc(avh)
      avh%istime=0
    endif
    !DBGINF0;    call flush(6);call SerRun

    IF(avh%wrwa/=0)THEN
      !OutDType=0;;
      !DBGINF0;    call flush(6);call SerRun
      call GetHistFn(avh%wafn,avh%tag,"wa",avh%ftimestr)
      !DBGINF0;    call flush(6);call SerRun
      if(mpi_id==0)then
        !DBGINF0;    call flush(6);call SerRun
        if(IsNewF(avh%wafn,newfile))then
          DBGO(0,*)'Create ',trim(avh%wafn)
          call flush(6)
          call HistFileDef(avh%wafn,ncid,0,dims)
          if(avh%LogWind/=0)call variable_define(ncid,ONAMEwindx,nf_float,Dims(1:3))
          if(avh%LogWind/=0)call variable_define(ncid,ONAMEwindy,nf_float,Dims(1:3))
          if(avh%Loghs  /=0)call variable_define(ncid,ONAMEhs   ,nf_float,Dims(1:3))
          if(avh%Logtz  /=0)call variable_define(ncid,ONAMEtz   ,nf_float,Dims(1:3))
          if(avh%Logtp  /=0)call variable_define(ncid,ONAMEtp   ,nf_float,Dims(1:3))
          if(avh%Logth  /=0)call variable_define(ncid,ONAMEth   ,nf_float,Dims(1:3))
          call end_define(ncid)
          call writexy(ncid)
          call close_nc(ncid)
        else
        end if
      endif

    end if
    IF(avh%wrbv/=0)THEN
      OutDType=1;;
      call GetHistFn(avh%bvfn,avh%tag,"bv",avh%ftimestr)
      if(mpi_id==0)then   
        if(IsNewF(avh%bvfn,newfile))then
          DBGO(0,*)'Create ',trim(avh%bvfn)
          call flush(6)
          call HistFileDef(avh%bvfn,ncid,1,dims)
          call variable_define(ncid,ONAMEbv,nf_float,Dims(1:4))
          call end_define(ncid)
          call writexy(ncid)
          call writenc(NcID,'dep',vdep(1:ndep))
          call close_nc(ncid)
        end if
      endif

    end if
  END subroutine OpenHistNc
  integer function VOpenHists(avh,bflg)
    type(AVHIST) avh
    integer bflg,ncid
    call NextIOP
    if(mpi_id/=IOCP)RETURN
    if(bflg==0)then
      VOpenHists=avh%ncidwa
      if(avh%ncidwa>=0)RETURN
      call open_nc(ncid,trim(avh%wafn),'write') 
      avh%ncidwa  =ncid
    else if(bflg==1)then
      VOpenHists=avh%ncidbv
      if(avh%ncidbv>=0)RETURN
      call open_nc(ncid,trim(avh%bvfn),'write') 
      avh%ncidbv  =ncid
    endif
    VOpenHists=ncid
  end function VOpenHists

  subroutine OutPutHistNc(timet,avh)  !{
    real(8) timet
    type(AVHIST) avh
    INTEGER ncid,vid
    INTEGER j,k
    character*(256) fn
    avh%istime=avh%istime+1
    IF(avh%wrwa>0)THEN
      OutDType=0;;
      ncid= VOpenHists(avh,0)
      if(mpi_id==IOCP)then
        call writenc(ncid,'time',timet,avh%istime)
      endif
      IF(avh%LogWind/=0)THEN
        call NCOUTV2D(avh,ONAMEWindx,wxyo(1,:),0)
        call NCOUTV2D(avh,ONAMEwindy,wxyo(2,:),0)
      end if  
      IF(avh%Loghs/=0)then
        call NCOUTV2d(avh,ONAMEhs,avh%h1_3,0)
      endif
      IF(avh%Logtp/=0)then
        call NCOUTV2d(avh,ONAMEtp,avh%tpf ,0)
      endif
      IF(avh%Logtz/=0)then
        call NCOUTV2d(avh,ONAMEtz,avh%ape ,0)
      endif
      IF(avh%Logth/=0)then
        call NCOUTV2d(avh,ONAMEth,avh%aet ,0) 
      endif
    endif
    call writegrd("hs_",avh%istime,avh%h1_3);
    IF(avh%wrbv>0)THEN
      OutDType=1;  
      ncid= VOpenHists(avh,1);
      if(mpi_id==IOCP)then
        call writenc(ncid,'time',timet,avh%istime)
      endif
      call NCOUTV3d(avh,ONAMEbv,avh%bv,ndep,1)  
    end if
    call WriteOut

    !  call DeinitOutPut;stop
  END subroutine OutPutHistNc !}

  subroutine NCOUTV3d(avh,vname,rv,nlvl,bflg)
    type(AVHIST) avh
    character*(*)vname
    INTEGER ncid,nlvl,bflg
    real(4) rv(:,:)
    INTEGER nstart(4)
    INTEGER il,it
    nstart=1
    ncid= VOpenHists(avh,bflg)     
    do il=1,nlvl
      call SetDbuf(ncid,vname,avh%istime,il)
      nstart(3)=il;

      work2d=nf_fill_real
      call GatherVar2d(work2d,rv(:,il),IOCP)

    enddo
  END subroutine NCOUTV3d !}

  subroutine NCOUTV2d(avh,vname,rv,bflg)
    type(AVHIST) avh
    character*(*)vname
    INTEGER ncid,bflg
    real(4) rv(:)
    real(4) vt
    INTEGER i,j

    ncid= VOpenHists(avh,bflg)
    call SetDbuf(ncid,vname,avh%istime,-1)
    work2d=nf_fill_real
    call GatherVar2d(work2d,rv,IOCP)
    !if(mpi_id==0)print*,'value of ',vname,work2d(1,giyl/4),rv(nwpc/4)

    !if(vname==ONAMEwindx)then
    ! DBGO0(0,*)'WINDX O check',runstate%timecur,mpi_id,maxval(rv),maxval(work2d)
    !endif

    !if(mpi_id==IOCP.and.vname==ONAMEhs)then
    ! vt=maxval(work2d)/2
    ! write(*,*)__FILE__,__LINE__,mpi_id,nTimeStep,vt
    ! do j=1,giyl
    !   do i=1,gixl
    !     if(work2d(i,j)>vt)then
    !       write(6,*)xcord(i),ycord(j),work2d(i,j)
    !     endif
    ! enddo
    ! enddo
    !endif
  END subroutine NCOUTV2d !}


  subroutine Write_ncV2d(fn,vname,rv,irec_)
    character*(*)fn,vname
    real(4) rv(:)
    INTEGER,optional:: irec_
    INTEGER ncid,vid,irec
    INTEGER dims(8)
    real*8 timet
    irec=1;
    if(present(irec_))irec=irec_
    if(mpi_id==0)then 
      call HistFileDef(fn,ncid,0,dims)
      call variable_define(ncid,vname,nf_float,Dims(1:3),vid)
      call end_define(ncid)
      call writenc(NcID,'lon',xcord)
      call writenc(NcID,'lat',ycord)
      work2d=nf_fill_real
    endif
    call GatherVar2d(work2d,rv,0)
    if(mpi_id==0)then
      timet=runstate%edaycur+runstate%timecur
      call writenc(ncid,'time',timet,irec)
      call writenc(ncID,vid,work2d,irec)
      call close_nc(ncID)
    endif
  END subroutine Write_ncV2d !}

#ifdef NO_MPI
#define OD "s/"
#else
#define OD "m/"
#endif
  subroutine testmpi(id)
    INTEGER,intent(in):: id
    INTEGER i
    real(4),allocatable:: rv(:)
    allocate(rv(0:nwpa));
    rv(0)=-2;
    do i=1,nwpc
      rv(i)=ieposc(1,i)+ieposc(2,i)*1000 +mpi_id*0.01
    enddo
    call writegrdi('ti',0,rv)
    do i=1,nwpc
      rv(i)=i +mpi_id*0.01
    enddo
    call writegrdi('td',0,rv)
    deallocate(rv);
  end subroutine 

  subroutine writegrdi(tag,id,rv)
    real(4),intent(in):: rv(:)
    character*(*),intent(in):: tag
    INTEGER,intent(in):: id
    INTEGER i,j,fid,ierr
    real(4):: nd,v0,v1,vd
    character*(256) fn
    !call WriteOut
    nd=0
    work2d=-1.
    call GatherVar2d(work2d,rv,0) 
    if(mpi_id==0)then
      write(fn,'(2a,i4.4,a)')OD,trim(tag),id,'.grd'
      open(newunit=fid, file=fn, status='replace', action='write', iostat=ierr)
      if (ierr /= 0) then
        print*,'open file:',trim(fn),'error'
        return
      endif
      v0=1e30;v1=-1e30
      do j=1,giyl
        do i=1,gixl
          vd=work2d(i,j)
          if(vd==nd)cycle
          v0=min(v0,vd);
          v1=max(v1,vd);
        enddo
      enddo

      write(fid,'(a)')"DSAA"
      write(fid,'(2i6)'  )gixl,giyl
      write(fid,'(2f8.2)')gxlon0,gxlon1
      write(fid,'(2f8.2)')gylat0,gylat1
      write(fid,'(2f10.2)')v0,v1
      do j=1,giyl
        write(fid,'(400f13.2)')work2d(:,j)
      enddo
      close(fid)
    endif
  END subroutine !}
  subroutine writegrd(tag,id,rv)
    real(4),intent(in):: rv(:)
    character*(*),intent(in):: tag
    INTEGER,intent(in):: id
    INTEGER i,j,fid,ierr
    real(4):: nd,v0,v1,vd
    character*(256) fn
    write(fn,'(2a,i4.4,a)')OD,trim(tag),id,'.grd'
    !call WriteOut
    nd=1.70141e+38
    work2d=nd
    call GatherVar2d(work2d,rv,0) 
    if(mpi_id==0)then
      open(newunit=fid, file=fn, status='replace', action='write', iostat=ierr)
      if (ierr /= 0) then
        print*,'open file:',trim(fn),'error'
        return
      endif
      v0=1e30;v1=-1e30
      do j=1,giyl
        do i=1,gixl
          vd=work2d(i,j)
          if(vd==nd)cycle
          v0=min(v0,vd);
          v1=max(v1,vd);
        enddo
      enddo

      write(fid,'(a)')"DSAA"
      write(fid,'(2i6)'  )gixl,giyl
      write(fid,'(2f8.2)')gxlon0,gxlon1
      write(fid,'(2f8.2)')gylat0,gylat1
      write(fid,'(2f10.2)')v0,v1
      do j=1,giyl
        write(fid,'(400g13.6)')work2d(:,j)
      enddo
      close(fid)
    endif
  END subroutine writegrd !}
  subroutine writeee(id,ee)
    integer id
    REALD ee(kl,jnthet)
    character*(256) fn
    INTEGER i,j,fid,ierr
    write(fn,'(2a,i4.4,i2.2,a)')OD,'e',id,mpi_id,'.txt'
    open(newunit=fid, file=fn, status='replace', action='write', iostat=ierr)
    do j=1,jnthet
      write(fid,'(32g13.4)')ee(:,j)
    enddo
    close(fid)
  end subroutine 
    
  subroutine writehs(id,rv)
    integer id
    real(4) rv(:)
    character*(256) fn
    INTEGER fid,ierr
    write(fn,'(2a,i4.4,i2.2,a)')OD,'v',id,mpi_id,'.txt'
    open(newunit=fid, file=fn, status='replace', action='write', iostat=ierr)
    write(fid,'(10g13.4)')rv(1:)
    close(fid)
  end subroutine 
  SUBROUTINE outhss(id,ee)
    REALD ee(kl,jnthet,0:*)
    integer id
    type (AVHIST),pointer::ph
    ph=>avhists(1)
    call mean1(ee,1 ,nwpc ,ph%ape,ph%tpf,ph%aet,ph%h1_3) 
    call writeee(id,ee(1,1,2));
    call writehs(id,ph%h1_3);
    call writegrd("v_",id,ph%h1_3);
  end SUBROUTINE 
END Module output_nc_mod
