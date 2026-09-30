#include "wavedeff.h"
Module output_loc_mod 
  use varcommon_mod
  use partition_mod
  use netcdf_mod
  IMPLICIT NONE
  !=======================
  private
  public::InitOutPutLoc,OpenHistLoc,OutPutHistLoc,closeoutpputloc,zftest
  INTEGER ::OutPutInited  =0
  integer(1) cstrbuf(16)
CONTAINS
  SUBROUTINE closeOutpputLoc(avh)
    type(AVHIST) avh
    if(avh%ncidwa>=0)call zf_close(avh%ncidwa)
    if(avh%ncidbv>=0)call zf_close(avh%ncidbv)
    avh%ncidwa=-1;avh%ncidbv=-1
  end SUBROUTINE closeOutpputLoc
  SUBROUTINE DeinitOutPut
    integer ko
    do ko=1,Navhists
      call closeOutpputLoc(avhists(ko))
    enddo
    OutPutInited=0
  end SUBROUTINE DeinitOutPut

  SUBROUTINE InitOutPutLoc(mode)!{
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
    call verifyLocDir
    cstrbuf=32 
    cstrbuf(1)=ichar('W');
    cstrbuf(2)=ichar('A');
    cstrbuf(3)=ichar('V');
    cstrbuf(4)=ichar('L');
    cstrbuf(5)=ichar('D');
    cstrbuf(6)=ichar('A');
    cstrbuf(7)=ichar(' ');

    do ko=1,Navhists
      avhists(ko)%ftimestr='tt'
      avhists(ko)%ncidwa=-1
      avhists(ko)%ncidbv=-1
    enddo
    NTYPE=4
    OutPutInited=1
  END SUBROUTINE InitOutPutLoc !}
  SUBROUTINE zftest
    integer fid
    integer*8 off,off1
    character*256 fn 
    write(fn,'(a,a,i6.6)')trim(pathLoc),"ttta.wav.h.wa.",mpi_id
    fid=zf_open(0,fn)
    off1=ZF_Size(fid)
    off=zf_write(fid,off1,mpi_id,4,1)
    off=zf_write(fid,off,mpi_id,4,1)
    DBGINF0,fid,off1,off
    call zf_close(fid)
  end SUBROUTINE zftest

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
  SUBROUTINE verifyDir(dir)
    character*(*) dir
    ! DBGO(0,*)'verifyDir',dir
    if(IsNewF(dir,0))then
      call zf_mkdir(dir)
    endif 
  end SUBROUTINE verifyDir

  SUBROUTINE verifyLocDir
    integer ip
    character*(256) dd,fn
    integer nsplit,nn

    if(mpi_id==0)then
      nsplit=1024
      nn=(mpi_npe+nsplit-1)/nsplit
      call verifyDir(pathLoc)
      do ip=0,nn
        write(dd,'(i8,"/")')ip
        dd=adjustl(dd)
        write(fn,'(a,a)')trim(pathLoc),trim(dd) 
        call verifyDir(fn)    
      enddo 
    endif
    call Barrier
  end SUBROUTINE verifyLocDir

  SUBROUTINE GetHistFn(fn,tag,vname,ftimestr)
    character*(*) fn,tag,vname,ftimestr
    character*(16) dd
    integer nn
    if(OutType==1)then
      write(fn,'(a,a,".wav.",a,".",a,".",a,".",i6.6)')trim(pathLoc),trim(caseid),trim(tag),trim(vname),trim(ftimestr),mpi_id
    else if(OutType==2)then
      write(dd,'(i8,"/")')mpi_id/1024
      dd=adjustl(dd)
      write(fn,'(a,a,a,".wav.",a,".",a,".",a,".",i6.6)')trim(pathLoc),trim(dd),trim(caseid),trim(tag),trim(vname),trim(ftimestr),mpi_id
    endif
  end SUBROUTINE GetHistFn
  subroutine OpenHistLoc(avh,newfile)
    type(AVHIST) avh
    INTEGER newfile
    DBGINF0
    call flush(6);call SerRun
    if(newfile/=0)then 
      call closeOutpputLoc(avh)
      DBGINF0
      call flush(6);call SerRun
      call GetHistFn(avh%wafn,avh%tag,"wa",avh%ftimestr)
      DBGINF0
      call flush(6);call SerRun
      call GetHistFn(avh%bvfn,avh%tag,"bv",avh%ftimestr)
      DBGINF0
      call flush(6);call SerRun
    endif
  END subroutine OpenHistLoc

  integer function VOpenHists(avh,bflg)
    type(AVHIST) avh
    integer bflg,ncid,nvar
    if(bflg==0)then
      VOpenHists=avh%ncidwa
      if(avh%ncidwa>=0)RETURN

      avh%ncidwa  =zf_open(0,avh%wafn)
      avh%waoff   =zf_size(avh%ncidwa)
      !print*,'open',avh%wafn;  call flush(6);call SerRun
      cstrbuf(7)=ichar('A')
      nvar=2
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,cstrbuf ,1,8)
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,mpi_npe ,4,1)
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,gnwpc   ,4,1)
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,nwpc    ,4,1)
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,nvar    ,4,1)
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,ieposc  ,4,nwpc*2)    
      VOpenHists=avh%ncidwa
    else if(bflg==1)then
      VOpenHists=avh%ncidbv
      if(avh%ncidbv>=0)RETURN
      avh%ncidbv  =zf_open(0,avh%bvfn)
      avh%bvoff   =zf_size(avh%ncidbv)
      cstrbuf(7)=ichar('A')
      nvar=2
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,cstrbuf ,1,8)
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,mpi_npe ,4,1)
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,gnwpc   ,4,1)
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,nwpc    ,4,1)
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,nvar    ,4,1)
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,ieposc  ,4,nwpc*2)    
      VOpenHists=avh%ncidbv
    endif
  end function VOpenHists

  subroutine OutPutHistLoc(timet,avh)  !{
    real(8) timet
    type(AVHIST) avh
    INTEGER ncid,vid,nvar
    INTEGER j,k
    character*(256) fn
    avh%istime=avh%istime+1
    IF(avh%wrwa>0)THEN
      OutDType=0;call starttimer(12+OutDType*NTYPE);
      ncid= VOpenHists(avh,0)
      !DBGINF;  call flush(6);call SerRun
      nvar=6
      cstrbuf(7)=ichar('B')
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,cstrbuf ,1,8)
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,gnwpc   ,4,1)
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,nwpc    ,4,1)
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,nvar    ,4,1)
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,runstate%cdatecur,4,1)
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,runstate%ctimecur,4,1)
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,timet   ,8,1)   
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,wxyo(1,:),4,nwpc)   
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,wxyo(2,:),4,nwpc)   
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,avh%h1_3    ,4,nwpc)    
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,avh%tpf     ,4,nwpc)    
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,avh%ape     ,4,nwpc)    
      avh%waoff=zf_write(avh%ncIDwa,avh%waoff,avh%aet     ,4,nwpc)    
      call  endtimer(12+OutDType*NTYPE)
    endif
    IF(avh%wrbv>0)THEN
      OutDType=1;call starttimer(12+OutDType*NTYPE);
      ncid= VOpenHists(avh,1);
      nvar=ndep
      cstrbuf(7)=ichar('C')
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,cstrbuf ,1,8)
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,gnwpc   ,4,1)
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,nwpc    ,4,1)
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,nvar    ,4,1)
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,runstate%cdatecur,4,1)
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,runstate%ctimecur,4,1)
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,timet   ,8,1)   
      avh%bvoff=zf_write(avh%ncIDbv,avh%bvoff,avh%bv      ,4,nwpc*ndep)   
      call endtimer(12+OutDType*NTYPE)
    end if
    ! call closeOutpput(avh);stop
  END subroutine OutPutHistLoc !}
END Module output_loc_mod
#ifdef NO_MPI
subroutine SerRun
end subroutine SerRun
#endif
