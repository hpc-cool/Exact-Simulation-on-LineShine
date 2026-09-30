#include "wavedeff.h"
! modify varcommon wavemdl surport restart Test
Module restart_mod
#ifndef NO_MPI
use wav_mpi_mod
#endif
use varcommon_mod
!use output_mod
use partition_mod
use windin_mod
use netcdf_mod
IMPLICIT NONE
   integer,parameter :: NB_RN = kind(1.0D0)              ! native real
   integer,parameter :: NB_IN = kind(1)                ! native integer
    integer ::restartinited=0,firstrestart=1
    integer irtag
    character*120 :: RestHead
    integer(8),allocatable::ioffs(:)
    integer(8)::ioffcur=-1
    integer ::aoneas(IMAXAVHIST)
    real(8)::aotimes(IMAXAVHIST)
CONTAINS
! use ZF_open ZF_READ ZF_WRITE writed in C language
! ZF_READ (ifile,off,var,NByte,nItem);
! ZF_WRITE(ifile,off,var,NByte,nItem);
SUBROUTINE InitRestart
    restartinited=1
end SUBROUTINE InitRestart
	SUBROUTINE CheckNeedRestart(flag)
		integer flag
		flag=0
    if(rest_N==0 .and. runstate%rest_now/=0)then
      flag=1
    else if(resttime<=runstate%timecur)then
      flag=1;
    endif
		
	end SUBROUTINE CheckNeedRestart
  subroutine CheckRestart(eec)
    REALD   eec(kl,jnthet,0:nwpa)
    integer ::rest
    rest=0;
    if(rest_N==0 .and. runstate%rest_now/=0)then
      rest=1
    else if(resttime<=runstate%timecur)then
      rest=1;call NextTime(resttime,rest_option,rest_N,0)
    endif
    if(rest/=0)then
      call Write_restart(eec)
    end if
  end subroutine CheckRestart

SUBROUTINE ErrStop(msg)
  character*(*) msg
  write(*,*)msg
  stop
end SUBROUTINE ErrStop
  SUBROUTINE TestRestart(iTestRestart,eec)
    REALD   eec(kl,jnthet,0:nwpa)
    integer iTestRestart
    integer intarr(16)
    integer iac ,ierr
    if(iTestRestart>1)then
      intarr( 1)=   Navhists
      intarr( 2)=   ndep
      intarr( 3)=   mpi_npe
      do iac=1,nwpc
        eec(:,:,iac) =iac*10+1
        wind0x (iac) =iac*10+2
        wind0y (iac) =iac*10+3
        windx  (iac) =iac*10+4
        windy  (iac) =iac*10+5
        !wiv    (iac) =iac*10+6
        !winc   (iac) =iac*10+7
      enddo
      call write_restart(eec)
    endif
    eec=0;wind0x=0;wind0y=0;windx=0;windy=0; 
    call read_restart(eec)
    ierr=0
    if(iTestRestart>1)then
      if(intarr( 1)/=   Navhists)then ;write(*,*)mpi_id,"Verify Restart Error Navhists  ";ierr= 1;endif
      if(intarr( 2)/=   ndep    )then ;write(*,*)mpi_id,"Verify Restart Error ndep      ";ierr= 2;endif
      if(intarr( 3)/=   mpi_npe )then ;write(*,*)mpi_id,"Verify Restart Error mpi_npe   ";ierr= 3;endif
    endif
    do iac=1,nwpc
      if(eec(1,1,iac) /=iac*10+1)then ;write(*,*)mpi_id,iac,"Verify Restart Error eec   ",eec(1,1,iac),iac*10+1;ierr=iac* 5;endif
      if(wind0x (iac) /=iac*10+2)then ;write(*,*)mpi_id,iac,"Verify Restart Error wind0x",wind0x (iac),iac*10+2;ierr=iac* 6;endif
      if(wind0y (iac) /=iac*10+3)then ;write(*,*)mpi_id,iac,"Verify Restart Error wind0y",wind0y (iac),iac*10+3;ierr=iac* 7;endif
      if(windx  (iac) /=iac*10+4)then ;write(*,*)mpi_id,iac,"Verify Restart Error windx ",windx  (iac),iac*10+4;ierr=iac* 8;endif
      if(windy  (iac) /=iac*10+5)then ;write(*,*)mpi_id,iac,"Verify Restart Error windy ",windy  (iac),iac*10+5;ierr=iac* 9;endif
      !if(wiv    (iac) /=iac*10+6)then ;write(*,*)mpi_id,iac,"Verify Restart Error wiv   ",wiv    (iac),iac*10+6;ierr=iac*10;endif
      !if(winc   (iac) /=iac*10+7)then ;write(*,*)mpi_id,iac,"Verify Restart Error winc  ",winc   (iac),iac*10+7;ierr=iac*11;endif
    enddo
    if(ierr/=0)then
      write(6,*)mpi_id,"check restart Err"
    else
      write(6,*)mpi_id,"check restart ErrOK"
    endif
  end SUBROUTINE TestRestart
  SUBROUTINE CalRestOff(off,nwpcst,mpi_npet,ioffst)
     integer(8) off
     integer mpi_npet
     integer nwpcst(0:mpi_npet)
     integer(8) ioffst(:)
     integer(8) sizet
     integer i,ko
     ioffst=0
     ioffst(1)=off ! Note nwpcs is (1:mpi_npe);ioffs is (0:mpi_npe-1)
     do i=2,mpi_npet
       sizet=(nwpcst(i-1))*NB_RN*kl*jnthet !e
       sizet=sizet+nwpcst(i-1)*6*NB_RN
       do ko=1,Navhists
         if(aoneas(ko)/=0)then
           sizet=sizet+(nwpcst(i-1))*NB_RN*kl*jnthet
         end if
       end do
       DBGO(1,*)mpi_id,sizet,'ioffst',nwpcst(i-1),ioffst
       ioffst(i)=ioffst(i-1)+sizet
     end do
     DBGO(1,*)mpi_id,'ioffst',ioffst
! 368
! 23724368
! 29654768 5930400
! 35587568
   end SUBROUTINE CalRestOff
   SUBROUTINE CalRestOffc(off)
     integer(8) off
     if(mpi_id==0)then
       if(ioffcur<=0)then
         allocate(ioffs(mpi_npe))
         call CalRestOff(off,nwpcs,mpi_npe,ioffs)
       end if
     endif
#ifndef NO_MPI
     if(ioffcur/=0)call wav_mpi_scatter(ioffs,ioffcur);
#endif
   end SUBROUTINE CalRestOffc

   SUBROUTINE write_restart(eec)
     REALD   eec(kl,jnthet,0:nwpa)
     integer ifile
     integer(8)::off
     integer i,ko
     integer intarr(16)
     character*256 :: RestFile
     write(RestFile,'(a,".wav.r.",i8.8,"-",i8.8,".restart")')trim(caseid),runstate%cdatecur,runstate%ctimecur
     if(mpi_id == 0) then
       DBGO(0,*) "Now create restart file:",trim(RestFile)
     end if
     if(mpi_id==0)then
       if(firstrestart==1)then
         firstrestart=0
         write(RestHead,'(a,".wav.r.",i8.8,"-",i8.8,".rhdr")')trim(caseid),runstate%cdatecur,runstate%ctimecur
         DBGO(0,*) "Now create restart Head file:",trim(RestHead)
         call SYSTEM_CLOCK(count=irtag)
         ifile=ZF_open(0,trim(pathrest)//RestHead)
         off=0;
         off=ZF_WRITE(ifile,off,irtag     ,NB_IN,1);
         off=ZF_WRITE(ifile,off,mpi_npe   ,NB_IN,1);
         off=ZF_WRITE(ifile,off,timeinit  ,8,1);
         off=ZF_WRITE(ifile,off,ndep      ,NB_IN,1);
         off=ZF_WRITE(ifile,off,vdep      ,4,ndep);
         off=ZF_WRITE(ifile,off,nwpcs     ,NB_IN,mpi_npe+1)
         off=ZF_WRITE(ifile,off,rects     ,NB_IN,4*mpi_npe)
         !off=ZF_WRITE(ifile,off,ieposg    ,2,2*nwpcs(0))
         call ZF_close(ifile);
       end if
       intarr=0
       intarr(1)=Navhists;intarr(2)=ndep;intarr(3)=mpi_npe;
       do ko=1,Navhists
           aoneas (ko)=avhists(ko)%nea;
           aotimes(ko)=avhists(ko)%nextTime
       enddo
       ifile=ZF_open(0,trim(pathrest)//RestFile)
       off=0
       off=ZF_WRITE(ifile,off,irtag   ,NB_IN,1);
       off=ZF_WRITE(ifile,off,timeinit ,8,1);
       off=ZF_WRITE(ifile,off,runstate%timecur ,8,1);
       off=ZF_WRITE(ifile,off,RestHead,1,len(RestHead));
       off=ZF_WRITE(ifile,off,intarr  ,NB_IN,16);
       off=ZF_WRITE(ifile,off,aoneas  ,NB_IN,Navhists);
       off=ZF_WRITE(ifile,off,aotimes ,8,Navhists);
       off=ZF_WRITE(ifile,off,vdep    ,4,ndep);
     end if

     call CalRestOffc(off)
     ! call WaitPrev
     if(mpi_id/=0)then
         ifile=ZF_open(0,trim(pathrest)//RestFile)
     end if
     off=ioffcur ;
     off=ZF_WRITE(ifile,off,eec(1,1,1)   ,NB_RN,nwpc*kl*jnthet);
     off=ZF_WRITE(ifile,off,wind0x,4,nwpc);
     off=ZF_WRITE(ifile,off,wind0y,4,nwpc);
     off=ZF_WRITE(ifile,off,windx ,4,nwpc);
     off=ZF_WRITE(ifile,off,windy ,4,nwpc);
     !off=ZF_WRITE(ifile,off,wiv   ,4,nwpc);
     !off=ZF_WRITE(ifile,off,winc  ,4,nwpc);
     !off=ZF_WRITE(ifile,off,bv    ,4,nwpc*ndep);
     do ko=1,Navhists
         if(aoneas(ko)/=0)then
            off=ZF_WRITE(ifile,off,avhists(ko)%ea,NB_RN,nwpc*kl*jnthet);
         end if
     enddo
     DBGO(1,*)mpi_id,"complete restart data out!"
     call ZF_close(ifile);
     if(mpi_id == 0) then
         open(32,file=trim(pathrest)//"../rpointer.wav",form='formatted',action='write')
         write(32,'(a)') trim(RestFile)
         write(32,*) trim(caseid)
         write(32,'(a)') trim(RestHead)
         close(32)
         DBGO(0,*) "create restart file End"
     end if
     RETURN
 END SUBROUTINE write_restart
! !******************************************************************************
 SUBROUTINE read_restart(eec)
    	REALD   eec(kl,jnthet,0:nwpa)
     integer ifile,ifileh,ko,i,j
     integer(8)::off,offt
     character(256) ::RestFile
     integer intarr(16)
     integer(NB_IN)::irtagt,ndept,mpi_npet
     real(4) vdept(160)
     real(4),allocatable::v2d(:),vbv(:,:)
     REALD,allocatable::vee(:,:,:)
     INTEGER(NB_IN),allocatable::nwpcst(:),rectst(:)
     INTEGER(8),allocatable::ioffst(:)
     INTEGER(2),allocatable::ieposgt(:,:)
     INTEGER nipos,samevdep,icwp,nicwp,ncwp,nepos;
     integer ::ix,iy,ia,ic,icvt,icv,ixlt,ixlrt,iac,iact,ixot,iyot,ixod,iyod,ixst,iyst,ixet,iyet,iip
     integer ::gnwpct,maxnwpc,ipoff
     ! Read rpointer
      intarr=0
     if(mpi_id == 0) then
       open(32,file=trim(pathrest)//"../rpointer.wav",form='formatted')
       read(32,'(a)') RestFile
       close(32)
     end if
#ifndef NO_MPI
     call wav_mpi_bcast(RestFile)
#endif
     if(mpi_id == 0) then
         DBGO(0,*) "Now Read restart file:",trim(RestFile)
     end if
     intarr=0
     if(mpi_id==0)then
         ifile=ZF_open(1,trim(pathrest)//RestFile)
         off=0
         off=ZF_READ(ifile,off,irtag  ,NB_IN,1);
         off=ZF_READ(ifile,off,timeinit ,8,1);
         off=ZF_READ(ifile,off,runstate%timecur ,8,1);
         off=ZF_READ(ifile,off,RestHead,1,len(RestHead));
         off=ZF_READ(ifile,off,intarr  ,NB_IN,16);
         Navhists   =intarr( 1);
         off=ZF_READ(ifile,off,aoneas  ,NB_IN,Navhists);
         off=ZF_READ(ifile,off,aotimes ,8,Navhists);
         off=ZF_READ(ifile,off,vdept   ,4,intarr( 2));
         samevdep=0
         if(ndep==intarr( 2))then
             samevdep=1
             do i=1,ndep
                 if(abs(vdep(i)-vdept(i))>1e10)then
                     samevdep=0;exit
                 endif
             enddo
         endif
     end if
#ifndef NO_MPI
     call wav_mpi_bcast(irtag  );
     call wav_mpi_bcast(runstate%timecur);
     call wav_mpi_bcast(intarr );
     call wav_mpi_bcast(aoneas );
     call wav_mpi_bcast(aotimes);
#endif
     Navhists   =intarr( 1);
     if(intarr(3)==mpi_npe)then
       call CalRestOffc(off)
       if(mpi_id/=0)then
           ifile=ZF_open(0,trim(pathrest)//RestFile)
       end if
       off=ioffcur
       off=ZF_READ(ifile,off,eec(1,1,1),NB_RN,nwpc*kl*jnthet);
       off=ZF_READ(ifile,off,wind0x,4,nwpc);
       off=ZF_READ(ifile,off,wind0y,4,nwpc);
       off=ZF_READ(ifile,off,windx ,4,nwpc);
       off=ZF_READ(ifile,off,windy ,4,nwpc);
       !off=ZF_READ(ifile,off,wiv   ,4,nwpc);
       !off=ZF_READ(ifile,off,winc  ,4,nwpc);
       do ko=1,Navhists
           if(aoneas(ko)/=0)then
              off=ZF_READ(ifile,off,avhists(ko)%ea,NB_RN,nwpc*kl*jnthet);
           end if
       enddo
       if(mpi_id.eq.mpi_npe-1)then
           print*,mpi_id,"complete restart data out!"
       end if
       call ZF_close(ifile);
     else
       DBGO(0,*) 'ERR NPE Changed not surport now ',mpi_id,mpi_npe,intarr
       stop
       offt=off
       if(mpi_id==0)then
           ifileh=ZF_open(0,trim(pathrest)//RestHead)
           off=0;
           off=ZF_READ(ifileh,off,irtagt,NB_IN,1);
           off=ZF_READ(ifileh,off,mpi_npet,NB_IN,1);
           off=off+8 !skip timeinit
           off=ZF_READ(ifileh,off,ndept,NB_IN,1);
           off=ZF_READ(ifileh,off,vdept,4,ndept);

           ALLOCATE(nwpcst(mpi_npet),rectst(4*mpi_npet))
           off=ZF_READ(ifileh,off,nwpcst,NB_IN,mpi_npet+1)
           off=ZF_READ(ifileh,off,rectst,NB_IN,4*mpi_npet)
           ALLOCATE(ieposgt(2,nwpcst(0)))
           off=ZF_READ(ifileh,off,ieposgt,2,2*nwpcst(0))
           call ZF_close(ifileh)
       end if
#ifndef NO_MPI
       call wav_mpi_bcast(irtagt);
       call wav_mpi_bcast(mpi_npet);
       call wav_mpi_bcast(nwpcst(0));
#endif
       if(irtagt/=irtag.or.mpi_npet/=intarr(3))then
           DBGO(0,*) "Restart File  ",RestFile,"& Restart HeadFile ",RestHead," Not Match"
           call wav_abort("Restart File  &Restart HeadFile Not Match")
       end if
       if(mpi_id/=0)then
           ALLOCATE(nwpcst(mpi_npet),rectst(4*mpi_npet))
           ALLOCATE(ieposgt(2,nwpcst(0)))
       endif
       allocate(ioffs (mpi_npet))
#ifndef NO_MPI
       call wav_mpi_bcast(ndept);
       call wav_mpi_bcast(nwpcst);
       call wav_mpi_bcast(rectst);
       call wav_mpi_bcast(ieposgt);
       call wav_mpi_bcast(offt );
#endif
       call CalRestOff(offt,nwpcst,mpi_npet,ioffst)

       maxnwpc=maxval (nwpcst(1:mpi_npet))
       allocate(v2d(maxnwpc),vee(kl,jnthet,0:999));
       ipoff=0
       do i=1,mpi_npet
         if(.not.intersect(rectst((i-1)*4+1)) )cycle
         off=ioffs(i-1);
         ncwp=nwpcst(i);icwp=0;nicwp=0
         do while(ncwp>0)
           icwp=icwp+nicwp;nicwp=min(ncwp,1000)
           off=ZF_READ(ifile,off,vee,NB_RN,nicwp*kl*jnthet);
           call cp4dvar(vee,nicwp,ipoff+icwp,ieposgt,eec)
           ncwp=ncwp-nicwp
         enddo
         ncwp=nwpcst(i);
         off=ZF_READ(ifile,off,v2d,4,   ncwp)
         call cp2dvar(v2d,ncwp,ipoff,ieposgt,wind0x)
         off=ZF_READ(ifile,off,v2d,4,   ncwp)
         call cp2dvar(v2d,ncwp,ipoff,ieposgt,wind0y)
         off=ZF_READ(ifile,off,v2d,4,   ncwp)
         call cp2dvar(v2d,ncwp,ipoff,ieposgt,windx)
         off=ZF_READ(ifile,off,v2d,4,   ncwp)
         call cp2dvar(v2d,ncwp,ipoff,ieposgt,windy)
         !off=ZF_READ(ifile,off,v2d,4,   ncwp)
         !call cp2dvar(v2d,ncwp,ipoff,ieposgt,wiv)
         !off=ZF_READ(ifile,off,v2d,4,   ncwp)
         !call cp2dvar(v2d,ncwp,ipoff,ieposgt,winc)
         do ko=1,Navhists
           if(aoneas(ko)/=0)then
             ncwp=nwpcst(i);icwp=0;nicwp=0
             do while(ncwp>0)
               icwp=icwp+nicwp;nicwp=min(ncwp,1000)
               off=ZF_READ(ifile,off,vee,NB_RN,nicwp*kl*jnthet);
               call cp4dvar(vee,nicwp,ipoff+icwp,ieposgt,avhists(ko)%ea)
               ncwp=ncwp-nicwp
             enddo
           end if
         enddo
         ipoff=ipoff+nwpcst(i)
       end do
       call ZF_close(ifile);
       deallocate(ioffs,nwpcst,v2d,vee)
     end if
 end SUBROUTINE read_restart
   SUBROUTINE cp2dvar(vs,nns,ioff,iepost,vd)
     real(4):: vs(:),vd(:)
     integer*2 iepost(2,*)
     integer nns,ioff
     integer ix,iy,i,n
     do i=1,nns
       ix=iepost(1,i+ioff);ix=iepost(2,i+ioff);
       n=pos2id(ix,iy)
       if(n/=0)vd(n)=vs(i)
     enddo
   end SUBROUTINE cp2dvar
   SUBROUTINE cp4dvar(ves,nns,ioff,iepost,ved)
     REALD ves(:,:,:),ved(:,:,:)
     integer*2 iepost(2,*)
     integer nns,ioff
     integer ix,iy,i,n
     do i=1,nns
       ix=iepost(1,i+ioff);ix=iepost(2,i+ioff);
       n=pos2id(ix,iy)
       if(n/=0)ved(:,:,n)=ves(:,:,i)
     enddo
   end SUBROUTINE cp4dvar

end Module restart_mod

