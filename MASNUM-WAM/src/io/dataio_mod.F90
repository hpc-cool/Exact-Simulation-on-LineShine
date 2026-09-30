#include "wavedeff.h"
Module DataIO_mod
  use varcommon_mod
  use netcdf_mod
  use partition_mod
  use wav_cal_mod
#ifdef NO_MPI
#define GINTERW
#else
#define GINTERW
#endif
IMPLICIT NONE
#ifdef GINTERW
    real(4) ,allocatable::gwind(:,:)
    real(4) ,allocatable::Interv(:,:)
#endif
    real(4) ,allocatable::windi(:,:)
    integer(2) ,allocatable::swindi(:,:)
    real(4) ,allocatable::lat(:),lon(:)
    real(8) ,allocatable::times(:)
    real(4) ,allocatable::ftimes(:)
    real(8) ::uoff,uscale,voff,vscale,timeoff,timescale
    real(8),allocatable::dlat(:),dlon(:)    
    integer ::nlat=0,nlon=0,ntime=0,itimec=1
    integer ::ncidu=-1,ncidv=-1,wuid,wvid,uwtype,vwtype,uos,vos
    INTEGER ::newwindforce=1
    character(8)::date_labels(6) = ['years   ','months  ','days    ','hours   ','minutes ','seconds ' ]
    real(8) ::datescales(6) = [365.,30.4,1.,1/24.,1/1440.,1/86400.]
    character*(1024) ::ncfnu,ncfnv,ncfnuo="NONE",ncfnvo
contains
subroutine InitWindIn
#ifdef GINTERW
   allocate(gwind(gixl,giyl))       
#endif   
   allocate(dlat(1),dlon(1))   
end subroutine InitWindIn
subroutine formfn(fn,form,y,m,d)
  INTEGER y,m,d
  character*(*) fn,form
  character*(16) tt,ft
  character c
  integer i,l,j,n,k,fl
  fn='';
  fl=len(fn)
  l=LEN_TRIM(form)
  i=1;j=1
  fn=''
  do while(i<=l) 
    c=form(i:i);
    if(c/='%')then
      fn(j:j)=c;j=j+1;
    else
      i=i+1;c=form(i:i)
      if(c=='%')then
        fn(j:j)=c;j=j+1
      else 
        ft='(i';k=3
        do while(i<=l) 
          c=form(i:i)         
          if((c>='0'.and.c<='9').or.c=='.')then
            ft(k:k)=c;k=k+1
          else
            exit
          endif
          i=i+1
        enddo
        ft(k:k)=')'
        if(c=='P')then
          n=LEN_TRIM(pathwind)
          ft(2:2)='a';fn(j:j+n-1)=trim(pathwind)
        else if(c=='Y')then
          ft(2:2)='i';write(fn(j:fl),ft)y
        else if(c=='M')then       
          ft(2:2)='i';write(fn(j:fl),ft)m
        else if(c=='D')then       
          ft(2:2)='i';write(fn(j:fl),ft)d
        endif
        j=LEN_TRIM(fn)+1
      endif
    endif
    i=i+1
  enddo
end subroutine formfn
subroutine GetDataFn(form,eday)
  character*(*) form
  integer eday
  integer y,m,d
  call wav_cal_eday2ymd(eday,y,m,d)
  !DBGINF,eday,y,m,d
  call formfn(ncfnu ,ufnform,y,m,d)
  if(vfnform=='')then
    ncfnv=ncfnu
  else
    call formfn(ncfnv ,vfnform,y,m,d)
  endif  
end subroutine GetDataFn
subroutine ReadTime
  integer ntimet,vid,ttype
  character(256) timeunits,stime
  integer status,n,y,m,d,h,mn,s,i,days
  real*8 dsec,vt,vts
  ntimet=get_dimension_len(ncidu,'time')
  if(ntimet/=ntime)then
    if(allocated(times))deallocate(times)
    if(allocated(ftimes))deallocate(ftimes)
    ntime=ntimet;allocate(times(ntime),ftimes(ntime))
  endif
  vid=GetVarId(ncidu,'time')
  ttype=GetVarType(ncidu,vid)
  
  if(ttype==nf_real)then    
    call readnc(ncidu,vid,ftimes)
    times=ftimes
  else if(ttype==nf_double)then
    call readnc(ncidu,vid,times)
  else
    
  endif
  itimec=1;timescale=1;timeoff=0
  timeunits='';
  status=get_attribute(ncidu,vid,"units"  ,timeunits)
  y=1950;m=1;d=1;h=0;m=0;s=0;
  !DBGO0(0,*)trim(timeunits),timescale
  if(status==NF_NOERR)then
    do i=1,6
      n=INDEX(timeunits,trim(date_labels(i)))
      if(n>0)then
        timescale=datescales(i)       
      endif
    enddo
    stime=timeunits(INDEX(timeunits,'since')+5:);
    n=LEN_TRIM(stime);
    if(stime(n:n)=='.')stime(n:n)=' '
    n=-1
    
    
    stime=stime(n+2:);n=INDEX(stime,'-'    )-1;call flush(6);read(stime(1:n),*)y
    stime=stime(n+2:);n=INDEX(stime,'-'    )-1;call flush(6);read(stime(1:n),*)m
    stime=stime(n+2:);n=INDEX(stime,' '    )-1;call flush(6);read(stime(1:n),*)d
    stime=stime(n+2:);n=INDEX(stime,':'    )-1;call flush(6);read(stime(1:n),*)h
    stime=stime(n+2:);n=INDEX(stime,':'    )-1;call flush(6);read(stime(1:n),*)mn
    stime=stime(n+2:);n=INDEX(stime,' '    )-1;call flush(6);read(stime(1:n),*)dsec
    s=dsec
  endif 
  call wav_cal_ymd2eday(y,m,d,days);
  days=days+WTIME_PATCH;
  timeoff=days+hms2day(h,mn,s)
  timescale=timescale
  vts=1./86400.
  do n=1,ntime
    vt=times(n)*timescale+timeoff+0.01*vts
    times(n)=vt-mod(vt,vts)
  enddo
  !times=times*timescale+timeoff
  !DBGO0(0,*)'read WIND time 1: ',trim(ncfnu),WTIME_PATCH,runstate%timecur,times(1),times(ntime);
  !stop
end subroutine ReadTime
subroutine OpenWind()
  integer nlatt,nlont,ttype,vid
  integer status
  if(ncfnu==ncfnuo)return
  if(ncidu>0)then
    if(ncidv/=ncidu)call close_nc(ncidv)
    call close_nc(ncidu) 
  endif
  !DBGO(1,*)' Open Nc for u or uv '//trim(ncfnu)
  call open_nc(ncidu,trim(ncfnu),'read')  
  ncfnuo=ncfnu
  if(ncfnu/=ncfnv)then
    !DBGO(1,*)' Open Nc for v       '//trim(ncfnv)
    call open_nc(ncidv,trim(ncfnv),'read')
  else
    ncidv=ncidu
  endif
  nlont=get_dimension_len(ncidu,'lon')
  nlatt=get_dimension_len(ncidu,'lat')
  wuid=GetVarID(ncidu,uwind_name)
  wvid=GetVarID(ncidv,vwind_name)
  uoff=0;uscale=1;voff=0;vscale=1;uos=0;vos=0
  uwtype=GetVarType(ncidu,wuid)
  vwtype=GetVarType(ncidv,wvid)
  status=get_attribute(ncidu,wuid,"add_offset"  ,uoff  )
  status=get_attribute(ncidu,wuid,"scale_factor",uscale)
  status=get_attribute(ncidv,wvid,"add_offset"  ,voff  )
  status=get_attribute(ncidv,wvid,"scale_factor",vscale)  
  if(abs(uoff)+abs(uscale-1)>1e-10)uos=1
  if(abs(voff)+abs(vscale-1)>1e-10)vos=1
  if(nlatt/=nlat.or.nlont/=nlon)then
    if(allocated(windi))deallocate(windi)
    if(allocated(swindi))deallocate(swindi)
    allocate(windi(nlont,nlatt))
    if(uwtype==nf_int2)then
      allocate(swindi(nlont,nlatt))
    endif
    newwindforce=1
  endif  
  if(nlatt/=nlat)then
    nlat=nlatt;    
    if(allocated(lat))deallocate(lat)
    if(allocated(dlat))deallocate(dlat)
    allocate(lat(nlat),dlat(nlat))
#ifdef GINTERW
    if(allocated(interv))deallocate(interv)
    allocate(interv(gixl,nlat))
#endif    
  endif
  
  if(nlont/=nlon)then
    if(allocated(lon))deallocate(lon)
    if(allocated(dlon))deallocate(dlon)
    nlon=nlont;allocate(lon(nlon),dlon(nlon))
  endif
    
  ttype=GetVarType(ncidu,'lon')
  if(ttype==nf_real)then
    call readnc(ncidu,'lon',lon)
    call readnc(ncidu,'lat',lat)
    dlon=lon;dlat=lat
    newwindforce=1
  else if(ttype==nf_double)then
    call readnc(ncidu,'lon',dlon)
    call readnc(ncidu,'lat',dlat)
    newwindforce=1
  endif
  call ReadTime
end subroutine OpenWind
integer function IFindtime(time)
    real(8),intent(in):: time
    real(8) timet
    integer i
    timet=time-1e-4
    !print*,'IFindtime=',time,times
    do i=1,ntime
      if(timet<=times(i))then
        exit
      endif
    enddo
    IFindtime=0;itimec=i
    if(i>ntime)return
    IFindtime=i;
end function IFindtime
integer function IvFindtime(time)
  real(8),intent(in):: time
  ncfnu=' '
  IvFindtime=IFindtime(time)
  if(IvFindtime>0)return
  call GetDataFn("",int(time))
  !DBGO(0,*)time,int(time),trim(ncfnu)
  call OpenWind
  IvFindtime=IFindtime(time)  
  !DBGO(0,*)'IvFindtime',IvFindtime
end function IvFindtime
#ifdef GINTERW
subroutine InterGrdreg(sv,dv)
  real(4)::sv(:,:),dv(:,:)
  integer i,j,ij,it,it1,jt,jt1,itc
  logical iflon
  real x,y,p
  itc=1;
  iflon=.false.; if(lon(1)>=0)iflon=.true.
  do i=1,gixl
    x=xcord(i,1)
    if(iflon)then
      if(x<0)x=x+360.
      if(itc>1)then
        if(x<lon(itc-1))itc=1
      endif
    endif
    do ij=itc,nlon
      if(x<lon(ij))exit
    enddo
    itc=ij;it1=ij;it=it1-1
    if(it<=0)then
      it=nlon; p=(x+360-lon(it))/(lon(it1)+360-lon(it))
    else if(it1>=nlon)then
      it1=1; p=(x-lon(it))/(lon(it1)+360-lon(it))
    else
      p=(x-lon(it))/(lon(it1)-lon(it))
    endif
    if(p<1e-8)then
      interv(i,:)=sv(it,:)
    else if(1-p<1e-8)then
      interv(i,:)=sv(it1,:)
    else
      interv(i,:)=sv(it,:)*(1-p)+sv(it1,:)*p    
      !write(6,'(a,4i4,4f8.2)')'A',i,j,it,it1,sv(it,j),sv(it1,j),p,Interv(i,j)
    endif
  enddo
  !print*,interv
  if(lat(1)<lat(nlat))then
    jt1=1
    do j=1,giyl
      y=ycord(1,j)
      do ij=jt1,nlat
        if(y<lat(ij))exit
      enddo
      jt1=ij;jt=jt1-1
      if(jt<=0)then
        jt=jt1;p=0
      else if(jt1>nlat)then
        jt=nlat;jt1=nlat;p=0
      else
        p=(y-lat(jt))/(lat(jt1)-lat(jt))
      endif
      do i=1,gixl
        dv(i,j)=interv(i,jt)*(1-p)+interv(i,jt1)*p      
      enddo
    enddo
  else
    jt1=nlat
    do j=1,giyl
      y=ycord(1,j)
      do ij=jt1,1,-1
        if(y<lat(ij))exit
      enddo
      jt1=ij;jt=jt1+1
      if(jt>nlat)then
        jt=jt1;p=0
      else if(jt1<1)then
        jt=1;jt1=1;p=0
      else
        p=(y-lat(jt))/(lat(jt1)-lat(jt))
      endif
      if(jt==jt1)then
        dv(:,j)=interv(:,jt)
      else
        do i=1,gixl
          dv(i,j)=interv(i,jt)*(1-p)+interv(i,jt1)*p      
        enddo
      endif
    enddo
  endif
end subroutine InterGrdreg
#endif
!======================================
! mode:
!   0: Init MPI msg etc
!   1: Init before startup
!   2: Get First Data for cold start
!   3: Get First Data for warm start
!   4: Get Data For for timeStep
!======================================
integer function IRDataIO (mode,time,windx_,windy_,timec,timep)
    integer,intent(in):: mode
    real(8),intent(in)::time
    real(8),intent(out)::timec,timep
    real(4),intent(out)::windx_(:),windy_(:)
    real(8) timet,dt
    integer it
    IRDataIO=0      
    !DBGO0(0,*)'RW',mode,time,timec
    if(mode==0)then
      IRDataIO=1
      return 
    else if(mode==1)then
      if(mpi_id==0)then
        call InitWindIn
      endif      
      IRDataIO=1
      return 
    else if(mode==-1)then
      IRDataIO=1
      return 
    endif 
    
  if(mpi_id==0)then     
    it=IvFindtime(time)
    timec=-1
    if(it>0)then
      timec=times(it)
      if(it>1)then
        timep=times(it-1)
      else if(it<ntime)then
        timep=timec*2-times(it+1)
      else 
        timep=timec-1/24.
      endif
    endif
    !DBGO0(0,*)'Got WindTime',timep,timec,time
  endif
#ifndef NO_MPI
  call wav_mpi_bcast(timec);call wav_mpi_bcast(timep)
#endif
  if(timec<0)   return
#ifndef NO_MPI
#ifndef GINTERW
  if(newwindforce/=0)then
    call scatter_wind_init( dlon, dlat,  nlon, nlat,nlon,nlat, 360.d0)
  endif
#endif  
#endif
  if(mpi_id==0)then    
    !read uwind
    if(uwtype==nf_int2)then
      call readnc(ncidu,uwind_name,swindi,it)
      windi=(swindi)*uscale+uoff
    else if(uwtype==nf_real)then
      call readnc(ncidu,uwind_name,windi,it)
      if(uos/=0)windi=(windi)*uscale+uoff
    endif
  endif
  !print*,'windu=',windi(nlon/3,nlat/3)
#ifdef GINTERW
  if(mpi_id==0)call InterGrdreg(windi,gwind)
  call ScatterVar2d(gwind,windx_)
#else  
  call scatter_wind(windi,windx_)
#endif   
  if(mpi_id==0)then   
    !read vwind
    if(uwtype==nf_int2)then
      call readnc(ncidv,vwind_name,swindi,it)
      windi=(swindi)*vscale+voff
    else if(uwtype==nf_real)then
      call readnc(ncidv,vwind_name,windi,it)
      if(vos/=0)windi=(windi)*vscale+voff
    endif
  endif
#ifdef GINTERW
  if(mpi_id==0)call InterGrdreg(windi,gwind)
  call ScatterVar2d(gwind,windy_)
#else  
  call scatter_wind(windi,windy_)
#endif    
  !DBGO0(0,*)'R',mpi_id,maxval(windx_),maxval(windy_)
  IRDataIO=1
end function IRDataIO
end Module DataIO_mod
