!BOP ===========================================================================
! !MODULE: wav_cal_mod -- calendar module, relates elapsed days to calendar date.
! !REMARKS:
!   Following are some internal assumptions.  These assumptions are somewhat
!   arbitrary -- they were chosen because they result in the simplest code given
!   the requirements of this module.  These assumptions can be relaxed as
!   necessary:
!   o the valid range of years is [0,9999]
!   o elapsed days = 0 <=> January 1st, year 0000
!   o all years have 365 days (no leap years)
!     This module is hard-coded to implement a 365-day calendar, ie. there
!     are no leap years.  This module can be modified to implement a calendar
!     with leap years if this becomes desireable.  This would make the internal
!     logic of this module more complex, but would not require any change to the
!     module API or the calling code because the module API hides these details
!     from all external routines.
! !INTERFACE: ------------------------------------------------------------------
module wav_cal_mod
! !USES:
   !use debughlp_mod
   implicit none
   private ! except
! !PUBLIC MEMBER FUNCTIONS:
   public :: SetBaseYear,ctime2day,day2ctime,day2cdatetime,hms2day
   public :: yearsday,IsLeapYear
   public :: dsm,leapdsm,baseDays
   public :: wav_cal_eday2date  ! converts elapsed days to coded-date
   public :: wav_cal_eday2ymd   ! converts elapsed days to yr,month,day
   public :: wav_cal_date2ymd   ! converts coded-date   to yr,month,day
   public :: wav_cal_date2eday  ! converts coded-date   to elapsed days
   public :: wav_cal_ymd2date   ! converts yr,month,day to coded-date
   public :: wav_cal_ymd2eday   ! converts yr,month,day to elapsed days
   public :: wav_cal_validDate  ! logical function: is coded-date valid?
   public :: wav_cal_validYMD   ! logical function: are yr,month,day valid?
! !PUBLIC DATA MEMBERS:
   integer,parameter :: dsm(13) = &     ! elapsed Days on Start of Month
   &                    (/ 0,31,59, 90,120,151, 181,212,243, 273,304,334,365/)
   integer,parameter :: dpm(12) = &     ! Days Per Month
   &                    (/31,28,31, 30, 31, 30,  31, 31, 30,  31, 30, 31/)
   integer ::leapdsm(13),leapdpm(12)
   integer ::Baseyear=1950,baseDays=711857
  real*8,parameter::DaysPerRyear=365.24219879D0
  real*8   ::DaysPeryear=365
  integer ::LEAPYEAR=0,maxyears=1000  
contains
subroutine SetBaseYear(UseLeap,year0,maxyear) 
  integer,optional::year0,UseLeap,maxyear
  integer i,yt,dt,d0
  if(present(year0))Baseyear=year0
  if(present(UseLeap))LEAPYEAR=UseLeap
  if(present(maxyear))maxyears=maxyear-year0  
  baseDays=yearsday(Baseyear);
  leapdsm=0;leapdpm=0
  if(LEAPYEAR/=0)then   
    DaysPeryear=DaysPerRyear
    leapdsm(3:13)=1
    leapdpm(2)=1
  else
    DaysPeryear=365
  endif
end subroutine SetBaseYear

  integer function IsLeapYear(year)
    integer year,yt
    !if((mod(year, 4)==0 .and. mod(year, 100) /= 0) .or. (mod(year, 400)==0))then
    ! IsLeapYear=1
    !endif
    yt=year
    IsLeapYear=yt/4-yt/100+yt/400
    yt=yt-1
    IsLeapYear=IsLeapYear-(yt/4-yt/100+yt/400)
 end function IsLeapYear

!BOP ===========================================================================
logical function wav_cal_validYMD(year,month,day)
   integer,intent(in ) :: year,month,day  ! calendar year,month,day
   wav_cal_validYMD = .true.
   if (year  <    0) wav_cal_validYMD = .false.
   if (year  >99999) wav_cal_validYMD = .false.
   if (month <    1) wav_cal_validYMD = .false.
   if (month >   12) wav_cal_validYMD = .false.
   if (day   <    1) wav_cal_validYMD = .false.
   if (wav_cal_validYMD) then
      if (day > dpm(month)+leapdpm(month)*IsLeapYear(year)) wav_cal_validYMD = .false.
   endif
end function wav_cal_validYMD
subroutine  wav_cal_ymd2eday(year,month,day,eday)
   integer,intent(in ) :: year,month,day  ! calendar year,month,day
   integer,intent(out) :: eday            ! number of elapsed days
   if (.not. wav_cal_validYMD(year,month,day)) stop 'wav_cal_validYMD'
   eday = yearsday(year)-baseDays + dsm(month)+leapdsm(month)*IsLeapYear(year) + day-1
end subroutine  wav_cal_ymd2eday

subroutine wav_cal_eday2ymd (eday,year,month,day)
  integer,intent(in)  :: eday             ! elapsed days
  integer,intent(out) :: year,month,day   ! calendar year,month,day
  !--- local ---
  integer :: k,md
  real dt,dt1
  dt1=eday+baseDays
  year=(dt1-10)/DaysPeryear+1
  ! if(eday>yearsday(year+1))year=year+1
  dt=dt1+1-yearsday(year+1)
  !DBGO0(0,*)'AA1',dt,year,dt1,dt1/DaysPeryear,(dt1-1)/DaysPeryear
  year=year+ifix((abs(dt)+dt)/(abs(dt)*2-1e-30))
  day  = dt1-yearsday(year)  ! elapsed days within current year
  month=day/31.001+1+1  ! make sure month is bigger
  md=dsm(month)+leapdsm(month)*IsLeapYear(year)
  dt=md-day
  month=month-ifix((abs(dt)+dt)/(abs(dt)*2-1e-30))
  md=dsm(month)+leapdsm(month)*IsLeapYear(year)
  day = day-md+1          ! calendar day
end subroutine wav_cal_eday2ymd

!============= Fore Date =============================
subroutine wav_cal_date2ymd (date,year,month,day)
   integer,intent(in)  :: date             ! coded-date (yyyymmdd)
   integer,intent(out) :: year,month,day   ! calendar year,month,day
   year =     date       /10000
   month= mod(date,10000)/  100
   day  = mod(date,  100)
   if (.not. wav_cal_validYMD(year,month,day)) then
      write(6,*) "(cal_date2ymd) ERROR: invalid date = ",date
   endif
end subroutine wav_cal_date2ymd
subroutine wav_cal_date2eday(date,eday)
   integer,intent(in ) :: date            ! coded (yyyymmdd) calendar date
   integer,intent(out) :: eday            ! number of elapsed days
   !--- local ---
   integer :: year,month,day
   call wav_cal_date2ymd(date,year,month,day)
   call wav_cal_ymd2eday(year,month,day,eday)
end subroutine wav_cal_date2eday
subroutine wav_cal_ymd2date(year,month,day,date)
   integer,intent(in ) :: year,month,day  ! calendar year,month,day
   integer,intent(out) :: date            ! coded (yyyymmdd) calendar date
! NOTE:
!   this calendar has a year zero (but no day or month zero)
   if (.not. wav_cal_validYMD(year,month,day)) stop 'wav_cal_validYMD'
   date = year*10000 + month*100 + day  ! coded calendar date
end subroutine wav_cal_ymd2date
subroutine wav_cal_eday2date(eday,date)
   integer,intent(in)  :: eday  ! number of elapsed days
   integer,intent(out) :: date  ! coded (yyyymmdd) calendar date
   !--- local ---
   integer :: k,year,month,day
! ASSUMPTIONS:
!   this calendar has a year zero (but no day or month zero)
  call wav_cal_eday2ymd (eday,year,month,day)
  date = year*10000 + month*100 + day  ! coded calendar date
end subroutine wav_cal_eday2date

logical function wav_cal_validDate(date)
   integer,intent(in ) :: date            ! coded (yyyymmdd) calendar date
   !--- local ---
   integer :: year,month,day
   call wav_cal_date2ymd(date,year,month,day)
   wav_cal_validDate = wav_cal_validYMD(year,month,day)
end function wav_cal_validDate
real*8 function ctime2day(ctime)
  integer ctime
  integer h,m,s,c
  c=ctime
  h=c/10000;c=c-h*10000
  m=c/100  ;s=c-m*100
  ctime2day=((h*60+m)*60+s)/86400.
end function ctime2day
real*8 function hms2day(h,m,s)
  integer ctime
  integer h,m,s,c
  hms2day=((h*60+m)*60+s)/86400.
end function hms2day

integer function day2ctime(daytime)
  real(8) ::daytime,dt
  integer ::h,m,s,c
  dt=(daytime-int(daytime))*24
  h=dt;dt=(dt-h)*60;
  m=dt;s=(dt-m)*60;
  day2ctime=h*10000+m*100+s
end function day2ctime

subroutine day2cdatetime(daytime,cdate,ctime)
  real(8) daytime,dt
  integer cdate,ctime
  dt=daytime+0.5/86400.
  call wav_cal_eday2date(int(dt),cdate)
  ctime=day2ctime(dt)  
end subroutine day2cdatetime
  integer function yearsday(year)
    integer year,yt
    yt=year-1
    yearsday=yt*365+(yt/4)-(yt/100)+(yt/400)
 end function yearsday
end module wav_cal_mod
