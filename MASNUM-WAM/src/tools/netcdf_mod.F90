!-------------------------------------------------------------------------------
!* netcdf_mod 2.6
  module netcdf_mod
!-------------------------------------------------------------------------------
! ******************************************************************************
!-------------------------------------------------------------------------------
!                                    Copyright (C) 20014 Zhaowei
!                                    MODULE NAME : netcdf_mod
!                                    History : netcdf_mod 2.6
!                                           change from xunqiang yin
!                                    History : netcdf_mod 2.1
!                                    History : netcdf_mod 2.0
!                                    History : NcFileData 2.0
!                                    Current VERSION : 2008/03/31
!
! --- USAGE : To output/input nc files in convenience.
! --- DEPEND: The package of netcdf for FORTRAN.
!
! --- NOTE for describing of subroutine / function :
!  A. The parameters bracketed with [], means optional parameter.
!  B. The describe for the parameters of subroutine / function, started with:
!   * It means input prameter;
!   # It means output prameter;
!   @ It means input and output prameter(it will be changed inside).
!
!-------------------------------------------------------------------------------
! ******************************************************************************
! ***                       INTERFACE DESCRIBE                               ***
! ******************************************************************************
!-------------------------------------------------------------------------------
!
!  1. subroutine open_nc(NcID, FileName, action)
!     # integer :: NcID = unit number of opened netcdf file.
!     * character(len=*) :: FileName = the name of netcdf file.
!     * character :: action = 'create', 'define', 'write' or 'read'.
!
!-------------------------------------------------------------------------------
!
!  2. subroutine close_nc(NcID)
!     * integer :: NcID = unit number of opened netcdf file.
!
!-------------------------------------------------------------------------------
!
!  3. subroutine end_define(NcID)
!     * integer :: NcID = unit number of opened netcdf file.
!
!-------------------------------------------------------------------------------
!
!  4. subroutine dimension_define 
!     - dimension_define(NcID, DimName, DimLen[,DimID,[ DimVarName, DimVarType [, &
!                              DimVarID]]])
!     * integer :: NcID = unit number of opened netcdf file.
!     * character :: Dimname = the name of this dimension.
!     * integer :: DimLen = the size of this dimension.
!     * character :: DimVarName = the name of dimension variable name.
!     * integer :: DimVarType = the type of dimension variable (use integer).
!          ( nf_int1 = 1, nf_char = 2, nf_int2 = 3, nf_int = 4,
!            nf_real = 5, nf_double = 6 )
!     # integer :: DimID = the record number of this dimension.
!     # integer :: DimVarID = the record number of dimension variable.
!     ( Note: for the unlimited dimension, 'DimLen = 0 '.)
!
!-------------------------------------------------------------------------------
!
!  5. function get_dimension_len(NcID, DimName)
!     * integer :: NcID = unit number of opened netcdf file.
!     * character(len=*) :: DimName = the name of dimension need to get length.
!     # integer :: get_dimension_len = the needed dimension length.
!
!-------------------------------------------------------------------------------
!
!  6. subroutine variable_define
!
!     --- 2 interfaces:
!
!    6.1 variable_define(NcID, VarName, VarType, VarDimsName [, VarID])
!     * integer :: NcID = unit number of opened netcdf file.
!     * character :: VarName = the name of this variable.
!     * integer :: VarType = the type of variable (use integer).
!          ( nf_int1 = 1, nf_char = 2, nf_int2 = 3, nf_int = 4,
!            nf_real = 5, nf_double = 6 )
!     * character(len=*) :: VarDimsName = the name list to which this variable
!                           related.
!     # integer :: VarID = the record number of this variable.
!
!    6.2 variable_define(NcID, VarName, VarType, VarDims [, VarID])
!     * integer :: NcID = unit number of opened netcdf file.
!     * character :: VarName = the name of this variable.
!     * integer :: VarType = the type of variable (use integer).
!          ( nf_int1 = 1, nf_char = 2, nf_int2 = 3, nf_int = 4,
!            nf_real = 5, nf_double = 6 )
!     * integer :: VarDims = the related DimID of this variable.
!     # integer :: VarID = the record number of this variable.
!
!-------------------------------------------------------------------------------
!
!  7. integer function set_attribute(NcID, AttName, Att [, VarName])
!     * integer :: NcID = unit number of opened netcdf file.
!     * character(len=*) :: AttName = the attribute name of netcdf file.
!     * character(len=*)/integer/integer*1/integer*2/real/double presion ::
!                           Att = attribute of opened netcdf file.
!     * character(len=*) :: VarName = the owner(variable) of attribute.
!     ( Note: for the globle attribute, VarName shouldn't be given.)
!
!-------------------------------------------------------------------------------
!
!  8. integer function get_attribute(NcID, AttName, Att [, VarName])
!     * integer :: NcID = unit number of opened netcdf file.
!     * character(len=*) :: AttName = the attribute name of netcdf file.
!     # character(len=*)/integer/integer*1/integer*2/real/double presion ::
!                           Att = attribute of opened netcdf file.
!     * character(len=*) :: VarName = the owner(variable) of attribute.
!     ( Note: for the globle attribute, VarName shouldn't be given.)
!
!-------------------------------------------------------------------------------
!
!  9. subroutine writenc(NcID, VarName, Var [, RecNum, locs])
!     * integer :: NcID = unit number of opened netcdf file.
!     * character(len=*) :: VarName = the name of variable in netcdf file.
!     * character(len=*)/char/integer/integer*1/integer*2/real/double presion ::
!                          Var = the variable needed to be output.
!     * integer :: RecNum = the record number at unlimited dimension.
!     * integer :: locs = the start locations for output.
!     ( Note: for the variable without unlimited dimension,
!             RecNum shouldn't be given.)
!
!-------------------------------------------------------------------------------
!
!  10. subroutine readnc(NcID, VarName, Var [, RecNum, locs])
!     * integer :: NcID = unit number of opened netcdf file.
!     * character(len=*) :: VarName = the name of variable in netcdf file.
!     # character(len=*)/char/integer/integer*1/integer*2/real/double presion ::
!                         Var = the variable needed to be input.
!     * integer :: RecNum = the record number at unlimited dimension.
!     * integer :: locs = the start locations for output.
!     ( Note: for the variable without unlimited dimension,
!             RecNum shouldn't be given.)
!
!-------------------------------------------------------------------------------
!
!                                                 E-Mail: zhaow@fio.org.cn
!
!-------------------------------------------------------------------------------
! ******************************************************************************
!-------------------------------------------------------------------------------
!  use debughlp_mod
  implicit none
!-------------------------------------------------------------------------------
  public ::open_nc, close_nc, end_define
  public ::readnc, writenc
  public ::GetVarid,GetVarType
  public ::variable_define, set_attribute, get_attribute
  public ::dimension_define, get_dimension_len, ncdump
  include 'netcdf.inc'
#ifndef GF13
  integer(8),external::  ZF_size   ! (ifile)
  integer(4),external::  ZF_open   ! (ifile,fn)
  integer(8),external::  ZF_skip   ! (ifile,off,NByte,nItem)
  integer(8),external::  ZF_READ   ! (ifile,off,var,NByte,nItem)
  integer(8),external::  ZF_WRITE  ! (ifile,off,var,NByte,nItem)
#else
#if 1
  integer(8),external::  ZF_READ   ! (ifile,off,var,NByte,nItem)
  integer(8),external::  ZF_WRITE  ! (ifile,off,var,NByte,nItem)
#else
  interface
    integer(8) function ZF_READ(ifile,off,var,NByte,nItem)
      integer, intent(in)    :: ifile,NByte,nItem
      integer(8),intent(in)    :: off
      type(*), dimension(*) ::var 
    end function
    integer(8) function ZF_WRITE(ifile,off,var,NByte,nItem)
      integer, intent(in)    :: ifile,NByte,nItem
      integer(8),intent(in)    :: off
      type(*), dimension(*) ::var 
    end function
  end interface
#endif
  interface
    SUBROUTINE zf_mkdir(fn)
      character*(*), intent(in) ::fn 
    end SUBROUTINE 
    SUBROUTINE zf_setpid (mpiid)
      integer, intent(in)    ::mpiid 
    end SUBROUTINE 
    integer(8) function ZF_size(ifile)
      integer, intent(in)    :: ifile
    end function 
    integer(4) function ZF_open(ifile,fn)
      integer, intent(in)    :: ifile
      character*(*), intent(in) ::fn 
    end function ZF_open
    integer(8) function ZF_skip(ifile,off,NByte,nItem)
      integer, intent(in)    :: ifile,NByte,nItem
      integer(8),intent(in)    :: off
    end function 
    SUBROUTINE zf_close(ifile)
      integer, intent(in)    :: ifile
    end SUBROUTINE
  end interface
#endif
!
!-------------------------------------------------------------------------------
  interface variable_define
    module procedure variable_define1, variable_define2
  end interface variable_define
!-------------------------------------------------------------------------------
  interface SetVARSC
      module procedure SetVARSC_TEXT   , &
                       SetVARSC_int1   , &
                       SetVARSC_int2   , &
                       SetVARSC_int    , &
                       SetVARSC_REAL   , &
                       SetVARSC_DOUBLE
  end interface SetVARSC
  interface set_attribute
    module procedure set_attribute_i_TEXT  ,set_attribute_n_TEXT  , &
                     set_attribute_i_int1  ,set_attribute_n_int1  , &
                     set_attribute_i_int2  ,set_attribute_n_int2  , &
                     set_attribute_i_int   ,set_attribute_n_int   , &
                     set_attribute_i_REAL  ,set_attribute_n_REAL  , &
                     set_attribute_i_DOUBLE,set_attribute_n_DOUBLE
  end interface set_attribute
!-------------------------------------------------------------------------------
  interface get_attribute
    module procedure get_attribute_i_TEXT  ,get_attribute_n_TEXT  , &
                     get_attribute_i_int1  ,get_attribute_n_int1  , &
                     get_attribute_i_int2  ,get_attribute_n_int2  , &
                     get_attribute_i_int   ,get_attribute_n_int   , &
                     get_attribute_i_REAL  ,get_attribute_n_REAL  , &
                     get_attribute_i_DOUBLE,get_attribute_n_DOUBLE
  end interface get_attribute
! interface AttSize
!   module procedure AttSize_TEXT, &
!                    AttSize_int1, &
!                    AttSize_int2, &
!                    AttSize_int4, &
!                    AttSize_real, &
!                    AttSize_double
! end interface AttSize
!-------------------------------------------------------------------------------
  interface GetVarType
    module procedure GetVarType_n,GetVarType_i
  end interface GetVarType
!-------------------------------------------------------------------------------
  interface writenc
    module procedure &
    writenc_i0d_int1, writenc_i0d_int2,   writenc_i0d_int,      &
    writenc_i0d_real, writenc_i0d_double, writenc_i0d_text,     &
    writenc_i1d_int1, writenc_i1d_int2,   writenc_i1d_int,      &
    writenc_i1d_real, writenc_i1d_double, writenc_i1d_text,     &
    writenc_i2d_int1, writenc_i2d_int2,   writenc_i2d_int,      &
    writenc_i2d_real, writenc_i2d_double, writenc_i2d_text,     &
    writenc_i3d_int1, writenc_i3d_int2,   writenc_i3d_int,      &
    writenc_i3d_real, writenc_i3d_double, writenc_i3d_text,     &
    writenc_i4d_int1, writenc_i4d_int2,   writenc_i4d_int,      &
    writenc_i4d_real, writenc_i4d_double, writenc_i4d_text,     &
    writenc_n0d_int1, writenc_n0d_int2,   writenc_n0d_int,      &
    writenc_n0d_real, writenc_n0d_double, writenc_n0d_text,     &
    writenc_n1d_int1, writenc_n1d_int2,   writenc_n1d_int,      &
    writenc_n1d_real, writenc_n1d_double, writenc_n1d_text,     &
    writenc_n2d_int1, writenc_n2d_int2,   writenc_n2d_int,      &
    writenc_n2d_real, writenc_n2d_double, writenc_n2d_text,     &
    writenc_n3d_int1, writenc_n3d_int2,   writenc_n3d_int,      &
    writenc_n3d_real, writenc_n3d_double, writenc_n3d_text,     &
    writenc_n4d_int1, writenc_n4d_int2,   writenc_n4d_int,      &
    writenc_n4d_real, writenc_n4d_double, writenc_n4d_text
  end interface writenc
!-------------------------------------------------------------------------------
  interface readnc
    module procedure &
    readnc_i0d_int1, readnc_i0d_int2,   readnc_i0d_int,        &
    readnc_i0d_real, readnc_i0d_double, readnc_i0d_text,       &
    readnc_i1d_int1, readnc_i1d_int2,   readnc_i1d_int,        &
    readnc_i1d_real, readnc_i1d_double, readnc_i1d_text,       &
    readnc_i2d_int1, readnc_i2d_int2,   readnc_i2d_int,        &
    readnc_i2d_real, readnc_i2d_double, readnc_i2d_text,       &
    readnc_i3d_int1, readnc_i3d_int2,   readnc_i3d_int,        &
    readnc_i3d_real, readnc_i3d_double, readnc_i3d_text,       &
    readnc_i4d_int1, readnc_i4d_int2,   readnc_i4d_int,        &
    readnc_i4d_real, readnc_i4d_double, readnc_i4d_text,       &
    readnc_n0d_int1, readnc_n0d_int2,   readnc_n0d_int,        &
    readnc_n0d_real, readnc_n0d_double, readnc_n0d_text,       &
    readnc_n1d_int1, readnc_n1d_int2,   readnc_n1d_int,        &
    readnc_n1d_real, readnc_n1d_double, readnc_n1d_text,       &
    readnc_n2d_int1, readnc_n2d_int2,   readnc_n2d_int,        &
    readnc_n2d_real, readnc_n2d_double, readnc_n2d_text,       &
    readnc_n3d_int1, readnc_n3d_int2,   readnc_n3d_int,        &
    readnc_n3d_real, readnc_n3d_double, readnc_n3d_text,       &
    readnc_n4d_int1, readnc_n4d_int2,   readnc_n4d_int,        &
    readnc_n4d_real, readnc_n4d_double, readnc_n4d_text
  end interface readnc
  character(len=80),private :: version = 'netcdf_mod 2.6,by zhaowei, 20014-3-1 .'
!-------------------------------------------------------------------------------
  contains
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
! --- To handle errors for NCFile input or output.
!* checknc_err
#define Mchecknc_err(stat) checknc_err(stat,__LINE__)
  subroutine checknc_err(status,cline)
  integer, intent(in) :: status,cline
  if (status .ne. NF_NOERR) then
    write(*, *)cline,nf_strerror(status)
    ! call wav_abort( NF_STRERROR(status))
    stop
  endif
  return
  end subroutine checknc_err
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
! --- To open a NcFile. (By create, if new; By redefine, if old.)
!* open_nc_def
  subroutine open_nc(NcID, FileName, action)
  character(len=*), intent(in) :: FileName, action
  integer, intent(inout) :: NcID
  logical alive
  integer status
! SDBGINF,"OPEN_NC "//trim(action)//"  "//trim(FileName)
  select case(action(1:1))
  case('c', 'C')
    status = NF_CREATE(trim(FileName), NF_CLOBBER+nf_share, ncID)
    IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
  case('d', 'D')
    inquire(file=trim(FileName), exist = alive)
    if(.not.alive)then
      write(*, *)trim(FileName), ' for define is not exist!'; stop
    endif
    status = NF_OPEN(trim(FileName), NF_WRITE+nf_share, NcID)
    IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
    status = NF_REDEF(NcID)
    IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
  case('w', 'W')
    inquire(file=trim(FileName), exist = alive)
    if(.not.alive)then
      write(*, *)trim(FileName), ' for write is not exist!'; stop
    endif
    status = NF_OPEN(trim(FileName), NF_WRITE+nf_share, NcID)
    IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
  case('r', 'R')
    inquire(file=trim(FileName), exist = alive)
    if(.not.alive)then
      write(*, *)trim(FileName), ' for read is not exist!'; stop
    endif
    status = NF_OPEN(trim(FileName), NF_NOWRITE+nf_share, NcID)
    IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
  end select
  return
  end subroutine open_nc
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
! --- Close a openning NcFile.
!* WriteNC1DInt2
  subroutine close_nc(NcID)
  integer, intent(inout) :: NcID
  integer status
  CALL Mchecknc_err(NF_CLOSE(NcID))
  NcID=-1
  end subroutine close_nc
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
! --- To end define.
!* end_define
  subroutine end_define(NcID)
  integer, intent(in) :: NcID
  character(len=80) :: lib_version
  character(len=10) :: dd, tt
  character(len=19) :: ts
  integer status
  call date_and_time(dd, tt)
  ts = dd(1:4)//'/'//dd(5:6)//'/'//dd(7:8)//' '//tt(1:2)//':'//tt(3:4)//':'//tt(5:10)
! --- Output sorftware versions.
  status=set_attribute(NcID,nf_global, 'CreatedTime', trim(ts))
  lib_version = 'netcdf '//trim(nf_inq_libvers())  
  status=set_attribute(NcID,nf_global, 'Sorftware1', lib_version)
  status=set_attribute(NcID,nf_global, 'Sorftware2', version)
! --- end of defination.
  status = NF_ENDDEF(NcID)
  IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
  return
  end subroutine end_define
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
! --- To define dimension for a NcFile.
!* dimension_define_notype
  subroutine dimension_define(NcID, DimName, DimLen,DimID, DimVarName, DimVarType, DimVarID)
  character(len=*), intent(in) :: DimName
  integer, intent(in) :: NcId, DimLen
  character(len=*), optional,intent(in) :: DimVarName
  integer, optional,intent(in) :: DimVarType
  integer, optional, intent(out) :: DimID, DimVarID
  integer :: DimID1, DimVarID1
  integer status
  status = NF_DEF_DIM(NcID, DimName, DimLen, DimID1)
  IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
  if(present(DimVarName).and. present(DimVarType))then
    status = nf_def_var(ncid, trim(DimVarName), DimVarType, 1, DimID1, DimVarID1)
    IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
  endif
  if(present(DimID))DimID = DimID1
  if(present(DimVarID))DimVarID = DimVarID1
  return
  end subroutine dimension_define
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
! --- Get the dimension length,
!* get_dimension_len
  integer function get_dimension_len(NcID, DimName)

  integer, intent(in) :: NcID
  character(len=*), intent(in) :: DimName
  integer :: DimID
  integer status
  status = NF_INQ_DIMID(NcID, trim(DimName), DimID)
  IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
  status = NF_INQ_DIMLEN(NcID, DimID, get_dimension_len)
  IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
  return
  end function get_dimension_len
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
! --- To define variables.
!* variable_define1
  subroutine variable_define1(NcID, VarName, VarType, VarDims, VarID)
  integer, intent(in) :: NcID
  character(len=*), intent(in) :: VarName
  integer, intent(in) :: VarType, VarDims(:)
  integer, optional, intent(out) :: VarID
  integer :: VarRank, VarID1
  integer status
  status=NF_INQ_VarID(NcID, trim(VarName), VarID1)
  if(status.ne.NF_NOERR)then
    VarRank = size(VarDims)
    status = NF_DEF_VAR(ncid, trim(VarName), VarType, VarRank, VarDims, VarID1)
    IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
  endif
  if(present(VarID))VarID = VarID1
  return
  end subroutine variable_define1
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
! --- To define variables.
!* variable_define2
  subroutine variable_define2(NcID, VarName, VarType, VarDimsName, VarID)
  integer, intent(in) :: NcID
  character(len=*), intent(in) :: VarName
  integer, intent(in) :: VarType
  character(len=*), intent(in) :: VarDimsName(:)
  integer, optional, intent(out) :: VarID
  integer :: VarRank, VarID1, i
  integer, allocatable :: VarDims(:)
  integer status
  status=NF_INQ_VarID(NcID, trim(VarName), VarID1)
  if(status.ne.NF_NOERR)then
    VarRank = size(VarDimsName); allocate(VarDims(VarRank))
    do i = 1, VarRank
      status = NF_INQ_DIMID(NcID, trim(VarDimsName(i)), VarDims(i))
      IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
    enddo
    status = NF_DEF_VAR(ncid, trim(VarName), VarType, VarRank, VarDims, VarID1)
    IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)
    deallocate(VarDims)
  endif
  if(present(VarID))VarID = VarID1
  return
  end subroutine variable_define2
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
! --- To define attribute: character
!* set_attribute_i_character
  integer function set_attribute_i_TEXT(NcID,VarID_, AttName, attribute)
  integer, intent(in) :: NcID,VarID_
  character(len=*), intent(in) :: AttName
  character(*), intent(in) :: attribute
  integer ::VarID;
  VarID=VarID_
  if(VarID<0)VarID=NF_GLOBAL
  set_attribute_i_TEXT = NF_PUT_ATT_TEXT(NcID, VarID, AttName,len_trim(attribute), trim(attribute))
  end function set_attribute_i_TEXT
  
#define TMPL_set_attribute_i(NA,VT)                                            \
  integer function set_attribute_i_##NA(NcID,VarID_, AttName, attribute)      ;\
  integer, intent(in) :: NcID,VarID_                                          ;\
  character(len=*), intent(in) :: AttName                                     ;\
  VT, intent(in) :: attribute                                                 ;\
  integer ::VarID;VarID=VarID_                                                ;\
  set_attribute_i_##NA = NF_PUT_ATT_##NA(NcID, VarID, AttName,nf_##NA,1, attribute)   ;\
  end function set_attribute_i_##NA
#define TMPL_get_attribute_i(NA,VT)                                            \
  integer function get_attribute_i_##NA(NcID,VarID_, AttName, attribute)      ;\
  integer, intent(in) :: NcID,VarID_                                          ;\
  character(len=*), intent(in) :: AttName                                     ;\
  VT, intent(out) :: attribute                                                 ;\
  integer ::VarID;VarID=VarID_                                                ;\
  get_attribute_i_##NA = NF_GET_ATT_##NA(NcID, VarID, AttName, attribute)     ;\
  end function get_attribute_i_##NA
  TMPL_set_attribute_i(INT1  ,integer(1)  )
  TMPL_set_attribute_i(INT2  ,integer(2)  )
  TMPL_set_attribute_i(INT   ,integer     )
  TMPL_set_attribute_i(REAL  ,real*4      )
  TMPL_set_attribute_i(double,real*8      )
  TMPL_get_attribute_i(TEXT  ,character(*))
  TMPL_get_attribute_i(INT1  ,integer(1)  )
  TMPL_get_attribute_i(INT2  ,integer(2)  )
  TMPL_get_attribute_i(INT   ,integer     )
  TMPL_get_attribute_i(REAL  ,real*4      )
  TMPL_get_attribute_i(double,real*8      )
#define TMPL_set_attribute_n(NA,VT)                                            \
  integer function set_attribute_n_##NA(NcID, AttName, attribute,VarName)     ;\
    integer, intent(in) :: NcID                                               ;\
    character(len=*), intent(in) :: AttName                                   ;\
    character(len=*), optional, intent(in) :: VarName                         ;\
    VT, intent(in) :: attribute                                               ;\
    integer :: VarID;VarID=GetVarid(NcID,VarName)                             ;\
    set_attribute_n_##NA=set_attribute(NcID, VarID, AttName, attribute)       ;\
  end function set_attribute_n_##NA
#define TMPL_get_attribute_n(NA,VT)                                            \
  integer function get_attribute_n_##NA(NcID, AttName, attribute, VarName)    ;\
    integer, intent(in) :: NcID                                               ;\
    character(len=*), intent(in) :: AttName                                   ;\
    character(len=*), optional, intent(in) :: VarName                         ;\
    VT, intent(out) :: attribute                                              ;\
    integer :: VarID;VarID=GetVarid(NcID,VarName)                             ;\
    get_attribute_n_##NA = get_attribute(NcID, VarID, AttName, attribute)     ;\
  end function get_attribute_n_##NA
  TMPL_set_attribute_n(TEXT  ,character(*))
  TMPL_set_attribute_n(INT1  ,integer(1)  )
  TMPL_set_attribute_n(INT2  ,integer(2)  )
  TMPL_set_attribute_n(INT   ,integer     )
  TMPL_set_attribute_n(REAL  ,real*4      )
  TMPL_set_attribute_n(double,real*8      )
  TMPL_get_attribute_n(TEXT  ,character(*))
  TMPL_get_attribute_n(INT1  ,integer(1)  )
  TMPL_get_attribute_n(INT2  ,integer(2)  )
  TMPL_get_attribute_n(INT   ,integer     )
  TMPL_get_attribute_n(REAL  ,real*4      )
  TMPL_get_attribute_n(double,real*8      )

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!* ncdump
  subroutine ncdump(filename)
  character(len=*), intent(in) :: filename
  !call systemqq('ncdump -h '//filename)
  return
  end subroutine ncdump
#define TMPL_readwritenc_i(ND,NA,VT,VSHAPE,VP0)                                \
  subroutine writenc_i##ND##d_##NA(ncid, VarID, Var, RecNum, locs)            ;\
    integer, intent(in) :: ncid,VarID                                         ;\
    VT, intent(in) :: var VSHAPE                                              ;\
    integer, optional, intent(in) :: RecNum,locs(:)                           ;\
    integer ::status,DimNum ,starts(8),counts(8)                              ;\
    DimNum = 1;starts=1;counts=1                                              ;\
    CALL SetVARSC(starts,counts,dimnum,shape(Var),var VP0,RecNum,locs)        ;\
    status = NF_PUT_VARA_##NA(ncID, varID, starts, counts, Var)               ;\
    IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)                          ;\
  end subroutine writenc_i##ND##d_##NA                                        ;\
  subroutine readnc_i##ND##d_##NA(ncid, VarID, Var, RecNum, locs)             ;\
    integer, intent(in) :: ncid,VarID                                         ;\
    VT, intent(out) :: var VSHAPE                                             ;\
    integer, optional, intent(in) :: RecNum,locs(:)                           ;\
    integer ::status,DimNum ,starts(8),counts(8)                              ;\
    DimNum = 1;starts=1;counts=1                                              ;\
    CALL SetVARSC(starts,counts,dimnum,shape(Var),var VP0,RecNum,locs)        ;\
    status = NF_GET_VARA_##NA(ncID, varID, starts, counts, Var)               ;\
    IF(status.NE.NF_NOERR) CALL Mchecknc_err(status)                          ;\
  end subroutine readnc_i##ND##d_##NA
  TMPL_readwritenc_i(0,int1  ,integer(1)  ,         ,         )
  TMPL_readwritenc_i(0,int2  ,integer(2)  ,         ,         )
  TMPL_readwritenc_i(0,int   ,integer(4)  ,         ,         )
  TMPL_readwritenc_i(0,real  ,real(4)     ,         ,         )
  TMPL_readwritenc_i(0,double,real(8)     ,         ,         )
  TMPL_readwritenc_i(0,text  ,character(*),         ,         )
  TMPL_readwritenc_i(1,int1  ,integer(1)  ,(:)      ,(1)      )
  TMPL_readwritenc_i(1,int2  ,integer(2)  ,(:)      ,(1)      )
  TMPL_readwritenc_i(1,int   ,integer(4)  ,(:)      ,(1)      )
  TMPL_readwritenc_i(1,real  ,real(4)     ,(:)      ,(1)      )
  TMPL_readwritenc_i(1,double,real(8)     ,(:)      ,(1)      )
  TMPL_readwritenc_i(1,text  ,character(*),(:)      ,(1)      )
  TMPL_readwritenc_i(2,int1  ,integer(1)  ,(:,:)    ,(1,1)    )
  TMPL_readwritenc_i(2,int2  ,integer(2)  ,(:,:)    ,(1,1)    )
  TMPL_readwritenc_i(2,int   ,integer(4)  ,(:,:)    ,(1,1)    )
  TMPL_readwritenc_i(2,real  ,real(4)     ,(:,:)    ,(1,1)    )
  TMPL_readwritenc_i(2,double,real(8)     ,(:,:)    ,(1,1)    )
  TMPL_readwritenc_i(2,text  ,character(*),(:,:)    ,(1,1)    )
  TMPL_readwritenc_i(3,int1  ,integer(1)  ,(:,:,:)  ,(1,1,1)  )
  TMPL_readwritenc_i(3,int2  ,integer(2)  ,(:,:,:)  ,(1,1,1)  )
  TMPL_readwritenc_i(3,int   ,integer(4)  ,(:,:,:)  ,(1,1,1)  )
  TMPL_readwritenc_i(3,real  ,real(4)     ,(:,:,:)  ,(1,1,1)  )
  TMPL_readwritenc_i(3,double,real(8)     ,(:,:,:)  ,(1,1,1)  )
  TMPL_readwritenc_i(3,text  ,character(*),(:,:,:)  ,(1,1,1)  )
  TMPL_readwritenc_i(4,int1  ,integer(1)  ,(:,:,:,:),(1,1,1,1))
  TMPL_readwritenc_i(4,int2  ,integer(2)  ,(:,:,:,:),(1,1,1,1))
  TMPL_readwritenc_i(4,int   ,integer(4)  ,(:,:,:,:),(1,1,1,1))
  TMPL_readwritenc_i(4,real  ,real(4)     ,(:,:,:,:),(1,1,1,1))
  TMPL_readwritenc_i(4,double,real(8)     ,(:,:,:,:),(1,1,1,1))
  TMPL_readwritenc_i(4,text  ,character(*),(:,:,:,:),(1,1,1,1))
#define TMPL_readwritenc_n(ND,NA,VT,VSHAPE)                                    \
  subroutine writenc_n##ND##d_##NA(ncid, VarName, Var, RecNum,locs)           ;\
    integer, intent(in) :: ncid                                               ;\
    character(len=*), intent(in) :: VarName                                   ;\
    VT, intent(in) :: var VSHAPE                                              ;\
    integer, intent(in), optional :: RecNum,locs(:)                           ;\
    call writenc(ncid, GetVarid(ncID,VarName), Var, RecNum,locs)              ;\
  end subroutine  writenc_n##ND##d_##NA                                       ;\
  subroutine readnc_n##ND##d_##NA(ncid, VarName, Var, RecNum,locs)            ;\
    integer, intent(in) :: NcID                                               ;\
    character(len=*), intent(in) :: VarName                                   ;\
    VT, intent(out) :: var VSHAPE                                             ;\
    integer, optional, intent(in) :: RecNum,locs(:)                           ;\
    call readnc(ncid, GetVarid(ncID,VarName), Var, RecNum,locs)               ;\
  end subroutine readnc_n##ND##d_##NA
  TMPL_readwritenc_n(0,int1  ,integer(1)  ,)
  TMPL_readwritenc_n(0,int2  ,integer(2)  ,)
  TMPL_readwritenc_n(0,int   ,integer(4)  ,)
  TMPL_readwritenc_n(0,real  ,real(4)     ,)
  TMPL_readwritenc_n(0,double,real(8)     ,)
  TMPL_readwritenc_n(0,text  ,character(*),)
  TMPL_readwritenc_n(1,int1  ,integer(1)  ,(:))
  TMPL_readwritenc_n(1,int2  ,integer(2)  ,(:))
  TMPL_readwritenc_n(1,int   ,integer(4)  ,(:))
  TMPL_readwritenc_n(1,real  ,real(4)     ,(:))
  TMPL_readwritenc_n(1,double,real(8)     ,(:))
  TMPL_readwritenc_n(1,text  ,character(*),(:))
  TMPL_readwritenc_n(2,int1  ,integer(1)  ,(:,:))
  TMPL_readwritenc_n(2,int2  ,integer(2)  ,(:,:))
  TMPL_readwritenc_n(2,int   ,integer(4)  ,(:,:))
  TMPL_readwritenc_n(2,real  ,real(4)     ,(:,:))
  TMPL_readwritenc_n(2,double,real(8)     ,(:,:))
  TMPL_readwritenc_n(2,text  ,character(*),(:,:))
  TMPL_readwritenc_n(3,int1  ,integer(1)  ,(:,:,:))
  TMPL_readwritenc_n(3,int2  ,integer(2)  ,(:,:,:))
  TMPL_readwritenc_n(3,int   ,integer(4)  ,(:,:,:))
  TMPL_readwritenc_n(3,real  ,real(4)     ,(:,:,:))
  TMPL_readwritenc_n(3,double,real(8)     ,(:,:,:))
  TMPL_readwritenc_n(3,text  ,character(*),(:,:,:))
  TMPL_readwritenc_n(4,int1  ,integer(1)  ,(:,:,:,:))
  TMPL_readwritenc_n(4,int2  ,integer(2)  ,(:,:,:,:))
  TMPL_readwritenc_n(4,int   ,integer(4)  ,(:,:,:,:))
  TMPL_readwritenc_n(4,real  ,real(4)     ,(:,:,:,:))
  TMPL_readwritenc_n(4,double,real(8)     ,(:,:,:,:))
  TMPL_readwritenc_n(4,text  ,character(*),(:,:,:,:))

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
integer function GetVarid(ncID,vname)
  INTEGER ncID
  character*(*),optional,intent(in)::vname
  if(present(vname)) then
    CALL Mchecknc_err(NF_INQ_VARID(ncID,vname,GetVarid))
  else
    GetVarid=NF_GLOBAL
  endif
end function GetVarid
integer function GetVarType_n(ncID,VName)
  INTEGER ncID
  character*(*) vname
  INTEGER vid,vtype
  CALL Mchecknc_err(NF_INQ_VARID(ncID,vname,vID))
  CALL Mchecknc_err(NF_INQ_VARTYPE(ncID,vID,vtype))
  GetVarType_n=vtype
end function GetVarType_n
integer function GetVarType_i(ncID,Vid)
  INTEGER ncID
  INTEGER vid,vtype
  CALL Mchecknc_err(NF_INQ_VARTYPE(ncID,vID,vtype))
  GetVarType_i=vtype
end function GetVarType_i
INTEGER function nfDefVar(ncid,name,logtype,Sc,Off,ndim,vdim)
  INTEGER ncid,logtype,varid,ndim
    INTEGER vdim(ndim)
  real(8) Sc,Off
  CHARACTER name*(*)
  IF(logtype==0)THEN
    nfDefVar=-1
    return
!  else IF(logtype==1)THEN
!    CALL Mchecknc_err(NF_DEF_VAR(ncID,name,NF_SHORT,ndim,vdim,varid))
!    CALL Mchecknc_err(NF_PUT_ATT_INT2(ncID,varid,'missing_value',NF_SHORT,1,32767))
!    CALL Mchecknc_err(NF_PUT_ATT_real(ncID,varid,'scale_factor',NF_real,1,sc))
!    CALL Mchecknc_err(NF_PUT_ATT_real(ncID,varid,'Offset',NF_real,1,Off))
  else
    CALL Mchecknc_err(NF_DEF_VAR(ncID,name,nf_float,ndim,vdim,varid))
    CALL Mchecknc_err(NF_PUT_ATT_real(ncID,varid,'missing_value',NF_SHORT,1,1.71e38))
  end if
  nfDefVar=varid
END function nfDefVar
  subroutine SetVARSC_TEXT(starts,counts,DimNum,vshape,var,RecNum,locs)
    integer,intent(inout)::DimNum
    integer,intent(out):: starts(:),counts(:)
    integer,intent(in) ::vshape(:)
    CHARACTER(*),intent(in)::var
    integer, optional, intent(in) :: RecNum,locs(:)
    DimNum=size(vshape);starts=1;counts=1;
    if(DimNum>0)counts(1:DimNum)=vshape
    counts(1) = len(var);counts(2:size(vshape)+1)=vshape
    if(present(locs))then
      DimNum=size(locs);
      starts(2:DimNum+1) = locs(1:DimNum)
    endif
    if(present(RecNum))starts(DimNum+2) = RecNum
  end subroutine SetVARSC_TEXT
#define TMPL_SetVARSC(NA,VT)                                                   \
  subroutine SetVARSC_##NA(starts,counts,DimNum,vshape,var,RecNum,locs)       ;\
    integer,intent(inout)::DimNum                                             ;\
    integer,intent(out):: starts(:),counts(:)                                 ;\
    integer,intent(in) ::vshape(:)                                            ;\
    VT,intent(in)::var                                                        ;\
    integer, optional, intent(in) :: RecNum,locs(:)                           ;\
    DimNum=size(vshape);starts=1;counts=1;                                    ;\
    if(DimNum>0)counts(1:DimNum)=vshape                                       ;\
    if(present(locs))then                                                     ;\
      DimNum=size(locs);                                                      ;\
      starts(1:DimNum) = locs(1:DimNum)                                       ;\
    endif                                                                     ;\
    if(present(RecNum))starts(DimNum+1) = RecNum                              ;\
  end subroutine SetVARSC_##NA
  TMPL_SetVARSC(INT1  ,integer(1)  )
  TMPL_SetVARSC(INT2  ,integer(2)  )
  TMPL_SetVARSC(INT   ,integer     )
  TMPL_SetVARSC(REAL  ,real*4      )
  TMPL_SetVARSC(double,real*8      )
end module netcdf_mod
!-------------------------------------------------------------------------------
!###############################################################################
