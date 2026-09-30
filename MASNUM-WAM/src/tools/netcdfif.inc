  !-------------------- Explicit Interfaces for nf routines ------------------
  !================= Misc functions first
  Interface nf_inq_libvers
   Character(LEN=80) Function nf_inq_libvers() 
    End Function nf_inq_libvers
  End Interface
  Interface nf_strerror
    Function nf_strerror(nerr) RESULT(errmsg)
      Integer, Intent(IN) :: nerr
      Character(LEN=80)   :: errmsg
    End Function nf_strerror
  End Interface
  Interface nf_issyserr
    Function nf_issyserr(nerr) RESULT(status)
      Integer, Intent(IN) :: nerr
      Logical             :: status
    End Function nf_issyserr
  End Interface
  !================= Control routines
  Interface nf_create
    Function nf_create(path, cmode, ncid) RESULT (status)
      Character(LEN=*), Intent(IN)  :: path
      Integer,          Intent(IN)  :: cmode
      Integer,          Intent(OUT) :: ncid
      Integer                       :: status
    End Function nf_create
  End Interface
  Interface nf_create_par
    Function nf_create_par (path, cmode, comm, info, ncid) RESULT(status)
      Integer,          Intent(IN)  :: cmode, comm, info
      Character(LEN=*), Intent(IN)  :: path
      Integer,          Intent(OUT) :: ncid
      Integer                       :: status
    End Function nf_create_par
  End Interface
  Interface nf__create
    Function nf__create(path, cmode, initialsz, chunksizehintp, ncid) &
        RESULT(status)
      Character(LEN=*), Intent(IN)  :: path
      Integer,          Intent(IN)  :: cmode, initialsz, chunksizehintp
      Integer,          Intent(OUT) :: ncid
      Integer                       :: status
    End Function nf__create
  End Interface
  Interface nf__create_mp
    Function nf__create_mp(path, cmode, initialsz, basepe, chunksizehintp, ncid) &
        RESULT(status)
      Character(LEN=*), Intent(IN)  :: path
      Integer,          Intent(IN)  :: cmode, initialsz, chunksizehintp, basepe
      Integer,          Intent(OUT) :: ncid
      Integer                       :: status
    End Function nf__create_mp
  End Interface
  Interface nf_open
    Function nf_open(path, mode, ncid) RESULT (status)
      Character(LEN=*), Intent(IN)    :: path
      Integer,          Intent(IN)    :: mode
      Integer,          Intent(INOUT) :: ncid
      Integer                         :: status
    End Function nf_open
  End Interface
  Interface nf_open_par
    Function nf_open_par (path, mode, comm, info, ncid) RESULT(status)
      Integer,          Intent(IN)  :: mode, comm, info
      Character(LEN=*), Intent(IN)  :: path
      Integer,          Intent(OUT) :: ncid
      Integer                       :: status
    End Function nf_open_par
  End Interface
  Interface nf__open
    Function nf__open(path, mode, chunksizehintp, ncid) RESULT (status)
      Character(LEN=*), Intent(IN)    :: path
      Integer,          Intent(IN)    :: mode, chunksizehintp
      Integer,          Intent(INOUT) :: ncid
      Integer                         :: status
    End Function nf__open
  End Interface
  Interface nf__open_mp
    Function nf__open_mp(path, mode, basepe, chunksizehintp, ncid) RESULT (status)
      Character(LEN=*), Intent(IN)    :: path
      Integer,          Intent(IN)    :: mode, chunksizehintp, basepe
      Integer,          Intent(INOUT) :: ncid
      Integer                         :: status
    End Function nf__open_mp
  End Interface
  Interface nf_open_mem
    Function nf_open_mem(path, mode, size, memory, ncid) RESULT (status)
      USE ISO_C_BINDING, ONLY : C_CHAR
      Character(LEN=*),       Intent(IN)           :: path
      Integer,                Intent(IN)           :: mode
      Integer,                Intent(IN)           :: size 
      Character(KIND=C_CHAR), Intent(IN),   TARGET :: memory(*)
      Integer,                Intent(INOUT)        :: ncid
      Integer                               :: status
    End Function nf_open_mem
  End Interface
  Interface nf_var_par_access
    Function nf_var_par_access( ncid, varid, iaccess) RESULT (status)
      Integer, Intent(IN) :: ncid, varid, iaccess
      Integer             :: status
    End Function nf_var_par_access
  End Interface
  Interface nf_inq_path
    Function nf_inq_path(ncid, pathlen, path) RESULT (status)
      Integer,          Intent(IN)    :: ncid
      Integer,          Intent(INOUT) :: pathlen
      Character(LEN=*), Intent(INOUT) :: path
      Integer                         :: status
    End Function nf_inq_path
  End Interface
  Interface nf_set_fill
    Function nf_set_fill(ncid, fillmode, old_mode) RESULT(status)
      Integer, Intent(IN)  :: ncid, fillmode
      Integer, Intent(OUT) :: old_mode
      Integer              :: status
    End Function nf_set_fill
  End Interface
  Interface nf_set_default_format
    Function nf_set_default_format(newform, old_format) RESULT(status)
      Integer, Intent(IN)  :: newform
      Integer, Intent(OUT) :: old_format
      Integer              :: status
    End Function nf_set_default_format
  End Interface
  Interface nf_redef
    Function nf_redef(ncid) RESULT(status)
      Integer, Intent(IN) :: ncid
      Integer             :: status
    End Function nf_redef
  End Interface
  Interface nf_enddef
    Function nf_enddef(ncid) RESULT(status)
      Integer, Intent(IN) :: ncid
      Integer             :: status
    End Function nf_enddef
  End Interface
  Interface nf__enddef
    Function nf__enddef(ncid, h_minfree, v_align, v_minfree, r_align) &
        RESULT(status)
      Integer, Intent(IN) :: ncid, h_minfree, v_align, v_minfree, r_align
      Integer             :: status
    End Function nf__enddef
  End Interface
  Interface nf_sync
    Function nf_sync(ncid) RESULT(status)
      Integer, Intent(IN) :: ncid
      Integer             :: status
    End Function nf_sync
  End Interface
  Interface nf_abort
    Function nf_abort(ncid) RESULT(status)
      Integer, Intent(IN) :: ncid
      Integer             :: status
    End Function nf_abort
  End Interface
  Interface nf_close
    Function nf_close(ncid) RESULT(status)
      Integer, Intent(IN) :: ncid
      Integer             :: status
    End Function nf_close
  End Interface
  Interface nf_delete
    Function nf_delete(path) RESULT(status)
      Character(LEN=*), Intent(IN) :: path
      Integer                      :: status
    End Function nf_delete
  End Interface
  Interface nf_delete_mp
    Function nf_delete_mp(path, pe) RESULT(status)
      Character(LEN=*), Intent(IN) :: path
      Integer,          Intent(IN) :: pe
      Integer                      :: status
    End Function nf_delete_mp
  End Interface
  Interface nf_set_base_pe
    Function nf_set_base_pe(ncid, pe) RESULT(status)
      Integer, Intent(IN) :: ncid, pe
      Integer             :: status
    End Function nf_set_base_pe
  End Interface
  Interface nf_inq_base_pe
    Function nf_inq_base_pe(ncid, pe) RESULT(status)
      Integer, Intent(IN)  :: ncid
      Integer, Intent(OUT) :: pe
      Integer              :: status
    End Function nf_inq_base_pe
  End Interface
  !================= Dimension definition and inquiry functions
  Interface nf_def_dim
    Function nf_def_dim(ncid, name, dlen, dimid) RESULT (status)
      Integer,          Intent(IN)  :: ncid, dlen
      Integer,          Intent(OUT) :: dimid
      Character(LEN=*), Intent(IN)  :: name
      Integer                       :: status
    End Function nf_def_dim
  End Interface
  Interface nf_inq_dim
    Function nf_inq_dim(ncid, dimid, name, dlen) RESULT (status)
      Integer,          Intent(IN)   :: ncid, dimid
      Integer,          Intent(OUT)  :: dlen
      Character(LEN=*), Intent(OUT)  :: name
      Integer                        :: status
    End Function nf_inq_dim
  End Interface
  Interface nf_inq_dimid
    Function nf_inq_dimid(ncid, name, dimid) RESULT (status)
      Integer,          Intent(IN)  :: ncid
      Integer,          Intent(OUT) :: dimid
      Character(LEN=*), Intent(IN)  :: name
      Integer                       :: status
    End Function nf_inq_dimid
  End Interface
  Interface nf_inq_dimlen
    Function nf_inq_dimlen(ncid, dimid, dlen) RESULT (status)
      Integer, Intent(IN)  :: ncid, dimid
      Integer, Intent(OUT) :: dlen
      Integer              :: status
    End Function nf_inq_dimlen
  End Interface
  Interface nf_inq_dimname
    Function nf_inq_dimname (ncid, dimid, name) RESULT (status)
      Integer,          Intent(IN)   :: ncid, dimid
      Character(LEN=*), Intent(OUT)  :: name
      Integer                        :: status
    End Function nf_inq_dimname
  End Interface
  Interface nf_rename_dim
    Function nf_rename_dim(ncid, dimid, name) RESULT (status)
      Integer,          Intent(IN)  :: ncid, dimid
      Character(LEN=*), Intent(IN)  :: name
      Integer                       :: status
    End Function nf_rename_dim
  End Interface
  !================= General inquiry functions
  Interface nf_inq
    Function nf_inq(ncid, ndims, nvars, ngatts, unlimdimid) RESULT(status)
      Integer, Intent(IN)  :: ncid
      Integer, Intent(OUT) :: ndims, nvars, ngatts, unlimdimid
      Integer              :: status
    End Function nf_inq
  End Interface
  Interface nf_inq_ndims
    Function nf_inq_ndims(ncid, ndims) RESULT(status)
      Integer, Intent(IN)  :: ncid
      Integer, Intent(OUT) :: ndims
      Integer              :: status
    End Function nf_inq_ndims
  End Interface
  Interface nf_inq_nvars
    Function nf_inq_nvars(ncid, nvars) RESULT(status)
      Integer, Intent(IN)  :: ncid
      Integer, Intent(OUT) :: nvars
      Integer              :: status
    End Function nf_inq_nvars
  End Interface
  Interface nf_inq_natts
    Function nf_inq_natts(ncid, ngatts) RESULT(status)
      Integer, Intent(IN)  :: ncid
      Integer, Intent(OUT) :: ngatts
      Integer              :: status
    End Function nf_inq_natts
  End Interface
  Interface nf_inq_unlimdim
    Function nf_inq_unlimdim(ncid, unlimdimid) RESULT(status)
      Integer, Intent(IN)  :: ncid
      Integer, Intent(OUT) :: unlimdimid
      Integer              :: status
    End Function nf_inq_unlimdim
  End Interface
  Interface nf_inq_format
    Function nf_inq_format(ncid, format_type) RESULT(status)
      Integer, Intent(IN)  :: ncid
      Integer, Intent(OUT) :: format_type
      Integer              :: status
    End Function nf_inq_format
  End Interface
  !================= General variable functions
  Interface nf_def_var
    Function nf_def_var(ncid, name, xtype, nvdims, vdims, varid) RESULT (status)
      Integer,          Intent(IN)  :: ncid, xtype, nvdims
      Integer,          Intent(IN)  :: vdims(*)
      Integer,          Intent(OUT) :: varid
      Character(LEN=*), Intent(IN)  :: name
      Integer                       :: status
    End Function nf_def_var
  End Interface
  Interface nf_inq_varndims
    Function nf_inq_varndims(ncid, varid, vndims) RESULT (status)
      Integer, Intent(IN)  :: ncid, varid
      Integer, Intent(OUT) :: vndims
      Integer              :: status
    End Function nf_inq_varndims
  End Interface
  Interface nf_inq_var
    Function nf_inq_var(ncid, varid, name, xtype, ndims, dimids, natts) &
        RESULT (status)
      Integer,          Intent(IN)  :: ncid, varid
      Character(LEN=*), Intent(OUT) :: name
      Integer,          Intent(OUT) :: dimids(*)
      Integer,          Intent(OUT) :: ndims, xtype, natts
      Integer                       :: status
    End Function nf_inq_var
  End Interface
  Interface nf_inq_vardimid
    Function nf_inq_vardimid(ncid, varid, dimids) RESULT (status)
      Integer, Intent(IN)  :: ncid, varid
      Integer, Intent(OUT) :: dimids(*)
      Integer              :: status
    End Function nf_inq_vardimid
  End Interface
  Interface nf_inq_varid
    Function nf_inq_varid(ncid, name, varid) RESULT (status)
      Integer,          Intent(IN)  :: ncid
      Integer,          Intent(OUT) :: varid
      Character(LEN=*), Intent(IN)  :: name
      Integer                       :: status
    End Function nf_inq_varid
  End Interface
  Interface nf_inq_varname
    Function nf_inq_varname (ncid, varid, name) RESULT (status)
      Integer,          Intent(IN)   :: ncid, varid
      Character(LEN=*), Intent(OUT)  :: name
      Integer                        :: status
    End Function nf_inq_varname
  End Interface
  Interface nf_inq_vartype
    Function nf_inq_vartype(ncid, varid, xtype) RESULT(status)
      Integer, Intent(IN)  :: ncid, varid
      Integer, Intent(OUT) :: xtype
      Integer              :: status
    End Function nf_inq_vartype
  End Interface
  Interface nf_inq_varnatts
    Function nf_inq_varnatts(ncid, varid, nvatts) RESULT(status)
      Integer, Intent(IN)  :: ncid, varid
      Integer, Intent(OUT) :: nvatts
      Integer              :: status
    End Function nf_inq_varnatts
  End Interface
  Interface nf_rename_var
    Function nf_rename_var(ncid, varid, name) RESULT (status)
      Integer,          Intent(IN) :: ncid, varid
      Character(LEN=*), Intent(IN) :: name
      Integer                      :: status
    End Function nf_rename_var
  End Interface
  Interface nf_copy_var
    Function nf_copy_var(ncid_in, varid, ncid_out) RESULT(status)
      Integer, Intent(IN) :: ncid_in, varid, ncid_out
      Integer             :: status
    End Function nf_copy_var
  End Interface
  !================= General attribute functions
  Interface nf_inq_att
    Function nf_inq_att(ncid, varid, name, xtype, nlen) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Integer,          Intent(OUT) :: nlen, xtype
      Character(LEN=*), Intent(IN)  :: name
      Integer                       :: status
    End Function nf_inq_att
  End Interface
  Interface nf_inq_atttype
    Function nf_inq_atttype(ncid, varid, name, xtype) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Integer,          Intent(OUT) :: xtype
      Character(LEN=*), Intent(IN)  :: name
      Integer                       :: status
    End Function nf_inq_atttype
  End Interface
  Interface nf_inq_attlen
    Function nf_inq_attlen(ncid, varid, name, nlen) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Integer,          Intent(OUT) :: nlen
      Character(LEN=*), Intent(IN)  :: name
      Integer                       :: status
    End Function nf_inq_attlen
  End Interface
  Interface nf_inq_attid
    Function nf_inq_attid(ncid, varid, name, attnum) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Integer,          Intent(OUT) :: attnum
      Character(LEN=*), Intent(IN)  :: name
      Integer                       :: status
    End Function nf_inq_attid
  End Interface
  Interface nf_inq_attname
    Function nf_inq_attname(ncid, varid, attnum, name) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid, attnum
      Character(LEN=*), Intent(OUT) :: name
      Integer                       :: status
    End Function nf_inq_attname
  End Interface
  Interface nf_copy_att
    Function nf_copy_att(ncid_in, varid_in, name, ncid_out, varid_out) &
        RESULT(status)
      Integer,          Intent(IN)  :: ncid_in, varid_in, ncid_out, varid_out
      Character(LEN=*), Intent(IN)  :: name
      Integer                       :: status
    End Function nf_copy_att
  End Interface
  Interface nf_rename_att
    Function nf_rename_att(ncid, varid, name, newname) RESULT(status)
      Integer,          Intent(IN) :: ncid, varid
      Character(LEN=*), Intent(IN) :: name, newname
      Integer                      :: status
    End Function nf_rename_att
  End Interface
  Interface nf_del_att
    Function nf_del_att(ncid, varid, name) RESULT(status)
      Integer,          Intent(IN) :: ncid, varid
      Character(LEN=*), Intent(IN) :: name
      Integer                      :: status
    End Function nf_del_att
  End Interface
  !================= var1 put and get functions
  Interface nf_put_var1_text
    Function nf_put_var1_text(ncid, varid, ndex, chval) RESULT(status)
      Integer,          Intent(IN) :: ncid, varid
      Integer,          Intent(IN) :: ndex(*)
      Character(LEN=1), Intent(IN) :: chval
      Integer                      :: status
    End Function nf_put_var1_text
  End Interface
  Interface nf_put_var1_int1
    Function nf_put_var1_int1(ncid, varid, ndex, ival) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT1
      Integer,         Intent(IN) :: ncid, varid
      Integer,         Intent(IN) :: ndex(*)
      Integer(NFINT1), Intent(IN) :: ival
      Integer                     :: status
    End Function nf_put_var1_int1
  End Interface
  Interface nf_put_var1_int2
    Function nf_put_var1_int2(ncid, varid, ndex, ival) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT2
      Integer,         Intent(IN) :: ncid, varid
      Integer,         Intent(IN) :: ndex(*)
      Integer(NFINT2), Intent(IN) :: ival
      Integer                     :: status
    End Function nf_put_var1_int2
  End Interface
  Interface nf_put_var1_int
    Function nf_put_var1_int(ncid, varid, ndex, ival) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT
      Integer,        Intent(IN) :: ncid, varid
      Integer,        Intent(IN) :: ndex(*)
      Integer(NFINT), Intent(IN) :: ival
      Integer                    :: status
    End Function nf_put_var1_int
  End Interface
  Interface nf_put_var1_real
    Function nf_put_var1_real(ncid, varid, ndex, rval) RESULT(status)
      USE netcdf_nf_data, ONLY: NFREAL
      Integer,      Intent(IN) :: ncid, varid
      Integer,      Intent(IN) :: ndex(*)
      Real(NFREAL), Intent(IN) :: rval
      Integer                  :: status
    End Function nf_put_var1_real
  End Interface
  Interface nf_put_var1_double
    Function nf_put_var1_double(ncid, varid, ndex, dval) RESULT(status)
      USE netcdf_nf_data, ONLY: RK8
      Integer,   Intent(IN) :: ncid, varid
      Integer,   Intent(IN) :: ndex(*)
      Real(RK8), Intent(IN) :: dval
      Integer               :: status
    End Function nf_put_var1_double
  End Interface
  Interface nf_get_var1_text
    Function nf_get_var1_text(ncid, varid, ndex, chval) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Integer,          Intent(IN)  :: ndex(*)
      Character(LEN=1), Intent(OUT) :: chval
      Integer                       :: status
    End Function nf_get_var1_text
  End Interface
  Interface nf_get_var1_int1
    Function nf_get_var1_int1(ncid, varid, ndex, ival) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT1
      Integer,         Intent(IN)  :: ncid, varid
      Integer,         Intent(IN)  :: ndex(*)
      Integer(NFINT1), Intent(OUT) :: ival
      Integer                      :: status
    End Function nf_get_var1_int1
  End Interface
  Interface nf_get_var1_int2
    Function nf_get_var1_int2(ncid, varid, ndex, ival) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT2
      Integer,         Intent(IN)  :: ncid, varid
      Integer,         Intent(IN)  :: ndex(*)
      Integer(NFINT2), Intent(OUT) :: ival
      Integer                      :: status
    End Function nf_get_var1_int2
  End Interface
  Interface nf_get_var1_int
    Function nf_get_var1_int(ncid, varid, ndex, ival) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT
      Integer,        Intent(IN)  :: ncid, varid
      Integer,        Intent(IN)  :: ndex(*)
      Integer(NFINT), Intent(OUT) :: ival
      Integer                     :: status
    End Function nf_get_var1_int
  End Interface
  Interface nf_get_var1_real
    Function nf_get_var1_real(ncid, varid, ndex, rval) RESULT(status)
      USE netcdf_nf_data, ONLY: NFREAL
      Integer,      Intent(IN)  :: ncid, varid
      Integer,      Intent(IN)  :: ndex(*)
      Real(NFREAL), Intent(OUT) :: rval
      Integer                   :: status
    End Function nf_get_var1_real
  End Interface
  Interface nf_get_var1_double
    Function nf_get_var1_double(ncid, varid, ndex, rval) RESULT(status)
      USE netcdf_nf_data, ONLY: RK8
      Integer,   Intent(IN)  :: ncid, varid
      Integer,   Intent(IN)  :: ndex(*)
      Real(RK8), Intent(OUT) :: rval
      Integer                :: status
    End Function nf_get_var1_double
  End Interface
  !================= var put and get functions
  Interface nf_put_var_text
    Function nf_put_var_text(ncid, varid, text) RESULT(status)
      Integer,          Intent(IN) :: ncid, varid
      Character(LEN=*), Intent(IN) :: text
      Integer                      :: status
    End Function nf_put_var_text
    !================= Array of characters
    Function nf_put_var_text_a(ncid, varid, text) RESULT(status)
      Integer,          Intent(IN) :: ncid, varid
      Character(LEN=1), Intent(IN) :: text(*)
      Integer                      :: status
    End Function nf_put_var_text_a
  End Interface
  Interface nf_put_var_int1
    Function nf_put_var_int1(ncid, varid, i1vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT1
      Integer,         Intent(IN) :: ncid, varid
      Integer(NFINT1), Intent(IN) :: i1vals(*)
      Integer                     :: status
    End Function nf_put_var_int1
  End Interface
  Interface nf_put_var_int2
    Function nf_put_var_int2(ncid, varid, i2vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT2
      Integer,         Intent(IN) :: ncid, varid
      Integer(NFINT2), Intent(IN) :: i2vals(*)
      Integer                     :: status
    End Function nf_put_var_int2
  End Interface
  Interface nf_put_var_int
    Function nf_put_var_int(ncid, varid, ivals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT
      Integer,        Intent(IN) :: ncid, varid
      Integer(NFINT), Intent(IN) :: ivals(*)
      Integer             :: status
    End Function nf_put_var_int
  End Interface
  Interface nf_put_var_real
    Function nf_put_var_real(ncid, varid, rvals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFREAL
      Integer,      Intent(IN) :: ncid, varid
      Real(NFREAL), Intent(IN) :: rvals(*)
      Integer                  :: status
    End Function nf_put_var_real
  End Interface
  Interface nf_put_var_double
    Function nf_put_var_double(ncid, varid, dvals) RESULT(status)
      USE netcdf_nf_data, ONLY: RK8
      Integer,   Intent(IN) :: ncid, varid
      Real(RK8), Intent(IN) :: dvals(*)
      Integer               :: status
    End Function nf_put_var_double
  End Interface
  Interface nf_get_var_text
    Function nf_get_var_text(ncid, varid, text) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Character(LEN=*), Intent(OUT) :: text
      Integer                       :: status
    End Function nf_get_var_text
    !================= array of characters
    Function nf_get_var_text_a(ncid, varid, text) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Character(LEN=1), Intent(OUT) :: text(*)
      Integer                       :: status
    End Function nf_get_var_text_a
  End Interface
  Interface nf_get_var_int1
    Function nf_get_var_int1(ncid, varid, i1vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT1
      Integer,         Intent(IN)  :: ncid, varid
      Integer(NFINT1), Intent(OUT) :: i1vals(*)
      Integer                      :: status
    End Function nf_get_var_int1
  End Interface
  Interface nf_get_var_int2
    Function nf_get_var_int2(ncid, varid, i2vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT2
      Integer,         Intent(IN)  :: ncid, varid
      Integer(NFINT2), Intent(OUT) :: i2vals(*)
      Integer                      :: status
    End Function nf_get_var_int2
  End Interface
  Interface nf_get_var_int
    Function nf_get_var_int(ncid, varid, ivals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT
      Integer,        Intent(IN)  :: ncid, varid
      Integer(NFINT), Intent(OUT) :: ivals(*)
      Integer                     :: status
    End Function nf_get_var_int
  End Interface
  Interface nf_get_var_real
    Function nf_get_var_real(ncid, varid, rvals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFREAL
      Integer,      Intent(IN)  :: ncid, varid
      Real(NFREAL), Intent(OUT) :: rvals(*)
      Integer                   :: status
    End Function nf_get_var_real
  End Interface
  Interface nf_get_var_double
    Function nf_get_var_double(ncid, varid, dvals) RESULT(status)
      USE netcdf_nf_data, ONLY: RK8
      Integer,   Intent(IN)  :: ncid, varid
      Real(RK8), Intent(OUT) :: dvals(*)
      Integer                :: status
    End Function nf_get_var_double
  End Interface
  !================= vars put and get functions
  Interface nf_put_vars_text
    Function nf_put_vars_text(ncid, varid, start, counts, strides, text) &
        RESULT(status)
      Integer,          Intent(IN) :: ncid, varid
      Integer,          Intent(IN) :: start(*), counts(*), strides(*)
      Character(LEN=*), Intent(IN) :: text
      Integer                      :: status
    End Function nf_put_vars_text
    !================= array of characters
    Function nf_put_vars_text_a(ncid, varid, start, counts, strides, text) &
        RESULT(status)
      Integer,          Intent(IN) :: ncid, varid
      Integer,          Intent(IN) :: start(*), counts(*), strides(*)
      Character(LEN=1), Intent(IN) :: text(*)
      Integer                      :: status
    End Function nf_put_vars_text_a
  End Interface
  Interface nf_put_vars_int1
    Function nf_put_vars_int1(ncid, varid, start, counts, strides, i1vals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT1
      Integer,         Intent(IN) :: ncid, varid
      Integer,         Intent(IN) :: start(*), counts(*), strides(*)
      Integer(NFINT1), Intent(IN) :: i1vals(*)
      Integer                     :: status
    End Function nf_put_vars_int1
  End Interface
  Interface nf_put_vars_int2
    Function nf_put_vars_int2(ncid, varid, start, counts, strides, i2vals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT2
      Integer,         Intent(IN) :: ncid, varid
      Integer,         Intent(IN) :: start(*), counts(*), strides(*)
      Integer(NFINT2), Intent(IN) :: i2vals(*)
      Integer                     :: status
    End Function nf_put_vars_int2
  End Interface
  Interface nf_put_vars_int
    Function nf_put_vars_int(ncid, varid, start, counts, strides, ivals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT
      Integer,        Intent(IN) :: ncid, varid
      Integer,        Intent(IN) :: start(*), counts(*), strides(*)
      Integer(NFINT), Intent(IN) :: ivals(*)
      Integer                    :: status
    End Function nf_put_vars_int
  End Interface
  Interface nf_put_vars_real
    Function nf_put_vars_real(ncid, varid, start, counts, strides, rvals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFREAL
      Integer,      Intent(IN) :: ncid, varid
      Integer,      Intent(IN) :: start(*), counts(*), strides(*)
      Real(NFREAL), Intent(IN) :: rvals(*)
      Integer                  :: status
    End Function nf_put_vars_real
  End Interface
  Interface nf_put_vars_double
    Function nf_put_vars_double(ncid, varid, start, counts, strides, dvals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: RK8
      Integer,   Intent(IN) :: ncid, varid
      Integer,   Intent(IN) :: start(*), counts(*), strides(*)
      Real(RK8), Intent(IN) :: dvals(*)
      Integer               :: status
    End Function nf_put_vars_double
  End Interface
  Interface nf_get_vars_text
    Function nf_get_vars_text(ncid, varid, start, counts, strides, text) &
        RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Integer,          Intent(IN)  :: start(*), counts(*), strides(*)
      Character(LEN=*), Intent(OUT) :: text
      Integer                       :: status
    End Function nf_get_vars_text
    !================= array of characters
    Function nf_get_vars_text_a(ncid, varid, start, counts, strides, text) &
        RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Integer,          Intent(IN)  :: start(*), counts(*), strides(*)
      Character(LEN=1), Intent(OUT) :: text(*)
      Integer                       :: status
    End Function nf_get_vars_text_a
  End Interface
  Interface nf_get_vars_int1
    Function nf_get_vars_int1(ncid, varid, start, counts, strides, i1vals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT1
      Integer,         Intent(IN)  :: ncid, varid
      Integer,         Intent(IN)  :: start(*), counts(*), strides(*)
      Integer(NFINT1), Intent(OUT) :: i1vals(*)
      Integer                      :: status
    End Function nf_get_vars_int1
  End Interface
  Interface nf_get_vars_int2
    Function nf_get_vars_int2(ncid, varid, start, counts, strides, i2vals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT2
      Integer,         Intent(IN)  :: ncid, varid
      Integer,         Intent(IN)  :: start(*), counts(*), strides(*)
      Integer(NFINT2), Intent(OUT) :: i2vals(*)
      Integer                      :: status
    End Function nf_get_vars_int2
  End Interface
  Interface nf_get_vars_int
    Function nf_get_vars_int(ncid, varid, start, counts, strides, ivals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT
      Integer,        Intent(IN)  :: ncid, varid
      Integer,        Intent(IN)  :: start(*), counts(*), strides(*)
      Integer(NFINT), Intent(OUT) :: ivals(*)
      Integer                     :: status
    End Function nf_get_vars_int
  End Interface
  Interface nf_get_vars_real
    Function nf_get_vars_real(ncid, varid, start, counts, strides, rvals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFREAL
      Integer,      Intent(IN)  :: ncid, varid
      Integer,      Intent(IN)  :: start(*), counts(*), strides(*)
      Real(NFREAL), Intent(OUT) :: rvals(*)
      Integer                   :: status
    End Function nf_get_vars_real
  End Interface
  Interface nf_get_vars_double
    Function nf_get_vars_double(ncid, varid, start, counts, strides, dvals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: RK8
      Integer,   Intent(IN)  :: ncid, varid
      Integer,   Intent(IN)  :: start(*), counts(*), strides(*)
      Real(RK8), Intent(OUT) :: dvals(*)
      Integer                :: status
    End Function nf_get_vars_double
  End Interface
  !================= varm put and get functions
  Interface nf_put_varm_text
    Function nf_put_varm_text(ncid, varid, start, counts, strides, maps, &
        text) RESULT(status)
      Integer,          Intent(IN) :: ncid, varid
      Integer,          Intent(IN) :: start(*), counts(*), strides(*), maps(*)
      Character(LEN=*), Intent(IN) :: text
      Integer                      :: status
    End Function nf_put_varm_text
    !================= array of characters
    Function nf_put_varm_text_a(ncid, varid, start, counts, strides, maps, &
        text) RESULT(status)
      Integer,          Intent(IN) :: ncid, varid
      Integer,          Intent(IN) :: start(*), counts(*), strides(*), maps(*)
      Character(LEN=1), Intent(IN) :: text(*)
      Integer                      :: status
    End Function nf_put_varm_text_a
  End Interface
  Interface nf_put_varm_int1
    Function nf_put_varm_int1(ncid, varid, start, counts, strides, maps, &
        i1vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT1
      Integer,         Intent(IN) :: ncid, varid
      Integer,         Intent(IN) :: start(*), counts(*), strides(*), maps(*)
      Integer(NFINT1), Intent(IN) :: i1vals(*)
      Integer                     :: status
    End Function nf_put_varm_int1
  End Interface
  Interface nf_put_varm_int2
    Function nf_put_varm_int2(ncid, varid, start, counts, strides, maps, &
        i2vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT2
      Integer,         Intent(IN) :: ncid, varid
      Integer,         Intent(IN) :: start(*), counts(*), strides(*), maps(*)
      Integer(NFINT2), Intent(IN) :: i2vals(*)
      Integer                     :: status
    End Function nf_put_varm_int2
  End Interface
  Interface nf_put_varm_int
    Function nf_put_varm_int(ncid, varid, start, counts, strides, maps, &
        ivals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT
      Integer,        Intent(IN) :: ncid, varid
      Integer,        Intent(IN) :: start(*), counts(*), strides(*), maps(*)
      Integer(NFINT), Intent(IN) :: ivals(*)
      Integer                    :: status
    End Function nf_put_varm_int
  End Interface
  Interface nf_put_varm_real
    Function nf_put_varm_real(ncid, varid, start, counts, strides, maps, &
        rvals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFREAL
      Integer,      Intent(IN) :: ncid, varid
      Integer,      Intent(IN) :: start(*), counts(*), strides(*), maps(*)
      Real(NFREAL), Intent(IN) :: rvals(*)
      Integer                  :: status
    End Function nf_put_varm_real
  End Interface
  Interface nf_put_varm_double
    Function nf_put_varm_double(ncid, varid, start, counts, strides, maps, &
        dvals) RESULT(status)
      USE netcdf_nf_data, ONLY: RK8
      Integer,   Intent(IN) :: ncid, varid
      Integer,   Intent(IN) :: start(*), counts(*), strides(*), maps(*)
      Real(RK8), Intent(IN) :: dvals(*)
      Integer               :: status
    End Function nf_put_varm_double
  End Interface
  Interface nf_get_varm_text
    Function nf_get_varm_text(ncid, varid, start, counts, strides, maps, &
        text) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Integer,          Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
      Character(LEN=*), Intent(OUT) :: text
      Integer                       :: status
    End Function nf_get_varm_text
    !================= array of characters
    Function nf_get_varm_text_a(ncid, varid, start, counts, strides, maps, &
        text) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Integer,          Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
      Character(LEN=1), Intent(OUT) :: text(*)
      Integer                       :: status
    End Function nf_get_varm_text_a
  End Interface
  Interface nf_get_varm_int1
    Function nf_get_varm_int1(ncid, varid, start, counts, strides, maps, &
        i1vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT1
      Integer,         Intent(IN)  :: ncid, varid
      Integer,         Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
      Integer(NFINT1), Intent(OUT) :: i1vals(*)
      Integer                      :: status
    End Function nf_get_varm_int1
  End Interface
  Interface nf_get_varm_int2
    Function nf_get_varm_int2(ncid, varid, start, counts, strides, maps, &
        i2vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT2
      Integer,         Intent(IN)  :: ncid, varid
      Integer,         Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
      Integer(NFINT2), Intent(OUT) :: i2vals(*)
      Integer                      :: status
    End Function nf_get_varm_int2
  End Interface
  Interface nf_get_varm_int
    Function nf_get_varm_int(ncid, varid, start, counts, strides, maps, &
        ivals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT
      Integer,        Intent(IN)  :: ncid, varid
      Integer,        Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
      Integer(NFINT), Intent(OUT) :: ivals(*)
      Integer                     :: status
    End Function nf_get_varm_int
  End Interface
  Interface nf_get_varm_real
    Function nf_get_varm_real(ncid, varid, start, counts, strides, maps, &
        rvals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFREAL
      Integer,      Intent(IN)  :: ncid, varid
      Integer,      Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
      Real(NFREAL), Intent(OUT) :: rvals(*)
      Integer                   :: status
    End Function nf_get_varm_real
  End Interface
  Interface nf_get_varm_double
    Function nf_get_varm_double(ncid, varid, start, counts, strides, maps, &
        dvals) RESULT(status)
      USE netcdf_nf_data, ONLY: RK8
      Integer,   Intent(IN)  :: ncid, varid
      Integer,   Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
      Real(RK8), Intent(OUT) :: dvals(*)
      Integer                :: status
    End Function nf_get_varm_double
  End Interface
  !================= vara put and get routines
  Interface nf_put_vara_text
    Function nf_put_vara_text(ncid, varid, start, counts, text) RESULT(status)
      Integer,          Intent(IN) :: ncid, varid
      Integer,          Intent(IN) :: start(*), counts(*)
      Character(LEN=*), Intent(IN) :: text(*)
      Integer                      :: status
    End Function nf_put_vara_text
  End Interface
  Interface nf_put_vara_int1
    Function nf_put_vara_int1(ncid, varid, start, counts, i1vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT1
      Integer,         Intent(IN) :: ncid, varid
      Integer,         Intent(IN) :: start(*), counts(*)
      Integer(NFINT1), Intent(IN) :: i1vals(*)
      Integer                     :: status
    End Function nf_put_vara_int1
  End Interface
  Interface nf_put_vara_int2
    Function nf_put_vara_int2(ncid, varid, start, counts, i2vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT2
      Integer,         Intent(IN) :: ncid, varid
      Integer,         Intent(IN) :: start(*), counts(*)
      Integer(NFINT2), Intent(IN) :: i2vals(*)
      Integer                     :: status
    End Function nf_put_vara_int2
  End Interface
  Interface nf_put_vara_int
    Function nf_put_vara_int(ncid, varid, start, counts, ivals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT
      Integer,        Intent(IN) :: ncid, varid
      Integer,        Intent(IN) :: start(*), counts(*)
      Integer(NFINT), Intent(IN) :: ivals(*)
      Integer                    :: status
    End Function nf_put_vara_int
  End Interface
  Interface nf_put_vara_int64
    Function nf_put_vara_int64(ncid, varid, start, counts, i8vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT8
      Integer,         Intent(IN) :: ncid, varid
      Integer,         Intent(IN) :: start(*), counts(*)
      Integer(NFINT8), Intent(IN) :: i8vals(*)
      Integer                     :: status
    End Function nf_put_vara_int64
  End Interface
  Interface nf_put_vara_real
    Function nf_put_vara_real(ncid, varid, start, counts, rvals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFREAL
      Integer,      Intent(IN) :: ncid, varid
      Integer,      Intent(IN) :: start(*), counts(*)
      Real(NFREAL), Intent(IN) :: rvals(*)
      Integer                  :: status
    End Function nf_put_vara_real
  End Interface
  Interface nf_put_vara_double
    Function nf_put_vara_double(ncid, varid, start, counts, dvals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: RK8
      Integer,   Intent(IN) :: ncid, varid
      Integer,   Intent(IN) :: start(*), counts(*)
      Real(RK8), Intent(IN) :: dvals(*)
      Integer               :: status
    End Function nf_put_vara_double
  End Interface
  Interface nf_get_vara_text
    Function nf_get_vara_text(ncid, varid, start, counts, text) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Integer,          Intent(IN)  :: start(*), counts(*)
      Character(LEN=*), Intent(OUT) :: text(*)
      Integer                       :: status
    End Function nf_get_vara_text
  End Interface
  Interface nf_get_vara_int1
    Function nf_get_vara_int1(ncid, varid, start, counts, i1vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT1
      Integer,         Intent(IN)  :: ncid, varid
      Integer,         Intent(IN)  :: start(*), counts(*)
      Integer(NFINT1), Intent(OUT) :: i1vals(*)
      Integer                      :: status
    End Function nf_get_vara_int1
  End Interface
  Interface nf_get_vara_int2
    Function nf_get_vara_int2(ncid, varid, start, counts, i2vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT2
      Integer,         Intent(IN)  :: ncid, varid
      Integer,         Intent(IN)  :: start(*), counts(*)
      Integer(NFINT2), Intent(OUT) :: i2vals(*)
      Integer                      :: status
    End Function nf_get_vara_int2
  End Interface
  Interface nf_get_vara_int
    Function nf_get_vara_int(ncid, varid, start, counts, ivals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT
      Integer,        Intent(IN)  :: ncid, varid
      Integer,        Intent(IN)  :: start(*), counts(*)
      Integer(NFINT), Intent(OUT) :: ivals(*)
      Integer                     :: status
    End Function nf_get_vara_int
  End Interface
  Interface nf_get_vara_int64
    Function nf_get_vara_int64(ncid, varid, start, counts, i8vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT8
      Integer,         Intent(IN)  :: ncid, varid
      Integer,         Intent(IN)  :: start(*), counts(*)
      Integer(NFINT8), Intent(OUT) :: i8vals(*)
      Integer                      :: status
    End Function nf_get_vara_int64
  End Interface
  Interface nf_get_vara_real
    Function nf_get_vara_real(ncid, varid, start, counts, rvals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFREAL
      Integer,      Intent(IN)  :: ncid, varid
      Integer,      Intent(IN)  :: start(*), counts(*)
      Real(NFREAL), Intent(OUT) :: rvals(*)
      Integer                   :: status
    End Function nf_get_vara_real
  End Interface
  Interface nf_get_vara_double
    Function nf_get_vara_double(ncid, varid, start, counts, dvals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: RK8
      Integer,   Intent(IN)  :: ncid, varid
      Integer,   Intent(IN)  :: start(*), counts(*)
      Real(RK8), Intent(OUT) :: dvals(*)
      Integer                :: status
    End Function nf_get_vara_double
  End Interface
  Interface nf_put_att_text
    Function nf_put_att_text(ncid, varid, name, nlen, text) RESULT(status)
      Integer,          Intent(IN) :: ncid, varid, nlen
      Character(LEN=*), Intent(IN) :: name, text
      Integer                      :: status
    End Function nf_put_att_text
    !================= array of characters
    Function nf_put_att_text_a(ncid, varid, name, nlen, text) RESULT(status)
      Integer,          Intent(IN) :: ncid, varid, nlen
      Character(LEN=*), Intent(IN) :: name
      Character(LEN=1), Intent(IN) :: text(*)
      Integer                      :: status
    End Function nf_put_att_text_a
  End Interface
  Interface nf_put_att_int1
    Function nf_put_att_int1(ncid, varid, name, xtype, nlen, i1vals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT1
      Integer,          Intent(IN) :: ncid, varid, nlen, xtype
      Character(LEN=*), Intent(IN) :: name
      Integer(NFINT1),  Intent(IN) :: i1vals(*)
      Integer                      :: status
    End Function nf_put_att_int1
  End Interface
  Interface nf_put_att_int2
    Function nf_put_att_int2(ncid, varid, name, xtype, nlen, i2vals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT2
      Integer,          Intent(IN) :: ncid, varid, nlen, xtype
      Character(LEN=*), Intent(IN) :: name
      Integer(NFINT2),  Intent(IN) :: i2vals(*)
      Integer                      :: status
    End Function nf_put_att_int2
  End Interface
  Interface nf_put_att_int
    Function nf_put_att_int(ncid, varid, name, xtype, nlen, ivals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT
      Integer,          Intent(IN) :: ncid, varid, nlen, xtype
      Character(LEN=*), Intent(IN) :: name
      Integer(NFINT),   Intent(IN) :: ivals(*)
      Integer                      :: status
    End Function nf_put_att_int
  End Interface
  Interface nf_put_att_int64
    Function nf_put_att_int64(ncid, varid, name, xtype, nlen, i8vals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT8
      Integer,          Intent(IN) :: ncid, varid, nlen, xtype
      Character(LEN=*), Intent(IN) :: name
      Integer(NFINT8),  Intent(IN) :: i8vals(*)
      Integer                      :: status
    End Function nf_put_att_int64
  End Interface
  Interface nf_put_att_real
    Function nf_put_att_real(ncid, varid, name, xtype, nlen, rvals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: NFREAL
      Integer,          Intent(IN) :: ncid, varid, nlen, xtype
      Character(LEN=*), Intent(IN) :: name
      Real(NFREAL),     Intent(IN) :: rvals(*)
      Integer                      :: status
    End Function nf_put_att_real
  End Interface
  Interface nf_put_att_double
    Function nf_put_att_double(ncid, varid, name, xtype, nlen, dvals) &
        RESULT(status)
      USE netcdf_nf_data, ONLY: RK8
      Integer,          Intent(IN) :: ncid, varid, nlen, xtype
      Character(LEN=*), Intent(IN) :: name
      Real(RK8),        Intent(IN) :: dvals(*)
      Integer                      :: status
    End Function nf_put_att_double
  End Interface
  Interface nf_get_att_text
    Function nf_get_att_text(ncid, varid, name, text) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Character(LEN=*), Intent(IN)  :: name
      Character(LEN=*), Intent(OUT) :: text
      Integer                       :: status
    End Function nf_get_att_text
    !================= Array of characters
    Function nf_get_att_text_a(ncid, varid, name, text) RESULT(status)
      Integer,          Intent(IN)  :: ncid, varid
      Character(LEN=*), Intent(IN)  :: name
      Character(LEN=1), Intent(OUT) :: text(*)
      Integer                       :: status
    End Function nf_get_att_text_a
  End Interface
  Interface nf_get_att_int1
    Function nf_get_att_int1(ncid, varid, name, i1vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT1
      Integer,          Intent(IN)  :: ncid, varid
      Character(LEN=*), Intent(IN)  :: name
      Integer(NFINT1),  Intent(OUT) :: i1vals(*)
      Integer                       :: status
    End Function nf_get_att_int1
  End Interface
  Interface nf_get_att_int2
    Function nf_get_att_int2(ncid, varid, name, i2vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT2
      Integer,          Intent(IN)  :: ncid, varid
      Character(LEN=*), Intent(IN)  :: name
      Integer(NFINT2),  Intent(OUT) :: i2vals(*)
      Integer                       :: status
    End Function nf_get_att_int2
  End Interface
  Interface nf_get_att_int
    Function nf_get_att_int(ncid, varid, name, ivals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT
      Integer,          Intent(IN)  :: ncid, varid
      Character(LEN=*), Intent(IN)  :: name
      Integer(NFINT),   Intent(OUT) :: ivals(*)
      Integer                       :: status
    End Function nf_get_att_int
  End Interface
  Interface nf_get_att_int64
    Function nf_get_att_int64(ncid, varid, name, i8vals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFINT8
      Integer,          Intent(IN)  :: ncid, varid
      Character(LEN=*), Intent(IN)  :: name
      Integer(NFINT8),  Intent(OUT) :: i8vals(*)
      Integer                       :: status
    End Function nf_get_att_int64
  End Interface
  Interface nf_get_att_real
    Function nf_get_att_real(ncid, varid, name, rvals) RESULT(status)
      USE netcdf_nf_data, ONLY: NFREAL
      Integer,          Intent(IN)  :: ncid, varid
      Character(LEN=*), Intent(IN)  :: name
      Real(NFREAL),     Intent(OUT) :: rvals(*)
      Integer                       :: status
    End Function nf_get_att_real
  End Interface
  Interface nf_get_att_double
    Function nf_get_att_double(ncid, varid, name, dvals) RESULT(status)
      USE netcdf_nf_data, ONLY: RK8
      Integer,          Intent(IN)  :: ncid, varid
      Character(LEN=*), Intent(IN)  :: name
      Real(RK8),        Intent(OUT) :: dvals(*)
      Integer                       :: status
    End Function nf_get_att_double
  End Interface
  !======Externals for functions that use C_CHAR arrays to pass data to void pointers
  Integer, External :: nf_put_var1
  Integer, External :: nf_get_var1
  Integer, External :: nf_put_vars
  Integer, External :: nf_get_vars
  Integer, External :: nf_put_vara
  Integer, External :: nf_get_vara
  !================= Integer, External :: nf_open_mem
