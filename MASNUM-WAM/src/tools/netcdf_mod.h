
typedef signed char int1;
typedef short       int2;
typedef int         int4;
typedef float       real4;
typedef double      real8;

typedef enum {
	NC_NAT    =	0,	/* NAT = 'Not A Type' (c.f. NaN) */
	NC_BYTE   =	1,	/* signed 1 byte integer */
	NC_CHAR   =	2,	/* ISO/ASCII character */
	NC_SHORT  =	3,	/* signed 2 byte integer */
	NC_INT    =	4,	/* signed 4 byte integer */
	NC_FLOAT  =	5,	/* single precision floating point number */
	NC_DOUBLE =	6,	/* double precision floating point number */
	NC_INT1   =	2,	/* ISO/ASCII character */
	NC_INT2   =	3,	/* signed 2 byte integer */
	NC_INT4   =	4,	/* signed 4 byte integer */
	NC_REAL4  =	5,	/* single precision floating point number */
	NC_REAL8  =	6	  /* double precision floating point number */
} nc_type;
/*
#define NC_FILL_BYTE	((signed char)-127)
#define NC_FILL_CHAR	((char)0)
#define NC_FILL_SHORT	((short)-32767)
#define NC_FILL_INT	(-2147483647L)
#define NC_FILL_FLOAT	(9.9692099683868690e+36f) /* near 15 * 2^119 */
#define NC_FILL_DOUBLE	(9.9692099683868690e+36)
#define _FillValue	"_FillValue"
#define NC_FILL		0	/* argument to ncsetfill to clear NC_NOFILL */
#define NC_NOFILL	0x100	/* Don't fill data section an records */
#define NC_NOWRITE	0	/* default is read only */
#define NC_WRITE    	0x1	/* read & write */
#define NC_CLOBBER	0
#define NC_NOCLOBBER	0x4	/* Don't destroy existing file on create */
#define NC_64BIT_OFFSET 0x0200  /* Use large (64-bit) file offsets */
#define NC_SHARE	0x0800	/* Share updates, limit cacheing */
#define NC_STRICT_NC3  (0x8)
#define NC_LOCK		0x0400	/* Use locking if available */
#define NC_FORMAT_CLASSIC (1)
#define NC_FORMAT_64BIT   (2)
#define NC_FORMAT_NETCDF4 (3)
#define NC_FORMAT_NETCDF4_CLASSIC  (4) /* create netcdf-4 files, with NC_STRICT_NC3. */
#define NC_SIZEHINT_DEFAULT 0
#define NC_ALIGN_CHUNK ((size_t)(-1))
#define NC_UNLIMITED 0L
#define NC_GLOBAL -1
#define NC_NOAXISTYPE 0
#define NC_LATITUDE 1
#define NC_LONGITUDE 2
#define NC_GEOX 3
#define NC_GEOY 4
#define NC_GEOZ 5
#define NC_HEIGHT_UP 6
#define NC_HEIGHT_DOWN 7
#define NC_PRESSURE 8
#define NC_TIME 9
#define NC_RADAZ 10
#define NC_RADEL 11
#define NC_RADDIST 12
#define NC_MAX_DIMS	1024	 /* max dimensions per file */
#define NC_MAX_ATTRS	8192	 /* max global or per variable attributes */
#define NC_MAX_VARS	8192	 /* max variables per file */
#define NC_MAX_NAME	256	 /* max length of a name */
#define NC_MAX_VAR_DIMS	NC_MAX_DIMS /* max per variable dimensions */
#define NC_ISSYSERR(err)	((err) > 0)
#define	NC_NOERR	0	/* No Error */
#define NC2_ERR         (-1)    /* Returned for all errors in the v2 API. */
#define	NC_EBADID	(-33)	/* Not a netcdf id */
#define	NC_ENFILE	(-34)	/* Too many netcdfs open */
#define	NC_EEXIST	(-35)	/* netcdf file exists && NC_NOCLOBBER */
#define	NC_EINVAL	(-36)	/* Invalid Argument */
#define	NC_EPERM	(-37)	/* Write to read only */
#define	NC_ENOTINDEFINE	(-38)	/* Operation not allowed in data mode */
#define	NC_EINDEFINE	(-39)	/* Operation not allowed in define mode */
#define	NC_EINVALCOORDS	(-40)	/* Index exceeds dimension bound */
#define	NC_EMAXDIMS	(-41)	/* NC_MAX_DIMS exceeded */
#define	NC_ENAMEINUSE	(-42)	/* String match to name in use */
#define NC_ENOTATT	(-43)	/* Attribute not found */
#define	NC_EMAXATTS	(-44)	/* NC_MAX_ATTRS exceeded */
#define NC_EBADTYPE	(-45)	/* Not a netcdf data type */
#define NC_EBADDIM	(-46)	/* Invalid dimension id or name */
#define NC_EUNLIMPOS	(-47)	/* NC_UNLIMITED in the wrong index */
#define	NC_EMAXVARS	(-48)	/* NC_MAX_VARS exceeded */
#define NC_ENOTVAR	(-49)	/* Variable not found */
#define NC_EGLOBAL	(-50)	/* Action prohibited on NC_GLOBAL varid */
#define NC_ENOTNC	(-51)	/* Not a netcdf file */
#define NC_ESTS        	(-52)	/* In Fortran, string too short */
#define NC_EMAXNAME    	(-53)	/* NC_MAX_NAME exceeded */
#define NC_EUNLIMIT    	(-54)	/* NC_UNLIMITED size already in use */
#define NC_ENORECVARS  	(-55)	/* nc_rec op when there are no record vars */
#define NC_ECHAR	(-56)	/* Attempt to convert between text & numbers */
#define NC_EEDGE	(-57)	/* Start+count exceeds dimension bound */
#define NC_ESTRIDE	(-58)	/* Illegal stride */
#define NC_EBADNAME	(-59)	/* Attribute or variable name
                                         contains illegal characters */
/* N.B. following must match value in ncx.h */
#define NC_ERANGE	(-60)	/* Math result not representable */
#define NC_ENOMEM	(-61)	/* Memory allocation (malloc) failure */
#define NC_EVARSIZE     (-62)   /* One or more variable sizes violate
				   format constraints */
#define NC_EDIMSIZE     (-63)   /* Invalid dimension size */
#define NC_ETRUNC       (-64)   /* File likely truncated or possibly corrupted */
/* Declaration modifiers for DLL support (MSC et al) */
#if defined(DLL_NETCDF) /* define when library is a DLL */
#   define MSC_EXTRA __declspec(dllimport)
#else
#define MSC_EXTRA
#endif	/* defined(DLL_NETCDF) */
# define EXTERNL extern MSC_EXTRA
/* When netCDF is built as a DLL, this will export ncerr and
 * ncopts. When it is used as a DLL, it will import them. */
MSC_EXTRA int ncerr;
MSC_EXTRA int ncopts;
EXTERNL const char *nc_inq_libvers(void);
EXTERNL const char *nc_strerror(int ncerr);
EXTERNL int nc__create(const char *path, int cmode, size_t initialsz,size_t *chunksizehintp, int *ncidp);
EXTERNL int nc_create(const char *path, int cmode, int *ncidp);
EXTERNL int nc__open(const char *path, int mode, size_t *chunksizehintp, int *ncidp);
EXTERNL int nc_open(const char *path, int mode, int *ncidp);
EXTERNL int nc_set_fill(int ncid, int fillmode, int *old_modep);
EXTERNL int nc_redef(int ncid);
EXTERNL int nc__enddef(int ncid, size_t h_minfree, size_t v_align,size_t v_minfree, size_t r_align);
EXTERNL int nc_enddef(int ncid);
EXTERNL int nc_sync(int ncid);
EXTERNL int nc_abort(int ncid);
EXTERNL int nc_close(int ncid);
EXTERNL int nc_inq(int ncid, int *ndimsp, int *nvarsp, int *nattsp, int *unlimdimidp);
EXTERNL int nc_inq_ndims(int ncid, int *ndimsp);
EXTERNL int nc_inq_nvars(int ncid, int *nvarsp);
EXTERNL int nc_inq_natts(int ncid, int *nattsp);
EXTERNL int nc_inq_unlimdim(int ncid, int *unlimdimidp);
EXTERNL int nc_set_default_format(int format, int *old_formatp);
EXTERNL int nc_inq_format(int ncid, int *formatp);
/* Begin _dim */
EXTERNL int nc_def_dim(int ncid, const char *name, size_t len, int *idp);
EXTERNL int nc_inq_dimid(int ncid, const char *name, int *idp);
EXTERNL int nc_inq_dim(int ncid, int dimid, char *name, size_t *lenp);
EXTERNL int nc_inq_dimname(int ncid, int dimid, char *name);
EXTERNL int nc_inq_dimlen(int ncid, int dimid, size_t *lenp);
EXTERNL int nc_rename_dim(int ncid, int dimid, const char *name);
/* End _dim */
/* Begin _att */
EXTERNL int nc_inq_att(int ncid, int varid, const char *name,	   nc_type *xtypep, size_t *lenp);
EXTERNL int nc_inq_attid(int ncid, int varid, const char *name, int *idp);
EXTERNL int nc_inq_atttype(int ncid, int varid, const char *name, nc_type *xtypep);
EXTERNL int nc_inq_attlen(int ncid, int varid, const char *name, size_t *lenp);
EXTERNL int nc_inq_attname(int ncid, int varid, int attnum, char *name);
EXTERNL int nc_copy_att(int ncid_in, int varid_in, const char *name, int ncid_out, int varid_out);
EXTERNL int nc_rename_att(int ncid, int varid, const char *name, const char *newname);
EXTERNL int nc_del_att(int ncid, int varid, const char *name);
/* End _att */
/* Begin {put,get}_att */
EXTERNL int nc_put_att(int ncid, int varid, const char *name, nc_type datatype,size_t len, const void *value);
EXTERNL int nc_get_att(int ncid, int varid, const char *name, void *value);
EXTERNL int nc_put_att_text(int ncid, int varid, const char *name,size_t len, const char *op);
EXTERNL int nc_get_att_text(int ncid, int varid, const char *name, char *ip);
EXTERNL int nc_put_att_uchar(int ncid, int varid, const char *name, nc_type xtype,size_t len, const unsigned char *op);
EXTERNL int nc_get_att_uchar(int ncid, int varid, const char *name, unsigned char *ip);
EXTERNL int nc_put_att_schar(int ncid, int varid, const char *name, nc_type xtype,size_t len, const signed char *op);
EXTERNL int nc_get_att_schar(int ncid, int varid, const char *name, signed char *ip);
EXTERNL int nc_put_att_short(int ncid, int varid, const char *name, nc_type xtype,size_t len, const short *op);
EXTERNL int nc_get_att_short(int ncid, int varid, const char *name, short *ip);
EXTERNL int nc_put_att_int(int ncid, int varid, const char *name, nc_type xtype,size_t len, const int *op);
EXTERNL int nc_get_att_int(int ncid, int varid, const char *name, int *ip);
EXTERNL int nc_put_att_long(int ncid, int varid, const char *name, nc_type xtype,	size_t len, const long *op);
EXTERNL int nc_get_att_long(int ncid, int varid, const char *name, long *ip);
EXTERNL int nc_put_att_float(int ncid, int varid, const char *name, nc_type xtype,size_t len, const float *op);
EXTERNL int nc_get_att_float(int ncid, int varid, const char *name, float *ip);
EXTERNL int nc_put_att_double(int ncid, int varid, const char *name, nc_type xtype,	size_t len, const double *op);
EXTERNL int nc_get_att_double(int ncid, int varid, const char *name, double *ip);
/* End {put,get}_att */
/* Begin _var */
EXTERNL int nc_def_var(int ncid, const char *name, nc_type xtype, int ndims,const int *dimidsp, int *varidp);
EXTERNL int nc_inq_var(int ncid, int varid, char *name, nc_type *xtypep,int *ndimsp, int *dimidsp, int *nattsp);
EXTERNL int nc_inq_varid(int ncid, const char *name, int *varidp);
EXTERNL int nc_inq_varname(int ncid, int varid, char *name);
EXTERNL int nc_inq_vartype(int ncid, int varid, nc_type *xtypep);
EXTERNL int nc_inq_varndims(int ncid, int varid, int *ndimsp);
EXTERNL int nc_inq_vardimid(int ncid, int varid, int *dimidsp);
EXTERNL int nc_inq_varnatts(int ncid, int varid, int *nattsp);
EXTERNL int nc_rename_var(int ncid, int varid, const char *name);
EXTERNL int nc_copy_var(int ncid_in, int varid, int ncid_out);
#ifndef ncvarcpy
/* support the old name for now */
#define ncvarcpy(ncid_in, varid, ncid_out) ncvarcopy((ncid_in), (varid), (ncid_out))
#endif
/* End _var */
/* Begin {put,get}_var1 */
EXTERNL int nc_put_var1(int ncid, int varid, const size_t *indexp, const void *value);
EXTERNL int nc_get_var1(int ncid, int varid, const size_t *indexp, void *value);
EXTERNL int nc_put_var1_text(int ncid, int varid, const size_t *indexp, const char *op);
EXTERNL int nc_get_var1_text(int ncid, int varid, const size_t *indexp, char *ip);
EXTERNL int nc_put_var1_uchar(int ncid, int varid, const size_t *indexp,const unsigned char *op);
EXTERNL int nc_get_var1_uchar(int ncid, int varid, const size_t *indexp,unsigned char *ip);
EXTERNL int nc_put_var1_schar(int ncid, int varid, const size_t *indexp,const signed char *op);
EXTERNL int nc_get_var1_schar(int ncid, int varid, const size_t *indexp,signed char *ip);
EXTERNL int nc_put_var1_short(int ncid, int varid, const size_t *indexp,const short *op);
EXTERNL int nc_get_var1_short(int ncid, int varid, const size_t *indexp,short *ip);
EXTERNL int nc_put_var1_int(int ncid, int varid, const size_t *indexp, const int *op);
EXTERNL int nc_get_var1_int(int ncid, int varid, const size_t *indexp, int *ip);
EXTERNL int nc_put_var1_long(int ncid, int varid, const size_t *indexp, const long *op);
EXTERNL int nc_get_var1_long(int ncid, int varid, const size_t *indexp, long *ip);
EXTERNL int nc_put_var1_float(int ncid, int varid, const size_t *indexp, const float *op);
EXTERNL int nc_get_var1_float(int ncid, int varid, const size_t *indexp, float *ip);
EXTERNL int nc_put_var1_double(int ncid, int varid, const size_t *indexp, const double *op);
EXTERNL int nc_get_var1_double(int ncid, int varid, const size_t *indexp, double *ip);
/* End {put,get}_var1 */
/* Begin {put,get}_vara */
EXTERNL int nc_put_vara(int ncid, int varid,  const size_t *start, const size_t *count, const void *value);
EXTERNL int nc_get_vara(int ncid, int varid,  const size_t *start, const size_t *count, void *value);
EXTERNL int nc_put_vara_text(int ncid, int varid, const size_t *startp, const size_t *countp, const char *op);
EXTERNL int nc_get_vara_text(int ncid, int varid, const size_t *startp, const size_t *countp, char *ip);
EXTERNL int nc_put_vara_uchar(int ncid, int varid, const size_t *startp, const size_t *countp, const unsigned char *op);
EXTERNL int nc_get_vara_uchar(int ncid, int varid, const size_t *startp, const size_t *countp, unsigned char *ip);
EXTERNL int nc_put_vara_schar(int ncid, int varid, const size_t *startp, const size_t *countp, const signed char *op);
EXTERNL int nc_get_vara_schar(int ncid, int varid, const size_t *startp, const size_t *countp, signed char *ip);
EXTERNL int nc_put_vara_short(int ncid, int varid, const size_t *startp, const size_t *countp, const short *op);
EXTERNL int nc_get_vara_short(int ncid, int varid, const size_t *startp, const size_t *countp, short *ip);
EXTERNL int nc_put_vara_int(int ncid, int varid, const size_t *startp, const size_t *countp, const int *op);
EXTERNL int nc_get_vara_int(int ncid, int varid, const size_t *startp, const size_t *countp, int *ip);
EXTERNL int nc_put_vara_long(int ncid, int varid, const size_t *startp, const size_t *countp, const long *op);
EXTERNL int nc_get_vara_long(int ncid, int varid, const size_t *startp, const size_t *countp, long *ip);
EXTERNL int nc_put_vara_float(int ncid, int varid, const size_t *startp, const size_t *countp, const float *op);
EXTERNL int nc_get_vara_float(int ncid, int varid, const size_t *startp, const size_t *countp, float *ip);
EXTERNL int nc_put_vara_double(int ncid, int varid, const size_t *startp, const size_t *countp, const double *op);
EXTERNL int nc_get_vara_double(int ncid, int varid, const size_t *startp, const size_t *countp, double *ip);
/* End {put,get}_vara */
/* Begin {put,get}_vars */
EXTERNL int nc_put_vars(int ncid, int varid,  const size_t *start, const size_t *count, const ptrdiff_t *stride,  const void * value);
EXTERNL int nc_get_vars(int ncid, int varid,  const size_t *start, const size_t *count, const ptrdiff_t *stride,  void * value);
EXTERNL int nc_put_vars_text(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const char *op);
EXTERNL int nc_get_vars_text(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, char *ip);
EXTERNL int nc_put_vars_uchar(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const unsigned char *op);
EXTERNL int nc_get_vars_uchar(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, unsigned char *ip);
EXTERNL int nc_put_vars_schar(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const signed char *op);
EXTERNL int nc_get_vars_schar(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, signed char *ip);
EXTERNL int nc_put_vars_short(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const short *op);
EXTERNL int nc_get_vars_short(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, short *ip);
EXTERNL int nc_put_vars_int(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const int *op);
EXTERNL int nc_get_vars_int(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, int *ip);
EXTERNL int nc_put_vars_long(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const long *op);
EXTERNL int nc_get_vars_long(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, long *ip);
EXTERNL int nc_put_vars_float(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const float *op);
EXTERNL int nc_get_vars_float(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, float *ip);
EXTERNL int nc_put_vars_double(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const double *op);
EXTERNL int nc_get_vars_double(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, double *ip);
/* End {put,get}_vars */
/* Begin {put,get}_varm */
EXTERNL int nc_put_varm(int ncid, int varid, const size_t *start, const size_t *count, const ptrdiff_t *stride, const ptrdiff_t *imapp,     const void *value);
EXTERNL int nc_get_varm(int ncid, int varid, const size_t *start, const size_t *count, const ptrdiff_t *stride, const ptrdiff_t *imapp, void *value);
EXTERNL int nc_put_varm_text(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, const char *op);
EXTERNL int nc_get_varm_text(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, char *ip);
EXTERNL int nc_put_varm_uchar(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, const unsigned char *op);
EXTERNL int nc_get_varm_uchar(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, unsigned char *ip);
EXTERNL int nc_put_varm_schar(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, const signed char *op);
EXTERNL int nc_get_varm_schar(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, signed char *ip);
EXTERNL int nc_put_varm_short(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, const short *op);
EXTERNL int nc_get_varm_short(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, short *ip);
EXTERNL int nc_put_varm_int(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, const int *op);
EXTERNL int nc_get_varm_int(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, int *ip);
EXTERNL int nc_put_varm_long(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, const long *op);
EXTERNL int nc_get_varm_long(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, long *ip);
EXTERNL int nc_put_varm_float(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, const float *op);
EXTERNL int nc_get_varm_float(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, float *ip);
EXTERNL int nc_put_varm_double(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t *imapp, const double *op);
EXTERNL int nc_get_varm_double(int ncid, int varid, const size_t *startp, const size_t *countp, const ptrdiff_t *stridep, const ptrdiff_t * imapp, double *ip);
/* End {put,get}_varm */
/* Begin {put,get}_var */
EXTERNL int nc_put_var_text(int ncid, int varid, const char *op);
EXTERNL int nc_get_var_text(int ncid, int varid, char *ip);
EXTERNL int nc_put_var_uchar(int ncid, int varid, const unsigned char *op);
EXTERNL int nc_get_var_uchar(int ncid, int varid, unsigned char *ip);
EXTERNL int nc_put_var_schar(int ncid, int varid, const signed char *op);
EXTERNL int nc_get_var_schar(int ncid, int varid, signed char *ip);
EXTERNL int nc_put_var_short(int ncid, int varid, const short *op);
EXTERNL int nc_get_var_short(int ncid, int varid, short *ip);
EXTERNL int nc_put_var_int(int ncid, int varid, const int *op);
EXTERNL int nc_get_var_int(int ncid, int varid, int *ip);
EXTERNL int nc_put_var_long(int ncid, int varid, const long *op);
EXTERNL int nc_get_var_long(int ncid, int varid, long *ip);
EXTERNL int nc_put_var_float(int ncid, int varid, const float *op);
EXTERNL int nc_get_var_float(int ncid, int varid, float *ip);
EXTERNL int nc_put_var_double(int ncid, int varid, const double *op);
EXTERNL int nc_get_var_double(int ncid, int varid, double *ip);
#ifdef LOGGING
#ifdef DEBUG
EXTERNL void nc_exit(void);
#endif
/* This is only defined for netcdf-4 apparently */
#ifdef USE_NETCDF4
EXTERNL void nc_set_log_level(int new_level);
#endif
#define NC_TURN_OFF_LOGGING (-1)
#else /* not LOGGING */
#define nc_set_log_level(e)
#endif
/* End {put,get}_var */
EXTERNL int nc__create_mp(const char *path, int cmode, size_t initialsz, int basepe,  size_t *chunksizehintp, int *ncidp);
EXTERNL int nc__open_mp(const char *path, int mode, int basepe, size_t *chunksizehintp, int *ncidp);
EXTERNL int nc_delete(const char * path);
EXTERNL int nc_delete_mp(const char * path, int basepe);
EXTERNL int nc_set_base_pe(int ncid, int pe);
EXTERNL int nc_inq_base_pe(int ncid, int *pe);
