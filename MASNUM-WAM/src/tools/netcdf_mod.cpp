#include "netcdf_mod.h"
  char* version = "netcdf_mod 2.6,by zhaowei, 20014-3-1 .";
//-------------------------------------------------------------------------------
  contains
//-------------------------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//-------------------------------------------------------------------------------
// --- To handle errors for NCFile input or output.
//* checknc_err
#define Mchecknc_err(stat) checknc_err(stat,__LINE__)
  void  checknc_err(int status,int cline){
	  if (status != NC_NOERR) {
	    printf("%d %s\n",cline,NC_strerror(status))
	    exit(0);
	  }
  } 
//-------------------------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//-------------------------------------------------------------------------------
// --- To open a NcFile. (By create, if new; By redefine, if old.)
//* open_nc_def
  int open_nc(char* FileName, char*action){
  int  NcID;
  int alive;
  int  status;
  switch(action[0]|0x40){
  case 'c':
    status = nc_create(FileName, NC_CLOBBER+NC_share, &ncID)
    Mchecknc_err(status);
  case 'd':
  	alive=access(FileName,0);
  	if(!alive){
  		printf("file %s for define is not exist\n",FileName);
  		exit(0);
  	}
    status = nc_open(FileName, NC_WRITE+NC_share, *NcID);
    Mchecknc_err(status);
    status = nc_redef(NcID);
    Mchecknc_err(status);
  case'w':
  	alive=access(FileName,0);
  	if(!alive){
  		printf("file %s for write is not exist\n",FileName);
  		exit(0);
  	}
    status = nc_open(FileName, NC_WRITE+NC_share, &NcID);
    Mchecknc_err(status)
  case'r':
  	alive=access(FileName,0);
  	if(!alive){
  		printf("file %s for read is not exist\n",FileName);
  		exit(0);
  	}
    status = nc_open((FileName), NC_NOWRITE+NC_share, NcID);
    Mchecknc_err(status);
  }
  return NcID;
  } //  open_nc
  void  close_nc(int &NcID){
	  Mchecknc_err(NC_CLOSE(*NcID));
	  *NcID=-1;
  } 
  void  end_define(int NcID){
  char lib_version[80];
  char *ts;
  int  status;
	__time64_t ltime;
  _time64( &ltime );  
  ts=_ctime64( &ltime )
  status=set_attribute(NcID,NC_GLOBAL, "CreatedTime", ts);
  sprintf(lib_version,"netcdf %s",nc_inq_libvers());  
  status=set_attribute(NcID,NC_GLOBAL, "Sorftware1", lib_version);
  status=set_attribute(NcID,NC_GLOBAL, "Sorftware2", version);
  status = nc_enddef(NcID);
  Mchecknc_err(status);
  } //  end_define
  void  dimension_define(int NcID, char* DimName, int DimLen,int *DimID, char*DimVarName, int DimVarType, int *DimVarID){
	  int  DimID1=-1, DimVarID1=-1;
	  int  status;
	  status = nc_def_dim(NcID, DimName, DimLen, &DimID1);
	  Mchecknc_err(status);
	  if(DimVarName&&DimVarType){
		  status = nc_def_var(ncid, DimVarName, DimVarType, 1, DimID1, &DimVarID1);
		  Mchecknc_err(status);
		}
	  if(DimID)*DimID = DimID1;
	  if(DimVarID)*DimVarID = DimVarID1;
  } 
  int  get_dimension_len(int NcID, char*DimName){
	  int   DimID,status,dimlen;
	  status = nc_inq_dimid(NcID, DimName, &DimID)
	  Mchecknc_err(status);
	  status = nc_inq_dimlen(NcID, DimID, &dimlen)
	  Mchecknc_err(status);
	  return dimlen;
	}
  void  variable_define(int NcID, char*VarName, int VarType, int VarRank,int *VarDims, int*VarID){
	  int   VarRank, VarID1,status;
	  status=nc_inq_varid(NcID, (VarName), &VarID1)
	  if(status != NC_NOERR){
	    status = nc_def_var(ncid, VarName, VarType, VarRank, VarDims, &VarID1);
	    Mchecknc_err(status);
	  }
	  if(VarID)*VarID = VarID1;
  } 
  void  variable_define(int NcID, char*VarName, int VarType, int   VarRank,char**VarDimsName, int*VarID){
	  int   VarRank, VarID1, i
	  int VarDims[NC_MAX_VAR_DIMS];
	  int  status;
	  status=nc_inq_varid(NcID, (VarName), &VarID1)
	  if(status != NC_NOERR){	    
	    for( i = 0;i<VarRank;i++){
	      status = nc_inq_dimid(NcID, VarDimsName(i), &VarDims(i));
	      Mchecknc_err(status);
	    }
	    status = nc_def_var(ncid, (VarName), VarType, VarRank, VarDims, &VarID1);
	    Mchecknc_err(status);
	  }
	  if(VarID)*VarID = VarID1;
  } 
  int  set_attribute_i_TEXT(int NcID,int VarID, char*AttName, char*attribute){
  	if(VarID<0)VarID=NC_GLOBAL;
  	return nc_put_att_text(NcID, VarID, AttName,strlen(attribute), attribute);
	}
  int  get_attribute_i_TEXT(int NcID,int VarID, char*AttName, char*attribute){
  	if(VarID<0)VarID=NC_GLOBAL;
  	return nc_gut_att_text(NcID, VarID, AttName,attribute);
	}
#define TMPL_set_attribute_i(NA,NAT,VT)                                \
  int  set_attribute(int NcID,int VarID, char*AttName, VT *attribute){ \
	  if(VarID<0)VarID=NC_GLOBAL;                                        \
	  return nc_put_att_##NA(NcID, VarID, AttName,NAT,1, attribute);     \
	}
#define TMPL_get_attribute_i(NA,NAT,VT)                                \
  int get_attribute(int NcID,int VarID, char*AttName, VT*attribute)    \
    if(VarID<0)VarID=NC_GLOBAL;                                        \
    return NC_GET_ATT_##NA(NcID, VarID, AttName, attribute);           \
	}
  TMPL_set_attribute_i(uchar ,INT1 ,int1 )
  TMPL_set_attribute_i(short ,INT2 ,int2 )
  TMPL_set_attribute_i(int   ,INT4 ,int4 )
  TMPL_set_attribute_i(float ,REAL4,real4)
  TMPL_set_attribute_i(double,REAL8,real8)
  TMPL_get_attribute_i(uchar ,INT1 ,int1 )
  TMPL_get_attribute_i(short ,INT2 ,int2 )
  TMPL_get_attribute_i(int   ,INT4 ,int4 )
  TMPL_get_attribute_i(float ,REAL4,real4)
  TMPL_get_attribute_i(double,REAL8,real8)
#define TMPL_set_attribute_n(VT)                                          \
  int  set_attribute(int NcID, char*AttName, VT*attribute,char*VarName){  \
    int   VarID=GetVarid(NcID,VarName);                                   \
    return set_attribute(NcID, VarID, AttName, attribute);                \
	}
#define TMPL_get_attribute_n(NA,VT)                                       \
  int  get_attribute(int NcID,char*AttName, VT*attribute,char*VarName){   \
  int   VarID=GetVarid(NcID,VarName)；                                    \
  return get_attribute(NcID, VarID, AttName, attribute);                  \
	}
  TMPL_set_attribute_i(char     )
  TMPL_set_attribute_i(int1     )
  TMPL_set_attribute_i(int2     )
  TMPL_set_attribute_i(int4     )
  TMPL_set_attribute_i(real4    )
  TMPL_set_attribute_i(real8    )
  TMPL_get_attribute_i(char     )
  TMPL_get_attribute_i(int1     )
  TMPL_get_attribute_i(int2     )
  TMPL_get_attribute_i(int4     )
  TMPL_get_attribute_i(real4    )
  TMPL_get_attribute_i(real8    )

#define TMPL_readwritenc_i(VT)                                      \
  void  writenc(int ncid, int VarID, VT*Var, int RecNum, int*locs){ \
    int  status,DimNum = 1,starts[8]={1}, counts[8]={1};            \
    SetVARSC(starts,counts,dimnum,shape(Var),var ,RecNum,locs);     \
    status = nc_put_vara(ncID, varID, starts, counts, Var);         \
    Mchecknc_err(status);                                           \
  }                                                                 \
  void  readnc(int ncid, int VarID, VT*Var, int RecNum, int*locs){  \
    int  status,DimNum = 1,starts[8]={1}, counts[8]={1};            \
    SetVARSC(starts,counts,dimnum,shape(Var),var ,RecNum,locs);     \
    status = nc_get_vara(ncID, varID, starts, counts, Var);         \
	  Mchecknc_err(status);                                           \
  } 
	TMPL_readwritenc_i(int1 )
	TMPL_readwritenc_i(int2 )
	TMPL_readwritenc_i(int4 )
	TMPL_readwritenc_i(real4)
	TMPL_readwritenc_i(real8)
	TMPL_readwritenc_i(text )

#define TMPL_readwritenc_n(ND,NA,VT,VSHAPE)                            \
  void  writenc(int ncid, char*VarName, VT*Var, int RecNum, int*locs){ \
	  writenc(ncid, GetVarid(ncID,VarName), Var, RecNum,locs);           \
  }                                                                    \
  void  readnc(int ncid, char*VarName, VT*Var, int RecNum, int*locs){  \
	  readnc(ncid, GetVarid(ncID,VarName), Var, RecNum,locs);            \
  } 
	TMPL_readwritenc_i(int1 )
	TMPL_readwritenc_i(int2 )
	TMPL_readwritenc_i(int4 )
	TMPL_readwritenc_i(real4)
	TMPL_readwritenc_i(real8)
	TMPL_readwritenc_i(text )

//-------------------------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//-------------------------------------------------------------------------------
int GetVarid(int ncID,char*vname){
	int  vid=NC_GLOBAL;
	if(vname)  Mchecknc_err(nc_inq_varid(ncID,vname,&vID));
	return vid;
}
int GetVarType_n(int ncID,char*VName){
	int  vid,vtype;
  Mchecknc_err(nc_inq_varid(ncID,vname,&vID));
  Mchecknc_err(NC_INQ_VARTYPE(ncID,vID,vtype));
  return vtype;
}
int  function GetVarType(int ncID,int Vid){
	int  vtype;
  Mchecknc_err(nc_inq_vartype(ncID,vID,vtype))
  return vtype;
}
int  nfDefVar(int ncid,char*name,int logtype,real8 Sc,real8 Off,int ndim,int*vdim){
  if(logtype==0){
    nfDefVar=-1;
  }else{
    Mchecknc_err(nc_def_var(ncID,name,NC_real4,ndim,vdim,&varid));
    Mchecknc_err(set_attribute(ncID,varid,"missing_value",1.71e38));
  }
  return varid;
}
	void  SetVARSC_TEXT(starts,counts,DimNum,vshape,var,RecNum,locs)
    int ,intent(inout)DimNum
    int ,intent(out) starts(:),counts(:)
    int ,intent(in) vshape(:)
    char(*),intent(in)var
	  int , optional RecNum,locs(:)
    DimNum=size(vshape);starts=1;counts=1;
    if(DimNum>0)counts(1:DimNum)=vshape
    counts(1) = len(var);counts(2:size(vshape)+1)=vshape
    if((locs)){
    	DimNum=size(locs);
    	starts(2:DimNum+1) = locs(1:DimNum)
    }
    if((RecNum))starts(DimNum+2) = RecNum
	} //  SetVARSC_TEXT
#define TMPL_SetVARSC(NA,VT)                                                   \
  void  SetVARSC_##NA(starts,counts,DimNum,vshape,var,RecNum,locs)       ;\
    int ,intent(inout)DimNum                                             ;\
    int ,intent(out) starts(:),counts(:)                                 ;\
    int ,intent(in) vshape(:)                                            ;\
    VT,intent(in)var                                                        ;\
	  int , optional RecNum,locs(:)                           ;\
    DimNum=size(vshape);starts=1;counts=1;                                    ;\
    if(DimNum>0)counts(1:DimNum)=vshape                          						  ;\
    if((locs)){                                                     ;\
    	DimNum=size(locs);                                                      ;\
    	starts(1:DimNum) = locs(1:DimNum)                                       ;\
    }                                                                     ;\
    if((RecNum))starts(DimNum+1) = RecNum                              ;\
  } //  SetVARSC_##NA
  TMPL_SetVARSC(INT1  ,int1  )
  TMPL_SetVARSC(INT2  ,int2  )
  TMPL_SetVARSC(INT   ,int      )
  TMPL_SetVARSC(REAL  ,real4      )
  TMPL_SetVARSC(real8,real8      )
