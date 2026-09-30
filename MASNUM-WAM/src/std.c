#define _GNU_SOURCE
#include "std.h"
#include "ctools.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <unistd.h>
prggpar gppar;
typedef struct _FC_PARA{
  int _kl,_jnthet,_kld,_NBVDEP,_SIGMALVL,_LOGSCURR;
  int gixl,giyl;
  int nwps,nwpc,nwpa;
  int LXB,LYB,LXN,LYN;
  int mpi_id,mpi_comm_wav,mpi_npe;
  float constwindx,constwindy;
  int mTimeStep;
}FC_PARA;
struct Iepos{
  int ix,iy;
};
void prtieind(struct Iepos*iepos,int*ieind,char*nsp);
// Set gppar for init,before memory allocated
// used by c_allocate_vee
void c_setgpar_(FC_PARA *fcp,struct Iepos*iepos,int*ieind,char*nsp,float*dep){
  int acwv;
  //printf("AAA %d %d\n",fcp->nwpc,fcp->nwpa);
  //pisp.mpi_id=fcp->mpi_id;
  gppar.mpi_comm_wav=fcp->mpi_comm_wav;
  gppar.nbvdep=fcp->_NBVDEP;
#define VD(v) gppar.v=fcp->v;
  VD(mpi_id ); VD(mpi_npe);
  VD(nwps   ); VD(nwpc   ); VD(nwpa   );
  VD(LXB    ); VD(LYB    );
  VD(LXN    ); VD(LYN    );
  VD(mTimeStep);
#undef VD
  gppar.nwpw  =0;
  gppar.nsp   =nsp;
  gppar.ieind =ieind;
  gppar.iepos =iepos;
  gppar.constwind.wx=fcp->constwindx;
  gppar.constwind.wy=fcp->constwindy;

  gppar.constwind.wv=gppar.constwind.wx*gppar.constwind.wx+gppar.constwind.wy*gppar.constwind.wy;
  gppar.constwind.wi=0;
  gppar.dep   =dep;
  acwv=gppar.constwind.wv*100;
  if(acwv>=10){
    gppar.uconstwind=2;
  }else if(acwv>=1){
    gppar.uconstwind=1;
  }else {
    gppar.uconstwind=0;
  }
  //prtieind(iepos,ieind,nsp);
}
int iwalltime();
void startinit_() {
	char hostname[256];
	gethostname(hostname, sizeof(hostname));
  char *p=getenv("OMPI_COMM_WORLD_RANK");
	printf("%s %s %8d BS:\n", hostname,p,iwalltime());
  fflush(stdout);
}
//#define USE_C_ALLOC
//#ifndef MMTHREAD
#ifdef USE_C_ALLOC //must define USE_C_ALLOC in platform_init_mod.F90
void f_setevar_  (double*,double*);
void f_setimpvar_(ImplschP*ip,double*pvws,windvs *pwvs);
void f_setprop_  (spc_interg*psis,geo_interg*pgis,double*pvkdp,propinf*ipos8,int *ipos12);
void f_setoutput_(double*cpvkdo,double*cpebdep);
static char *fbuff=NULL;
void c_allocate_vee_(int *res){
  char*buf;
  size_t size=0;
  //printf("c_alloc %d %d\n",gppar.nwpc,gppar.nwpa);
#define IMDPSP  2
#define IMDPSO  4
#define BVIMDPS 1
#define VAVARS \
  VVAR(pebdep ,double  ,(size_t)(gppar.nwpc+1)*kld*(BVIMDPS+gppar.nbvdep));\
  VVAR(wav_ec ,double  ,(size_t)(gppar.nwpa+1)*mkj);\
  VVAR(wav_et ,double  ,(size_t)(gppar.nwpa+1)*mkj);\
  VVAR(ip     ,ImplschP,(size_t)(gppar.nwpc+1)    );\
  VVAR(pvws   ,double  ,(size_t)(gppar.nwpc+1)*kl );\
  VVAR(pvkdp  ,double  ,(size_t)(gppar.nwpc+1)*kl *IMDPSP);\
  VVAR(ipos12 ,int     ,(size_t)(gppar.nwpc+1)*12 );\
  VVAR(ipos8  ,propinf ,(size_t)(gppar.nwpc+1)    );\
  VVAR(psis   ,spc_interg,(size_t)(gppar.nwpc+1)*mkj);\
  VVAR(pgis   ,geo_interg,(size_t)(gppar.nwpc+1)*mkj);\
  VVAR(pwvs   ,windvs  ,(size_t)(gppar.nwpc+1)    );\
  VVAR(pvkdo  ,double  ,(size_t)(gppar.nwpc+1)*kld*IMDPSO);
#define  VVAR(vn,T,N) size+=sizeof(T)*(size_t)(N)+256
  VAVARS
#undef VVAR
#ifdef MMTHREAD
  buf=fbuff=PNEWN(char,size);
#else
  buf=fbuff=HNEWN(char,size);
#endif
  if(buf==NULL){
    printf(" alloc Fortran Mem error %ldMB\n",size>>20);
    *res=-1;
    return ;
  }
  for(size_t i=0;i<size;i+=4096) *(long*)(buf+i)=0;
  //memset(buf,0,size);
  buf=(char*)(( ((long)buf)+255)&(-256) );
#define  VVAR(vn,T,N) gppar.vn=(T*)buf;buf=(char*)( (((long)buf)+sizeof(T)*(size_t)(N)+255) & (-256) )
  VAVARS
#undef VVAR
#undef VAVARS
  f_setevar_(gppar.wav_ec,gppar.wav_et);
  f_setimpvar_(gppar.ip ,gppar.pvws ,gppar.pwvs);
  f_setprop_(gppar.psis,gppar.pgis,gppar.pvkdp  ,gppar.ipos8,gppar.ipos12);
  f_setoutput_(gppar.pvkdo  ,gppar.pebdep );
  *res= 0;
}
#endif
//set vars pointer in gppar
void c_setpointers_(double*eec,double*eet,ImplschP*ip,spc_interg*ps_vs,geo_interg*pg_vs,int*ipos12,propinf*ipos8,windvs*pwvs,double*pvkdo,double*pebdep){
  gppar.wav_ec =eec   ;
  gppar.wav_et =eet   ;
  gppar.ipos12 =ipos12;
  gppar.ipos8  =ipos8 ;
  gppar.pwvs   =pwvs  ;
  gppar.ip     =ip    ;
  gppar.psis   =ps_vs ;
  gppar.pgis   =pg_vs ;
  gppar.pvkdo  =pvkdo ;
  gppar.pebdep =pebdep;
}
void c_init_distinct_(){
  register_segfault_handler_();
#ifdef C_CALCULATE
  InitGPar();
#endif
#ifdef MMTHREAD
  InitThreadsGPar();
#endif
#ifndef MMTHREAD
#ifdef C_CALCULATE
  //NO MMTHREAD can init here
  InitCPropgats();
#endif
#endif
}
//#endif
int GetVInt(int volatile *volatile p){
	return *p;
}
int getvint_(int volatile *volatile p){
	return *p;
}
