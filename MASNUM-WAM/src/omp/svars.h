#ifndef SVARS_H_INCLUDED
#define SVARS_H_INCLUDED
#include "wavedefc.h"
#include "std.h"
#ifdef ARM_VEC
# include "zarm.h"
#endif
#ifdef USEARM_SME
#define ARM_STREAMING __arm_streaming  
#else
#define ARM_STREAMING
#endif
#ifndef DOUBLE 
# define DOUBLE double
#endif
#define DBGLN printf("%s %d",__FILE__,__LINE__)
#define DBGL printf("%s %d\n",__FILE__,__LINE__)

typedef union _prgstate{
  double vlt;
  struct {
    long   state;
    long   ntimestep;
    long   slave_state;
    long   tstate;
  };
}prgstate;

#define ZBYTE signed char
#define  IMDPS   4
typedef struct _ImplschPar{
  double  wkal[4*12]  ;      //0x180
  long    kjps[CJNTHET*24];  //0x900 (  2,12,jn)
  ZBYTE    kjs[CJNTHET*48];  //0x900 (2,2,12,jn)
  double cwks17[CKL];        //0x100
  double grolim[CKL];        //0x100
  double wk[CKL],wkh[CKL];   //0x200
  double dwk[CKL],wkdk[CKL],wkibdk[CKL];     //0x300
#if SIGMALVL==2
  double expdep[CKL*(MBVDEP)];
#endif
  double  cosths[CJNTHET],sinths[CJNTHET]  ;     //0xc0
  double  lmdpd,lmdmd,lmdpd20,lmdmd20,Dm2,D1;
  float   beta10,beta11;        //0x88 beta11
  float   ads,brkd1,brkd2,deltts;
  float spdeltts;
  int   ntsplit;

  double  sds,awk,enh;
  double  zpi;


  float   D1m30,F44_24778,F16_0,F0_467;//mean3
  float   F18_20832,F0_001930368,F0_5,F1_0;//mean2
  float   F0_75,F5_5,F0_833,Fm1_25;//enh
  float   F28_,F0_80,F0_065,F0_001;
  float   tztp,tzzpi;

  ImplschP*pip;
  int iacb,iace;
  int mpi_id,tid,dbglvl,dmalog;     //0x138
  int nbvlvl,nwps,nwpc,nwpa;     //0x148
}ImplschPar;

//void checksp(gdata *pgpp,int ll);
#define CHECKESP  //checksp(pgpp,__LINE__);
#define CHECKMEM  //CheckMem(pgpp,__LINE__);

typedef struct _Boundary{
  double cosths[CJNTHET],sinths[CJNTHET];
  double rwk[CKL];
  int nwpc;
  double rg,g,gama,zpi,windfield,xj0,xj_,arlfa_,wsj_;
  double rsigma1,rsigma2,wv0;
  double F_1,F_m1d25,F_m0d5,F_m0d4,F_m0d33;
  double* pvws;
  char*nsp;
}Boundary;
typedef struct _Stencil_Smooth{
  double em[mkj*9];
  double stc[9], stcv[8];
  int ind[9],old_neb[9], all[9];
  short clcy; 
  char inited_neb, inited_upem;
} SSmooth;
typedef struct TData{ 
  double *wav_ec, *wav_et;
  int iacb1,iace1,iacb2,iace2;
  int volatile nstep,nostep;
  SSmooth smt;
}TData;
#ifdef MTHREAD
#include "mthread.h"
#ifdef MMTHREAD
#include "gpart.h"
#endif
#define __THREAD __thread 
#else
# define __THREAD 
#endif
typedef struct GData{ 
  int mpi_id,mpi_npe,mpi_comm_wav;
  int volatile hist_eot,stop_now,rest_eot;
  int nstep,nostep;
  int nwps,nwpc,nwpa;
  int LXB,LYB,LXN,LYN;
  int mTimeStep;
  windvs constwind;
  int   uconstwind;
#ifdef MMTHREAD
  int NPART,block_id,nwpo,nwpt;//GD
  int *l2gind; //GD
  sendcomm s_segs;//GD
  Rect recti,recto; //GD
#endif
  char*fbuff;//HBM
  //struct g2lInd *g2l;
  int*gg2l;
  double *wav_ec, *wav_et;
  char *nsp;    //g_initSetspec
  Boundary bdy;
  float *dep;
  spc_interg *psis; // g_initPropgats 
  geo_interg *pgis; // g_initPropgats 
#ifdef ARM_VEC
  struct vspc_interg *psvis;// g_initPropgats 
  struct vgeo_interg *pgvis;// g_initPropgats 
#endif
  int *ipos12; 
  propinf*ipos8;

  int *ieind;
  struct xnlinf*xi;
  struct Iepos*iepos;
  windvs *pwvs;
  windvs *pwvsg;
  double*tee;
  /*inited at g_initImplsch OK*/
  ImplschPar pisp;
  /*inited at g_initOutput_cal OK */
  int ndep,npdep;
  double *pvkd,*pebdep;
  float *ape , *tpf , *aet , *h1_3, *bv;
}GData;
__BEGIN_DECLS
extern __THREAD struct TData *td;
extern __THREAD struct GData *gd;
extern GData pd;
#ifdef MMTHREAD
void c_fillsendcomm_();
// in propagat.cpp
extern void f_initGPar();
extern void g_initPropgats();
extern void g_initImplsch();
extern void g_initSetspec();
// in svars.c
void g_setoutput_();
#endif
extern ImplschPar pisp;
// in propagat.cpp
extern void InitCPropgats();
extern void InitCImplsch();
void InitGPar();
void c_initimplsch_(ImplschPar *pisp_,ImplschP*pip);
void c_setpropinterg_(int *iac_,int*j_,int*k_,int*iquad_,int*kpos,float*pab,float*pabp);
void c_setwind_();
void c_setwind();
void setflag_(int*nostep,int *hist_eot,int *stop_now,int *rest_eot);
extern void setspec(int mode,double *e,int nwpb,int nwpe,int n) ;
extern void implschs_(double *eec,double *eet,int*iacb,int*iace) ;
extern void implschs (double *eec,double *eet,int iacb,int iace) ;
extern void setspec_(int*mode,double *eet,int *iacb,int *iace,int *type) ;
extern void accumea_(int*ind,double*ee,int*nwpb,int *nwpe);
extern void checkoutput_(int*igrp,int*ind,int*RFB,double*ee,int*iacb,int*iace,int*iacb2,int*iace2,int*nout);
extern void t_endoutput_cal_(int*iacb1,int*iace1);
extern void g_initOutput_cal();
extern void propagats_spec(double *ed,double *es,int iacb,int iace) ;
extern void propagats_geo(double *ed,double *es,int iacb,int iace) ;
extern void propagats(double *wav_e,double *wav_ee,int iacb,int iace);
extern void propagats_spec_(double *ed,double *es,int *iacb,int *iace) ;
extern void propagats_geo_(double *ed,double *es,int *iacb,int *iace) ;
extern void propagats_(double *wav_e,double *wav_ee,int *iacb,int *iace);
void addea(double*ea,double*ee,int iacb,int iace);
void addea_(double*ea,double*ee,int*iacb,int*iace);
void averageea(double*ea,int nea,int iacb,int iace);
void averageea_(double*ea,int*nea,int*iacb,int*iace);
void zeroea_(double*ea,int*iacb,int*iace);
void zeroea(double*ea,int iacb,int iace);
void init_xnlg_();
void init_xnl_();
void init_dia_();
int xnlcal(double*ee,int ks,double*xnl,double*diag);
void xnlcal_(double*nspec,int*ks,double*xnl,double*diag, int*ierror);
struct wws_type;
void xnlsetgrd_(int*naq,int*nkq,int*mkq,int*maq,int*klocus,int*qf_dn,int*qf_kn,
             struct wws_type*quads,int*quadnl,
             double*q_k2,double*q_ka,double*q_sig,double*q_sigr,
             double*q_kpow,double*q_tail,double*q_lamd);
void xnlsetgrd(int naq,int nkq,int mkq,int maq,int klocus,int qf_dn,int qf_kn,
            struct wws_type*quads,int*quadnl,
            double*q_k2,double*q_ka,double*q_sig,double*q_sigr,
            double*q_kpow,double*q_tail,double*q_lamd);
extern void diffract_smooth_(double *wav_e,double *wav_ee,int *iacb,int *iace);
extern void diffract_smooth(double *et, double *ee, int iacb, int iace) ARM_STREAMING;
extern void init_SSmooth();

void initsetspec_  (int*iacb,int*iace);
void initimplschs_ (int*iacb,int*iace);
void initpropagat_ (int*iacb,int*iace);
void initoutputcal_(int*iacb,int*iace,float *dep);
__END_DECLS
#endif //SVARS_H_INCLUDED
