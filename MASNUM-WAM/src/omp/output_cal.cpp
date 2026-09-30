#include "svars.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#ifndef CONSTKGRID
#error Must use CONSTKGRID
#endif
//!-----------------------------------------------
//!计算波高周期波向等
//!-----------------------------------------------
#define  IWFDK   0
#define  IIWF2DK 1
#define  IWS2DK  2
#define  IWF     3
#define  IMDPS   4
#define  BVIWS2  0
#define  BVIMDPS 1
__BEGIN_DECLS
void initc_output_cal_(int *nwpc_,int *ndep_,double*pvkd_,double*pebdep_);
void mean1_(double *ea,int *iacb,int *iace, float *ape ,float *tpf ,float *aet ,float *h1_3);
void mixture_(double *ea,int *iacb,int *iace,float*bv);
void mean1t_(double *ea,int *iac_,float*apet,float*tpft,float*aett,float*h1_3t);
__END_DECLS
void initc_output_cal_(int *nwpc_,int *ndep_,double*pvkd_,double*pebdep_){
  pd.nwpc=*nwpc_;
  pd.pvkd=pvkd_;
  pd.ndep =*ndep_;
  pd.npdep=BVIMDPS+pd.ndep;
  pd.pebdep=pebdep_;
}
#ifdef MMTHREAD
void g_initOutput_cal(){
  gd->ndep=pd.ndep;
  gd->npdep=pd.npdep;
  long sa=(gd->nwpc+1)*(size_t)gd->npdep*kld;
  long sb=(gd->nwpc+1)*(size_t)IMDPS*kl;
  double *p=HNEWN(double,sa+sb+64);
  gd->pebdep=p;p+=sa+32;
  gd->pvkd  =p;
  //write through ape,tpf,aet,h1_3,bv; 
  for(int iac=1;iac<=gd->nwpc;iac++){
    double *d,*s;
    d=gd->pvkd+iac*(size_t)IMDPS*kl;
    s=pd.pvkd+gd->l2gind[iac]*(size_t)IMDPS*kl;
    for(int j=0;j<IMDPS*kl;j+=8){
      VVD(d+j) =VVD(s+j);
    }
    d=gd->pebdep+iac*(size_t)gd->npdep*kl;
    s=pd.pebdep+gd->l2gind[iac]*(size_t)gd->npdep*kl;
    for(int j=0;j<gd->npdep*kl;j+=8){
      VVD(d+j) =VVD(s+j);
    }
  }
}
#endif
void mean1_1(double *eei,double*pvkdt, float *ape ,float *tpf ,float *aet ,float *h1_3){
  ImplschPar *lpp=&gd->pisp;
  double tztp=1.2,tztz=1.099314;
  double zpi=2*3.1415926535897932384626433832795;
  double pi2d=360/zpi;
#ifdef ARM_VEC
  svbool_t pm=svptrue_b64();
  VDOUBLE vaess=VDUP(0.),vawfss=VDUP(0.),vasiss=VDUP(0.),vapess=VDUP(0.),vaets=VDUP(0.),vaetc=VDUP(0.);
  for(int  k=0;k<kl ;k+=NKP){
    VDOUBLE vekjs, vaetst, vaetct;
    vekjs=vaetst=vaetct=VDUP(0.);
    for(int j=0,kj=k;j<jnthet;kj+=kl,j++){
      VDOUBLE vekj=VVD(eei+kj);
      vekjs+=vekj;
      vaetst+=vekj*VVD(lpp->sinths+j);
      vaetct+=vekj*VVD(lpp->cosths+j);
    }
    vaess +=vekjs  *VLD(lpp->dwk          +k);
    vapess+=vekjs  *VLD((pvkdt+kl*IWS2DK )+k);
    vasiss+=vekjs  *VLD((pvkdt+kl*IIWF2DK)+k);
    vawfss+=vekjs  *VLD((pvkdt+kl*IWFDK  )+k);
    VDOUBLE vwkdkk= VLD(lpp->wkdk+k);
    vaets +=vaetst *vwkdkk;
    vaetc +=vaetct *vwkdkk;
  }
  double aets=VSUM(vaets);
  double aetc=VSUM(vaetc);
  double awfss=VSUM(vawfss);
  double asiss=VSUM(vasiss);
  double apess=VSUM(vapess);
  double aess=VSUM(vaess);
#else
  double aess=0,awfss=0,asiss=0,apess=0,aets=0,aetc=0;
  for(int k=0;k<kl ;k++){
    double ekjs=0, aetst=0.,aetct=0.;
    for(int j=0,kj=k;j<jnthet;kj+=kl,j++){
      double ekj=eei[kj];
      ekjs +=ekj;
      aetst+=ekj*lpp->sinths[j];
      aetct+=ekj*lpp->cosths[j];
    }
    aess +=ekjs  *lpp->dwk[k];
    apess+=ekjs  *pvkdt[k+kl*IWS2DK ];
    asiss+=ekjs  *pvkdt[k+kl*IIWF2DK];
    awfss+=ekjs  *pvkdt[k+kl*IWFDK  ];
    double wkdkk=lpp->wkdk[k];
    aets +=aetst *wkdkk;
    aetc =aetc +aetct *wkdkk;
  }
#endif
  //!^^^^^^^^^^^^^^^^^^
  //!ape: tz; tpf: tp; aet: th; h1_3:hs
  *aet=atan2(aets,aetc)*pi2d; //MATH atan2
  if (*aet<0.)*aet=360.+*aet;
  *h1_3=4.*sqrt(aess);
  if(aess>1e-30){
    *ape=tztz*zpi*sqrt(aess/apess);//MATH sqrt
    *tpf=asiss*awfss/(aess*aess);
  }else{
    *ape=0;
    *tpf=0;
  }
}
void mean1(double *ea,int iacb,int iace, float *ape ,float *tpf ,float *aet ,float *h1_3){
  int iac;
  for(iac=iacb;iac<=iace;iac++){
#ifdef MMTHREAD
    int iacg=gd->l2gind[iac];
#else
    int iacg=iac;
#endif
    mean1_1(ea+iac*mkj,gd->pvkd+iac*IMDPS*kld,ape+iacg,tpf+iacg,aet+iacg,h1_3+iacg);
  }
}
void mean1_(double *ea,int *iacb,int *iace, float *ape ,float *tpf ,float *aet ,float *h1_3){
  mean1(ea,*iacb,*iace,ape,tpf,aet,h1_3);
}
void mean1t_(double *ea,int *iac_,float*apet,float*tpft,float*aett,float*h1_3t){
  int iac=*iac_;
  mean1_1(ea+iac*mkj,gd->pvkd+iac*IMDPS*kld,apet,tpft,aett,h1_3t);
}
void mean1t(double *ea,int iac,float*apet,float*tpft,float*aett,float*h1_3t){
  mean1_1(ea+iac*mkj,gd->pvkd+iac*IMDPS*kld,apet,tpft,aett,h1_3t);
}
void mixture_1(double *eei,double *pebdept,float *bvo,int nskip){
  ImplschPar *lpp=&pisp;
  for(int kh=0;kh<gd->ndep;kh++){
#ifdef ARM_VEC
    VDOUBLE vbv1,vbv2,vbv3;
    vbv1=vbv2=vbv3=VDUP(0.);
    svbool_t pm=svptrue_b64();
    for(int k=0;k<kl ;k+=8){
      VDOUBLE vekjs=VDUP(0.0);
      for(int kj=k;kj<mkj;kj+=kl){
        vekjs +=VLD(eei+kj);
      }
      VDOUBLE vebdep=VVD(pebdept+k+(BVIMDPS+kh)*kl);
      VDOUBLE vwsk2=VVD(pebdept+k+BVIWS2*kl);
      VDOUBLE vwkk=VVD(lpp->wk+k);
      VDOUBLE vbvt =vekjs*vebdep;vbv1+=vbvt;
      vbvt*=      vwsk2 ;vbv2+=vbvt;
      vbvt*=      vwkk  ;vbv3+=vbvt;
    }
    double bv1=VSUM(vbv1);
    double bv2=VSUM(vbv2);
    double bv3=VSUM(vbv3);
#else
    double bv1=0,bv2=0,bv3=0;
    for(int k=0;k<kl ;k++){
      double ekjs=0;
      for(int kj=mkj+k-kl;kj>0;kj-=kl){
        ekjs =ekjs +eei[kj];
      }
      double wkk=lpp->wk[k];
      double wsk2=pebdept[k+BVIWS2  *kl];
      double ebdep=pebdept[k+(BVIMDPS+kh)*kl];
      double bvt=ekjs*ebdep  ;bv1=bv1+bvt;
      bvt=bvt*wsk2    ;bv2=bv2+bvt;
      bvt=bvt*wkk     ;bv3=bv3+bvt;
    }
#endif
    if(bv2>1e-30){
      bvo[kh*nskip]=bv1/sqrt(bv2)*bv3;
    }else{
      bvo[kh*nskip]=0;
    }
  }
}
void mixture_(double *ea,int *iacb,int *iace,float*bv){
  int iac;
  for(iac=*iacb;iac<=*iace;iac++){
#ifdef MMTHREAD
    int iacg=gd->l2gind[iac];
#else
    int iacg=iac;
#endif
    mixture_1(ea+iac*mkj,gd->pebdep+iac*gd->npdep*kld,bv+iacg,pd.nwpc+1);
  }
}
void addea(double*ea,double*ee,int iacb,int iace){
  for(int iac=iacb;iac<=iace;iac++){
    double *pea=ea+mkj*iac;
    double *pee=ea+mkj*iac;
#ifdef ARM_VEC
    for(int kj=0;kj<mkj;kj+=8){
      VVD(pea+kj)+=VVD(pee+kj);
    }
#else
    for(int kj=0;kj<mkj;kj++){
      pea[kj]+=pee[kj];
    }
#endif
  }
}
void addea_(double*ea,double*ee,int*iacb,int*iace){ addea(ea,ee,*iacb,*iace); }
void averageea(double*ea,int nea,int iacb,int iace){
  double rv=1./ nea;
#ifdef ARM_VEC
  VDOUBLE vrv=VDUP(rv);
  for(int iac=iacb;iac<=iace;iac++){
    double *pe=ea+mkj*iac;
    for(int kj=0;kj<mkj;kj+=8)VVD(pe+kj)*=vrv;
  }
#else
  for(int iac=iacb;iac<=iace;iac++){
    double *pe=ea+mkj*iac;
    for(int kj=0;kj<mkj;kj++)pe[kj]*=rv;
  }
#endif
}
void averageea_(double*ea,int*nea,int*iacb,int*iace){ averageea(ea,*nea,*iacb,*iace); }
void zeroea_(double*ea,int*iacb,int*iace){
  memset(ea,0,kl*jnthet*(*iace-*iacb+1)*sizeof(double));
}
void zeroea(double*ea,int iacb,int iace){
  memset(ea,0,kl*jnthet*(iace-iacb+1)*sizeof(double));
}
