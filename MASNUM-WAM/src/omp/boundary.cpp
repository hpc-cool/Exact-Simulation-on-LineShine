#include "svars.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
typedef struct _FBoundary{
  double cosths[CJNTHET],sinths[CJNTHET];
  double wk[CKL];
  double windfield;
  int nwpc;
}FBoundary;
__BEGIN_DECLS
void initc_boundary_(FBoundary *bdy_,double*pvws,char *nsp);
void setspec_(int *mode,double *eet,int *iacb,int *iace,int *type) ;
void setspec(int mode,double *e,int nwpb,int nwpe,int n) ;
__END_DECLS

void initc_boundary_(FBoundary *bdy_,double*pvws,char *nsp){
  int j;
  for(j=0;j<CJNTHET;j++)pd.bdy.cosths[j]=bdy_->cosths[j];
  for(j=0;j<CJNTHET;j++)pd.bdy.sinths[j]=bdy_->sinths[j];
  for(j=0;j<CKL    ;j++)pd.bdy.rwk[j]=1./bdy_->wk[j];
  pd.bdy.windfield=bdy_->windfield;
  pd.bdy.nwpc=bdy_->nwpc;
  pd.bdy.g=9.81;
  pd.bdy.rg=1./9.81;
  pd.bdy.gama=3.3;
  pd.bdy.zpi=2*3.1415926535897932384626433832795;
  pd.bdy.xj0=pd.bdy.windfield*1000.;
  pd.bdy.xj_=pd.bdy.g*pd.bdy.xj0;
  pd.bdy.arlfa_=0.076*2/pd.bdy.zpi;
  pd.bdy.wsj_=22.*pd.bdy.g;
  pd.bdy.F_1=1.;
  pd.bdy.F_m1d25=-1.25;
  pd.bdy.F_m0d5=-0.5;
  pd.bdy.rsigma1=1/0.07;
  pd.bdy.rsigma2=1./0.09;
  pd.bdy.F_m0d4=-0.4;
  pd.bdy.F_m0d33=-0.33;
  pd.bdy.wv0=0.9;
  pd.bdy.pvws=pvws;
}
#ifdef MMTHREAD
void g_initSetspec(){
  memcpy(&gd->bdy,&pd.bdy,sizeof(gd->bdy));
  gd->bdy.pvws=(double*)HNEWN(double,kl*(gd->nwpc+1));
  gd->bdy.nwpc=gd->nwpc;
  double *d,*s;
  for(int i=1;i<=gd->nwpc;i++){
    int l2g=gd->l2gind[i];
    s=pd.bdy.pvws+kl*l2g;
    d=gd->bdy.pvws+kl*i;
    VVD(d+ 0)=VVD(s+ 0);
    VVD(d+ 8)=VVD(s+ 8);
    VVD(d+16)=VVD(s+16);
  }				
}
#endif
#ifdef ARM_VEC
#include "ksvml.h"
#if 0
__attribute__ ((optnone)) void calwl(DOUBLE vx,DOUBLE vy,DOUBLE rwv,double *wl){
  VDOUBLE vvx,vvy,vvt,vc,vs,vrwv;
  svbool_t pm=svptrue_b64();
  vvx=VDUP(vx);       vvy=VDUP(vy); 
  vrwv=VDUP(rwv);         vrwv=vrwv*vrwv;
  vc=VLD(gd->bdy.cosths);     vs=VLD(gd->bdy.sinths);
  vvt=(vc*vvx+vs*vvy);    
  VST(wl,vvt*vvt*vrwv);
  vc=VLD(gd->bdy.cosths+NKP); vs=VLD(gd->bdy.sinths+NKP);
  vvt=(vc*vvx+vs*vvy);    
  VST(wl+NKP,vvt*vvt*vrwv);
}
#endif
//风浪谱，n=1 全场初始化，n=2 边界
#ifdef ALLLOG
#define RINFO(TAG) write(*,'(i8," ",a,f8.3," ",2i8)')iwalltime(),TAG,Difftimer(2),mpi_id,mpi_npe;  call flush(6);
#else
#define RINFO(TAG) if(mpi_id==0)then;write(*,'(i8," ",a,f8.3," ",i8)')iwalltime(),TAG,Difftimer(2),mpi_npe;  call flush(6);endif
#endif
void setspec(int mode,double *e,int nwpb,int nwpe,int n) {
  int iac;
  GData *ppd=&pd;
#ifdef MMTHREAD
  if(mode)ppd=gd;
#endif
  Boundary *pbdy=&ppd->bdy;
  for(iac=nwpb;iac<=nwpe;iac++){     //CALA NP_BOUNDARY*(pm NK*14+NJ*5+NK*NJ*1 exp NK*2 pow MK*1)
    if(ppd->nsp[iac]<n)continue;
    windvs*wvs=ppd->pwvs+iac;
    double wv=wvs->wv;
    double*pvws=pbdy->pvws+iac*kl;

    if (wv<=0) wv=pbdy->wv0;
    double rwv=pbdy->F_1/wv;          //CAL pm 0 div  1 
    double xj=pbdy->xj_*rwv*rwv;      //CAL pm 2 
    double arlfa=pbdy->arlfa_*pow(xj,pbdy->F_m0d4);//CAL pm 1 pow 1 //MATH pow
    double wsj=pbdy->wsj_*pow(xj,pbdy->F_m0d33)*rwv;//CAL pm 2 pow 1 //MATH pow
    double wkj=wsj*wsj*pbdy->rg;//CAL pm 2
    double rwsj=pbdy->F_1/wsj;//CAL pm 0 div1
    int j,k,kj;
    double alpha,rwk0,wsk,rsigma,vt,vt1,vt2;
    double wl[jnthet+NKP];
    svbool_t pm=svptrue_b64();
    VDOUBLE vrwsj=VDUP(rwsj);
#if 0
    calwl(wvx->wx,wvx->wy,rwv,wl);
#else
    {
      VDOUBLE vvx,vvy,vvt,vc,vs,vrwv;
      vvx=VDUP(wvs->wx);       vvy=VDUP(wvs->wy); 
      vrwv=VDUP(rwv);         vrwv=vrwv*vrwv;
      for(j=0;j<jnthet;j+=NKP){
        vc=VLD(pbdy->cosths+j);     vs=VLD(pbdy->sinths+j);
        vvt=(vc*vvx+vs*vvy);   //CAL vpm 3  
        VST(wl+j,vvt*vvt*vrwv);   //CAL vpm 2
      }
    }
#endif

    double *eei=e+mkj*iac;
    for (k=0;k<kl;k+=NKP){

      VDOUBLE vrwk0,vwsk,va1,vvt,vvt1,vvt2,vrsigma;
      vrwk0=VLD(pbdy->rwk+k);
      vwsk=VLD(pvws+k);
      svbool_t pq=svcmple(pm,vwsk,VDUP(wsj));
      vrsigma=svsel_f64(pq,VDUP(pbdy->rsigma1),VDUP(pbdy->rsigma2));
      //alpha=arlfa/wk0**4 * exp(-1.25*(wkj/wk0)**2) * gama**(exp(-0.5*((1.-wsk/wsj)/sigma)**2))*(wl/wv)**2
      vvt= VDUP(wkj)*vrwk0;
      vvt=VEXP(pbdy->F_m1d25*vvt*vvt);//MATH vexp
      vrwk0=vrwk0*vrwk0*vrwk0*vrwk0;
      vvt1=vrsigma-vrsigma*vwsk*vrwsj;
      vvt1=VEXP(pbdy->F_m0d5*vvt1*vvt1);//MATH vexp
      va1=VDUP(arlfa)*(vrwk0)*vvt*VPOW(VDUP(pbdy->gama),vvt1);//MATH vpow
      for (j=0,kj=k;j<jnthet;j++,kj+=kl){
        VST(eei+kj,VDUP(wl[j])*va1);
      } 
    }
  }
}
#else
//风浪谱，n=1 全场初始化，n=2 边界
void setspec(int mode,double *e,int nwpb,int nwpe,int n) {
  int iac;
  GData *ppd=&pd;
#ifdef MMTHREAD
  if(mode)ppd=gd;
#endif
  Boundary *pbdy=&ppd->bdy;
  for(iac=nwpb;iac<=nwpe;iac++){     //CALA NP_BOUNDARY*(pm NK*14+NJ*5+NK*NJ*1 exp NK*2 pow MK*1)
    if(gd->nsp[iac]<n)continue;
    windvs*wvs=gd->pwvs+iac;
    double vx=wvs->wx;
    double vy=wvs->wy;
    double wv=wvs->wv;
    double*pvws=pbdy->pvws+iac*kl;

    if (wv<=0) wv=pbdy->wv0;
    double rwv=pbdy->F_1/wv;          //CAL pm 0 div  1 
    double xj=pbdy->xj_*rwv*rwv;      //CAL pm 2 
    double arlfa=pbdy->arlfa_*pow(xj,pbdy->F_m0d4);//CAL pm 1 pow 1
    double wsj=pbdy->wsj_*pow(xj,pbdy->F_m0d33)*rwv;//CAL pm 2 pow 1
    double wkj=wsj*wsj*pbdy->rg;//CAL pm 2
    double rwsj=pbdy->F_1/wsj;//CAL pm 0 div1
    int j,k,kj;
    double alpha,rwk0,rsigma,vt,vt1,vt2;
    double wl[jnthet];
    double *eei=e+mkj*iac;
    for (j=0;j<jnthet;j++){
      double vt;
      vt=(vx*pbdy->cosths[j]+vy*pbdy->sinths[j])*rwv;//CAL pm NJ*5
      wl[j]=vt*vt;
    }
    for (k=0;k<kl;k++){
      double vrwk0,wsk,va1,vvt,vvt1,vvt2,vrsigma;
      vrwk0=pbdy->rwk[k];
      wsk =pvws[k];
      if (wsk<=wsj) vrsigma=pbdy->rsigma1; // 14.2857142857;//1/0.07;
      else vrsigma=pbdy->rsigma2; // 11.1111111111;//1/0.09;
      //alpha=arlfa/wk0**4 * exp(-1.25*(wkj/wk0)**2) * gama**(exp(-0.5*((1.-wsk/wsj)/sigma)**2))*(wl/wv)**2
      vvt= wkj*vrwk0;//CAL pm NK*1 
      vvt=exp(pbdy->F_m1d25*vvt*vvt);//CAL pm NK*2 exp NK*1
      vrwk0=vrwk0*vrwk0*vrwk0*vrwk0;//CAL pm NK*3
      vvt1=vrsigma-vrsigma*wsk*rwsj;//CAL pm NK*3
      vvt1=exp(pbdy->F_m0d5*vvt1*vvt1);//CAL pm NK*2 exp NK*1
      va1=arlfa*vrwk0*vvt*pow(pbdy->gama,vvt1);//CAL pm NK*3 pow NK*1
      for (j=0,kj=k;j<jnthet;j++,kj+=kl){
        eei[kj]=wl[j]*va1;//CAL pm NK*NJ*1
      } 
    }
  }
}
#endif
void setspec_(int*mode,double *e,int *nwpb,int*nwpe,int*n) {
  setspec(*mode,e,*nwpb,*nwpe,*n) ;
}
