#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#ifndef __THREAD
#include "svars.h"
#endif
#ifdef ARM_VEC
#define VEC_XNL
#endif
typedef struct wws_type {
  int ik2,ia2;      //lower wave number,direction  index of k2
  int ik4,ia4;      //lower wave number,direction  index of k2
  float wr_k2,wr_k4;
  float w1k2,w2k2,w3k2,w4k2; //weights  of k2
  float w1k4,w2k4,w3k4,w4k4; //weights  of k4
  float zz;         //!base on depth, compound product of cple*ds*sym/jac
} wws_type ;
typedef struct wws_typev {
  short ik2,ia2;      //lower wave number,direction  index of k2
  short ik4,ia4;      //lower wave number,direction  index of k2
#ifndef VEC_XNL
  float wr_k2,wr_k4;
#endif
  float w1k2,w2k2,w3k2,w4k2; //weights  of k2
  float w1k4,w2k4,w3k4,w4k4; //weights  of k4
  float zz;         //!base on depth, compound product of cple*ds*sym/jac
} wws_typev ;
#define min(a,b) a>b?b:a;
typedef struct xnlinf{
  int naq,nkq,nekq;
  int mkq,maq,klocus;
  int qf_dn,qf_kn;
  int  *quadnl;
  wws_typev*quads;
  double*q_k2,*q_ka,*q_sig,*q_sigr;
  double*q_kpow,*q_tail,*q_lamd;
}xnlinf;
xnlinf xi;
#ifdef DEBUG
#define CHKEE(ee,iat) Checkee(ee,iat,__FILE__,__LINE__)
#else
#define CHKEE(ee,iat) 
#endif
#ifndef __THREAD
#define __THREAD
#endif
__THREAD xnlinf *ppxi=&xi;
__THREAD double *nspec=NULL;
#ifndef HMALLOC
#define HMALLOC malloc
#define HNEWZN(p,N) p=(typeof(p))HMALLOC(((long)sizeof(*p))*(N));if(p)memset(p,0,sizeof(*p)*(N))
#endif
#define NXNLKPL NKP
#define NXNLKPR (NKP*2)
#define NXNLKP (NXNLKPL+NXNLKPR)
int getspecsize(){
#ifdef VEC_XNL
  return (xi.naq*(xi.nkq+NXNLKP )+NXNLKPL);
#else
  return (xi.naq*(xi.nkq));
#endif
}
#ifdef MMTHREAD
void init_xnlg_(){
  xnlinf *pxi;
  size_t st=getspecsize();
#ifndef USEHBM
  pxi=ppxi=gd->xi=(xnlinf*)HMALLOC(sizeof(*pxi));
  memcpy(pxi,&xi,sizeof(xi));
#else
  size_t size=0;
  size+=sizeof(xnlinf)+64;
#define VC(v,N) size+=((long)sizeof(*xi.v))*(N)+64
  VC(quads ,xi.klocus*xi.mkq*xi.maq);
  VC(quadnl,xi.mkq*xi.maq);
  VC(q_k2  ,xi.nkq    );VC(q_ka  ,xi.nkq);
  VC(q_sig ,xi.nkq    );VC(q_sigr,xi.nkq);
  VC(q_kpow,xi.nkq    );VC(q_tail,xi.nkq);
  VC(q_lamd,xi.nkq    );
#undef VC
  size+=gi->Nthreads*((st+64)*sizeof(*gd->tee));
  char*p=(char*)HMALLOC(size+64);
  pxi=ppxi=gd->xi=(xnlinf*)p;p+=sizeof(xnlinf)+64;
  memcpy(pxi,&xi,sizeof(xi));
#define VC(v,N) pxi->v=(typeof(xi.v))p; p+=((long)sizeof(*xi.v))*(N)+64; if(pxi->v&&xi.v)memcpy(pxi->v,xi.v,sizeof(*xi.v)*(N))
  VC(quads ,xi.klocus*xi.mkq*xi.maq);
  VC(quadnl,xi.mkq*xi.maq);
  VC(q_k2  ,xi.nkq    );VC(q_ka  ,xi.nkq);
  VC(q_sig ,xi.nkq    );VC(q_sigr,xi.nkq);
  VC(q_kpow,xi.nkq    );VC(q_tail,xi.nkq);
  VC(q_lamd,xi.nkq    );
#undef VC
  gd->tee=(typeof(gd->tee))p;
#endif
}
#endif
/*-3,9 align to vector*/
void init_xnl_(){
  size_t st=getspecsize();
#ifdef MMTHREAD
  ppxi=gd->xi;
  nspec=gd->tee+ti->ind*(st+64);
#else
  ppxi=&xi;
  nspec=(double*)HMALLOC((st+64)*sizeof(double));
#endif
#ifdef VEC_XNL
  nspec+=NXNLKPL; 
#endif
}
void cpwws(wws_typev*qv,wws_type*q){
#define VS(n) qv->n=q->n
  VS(ik2);VS(ia2);
  VS(ik4);VS(ia4);
  VS(w1k2);VS(w2k2);VS(w3k2);VS(w4k2);
  VS(w1k4);VS(w2k4);VS(w3k4);VS(w4k4);
  VS(zz);
#ifndef VEC_XNL
  VS(wr_k2);VS(wr_k4);
#endif
}
void ckiks(wws_type*quads,int*quadnl,int mkq,int maq,int klocus){
  int m0=1000,m1=-1000;
  for(int i=0;i<mkq*maq;i++){
    int nl=quadnl[i];
    wws_type*q=quads+i*klocus;
    for(int l=0;i<nl;i++,q++){
      if(m0>q->ik2)m0=q->ik2;
      if(m0>q->ik4)m0=q->ik4;
      if(m1<q->ik2)m1=q->ik2;
      if(m1<q->ik4)m1=q->ik4;
    }
  }
  if(m0<-NXNLKPL||m1>NXNLKPR){
    printf("xnl vec error:m0=%d < -NXNLKPL=%d || m1:%d> NXNLKPR=%d\n",m0,NXNLKPL,m1,NXNLKPR);
    exit(-1);
  }
  if((NXNLKPL%NKP!=0)||(NXNLKPR%NKP)){
    printf("xnl vec error:NXNLKPL=%d ||  NXNLKPR%d is not aligned NKP\n",NXNLKPL,NXNLKPR);
    exit(-2);
  }
  if(m0+NXNLKPL>NKP||NXNLKPR-m1>NKP){
    printf("xnl vec error:m0=%d NXNLKPL=%d || m1:%d NXNLKPR=%d is too big\n",m0,NXNLKPL,m1,NXNLKPR);
    exit(-2);
  }
}
int xnlcals(xnlinf *pxi);
void xnlsetgrd(int naq,int nkq,int mkq,int maq,int klocus,int qf_dn,int qf_kn,
            struct wws_type*quads,int*quadnl,
            double*q_k2,double*q_ka,double*q_sig,double*q_sigr,
            double*q_kpow,double*q_tail,double*q_lamd){
  //ckiks(quads,quadnl,mkq,maq,klocus);
#define VP(v) xi.v=v
  VP(naq   );VP(nkq   );VP(mkq   );VP(maq   );VP(klocus);
  VP(qf_dn );VP(qf_kn );
#undef VP
  int size=0;
  char *p=NULL;
  size+=klocus*mkq*maq*sizeof(*xi.quads)+64;

#define VC(v,N) size+=((long)sizeof(*xi.v))*(N)+64
  VC(quadnl,xi.mkq*xi.maq);
  VC(q_k2  ,xi.nkq    );VC(q_ka  ,xi.nkq);
  VC(q_sig ,xi.nkq    );VC(q_sigr,xi.nkq);
  VC(q_kpow,xi.nkq    );VC(q_tail,xi.nkq);
  VC(q_lamd,xi.nkq    );
#undef VC
  p=(char*)HMALLOC(size+64);
  xi.quads=(typeof(xi.quads))p;p+=klocus*mkq*maq*sizeof(*xi.quads)+64;
  for(int i=0;i<klocus*mkq*maq;i++)cpwws(xi.quads+i,quads+i);
#define VC(v,N) xi.v=(typeof(xi.v))p; p+=((long)sizeof(*xi.v))*(N)+64;for(int _i=0;_i<N;_i++)xi.v[_i]=v[_i]
  VC(quadnl,xi.mkq*xi.maq);
  VC(q_k2  ,xi.nkq    );VC(q_ka  ,xi.nkq);
  VC(q_sig ,xi.nkq    );VC(q_sigr,xi.nkq);
  VC(q_kpow,xi.nkq    );VC(q_tail,xi.nkq);
  VC(q_lamd,xi.nkq    );
#undef VC
  ppxi=&xi;
  //xnlcals(ppxi);
}
void xnlsetgrd_(int*naq,int*nkq,int*mkq,int*maq,int*klocus,int*qf_dn,int*qf_kn,
             struct wws_type*quads,int*quadnl,
             double*q_k2,double*q_ka,double*q_sig,double*q_sigr,
             double*q_kpow,double*q_tail,double*q_lamd){
  xnlsetgrd(*naq,*nkq,*mkq,*maq,*klocus,*qf_dn,*qf_kn,
         quads,quadnl,q_k2,q_ka,q_sig,q_sigr,q_kpow,q_tail,q_lamd);
}
void DBGH(){
  printf("DBGH\n");
}
static int q_getlocus(xnlinf*pxi,wws_typev**tqs,int kdif,int adif){
  int kmem = abs(kdif) ;
  int amem = adif;
  if(adif>=0) amem=adif;
  else        amem=pxi->maq+adif;
  *tqs=pxi->quads+(kmem+amem*pxi->mkq)*pxi->klocus;
  return pxi->quadnl[kmem+amem*pxi->mkq];
}
#ifdef VEC_XNL
#define INDB(k,j) (k+j*pxi->nkq)
#define NSD(k,j) nspec[k+j*nkqp]
#define NSV(k,j) VVD(nspec+(k+j*nkqp))
int xnlcal(double*ee,int ks,double*xnl,double*diag){
  int ia1,ia3,la3 ; // limits for directional loops
  int ik1,ik3,lk3 ; // counters for wave number loop
  int ibeta, nloc,iloc;
  wws_typev* tq;
  wws_typev*tqs;
  int ja2,ja2p,ja4,ja4p;      // direction indices for interpolation of k2,k4
  int jk2,jk2p,jk4,jk4p;      // wave number indices for interpolation of k2,k4
  VDOUBLE qn1,qn3,qn2,qn4;     // action densities at wave numbers k1, k2, k3 and k4
  VDOUBLE qn13p,qn13d; // N1*N3,N3-N1
  VDOUBLE t13    ;     // value of sub-integral
#ifdef USEDIAG
  VDOUBLE qn24p,qn24d,diagk1,diagk3;     // diagonal term for k1,k3
#endif
  VDOUBLE zz,t2,t4;    // tail factors for k2 and k4
  VDOUBLE z_lambda;
  xnlinf *pxi=ppxi;
  int nkqp=pxi->nkq+NXNLKP;
  int *piac=(int*)diag;
  int iact=*piac;
  for(ia1=0;ia1<pxi->naq;ia1++){
    int ilo=ia1*pxi->nkq;
    for(ik1=0;ik1<pxi->nkq;ik1+=NKP){
      VVD(xnl +ik1+ilo) =VDUP(0.);
#ifdef USEDIAG
      VVD(diag+ik1+ilo) =VDUP(0.);
#endif
    }
    // wr_k2 = k2m/kqmin
    // wr_k2 = q_kfac**(-3.5*(wk_k2-nkq))
    // q_tail(ik1) =q_kfac**(-3.5*ik1)
    // q_kpow(ik1) =q_kfac**((ik1-1))
    // t2 = q_kfac**(-3.5*(wr_k2-nkq))
    // t2 = q_kfac**(-3.5*(jk2-nkq))=q_tail(jk2-nkq+1)
    for(ik1=0;ik1<pxi->nkq;ik1+=NKP){
      NSV(ik1,ia1)=VVD(ee+ik1+ilo)*VVD(pxi->q_sigr+ik1);
    }
    for(ik1=0;ik1<NXNLKPR;ik1+=NKP){ // q_tail(ik1) =q_kfac**(-3.5*ik1)
      NSV(ik1+pxi->nkq,ia1)=NSD(pxi->nkq-1,ia1)*VVD(pxi->q_tail+ik1+1);
    }
    for(ik1=-NXNLKPL;ik1<0;ik1+=NKP){ 
      // q_kpow(ik1) =q_kfac**((ik1-1))
      // h2 = exp((wr_k2-1)*alog(q_kfac))=q_kfac**(wr_k2-1)
      // h2 = q_kfac**((ik2-1))*q_kpow[ik1]=q_kfac**(jk2-1)
      NSV(ik1,ia1)=NSD(0,ia1)*VVD(pxi->q_kpow-ik1+2);
    }
  }
  int ssloc=0;
  for(lk3=0;lk3<pxi->mkq;lk3++){
    for(la3 = -pxi->qf_dn;la3<=pxi->qf_dn;la3++){ // loop over all possible wave directions
      nloc=q_getlocus(pxi,&tqs,lk3,la3);
      if(nloc<=0)continue;
      if(lk3==0 && la3 ==0)continue;  // skip routine if k1=k3
      for(ia1= 0;ia1<pxi->naq;ia1++){      // loop over selected part of grid, set in q_init
        ia3=(ia1+la3+pxi->naq)%pxi->naq;        // cycle
        ibeta=ia1; //if(lk3<0) ibeta=ia3;
        for(ik1 = 0;ik1<pxi->nkq-lk3;ik1+=NKP){
          ssloc+=nloc;
          ik3=ik1+lk3;
          // assume ik1<ik3
          z_lambda=VVD(pxi->q_kpow+ik1);/*= q_kfac**((ik1)*7.5)*/
          // 11 33 ; 13 31 ; 31 13 ; 33 11
          // perform integration along locus
          //call q_t13v4(nspec,ik1,ia1,ik3,ia3,t13,diagk1,diagk3)   //jiangxj:diag
          t13    = VDUP(0.); 
#ifdef USEDIAG
          diagk1 = VDUP(0.); diagk3 = VDUP(0.);
#endif
          qn1 = NSV(ik1,ia1);
          qn3 = NSV(ik3,ia3);
          qn13p = qn1*qn3; qn13d = qn3-qn1;

          //-----------------------------------------------------------------------------------
          //  loop over the locus
          for(iloc=0;iloc<nloc;iloc++){
            tq=tqs+iloc; zz=VDUP(tq->zz*2);
            ja2 = (tq->ia2+ibeta  +pxi->naq)%pxi->naq;
            ja2p= (tq->ia2+ibeta+1+pxi->naq)%pxi->naq;
            ja4 = (tq->ia4+ibeta  +pxi->naq)%pxi->naq;
            ja4p= (tq->ia4+ibeta+1+pxi->naq)%pxi->naq;
            jk2 = tq->ik2+ik1; jk4 = tq->ik4+ik1;//-3<= ik2,ik4<=9
            qn2 = (tq->w1k2*NSV(jk2,ja2 )+tq->w2k2*NSV(jk2+1,ja2 )+
                   tq->w3k2*NSV(jk2,ja2p)+tq->w4k2*NSV(jk2+1,ja2p) );
            qn4 = (tq->w1k4*NSV(jk4,ja4 )+tq->w2k4*NSV(jk4+1,ja4 )+
                   tq->w3k4*NSV(jk4,ja4p)+tq->w4k4*NSV(jk4+1,ja4p) );
            t13+= zz*z_lambda*(qn13p*(qn4-qn2) + qn2*qn4*qn13d);
#ifdef USEDIAG
            qn24d=qn4-qn2;qn24p=qn2*qn4;
            diagk1+= (qn3*qn24d - qn24p)*zz;
            diagk3+= (qn1*qn24d + qn24p)*zz;
#endif
          }
          //end q_t13v4
          VVD(xnl +ik1+ia1*pxi->nkq) +=VVD(pxi->q_k2+ik3)*t13   ;//user mannual (2.33)
          VVD(xnl +ik3+ia3*pxi->nkq) -=VVD(pxi->q_k2+ik1)*t13   ;//move *2 to zz

#ifdef USEDIAG                                                  
          VVD(diag+ik1+ia1*pxi->nkq) +=VVD(pxi->q_ka+ik3)*diagk1;
          VVD(diag+ik3+ia3*pxi->nkq) -=VVD(pxi->q_ka+ik1)*diagk3;
#endif
        }
      }
    }
  }
  for(ia1=0;ia1<pxi->naq;ia1++){
    for(ik1=0;ik1<pxi->nkq;ik1+=NKP){
      VVD(xnl+ik1+ia1*pxi->nkq)=VVD(xnl+ik1+ia1*pxi->nkq)*VVD(pxi->q_sig+ik1);//q_dfac*
    }
  }
  return 0;
}
#else
#define INDB(k,j) (k+j*pxi->nkq)
#define NSV(k,j) nspec[(k+j*pxi->nkq)]
void Checkee(DOUBLE*ee,int iac,const char*fn,int ln);
#ifdef DEBUG
#define EC(ae) if(isnan(ae)){ printf("ee is nan %s %d %d\n",__FILE__,__LINE__,iact); exit(-121); } 
#else
#define EC(ae) 
#endif
extern "C" int iact;
int xnlcal(double*ee,int ks,double*xnl,double*diag){
  int ia1,ia3,la3 ; // limits for directional loops
  int ik1,ik3,lk3 ; // counters for wave number loop
  int ibeta, nloc,iloc;
  wws_typev* tq;
  wws_typev*tqs;
  int ja2,ja2p,ja4,ja4p;      // direction indices for interpolation of k2,k4
  int jk2,jk2p,jk4,jk4p;      // wave number indices for interpolation of k2,k4
  double qn1,qn3,qn2,qn4;     // action densities at wave numbers k1, k2, k3 and k4
  double qn13p,qn13d; // N1*N3,N3-N1
  double t13    ;     // value of sub-integral
#ifdef USEDIAG
  double qn24p,qn24d,diagk1,diagk3;     // diagonal term for k1,k3
#endif
  double zz,t2,t4;    // tail factors for k2 and k4
  double z_lambda;
  xnlinf *pxi=ppxi;
  for(ia1=0;ia1<pxi->naq;ia1++){
    for(ik1=0;ik1<pxi->nkq;ik1++){
      NSV(ik1,ia1)=ee[ik1+ia1*pxi->nkq]*pxi->q_sigr[ik1];
    }
  }
  memset(xnl,0,sizeof(xnl[0])*pxi->naq*pxi->nkq);
  for(lk3=0;lk3<pxi->mkq;lk3++){
    for(la3 = -pxi->qf_dn;la3<=pxi->qf_dn;la3++){ // loop over all possible wave directions
      nloc=q_getlocus(pxi,&tqs,lk3,la3);
      if(nloc<=0)continue;
      if(lk3==0 && la3 ==0)continue;  // skip routine if k1=k3
      for(ia1= 0;ia1<pxi->naq;ia1++){      // loop over selected part of grid, set in q_init
        ia3=(ia1+la3+pxi->naq)%pxi->naq;        // cycle
        ibeta=ia1; //if(lk3<0) ibeta=ia3;
        for(ik1 = 0;ik1<pxi->nkq-lk3;ik1++){
          ik3=ik1+lk3;
          // assume ik1<ik3
          z_lambda=pxi->q_kpow[ik1];/*= q_kfac**((ik1)*7.5)*/
          // 11 33 ; 13 31 ; 31 13 ; 33 11
          // perform integration along locus
          //call q_t13v4(nspec,ik1,ia1,ik3,ia3,t13,diagk1,diagk3)   //jiangxj:diag
          t13    = 0.; 
#ifdef USEDIAG
          diagk1 = 0.; diagk3 = 0.;
#endif
          qn1 = NSV(ik1,ia1);
          qn3 = NSV(ik3,ia3);
          qn13p = qn1*qn3; qn13d = qn3-qn1;
          //-----------------------------------------------------------------------------------
          //  loop over the locus
          for(iloc=0;iloc<nloc;iloc++){
            tq=tqs+iloc;
            zz=tq->zz;
            ja2 = tq->ia2+ibeta;ja2p = (ja2+1+pxi->naq)%pxi->naq;ja2  = (ja2+pxi->naq)%pxi->naq;
            ja4 = tq->ia4+ibeta;ja4p = (ja4+1+pxi->naq)%pxi->naq;ja4  = (ja4+pxi->naq)%pxi->naq;
            jk2 = tq->ik2+ik1;
            if(jk2>=pxi->nkq-1){
              t2 = tq->wr_k2* pxi->q_tail[ik1];//wr_k2 = q_kfac**(-3.5*(pg%wk_k2-pxi->nkq))
              qn2 =((tq->w1k2+tq->w2k2)*NSV(pxi->nkq-1,ja2)+
                    (tq->w3k2+tq->w4k2)*NSV(pxi->nkq-1,ja2p))*t2;
            }else if(jk2>=0){
              jk2p = min(jk2+1,pxi->nkq);
              qn2 = (tq->w1k2*NSV(jk2,ja2 )+tq->w2k2*NSV(jk2p,ja2 )+
                     tq->w3k2*NSV(jk2,ja2p)+tq->w4k2*NSV(jk2p,ja2p) );
            }else{
              t2  = tq->wr_k2*pxi->q_kpow[ik1]; //wr_k2 = k2m/kqmin
              qn2 =((tq->w1k2+tq->w2k2)*NSV(0,ja2)+(tq->w3k2+tq->w4k2)*NSV(0,ja2p))*t2;
            }
            jk4 = tq->ik4+ik1;
            if(jk4>=pxi->nkq-1){
              t4  = tq->wr_k4* pxi->q_tail[ik1];
              qn4 =((tq->w1k4 +tq->w2k4)*NSV(pxi->nkq-1,ja4)+
                    (tq->w3k4+tq->w4k4)*NSV(pxi->nkq-1,ja4p))*t4;
            }else if(jk4>=0){
              jk4p = jk4+1;
              qn4 = (tq->w1k4*NSV(jk4,ja4 )+tq->w2k4*NSV(jk4p,ja4 )+
                     tq->w3k4*NSV(jk4,ja4p)+tq->w4k4*NSV(jk4p,ja4p) );
            }else{
              t4  = tq->wr_k4*pxi->q_kpow[ik1];
              qn4 = ((tq->w1k4+tq->w2k4)*NSV(0,ja4) + (tq->w3k4 + tq->w4k4)*NSV(0,ja4p))*t4;
            }
            t13    = t13 + zz*z_lambda*(qn13p*(qn4-qn2) + qn2*qn4*qn13d);
#ifdef USEDIAG
            qn24d=qn4-qn2;qn24p=qn2*qn4;
            diagk1 = diagk1 + (qn3*qn24d - qn24p)*zz;
            diagk3 = diagk3 + (qn1*qn24d + qn24p)*zz;
#endif
          }
          //end q_t13v4
          xnl [INDB(ik1,ia1)] +=2*t13*pxi->q_k2[ik3] ;         //user mannual (2.33)
          xnl [INDB(ik3,ia3)] -=2*t13*pxi->q_k2[ik1] ;
#ifdef USEDIAG
          diag[INDB(ik1,ia1)] +=2*diagk1*pxi->q_ka[ik3];
          diag[INDB(ik3,ia3)] -=2*diagk3*pxi->q_ka[ik1];
#endif
        }
      }
    }
  }
  //for(ia1=0;ia1<pxi->naq;ia1++){
  //  for(ik1=0;ik1<pxi->nkq;ik1++){
  //    printf("%10.2e", xnl[ik1+ia1*pxi->nkq]);
  //  }
  //  printf("\n");
  //}
  //printf("\n");
  for(ia1=0;ia1<pxi->naq;ia1++){
    for(ik1=0;ik1<pxi->nkq;ik1++){
      xnl[ik1+ia1*pxi->nkq]=xnl[ik1+ia1*pxi->nkq]*pxi->q_sig[ik1];//q_dfac*
    }
  }
  return 0;
}
#endif
void xnlcal_(double*nspec,int*ks,double*xnl,double*diag, int*ierror){
  *ierror=xnlcal(nspec,*ks,xnl,diag);
}
int xnlcals(xnlinf *pxi){
  int ia1,ia3,la3 ; // limits for directional loops
  int ik1,ik3,lk3 ; // counters for wave number loop
  int ibeta, nloc,iloc;
  wws_typev* tq;
  wws_typev*tqs;
  int ssloc=0;
  for(lk3=0;lk3<pxi->mkq;lk3++){
    for(la3 = -pxi->qf_dn;la3<=pxi->qf_dn;la3++){ // loop over all possible wave directions
      nloc=q_getlocus(pxi,&tqs,lk3,la3);
      printf("VVA %2d %2d %3d %3d %3d\n",la3,lk3,nloc,pxi->mkq,pxi->qf_dn);
      if(nloc<=0)continue;
      if(lk3==0 && la3 ==0)continue;  // skip routine if k1=k3
      for(ia1= 0;ia1<pxi->naq;ia1++){      // loop over selected part of grid, set in q_init
        ia3=(ia1+la3+pxi->naq)%pxi->naq;        // cycle
        for(ik1 = 0;ik1<pxi->nkq-lk3;ik1+=NKP){
          ssloc+=nloc;
          printf("VVB %2d %2d %2d %2d %3d %8d\n",ik1,ia1,la3,lk3,nloc,ssloc);
        }
      }
    }
  }
  printf("VVC %8d\n",ssloc);
  exit(0);
}
