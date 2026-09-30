#include <math.h>
#include <string.h>
#include <stdio.h>
#include "svars.h"
//#define NODELTEE
#define LOGID 0
#define DBGINFN printf("%s %d %d %d:",__FILE__,__LINE__,itt,iat);
void Checkee(double*ee,int iac,const char*fn,int ln);
#ifdef DEBUG
#define CHKEE(ee,iat) Checkee(ee,iat,__FILE__,__LINE__)
#else
#define CHKEE(ee,iat) 
#endif
#define CNAN(v,i1,i2,i3,i4) if(isnan(v)){printf("NAN %s %d %d %d %d:%s %d\n",#v ,i1,i2,i3,i4,__FILE__,__LINE__);}
void mean2(ImplschPar *lpp,DOUBLE*ee,ImplschP*pi,DOUBLE *sds,DOUBLE*awk) ;
void mean3(ImplschPar *lpp,DOUBLE*ee,ImplschP*pi) ;
void InitCImplsch(){
}
#ifdef MMTHREAD
void g_initImplsch(){
  // ZZZ fill groups gd->pisp
  memcpy(&gd->pisp,&pd.pisp,sizeof(gd->pisp));
  gd->pisp.pip=HNEWN(ImplschP,(gd->nwpc+1));
  ImplschP*pip=gd->pisp.pip;
  ImplschP*pipn=pd.pisp.pip;
  for(int i=1;i<=gd->nwpc;i++){
    pip[i]=pipn[gd->l2gind[i]];
  }
}
#endif
#define register 
#define NKP 8
#define KLP (CKL+NKP)
__THREAD double*tee=NULL;
void init_dia_(){
  tee=HNEWZN(tee, KLP*CJNTHET+NKP+64); 
#ifndef ARM_VEC
  tee+=NKP;
#endif
}
#ifndef ARM_VEC
#define  PWK(I0)   \
    { wka00=lpp->wkal[I0*4  ];wka01=lpp->wkal[I0*4+1];wka02=lpp->wkal[I0*4+2];wka03=lpp->wkal[I0*4+3];\
      wka10=lpp->wkal[I0*4+4];wka11=lpp->wkal[I0*4+5];wka12=lpp->wkal[I0*4+6];wka13=lpp->wkal[I0*4+7];\
      register int lt0,lt1,ii,kt,jt; ZBYTE *kjs=lpp->kjs+j*48+I0*4;\
      kt=kjs[0];jt=kjs[1]; lt0=jt*KLP+kt;\
      kt=kjs[2];jt=kjs[3]; lt1=jt*KLP+kt;\
      pe00=tee+NKP+lt0;pe01=tee+NKP+lt0+1;pe02=tee+NKP+lt1;pe03=tee+NKP+lt1+1;\
      kjs+=4;\
      kt=kjs[0];jt=kjs[1]; lt0=jt*KLP+kt;\
      kt=kjs[2];jt=kjs[3]; lt1=jt*KLP+kt;\
      pe10=tee+NKP+lt0;pe11=tee+NKP+lt0+1;pe12=tee+NKP+lt1;pe13=tee+NKP+lt1+1;\
    }
#define INTERE12(coff1,coff2) \
      e1=coff1*( wka00*pe00[k]+wka01*pe01[k]+wka02*pe02[k]+wka03*pe03[k]); \
      e2=coff2*( wka10*pe10[k]+wka11*pe11[k]+wka12*pe12[k]+wka13*pe13[k]);
//CAL pm 8+16
#define CALC12(I0)                                               \
    PWK(I0);                                                     \
    for(k=0;k<CKL;k++){                                         \
      register DOUBLE e0,e1,e2;e0=pe0[k]; INTERE12(lmdpd,lmdmd);  \
      pse[k]=pse[k]+LIV_m2*e0*(e0*e1+e0*e2+e1*e2);               \
    }
//CAL pm 9+16
#define CALC34(I0)                                               \
    PWK(I0);                                                     \
    for(k=0;k<CKL;k++){                                         \
      register DOUBLE e0,e1,e2;e0=lmdpd*pe0[k]; INTERE12(LIV_1,lmdmd);  \
      pse[k]=pse[k]+lmdpd20*e1*(e1*e0+e1*e2+e0*e2);   \
    }
//CAL pm 9+16
#define CALC56(I0)                                               \
    PWK(I0);                                                     \
    for(k=0;k<CKL;k++){                                         \
      register DOUBLE e0,e1,e2;e0=lmdmd*pe0[k]; INTERE12(lmdpd,LIV_1);  \
      pse[k]=pse[k]+lmdmd20*e2*(e2*e0+e2*e1+e0*e1);   \
    }
//CALC pm NJ*NK*150
void winter(ImplschPar *lpp,DOUBLE *ee,DOUBLE *se){
  register DOUBLE lmdpd,lmdmd,lmdpd20,lmdmd20,LIV_m2,LIV_1;
  int j;
  if(1){
    DOUBLE *ped=tee,*pes=ee;
    int k,j;
    for(k=0;k<NKP;k++) *ped++=0;
    for(j=0;j<CJNTHET;j++){
      for(k=0;k<CKL;k++) *ped++=*pes++;
      for(k=0;k<NKP;k++) *ped++=0;
    }
  }
  LIV_m2=-2.;LIV_1=1.;
  lmdpd=lpp->lmdpd;lmdmd=lpp->lmdmd;lmdpd20=lpp->lmdpd20;lmdmd20=lpp->lmdmd20;
  memset(se,0,sizeof(DOUBLE)*mkj);
  for(j=jnthet-1;j>=0;j--){
    register DOUBLE wka00,wka01,wka02,wka03;
    register DOUBLE wka10,wka11,wka12,wka13;
    register DOUBLE*pe00,*pe01,*pe02,*pe03;
    register DOUBLE*pe10,*pe11,*pe12,*pe13;
    register DOUBLE*pse;
    DOUBLE *pe0=tee+NKP+j*KLP;
    register int k;
    pse=se+j*(CKL);
    //CAL pm 8*2+9*4+16*6=148*NJ*NK
    CALC12( 0);//* resonate 1 // kjoff : -2,-4   (16+8)*kl
    CALC12( 2);//* resonate 2 // kjoff : -2,-4   (16+8)*kl
    CALC34( 4);//* resonate 3 // kjoff : -3,-6   (16+9)*kl
    CALC34( 6);//* resonate 4 // kjoff : -3,-6   (16+9)*kl
    CALC56( 8);//* resonate 5 // kjoff :  3, 5   (16+9)*kl
    CALC56(10);//* resonate 6 // kjoff :  3, 5   (16+9)*kl
  }
}
#else
#define LDO(P,l) VLD(pm,P+kpos[l*8])
#define LDWK(O)  VDUP(lpp->wkal[O])
#define  PWK(I0){\
    wka00=LDWK(I0*4+0);wka01=LDWK(I0*4+1);wka02=LDWK(I0*4+2);wka03=LDWK(I0*4+3);\
    wka10=LDWK(I0*4+4);wka11=LDWK(I0*4+5);wka12=LDWK(I0*4+6);wka13=LDWK(I0*4+7);\
    register int lt0,lt1,ii,kt,jt; \
    ZBYTE *kjs=lpp->kjs+j*48+I0*4;\
    kt=kjs[0];jt=kjs[1]; lt0=jt*KLP+kt;\
    kt=kjs[2];jt=kjs[3]; lt1=jt*KLP+kt;\
    pe00=tee+NKP+lt0;pe01=tee+NKP+lt0+1;\
    pe02=tee+NKP+lt1;pe03=tee+NKP+lt1+1;\
    kjs+=4;\
    kt=kjs[0];jt=kjs[1]; lt0=jt*KLP+kt;\
    kt=kjs[2];jt=kjs[3]; lt1=jt*KLP+kt;\
    pe10=tee+NKP+lt0;pe11=tee+NKP+lt0+1;\
    pe12=tee+NKP+lt1;pe13=tee+NKP+lt1+1;\
  }
// e1=MUL(coff1,MLA(MLA(MLA(MUL( wka00,LD(pe00+k)),wka01,LD(pe01+k)),wka02,LD(pe02+k)),wka03,LD(pe03+k))); \
// e2=MUL(coff2,MLA(MLA(MLA(MUL( wka10,LD(pe10+k)),wka11,LD(pe11+k)),wka12,LD(pe12+k)),wka13,LD(pe13+k)));
#define INTERE12(coff1,coff2) \
      e1=coff1*( wka00*VVD(pe00+k)+wka01*VVD(pe01+k)+wka02*VVD(pe02+k)+wka03*VVD(pe03+k)); \
      e2=coff2*( wka10*VVD(pe10+k)+wka11*VVD(pe11+k)+wka12*VVD(pe12+k)+wka13*VVD(pe13+k));
#define CALC12(I0)                                                 \
    PWK(I0);                                                       \
    for(k=0;k<CKL;k+=NKP){                                          \
      register VDOUBLE e0,e1,e2;e0=VVD(pe0+k);INTERE12(lmdpd,lmdmd);\
      VVD(pse+k)+=LIV_m2*e0*(e1*e2+e0*(e1+e2));                     \
    }
#define CALC34(I0)                                                  \
    PWK(I0);                                                        \
    for(k=0;k<CKL;k+=NKP){                                          \
      register VDOUBLE e0,e1,e2;e0=VVD(pe0+k);INTERE12(LIV_1,lmdmd);\
      VVD(pse+k)+=lmdpd20*e1*(e1*e2+lmdpd*e0*(e1+e2));              \
    }
#define CALC56(I0)                                                  \
    PWK(I0);                                                        \
    for(k=0;k<CKL;k+=NKP){                                          \
      register VDOUBLE e0,e1,e2;e0=VVD(pe0+k);INTERE12(lmdpd,LIV_1);\
      VVD(pse+k)+=lmdmd20*e2*(e2*e1+lmdmd*e0*(e1+e2));              \
    }

void winter(ImplschPar *lpp,DOUBLE *ee,DOUBLE *se){
  register VDOUBLE lmdpd,lmdmd,lmdpd20,lmdmd20,LIV_m2,LIV_1;
  int j;
  svbool_t pm=svptrue_b64();
  if(1){
    DOUBLE *ped=tee,*pes=ee;
    int k,j;
#if 1
    VVD(tee)=VDUP(0.);
    for(j=0;j<CJNTHET;j++){
      VVD(tee+j*(CKL+NKP)+NKP)=VVD(ee+j*CKL);
      VVD(tee+j*(CKL+NKP)+2*NKP)=VVD(ee+j*CKL+NKP);
      VVD(tee+j*(CKL+NKP)+3*NKP)=VVD(ee+j*CKL+2*NKP);
      VVD(tee+j*(CKL+NKP)+4*NKP)=VVD(ee+j*CKL+3*NKP);
      VVD(tee+j*(CKL+NKP)+5*NKP)=VDUP(0.);
    }
#else
    for(k=0;k<NKP;k++) *ped++=0;
    for(j=0;j<CJNTHET;j++){
      for(k=0;k<CKL;k++) *ped++=*pes++;
      for(k=0;k<NKP;k++) *ped++=0;
    }
#endif
  }
  LIV_m2=VDUP(-2.);LIV_1=VDUP(1.);
  lmdpd=VDUP(lpp->lmdpd);lmdmd=VDUP(lpp->lmdmd);lmdpd20=VDUP(lpp->lmdpd20);lmdmd20=VDUP(lpp->lmdmd20);
  memset(se,0,sizeof(DOUBLE)*mkj);
  for(j=jnthet-1;j>=0;j--){
    register VDOUBLE wka00,wka01,wka02,wka03;
    register VDOUBLE wka10,wka11,wka12,wka13;
    register DOUBLE*pe00,*pe01,*pe02,*pe03;
    register DOUBLE*pe10,*pe11,*pe12,*pe13;
    register DOUBLE*pse;
    DOUBLE *pe0=tee+NKP+j*KLP;
    register int k;
    pse=se+j*(CKL);
    CALC12( 0);//* resonate 1 // kjoff : -2,-4   (16+8)*kl
    CALC12( 2);//* resonate 2 // kjoff : -2,-4   (16+8)*kl
    CALC34( 4);//* resonate 3 // kjoff : -3,-6   (16+9)*kl
    CALC34( 6);//* resonate 4 // kjoff : -3,-6   (16+9)*kl
    CALC56( 8);//* resonate 5 // kjoff :  3, 5   (16+9)*kl
    CALC56(10);//* resonate 6 // kjoff :  3, 5   (16+9)*kl
  }
}
#endif
#undef PWKA
#undef PWKB
#undef INTERE12
#undef CALC12
#undef CALC34
#undef CALC56
void implschs_(DOUBLE *eec,DOUBLE *eet,int *piacb,int *piace) {
  implschs(eec,eet,*piacb,*piace) ;
}
//CALC pm (8+NJ)*NK+9 div 5 sqrt 1 exp 1
//CALC pm 150*NJ*NK+7 div 1 exp 1
//CALC pm (4+NJ)*NK+5+0.1*NK*NJ div 2 sqrt 1 tanh 1
//CAL  pm (6+11*NK)*NJ+7 div 0 sqrt 1
//CALA pm NWPC*(6*NJ+12*NK+163*NK*NJ+28+0.1*NK*NJ div 8 sqrt 3 exp 2 tanh 1)
extern "C" void cwriteee(const char*fn,double*ee);
void cwriteee(const char*fn,double*ee){
  char fnt[256];
  int i,j,fid,ierr;
#ifdef ARM_VEC
#ifndef MTHREAD
  sprintf(fnt,"v/%s",fn);
#endif
#else
  sprintf(fnt,"s/%s",fn);
#endif
  FILE*fo=fopen(fnt,"wt");
  for(j=0;j<jnthet;j++){
    for(i=0;i<kl;i++){
      fprintf(fo,"%13.4g",ee[i+j*kl]);
      if(i%16==15)fprintf(fo,"\n");
    }
  }
  fclose(fo);
}
void implschs(DOUBLE *eec,DOUBLE *eet,int iacb,int iace) {
  int iac;
  int k,j,kj;
  DOUBLE se[mkj+kl];
  DOUBLE sds,awk,enh;
  ImplschPar *lpp=&gd->pisp;
  DOUBLE*ee=eet+mkj*iacb;
  DOUBLE*ec=eec+mkj*iacb;
#ifdef ARM_VEC
  svbool_t pm=svptrue_b64();
#endif
  for(iac=iacb;iac<=iace;iac++,ee+=mkj,ec+=mkj){
    DOUBLE enh;
    if(gd->nsp[iac]!=1)continue;
    ImplschP*pi=&lpp->pip[iac];
    windvs *wvs=gd->pwvs+iac;

    //CALC pm (8+NJ)*NK+ 9 div 5 sqrt 1 exp 1
    mean2(lpp,ee,pi,&sds,&awk);
    //!*************************************
    //! snonlin(e)
    //!*************************************
    //CALC pm NJ*NK*150+7 div 1 exp 1
    DOUBLE xx=0.75*pi->depth*awk;
    if (xx<0.5)xx=0.5;
    enh=1.+(5.5/xx)*(1.-0.833*xx)*exp(-1.25*xx); //MATH exp
#ifdef USEXNL
    xnlcal(ee,kl,se,(double*)&iac);
#else 
    winter(lpp,ee,se);
#endif
    //if(iac==iacb) cwriteee("se.txt",se);
    DOUBLE *cwks17=lpp->cwks17;
#ifdef ARM_VEC
    VDOUBLE venh=VDUP(enh);
      for(j=jnthet-1;j>=0;j--){
        double *pse=se+j*CKL;
        for(k=0;k<CKL;k+=NKP){// 9 click 3*4 op
          VVD(pse+k)*=venh*VVD(cwks17+k);
        }
      }
#else
      for(j=jnthet-1;j>=0;j--){
        double *pse=se+j*CKL;
        for(k=0;k<CKL;k++){// 9 click 3*4 op
          pse[k]*=enh*cwks17[k];
        }
      }
#endif
    //CAL pm 7+(6+11*NK)*NJ div 0 sqrt 1
    {   
      register DOUBLE deltts;
      DOUBLE scd=sqrt((0.80+0.065*wvs->wv)*0.001);
      DOUBLE bett =lpp->beta10*(1.+lpp->beta11*wvs->wi);
#ifndef NODELTEE
      DOUBLE wstarm=wvs->wv*scd;
#endif  
      //Wind Input

      deltts=lpp->deltts;
      for(j=0;j<jnthet;j++) {
        //!sinput(e)
        DOUBLE *pe =ee+j*kl;
        DOUBLE *pc =ec+j*kl;
        DOUBLE *pse=se+j*kl;
#ifdef ARM_VEC
        VDOUBLE betta;
        {
          DOUBLE wl=wvs->wx*lpp->cosths[j]+wvs->wy*lpp->sinths[j];
          betta=VDUP(bett*28.*wl*scd);
        }
        for(k=0;k<kl;k+=NKP){ // vector
          //! sinput(e)
          VDOUBLE wk=VVD(lpp->wk+k);
          VDOUBLE beta=betta*wk-bett*VVD(pi->ws+k);
          //if(beta<0)beta=0;
          //!*************************************
          //!      source terms end
          //!beta wind input
          //!ssds:sdissip
          //!ISSBO:sbottom
          //!cg*cosths[j]*Rsd_tanLat spherical coords
          //! ZZP cg=pi->CCG[k]  ! Big Circle
          VDOUBLE dset=beta-sds*wk+VVD(pi->ssbo+k); //! +pi->CCG[k]*cosths[j]*tRsd_tanLat
          {
            VDOUBLE eev=VVD(pe+k);
            VDOUBLE set=VVD(pse+k)+dset*eev;
            svbool_t pq;
#ifndef NODELTEE
            {
              VDOUBLE deltee=wstarm*VVD(lpp->grolim+k);
              VDOUBLE delteem=-deltee;
              //if(set<-deltee) set=-deltee;
              //else if(set>deltee) set=deltee;
              pq=svcmpgt(pm,set,deltee);
              set=svsel(pq,deltee,set);
              pq=svcmplt(pm,set,delteem);
              set=svsel(pq,delteem,set);
            }
#endif          
            eev=eev+set*deltts;
            //if(eev<0)eev=0;
            VDOUBLE vzero=VDUP(0.);
            pq=svcmplt(pm,eev,vzero);
            eev=svsel(pq,vzero,eev);
            VVD(pc+k)=eev;
          }
        }
#else
        DOUBLE betta;
        {
          DOUBLE wl=wvs->wx*lpp->cosths[j]+wvs->wy*lpp->sinths[j];
          betta=bett*28.*wl*scd;
        }
        for(k=0;k<kl;k++){ // vector
          //! sinput(e)

          DOUBLE beta=betta*lpp->wk[k]-bett*pi->ws[k];
          //if(beta<0)beta=0;
          //!*************************************
          //!      source terms end
          //!beta wind input
          //!ssds:sdissip
          //!ISSBO:sbottom
          //!cg*cosths[j]*Rsd_tanLat spherical coords
          //! ZZP cg=pi->CCG[k]  ! Big Circle
          DOUBLE dset=beta-sds*lpp->wk[k]+pi->ssbo[k]; //! +pi->CCG[k]*cosths[j]*tRsd_tanLat
          {
            DOUBLE eev=pe[k];
            DOUBLE set=pse[k]+dset*eev;
#ifndef NODELTEE
            {
              DOUBLE deltee=wstarm*lpp->grolim[k];
              if(set<-deltee){
                set=-deltee;
              }else if(set>deltee){
                set=deltee;
              }
            }
#endif          
            eev=eev+set*deltts;
            if(eev<0)eev=0;
            pc[k]=eev;
          }
        }
#endif
      }
    }
    //CALC pm (4+NJ)*NK+0.1*NK*NJ+ 5 div 2 sqrt 1 tanh 1
    mean3(lpp,ec,pi);
  }
}
//CALC pm (8+NJ)*NK+ 9 div 5 sqrt 1 exp 1
void mean2(ImplschPar *lpp,DOUBLE*ee,ImplschP*pi,DOUBLE *psds,DOUBLE*pawk) {
  DOUBLE ae,asi,ark,ekspm,awk;
  int k,j,kj;
#ifdef ARM_VEC
  VDOUBLE vae,vasi,vark,vawk,vpes;
  vae=vasi=vawk=vark=vpes=VDUP(0.);
  svbool_t pm=svptrue_b64();
  for(k=0;k<kl;k+=8){
    vpes=VDUP(0.);
    for(j=0,kj=k;j<jnthet;kj+=kl,j++){//vector
      VDOUBLE vpe=VLD(ee+kj);
      vpes+=vpe;
    }
    vae +=vpes*VLD(lpp->dwk+k);
    vawk+=vpes*VLD(lpp->wkdk+k);
    vark+=vpes*VLD(lpp->wkibdk+k);
    vasi+=vpes*VLD(pi->iwsdk+k);
  }
  ae =VSUM(vae ); awk=VSUM(vawk);
  ark=VSUM(vark); asi=VSUM(vasi);
#else
  ae=asi=awk=ark=0.;
  for(k=0;k<kl;k++){
    DOUBLE pes=0.;
    for(j=0,kj=k;j<jnthet;kj+=kl,j++){//vector
      pes+=ee[kj];
    }
    ae +=pes*lpp->dwk[k];
    awk+=pes*lpp->wkdk[k];
    ark+=pes*lpp->wkibdk[k];
    asi+=pes*pi->iwsdk[k];
  }
#endif
  {
    if(ae>1.e-100){
      awk=awk/ae;
      ark=(ae/ark);
      ark=ark*ark;
      asi=ae/(asi*ark);
      //!     sds=2.36e-5*asi*ark**3*ae**2/alpm2
      //!         2.36e-5/alpm2=2.587605807
      ekspm=ae*ark*ark/0.0030162;
      *psds=lpp->ads*lpp->brkd1*asi*sqrt(ekspm)*exp(-lpp->brkd2*0.64/ekspm); //MATH sqrt exp
    }else{
      *psds=0;
    }
    *pawk=awk;
  }
}
//CALC pm (4+NJ)*NK+0.1*NK*NJ+ 5 div 2 sqrt 1 tanh 1
void mean3(ImplschPar *lpp,DOUBLE*ee,ImplschP*pi) {
  int    k,j,kj;
  register DOUBLE ae,awk;
#ifdef ARM_VEC
  VDOUBLE vae,vawk,vev;
  vae=vawk=VDUP(0.);
  svbool_t pm=svptrue_b64();
  for(k=0;k<kl;k+=8){
    vev=VDUP(0.);
    for(kj=k,j=0;j<jnthet ;kj+=kl,j++){//vector
      vev+=VLD(ee+kj);
    }
    vae +=vev*VLD(lpp->dwk+k);
    vawk+=vev*VLD(lpp->wkdk+k);
  }
  ae =VSUM(vae ); awk=VSUM(vawk); 
#else
  ae=0.;    awk=0.;
  for(k=0;k<kl;k++){
    DOUBLE vev=0.;
    for(kj=k,j=0;j<jnthet ;kj+=kl,j++){//vector
      vev+=ee[kj];
    }
    ae +=vev*lpp->dwk[k];
    awk+=vev*lpp->wkdk[k];
  }
#endif
  {
    if(ae>1e-30){
      DOUBLE hs,hb,hbb;
      awk=awk/ae;
      hs=4.*sqrt(ae); //MATH sqrt 
      hb=lpp->zpi/awk*0.142*tanh(pi->depth*awk); //!hb=zpi/awk*0.12*tanh(dep(ia,ic)*awk)/1.6726 //MATH tanh
      hbb =pi->depth*(0.78125/1.6726)     ; //!hbb=0.78125*dep(ia,ic)/1.5864792
      if(hb>hbb) hb=hbb;
      if(hs>hb) {
        DOUBLE vt;
        DOUBLE chbh;
        int kj;
        vt=(hb/hs);chbh=vt;
        chbh=chbh*chbh;
#ifdef ARM_VEC
        VDOUBLE vchbh=VDUP(chbh);
        for(kj=0;kj<mkj;kj+=8){//vector          
          VVD(ee+kj)*=vchbh;
        }
#else
        for(kj=0;kj<mkj;kj++){//vector          
          ee[kj]*=chbh;
        }
#endif
      }
    }
  }
}
void c_checkee_(DOUBLE*ee,int *iac,const char*fn,int *ln){
  Checkee(ee,*iac,fn,*ln);
}
extern "C" void wav_abort();
void Checkee(DOUBLE*ee,int iac,const char*fn,int ln){
  DOUBLE ae=0,be,me=0;
  int k,j;
  for(j=0;j<jnthet;j++){
    DOUBLE*pe=ee+j*kl;
    for(  k=0;k<kl ;k++){//vector
      be=fabs(pe[k]);
      if(be>me)me=be;
      ae =ae +be;
    }
  }
  if(isnan(ae)){
    printf("ee is nan %s %d %d \n",fn,ln,iac);
    wav_abort();
  } 
  if(isinf(ae)){
    printf("ee is inf %s %d %d \n",fn,ln,iac);
    wav_abort();
  } 
  if(fabs(me)>1e80){
    printf("ee is big %s %d %d %e\n",fn,ln,iac,me);
    wav_abort();
  } 
}
__BEGIN_DECLS
void prtee(int id,DOUBLE*ee);
void prtee_(int *id,DOUBLE*ee);
__END_DECLS
void prtee(int id,DOUBLE*ee){
  ee+=mkj;
  printf("AA %d:",id);
  for(int j=0;j<1;j++) {
    for(int k=0;k<10;k++){
      printf(" %9.2e ",ee[j*kl+k]);
    }
    printf("\n");
  }
}

void prtee_(int *id,DOUBLE*ee){
  prtee(*id,ee);
}
