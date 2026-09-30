#include "svars.h"
extern prgstate pstates;
//#define NODELTEE
#define LOGID 0
#define DBGINFN printf("%s %d %d %d:",__FILE__,__LINE__,itt,iat);
void Checkee(double*ee,int iac,const char*fn,int ln);
#define CHKEE(ee,iat) //Checkee(ee,iat,__FILE__,__LINE__)
#define CNAN(v,i1,i2,i3,i4) if(isnan(v)){printf("NAN %s %d %d %d %d:%s %d\n",#v ,i1,i2,i3,i4,__FILE__,__LINE__);}
void mean2(ImplschPar *lpp,double*ee,ImplschP*pi,double *sds,double*awk) ;
void mean3(ImplschPar *lpp,double*ee,ImplschP*pi) ;
extern "C" int PT_openw(const char*fn);
extern "C" int PT_write(int ifile,int ind,double *ee,double*se);
extern "C" int PT_close_(int ifile);
int iact=0;
static int ofid=-1,eeind=0;
static void ofai4sinit(){
  if(ofid<0){
    ofid=PT_openw("winter.dat");
    eeind=0;
  }
}
static void ofai4s(double *ee,double*se){
  ofai4sinit();
  PT_write(ofid,eeind,ee,se);
  eeind++;
}
void InitCImplsch(){
}
#define register 
#define NKP 8
#define KLP (CKL+NKP)
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
      register double e0,e1,e2;e0=pe0[k]; INTERE12(lmdpd,lmdmd);  \
      pse[k]=pse[k]+LIV_m2*e0*(e0*e1+e0*e2+e1*e2);               \
    }
//CAL pm 9+16
#define CALC34(I0)                                               \
    PWK(I0);                                                     \
    for(k=0;k<CKL;k++){                                         \
      register double e0,e1,e2;e0=lmdpd*pe0[k]; INTERE12(LIV_1,lmdmd);  \
      pse[k]=pse[k]+lmdpd20*e1*(e1*e0+e1*e2+e0*e2);   \
    }
//CAL pm 9+16
#define CALC56(I0)                                               \
    PWK(I0);                                                     \
    for(k=0;k<CKL;k++){                                         \
      register double e0,e1,e2;e0=lmdmd*pe0[k]; INTERE12(lmdpd,LIV_1);  \
      pse[k]=pse[k]+lmdmd20*e2*(e2*e0+e2*e1+e0*e1);   \
    }
//CALC pm NJ*NK*150
void winter(ImplschPar *lpp,double *ee,double *se){
  register double lmdpd,lmdmd,lmdpd20,lmdmd20,LIV_m2,LIV_1;
  int j;
  double tee[KLP*CJNTHET+NKP]; 
  if(1){
    double *ped=tee,*pes=ee;
    int k,j;
    for(k=0;k<NKP;k++) *ped++=0;
    for(j=0;j<CJNTHET;j++){
      for(k=0;k<CKL;k++) *ped++=*pes++;
      for(k=0;k<NKP;k++) *ped++=0;
    }
  }
  LIV_m2=-2.;LIV_1=1.;
  lmdpd=lpp->lmdpd;lmdmd=lpp->lmdmd;lmdpd20=lpp->lmdpd20;lmdmd20=lpp->lmdmd20;
  memset(se,0,sizeof(double)*mkj);
  for(j=jnthet-1;j>=0;j--){
    register double wka00,wka01,wka02,wka03;
    register double wka10,wka11,wka12,wka13;
    register double*pe00,*pe01,*pe02,*pe03;
    register double*pe10,*pe11,*pe12,*pe13;
    register double*pse;
    double *pe0=tee+NKP+j*KLP;
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

#undef PWKA
#undef PWKB
#undef INTERE12
#undef CALC12
#undef CALC34
#undef CALC56
__BEGIN_DECLS
void prtee(int id,double*ee,double*eeo);
void prtee_(int *id,double*ee,double*eeo);
__END_DECLS
void prtee(int id,double*ee,double*eeo){
  int j=1;
  double ae=0.,me=0.,hs,hb,awk=0;
  int jt=1; 
  double mdes[12],mde,de;
  int mks[12];
  ImplschPar *lpp=&pisp;
  ImplschP*pi=&lpp->pip[1];
  printf("BB :");
  for(int j=0;j<jnthet ;j++){//vector
    double vev=0.,vew=0;
    for(int k=0;k<kl;k++){
      vev+=ee[j*kl+k]*lpp->dwk[k];
      vew+=ee[j*kl+k]*lpp->wkdk[k];
    }
    if(vev>me){me=vev;jt=j;}
    if(vev>1e-1) printf("%2d: %10.2e %10.2e ",j,vev,vew);
    ae +=vev;
    awk+=vew;
  }
  printf("\n");
  printf("CC :");
  for(int k=0;k<kl;k++){
    double vev=0.,vew=0;
    for(int j=0;j<jnthet ;j++){//vector
      vev+=ee[j*kl+k];
    }
    if(vev*lpp->dwk[k]>1e-1) printf("%2d: %10.2e %10.2e ",k,vev*lpp->dwk[k],vev*lpp->wkdk[k]);
  }
  printf("\n");
  awk=awk/ae;
  hs=4.*sqrt(ae);
  hb=lpp->zpi/awk*0.142*tanh(pi->depth*awk); //!hb=zpi/awk*0.12*tanh(dep(ia,ic)*awk)/1.6726
  printf("AA %2d %4ld %5.2f %5.2f %10.4e %10.3e:\n",id,pstates.ntimestep,hs,hb,awk,ae);
  for(j=0;j<12;j++){
    printf("%2d:",j);
    for(int k=0;k<32;k++){
      printf(" %9.2e",ee[j*kl+k]);
    }
    printf("\n");
  }
  printf("  :");
  for(int k=0;k<32;k++){
    printf(" %9d",k);
  }
  printf("\n");
  for(j=0;j<12;j++){
    printf("%2d:",j);
    mde=0;
    for(int k=0;k<32;k++){
      de=ee[j*kl+k]-eeo[j*kl+k];
      if(abs(de)>mde){mde=abs(de);mks[j]=k;}
      printf(" %9.2e",de);
    }
    printf("\n");
    mdes[j]=mde;
  }
  printf("MMM :");
  for(j=0;j<12;j++){
    printf("%2d:%2d %9.2e ",j,mks[j],mdes[j]);
  }
  printf("\n");
  if(pstates.ntimestep>4)exit(0);
}

void prtee_(int *id,double*ee,double*eeo){
  prtee(*id,ee,eeo);
}

//CALC pm (8+NJ)*NK+9 div 5 sqrt 1 exp 1
//CALC pm 150*NJ*NK+7 div 1 exp 1
//CALC pm (4+NJ)*NK+5+0.1*NK*NJ div 2 sqrt 1 tanh 1
//CAL  pm (6+11*NK)*NJ+7 div 0 sqrt 1
//CALA pm NWPC*(6*NJ+12*NK+163*NK*NJ+28+0.1*NK*NJ div 8 sqrt 3 exp 2 tanh 1)
void implschs(double *eec,double *eet,int iacb,int iace) {
  int iac;
  int k,j,kj;
  double se[mkj];
  double sds,awk;
  ImplschPar *lpp=&pisp;
  //iacb=1; iace=gppar.nwpc;
  double*ee=eet+mkj*iacb;
  double*ec=eec+mkj*iacb;
  double eo[mkj];
  pstates.ntimestep++;
  for(iac=iacb;iac<=iace;iac++,ee+=mkj,ec+=mkj){
    double enh;
    if(gppar.nsp[iac]!=1)continue;
    ImplschP*pi=&lpp->pip[iac];
    iact=iac;

    //CALC pm (8+NJ)*NK+ 9 div 5 sqrt 1 exp 1
    mean2(lpp,ee,pi,&sds,&awk);
    //!*************************************
    //! snonlin(e)
    //!*************************************
    //CALC pm NJ*NK*150+7 div 1 exp 1
    winter(lpp,ee,se);
    {
      double xx=0.75*pi->depth*awk;
      double enh;
      if (xx<0.5)xx=0.5;
      enh=1.+(5.5/xx)*(1.-0.833*xx)*exp(-1.25*xx);
      for(j=jnthet-1;j>=0;j--){
        double *pse=se+j*(CKL);
        double *cwks17=lpp->cwks17;
        for(k=0;k<CKL;k++){// 9 click 3*4 op
          pse[k]=pse[k]*enh*cwks17[k];
        }
      }
    }
    ofai4s(ee,se);
    //CAL pm 7+(6+11*NK)*NJ div 0 sqrt 1
    {   
      register double deltts;
      double scd=sqrt((0.80+0.065*pi->wvs.wv)*0.001);
      double bett =lpp->beta10*(1.+lpp->beta11*pi->wvs.wi);
#ifndef NODELTEE
      double wstarm=pi->wvs.wv*scd;
#endif  
      //Wind Input

      deltts=lpp->deltts;
      for(j=0;j<jnthet;j++) {
        //!sinput(e)
        double *pe =ee+j*kl;
        double *pc =ec+j*kl;
        double *pse=se+j*kl;
        double betta;
        {
          double wl=pi->wvs.wx*lpp->cosths[j]+pi->wvs.wy*lpp->sinths[j];
          betta=bett*28.*wl*scd;
        }
        for(k=0;k<kl;k++){ // vector
          //! sinput(e)

          double beta=betta*lpp->wk[k]-bett*pi->ws[k];
          //if(beta<0)beta=0;
          //!*************************************
          //!      source terms end
          //!beta wind input
          //!ssds:sdissip
          //!ISSBO:sbottom
          //!cg*cosths[j]*Rsd_tanLat spherical coords
          //! ZZP cg=pi->CCG[k]  ! Big Circle
          double dset=beta-sds*lpp->wk[k]+pi->ssbo[k]; //! +pi->CCG[k]*cosths[j]*tRsd_tanLat
          {
            double eev=pe[k];
            double set=pse[k]+dset*eev;
#ifndef NODELTEE
            {
              double deltee=wstarm*lpp->grolim[k];
              if(set<-deltee){
                set=-deltee;
              }else if(set>deltee){
                set=deltee;
              }
            }
#endif          
            if(iac==1&& pstates.ntimestep==4){
              printf("DD:%2d %2d %10.2e %10.2e %10.2e %10.2e\n",j,k,eev,pse[k]*deltts,set*deltts,eev+set*deltts);
            }
            eev=eev+set*deltts;
            if(eev<0)eev=0;
            pc[k]=eev;
          }
        }
      }
    }
    //CALC pm (4+NJ)*NK+0.1*NK*NJ+ 5 div 2 sqrt 1 tanh 1
    if(iac==1) {
      prtee(3,ec,ee);
      printf("WV %f ",pi->wvs.wv);
      memcpy(eo,ec,mkj*sizeof(double));
    }
    mean3(lpp,ec,pi);
    if(iac==1) prtee(5,ec,eo);
  }
}

//CALC pm (8+NJ)*NK+ 9 div 5 sqrt 1 exp 1
void mean2(ImplschPar *lpp,double*ee,ImplschP*pi,double *psds,double*pawk) {
  double ae,asi,ark,ekspm,awk;
  int k,j,kj;
  ae=asi=awk=ark=0.;
  for(k=0;k<kl;k++){
    double pes=0.;
    for(j=0,kj=k;j<jnthet;kj+=kl,j++){//vector
      pes+=ee[kj];
    }
    ae +=pes*lpp->dwk[k];
    awk+=pes*lpp->wkdk[k];
    ark+=pes*lpp->wkibdk[k];
    asi+=pes*pi->iwsdk[k];
  }
  {
    if(ae>1.e-100){
      awk=awk/ae;
      asi=ae/asi;
      ark=(ae/ark);
      ark=ark*ark;
      //!     sds=2.36e-5*asi*ark**3*ae**2/alpm2
      //!         2.36e-5/alpm2=2.587605807
      ekspm=ae*ark*ark/0.0030162;
      *psds=lpp->ads*lpp->brkd1*asi/ark*sqrt(ekspm)*exp(-lpp->brkd2*0.64/ekspm);
    }else{
      *psds=0;
    }
    *pawk=awk;
  }
}
//CALC pm (4+NJ)*NK+0.1*NK*NJ+ 5 div 2 sqrt 1 tanh 1
void mean3(ImplschPar *lpp,double*ee,ImplschP*pi) {
  int    k,j,kj;
  register double ae,awk;
  ae=0.;    awk=0.;
  for(k=0;k<kl;k++){
    double vev=0.;
    for(kj=k,j=0;j<jnthet ;kj+=kl,j++){//vector
      vev+=ee[kj];
    }
    ae +=vev*lpp->dwk[k];
    awk+=vev*lpp->wkdk[k];
  }
  {
    if(ae>1e-30){
      double hs,hb,hbb;
      awk=awk/ae;
      hs=4.*sqrt(ae);
      hb=lpp->zpi/awk*0.142*tanh(pi->depth*awk); //!hb=zpi/awk*0.12*tanh(dep(ia,ic)*awk)/1.6726
      hbb =pi->depth*(0.78125/1.6726)     ; //!hbb=0.78125*dep(ia,ic)/1.5864792
      if(hb>hbb) hb=hbb;
      if(iact==1) printf(" :  %f %f %f %f %f ",ae,awk,hs,hb,hbb);
      if(hs>hb) {
        double vt;
        double chbh;
        int kj;
        vt=(hb/hs);chbh=vt;
        chbh=chbh*chbh;
        if(iact==1) printf(":   %f %f",vt,chbh);

        for(kj=0;kj<mkj;kj++){//vector          
          ee[kj]*=chbh;
        }
      }
      if(iact==1) printf("\n");
    }
  }
}
void c_checkee_(double*ee,int *iac,const char*fn,int *ln){
  Checkee(ee,*iac,fn,*ln);
}
extern "C" void wav_abort();
void Checkee(double*ee,int iac,const char*fn,int ln){
  double ae=0,be,me=0;
  int k,j;
  for(j=0;j<jnthet;j++){
    double*pe=ee+j*kl;
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
void implschs_(double *eec,double *eet,int *piacb,int *piace) {
  implschs(eec,eet,*piacb,*piace) ;
}
