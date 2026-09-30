//#define DTIMES
#include "svars.h"
#include "ctools.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <unistd.h>
double dclktime(void);
ImplschPar pisp;
prgstate pstates;
GData pd;
#ifndef MTHREAD
TData tdv;
__THREAD GData *gd=&pd;
__THREAD TData *td=&tdv;
#else
__thread struct TData*td=NULL;
#endif
//save all vars to pd from gppar
//same struct with gd
void InitGPar(){
  dclktime();
  pd.mpi_id=gppar.mpi_id;
  pd.mpi_comm_wav=gppar.mpi_comm_wav;
  pd.mpi_npe=gppar.mpi_npe;
  pd.pisp.mpi_id=gppar.mpi_id;
  pd.nwps  =gppar.nwps  ;
  pd.nwpc  =gppar.nwpc  ;
  pd.nwpa  =gppar.nwpa  ;
  pd.dep   =gppar.dep   ;
  //pd.nwpw  =0;
  pd.LXB   =gppar.LXB  ;
  pd.LYB   =gppar.LYB  ;
  pd.LXN   =gppar.LXN  ;
  pd.LYN   =gppar.LYN  ;
  pd.nsp   =gppar.nsp;
  pd.mTimeStep= gppar.mTimeStep;
  pd.ieind =gppar.ieind;
  pd.iepos =gppar.iepos;
  pd.wav_et=gppar.wav_et;//nwpc
  pd.wav_ec=gppar.wav_ec;//nwpc
  pd.psis  =gppar.psis  ;
  pd.pgis  =gppar.pgis  ;
  pd.ipos12=gppar.ipos12;
  pd.ipos8 =gppar.ipos8;
  pd.constwind=gppar.constwind;
 	pd.constwind.wv=pd.constwind.wx*pd.constwind.wx+pd.constwind.wy*pd.constwind.wy;
  pd.constwind.wi=0;
  pd.pwvsg  =gppar.pwvs;
#ifdef MMTHREAD
  pd.pwvs  =gppar.pwvs;
#elif defined(MTHREAD)
  pd.pwvs  =HNEWN(windvs ,pd.nwpc+1);
#else
  pd.pwvs  =gppar.pwvs;
#endif
 	int acwv=pd.constwind.wv*100;
  if(acwv>=10){
  	pd.uconstwind=2;
  }else if(acwv>=1){
  	pd.uconstwind=1;  	
  }else {
  	pd.uconstwind=0;  	
  }
  gd=&pd; //for MTHREAD InitCPropgats
}
void c_initimplsch_(ImplschPar *pisp_,ImplschP*pip){
  int isize;
  pd.pisp.pip=pip;
  isize=(long)&(((ImplschPar *)0)->sds);
  memcpy((void*)&pd.pisp,(void*)pisp_,isize );//for align so copy it
  //pd.pisp.pis=pd.pis;
  //pd.pisp.pis=pd.pis;
  //pd.pisp.ipos8=pd.ipos8; // nwpc
  pd.pisp.Dm2=-2;pd.pisp.D1=1.;pd.pisp.enh=1.;
  pd.pisp.zpi=3.1415926535897932384626433832795*2.;
  pd.pisp.D1m30=1e-30;
  pd.pisp.F44_24778=44.24778;
  pd.pisp.F16_0=16.;
  pd.pisp.F0_467=0.467;//mean3			
  pd.pisp.F18_20832=18.20832;
  pd.pisp.F0_001930368=0.001930368;
  pd.pisp.F0_5=0.5;
  pd.pisp.F1_0       = 1        ;
  pd.pisp.F0_75      = 0.75     ;
  pd.pisp.F5_5       = 5.5      ;
  pd.pisp.F0_833     = 0.833    ;
  pd.pisp.Fm1_25     =-1.25     ;
  pd.pisp.F28_       = 28       ;
  pd.pisp.F0_80      = 0.80     ;
  pd.pisp.F0_065     = 0.065    ;
  pd.pisp.F0_001     = 0.001    ;
  pd.pisp.tztp       = 1.2      ;
  pd.pisp.tzzpi      = pd.pisp.zpi*1.099314;//zpi*1.099314;
}
void getwind(){ // ZZZ should opt to sdma 
  static int first=1;
#ifdef MMTHREAD
	if(pd.uconstwind==0){
		windvs * wvs=gd->pwvs;
		windvs * wvsg=pd.pwvsg;
    for(int i=1;i<=gd->nwpc;i++){
      wvs[i]=wvsg[gd->l2gind[i]];
    }
	}else{
    if(!first)return;
    first=0;			
		windvs * wvs=gd->pwvs;
    if(pd.uconstwind==2){
      for(int i=1;i<=gd->nwpc;i++){
        wvs[i]  =pd.constwind;
      }
    }else{
      windvs * wvsg=pd.pwvsg;
      for(int i=1;i<=gd->nwpc;i++){
        wvs[i]=wvsg[gd->l2gind[i]];
      }				
    }
	}
#elif defined( MTHREAD)
	if(pd.uconstwind==0){
    memcpy(pd.pwvs+1,pd.pwvsg+1,sizeof(windvs)*(pd.nwpc));
  }else{
    if(!first)return;
    first=0;
    windvs *wvs=pd.pwvs;
    if(pd.uconstwind==2){
      for(int i=1;i<=pd.nwpc;i++){
        wvs[i]  =pd.constwind;
      }
      return;
    }
    memcpy(pd.pwvs+1,pd.pwvsg+1,sizeof(windvs)*(pd.nwpc));
  }
#endif
}
void setflag_(int*nostep,int *hist_eot,int *stop_now,int *rest_eot){
  pd.hist_eot=*hist_eot;
  pd.stop_now=*stop_now;
  pd.rest_eot=*rest_eot;
  pd.nostep=*nostep;
#ifdef MMTHREAD
  for(int i=0;i<md.ngrp;i++){
    struct GData*pgd=md.grps[i]->gd;
    pgd->nostep=*nostep;
    pgd->hist_eot=*hist_eot;
    pgd->stop_now=*stop_now;
    pgd->rest_eot=*rest_eot;
  }
#endif
}
#ifdef MTHREAD
//NGrpPProc,NThPGrp;
/*
 *iepos[ 0:nwpa ]
 *ieind[0,LXN*LYN-1],ind=ix*LXN
 *nsp [ 0:nwpa ],flags,0:land,1:water,2:open boundary
 * */
void f_initGPar(){
  //init in ouput_cal
  //double *pvkd,*pebdep;
  //float *ape , *tpf , *aet , *h1_3, *bv ;
}
#ifdef MMTHREAD
void fillsendcomm();
Part *spart;
struct Iepos{
  int ix,iy;
};
void psegs(int tag,int id,sendcomm* sc,int itemsize,int m){
  for (int j=0;j<=m;j++){
    struct sendsegs *cur_segs=&sc->ssegs[j];
    if (cur_segs->dst_id==-1)continue;//break;
		printf("VV %4d %d %2.2d %2d %2d:%2d %p\n",tag,pd.mpi_id,id,j ,cur_segs->dst_id,cur_segs->nsegs,cur_segs->pflag);
	}
}
static int NCSEP[16*2]={
  1,1, 2,1, 3,1, 2,2, //1   2  3  4
  5,1, 3,2, 7,1, 4,2, //5   6  7  8
  3,3, 5,2,11,1, 4,3, //9  10 11 12
 13,1, 7,2, 5,3, 4,4, //13 14 15 16
};
void InitThreadsGPar(){
  int NCX=4,NCY=4;
  if(cpuinf.NGrpPProc<=16&&cpuinf.NGrpPProc>0){
    NCX=NCSEP[cpuinf.NGrpPProc*2-2];
    NCY=NCSEP[cpuinf.NGrpPProc*2-1];
  }else{
    NCX=cpuinf.NGrpPProc;NCY=1;
  }
  NEWZN(spart,cpuinf.NGrpPProc +1);
  GPart(spart,pd.ieind,NCX,NCY,pd.nwps,pd.nwpc,pd.nwpa,pd.LXN,pd.LYN,-(int)(1+100*1/38));
  /* partion here,alloc mems
   * ZZZ fill  group info  */
  pd.s_segs=spart[0].send_segs;
#if 0
  char fn[256];
  sprintf("TTA_%2.2d.txt",mpi_id);
  FILE*fo=fopen(fn,"wt");
  for(int ii=0;ii<pd.nwpa;ii++){
    fprintf(fo,"%d %d %d %d\n",ii,pd.g2l[ii].ipart,pd.g2l[ii].iac);
  }
  fclose(fo);
#endif
  int gbase=(mpi_id%cpuinf.NProcPNode)*cpuinf.NGrpPProc;
  for(int p=0;p<cpuinf.NGrpPProc;p++){
    GData *pgd=(GData*)md.grps[p]->gd;
    bindcpu((p+gbase)*cpuinf.NCorePClu );//FFA
    Part *pt=spart+p+1;
    //printf("GGGGGGEEE %d %d %5d \n",mpi_id,p,pt->nwpa);
    pgd->block_id=pt->block_id;
    pgd->NPART   =pt->NPART;
    pgd->LXB =pt->glb_ixb;
    pgd->LYB =pt->glb_jyb;
    pgd->LXN     =pt->GRXN;
    pgd->LYN     =pt->GRYN;
    pgd->recti   =pt->recti;
    pgd->recto   =pt->recto;
    pgd->nwpt  =pt->nwpt;   //capital
    pgd->nwps  =pt->nwps;   //neighbor
    pgd->nwpc  =pt->nwpc;   //inside
    pgd->nwpo  =pt->nwpo;   //outside
    pgd->nwpa  =pt->nwpa;   //Z
    pgd->s_segs=pt->send_segs;
    size_t size=0;

#define VAVARS \
    VVAR(wav_ec,double,(size_t)(pgd->nwpa+1)*mkj);\
    VVAR(wav_et,double,(size_t)(pgd->nwpa+1)*mkj);\
    VVAR(l2gind,int   ,(size_t)(pgd->nwpa+1)    );\
    VVAR(gg2l  ,int   ,(size_t)(pd  .nwpa+1)    );\
    VVAR(nsp   ,char  ,(size_t)(pgd->nwpc+1)    );\
    VVAR(pwvs  ,windvs,(size_t)(pgd->nwpc+1)    );\
    VVAR(dep   ,float ,(size_t)(pgd->nwpc+1)    );

#define  VVAR(vn,T,N) size+=sizeof(T)*(size_t)(N)+256
    VAVARS;
#undef VVAR
    char*buf=pgd->fbuff=HNEWN(char,size);
    memset(buf,0,size);
#define  VVAR(vn,T,N) pgd->vn=(T*)buf;buf=(char*)( (((long)buf)+sizeof(T)*(size_t)(N)+255) & (-256) )
    VAVARS;
#undef VVAR
#undef VAVARS
    CPYN(pgd->l2gind,pt->grp2glb,pgd->nwpa+1);
    for(int i=1;i<=pgd->nwpc;i++){
      int l2g=pgd->l2gind[i];
      pgd->nsp[i]=pd.nsp[l2g];
      pgd->dep[i]=pd.dep[l2g];
    }
    for(int i=1;i<=pgd->nwpa;i++){
      pgd->gg2l[pt->grp2glb[i]]=i;
    }
    /* not used:
     * recv_segs  pemask
     * iepos  grp_ieind  */
  }
  pd.wav_et=gppar.wav_et;
  pd.wav_ec=gppar.wav_ec;
  fillsendcomm();
  bindcpu(0 );
  for(int p=0;p<cpuinf.NGrpPProc;p++){
    gp_freepart(spart+p);
  }
  free(spart);
}
void fillsendcomm(){
  sendcomm**scs;
  sendcomm *sc;
  scs=NEWN(sendcomm*,cpuinf.NGrpPProc+1);
  scs[0]=sc=&pd.s_segs;
  pd.wav_et=gppar.wav_et;
  pd.wav_ec=gppar.wav_ec;
  sc->src_ee=pd.wav_ec;
  //printf("FILL %d  %p %p\n",mpi_id,sc->src_ee,pd.wav_et); fflush(stdout);
  for(int p=0;p<cpuinf.NGrpPProc;p++){
    GData *pgd=(GData*)md.grps[p]->gd;
    scs[p+1]=sc=&pgd->s_segs;
    //*sc=spart[p].send_segs;
    sc->src_ee=pgd->wav_ec;
  }
  gp_filldst(scs,cpuinf.NGrpPProc+1);

	//psegs(125,0,&pd.s_segs,0,cpuinf.NGrpPProc);
  //for(int p=0;p<cpuinf.NGrpPProc;p++){
  //  GData *pgd=(GData*)md.grps[p]->gd;
	//	psegs(125,p+1,&pgd->s_segs,0,cpuinf.NGrpPProc);
	//}	
  free(scs);
}
#endif

int _getgdsize_(){return sizeof(GData);}
int _gettdsize_(){return sizeof(TData);}
int CheckStop(){
  return pd.stop_now;
}
void caliacbe(){
  int nt,n1,n2,nm,nn1,nwps,nwpc;
  float sc,ntsc;
  int id=ti->ind;
  nt=ti->Nthreads;  //ti->NThreads =gi->NThreads
  sc=0.3;ntsc=nt-sc;
  nwpc=gd->nwpc;// if NO MMTHREAD gd => pd
  nwps=gd->nwps;

  n2=nwps/ntsc+1;
  if(n2*(nt-1)>nwps)nwps=n2*(nt-1);
  td->iacb2=nwps-(nt-id  )*n2 ;
  td->iace2=nwps-(nt-id-1)*n2 ;
  if(td->iacb2<=0)td->iacb2=1;
  if(td->iace2<=td->iacb2)td->iace2=td->iacb2-1;

  n1=(nwpc-nwps)/ntsc+1;
  td->iacb1=nwpc-(nt-id  )*n1 ;
  td->iace1=nwpc-(nt-id-1)*n1 ;
  if(td->iacb1<=nwps)td->iacb1=nwps+1;
  if(td->iace1<td->iacb1)td->iace1=td->iacb1-1;
#if 0
  FILE*fo=ti->fo;
  if(!fo)fo=stdout;
  if(fo){
    fprintf(fo,"VGT %2.2d %2.2d:%4d %4d,%4d %4d %4d %4d:%4d %4d %4d: "
            "%4d %4d %4d %4d %4d %4d %f\n",
            ti->igrp,ti->ind,
            td->iace1-td->iacb1,td->iace2-td->iacb2,
            td->iacb1,td->iace1,td->iacb2,td->iace2,
            td->iace1-td->iacb1+1,td->iace2-td->iacb2+1,
            gd->nwps,nwps,nwpc,n1,n2,nt,id,ntsc
           );
    fflush(fo);
  }
#endif
}
#ifdef MMTHREAD
void c_sendboundary_(int *RFG){
  //printf("main set FLG %d\n",*RFG);
  gp_send_data(&pd.s_segs,sizeof(double)*mkj,cpuinf.NGrpPProc,*RFG);
  //gp_sdma_send_data(&pd.s_segs,sizeof(double)*mkj,cpuinf.NGrpPProc,*RFG);
}
void c_checkboundary_(int *RFG){
  //printf("main check FLG %d\n",*RFG);
  gp_check_recv(&pd.s_segs,cpuinf.NGrpPProc,*RFG);
}
#endif
/* M:main mpi thread
 * G:group main thread
 * S:group sub thread 
 * W:wait
 * S:set * */
#define DBG  //printf("AA %d %d %2.2d %2.2d %d %d %d\n",__LINE__,mpi_id,td->nstep,ti->igrp,ti->ind,RFB,RFG);
void prth(){
  char *hdm= "MT(s):id ig ind|   RWIND |  checkb |Exchange |setbbound|      NUL|      MWG|SetWind  |implsch 2| setspec2|    SSG  |      NUL|  accum 2|  SWG    |";       
  fprintf(ti->fo,"%s\n",hdm);
  char *hds= "ST(s):id ig ind|   prop 1|      SWG|implsch 1|setspec 1|  accum 1|      SWG|propagat2|implsch 2| setspec2|    SSG  |      NUL|  accum 2|  SWG    |";       
  fprintf(ti->fo,"%s\n",hds);
#ifdef MMTHREAD
  char *hdg= "GT(s):id ig ind|   prop 1|      GWM|implsch 1|setspec 1|  accum 1|  ckrecv1|propagat2|implsch 2| setspec2|    GWS  | senddata|  accum 2|  GWS    |";       
  fprintf(ti->fo,"%s\n",hdg);
#endif
}
#define DBGA if(ti->MainThread){printf("VV %4d %4d %4d %4d %d\n",pd.mpi_id,ti->igrp,ti->ind,__LINE__,RFB);fflush(stdout);}
void initgp(){
  int RFB=111;
  int iacb,iace;
  float anp=((float)pd.nwpc)/md.Nthreads;
  iacb=anp*ti->indp+1.3;
  iace=anp*ti->indp+anp+1.3-1;
  //printf("initgp %2.2d %2.2d %2.2d %3.3d %8.6d %8.6d %f %d\n",
  //       pd.mpi_id,ti->igrp,ti->ind,ti->indp,iacb,iace,anp,pd.nwpc);
  initsetspec_ (&iacb,&iace);
  initimplschs_(&iacb,&iace);
  initpropagat_(&iacb,&iace);
  initoutputcal_(&iacb,&iace,pd.dep);//fixme
}
void DBGHA(){
  printf("DBGHERE\n");
}

typedef struct xnlinf{
  int naq,nkq,nekq;
  int mkq,maq,klocus;
  int qf_dn,qf_kn;
  int  *quadnl;
  struct wws_typev*quads;
  double*q_k2,*q_ka,*q_sig,*q_sigr;
  double*q_kpow,*q_tail,*q_lamd;
}xnlinf;
void ckinf_(int*tag,int*ind,int*rfb){
  zlock_();
  printf("CKOUT %15.5f %5d %2.2d %2.2d %2.2d :%2.2d %2.2d \n",
         dclktime(),*tag,pd.mpi_id,ti->igrp,ti->ind,*ind,*rfb);
  fflush(stdout);
  zunlock_();
}
extern __THREAD xnlinf *ppxi;
int getcpuid();
int GetVInt(int volatile *volatile p);
#ifdef MMTHREAD
void _threadmain_(HTHREADINFO ti){
  int istep=0, RFB=1,RFG=1;
  int ind=ti->ind,indg=ti->indg,igrp=ti->igrp;
  double t0,tdd;
  if(!gd)gd=&pd; //for MTHREAD
  caliacbe();	  //not need wait	, ready here
  RFB=1;
  if(ti->MainThread){
    GSM(RFB); GWM(RFB); GSS(RFB);
  }else{
    SWG(RFB);
  }
  initgp();
  RFB=3;
  td->nstep=-1;
#define ILOG //if(istep<3){zlock_();printf("VV %2.2d %2.2d %2.2d %4d %4d\n",pd.mpi_id,ti->igrp,ti->ind,__LINE__,istep);fflush(stdout);zunlock_();}
  if(ti->MainThread){
    gd->nstep=-1;
    GWS(RFB); GSM(RFB);GWM(RFB);
    //usleep(5*1000000);
    g_initPropgats(); 
    g_initImplsch();
#ifdef USEXNL
    init_xnlg_();
#endif
    g_initSetspec();
    g_initOutput_cal();
    GSS(RFB); GSM(RFB);
  }else{ 
    SSG(RFB);SWG(RFB); 
  }
#ifdef USEXNL
  init_xnl_();
#else
  init_dia_();
#endif
  td->wav_et=gd->wav_et;
  td->wav_ec=gd->wav_ec;
  tscinit();
  RFB=9;
  if(ti->MainThread){ // forced sync
    GSM(RFB); GWM(RFB); 
    getwind();
    GSS(RFB); 
  }else{
    SSG(RFB); SWG(RFB);
  }
  setspec (1,td->wav_ec,td->iacb1,td->iace1,1)  ;
  setspec (1,td->wav_ec,td->iacb2,td->iace2,1)  ;
  RFB++;
  if(ti->MainThread){ // forced sync
    GWS(RFB); GSM(RFB);GWM(RFB);GSS(RFB); 
  }else{
    SSG(RFB);SWG(RFB); 
  }
#define CKCPU //if(ti->indg!=getcpuid()) printf("BZ:%d %2.2d %2.2d %2.2d %2.2d: %d %d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid(),__LINE__,td->nstep);
              //printf("BD:%d %2.2d %2.2d %2.2d %2.2d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid());
  if(ti->MainThread){ // forced sync
    RFG=1; gp_send_data(&gd->s_segs,sizeof(double)*mkj,cpuinf.NGrpPProc,RFG);
    t0=dclktime();
    //printf("BE:%d %2.2d %2.2d %2.2d %2.2d:%d %d %d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid(),td->nstep,pd.mTimeStep,pd.stop_now);
    for(;td->nstep<pd.mTimeStep;){
      CKCPU;
      gd->nstep++;
      td->nstep++;
      istep++;
      if(RFB>MAXRFB ){// reset Sync ID
        RFB=10; GWSR(RFB); GSM(RFB); GWMR(RFB); GSS(RFB);
      }
      RFB++; GSS(RFB);GSM(RFB);//S1 wind cpy end
      TBT(0);propagats_spec(td->wav_et,td->wav_ec,td->iacb1,td->iace1 );
      propagats_geo(td->wav_ec,td->wav_et,td->iacb1,td->iace1 );
      diffract_smooth(td->wav_et,td->wav_ec,td->iacb1,td->iace1 );
      TEBT(2);implschs(td->wav_ec,td->wav_et,td->iacb1,td->iace1 );
      TEBT(3);setspec (1,td->wav_ec,td->iacb1,td->iace1,2)  ;
      TEBT(4);accumea_(&ind,td->wav_ec,&td->iacb1,&td->iace1);
      TEBT(5);gp_check_recv(&gd->s_segs,cpuinf.NGrpPProc,RFG); //C0 wait boudary ok
      //printf("VW %2.2d %2.2d %4d %4d %4d %4d %4d\n",pd.mpi_id,ti->igrp,ti->ind,__LINE__,td->nstep,RFG,RFB);
      RFB++;  //S2 wait by gp_check_recv
#ifdef NO_MPI
      GSM(RFB);GWM(RFB);//if no mpi ,mainthread have no send/recv 
#endif
      GSS(RFB);
      TEBT(6);propagats_spec(td->wav_et,td->wav_ec,td->iacb2,td->iace2 );
      propagats_geo(td->wav_ec,td->wav_et,td->iacb2,td->iace2 );
      diffract_smooth(td->wav_et,td->wav_ec,td->iacb2,td->iace2 );
      TEBT(7);implschs (td->wav_ec,td->wav_et,td->iacb2,td->iace2 );
      TEBT(8);setspec  (1,td->wav_ec,td->iacb2,td->iace2,2);
      TEBT(9);accumea_ (&ind,td->wav_ec,&td->iacb2,&td->iace2);
      TEBT(10); 
      RFB++; //S3 //wait subthread compute end
      GWS(RFB);
      RFG++;
      TEBT(11);gp_send_data(&gd->s_segs,sizeof(double)*mkj,cpuinf.NGrpPProc,RFG);//C1 fixme:RFG
      GSS(RFB);
      if(gd->nostep==td->nstep){
        checkoutput_(&igrp,&ind,&RFB,td->wav_ec,&td->iacb1,&td->iace1,&td->iacb2,&td->iace2,(void*)&gd->hist_eot);
      }
      getwind();//double buff not need wait
      if(td->nstep==0){ TPT("gt"); }
      //gp_sdma_send_data(&gd->s_segs,sizeof(double)*mkj,cpuinf.NGrpPProc,RFG);
      TET(11);
      //must at after check flag
      if(CheckStop())break;
      //printf("AAE %2.2d %d\n",ti->igrp,RFB);
      //RFB++;GWS(RFB);GSM(RFB);GWM(RFB);GSS(RFB); //S4 //wait boudary ok
    }
    //printf("BF:%d %2.2d %2.2d %2.2d %2.2d:%d %d %d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid(),td->nstep,pd.mTimeStep,pd.stop_now);
    RFG++; gp_send_data(&gd->s_segs,sizeof(double)*mkj,cpuinf.NGrpPProc,RFG);
    TPH(); TPT("GT");
    tdd=dclktime();
    if(0)printf("RG%2.2d %2.2d %4d %f %d:%d %d end\n",ti->igrp,ti->ind,RFB,tdd-t0,pd.mTimeStep,
                td->iace1-td->iacb1,td->iace2-td->iacb2);
    RFB+=1000; //RFB-=10;GWS(RFB);RFB+=10;
    GSM(RFB); //RFB-=10;GWM(RFB); //if NEED Clean
    RFB+=10;GSS(RFB);
  }else{
    t0=dclktime();
    //printf("BE:%d %2.2d %2.2d %2.2d %2.2d:%d %d %d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid(),td->nstep,pd.mTimeStep,pd.stop_now);
    //if(ti->ind==1) printf("BB  %2.2d %d\n",ti->igrp,RFB);
    for(;td->nstep<pd.mTimeStep;){
      CKCPU;
      td->nstep++;
      istep++;
      if(RFB>MAXRFB ){// reset Sync ID
        RFB=10; SSG(RFB); SWGR(RFB);
      } 
      TBT(0);propagats_spec(td->wav_et,td->wav_ec,td->iacb1,td->iace1 );
      propagats_geo(td->wav_ec,td->wav_et,td->iacb1,td->iace1 );
      TEBT(1);RFB++; //SWG(RFB); //S1,not needed,wait wind ok,must be before implschs, will beready 
      TEBT(2);implschs(td->wav_ec,td->wav_et,td->iacb1,td->iace1 );
      TEBT(3);setspec (1,td->wav_ec,td->iacb1,td->iace1,2)  ;
      TEBT(4);accumea_(&ind,td->wav_ec,&td->iacb1,&td->iace1);
      TEBT(5);RFB++;SWG(RFB); //S2 //wait boudary ok
      TEBT(6);propagats_spec(td->wav_et,td->wav_ec,td->iacb2,td->iace2 );
      propagats_geo(td->wav_ec,td->wav_et,td->iacb2,td->iace2 );
      diffract_smooth(td->wav_et,td->wav_ec,td->iacb2,td->iace2 );
      TEBT(7);implschs (td->wav_ec,td->wav_et,td->iacb2,td->iace2 );
      TEBT(8);setspec  (1,td->wav_ec,td->iacb2,td->iace2,2);
      TEBT(9);RFB++;SSG(RFB);  //S3 compute end
      TEBT(10);
      TEBT(11);accumea_(&ind,td->wav_ec,&td->iacb2,&td->iace2);
      TET(11);
      SWG(RFB);
      if(gd->nostep==td->nstep){
        checkoutput_(&igrp,&ind,&RFB,td->wav_ec,&td->iacb1,&td->iace1,&td->iacb2,&td->iace2,(void*)&gd->hist_eot);
      }
      if(td->nstep==0){ TPT("st"); }
      if(CheckStop())break;
      //if(ti->ind==1) printf("BBE %2.2d %d\n",ti->igrp,RFB);
      //RFB++;SSG(RFB);GWS(RFB);  //S4 compute end
    }
    //printf("BF:%d %2.2d %2.2d %2.2d %2.2d:%d %d %d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid(),td->nstep,pd.mTimeStep,pd.stop_now);
    tdd=dclktime();
    if(0)printf("RG%2.2d %2.2d %4d %f %d end\n",ti->igrp,ti->ind,RFB,tdd-t0,pd.mTimeStep);
    TPT("ST");
    RFB+=1000; SSG(RFB);
    //RFB-=10;SWG(RFB);//not need wait
  }
}
void o_threadmain_(HTHREADINFO ti){
  int istep=0, RFB=1,RFG=1;
  int ind=ti->ind,indg=ti->indg,igrp=ti->igrp;
  double t0,tdd;
  if(!gd)gd=&pd; //for MTHREAD
  caliacbe();	  //not need wait	, ready here
  RFB=1;
  if(ti->MainThread){
    GSM(RFB); GWM(RFB); GSS(RFB);
  }else{
    SWG(RFB);
  }
  initgp();
  RFB=3;
  td->nstep=-1;
#define ILOG //if(istep<3){zlock_();printf("VV %2.2d %2.2d %2.2d %4d %4d\n",pd.mpi_id,ti->igrp,ti->ind,__LINE__,istep);fflush(stdout);zunlock_();}
  if(ti->MainThread){
    gd->nstep=-1;
    GWS(RFB); GSM(RFB);GWM(RFB);
    //usleep(5*1000000);
    g_initPropgats(); 
    g_initImplsch();
#ifdef USEXNL
    init_xnlg_();
#endif
    g_initSetspec();
    g_initOutput_cal();
    GSS(RFB); GSM(RFB);
  }else{ 
    SSG(RFB);SWG(RFB); 
  }
#ifdef USEXNL
  init_xnl_();
#else
  init_dia_();
#endif
  td->wav_et=gd->wav_et;
  td->wav_ec=gd->wav_ec;
  tscinit();
  RFB=9;
  if(ti->MainThread){ // forced sync
    GSM(RFB); GWM(RFB); 
    getwind();
    GSS(RFB); 
  }else{
    SSG(RFB); SWG(RFB);
  }
  RFB++;
  setspec (1,td->wav_ec,td->iacb1,td->iace1,1)  ;
  setspec (1,td->wav_ec,td->iacb2,td->iace2,1)  ;
  if(ti->MainThread){ // forced sync
    GWS(RFB); GSM(RFB); 
  }else{
    SSG(RFB); 
  }
#define CKCPU //if(ti->indg!=getcpuid()) printf("BZ:%d %2.2d %2.2d %2.2d %2.2d: %d %d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid(),__LINE__,td->nstep);
              //printf("BD:%d %2.2d %2.2d %2.2d %2.2d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid());
  if(ti->MainThread){ // forced sync
    RFG=1; gp_send_data(&gd->s_segs,sizeof(double)*mkj,cpuinf.NGrpPProc,RFG);
    t0=dclktime();
    //printf("BE:%d %2.2d %2.2d %2.2d %2.2d:%d %d %d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid(),td->nstep,pd.mTimeStep,pd.stop_now);
    for(;td->nstep<pd.mTimeStep;){
      CKCPU;
      gd->nstep++;
      td->nstep++;
      istep++;
      if(RFB>MAXRFB ){// reset Sync ID
        RFB=10; GWSR(RFB); GSM(RFB); GWMR(RFB); GSS(RFB);
      }
      RFB++; GSS(RFB);GSM(RFB);//S1 wind cpy end
      TBT(0);propagats_spec(td->wav_et,td->wav_ec,td->iacb1,td->iace1 );
      propagats_geo(td->wav_ec,td->wav_et,td->iacb1,td->iace1 );
      diffract_smooth(td->wav_et,td->wav_ec,td->iacb1,td->iace1 );
      TEBT(2);implschs(td->wav_ec,td->wav_et,td->iacb1,td->iace1 );
      TEBT(3);setspec (1,td->wav_ec,td->iacb1,td->iace1,2)  ;
      TEBT(4);accumea_(&ind,td->wav_ec,&td->iacb1,&td->iace1);
      TEBT(5);gp_check_recv(&gd->s_segs,cpuinf.NGrpPProc,RFG); //C0 wait boudary ok
      //printf("VW %2.2d %2.2d %4d %4d %4d %4d %4d\n",pd.mpi_id,ti->igrp,ti->ind,__LINE__,td->nstep,RFG,RFB);
      RFB++; GSS(RFB); //S2 wait by gp_check_recv
#ifdef NO_MPI
      GSM(RFB);GWM(RFB);//if no mpi ,mainthread have no send/recv 
#endif
      TEBT(6);propagats_spec(td->wav_et,td->wav_ec,td->iacb2,td->iace2 );
      propagats_geo(td->wav_ec,td->wav_et,td->iacb2,td->iace2 );
      diffract_smooth(td->wav_et,td->wav_ec,td->iacb2,td->iace2 );
      TEBT(7);implschs (td->wav_ec,td->wav_et,td->iacb2,td->iace2 );
      TEBT(8);setspec  (1,td->wav_ec,td->iacb2,td->iace2,2);
      TEBT(9);accumea_ (&ind,td->wav_ec,&td->iacb2,&td->iace2);
      TEBT(10); 
      RFB++; //S3 //wait subthread compute end
      GWS(RFB);
      if(gd->nostep==td->nstep){
        checkoutput_(&igrp,&ind,&RFB,td->wav_ec,&td->iacb1,&td->iace1,&td->iacb2,&td->iace2,(void*)&gd->hist_eot);
      }
      getwind();//double buff not need wait
      if(td->nstep==0){ TPT("gt"); }
      RFG++;
      TEBT(11);gp_send_data(&gd->s_segs,sizeof(double)*mkj,cpuinf.NGrpPProc,RFG);//C1 fixme:RFG
      //gp_sdma_send_data(&gd->s_segs,sizeof(double)*mkj,cpuinf.NGrpPProc,RFG);
      TET(11);
      //must at after check flag
      if(CheckStop())break;
    }
    //printf("BF:%d %2.2d %2.2d %2.2d %2.2d:%d %d %d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid(),td->nstep,pd.mTimeStep,pd.stop_now);
    RFG++; gp_send_data(&gd->s_segs,sizeof(double)*mkj,cpuinf.NGrpPProc,RFG);
    TPH(); TPT("GT");
    tdd=dclktime();
    if(0)printf("RG%2.2d %2.2d %4d %f %d:%d %d end\n",ti->igrp,ti->ind,RFB,tdd-t0,pd.mTimeStep,
                td->iace1-td->iacb1,td->iace2-td->iacb2);
    RFB+=1000; //RFB-=10;GWS(RFB);RFB+=10;
    GSM(RFB); //RFB-=10;GWM(RFB); //if NEED Clean
    RFB+=10;GSS(RFB);
  }else{
    t0=dclktime();
    //printf("BE:%d %2.2d %2.2d %2.2d %2.2d:%d %d %d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid(),td->nstep,pd.mTimeStep,pd.stop_now);
    for(;td->nstep<pd.mTimeStep;){
      CKCPU;
      td->nstep++;
      istep++;
      if(RFB>MAXRFB ){// reset Sync ID
        RFB=10; SSG(RFB); SWGR(RFB);
      } 
      TBT(0);propagats_spec(td->wav_et,td->wav_ec,td->iacb1,td->iace1 );
      propagats_geo(td->wav_ec,td->wav_et,td->iacb1,td->iace1 );
      TEBT(1);RFB++; //SWG(RFB); //S1,not needed,wait wind ok,must be before implschs, will beready 
      TEBT(2);implschs(td->wav_ec,td->wav_et,td->iacb1,td->iace1 );
      TEBT(3);setspec (1,td->wav_ec,td->iacb1,td->iace1,2)  ;
      TEBT(4);accumea_(&ind,td->wav_ec,&td->iacb1,&td->iace1);
      TEBT(5);RFB++;SWG(RFB); //S2 //wait boudary ok
      TEBT(6);propagats_spec(td->wav_et,td->wav_ec,td->iacb2,td->iace2 );
      propagats_geo(td->wav_ec,td->wav_et,td->iacb2,td->iace2 );
      diffract_smooth(td->wav_et,td->wav_ec,td->iacb2,td->iace2 );
      TEBT(7);implschs (td->wav_ec,td->wav_et,td->iacb2,td->iace2 );
      TEBT(8);setspec  (1,td->wav_ec,td->iacb2,td->iace2,2);
      TEBT(9);RFB++;SSG(RFB);  //S3 compute end
      TEBT(10);
      TEBT(11);accumea_(&ind,td->wav_ec,&td->iacb2,&td->iace2);
      TET(11);
      if(gd->nostep==td->nstep){
        checkoutput_(&igrp,&ind,&RFB,td->wav_ec,&td->iacb1,&td->iace1,&td->iacb2,&td->iace2,(void*)&gd->hist_eot);
      }
      if(td->nstep==0){ TPT("st"); }
      if(CheckStop())break;
    }
    //printf("BF:%d %2.2d %2.2d %2.2d %2.2d:%d %d %d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid(),td->nstep,pd.mTimeStep,pd.stop_now);
    tdd=dclktime();
    if(0)printf("RG%2.2d %2.2d %4d %f %d end\n",ti->igrp,ti->ind,RFB,tdd-t0,pd.mTimeStep);
    TPT("ST");
    RFB+=1000; SSG(RFB);
    //RFB-=10;SWG(RFB);//not need wait
  }
}
#else
//MTHREAD
void _threadmain_(HTHREADINFO ti){
  int istep=0, RFB=1;
  int ind=ti->ind,indg=ti->indg,igrp=ti->igrp;
  double t0,tdd;
  gd=&pd; //for MTHREAD
  //printf("BC:%d %2.2d %2.2d %2.2d %2.2d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid());
  caliacbe();	  //not need wait	, ready here
  RFB=1; SSM(RFB);SWM(RFB);
  initgp();
  RFB=3;
  td->nstep=-1;
#define ILOG //if(istep<3){zlock_();printf("VV %2.2d %2.2d %2.2d %4d %4d\n",pd.mpi_id,ti->igrp,ti->ind,__LINE__,istep);fflush(stdout);zunlock_();}
  SSM(RFB);SWM(RFB);
  getwind();//fixme
  td->wav_ec=gppar.wav_ec;
  td->wav_et=gppar.wav_et;
  setspec  (1,td->wav_ec,td->iacb1,td->iace1,1)  ;
  setspec  (1,td->wav_ec,td->iacb2,td->iace2,1)  ;
  RFB=10; SSM(RFB);SWM(RFB);
#ifdef USEXNL
  init_xnl_();
#else
  init_dia_();
#endif
  for(;td->nstep<pd.mTimeStep;){
    td->nstep++;
    if(RFB>MAXRFB ){
      RFB=10;SSM(RFB);SWMR(RFB);
    }
    RFB++; SSM(RFB); //S1
    //wait wind ok
    getwind();
    //wait boudary ok
    RFB++;SSM(RFB);SWM(RFB); //S2
    //RFB++;SSM(RFB);SWM(RFB);//o0
    TBT(4);propagats_spec(td->wav_et,td->wav_ec,td->iacb1,td->iace1 );
    TEBT(10);propagats_spec(td->wav_et,td->wav_ec,td->iacb2,td->iace2 );
    RFB++;SSM(RFB);SWM(RFB);//o1
    propagats_geo(td->wav_ec,td->wav_et,td->iacb1,td->iace1 );
    propagats_geo(td->wav_ec,td->wav_et,td->iacb2,td->iace2 );
    //RFB++;SSM(RFB);SWM(RFB);//o2
    diffract_smooth(td->wav_et,td->wav_ec,td->iacb1,td->iace1 );
    diffract_smooth(td->wav_et,td->wav_ec,td->iacb2,td->iace2 );
    //RFB++;SSM(RFB);SWM(RFB);//o3
    TEBT( 7);    implschs (td->wav_ec,td->wav_et,td->iacb1,td->iace1 );
    TEBT(11);    implschs (td->wav_ec,td->wav_et,td->iacb2,td->iace2);
    //RFB++;SSM(RFB);SWM(RFB);//o4
    TEBT(12);    setspec  (1,td->wav_ec,td->iacb1,td->iace1,2)  ;
    TEBT(12);    setspec  (1,td->wav_ec,td->iacb2,td->iace2,2);
    //RFB++;SSM(RFB);SWM(RFB);//o5
    TEBT(12);    accumea_ (&ind,td->wav_ec,&td->iacb1,&td->iace1);
    TEBT(12);    accumea_ (&ind,td->wav_ec,&td->iacb2,&td->iace2);
    RFB++;SSM(RFB);//S3
    if(gd->nostep==td->nstep){
      checkoutput_(&igrp,&ind,&RFB,td->wav_ec,&td->iacb1,&td->iace1,&td->iacb2,&td->iace2,(void*)&gd->hist_eot);
    }
    if(CheckStop())break;
  }
  RFB+=1000;SSM(RFB);usleep(100);SWM(RFB);
}

void o_threadmain_(HTHREADINFO ti){
  int istep=0, RFB=1;
  int ind=ti->ind,indg=ti->indg,igrp=ti->igrp;
  double t0,tdd;
  gd=&pd; //for MTHREAD
  //printf("BC:%d %2.2d %2.2d %2.2d %2.2d\n",pd.mpi_id,ti->igrp,ti->ind,ti->indg,getcpuid());
  caliacbe();	  //not need wait	, ready here
  RFB=1; SSM(RFB);SWM(RFB);
  initgp();
  RFB=3;
  td->nstep=-1;
#define ILOG //if(istep<3){zlock_();printf("VV %2.2d %2.2d %2.2d %4d %4d\n",pd.mpi_id,ti->igrp,ti->ind,__LINE__,istep);fflush(stdout);zunlock_();}
  SSM(RFB);SWM(RFB);
  getwind();//fixme
  td->wav_ec=gppar.wav_ec;
  td->wav_et=gppar.wav_et;
  setspec  (1,td->wav_ec,td->iacb1,td->iace1,1)  ;
  setspec  (1,td->wav_ec,td->iacb2,td->iace2,1)  ;
  RFB=10; SSM(RFB);SWM(RFB);
#ifdef USEXNL
  init_xnl_();
#else
  init_dia_();
#endif
  for(;td->nstep<pd.mTimeStep;){
    td->nstep++;
    if(RFB>MAXRFB ){
      RFB=10;SSM(RFB);SWMR(RFB);
    }
    RFB++; SSM(RFB); //1
    TBT(4);propagats_spec(td->wav_et,td->wav_ec,td->iacb1,td->iace1 );
    propagats_geo(td->wav_ec,td->wav_et,td->iacb1,td->iace1 );
    diffract_smooth(td->wav_et,td->wav_ec,td->iacb1,td->iace1 );
    //wait wind ok
    getwind();
    SWM(RFB);
    TEBT( 7);    implschs (td->wav_ec,td->wav_et,td->iacb1,td->iace1 );
    TEBT(12);    setspec  (1,td->wav_ec,td->iacb1,td->iace1,2)  ;
    TEBT(12);    accumea_ (&ind,td->wav_ec,&td->iacb1,&td->iace1);
    //wait boudary ok
    RFB++;SSM(RFB);SWM(RFB); //2
    TEBT(10);propagats_spec(td->wav_et,td->wav_ec,td->iacb2,td->iace2 );
    propagats_geo(td->wav_ec,td->wav_et,td->iacb2,td->iace2 );
    diffract_smooth(td->wav_et,td->wav_ec,td->iacb2,td->iace2 );
    TEBT(11);    implschs (td->wav_ec,td->wav_et,td->iacb2,td->iace2);
    TEBT(12);    setspec  (1,td->wav_ec,td->iacb2,td->iace2,2);
    TEBT(12);    accumea_ (&ind,td->wav_ec,&td->iacb2,&td->iace2);
    RFB++;SSM(RFB);//S3
    if(gd->nostep==td->nstep){
      checkoutput_(&igrp,&ind,&RFB,td->wav_ec,&td->iacb1,&td->iace1,&td->iacb2,&td->iace2,(void*)&gd->hist_eot);
    }
    if(CheckStop())break;
  }
  RFB+=1000;SSM(RFB);usleep(100);SWM(RFB);
}
#endif
#endif
