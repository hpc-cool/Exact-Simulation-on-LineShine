#define _GNU_SOURCE
#include <sched.h>
#include <pthread.h>
#include <stdio.h>
#include "mthread.h" //
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#ifdef DEBUG
//#define DBGSYNC
#endif

#ifdef USEHBM
// /pacific_fs/HPCKit/third/include/hbwmalloc.h
#include "hmalloc.h"
#define HMALLOC(n) hmalloc((n))
#define HREALLOC(p,n) hrealloc(p,(n))
#define HFREE(n) hfree((n))
#else
#define HREALLOC(p,n) realloc(p,(n))
#define HMALLOC(n) malloc((n))
#define HFREE(n) free((n))
#endif
#define HNEW(T)     (T*)HMALLOC(sizeof(T))
#define HNEWN(T,N)  (T*)HMALLOC(((long)sizeof(T))*(N))
#define HNEWPN(p,N) (  typeof(p))HMALLOC(((long)sizeof(*p))*(N))
#define HNEWZN(p,N) p=(typeof(p))HMALLOC(((long)sizeof(*p))*(N));if(p)memset((void*)p,0,sizeof(*p)*(N))
#ifndef SOLIDINFO
#define MSB     16
#define MSBG    1
#endif
//NGG >1: 2 level sync
//NGG =1: 1 level sync
#define NGG 6
//
#define MSWG
//#define USECSTATE
//MSB : group state distance,16*4 is a block
//MSBG: thread state in group distance,16*4 is a block
//thread set state for gmain ,ti->state
#ifdef USECSTATE
//sub thread state is in group data
#define PSSTATE(i) &gi->sstate[i*MSBG]
#define SSSTATE(i)  gi->sstate[i*MSBG]
#else
//sub thread state is in thread data
#define PSSTATE(i) &ti->threads[i].state
#define SSSTATE(i)  ti->state
#endif
struct TData;
__thread THREADINFO *ti=NULL;
__thread threadGroup *gi=NULL;
__thread extern struct TData*td;
__thread struct GData*gd=NULL;
int ThreadG=-1; //分组标志
#define DBG printf("%s %d %3.3d\n",__FILE__,__LINE__,mpi_id);fflush(stdout);

struct CPUINFO cpuinf={
  .OffClu=0,.SkipClu=1,.OffCore=0,.SkipCore=1,
  .NCorePClu=38,.MCorePClu=37,.NCluPNode =16,.NManageCore=4,
  .NThPGrp  =38,.NGrpPProc=16,.NProcPNode=1, .ManageCoreId=-1
};
int GetVInt(int volatile *volatile p);
void print_call_stack(int lock) ;
/* //这个函数需要再其它文件中实现,在本文件中可能会被优化
int __attribute__((noinline, optimize("O0")))
GetVInt(int volatile *volatile p){
  return *p;
}*/
int ntdelay(int n,int vt){
  for(int i=0;i<n;i++) vt=(vt*1357+2581);
  return vt;
}
static inline uint64_t atsc(void) {
#ifdef __ARM_ARCH
  uint64_t tsc;
  asm volatile("mrs %0, cntvct_el0" : "=r" (tsc));    //读取系统时间戳
  return tsc;//most tsc is 56bits
  //return tsc&0xFFFFFFFFFFFFFFL;//most tsc is 56bits
#else
  uint64_t a,d;
  asm volatile("rdtsc " : "=a" (a),"=d"(d));    //读取系统时间戳
  return (d<<32)|a;
#endif
}
#define MAXCHECK 10000
#define VLINE ,line
//#define PLINE ,int line
#define FVLINE ,*line
//ntdelay(10) 2ns
//1000*2ns=2us 10000*2ns*10 200us ,100 000 *1us*1000 =100s
static void DBGH(){
  printf("DBGHERE\n");
}
void zunlock_();
void zlock_();
#define FWAIT(NM,COP) int Wait_##NM(volatile int *pv,int cv,int nt PLINE){ \
  int i=0,vt=cv,val,delay=20*nt; uint64_t t0=atsc();\
  for(;((val=GetVInt(pv)) COP cv);){ \
    vt=ntdelay(delay,vt);\
    if(i++>1000){\
      for(;GetVInt(pv) COP cv;){\
        vt=ntdelay(delay*10,vt);\
        if(i++>10000){\
          for(delay=10*nt;GetVInt(pv) COP cv;){\
            if(i++>=1000000){\
              uint64_t t1=atsc()-t0;\
              zlock_();\
              printf("wait state timeout %d %s %d %d %x %10.6fs %d\n",val,#COP,cv,nt,vt,t1/100000000. VLINE);\
							print_call_stack(0); fflush(stdout);\
              zunlock_();usleep(1000*1000*10);\
              exit(0);\
            }\
            usleep(delay);if(((i&0x3f)==0)&&(delay<10000))delay++;\
          }\
        }\
      }\
    }\
  }\
  return vt;\
}
int Wait_LGE( volatile int *pv,int cv,int nt ,int line);
FWAIT(LNE,==);
FWAIT(LEQ,!=);
FWAIT(LGE,<);
FWAIT(LLE,>);
static void zStartThreads(TFunc tfun,void*para,int detach,int clear);
void initthreads_(int *mpi_id_, struct CPUINFO*ci ,int *err){
  *err=InitThreads(*mpi_id_,ci);
}
void startthreads_(){
  StartThreads();
}
void endthreads_(){
  EndThreads();
}
void SetLocV(int typ,int ind,void*p){
  struct LocVar *pl;
  if(typ==0) pl=&ti->locv;
  else if(typ==1) pl=&gi->locv;
  else if(typ==2) pl=&md.locv;
  if(ind<0||ind>100) ind=0;
  if(ind>=pl->nlocv)pl->nlocv=ind+1;
  if(ind<8){
    pl->locv[ind]=p;
  }else{
    void ***pp=(void***)&pl->locv[8];
    if(pl->mlocv<=ind){
      int ne=ind-pl->mlocv;
      ne=(ne+7)&(~7);

      pl->mlocv=ind+ne;
      *pp=HREALLOC(*pp,pl->mlocv-8);
      memset(*pp+pl->mlocv-8-ne,0,(ne)*sizeof(void*));
    }
    pp[ind-8]=p;
  }
}
void *GetLocV(int typ,int ind){
  struct LocVar *pl;
  if(typ==0) pl=&ti->locv;
  else if(typ==1) pl=&gi->locv;
  else if(typ==2) pl=&md.locv;
  if(ind<0||ind>100) return pl->locv[0];
  if(pl->nlocv<=ind) return NULL;
  if(ind<8) return pl->locv[ind];
  void**pp=(void**)pl->locv[8];
  return pp[ind-8];
}
int InitThreads(int mpi_id_, struct CPUINFO*ci){
  mpi_id=mpi_id_;
  cpuinf=*ci;
  if(cpuinf.NGrpPProc>1)ThreadG=1;
  else {ThreadG=0;cpuinf.NGrpPProc=1;}
  if(cpuinf.NThPGrp>cpuinf.MCorePClu){
    printf("NThPGrp %d is biger than MCorePClu %d,reset it\n",cpuinf.NThPGrp,cpuinf.MCorePClu);
    cpuinf.NThPGrp=cpuinf.MCorePClu;
  }
#ifdef SOLIDINFO
  if(cpuinf.NThPGrp>MCOREPC){
    printf("NThPGrp %d biger than surported MCOREPC %d\n",cpuinf.NThPGrp,MCOREPC);
    return (-1);
  }
  if(cpuinf.NGrpPProc>MCLUST){
    printf("NGrpPProc %d biger than surported MCLUST %d\n",cpuinf.NGrpPProc,MCLUST);
    return (-3);
  }
  if(cpuinf.NProcPNode > MCLUST){
    printf("NProcPNode %d biger than surported MCLUST %d\n", cpuinf.NProcPNode , MCLUST);
    return (-5);
  }
#endif
  if(mpi_id==0){
    printf("NGrpPProc=%d,NThPGrp=%d,NProcPNode=%d,ManageCoreId=%d\n",
           cpuinf.NGrpPProc,cpuinf.NThPGrp,cpuinf.NProcPNode,cpuinf.ManageCoreId);
    printf("OffClu  =%d, SkipClu =%d, OffCore =%d, SkipCore=%d\n",
           cpuinf.OffClu  , cpuinf.SkipClu , cpuinf.OffCore , cpuinf.SkipCore);
  }
  initmd();
  return 0;
}
static void clearlocv(struct LocVar *pl);
void threadMain(HTHREADINFO pti){
  ti=pti;td=pti->td;
  gi=ti->pg;
  if(gi)gd=gi->gd;
  bindthread();
  opentf();
  _threadmain_(pti);
  ti->threadid=0;
  clearlocv(&ti->locv);
}
void StartThreads(){
  zStartThreads(threadMain,0,1,1);
}
void EndThreads(){
  zStartThreads(NULL,0,1,1);
}
int bindcpu(int id){
  cpu_set_t mask;  //CPU核的集合
  CPU_ZERO(&mask);    //置空
  CPU_SET(id,&mask);   //设置亲和力值
  if (sched_setaffinity(0, sizeof(mask), &mask) == -1){//设置线程CPU亲和力
    printf("warning: could not set CPU affinity, continuing...\n");
    return -1;
  }
  return 0;
}
void bindthread(){
  if(ti) bindcpu(ti->indg);
}
void bindthread_(){
  bindthread();
}
threadProc md={0};
//#define DBGSYNC
//#define THLOG
#ifdef DBGSYNC
#define THLOG
#endif
void opentf(){
#ifdef THLOG
  if(!ti->fo){
    char fn[256];
    if(ThreadG) sprintf(fn,"MTW_%2.2d_%2.2d_%2.2d.log",mpi_id,ti->igrp,ti->ind);
    else sprintf(fn,"MTW_%d.log",ti->ind);
    ti->fo=fopen(fn,"wt");
  }
#endif
}
//SM
void sWaitMain(int state,int nt PLINE){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"TWT:3.2d %3.2d %3.2d ... ",md.nstep,state,ti->ind);
    fflush(ti->fo);
  }
#endif
  Wait_LGE(&md.state,state,nt VLINE);
#ifdef DBGSYNC
  fprintf(ti->fo,"twt end\n");fflush(ti->fo);
#endif
}
void sWaitMainr(int state,int nt PLINE){
  Wait_LLE(&md.state,state,nt VLINE);
}
void sSetMain(int state PLINE){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"TST:%3.2d %3.2d %3.2d %3.2d ...",md.nstep,ti->state,state,ti->ind);
    fflush(ti->fo);
  }
#endif
  ti->state=state;
#ifdef DBGSYNC
  fprintf(ti->fo,"tst end\n");fflush(ti->fo);
#endif
}
//MS
void mWaitSubs(int state,int nt PLINE){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"MWT:%3.2d %3.2d:",md.nstep,state);
    for(int i=0;i<md.Nthreads;i++){
      fprintf(ti->fo,"%3.2d ",md.threads[i].state);
    }
    fprintf(ti->fo," ...  "); fflush(ti->fo);
  }
#endif
  for(int i=0;i<md.Nthreads;i++){
#ifdef DBGSYNC
    fprintf(ti->fo,"V %d:%3.2d ",i,md.threads[i].state); fflush(ti->fo);
#endif
    Wait_LGE(&md.threads[i].state,state,nt VLINE);
  }
#ifdef DBGSYNC
  fprintf(ti->fo,"mws end\n");fflush(ti->fo);
#endif
}
void mWaitSubsr(int state,int nt PLINE){
  for(int i=0;i<md.Nthreads;i++){
    Wait_LLE(&md.threads[i].state,state,nt VLINE);
  }
}
void mSetSubs(int state PLINE){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"MST:%3.2d %3.2d %3.2d...",md.nstep,md.state,*state);
    fflush(ti->fo);
  }
#endif
  md.state=state;
#ifdef DBGSYNC
  fprintf(ti->fo,"mss end\n");fflush(ti->fo);
#endif
}
//SG thread wait gmain state,pg->state
void sWaitGrp(int state,int nt PLINE){
#ifdef DBGSYNC
  fprintf(ti->fo,"SWG:%3.2d %3.2d %3.2d ...",md.nstep,state,ti->ind);
  fflush(ti->fo);
#endif
#ifdef MSWG
  Wait_LGE(&ti->gstate,state,nt VLINE);
#else
  Wait_LGE(&gi->state,state,nt VLINE);
#endif
#ifdef DBGSYNC
  fprintf(ti->fo,"swg end\n");fflush(ti->fo);
#endif
}
void sWaitGrpr(int state,int nt PLINE){
#ifdef MSWG
  Wait_LLE(&ti->gstate,state,nt VLINE);
#else
  Wait_LLE(&gi->state,state,nt VLINE);
#endif
}
void sSetGrp(int state PLINE){
#ifdef DBGSYNC
  fprintf(ti->fo,"SSG:%3.2d %3.2d %3.2d %3.2d :",md.nstep,ti->pg->sstate[ti->ind*MSBG],state,ti->ind);
#endif
  if(NGG>1){
#ifdef DBGSYNC
    fprintf(ti->fo,"gWS S:%3.2d %3.2d:",md.nstep,state);
    for(int i= ti->sib;i<ti->sie;i++){
      fprintf(ti->fo,"%3.2d ",SSSTATE(i));
    }
    fprintf(ti->fo," ..."); fflush(ti->fo);
#endif
    for(int i= ti->sib;i<ti->sie;i++){
      Wait_LGE(PSSTATE(i),state,1 VLINE);
    }
  }
  SSSTATE(ti->ind)=state;
#ifdef DBGSYNC
  fprintf(ti->fo,"ssg end\n");fflush(ti->fo);
#endif
}
//GS gmain thread wait sub threads,threads[x].state
void gWaitSubs(int state,int nt PLINE){//mt
#ifdef DBGSYNC
  fprintf(ti->fo,"GWS S:%3.2d %3.2d:",md.nstep,state);
#endif
  if(NGG>1){
#ifdef DBGSYNC
    for(int i= ti->sib;i<ti->sie;i++){
      fprintf(ti->fo,"%3.2d ",SSSTATE(i));
    }
    fprintf(ti->fo," ... "); fflush(ti->fo);
#endif
    for(int i= ti->sib;i<ti->sie;i++){
      Wait_LGE(PSSTATE(i),state,nt VLINE);
    }
  }
#ifdef DBGSYNC
  fprintf(ti->fo,"gws s end\nGWS:%3.2d %3.2d:",md.nstep,state);
  fflush(ti->fo);
  for(int i=NGG;i<ti->Nthreads;i+=NGG){
    fprintf(ti->fo,"%3.2d ",SSSTATE(i));
  }
  fprintf(ti->fo,"... "); fflush(ti->fo);
#endif
  for(int i=NGG;i<ti->Nthreads;i+=NGG){
    Wait_LGE(PSSTATE(i),state,nt VLINE);
  }
#ifdef DBGSYNC
  fprintf(ti->fo,"gws end\n");fflush(ti->fo);
#endif
}
void gWaitSubsr(int state,int nt PLINE){//mt
  for(int i= ti->sib;i<ti->sie;i++){
    Wait_LLE(PSSTATE(i),state,nt VLINE);
  }
  for(int i=NGG;i<ti->Nthreads;i+=NGG){
    Wait_LLE(PSSTATE(i),state,nt VLINE);
  }
}
//gmain thread set for subthread ,pg->state
void gSetSubs (int state PLINE){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"GSS:%3.2d %3.2d %3.2d %3.2d\n",md.nstep,gi->state,state,ti->ind);
    fflush(ti->fo);
  }
#endif
#ifdef MSWG
  for(int i=0;i<ti->Nthreads;i++){
    ti->threads[i].gstate=state;
  }
#else
  gi->state=state;
#endif
}

//GM gmain thread wait mmthreads,pg->gstate
void gWaitMain(int state,int nt PLINE){//mt
  int i;
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"GWM:%3.2d %3.2d %3.2d:... ",md.nstep,state,ti->igrp); fflush(ti->fo);
  }
#endif
  Wait_LGE(&gi->gstate,state,nt VLINE);
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"gwm:%3.2d %3.2d %d %3.2d end\n",md.nstep,state,ti->pg->gstate,ti->igrp); fflush(ti->fo);
  }
#endif
}
void gWaitMainr(int state,int nt PLINE){//mt
  Wait_LLE(&gi->gstate,state,nt VLINE);
}
//gmain thread set for mmthread ,md.gstate[(x)*MSB]
void gSetMain(int state PLINE){
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"GSM:%3.2d %3.2d %3.2d %3.2d ...",md.nstep,md.gstate[(ti->igrp)*MSB],state,ti->igrp);
    fflush(ti->fo);
  }
#endif
  md.gstate[(ti->igrp)*MSB]=state;
#ifdef DBGSYNC
  if(ti->fo) {
    printf(ti->fo,"gsm:%3.2d %3.2d %3.2d %3.2d end\n",md.nstep,md.gstate[(ti->igrp)*MSB],state,ti->igrp);
    fflush(ti->fo);
  }
#endif
}
//called by main main threads
//mmthread wait gmain thread,md.gstate[(x)*MSB]
void mWaitGrps(int state,int nt PLINE){//mmt
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"MWG:%d %d %d:",md.nstep,state,md.ngrp);
    for(int i=0;i<md.ngrp;i++){ fprintf(ti->fo,"%d ",md.gstate[(i)*MSB]); }
    fprintf("\n ");
    fflush(ti->fo);
  }
#endif
  for(int i=0;i<md.ngrp;i++){
    Wait_LGE(&md.gstate[(i)*MSB],state,nt VLINE);
  }
}
void mWaitGrpsr(int state,int nt PLINE){//mmt
  for(int i=0;i<md.ngrp;i++){
    Wait_LLE(&md.gstate[(i)*MSB],state,nt VLINE);
  }
}
//mmthread set for gmain threads,md.grps[x].gstate
void mSetGrps(int state PLINE){// mmt
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"MSG:%3.2d %3.2d %3.2d ...",md.nstep,md.state,state);
    fflush(ti->fo);
  }
#endif
  md.state=state;
  for(int i=0;i<md.ngrp;i++){
    md.grps[i]->gstate=state;
  }
#ifdef DBGSYNC
  if(ti->fo) {
    fprintf(ti->fo,"MSG:%3.2d %3.2d ...",md.nstep,state);
    for(int i=0;i<md.ngrp;i++){
      fprintf(ti->fo,"%d:%3.2d ",i,md.grps[i]->gstate);
    }
    fprintf(ti->fo,"msg end\n");fflush(ti->fo);
  }
#endif
}

// Fortran interface
//SM
void swaitmain_ (int *state,int*nt FPLINE){ sWaitMain (*state,*nt FVLINE); }
void swaitmainr_(int *state,int*nt FPLINE){ sWaitMainr(*state,*nt FVLINE); }
void ssetmain_  (int *state        FPLINE){ sSetMain  (*state     FVLINE); }
//MS
void mwaitsubs_ (int *state,int*nt FPLINE){ mWaitSubs (*state,*nt FVLINE); }
void mwaitsubsr_(int *state,int*nt FPLINE){ mWaitSubsr(*state,*nt FVLINE); }
void msetsubs_  (int *state        FPLINE){ mSetSubs  (*state     FVLINE); }
//SG
void swaitgrp_  (int *state,int*nt FPLINE){ sWaitGrp  (*state,*nt FVLINE); }
void swaitgrpr_ (int *state,int*nt FPLINE){ sWaitGrpr (*state,*nt FVLINE); }
void ssetgrp_   (int *state        FPLINE){ sSetGrp   (*state     FVLINE); }
//GS
void gwaitsubs_ (int *state,int*nt FPLINE){ gWaitSubs (*state,*nt FVLINE); }
void gwaitsubsr_(int *state,int*nt FPLINE){ gWaitSubsr(*state,*nt FVLINE); }
void gsetsubs_  (int *state        FPLINE){ gSetSubs  (*state     FVLINE); }
//GM
void gwaitmain_ (int *state,int*nt FPLINE){ gWaitMain (*state,*nt FVLINE); }
void gwaitmainr_(int *state,int*nt FPLINE){ gWaitMainr(*state,*nt FVLINE); }
void gsetmain_  (int *state        FPLINE){ gSetMain  (*state     FVLINE); }
//MG
void mwaitgrps_ (int *state,int*nt FPLINE){ mWaitGrps (*state,*nt FVLINE); }
void mwaitgrpsr_(int *state,int*nt FPLINE){ mWaitGrpsr(*state,*nt FVLINE); }
void msetgrps_  (int *state        FPLINE){ mSetGrps  (*state     FVLINE); }

void tscinit(){
  memset(ti->tscds,0,sizeof(ti->tscds));
}
void tscb(int id){
  ti->tscbs[id]=atsc();
}
void tsce(int id){
  uint64_t d=atsc();
  ti->tscds[id]+=((d-ti->tscbs[id]));
}
void tsceb(int id){
  uint64_t d=atsc();
  ti->tscds[id-1]+=((d-ti->tscbs[id-1]));
  ti->tscbs[id]=d;
}
void tscb_(int *id_){
  int id=*id_;
  ti->tscbs[id]=atsc();
}
void tsce_(int *id_){
  int id=*id_;
  uint64_t d=atsc();
  ti->tscds[id]+=((d-ti->tscbs[id]));
}
void tsceb_(int *id_){
  int id=*id_;
  uint64_t d=atsc();
  ti->tscds[id-1]+=((d-ti->tscbs[id-1]));
  ti->tscbs[id]=d;
}
void prtsc(const char*tag){
  char buf[1024]={0};
  for(int i=0;i<16;i++){
    sprintf(buf+i*10,"|%9.5f     Z",ti->tscds[i]*(1./(100*1000*1000)));
  }
  fprintf(ti->fo,"%2s(s):%2.2d %2.2d %2d:%s\n",tag,md.mpi_id,ti->igrp,ti->ind,buf);
  fflush(ti->fo);
}
void prtsc_(){
  prtsc("MT");
}
#define VFREE(p) if(p)free(p)
  
static void clearlocv(struct LocVar *pl){
  if(pl->nlocv){
    int i;
    if(pl->nlocv>8){
      void**p=(void**)pl->locv[8];
      VFREE(p);pl->locv[8]=NULL;
      pl->nlocv=8;
    }
  }
  pl->nlocv=0;
}
static void ClearThread(HTHREADINFO pti){
  if(pti->threadid){
    void*vres;
    pthread_cancel(pti->threadid);
    pthread_join(pti->threadid,&vres);
    //pti->pg->sstate[pti->ind*MSBG]=0;
    pti->threadid=0;
    clearlocv(&pti->locv);
  }
}
void initmd(){
  if(md.PThreadInited) return;
  int ncores=cpuinf.NCorePClu*cpuinf.NCluPNode;
  memset(&md,0,sizeof(md));
  md.PThreadInited=1;
  ti=&md.tm;
  if(ThreadG){
    md.ngrp=cpuinf.NGrpPProc;
    int mg=0;
    int mores=cpuinf.NCorePClu*cpuinf.NGrpPProc;
    if(cpuinf.NProcPNode<1)cpuinf.NProcPNode=1;
    int ipr=mpi_id%cpuinf.NProcPNode;
    int bclu=ipr*cpuinf.NGrpPProc+cpuinf.OffClu;//begin of cluster
    int nclupm=0;//clusters per MangeCore
    if(cpuinf.NManageCore)nclupm=cpuinf.NCluPNode/cpuinf.NManageCore;
    //set to managecore
    if(cpuinf.ManageCoreId <=-2)cpuinf.ManageCoreId =ncores;
    else if(cpuinf.ManageCoreId==-1)cpuinf.ManageCoreId=0;
    //manage core ind at cluster
    int ic=cpuinf.ManageCoreId%cpuinf.NCorePClu;
    if(cpuinf.ManageCoreId>=cpuinf.NCorePClu){
      // Manage core is not in first cluster
      if(cpuinf.ManageCoreId>=ncores){
        if(nclupm){
          //use Managecore
          mg=-1; cpuinf.ManageCoreId = ncores+bclu/nclupm;
        }else cpuinf.ManageCoreId =0; //no Managecore
      }else {
        if(ic>=cpuinf.NThPGrp){//cpuinf.MCorePClu
          // Managecore is reserved core
          mg=-1;
        }else {
          // index of cluster with manage core
          mg=cpuinf.ManageCoreId/cpuinf.NCorePClu;
          //not allow ManageCore inner groupcore
          //cpuinf.ManageCoreId=mg*cpuinf.NCorePClu;
        }
      }
    }else {//manage core in first cluster
      if(ic>=cpuinf.NThPGrp){//fixme not surport more than 1 proc in 1 cluster
        mg=-1; // Managecore is reserved core
      }else {
        mg=0;
        //cpuinf.ManageCoreId=0; //not allow ManageCore inner groupcore
      }
    }
    //ManageCoreId from ind in proc to node
    cpuinf.ManageCoreId+=bclu*cpuinf.NCorePClu;
    bindcpu(cpuinf.ManageCoreId);
#ifndef SOLIDINFO
    HNEWZN(md.grps,md.ngrp);
    HNEWZN(md.gstate,md.ngrp*MSB);
#endif
    int nths=0;
    for(int i=0;i<cpuinf.NGrpPProc;i++){
      int iclu=(bclu+i*cpuinf.SkipClu)%cpuinf.NCluPNode;
      bindcpu(iclu*cpuinf.NCorePClu);// bind to core in grp for hbm
      threadGroup *pgi;
      HNEWZN(pgi,1); md.grps[i]=pgi;
      pgi->pmd=&md;
      pgi->ngrp=cpuinf.NGrpPProc;
      pgi->Nthreads=cpuinf.NThPGrp;
      pgi->idMainThread=0;
      pgi->iclu=iclu;
      pgi->igrp=i;
      int im=0;
      if(mg==i){
        if(pgi->Nthreads>=cpuinf.MCorePClu)pgi->Nthreads--;
      }
      int cbase=iclu*cpuinf.NCorePClu;
      pgi->NSpecial=0;
      int nsize=_getgdsize_();
      if(nsize){
        pgi->gd=(struct GData*)HNEWN(char,nsize);
        memset(pgi->gd,0,nsize);
      }
#ifndef SOLIDINFO
      HNEWZN(pgi->threads,pgi->Nthreads);
      HNEWZN(pgi->sstate,pgi->Nthreads*MSBG);
#endif
      int icore=cpuinf.OffCore;
      for(int j=0;j<pgi->Nthreads;j++,icore+=cpuinf.SkipCore){
        struct _THREADINFO *pti=pgi->threads+j;
        int nsize=_gettdsize_();
        if(nsize)pti->td=HNEWN(char,nsize);
        icore=icore%cpuinf.MCorePClu;
        if(cbase+icore==cpuinf.ManageCoreId)icore=(icore+1)%cpuinf.MCorePClu;
        pti->pg=pgi;
        pti->threads=pgi->threads;
        pti->Nthreads=pgi->Nthreads;
        pti->MainThread=(j==0);//pgi->idMainThread);
        pti->ind=j;
        pti->indg=cbase+icore;
        pti->igrp=pgi->igrp;
        pti->iclu=pgi->iclu;
        pti->indp=nths++;
        if(pti->ind%NGG==0){
          pti->sib=pti->ind+1;
          pti->sie=pti->sib+NGG-1;
          if(pti->sie>pti->Nthreads)pti->sie=pti->Nthreads;
        } else {
          pti->sib=0; pti->sie=-1;
        }
      }
    }
    md.Nthreads=nths;
    bindcpu(cpuinf.ManageCoreId);
    int nsize=_gettdsize_();
    if(nsize)td=ti->td=HNEWN(char,nsize);
  }else{
    md.idMainThread=1;
    //proc ind in node
    int ipr=mpi_id%cpuinf.NProcPNode;
    //groups per cluster
    int ng=cpuinf.MCorePClu/cpuinf.NThPGrp;
    // 1 group per cluster only
    if(cpuinf.NProcPNode<=cpuinf.NCluPNode)ng=1;
    //max core per group
    int mcore=cpuinf.MCorePClu/ng;
    //`cluster ind
    int iclu=((ipr/ng)*cpuinf.SkipClu+cpuinf.OffClu)%cpuinf.NCluPNode;
    //can not use managecore,
    int cbase=iclu*cpuinf.NCorePClu+(ipr%ng)*mcore;
    if(cpuinf.ManageCoreId<0||cpuinf.ManageCoreId>=cpuinf.NCorePClu)cpuinf.ManageCoreId=0;
    if(ng>1||cpuinf.ManageCoreId<cpuinf.MCorePClu)cpuinf.ManageCoreId=0;
    if(cpuinf.ManageCoreId==0){
      cpuinf.ManageCoreId=cbase;
    }
    int nth=cpuinf.NThPGrp;
    int icore=cbase;
    if(cpuinf.ManageCoreId==cbase&&nth==mcore)nth--;
    md.Nthreads=nth;
    bindcpu(cpuinf.ManageCoreId);
#ifndef SOLIDINFO
    HNEWZN(md.threads,nth);
#endif
    for(int j=0;j<nth;j++,icore++){
      struct _THREADINFO *pti=md.threads+j;
      int nsize=_gettdsize_();
      if(nsize)pti->td=HNEWN(char,nsize);
      if(icore==cpuinf.ManageCoreId) icore++;
      pti->ind=j;
      pti->indp=j;
      pti->indg=icore;
      pti->Nthreads=nth;
      pti->MainThread=(j==0);
      if(pti->ind%NGG==0){
        pti->sib=pti->ind+1;
        pti->sie=pti->sib+NGG-1;
        if(pti->sie>pti->Nthreads)pti->sie=pti->Nthreads;
      } else {
        pti->sib=0; pti->sie=-1;
      }
    }
    int nsize=_gettdsize_();
    if(nsize)ti->td=HNEWN(char,nsize);
  }
}
static void InitThread(HTHREADINFO pti,TFunc tfun,void*para,int detach){
  ClearThread(pti);
  pti->tfun=tfun;
  pti->para=para;
  pti->detach=detach;
  if(tfun==NULL)return;
}

static int netis=0,metis=0;
static THREADINFO *etis=NULL;
void eThreadMain(HTHREADINFO pti){
  ti=pti;td=pti->td;
  gi=ti->pg;gd=gi->gd;
  bindthread();
  ti->tfun(ti->para);
  ti->threadid=0;
  clearlocv(&ti->locv);
}
void zStartEThread(TFunc tfun,void*para,int id){
  int eid=netis;
  if(eid>=metis){
    metis+=4;
    etis=(THREADINFO*)realloc(etis,sizeof(THREADINFO )*metis);
    memset(etis+netis,0,sizeof(THREADINFO )*(metis-netis));
  }
  netis++;
  THREADINFO*pti= &etis[eid];
  if(id>=0) pti->indg=id;
  else pti->indg=cpuinf.ManageCoreId;
  pti->tfun=tfun;
  pti->para=para;
  pthread_create(&pti->threadid,NULL,(TSFunc)eThreadMain,pti);
}
void zstartethread_(TFunc tfun,void*para,int *id){
  zStartEThread(tfun,para,*id);
}
void zStartThreads(TFunc tfun,void*para,int detach,int clear){
  int ib,ie;
  HTHREADINFO tis;
  if(tfun){
    //printf("%d start threads %p\n",mpi_id,tfun);
    ti=&md.tm;
    ti->threadid=0;
    ti->igrp=-1;
    ti->indg=cpuinf.ManageCoreId;
    InitThread(&md.tm,0,0,0);
    if(ThreadG){
      for(int j=0;j<cpuinf.NGrpPProc;j++){
        threadGroup *pgi=md.grps[j];
        tis=pgi->threads;
        for(int i=0;i<pgi->Nthreads;i++){
          InitThread(tis+i,tfun,para,detach);
          pthread_create(&tis[i].threadid,NULL,(TSFunc)tfun,&tis[i]);
          if(detach)pthread_detach(tis[i].threadid);
        }
      }
    }else{
      tis=md.threads;
      md.NSpecial=0;
      for(int i=0;i<md.Nthreads;i++){
        tis[i].Nthreads=md.Nthreads;
      }
      for(int i=0;i<md.Nthreads;i++){
        InitThread(tis+i,tfun,para,detach);
        pthread_create(&tis[i].threadid,NULL,(TSFunc)tfun,&tis[i]);
        if(detach)pthread_detach(tis[i].threadid);
      }
    }
  }else{
    usleep(100*1000);//100ms
    if(ThreadG){
      for(int j=0;j<cpuinf.NGrpPProc;j++){
        threadGroup *pgi=md.grps[j];
        tis=pgi->threads;
        for(int i=0;i<pgi->Nthreads;i++){
          ClearThread(tis+i);
        }
      }
    }else{
      tis=md.threads;
      for(int i=0;i<md.Nthreads;i++){
        ClearThread(tis+i);
      }
    }
  }
}
void threads_abort(){
	EndThreads();
}
//=================================
int getgid_(){ return ti->igrp; }
int gettid_(){ return ti->ind; }
int getnthreads_(){ return ti->Nthreads+1; }
