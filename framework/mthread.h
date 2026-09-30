#ifndef MTHREAD_H_INCLUDED
#define MTHREAD_H_INCLUDED
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#if 1
struct CPUINFO{
  int NCorePClu   ;// 每簇(NUMA)核数(38)
  int MCorePClu   ;// 每簇(NUMA)可用核数(38),假设不可用核在最后
  int NCluPNode   ;// 每节点簇数(16)
  int NManageCore ;// 管理核核数，编号>NCorePClu*NCluPNode    

  int NThPGrp     ;// 每组启动线程数(38)
  int NGrpPProc   ;// 每进程组数(16)
  int NProcPNode  ;// 每节点进程数(38)
  int ManageCoreId;// 管理核(-1)

  int OffClu      ;// 第一组使用的簇(0)
  int SkipClu     ;// 几个簇启动1个组间(1)
  int OffCore     ;// 每个组使用的第一个核(0)
  int SkipCore    ;// 几个核启动一个线程(1)
};
#endif
extern int mpi_id,NThreads;
extern struct CPUINFO cpuinf;
//#define SOLIDINFO
#ifdef SOLIDINFO
#ifndef MCOREPC 
#define MCOREPC 64
#endif
#ifndef MCLUST  
#define MCLUST  20
#endif
#ifndef MSB
#define MSB     16
#endif
#ifndef MSBG
#define MSBG    1
#endif
#endif 
typedef struct _THREADINFO THREADINFO,*HTHREADINFO;
typedef void (*TFunc)(HTHREADINFO);
typedef void *(*TSFunc)(void*);

void initmd();
struct LocVar{
  int nlocv,mlocv;
  void*locv[9];
};
typedef struct _THREADINFO{
  short ind,indp,Nthreads,MainThread,igrp,iclu;      // 4 * 2 bytes
  short detach,indg,sib,sie;   // 3 * 2 bytes
  int nstep;
  int volatile state;     //MWS SSM //GWS SSG
  int volatile gstate;     //GSS SWG
  uint64_t tscbs[16];
  uint64_t tscds[16];
  struct _threadGroup *pg;
  struct _THREADINFO  *threads;
  TFunc tfun;
  void *para;
  pthread_t threadid;
  FILE *fo;
  void*td;
  struct LocVar locv;
}THREADINFO;
struct GData;
typedef struct _threadGroup{
  short igrp,iclu,ngrp,mainGroup,Nthreads,idMainThread,NSpecial;
  short PThreadInited;
  int nstep;
  int volatile state;         //SWG GSS 
  int volatile gstate;        //GWM MSG 
  struct _threadProc *pmd;
  struct GData*gd;
#ifdef SOLIDINFO
  int volatile sstate[(MCOREPC)*MSBG]; //GWS GSM 
  //int volatile gstates[(MCOREPC)*MSBG]; //GWM MSG 
  THREADINFO threads[MCOREPC];
#else
  int volatile *sstate;//[(MCOREPC)*MSBG]; //GWS GSM 
  //int volatile gstates[(MCOREPC)*MSBG]; //GWM MSG 
  THREADINFO *threads;//[MCOREPC];
#endif
  struct LocVar locv;
}threadGroup;

typedef struct _threadProc{
  short Nthreads,NSpecial,idMainThread,ngrp;
  short PThreadInited;
  int nstep;
  int mpi_id,mpi_npe;
  int volatile state;       //SWM MSS 
  THREADINFO tm;
#ifdef SOLIDINFO
  THREADINFO threads[MCOREPC];//for 2 level mthread
  threadGroup *grps[MCLUST];
  int gstate[(MCLUST)*MSB]; //MWG GSM 
#else
  THREADINFO *threads;//for 2 level mthread
  threadGroup **grps;//[MCLUST];
  int *gstate;//[(MCLUST)*MSB]; //MWG GSM 
#endif
  struct LocVar locv;
}threadProc;

__BEGIN_DECLS

extern threadProc md;
extern __thread THREADINFO *ti;
extern __thread threadGroup *gi;
//extern __thread void*td;
//extern __thread void*gd;
extern int ThreadG; //分组标志
//tread main  function
#define PTI HTHREADINFO ti
#define PTIM HTHREADINFO ti,
#define VTI ti
#define VTIM ti,
extern int _getgdsize_();
extern int _gettdsize_();
extern void _threadmain_(HTHREADINFO ti);
//tread ctrl 
int InitThreads(int mpi_id_, struct CPUINFO*ci);
void StartThreads();
void EndThreads();
void initthreads_(int *mpi_id_, struct CPUINFO*ci ,int *err);
void startthreads_();
void zStartEThread(TFunc tfun,void*para,int id);
enum {
  LOC_THREAD=0,LOC_GROUP,LOC_PROC
};
void SetLocV(int typ,int ind,void*p);
void *GetLocV(int typ,int ind);
void endthreads_();
void bindthread();
void bindthread_();
int bindcpu(int id);
void opentf();
int getgid_();
int gettid_();
int getnthreads_();

void tscinit();
void tscb(int id);
void tsce(int id);
void tsceb(int id);
void prtsc(const char*tag);

//sync function
#define PLINE ,int line
#define FPLINE ,int *line
//SG
void sWaitGrp(int state,int nt PLINE);
void sWaitGrpr(int state,int nt PLINE);
void sSetGrp(int state PLINE);

void gWaitSubs (int state,int nt PLINE);
void gWaitSubsr(int state,int nt PLINE);
void gSetSubs  (int state PLINE);
//GM
void gWaitMain(int state,int nt PLINE);
void gWaitMainr(int state,int nt PLINE);
void gSetMain(int state PLINE);

void mWaitGrps(int state,int nt PLINE);
void mWaitGrpsr(int state,int nt PLINE);
void mSetGrps(int state PLINE);
//SM 
void sWaitMain(int state,int nt PLINE);
void sWaitMainr(int state,int nt PLINE);
void sSetMain(int state PLINE);
// Fortran interface
//SG
void swaitgrp_  (int *state,int*nt FPLINE);
void swaitgrpr_ (int *state,int*nt FPLINE);
void ssetgrp_   (int *state FPLINE);
//GS
void gwaitsubs_ (int *state,int*nt FPLINE);
void gwaitsubsr_(int *state,int*nt FPLINE);
void gsetsubs_  (int *state FPLINE);
//GM
void gwaitmain_ (int *state,int*nt FPLINE);
void gwaitmainr_(int *state,int*nt FPLINE);
void gsetmain_  (int *state FPLINE);
//MS
void mwaitsubs_ (int *state,int*nt FPLINE);
void mwaitsubsr_(int *state,int*nt FPLINE);
void msetsubs_  (int *state FPLINE);
//MG
void mwaitgrps_ (int *state,int*nt FPLINE);
void mwaitgrpsr_(int *state,int*nt FPLINE);
void msetgrps_  (int *state FPLINE);

#define MLINE ,__LINE__
//SM
#define SSM(RFB)      sSetMain  (RFB MLINE) /*set sub  ready */
#define NSWM(RFB,nt)  sWaitMain (RFB,nt MLINE) /*wait grp ready */
#define NSWMR(RFB,nt) sWaitMainr(RFB,nt MLINE) /*wait grp ready */
//MS
#define MSS(RFB)      mSetSubs  (RFB MLINE) /*mmt */
#define NMWS(RFB,nt)  mWaitSubs (RFB,nt MLINE) /*mmt */
#define NMWSR(RFB,nt) mWaitSubsr(RFB,nt MLINE) /*mmt */
//SG
#define SSG(RFB)      sSetGrp   (RFB MLINE) /*set sub  ready */
#define NSWG(RFB,nt)  sWaitGrp  (RFB,nt MLINE) /*wait grp ready */
#define NSWGR(RFB,nt) sWaitGrpr (RFB,nt MLINE) /*wait grp ready */
//GS
#define GSS(RFB)      gSetSubs  (RFB MLINE) /*set gmthread ok*/
#define NGWS(RFB,nt)  gWaitSubs (RFB,nt MLINE) /*wait sub thread*/
#define NGWSR(RFB,nt) gWaitSubsr(RFB,nt MLINE) /*wait sub thread*/
//MG
#define MSG(RFB)      mSetGrps  (RFB MLINE) /*mmt */
#define NMWG(RFB,nt)  mWaitGrps (RFB,nt MLINE) /*mmt */
#define NMWGR(RFB,nt) mWaitGrpsr(RFB,nt MLINE) /*mmt */
//GM
#define GSM(RFB)      gSetMain  (RFB MLINE) /*set grp ok     */
#define NGWM(RFB,nt)  gWaitMain (RFB,nt MLINE) /*wait mmthread  */
#define NGWMR(RFB,nt) gWaitMainr(RFB,nt MLINE) /*wait mmthread  */
//SM
#define SSM(RFB)      sSetMain  (RFB MLINE) /*set sub  ready */
#define NSWM(RFB,nt)  sWaitMain (RFB,nt MLINE) /*wait grp ready */
#define NSWMR(RFB,nt) sWaitMainr(RFB,nt MLINE) /*wait grp ready */
//MS
#define MSS(RFB)      mSetSubs  (RFB MLINE) /*mmt */
#define NMWS(RFB,nt)  mWaitSubs (RFB,nt MLINE) /*mmt */
#define NMWSR(RFB,nt) mWaitSubsr(RFB,nt MLINE) /*mmt */
//SG
#define SSG(RFB)      sSetGrp   (RFB MLINE) /*set sub  ready */
#define NSWG(RFB,nt)  sWaitGrp  (RFB,nt MLINE) /*wait grp ready */
#define NSWGR(RFB,nt) sWaitGrpr (RFB,nt MLINE) /*wait grp ready */
//GS
#define GSS(RFB)      gSetSubs  (RFB MLINE) /*set gmthread ok*/
#define NGWS(RFB,nt)  gWaitSubs (RFB,nt MLINE) /*wait sub thread*/
#define NGWSR(RFB,nt) gWaitSubsr(RFB,nt MLINE) /*wait sub thread*/
//MG
#define MSG(RFB)      mSetGrps  (RFB MLINE) /*mmt */
#define NMWG(RFB,nt)  mWaitGrps (RFB,nt MLINE) /*mmt */
#define NMWGR(RFB,nt) mWaitGrpsr(RFB,nt MLINE) /*mmt */
//GM
#define GSM(RFB)      gSetMain  (RFB MLINE) /*set grp ok     */
#define NGWM(RFB,nt)  gWaitMain (RFB,nt MLINE) /*wait mmthread  */
#define NGWMR(RFB,nt) gWaitMainr(RFB,nt MLINE) /*wait mmthread  */


#define SWM(RFB)  NSWM(RFB,1)  
#define SWMR(RFB) NSWMR(RFB,1) 
#define MWS(RFB)  NMWS(RFB,1)  
#define MWSR(RFB) NMWSR(RFB,1) 
#define SWG(RFB)  NSWG(RFB,1)  
#define SWGR(RFB) NSWGR(RFB,1) 
#define GWS(RFB)  NGWS(RFB,1)  
#define GWSR(RFB) NGWSR(RFB,1) 
#define MWG(RFB)  NMWG(RFB,1)  
#define MWGR(RFB) NMWGR(RFB,1) 
#define GWM(RFB)  NGWM(RFB,1)  
#define GWMR(RFB) NGWMR(RFB,1) 
#define SWM(RFB)  NSWM(RFB,1)  
#define SWMR(RFB) NSWMR(RFB,1) 
#define MWS(RFB)  NMWS(RFB,1)  
#define MWSR(RFB) NMWSR(RFB,1) 
#define SWG(RFB)  NSWG(RFB,1)  
#define SWGR(RFB) NSWGR(RFB,1) 
#define GWS(RFB)  NGWS(RFB,1)  
#define GWSR(RFB) NGWSR(RFB,1) 
#define MWG(RFB)  NMWG(RFB,1)  
#define MWGR(RFB) NMWGR(RFB,1) 
#define GWM(RFB)  NGWM(RFB,1)  
#define GWMR(RFB) NGWMR(RFB,1) 

//#define DTIMES
#ifdef DTIMES
#define TBT(i)  tscb(i)
#define TET(i)  tsce(i)
#define TEBT(i) tsceb(i)
#define TPT(t)   prtsc(t)
#define TPH()   prth()
#else
#define TBT(i)
#define TET(i)
#define TEBT(i)
#define TPT(t)
#define TPH()
#endif
__END_DECLS
#endif
