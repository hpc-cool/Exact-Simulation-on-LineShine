#ifndef MTHREAD_H_INCLUDED
#include <pthread.h>
#include <stdint.h>
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
#define MSBG    16
#endif
typedef struct _THREADINFO THREADINFO,*HTHREADINFO;
typedef void (*TFunc)(HTHREADINFO);
typedef void *(*TSFunc)(void*);
void initmd();
typedef struct _THREADINFO{
  short Nthreads,MainThread,igrp;      // 4 * 2 bytes
  short nlocv,mlocv,detach,indg,sib,sie;   // 3 * 2 bytes
  int ind;
  int sync_flag;
  int nstep;
  int volatile state;     //MWS SSM //GWS SSG
  uint64_t tscbs[16];
  uint64_t tscds[16];
  struct _threadGroup *pg;
  struct _THREADINFO  *threads;
  TFunc tfun;
  void *para;
  pthread_t threadid;
  FILE *fo;
  void*td;
  void*locv[8];
}THREADINFO;
typedef struct _threadGroup{
  short igrp,ngrp,mainGroup,Nthreads,idMainThread,NSpecial;
  short nlocv,mlocv,PThreadInited;
  int nstep;
  int volatile state;         //SWG GSS 
  int volatile gstate;        //GWM MSG 
  struct _threadProc *pmd;
  int volatile sstate[(MCOREPC)*MSBG]; //GWS GSM 
  THREADINFO threads[MCOREPC];
  void*gd;
  void*locv[8];
}threadGroup;

typedef struct _threadProc{
  short Nthreads,NSpecial,idMainThread,ngrp;
  short nlocv,mlocv,PThreadInited;
  int nstep;
  int mpi_id,mpi_npe;
  THREADINFO tm;
  THREADINFO threads[MCOREPC];
  threadGroup *grps[MCLUST];
  int volatile state;       //SWM MSS 
  int gstate[(MCLUST)*MSB]; //MWG GSM 
  void*locv[8];
}threadProc;

__BEGIN_DECLS
extern threadProc md;
extern int mpi_id, NCorePClu,NCluPNode, NCorePGrp,NGrpPProc,NThPGrp,NProcPNode,ManageCoreId,NThreads;
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
//extern int _getgdsize_();
//extern int _gettdsize_();
#ifndef USEHBM 
#define hmalloc malloc
#define hrealloc realloc
#else
extern void*hmalloc(size_t size);
extern void*hrealloc(void*p,size_t size);
#endif
extern void _threadmain_(HTHREADINFO ti);



//tread ctrl 
void setthread(int NCorePClu_ ,int NThPClu_,int NGrpPProc_,int NProcPNode_,int ManageCoreId_);
int InitThreads(int mpi_id_,int NCorePClu_ ,int NCluPNode_,int NCorePProc_,int NThPGrp_,int NGrpPProc_,int NProcPNode_,int *ManageCoreId_);
void StartThreads();
void EndThreads();
void inithreads_(int *mpi_id_,int *NCorePClu_ ,int *NCluPNode_,int *NCorePProc_,int *NThPGrp_,int *NGrpPProc_,int *NProcPNode_,int *ManageCoreId_,int *err);
void startthreads_();
void endthreads_();
void bindthread();
void bindthread_();
int bindcpu(int id);
void opentf();
void ntdelay_(int n);

void tscinit();
void tscb(int id);
void tsce(int id);
void tsceb(int id);
void prtsc(const char*tag);

//sync function
//SG
void sWaitGrp(int state);
void sWaitGrpr(int state);
void sSetGrp(int state);

void gWaitSubs (int state);
void gWaitSubsr(int state);
void gSetSubs  (int state);
//GM
void gWaitMain(int state);
void gWaitMainr(int state);
void gSetMain(int state);

void mWaitGrps(int state);
void mWaitGrpsr(int state);
void mSetGrps(int state);
//SM 
void sSetState(int state);
void sWaitState(int state);
void sWaitStater(int state);
void mSetSubs(int state);
void mWaitSubs(int state);
void mWaitSubsr(int state);
// Fortran interface
//SG
void swaitgrp_   (int *state);
void swaitgrpr_  (int *state);
void ssetgrp_    (int *state);
//k
//m
//GS
void gwaitsubs_  (int *state);
void gwaitsubsr_ (int *state);
void gsetsubs_   (int *state);
//GM
void gwaitmain_  (int *state);
void gwaitmainr_ (int *state);
void gsetmain_   (int *state);
//MS
void mwaitsubs_  (int *state);
void mwaitsubsr_ (int *state);
void msetsubs_   (int *state);
//MG
void mwaitgrps_ (int *state);
void mwaitgrpsr_(int *state);
void msetgrps_  (int *state);

//SM
#define SSM  sSetState   (RFB) /*set sub  ready */
#define SWM  sWaitState  (RFB) /*wait grp ready */
#define SWMR sWaitStater (RFB) /*wait grp ready */
//MS
#define MSS  mSetSubs    (RFB) /*mmt */
#define MWS  mWaitSubs   (RFB) /*mmt */
#define MWSR mWaitSubsr  (RFB) /*mmt */
//SG
#define SSG  sSetGrp     (RFB) /*set sub  ready */
#define SWG  sWaitGrp    (RFB) /*wait grp ready */
#define SWGR sWaitGrpr   (RFB) /*wait grp ready */
//GS
#define GSS  gSetSubs    (RFB) /*set gmthread ok*/
#define GWS  gWaitSubs   (RFB) /*wait sub thread*/
#define GWSR gWaitSubsr  (RFB) /*wait sub thread*/
//MG
#define MSG  mSetGrps    (RFB) /*mmt */
#define MWG  mWaitGrps   (RFB) /*mmt */
#define MWGR mWaitGrpsr  (RFB) /*mmt */
//GM
#define GSM  gSetMain    (RFB) /*set grp ok     */
#define GWM  gWaitMain   (RFB) /*wait mmthread  */
#define GWMR gWaitMainr  (RFB) /*wait mmthread  */

__END_DECLS
#endif