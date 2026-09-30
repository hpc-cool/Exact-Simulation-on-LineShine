#ifndef HMALLOC_H_INCLUDED
#define HMALLOC_H_INCLUDED
#include <sys/types.h>
__BEGIN_DECLS
struct CPUINFO;
#ifdef MMTHREAD
#ifndef MTHREAD
#define MTHREAD
#endif
#endif

#ifndef MTHREAD
extern int mpi_id,NThreads;
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
extern struct CPUINFO cpuinf;
#else
#include "mthread.h"
#endif
#ifndef NO_MPI
#include <mpi.h>
#define FCOMM2C(comm)  MPI_Comm_f2c(comm)
#else
#define FCOMM2C(comm)  (comm)
#define MPI_Comm  int
#endif
int init_cpu(int mode,MPI_Comm comm,int mpi_id_,int mpi_npe_,struct CPUINFO *ci);
int bindcpu(int id);
void inithbms(int nnodes_,MPI_Comm comm);
void*hmalloc_onnode(size_t size,int nodeid,int bind);
void*hmalloc(size_t size);
void*hpmalloc(size_t size,int bind);
void*hrealloc(void*p,size_t size);
int hfree(void*p);

void*hvmalloc_onnode(size_t size,int nodeid,int bind);
void*hvmalloc(size_t size);
void*hvrealloc(void*p,size_t size);
void hvfree(void*p);

void init_cpu_(int *mode,int *icomm,int *mpi_id_,int *mpi_npe_,struct CPUINFO *ci,int *err);
void inithbms_(int *nnodes_,int *icomm);
void*hmalloc_onnode_(size_t*size,int*nodeid,int *bind);
void*hmalloc_(size_t *size);
void*hrealloc_(void*p,size_t *size);
void hfree_(void*p);
void*hvmalloc_onnode_(size_t*size,int*nodeid,int *bind);
void*hvmalloc_(size_t *size);
void*hvrealloc_(void*p,size_t *size);
void hvfree_(void*p);
__END_DECLS
#endif
