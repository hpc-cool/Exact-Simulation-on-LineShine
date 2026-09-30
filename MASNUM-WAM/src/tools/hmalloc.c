#define _GNU_SOURCE
#include <sys/types.h>
#include <string.h>
#include <stdlib.h>

#include <stdio.h>

#include <malloc.h>

#include <errno.h>
//mmap
#include <unistd.h>
#include <asm/unistd.h>
#include <linux/mempolicy.h>
#include <sys/mman.h>

//dir
#include <dirent.h>
//bindcpu
#include <sched.h>
#include <pthread.h>
//shm 
#include <sys/ipc.h>
#include <sys/shm.h>

#include <hmalloc.h>

#ifndef MAXCPUS
#define MAXCPUS 4096
#endif
#ifndef MAXNODES 
#define MAXNODES 64
#define HBMBASE 16
#endif
#define BITPL (sizeof(long)*8)
#define MNL ((MAXNODES-1)/BITPL+1)
int mpi_id=0,mpi_npe=1,NThreads=1;
#ifndef MTHREAD
struct CPUINFO cpuinf={
  .OffClu=0,.SkipClu=1,.OffCore=0,.SkipCore=1,
  .NCorePClu=38,.MCorePClu=38,.NCluPNode =16,.NManageCore=4,
  .NThPGrp  =38,.NGrpPProc=16,.NProcPNode=1, .ManageCoreId=-1
};
#endif
struct memha{//4byte align，used size<1k，free size<4G
  unsigned int size;//BIT 0-1; xx11,used memha,xx01,free memha
                    //size=(size&(-8))
};
struct memhb{//16byte align，1KB<used size<64MB
  unsigned int sizeb;//BIT 0-3; 0010:used memhb,0000:free memb
  int fill[2];
  unsigned int size;//BIT 0-3; 0010:used memhb,0000:free memb
                    //size=(size&(-16))
};
struct heap{
  void*p;
  unsigned int size, used;
  int node,align;//0:4B,1:16B;2:64B,3:256B,4:1KB,5:4KB;
  struct heap *next;
};
struct hbmheap{
  struct heap h0,h1,h2,h3,h4,h5; 
};
struct hbmheap hbmheaps[MAXNODES]={0};

struct hmallocstate{
  void*p;
  int size;
  short node,pgtype;
};
struct hugepageinf{
  int hpsize,hpUsed, hpTotal,hpFree;
};
struct hbmstate{//32
  int lock;
  size_t memUsed,memTotal,memFree;
  struct hugepageinf hps[4];
};
struct hbmstates{//32+32*MAXNODES
  int jobid, lock;
  int nnodes,nodemasklen;
  int pagesize,pgbits,pgmask,pgmaskh;
  long nodemask[MNL];
  struct hbmstate hbms[MAXNODES];
};
#define HPAGESIZE 2*1024*1024
#ifndef MTHREAD
int bindcpu(int id){
  cpu_set_t mask;  //CPU核的集合
  CPU_ZERO(&mask);    //置空
  CPU_SET(id,&mask);   //设置亲和力值
  if (sched_setaffinity(0, sizeof(mask), &mask) == -1){//设置线程CPU亲和力
    printf("warning: could not set CPU affinity %d, %d ,continuing...\n",mpi_id,id);
    return -1;
  }
  return 0;
}
#endif
static void zunlock(int *p){ *p=0;}
static void zlock(int *p){ while(__sync_val_compare_and_swap(p,0,1))usleep(1);}
int init_cpu(int mode,MPI_Comm comm,int mpi_id_,int mpi_npe_,struct CPUINFO *ci){
  mpi_id      =mpi_id_;
  mpi_npe     =mpi_npe_;
  cpuinf=*ci;
  if(0&&mpi_id==0){
    printf("cpuset mpi_npe=%d,NGrpPProc=%d,NThPGrp=%d,NProcPNode=%d,ManageCoreId=%d\n",
           mpi_npe,cpuinf.NGrpPProc,cpuinf.NThPGrp,cpuinf.NProcPNode,cpuinf.ManageCoreId);
  }
  if(mode==0){
    int ind=mpi_id%cpuinf.NProcPNode;
    int ncpg=(cpuinf.NProcPNode+cpuinf.NCluPNode-1)/cpuinf.NCluPNode;
    int indg=(ind/ncpg)*cpuinf.NCorePClu + (ind%ncpg); 
    bindcpu(indg);
  }
  inithbms( 0,comm);
  return 0;
}
void init_cpu_(int *mode,int*icomm,int *mpi_id_,int *mpi_npe_,struct CPUINFO *ci,int *err){
  MPI_Comm comm=FCOMM2C(*icomm);
  *err=init_cpu(*mode,comm,*mpi_id_,*mpi_npe_,ci);
}
static struct hbmstates *hbms=NULL;
static struct hmallocstate*hms=NULL;
static int nhmalloc=0,mhmalloc=0;
static int lock=0;
static long getmeminfo(struct hbmstate*ph,int node,int prt){
  char fn[256];
  if(ph){ph->memUsed=ph->memTotal=ph->memFree=0;}
  if(node>=hbms->nnodes)return 0;
  sprintf(fn,"/sys/devices/system/node/node%d/meminfo",node);
  FILE*fi=fopen(fn,"rt");
  if(!fi)return 0;
  size_t total=0,free=0,used=-1,vt,vtt;
  int nr=3,nd;
  char s[64],tag[64]="",ad[64];
  while(!feof(fi)){
    ad[0]=0;
    fscanf(fi,"%s %d %s %ld %s\n",s,&nd,tag,&vt,ad);
    char c=(*ad)|0x20;
    vtt=vt;
    if(c=='k')vt*=1024;
    else if(c=='m')vt*=1024*1024;
    else if(c=='g')vt*=1024*1024*1024;
    if(strcmp(tag,"MemTotal:")==0){
      total=vt;nr--;
    }else if(strcmp(tag,"MemFree:")==0){
      free=vt;nr--;
    }else if(strcmp(tag,"MemUsed:")==0){
      used=vt;nr--;
    }
    if(prt)printf("A:%s %d %s %ld %s\n",s,nd,tag,vtt,ad);
    if(nr<=0)break;
  }
  if(used>total){
    if(total<free)used=40*1024*1024;
    else used=total-free;
  }
  if(ph){
    ph->memUsed=used;
    ph->memTotal=total;
    ph->memFree=free;
  }
  return 1;
}
static int gethugepageinfo(struct hugepageinf*hp,int node,int pgsize, int prt) {
  //      1GB ,    2MB ,    32MB
  char fn[256];
  int total=0,free=0,used=-1;
  if(hp){
    hp->hpUsed =hp->hpTotal=hp->hpFree =0;
    hp->hpsize=pgsize<<10;
  }
  // 2. 读取总大页数量
  sprintf(fn,"/sys/devices/system/node/node%d/hugepages/hugepages-%dkB/nr_hugepages", node,pgsize);
  FILE* fi = fopen(fn, "rt");
  if (fi){
    fscanf(fi, "%d",  &total);
    fclose(fi);
  }
  // 2. 读取空闲大页数量
  sprintf(fn,"/sys/devices/system/node/node%d/hugepages/hugepages-%dkB/free_hugepages", node,pgsize);
  fi = fopen(fn, "rt");
  if (fi){
    fscanf(fi, "%d", &free);
    fclose(fi);
  }
  if(prt)printf("HUGEPAGE:: %d %8d: %d :%d\n", node, pgsize,total, free);
  used=total-free;
  if(free>total){ free=total; used=0; }
  if(hp){
    hp->hpUsed =used;
    hp->hpTotal=total;
    hp->hpFree =free;
    hp->hpsize=pgsize<<10;
  }
  return 1;
}

static void setbit(long*nm,int n);
static int getNumNodes(long*nodemask){
  DIR *d;
  struct dirent*de;
  d=opendir("/sys/devices/system/node");
  if(!d)return 0;
  int md=0;
  while((de=readdir(d))){
    if(strncmp(de->d_name,"node",4))continue;
    int nd=atoi(de->d_name+4);
    if(nodemask)setbit(nodemask,nd);
    if(md<nd)md=nd;
  }
  closedir(d);
  return md+1;
}
static void *shm(int id,int size){
  int shmid,key;
  char*pshm;
  if((key=ftok("/bin",id))==-1){ perror("ftok");return NULL; }
  if(mpi_id % cpuinf.NProcPNode==0){
    if((shmid=shmget(key,size,IPC_CREAT|0666))==-1){ //|IPC_EXCL
      perror("shmgeta");return NULL; 
    }
    if((pshm=shmat(shmid,NULL,0))==(char*)-1){perror("shmat");return NULL;}
    memset(pshm,0,size);
  }else{
    if((shmid=shmget(key,size,0666))==-1){
      perror("shmgetb");return NULL; 
    }
    if((pshm=shmat(shmid,NULL,0))==(char*)-1){perror("shmat");return NULL;}
  }
  return pshm;
}
static void lockhbm(){
  if(hbms)zlock(&hbms->lock);
}
static void unlockhbm(){
  if(hbms)zunlock(&hbms->lock);
}
static void cleanup_hbms(){
  struct hmallocstate*pm=hms;
  if(pm){
    for(int i=0;i<mhmalloc;i++,pm++){
      munmap(pm->p,pm->size);
      pm->p=NULL;
      struct hbmstate *ph=&hbms->hbms[pm->node];
      lockhbm();
      nhmalloc--;
      ph->memUsed-=pm->size;
      unlockhbm();
    }
    free(hms);
  }
  shmdt(hbms);
}
static int hugepagesizes[4]={1048576 , 32768,2048 ,   64};// 30 25 21 16
#define nhugepagesizes sizeof(hugepagesizes)/sizeof(hugepagesizes[0])
//called by Main thread
void inithbms(int nnodes_,MPI_Comm comm){
  int onn=0;
  int pagesize=getpagesize();
  int ssize=sizeof(*hbms)-sizeof(hbms->hbms)+sizeof(hbms->hbms[0])*MAXNODES;
  int size=pagesize;
  for(;size<ssize;size+=pagesize);
  if(mpi_id % cpuinf.NProcPNode==0){
    if(!hbms){
      int nnodes=MAXNODES;
      if(nnodes_)nnodes=nnodes_;
      hbms=shm(41528,size);
      if(!hbms){
        printf("hbms Init err 1 %d %p\n",mpi_id,hbms);
        exit(-1);
      }
      atexit(cleanup_hbms);
      zlock(&hbms->lock);
      if(hbms->pagesize==0){
        hbms->pagesize=0;
        hbms->pgmask=pagesize-1; 
        hbms->pgmaskh=-pagesize;
        for(int i=0,ps=pagesize>>1;;i++,ps>>=1){ 
          if(!ps){ hbms->pgbits=i;break;} 
        }
        hbms->nnodes=nnodes;
        memset(hbms->hbms,0,sizeof(*hbms->hbms)*MAXNODES);
        hbms->nodemasklen=(hbms->nnodes-1)/BITPL+1;
        //printf("HBMF: %d %d %d :\n",mpi_id,ti->igrp,ti->ind);
        //fflush(stdout);
        for(int i=0;i<hbms->nnodes;i++){
          struct hbmstate*pm=hbms->hbms+i;
          if(getmeminfo(pm,i,0)==0)continue;
          for(int j=0;j<nhugepagesizes;j++){
            gethugepageinfo(pm->hps+j,i,hugepagesizes[j],0);
          }
          if(pm->memTotal==0) {hbms->nnodes=i;break;}
          //if(i>=16)printf("HBMG:%d %d/%d: %ld %ld %ld\n",mpi_id,i,hbms->nnodes,pm->memFree/1024,pm->memTotal/1024,pm->memUsed/1024);
#if 0
          if(pm->memTotal-pm->memUsed<512*1024*1024){
            printf("Hbm error node:%d %d total:%ldKB Used:%ldKB\n",mpi_id,i,pm->memTotal/1024,pm->memUsed/1024);
            exit(-10);
          }
#endif
        }
        hbms->pagesize=pagesize;
      }
      zunlock(&hbms->lock);
    }
  }
#ifndef NO_MPI
  MPI_Barrier(comm);
  if(mpi_id % cpuinf.NProcPNode>0){
    hbms=shm(41528,size);
    atexit(cleanup_hbms);
    if(!hbms){
      printf("hbms Init err 2 %d \n",mpi_id);
      exit(-1);
    }
  }
#endif
}
int getcpus(){
  return sysconf(_SC_NPROCESSORS_CONF);

}
int getcpuid(){
  cpu_set_t mask;  //CPU核的集合
  CPU_ZERO(&mask);    //置空
  if (sched_getaffinity(0, sizeof(mask), &mask) == -1){//设置线程CPU亲和力
    printf("warning: could not get CPU affinity , continuing...\n");
    return -1;
  }
  for(int i=0;i<CPU_SETSIZE;i++){
    if(CPU_ISSET(i,&mask)){
      return i;
    }
  }
  return 0;
}
static void clearbits(long*nm){
  for(int i=0;i<MNL;i++)nm[i]=0;
}
static void setbit(long*nm,int n){
  if(n<hbms->nnodes){ int i=n/BITPL; nm[i]|=1UL<<(n%BITPL); }
}
static void clearbit(long*nm,int n){
  if(n<hbms->nnodes){ int i=n/BITPL; nm[i]&=~(1UL<<(n%BITPL)); }
}
static void setnode(long*nm,int n){
  int it=n/BITPL; 
  for(int i=0;i<hbms->nodemasklen;i++){
    if(i!=it) nm[i]=0;
    else nm[i]=1UL<<(n%BITPL); 
  }
}
void zunlock_();
void zlock_();

void print_call_stack(int lock) ;
void*hmalloc_onnode(size_t size,int nodeid,int bind) {
  void*p=NULL;
  size_t sizea=0;
  if(!hbms){
    printf(" HBM Not Init %d \n",mpi_id);
    exit(-2);
  }
  struct hbmstate *ph=&hbms->hbms[nodeid];
  int ipgs=hbms->pagesize/sizeof(int);
  int isize=size/sizeof(int);
  int pgtype=0;
  zlock(&ph->lock);
#if 0
  printf("HM:%d %2.2d %2.2d %2.2d %ldKB\n",mpi_id,bind,nodeid,getcpuid(),size>>10);
  if(mpi_id%2){
    if(bind&&nodeid<24){
      //print_call_stack(1); 
    }
  }
#endif
  for(int j=0;j<nhugepagesizes;j++){
    struct hugepageinf *hp=ph->hps+j;
    int hpsize=hp->hpsize;
    if(size>hpsize){
      int npgs=(size-1)/hpsize+1;
      if(npgs<hp->hpFree){
        sizea = (size+hpsize-1)&(-hpsize);
        p = mmap(NULL, sizea, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS|MAP_HUGETLB, -1, 0);
        if (p == MAP_FAILED) p=NULL;
        else {
          pgtype=j+1; hp->hpUsed+=npgs;hp->hpFree-=npgs;
          //printf("mapped hupgepage %d %d %d %d %ldKB %ldKB\n",mpi_id,nodeid,hpsize,npgs,size>>10,sizea>>10);
          break;
        }
      }
    }
  }
  if(!p){
    if(!bind){
      zunlock(&ph->lock);
      return malloc(size);
    }
    if(ph->memUsed+size>ph->memTotal){//fixme
      printf("HBM warning:nospace mpi_id:%d %d %ldMB %ldMB %ldMB\n",mpi_id,nodeid,ph->memUsed>>20,ph->memFree>>20,size>>20);
      zunlock(&ph->lock);
      return malloc(size);
      //getmeminfo(0,nodeid,1);
      for(int j=0;j<nhugepagesizes;j++){
        gethugepageinfo(NULL,nodeid,hugepagesizes[j],0);
      }
      return NULL;
    }
    sizea = (size+hbms->pagesize-1)&hbms->pgmaskh;
    p = mmap(NULL, sizea, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (p == MAP_FAILED) {
      printf("mmap failed  on numa-%d, length=%ld. \n", nodeid, sizea); return NULL;
    }
  }
  //if(sizea>hbms->hpsize) madvise(p, sizea, MADV_HUGEPAGE);
  if(bind){
    long nm[MNL];
    setnode(nm, nodeid);
    //int ret = mbind(p, sizea, MPOL_BIND, nm, MAXNODES, 0);
    int ret=syscall(__NR_mbind,p, sizea, MPOL_BIND, nm, MAXNODES, 0);
    if (ret != 0) {
      printf( "mbind failed  on numa-%d. \n", nodeid );
      munmap(p,sizea);
      zunlock(&ph->lock);
      return NULL;
    }
  }

  if(pgtype==0) {
    ph->memUsed+=sizea;
    ph->memFree-=sizea;
  }
  //printf("mmap ok %d, %ldK %ldK %ldM %ldM %ldM.end \n", nodeid, sizea>>10,size>>10,ph->memUsed>>20,ph->memTotal>>20,ph->memFree>>20); fflush(stdout);
  zunlock(&ph->lock);
  zlock(&lock);
#if 0
  zlock_();
  printf("mmap ok %d, %ld %ld %ldM %ldM %ldM.checking... \n", nodeid, sizea,size,ph->memUsed>>20,ph->memTotal>>20,ph->memFree>>20); fflush(stdout);
  for(int i=0;i<isize;i+=ipgs){ *(((int*)p)+i)=0; }
  printf("mmap ok %d, %ld %ld %ldM %ldM %ldM.check end \n", nodeid, sizea,size,ph->memUsed>>20,ph->memTotal>>20,ph->memFree>>20); fflush(stdout);
  zunlock_();
#endif
  if(mhmalloc<=nhmalloc){
    mhmalloc+=16;
    hms=realloc(hms,sizeof(*hms)*mhmalloc);
    memset(hms+mhmalloc-16,0,16*sizeof(*hms));
  }
  struct hmallocstate*pm=hms;
  for(int i=0;i<mhmalloc;i++,pm++)if(pm->p==NULL)break;
  pm->p=p;pm->size=sizea;pm->node=nodeid;pm->pgtype=pgtype;nhmalloc++;
  zunlock(&lock);
  return p;
}
void*hmalloc(size_t size){
  int cpuid=getcpuid();
  int node=cpuid/cpuinf.NCorePClu+HBMBASE;
  return hmalloc_onnode(size,node,1);
}
void*hpmalloc(size_t size,int bind){
  if(size<4096&&bind==0)return malloc(size);
  int cpuid=getcpuid();
  int node=cpuid/cpuinf.NCorePClu;
  return hmalloc_onnode(size,node,bind);
}
void*hrealloc(void*p,size_t size){
  if(!p) return NULL;
  struct hmallocstate*pm=hms;
  for(int i=0;i<mhmalloc;i++,pm++){
    if(pm->p==p){
      struct hbmstate *ph=&hbms->hbms[pm->node];
      size_t sizea = (size+hbms->pagesize-1)&hbms->pgmaskh;
      void*p=mremap(pm->p,pm->size,sizea,MREMAP_MAYMOVE);
      if(p==MAP_FAILED)return NULL;
      pm->p=p;pm->size=sizea;
      zlock(&ph->lock);
      ph->memUsed-=sizea-pm->size;
      zunlock(&ph->lock);
      return p;
    }
  }
  return NULL;
}
void*hvrealloc(void*p,size_t size){
  if(!p)return malloc(size);
  void *pn=hrealloc(p,size);
  if(pn)return pn;
  return realloc(p,size);
}
void*hvmalloc_onnode(size_t size, const int nodeid,int bind){
  void *p=hmalloc_onnode(size,nodeid,bind);
  if(p)return p;
  return malloc(size);
}
void*hvmalloc(size_t size){
  void *p=hmalloc(size);
  if(p)return p;
  return malloc(size);
}
int hfree(void*p){
  if(!p)return 1;
  struct hmallocstate*pm=hms;
  for(int i=0;i<mhmalloc;i++,pm++){
    if(pm->p==p){
      munmap(pm->p,pm->size);
      struct hbmstate *ph=&hbms->hbms[pm->node];
      zlock(&ph->lock);
      ph->memUsed-=pm->size;
      zunlock(&ph->lock);
      zlock(&lock);
      pm->p=NULL; nhmalloc--;
      zunlock(&lock);
      return 1;
    }
  }
  return 0;
}
void hvfree(void*p){
  if(p){
    if(!hfree(p))free(p);
  }
}
void inithbms_(int *nnodes_,int *icomm){ 
  MPI_Comm comm=FCOMM2C(*icomm);
  inithbms(*nnodes_,comm); 
}
void*hmalloc_onnode_(size_t*size,int*nodeid,int *bind){
return hmalloc_onnode(*size,*nodeid,bind?*bind:1);
}
void*hvmalloc_onnode_(size_t*size,int*nodeid,int *bind){
return hvmalloc_onnode(*size,*nodeid,bind?*bind:1);
}
void*hmalloc_(size_t *size){ return hmalloc(*size); }
void*hrealloc_(void*p,size_t *size){ return hrealloc(p,*size);}
void*hvmalloc_(size_t *size){ return hvmalloc(*size); }
void*hvrealloc_(void*p,size_t *size){ return hvrealloc(p,*size);}
void hfree_(void*p){ hfree(p); }
void hvfree_(void*p){ hvfree(p); }
#ifdef TESTHMALLOC
int main(int ac,char**av){
  int size=1024*1024*1024;
  void*p=hmalloc_onnode(size,17);
  printf("hmalloc %p size 0x%x %8.2fM\n",p,size,size/1024./1024.);
  for(int i=0;i<10;i++){
    printf("\r %4d ",i);fflush(stdout);
    sleep(1);
  }
  printf("\n");
  hfree(p);
  int s=4096;
  printf("%x %x %x %x\n",s,!s,-s,~s);

}
#endif
#if 0
#ifdef USEMEMKIND
void*hmalloc(size_t size){
  int cpuid=getcpuid();
  lockhbm();
  struct hbmstate *ph=&hbms->hbms[cpuid/NCorePClu];
  void*pr=NULL;
  if(ph->HbmMax==0){
    ph->HbmMax=0xF0000000L;
    ph->hbmused=0;
  }
  if(ph->hbmused+size>ph->HbmMax){//fixme
    printf("HBM warning:nospace mpi_id:%d %lX %lX\n",mpi_id,ph->hbmused,size);
    pr=malloc(size);
  }else{
    pr=memkind_malloc(MEMKIND_HBW,size);
    if(!pr){
      printf("HBM warning:alloc failed pi_id:%d %d %lx %lx\n",mpi_id,cpuid/NCorePClu,ph->hbmused,size);
      pr=malloc(size);
    }else{
      size_t st=memkind_malloc_usable_size(MEMKIND_HBW, pr);
      ph->hbmused+=st;
    }
  }
  unlockhbm();
  return pr;
}
void hfree(void*p){
  if(p){
    if(hbw_verify_memory_region(p,16,HBW_TOUCH_PAGES)==1)
      hbw_free(p);
    else free(p);
  }
}
void*hrealloc(void*p,size_t size){
  void*pr=hmalloc(size);
  if(p&&pr){
    memcpy(pr,p,size);
    hfree(p);
  }
  return pr;
}
#elif defined(USEKUPL)
void hfree(void*p){
  if(p){
    if(kupl_hbw_verify(p,16,KUPL_HBW_TOUCH_PAGES)==1)
      kupl_hbw_free(p);
    else free(p);
  }
}
void*hmalloc(size_t size){
  void *pr=kupl_hbw_malloc(size);
  if(!pr){
    printf("HBM warning:alloc failed pi_id:%d %lx\n",mpi_id,size);
    pr=malloc(size);
  }else{
    printf("HBM INFO:alloc OK pi_id:%d %lx\n",mpi_id,size);
  }
  return pr;
}
void*hrealloc(void*p,size_t size){
  void *pr=kupl_hbw_malloc(size);
  if(!pr){
    printf("HBM warning:alloc failed pi_id:%d %lx\n",mpi_id,size);
    if(p&&kupl_hbw_verify(p,16,KUPL_HBW_TOUCH_PAGES)!=1)
      return realloc(p,size);
    pr=malloc(size);
  }
  if(p&&pr){
    memcpy(pr,p,size);
    hfree(p);
  }
  return pr;
}
#endif
#endif
