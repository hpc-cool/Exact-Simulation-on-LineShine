/*
	some tools write in c code,used by fortran code
*/
#define _GNU_SOURCE
#include <stdio.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#define DBGPINT(a) //printf a
#define PENV(name) printf("ENV %s=%s\n",#name,getenv(#name))
extern char **environ;
#define MAXTIMER 30
#ifdef MTHREAD
#define __THREADV __thread 
#else
#define __THREADV 
#endif
static __THREADV double rtcb[MAXTIMER+1],rtcl[MAXTIMER+1];
static __THREADV long dclkb=0;
int iwalltime() {
  time_t t=time(NULL);
  struct tm *tm;
  tm=localtime(&t);
  //tm->tm_year*10000+tm->tm_mon*100+tm->tm_nday;
  return tm->tm_hour*10000+tm->tm_min*100+tm->tm_sec;
}
int iwalltime_() { return iwalltime() ; }
double dclktime(){
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME,&ts);
  if(dclkb==0)dclkb=ts.tv_sec;
  return (ts.tv_sec-dclkb)+ts.tv_nsec*1.e-9;
}
double dclktime_(){ return dclktime(); }
void  resettimer(int iid){
  double rtime=dclktime();
  if( iid>=0){
    rtcl[iid]=0;rtcb[iid]=rtime;
  }else{
    for(int i=0;i<=MAXTIMER;i++){
      rtcl[i]=0;rtcb[i]=rtime;
    }
  }
}
void  resettimer_(int *iid){ resettimer(iid?*iid:-1); }
void endstarttimer(int iide,int iids){
  double rtime=dclktime();
  if(iide>=0)rtcl[iide]+=rtime-rtcb[iide];
  if(iids>=0)rtcb[iids] =rtime;
}
void endstarttimer_(int *iide,int *iids){ endstarttimer(*iide,*iids); }
void starttimer(int iid){ if(iid>=0)rtcb[iid]=dclktime(); }
void starttimer_(int *iid){ if(iid&&*iid>=0)rtcb[*iid]=dclktime(); }
void endtimer(int iid){ if(iid>=0)rtcl[iid]+=(dclktime()-rtcb[iid]); }
void endtimer_(int *iid){ if(iid&&*iid>=0)rtcl[*iid]+=(dclktime()-rtcb[*iid]); }
double difftimer(int iid){ return (dclktime()-rtcb[iid]); /*+rtcl(iid)*/ }
double difftimer_(int *iid){ return (dclktime()-rtcb[*iid]); /*+rtcl(iid)*/ }
double gettimer(int iid){ return rtcl[iid]; }
double gettimer_(int *iid){ return rtcl[*iid]; }
void msleep(int ms){ usleep(ms*1000); }
void msleep_(int *ms){ usleep(*ms*1000); }
void checkenv_(){
	char *s;
	int i;
	PENV(TM_API);
}
void setmask1_(int *p,int *nn){
	int i;
	for(i=0;i<*nn;i++)p[i]=p[i]>0?1:0;
}
void setgrid_(int *p,int *nn,float*pvb,float *pv0,float *pv1){
	int i;
	float vb,v0,v1,v,vz;
	vb=*pvb;v0=*pv0;v1=*pv1;vz=0.;
	//printf("%p %p %d\n",p,p+1,*nn)
	for(i=0;i<*nn;i++){
		v=p[i];
		p[i]=v>vb?(v<v0?v0:(v>v1?v1:v)):vz;
	}
}
void countrect_(int*p,int *pm,int *pim,int*pjm,int *prect,int*pnc){
	int i,j;
	int im,jm,nc,m,i1,i2,j1,j2,i0,j0;
	im=*pim;jm=*pjm;m=*pm;nc=0;
	i1=i2=j1=j2=-1;
	for(j=1;j<=jm;j++){
		for(i=1;i<=im;i++,p++){
			if(*p!=m)continue;
			i0=i1=i2=i;j0=j1=j2=j;			
			nc++;
			break;
		}
		if(i1>0)break;
	}
	for(j=j0;j<=jm;j++){
		for(i=i0;i<=im;i++,p++){
			if(*p!=m)continue;
			nc++;
			i1=i1<=i?i1:i;
			i2=i2>=i?i2:i;
			j1=j1<=j?j1:j;
			j2=j2>=j?j2:j;
		}
		i0=1;
	}
	prect[0]=i1;	prect[1]=j1;	prect[2]=i2;	prect[3]=j2;
	*pnc=nc;
}
/*
	fortran interface :
	return int value at given pointer,
	for volatile variable,avoid optimize
*/
int IntVolatile(void*p){
	return *(int *)p;
}
int intvolatile_(void*p){
	return *(int *)p;
}

//*====================================================================*//
#ifdef CP_SWF90
/*  in SW set cpu mode,for underflow
*/
union long_double {
        unsigned long long ldata;        double ddata;
};
void set_fpcr(unsigned long long a){
        union long_double fpcr;fpcr.ldata = a; asm("wfpcr %0"::"f"(fpcr.ddata));
}

unsigned long long get_fpcr(){
        union long_double fpcr; asm("rfpcr %0":"=f"(fpcr.ddata));return fpcr.ldata;
}
void syslimit_(long value){
        struct rlimit rlim;
        if (value == 0) rlim.rlim_cur = RLIM_INFINITY;
        else    rlim.rlim_cur = value;
        rlim.rlim_max = RLIM_INFINITY;
        if (setrlimit(RLIMIT_STACK, &rlim) == -1) {
                perror("setrlimit A:");
                exit(-1);
        }
        rlim.rlim_cur = RLIM_INFINITY;
        if (setrlimit(RLIMIT_CORE, &rlim) == -1) {
                perror("setrlimit B:");
                exit(-1);
        }
        return;
}

void set_rndmode_(){
        union long_double fpcr; fpcr.ldata = get_fpcr();fpcr.ldata |= (0x1ULL<<48)   ;set_fpcr(fpcr.ldata);
}
void reset_rndmode_(){
        union long_double fpcr; fpcr.ldata = get_fpcr();fpcr.ldata &= (~(0x1ULL<<48));set_fpcr(fpcr.ldata);
}

void initcpu_sw_(){  syslimit_(0);  set_rndmode_();}
void initcpu_(){  syslimit_(0);  set_rndmode_();}

void rpcc_(unsigned long *c)
{
  unsigned long a;
  asm("rtc %0": "=r" (a) : );
  *c=a;
};

#elif defined( INTEL)
void initcpu_(){  }
void rpcc_(unsigned long *c)
{
unsigned long a, d;
 __asm__ __volatile__ ("rdtsc" : "=a" (a), "=d" (d));
*c = a | (d << 32);
};

#else
void initcpu_(){  }
void rpcc_(unsigned long *c)
{
unsigned long a, d;
printf("unknown CPU not surport rpcc\not");
*c=0;
};
#endif
void meminf(int mpiid){
	FILE*fi;
	char fn[256];
	long rss;
	long pid=getpid();
	sprintf(fn,"/proc/%ld/status",pid);
	fi=fopen(fn,"rt");
	rss=0;
	while(!feof(fi)){
		char*ir=fgets(fn,256,fi);
		if(strncmp(fn,"VmRSS",5)==0){
			//printf("MEMINFA: %d %d %s AAA\n",mpiid,pid,fn);
			rss=atol(fn+6);break;
		}
		//if(strncmp(fn,"VMRSS")==0)rss=atol(fn+6);
	}
	fclose(fi);
	printf("MEMINFB: %d %ld %ld kb\n",mpiid,pid,rss);
}
void meminf_(int *mpiid){
	meminf(*mpiid);
}

#include <stdarg.h>
#define ferrrpt stderr
#define finfrpt stderr


#ifndef finfrpt
INMIC FILE* finfrpt=NULL;
static void vopen(FILE**f,char*fnt){
	char fn[128];
	if(*f)return;
	sprintf(fn,"%s_%d.txt",fnt,md.si.selfid);	
	*f=fopen(fn,"wt");
}
#endif
//#define MID md.si.selfid
#define MID 0
void errrpt(const char*s,int id,const char *fmt,...) {
	va_list args;
	int err;
	err=errno;
  va_start(args, fmt);
  fflush(stdout);
  fprintf(ferrrpt," %s %d %d",s,id,MID );
  vfprintf(ferrrpt,fmt, args);
  fprintf(ferrrpt," %d %s\n",err,strerror(err));
  va_end(args);	
  fflush(ferrrpt);
}

void infrpt(const char*s,int id,const char *fmt,...) {
	va_list args;
  va_start(args, fmt);
  fflush(stdout);
#ifndef finfrpt
  vopen(&ferrrpt,"Infrpt");
#endif
  fprintf(finfrpt," %s %d %d",s,id,MID);
  vfprintf(finfrpt,fmt, args);
  fprintf(finfrpt,"\n");
  va_end(args);	
  fflush(finfrpt);
}
#include <math.h>
int ccheckee(double *ee,int n){
  int i ,err=0;
  for(i=0;i<n;i++){
    if(isnan(ee[i])||isinf(ee[i])){
      printf("ee err:%d %lf \n",i,ee[i]);
      err=i+1;
      return err;
    }
  }
  return 0;
}
void ccheckee_(double *ee,int *n,int *err){
  *err=ccheckee(ee,*n);
}

#include <mpi.h>
void wav_abort(){
#ifdef NO_MPI
  exit(0);
#else
  MPI_Abort(MPI_COMM_WORLD,0);
#endif

}

#include <execinfo.h>
// 解析地址为「函数名+行号+文件名」
static char*exec_cmd_get_output(const char *cmd, char *buf, size_t buf_size) {
	buf[0]=0;
	if (!(cmd == NULL || buf== NULL || buf_size == 0) ){
		FILE *fp = popen(cmd, "r");
		if (fp ) {
			int len=fread(buf, 1, buf_size - 1, fp);
			pclose(fp);
      buf[len]=0;
		}
	}
	return buf;
}
extern int mpi_id,mpi_npe,NThreads;
static int iolock=0;
static int psid=0;
void zunlock(int *p){ 
  *p=0;
}
void zlock(int *p){ 
  while(__sync_val_compare_and_swap(p,0,1))usleep(1);
}
void zunlock_(){ 
  iolock=0;
}
void zlock_(){ 
  while(__sync_val_compare_and_swap(&iolock,0,1))usleep(1);
}
void print_call_stack(int lock) {
	char prgpath[256],buf[256], cmd[512];
	void *callstack[100];
	int frame_num = backtrace(callstack, 100);
	ssize_t len = readlink("/proc/self/exe", prgpath,sizeof(prgpath) - 1);
	prgpath[len]=0;
  //if(frame_num>4)frame_num=4;
  if(lock)zlock_();
	//printf("%2d %3.3d Frames : %d \n",mpi_id,psid, frame_num);
	for (int i = 1; i < frame_num; i++) {
		snprintf(cmd, sizeof(cmd), "addr2line -e %s -f -C %p", prgpath, callstack[i]);
    exec_cmd_get_output(cmd,buf,256);
    for(char*p=buf;*p;p++){if(*p=='\r'||*p=='\n')*p=' ';}
		printf("%2d %3.3d Frame %d:%p %s\n",mpi_id,psid, i-1,callstack[i],buf);
	}
  psid++;
  if(lock)zunlock_();
  msleep(1);
}

void print_call_stack_(int *lock) {
  print_call_stack(*lock);
}
#include <sys/mman.h>

int is_ptr_valid(const void *ptr, size_t size) {
    if (ptr == NULL || size == 0) {
        errno = EINVAL;
        return -1;
    }
    unsigned char vec;
    void *page_aligned_ptr = (void *)((unsigned long)ptr & ~(sysconf(_SC_PAGESIZE) - 1));
    int ret = mincore(page_aligned_ptr, size, &vec);
    if (ret == 0) return 1;
    if (errno == EFAULT || errno == ENOMEM) return 0;
    return -1;
}
int is_ptr_valid_(const void *ptr, int* size) {
  return is_ptr_valid(ptr, *size) ;
}
#define IS_PTR_VALID(ptr) is_ptr_valid_(ptr, 1)

#include <signal.h>
/**
 * @brief 段错误信号处理函数
 * @param sig 信号值（SIGSEGV=11）
 * @param info 信号详情（包含错误地址）
 * @param ctx 上下文（寄存器信息）
 */
void perr();
static void segfault_handler(int sig, siginfo_t *info, void *ctx) {
  zlock_();
  printf("===== 段错误捕获 =====\n");
  //printf("threadinfo：%d %d %d %d\n", mpi_id,ti->igrp,ti->ind,ti->nstep);
  printf("信号：%d (SIGSEGV)\n", sig);
  printf("错误地址：%p\n", info->si_addr); // 触发段错误的指针地址
  printf("错误原因：%s\n", info->si_code == SEGV_MAPERR ? "地址未映射" : "权限不足");
  //perr();
  print_call_stack(0);
  fflush(stdout);
  zunlock_();
  // 可选：记录日志、释放资源、优雅退出
  usleep(10*1000*1000);
  exit(EXIT_FAILURE); // 必须退出，否则会无限触发信号
}
#include <fenv.h>   // 浮点异常头文件
// 浮点异常处理函数
static void sigfpe_handler(int sig, siginfo_t *info, void *ucontext) {
    // 判断异常类型
    if (fetestexcept(FE_DIVBYZERO)) {
        printf("异常原因：除零错误\n");
    }
    else if (fetestexcept(FE_INVALID)) {
        printf("异常原因：无效浮点运算（如 sqrt(-1)）\n");
    }
    else if (fetestexcept(FE_OVERFLOW)) {
        printf("异常原因：数值溢出\n");
    }
    else if (fetestexcept(FE_UNDERFLOW)) {
        printf("异常原因：数值下溢\n");
    }
    else if (fetestexcept(FE_INEXACT)) {
        printf("异常原因：精度丢失\n");
    }
		printf("错误地址：%p\n", info->si_addr); 
		print_call_stack(0);
    // 清除异常标志
    feclearexcept(FE_ALL_EXCEPT);

}

// 注册段错误信号处理函数
void register_segfault_handler_() {
    struct sigaction sa;
    memset(&sa, 0, sizeof(sa));
    sa.sa_sigaction = segfault_handler; // 自定义处理函数
    sa.sa_flags = SA_SIGINFO;           // 启用详细信息（si_addr）
    // 注册 SIGSEGV 信号
    if (sigaction(SIGSEGV, &sa, NULL) == -1) {
        perror("sigaction failed");
        exit(EXIT_FAILURE);
    }
    // 2. 启用浮点异常捕获（必须加！）
    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    sa.sa_sigaction = sigfpe_handler;
    sa.sa_flags = SA_SIGINFO;
    sigemptyset(&sa.sa_mask);
    sigaction(SIGFPE, &sa, NULL);
}


