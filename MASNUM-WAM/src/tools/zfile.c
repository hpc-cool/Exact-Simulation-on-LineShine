#define _FILE_OFFSET_BITS 64
#define __USE_LARGEFILE64
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include <string.h>
#include<stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include<errno.h>
#define lseek64 lseek
//#define FORTRANUNDERSCORE
#if ( defined FORTRANCAPS )
#define zf_mkdir  ZF_MKDIR
#define zf_setpid ZF_SETPID
#define zf_size   ZF_SIZE
#define zf_open   ZF_OPEN
#define zf_skip   ZF_SKIP
#define zf_read   ZF_READ
#define zf_write  ZF_WRITE
#define zf_close  ZF_CLOSE
#elif ( defined FORTRANDOUBLEUNDERSCORE )
#define zf_mkdir  zf_mkdir__
#define zf_setpid zf_setpid__
#define zf_size   zf_size__
#define zf_open   zf_open__
#define zf_skip   zf_skip__
#define zf_read   zf_read__
#define zf_write  zf_write__
#define zf_close  zf_close__
#else
//default FORTRANUNDERSCORE 
#define zf_mkdir  zf_mkdir_
#define zf_setpid zf_setpid_
#define zf_size   zf_size_
#define zf_open   zf_open_
#define zf_skip   zf_skip_
#define zf_read   zf_read_
#define zf_write  zf_write_
#define zf_close  zf_close_
#endif
#define off_t long long
typedef struct _ZFILE{
  char fn[1024];
  int fd;
  short swap,ftype;
}ZFILE;

#define NFILE 256
ZFILE ZFiles[NFILE]={0};
static int ZFInited=0;
static int tag=0x12345678;
static FILE*flog=NULL;
static int mypid=0;
extern int zmsg (__const char *__restrict form, ...);
int zmsg (__const char *__restrict form, ...){
  int n=0;
  va_list arglist;
  va_start(arglist, form);
  n=vfprintf(stdout,form,arglist);
  fflush(stdout);
  return n;
}
ZFILE*GetZfile(int ifi){
  if(ifi<=0||ifi>NFILE)return NULL;
  return &ZFiles[ifi-1];
}

void swap2(char *b,int n){
  int i;
  for(i=0;i<n;i++,b+=4){
    char c;
    c=b[0];b[0]=b[1];b[0]=c;
  }
}
void swap4(char *b,int n){
  int i;
  for(i=0;i<n;i++,b+=4){
    char c;
    c=b[0];b[0]=b[3];b[3]=c;
    c=b[1];b[1]=b[2];b[2]=c;
  }
}
void swap8(char *b,int n){
  int i;
  for(i=0;i<n;i++,b+=8){
    char c;
    c=b[0];b[0]=b[7];b[7]=c;
    c=b[1];b[1]=b[6];b[6]=c;
    c=b[2];b[2]=b[5];b[5]=c;
    c=b[3];b[3]=b[4];b[4]=c;
  }
}
static char*strf2c(char*dfn,const char*fn,int nfn_){
  int nfn;
  for(nfn=nfn_-1;nfn;nfn--){
    char c=fn[nfn];
  	if(c!=' '&&c!=0)break;
  }
  nfn++;
  strncpy(dfn,fn,nfn);
  dfn[nfn]=0;
  return dfn;
}
int zopen(char*fn,int flg){
//  struct ffsw stat;
  return open(fn,flg,0666);//,0,&stat);
}
int checktag(int fd){
	int ttag;
  int ir;
  lseek(fd,0,SEEK_SET);
  ir=read(fd,&ttag,4);
  if(ttag!=tag){
    swap4((char*)&ttag,1);
    if(ttag!=tag){
      printf("zropen Tag Error\n");
      return 0;
    }
    return 2;
  }
  return 1;
}
void zf_setpid(int pid){
	mypid=pid;
}
void zf_mkdir(char*dir_,int nfn_){
  char dir[256];
  mkdir(strf2c(dir,dir_,nfn_),0755);
}
int ZF_open(int mode,const char*fn){
  int fd;
  char c,*pfn;
  ZFILE*zfi=NULL;
  int i,ifile=0,ftype;
  int ttag,swap=0;
  int oflg;
  off_t off;
  // zmsg("sizeof(off_t)=%d\n",sizeof(off_t));
  if(ZFInited==0){
    ZFInited=1;
    memset(ZFiles,0,NFILE*sizeof(ZFiles[0]));
    for(i=0;i<NFILE;i++){
      ZFiles[i].fd=-1;
    }
    ifile=1;
  }else{
    for(i=0;i<NFILE;i++){
      if(ZFiles[i].fd==-1){
        ifile=i+1;break;
      }
    }
    if(ifile==0){
      zmsg("zfopen NoMore File\n");
      return ifile;
    }
  }
  zfi=GetZfile(ifile);
  strcpy(pfn=zfi->fn,fn);
  // zmsg("zf_open A %d %d %d @@%s@@\n",ifile,*mode,nfn,pfn);
  swap=0;
  // fflag=O_LARGEFILE;
  if((mode&1)==0)oflg=O_RDWR|O_CREAT;
	else  				oflg=O_RDONLY;
	if((mode&5)==4)oflg|=O_TRUNC;
  // fd=zopen(pfn,oflg);
  fd=open(pfn,oflg,0666);
  if(fd<0){
  	perror("Open Err:");
  	zmsg("openfile %d %d %x %x %x fd=%x %s\n",mypid,mode,oflg,O_RDONLY,O_RDWR|O_CREAT,fd,pfn);
  	return 0;
  }
 	off=lseek64(fd,0,SEEK_END);
  int ir;
  if((mode&2)==0){
    ftype=1;
	  if((mode&1)==0){
  		if(off==0){
	      ir=write(fd,&tag,4);
	    }
	  }else if(mode==1){
      lseek(fd,0,SEEK_SET);
      ir=read(fd,&ttag,4);
      if((swap=checktag(fd))==0){
          printf("zropen Tag Error\n");
          ifile=0;close(fd);fd=-1;
          return 0;
      }
      if(swap>1)swap=1;
	  }
  }else{
    ftype=0;
  }
  zfi->fd=fd;
  zfi->ftype=ftype;
  zfi->swap=swap;
  return ifile;
}
int zf_open(int *mode,char*fn,int nfn_){
  char fnt[1024];
  return ZF_open(*mode,strf2c(fnt,fn,nfn_));
}
void dump(void*b,int size){
  int i,n;
  n=64;if(n>size)n=size;
  for(i=0;i<n;i++)printf("%2.2x ",((char*)b)[i]);
}
off_t zf_size(int *ifile){
  ZFILE*zfi;
  if((zfi=GetZfile(*ifile))==NULL)return 0;
  return lseek64(zfi->fd,0,SEEK_END);
}
off_t zf_skip(int *ifile,off_t *off,int *size,int *count){
  return *off+=(*size)*(*count);
}
off_t zf_write(int *ifile,off_t *off,void*b,int *size,int *count){
  ZFILE*zfi;
  // zmsg("WA %d %ld\n",*ifile,*off);
  if((zfi=GetZfile(*ifile))==NULL)return 0;
  if(*off<4 && zfi->ftype)*off=4;
  // zmsg("WA %d %ld\n",*ifile,*off);
  lseek64(zfi->fd,*off,SEEK_SET);
  int ir=write(zfi->fd,b,*size*(*count));
  return lseek64(zfi->fd,0,SEEK_CUR);
}
off_t zf_read(int *ifile,off_t *off,void*b,int *size,int *count){
  ZFILE*zfi;
  if((zfi=GetZfile(*ifile))==NULL)return 0;
  if(*off<4 && zfi->ftype) *off=4;
  lseek64(zfi->fd,*off,SEEK_SET);
  int ir=read(zfi->fd,b,*size*(*count));
  if(zfi->swap){
    if(*size==4){
      swap4(b,*count);
    }else if(*size==8){
      swap8(b,*count);
    }else if(*size==2){
      swap2(b,*count);
    }
  }
  return lseek64(zfi->fd,0,SEEK_CUR);
}
int ZF_close(int ifile){
  ZFILE*zfi;
  if((zfi=GetZfile(ifile))==NULL)return 0;
  close(zfi->fd);
  zfi->fd=-1;
  zfi->swap=0;
  return 0;
}
int zf_close(int *ifile){
  return ZF_close(*ifile);
}
int PT_openw(const char*fn){
  int mode=4;
  return ZF_open(mode,fn);
}
int pt_openw_(const char*fn,int nfn_){
  char fnt[1024];
  return PT_openw(strf2c(fnt,fn,nfn_));
}
int PT_write(int ifile,int ind,double *ee,double*se){
  ZFILE*zfi;
  int lr;
  if((zfi=GetZfile(ifile))==NULL)return 0;
  lseek64(zfi->fd,sizeof(double)*12*32*2*ind,SEEK_SET);
  lr=write(zfi->fd,ee,sizeof(double)*12*32);
  lr=write(zfi->fd,se,sizeof(double)*12*32);
  return lseek64(zfi->fd,0,SEEK_CUR);
}
int pt_write_(int *ifile,int *ind,double *ee,double*se){
  return PT_write(*ifile,*ind,ee,se);
}
int PT_close_(int ifile){
  return ZF_close(ifile);
}
int pt_close_(int*ifile){
  return zf_close(ifile);
}
#if 0
void  test(const std::string &filename){
  std::ifstream input(filename, std::ios::binary);
  torch::Tensor ee=load_ee(intput,1);
  torch::Tensor se=load_se(intput,1);
}
torch::Tensor load_ee(std::ifstream input, int ind) {
  std::vector<double> buffer(12*32);
  input.seekg(12*32*sizeof(double)*(ind*2),input.beg);
  if (!input.read(reinterpret_cast<char*>(buffer.data()), 12*32 * sizeof(double))) {
    throw std::runtime_error("Error reading binary file");
  }
  return torch::from_blob(buffer.data(), {length}, torch::kFloat64).clone();
}
torch::Tensor load_se(std::ifstream input, int ind) {
  std::vector<double> buffer(12*32);
  input.seekg(12*32*sizeof(double)*(ind*2),input.beg);
  if (!input.read(reinterpret_cast<char*>(buffer.data()), 12*32 * sizeof(double))) {
    throw std::runtime_error("Error reading binary file");
  }
  return torch::from_blob(buffer.data(), {length}, torch::kFloat64).clone();
}
#endif
