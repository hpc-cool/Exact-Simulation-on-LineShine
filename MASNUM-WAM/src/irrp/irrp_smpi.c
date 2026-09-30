#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
//#include "irrp_inner.h"
#define NEWPN(p,N)   *((void**)&(p))=malloc(sizeof(*(p))*(N))
#define RENEWN(p,N)  *((void**)&(p))=realloc(p,sizeof(*(p))*(N))
#define VFREE(p)      if(p)free(p);p=NULL
typedef struct mpipacket{
  int bsize ,pos ,dsize ,*buf;
}mpipacket;
void irrp_abort(const char *file, int line,const char*msg){
    printf("IRRP Error %s %d %s \n", file,line,msg);
#ifndef SERIAL_TEST
    MPI_Finalize(); 
#endif
    exit(-101);
}
#ifndef SERIAL_TEST
//-------------------------------------------------------------------------------------------------
// Set dsize and pos equal to zero as initial of mpipacket
// If lsize is given, allocate memory for the buffer.
void VMpiPacket(mpipacket*pk,int lsize){
  if(pk->bsize<pk->pos+lsize){
    pk->bsize += lsize+1024;  // set it as bsize and use it to allocate buf.
    RENEWN(pk->buf,(pk->bsize+20)/4);
  }
}
void InitMpiPacket(mpipacket*pk,int lsize){
  pk->dsize = 0; pk->pos = 0;  // Initialize dsize and pos for the begining of packet.
  if(lsize == 0){
    VFREE(pk->buf); // If buf is not empty, deallocate it first.
    pk->bsize = 0 ;      // And set bsize equal to zero.
  }else if(pk->bsize<lsize){
    pk->bsize = lsize;  // set it as bsize and use it to allocate buf.
    RENEWN(pk->buf,(pk->bsize+20)/4);
  }
}
// Get dsize before using it (i.e. send, isend, bcast, unpack, recv, irecv, ...)
void SetMpiPacketDSize(mpipacket*pk, int dsize){
  if(dsize>0) pk->dsize = dsize;
  else pk->dsize = pk->pos;
}

void bcast_packet(mpipacket*pk,int root,int pid,MPI_Comm mpicomm){
  int lsize, ierr;
  if(pid == root){
    SetMpiPacketDSize(pk,0);
    lsize = pk->dsize;
  }
  ierr=MPI_Bcast(&lsize, 1, MPI_INTEGER4, root, mpicomm);
  if(pid != root){
    InitMpiPacket(pk, lsize+100);
  }
  ierr=MPI_Bcast(pk->buf, lsize, MPI_PACKED, root, mpicomm);
  if(pid != root){
    SetMpiPacketDSize(pk, lsize);
  }
}

#endif

//-------------------------------------------------------------------------------------------------

#ifndef SERIAL_TEST

void next(int flg,int pid,int npe,MPI_Comm mpicomm){
  int ierr,v;
  fflush(stdout);
  if(flg==2){
    if(pid==npe-1)ierr=MPI_Send(&v,1,MPI_INTEGER,0,1000,mpicomm);
  }else{
    if(pid<npe-1){
      ierr=MPI_Send(&v,1,MPI_INTEGER,pid+1,1000,mpicomm);
    }else if(flg!=0){
      ierr=MPI_Send(&v,1,MPI_INTEGER,0,1000,mpicomm);
    }
  }
}

void prev(int flg,int pid,int npe,MPI_Comm mpicomm){
    int ierr,v;
    MPI_Status status;
    if(flg==2){
      if(pid==0)
        ierr=MPI_Recv(&v,1,MPI_INTEGER,npe-1,1000,mpicomm,&status);
    }else{
      if(pid>0){
        ierr=MPI_Recv(&v,1,MPI_INTEGER,pid-1,1000,mpicomm,&status);
      }else if(flg){
        ierr=MPI_Recv(&v,1,MPI_INTEGER,npe-1,1000,mpicomm,&status);
      }
    }
}

void prevnext (int flg,int pid,int npe,MPI_Comm mpicomm){
    prev(flg, pid, npe, mpicomm);
    next(flg, pid, npe, mpicomm);
}

#endif

//-------------------------------------------------------------------------------------------------

#ifndef SERIAL_TEST

void barr(int pid,int npe,MPI_Comm mpicomm){
    int it=1,ierr;
    int *vi;
    NEWPN(vi,npe);
   	ierr=MPI_Gather(&it,1,MPI_INTEGER,vi,1,MPI_INTEGER,0,mpicomm);
    VFREE(vi);
}

#endif
