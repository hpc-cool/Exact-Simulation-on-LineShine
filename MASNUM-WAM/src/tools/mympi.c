#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#define NEWN(T,N)  (T*)malloc(sizeof(T)*(N))
#define max(a,b) (a>b?a:b)
#define MYMPI_GATHER   mympi_gather_
#define MYMPI_GATHERV  mympi_gatherv_
#define MYMPI_SCATTER  mympi_scatter_
#define MYMPI_SCATTERV  mympi_scatterv_
#define HMPI
#ifdef HMPI
#define FCOMMT MPI_Fint
#define FCOMM2C(comm)  MPI_Comm_f2c(*comm)
#define FTYPE MPI_Fint
#define FTYPE2C(typ) MPI_Type_f2c(*typ)
#else
#define FCOMMT MPI_Comm
#define FCOMM2C(comm)  (*comm)
#define FTYPE MPI_Datatype 
#define FTYPE2C(typ) (*typ)
#endif
int MyMPI_Gather(void * sendbuf, int sendcount, MPI_Datatype sendtype,
                 void * recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MyMPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                  void* recvbuf, int recvcount, MPI_Datatype recvtype,int root, MPI_Comm comm);
int MyMPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,int *recvcounts, int *displs, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MyMPI_Scatterv(void* sendbuf, int *sendcounts, int *displs, MPI_Datatype sendtype,
                 void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,MPI_Comm comm);

void MYMPI_GATHER(void * sendbuf, int *sendcount, FTYPE*sendtype,
                  void * recvbuf, int *recvcount, FTYPE*recvtype, int *root, FCOMMT * comm,int*ierr){
    *ierr=MyMPI_Gather(sendbuf,*sendcount,FTYPE2C(sendtype),recvbuf,*recvcount,FTYPE2C(recvtype),*root,FCOMM2C(comm));
}
void MYMPI_SCATTER(void* sendbuf, int *sendcount, FTYPE*sendtype,
                void* recvbuf, int *recvcount, FTYPE*recvtype,int *root, FCOMMT * comm,int *ierr){
  *ierr=MyMPI_Scatter(sendbuf, *sendcount,FTYPE2C(sendtype),recvbuf, *recvcount,FTYPE2C(recvtype),*root,FCOMM2C(comm));
}
void MYMPI_GATHERV(void* sendbuf, int *sendcount, FTYPE*sendtype, void* recvbuf,
      int *recvcounts, int *displs,FTYPE*recvtype, int *root, FCOMMT * comm,int*ierr){
  *ierr=MyMPI_Gatherv(sendbuf, *sendcount,FTYPE2C(sendtype), recvbuf,recvcounts, displs,FTYPE2C(recvtype), *root,FCOMM2C(comm));
}
void MYMPI_SCATTERV(void* sendbuf, int *sendcounts, int *displs, FTYPE*sendtype,
                 void* recvbuf, int *recvcount, FTYPE*recvtype, int *root,FCOMMT * comm,int *ierr){
  *ierr=MyMPI_Scatterv(sendbuf, sendcounts, displs,FTYPE2C(sendtype),recvbuf, *recvcount,FTYPE2C(recvtype), *root,FCOMM2C(comm));
}

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#define DEBUGMPI
#ifdef DEBUGMPI
static FILE*flog=NULL;
static char fn[256];
static char NodeName[256];
void OpenLog(int pid){
  if(pid>=0){
    if(!flog){
      int len;
      MPI_Get_processor_name(NodeName,&len);
      sprintf(fn,"clogs/%s_%4.4d.txt",NodeName,pid);
      flog=fopen(fn,"wt");
    }
  }else{
    if(flog){
      fprintf(flog,"End Log\n");
      fclose(flog);flog=NULL;
    }
  }
}
static void my_log(const char *format, ...) {
	if(!flog)return;
	va_list ap;
	va_start(ap, format);
	int len = vfprintf(flog, format, ap);
	va_end(ap);
	fflush(flog);
}
#define MYLOG(a) my_log a
#else
void OpenLog(int pid){}
#define MYLOG(a) 
#endif
static int myid,myidr=-1,numprocs=0,_myidr=-1;
MPI_Comm _comm=NULL;
void InitMyMpi(MPI_Comm comm,int root){
  if(_comm!=comm){
    _comm=comm;
    MPI_Comm_rank( comm, &myidr);
    MPI_Comm_size( comm, &numprocs);
    if(_myidr!=myidr){
      _myidr=myidr;
      OpenLog(myidr);
    }
  }
  myid=myidr-root;
  if(myid<0)myid+=numprocs;
}
int MyMPI_Gather(void * sendbuf, int sendcount, MPI_Datatype sendtype,
                 void * recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm){
    // assume recvcount==SUM(sendcount)
    // assume root is 0 now
    // assume simple datatype
    int isize,dsize,msize;
    char *pbuf=NULL,*rbuf;
    MPI_Status status;
    int deltaid;
    int inr,nnr,nns,i,*rcs;
    if(!sendbuf)return -101;
    InitMyMpi(comm,root);
    if(myid==0&&!recvbuf)return -102;
    MPI_Type_size(sendtype,&isize);
    for(nns=1,nnr=0;nns<numprocs;nnr++,nns*=2);
    if(myid==0){
      msize=numprocs*recvcount*isize;
      if(root==0)pbuf=recvbuf;
      else pbuf=NEWN(char,msize);
    }
    for(deltaid=1,nnr=0;deltaid<numprocs;nnr++,deltaid*=2);
    for(inr=0;inr<nnr;inr++){
    	int nid;
      deltaid/=2;
      if(myid%deltaid)continue;
      nid=myid/deltaid;
      if(nid%2){
        msize=deltaid*recvcount*isize;
        pbuf=NEWN(char,msize);
        break;
      }
    }
    dsize=isize*sendcount;
    memcpy(pbuf,sendbuf,dsize); //copy self Data
    deltaid=1;
    for(inr=0;inr<nnr;inr++){
    	int nid=myid/deltaid;
      if(nid%2==0){ //recv
      	int sid=myid+deltaid; //cal Src ID
        if(sid<numprocs){
          int rsize;
          rbuf=pbuf+dsize;    //recvbuf
          sid+=root;
          if(sid>=numprocs)sid-=numprocs;
          MYLOG(("gather recv %d <= %d %d %d ",myidr,sid,inr,msize));
          MPI_Recv(rbuf,msize,MPI_CHAR,sid,10001,comm,&status);
					MYLOG(("end\n"));
          MPI_Get_count(&status,MPI_CHAR,&rsize);
          dsize+=rsize;//have data
        }
      }else{ //send
        int did=myid-deltaid; //cal Dst Id
        did+=root;
        if(did>=numprocs)did-=numprocs;
        MYLOG(("gather send %d => %d %d %d ",myidr,did,inr,dsize));
        MPI_Send(pbuf,dsize,MPI_CHAR,did,10001,comm);
				MYLOG(("end\n"));
        break;
      }
      deltaid*=2;
    }
    if(myid==0&&root!=0){
      //exchange Data Pos Root & 0
      int roff,doff;
      roff=recvcount*root;
      doff=recvcount*(numprocs-root);
      memmove(recvbuf           ,pbuf+isize*doff,isize*roff);
      memcpy (((char*)recvbuf)+isize*roff,pbuf,isize*doff);
    }
    if(pbuf!=recvbuf)free(pbuf);
    return 0;
}

int MyMPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                  void* recvbuf, int recvcount, MPI_Datatype recvtype,int root, MPI_Comm comm){
    // assume SUM(recvcount)==sendcount
    // assume root is 0 now
    // assume simple datatype
    int myid,isize,dsize,msize;
    char *pbuf;
    MPI_Status status;
    int deltaid;
    int inr,nnr,nns,i,soff;
    //int dnb=recvcounts[root]-recvcounts[0];
    if(!sendbuf)return -101;
    InitMyMpi(comm,root);
    if(myid==0&&!recvbuf)return -102;
    MPI_Type_size(sendtype,&isize);
    for(nns=1,nnr=0;nns<numprocs;nnr++,nns*=2);
    if(myidr==root){
      msize=isize*numprocs*sendcount;
      if(root==0) pbuf=sendbuf;
      else {
        int roff,doff;
        pbuf=NEWN(char,msize);
        //exchange Data Pos Root & 0
        roff=sendcount*root;
        doff=sendcount*(numprocs-root);
        memmove(pbuf           ,((char*)sendbuf+isize*roff),isize*doff);
        memcpy (pbuf+isize*doff,sendbuf     ,isize*roff);
      }
    }

    for(deltaid=1,nnr=0;deltaid<numprocs;nnr++,deltaid*=2);
    for(inr=0;inr<nnr;inr++){
    	int nid;
      deltaid/=2;
      if(myid%deltaid)continue;
      nid=myid/deltaid;

      if(nid%2==0){//send
        int did=myid+deltaid; //cal Dst Id
        if(did<numprocs){
          soff=deltaid*recvcount*isize;
          dsize=deltaid*recvcount*isize;
          if(myid+deltaid*2>=numprocs)
            		dsize=(numprocs-myid-deltaid)*recvcount*isize;
          else  dsize=(deltaid)*recvcount*isize;
          did+=root;
          if(did>=numprocs)did-=numprocs;
					MYLOG(("scatter send %d => %d %d %d %d ",myidr,did,inr,dsize,soff));
          MPI_Send(pbuf+soff,dsize,MPI_CHAR,did,10001,comm);
          MYLOG(("end\n"));
        }
      }else{//recv
        int sid=myid-deltaid; //cal Src ID
        if(myid+deltaid>=numprocs)
          		dsize=(numprocs-myid)*recvcount*isize;
        else  dsize=(deltaid)*recvcount*isize;
        dsize*=isize;
        pbuf=NEWN(char,msize);
        sid+=root;
        if(sid>=numprocs)sid-=numprocs;
        MYLOG(("scatter recv %d <= %d %d %d ",myidr,sid,inr,dsize));
        MPI_Recv(pbuf,dsize,MPI_CHAR,sid,10001,comm,&status);
        MYLOG(("end\n"));
      }
    }
    memcpy(recvbuf,pbuf,isize*recvcount); //copy self Data
    if(pbuf!=sendbuf)free(pbuf);
    //OpenLog(-1);
    return 0;
}
int MyMPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,int *recvcounts, int *displs, MPI_Datatype recvtype, int root, MPI_Comm comm){
    // assume recvcounts & displs have value in all procs
    // assume recvcounts[mpi_id]==sendcount_of_mpi_id
    // assume root is 0 now
    // assume simple datatype
    int isize,dsize,msize;
    char *pbuf=NULL,*rbuf;
    MPI_Status status;
    int deltaid,inr,nnr,nns,i,*rcs;
    int NeedSort=0;

    if(!sendbuf)return -101;
    if(!recvcounts)return -103;
    InitMyMpi(comm,root);
    if(myidr==root&&!recvbuf)return -102;
    MPI_Type_size(sendtype,&isize);

    if(myidr==root){
      pbuf=recvbuf;
      if(root!=0)NeedSort=1;
      if(displs){
        for(msize=0,i=0;i<numprocs;i++){
          if(displs[i]!=msize)NeedSort=2;
          msize+=recvcounts[i];
        }
      }else if(NeedSort){
        for(msize=0,i=0;i<numprocs;i++){
          msize+=recvcounts[i];
        }
      }
      msize*=isize;
      if(NeedSort){
        pbuf=NEWN(char,msize);
      }
    }

    for(nns=1,nnr=0;nns<numprocs;nnr++,nns*=2);
    rcs=NEWN(int,nns);
    memcpy(rcs              ,recvcounts+root,sizeof(int)*(numprocs-root));
    memcpy(rcs+numprocs-root,recvcounts     ,sizeof(int)*root);
    memset(rcs+numprocs,0,sizeof(int)*(nns-numprocs));

    deltaid=nns;
    for(inr=0;inr<nnr;inr++){
    	int nid;
      deltaid/=2;
      if(myid%deltaid)continue;
      nid=myid/deltaid;
      if(nid%2!=0){
        for(msize=0,i=0;i<deltaid;i++)msize+=rcs[myid+i];
        msize*=isize;
        pbuf=NEWN(char,msize);
        break;
      }
    }
    dsize=isize*sendcount;
    memcpy(pbuf,sendbuf,dsize); //copy self Data
    deltaid=1;
    for(inr=0;inr<nnr;inr++){
    	int nid=myid/deltaid;
      if(nid%2==0){ //recv
      	int sid=myid+deltaid; //cal Src ID
        if(sid<numprocs){
          int rsize;
          rbuf=pbuf+dsize;    //recvbuf
          sid+=root;
          if(sid>=numprocs)sid-=numprocs;
          //fprintf(flog,"RECV %d %d <= %d %d\n",myidr,__LINE__,sid,msize);fflush(flog);
          MYLOG(("gatherv recv %d <= %d %d %d ",myidr,sid,inr,msize));
          MPI_Recv(rbuf,msize,MPI_CHAR,sid,10001,comm,&status);
					MYLOG(("end\n"));
          MPI_Get_count(&status,MPI_CHAR,&rsize);
          dsize+=rsize;//have data
        }
      }else{ //send
        int did=myid-deltaid; //cal Dst Id
        did+=root;
        if(did>=numprocs)did-=numprocs;
        MYLOG(("gather send %d => %d %d %d ",myidr,did,inr,dsize));
        MPI_Send(pbuf,dsize,MPI_CHAR,did,10001,comm);
				MYLOG(("end\n"));
        break;
      }
      deltaid*=2;
    }
    //exchange Data Pos Root & 0
    if(NeedSort==1){
      int roff,doff;
      for(roff=0,i=0;i<root;i++)roff+=recvcounts[i];
      for(doff=0,i=root;i<numprocs;i++)doff+=recvcounts[i];
      memcpy(recvbuf          ,pbuf+isize*doff,isize*roff);
      memcpy(((char*)recvbuf+isize*roff),pbuf           ,isize*doff);
    }else if(NeedSort==2){
      int roff;
      for(roff=0,i=0;i<numprocs;i++){
        int n,it=i+root;
        if(it>=numprocs)it-=numprocs;
        n=recvcounts[it]*isize;
        //fprintf(flog,"%d %d %d %d\n",myidr,__LINE__,i,it,roff/isize,displs[it],n/isize);fflush(flog);
        memcpy(((char*)recvbuf+displs[it]*isize),pbuf+roff,n);
        roff+=n;
      }
    }
    if(pbuf!=recvbuf)free(pbuf);
    if(rcs!=recvcounts)free(rcs);
    //OpenLog(-1);
    return 0;
}

int MyMPI_Scatterv(void* sendbuf, int *sendcounts, int *displs, MPI_Datatype sendtype,
                 void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,MPI_Comm comm){
    // assume recvcounts & displs have value in all procs
    // assume recvcounts[mpi_id]==sendcount_of_mpi_id
    // not assume root is 0 now
    int isize,dsize,msize;
    char *pbuf,*rbuf;
    MPI_Status status;
    int deltaid;
    int inr,nnr,nns,i,dd,soff,it;
    int *scs;
    int NeedSort=0;
    //Check inputs
    if(!recvbuf)return -101;
    if(!sendcounts)return -103;
    InitMyMpi(comm,root);
    if(myidr==root&&!sendbuf)return -102;
    MPI_Type_size(sendtype,&isize);
    if(myidr==root){
      NeedSort=0;
      pbuf=sendbuf;
      if(root!=0)NeedSort=1;
      if(displs){
        for(msize=0,i=0;i<numprocs;i++){
          if(displs[i]!=msize)NeedSort=2;
          msize+=sendcounts[i];
        }
      }else if(NeedSort){
        for(msize=0,i=0;i<numprocs;i++){
          msize+=sendcounts[i];
        }
      }
      msize*=isize;
      if(NeedSort){
        pbuf=NEWN(char,msize);
        if(NeedSort==1){
          int roff,doff;
          //exchange Data Pos Root & 0
          for(roff=0,i=0;i<root;i++)roff+=sendcounts[i];
          for(doff=0,i=root;i<numprocs;i++)doff+=sendcounts[i];
          memcpy(pbuf           ,((char*)sendbuf+isize*roff),isize*doff);
          memcpy(pbuf+isize*doff,sendbuf           ,isize*roff);

        }else{
          int roff;
          for(roff=0,i=0;i<numprocs;i++){
            int n,it=i+root;
            if(it>=numprocs)it-=numprocs;
            n=sendcounts[it]*isize;
            //            fprintf(flog,"%d %d %d %d\n",myidr,__LINE__,i,it,roff/isize,displs[it],n/isize);fflush(flog);
            memcpy(pbuf+roff,((char*)sendbuf+displs[it]*isize),n);
            roff+=n;
          }
        }
      }
    }
    for(nns=1,nnr=0;nns<numprocs;nnr++,nns*=2);
    scs=NEWN(int,nns);
    memcpy(scs              ,sendcounts+root,sizeof(int)*(numprocs-root));
    memcpy(scs+numprocs-root,sendcounts     ,sizeof(int)*root);
    memset(scs+numprocs,0,sizeof(int)*(nns-numprocs));
    deltaid=nns;
    for(inr=0;inr<nnr;inr++){
    	int nid;
      deltaid/=2;
      if(myid%deltaid)continue;
      nid=myid/deltaid;
      if(nid%2==0){        //send
        int did=myid+deltaid; //cal Dst Id
        if(did<numprocs){
          for(dsize=0,soff=0 ,i=0;i<deltaid;i++){
            soff +=scs[myid+i];
            dsize+=scs[ did+i];
          }
          did+=root;if(did>=numprocs)did-=numprocs;
					MYLOG(("scatterv send %d => %d %d %d %d ",myidr,did,inr,dsize,soff));
          MPI_Send(pbuf+soff*isize,dsize*isize,MPI_CHAR,did,10001,comm);
          MYLOG(("end\n"));
        }
      }else{ //recv
        int sid=myid-deltaid; //cal Src ID
        for(dsize=0,i=0;i<deltaid;i++)dsize+=scs[myid+i];
        dsize*=isize;
        pbuf=NEWN(char,dsize);
        sid+=root;if(sid>=numprocs)sid-=numprocs;
        MYLOG(("scatterv recv %d <= %d %d %d ",myidr,sid,inr,dsize));
        MPI_Recv(pbuf,dsize,MPI_CHAR,sid,10001,comm,&status);
        MYLOG(("end\n"));
      }
    }
    memcpy(recvbuf,pbuf,isize*recvcount); //copy self Data
    free(scs);
    return 0;
}
