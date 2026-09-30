#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <mpi.h>
#include "irrp_inner.h"
#define DBGINF   printf("%s %d %d\n",__FILE__,__LINE__,gsi->pid)
#define DBGB0(msg); //if(gsi->pid>=0)  printf("%5d %5d %8.3fs %s\n", __LINE__,gsi->pid,Difftimer(),msg);fflush(stdout)
//#define MPI_GATHER   MYMPI_GATHER
//#define MPI_GATHERV  MYMPI_GATHERV
//#define MPI_SCATTER  MYMPI_SCATTER
//#define MPI_SCATTERV  MYMPI_SCATTERV
enum {
  VT_NULL=0,VT_FLOAT=1, VT_DOUBLE,VT_INT,  VT_SHORT ,
  VT_MAX
};
static int vtsize[]={1,4,8,2,4};
static MPI_Datatype mpidt[]={
  MPI_CHAR,MPI_FLOAT , MPI_DOUBLE, MPI_SHORT , MPI_INT   
};

//-------------------------------------------------------------------------------------------------
typedef struct vars_def{                        // Variables used for exchange.
  char  vname[16];        // Name of variable.
  int   vartype ;         // Type of variable:0=null, 1 = r4, 2 = r8, 3 = i2, 4 = i4
  MPI_Datatype mpivt;     //
  int   km    ;           // Number of Vertical layers.
                          // This should be in front, i.e. (km, im, jm)
                          // If it is not the case, it can be done one by one.
                          // The follows are 4 kinds of ponter variables.
                          // v**: the variable need to be exchanged.
                          // They need to be setted by irrp_exg_setvar, and use
                          // nullify to clean up.
  int   itemsize;
  BYTE* vp;
  //short *vi2;        // Type of int (kind=2), dims of (km,snpc).
  //int   *vi4;        // Type of int (kind=4), dims of (km,snpc).
  //float *vr4;        // Type of real(kind=4), dims of (km,snpc).
  //double*vr8;        // Type of real(kind=8), dims of (km,snpc).
}vars_def;
typedef struct exchange_group{
  int   nv ;  // Number of variables.
  int   mv ;  // Max-length of exchange buffer.
  int   state ;
  // 0:for not initting . 0 means not.
  // 1:inited
  // 2:sending
  // 3:recving & sending
  // 4:check end
  // 5:wait end
  // 6:exchange end
  vars_def*vardef;
  // Variables defining. (nv)
  mpipacket*sbuf; // Buffer for sending info/data. (mv)
  mpipacket*rbuf; // Buffer for receiving info/data. (mv)
}exchange_group;
typedef struct force_info{   // needed for scatter_ext, scatter_force
  int nx, ny, nxy;    // Number of input dim-size of forcing data.
  int nfp   ;         // Number points of input forcing data
  int infp  ;         // input pnts for cpart
  int infpg ;         // input pnts for all part infpg = sum(infp)
  int*fpos;           // ????(infp)
  int*fposg;          // Start position of each part (infpg), only master id
  int*infps;          // Input points of each part (npe).
  int*fdispls;        // Displacements of each PEs (npe).
                      // The follows are 4 kinds of ponter variables.
                      // gsv** variables continuely stored in global. (gnpc)
                      // lv**  local variables in current partition. (snpc)
  int itemsize;
  BYTE*gsvp,*lvp;  // 
  int *idx;           // Indeces of each related points (4 * nfp) for interp.
  double*w;             // Weighting of each related points (4 * nfp)
}force_info;
typedef struct globle_buf{                      // --- Variables related with global buffer.
  PI_Ginfo  gsi;
  PI_Part  cpart;
  int gnpc;            // Number of global calculation points.
  int*gdispls;         // Displacements of each partition (npe).
                       // Being used together with gsi.npcs in functions of
                       // MPI_Gatherv & MPI_Scatterv (npe).
  int*pposg;           // Relationship between global storage and partition
                       // storage (gnpc).
                       // The follows are 4 kinds of ponter variables.
                       // gsv** variables continuely stored in global. (gnpc)
                       // lv**  local variables in current partition. (snpc)
  int itemsize;
  BYTE*gsvp,*lvp;
  int  maxgroups ,maxforces ;
  struct exchange_group*exgroups;
  struct force_info*forces;
  struct force_info*scatterext;
}globle_buf;
//-------------------------------------------------------------------------------------------------
globle_buf gsd0={0};
globle_buf *gsds[20]={&gsd0,0};
globle_buf *gsd=&gsd0;
PI_Ginfo *gsi=&gsd0.gsi;
PI_Part  *cpart=&gsd0.cpart;
static int curgpid=0;
int setgpid(int gpid){
  if(gpid<0)return curgpid;
  if(gpid==0){
    gsd=&gsd0;
  }
  gsi=&gsd->gsi;
  cpart=&gsd->cpart;
  curgpid=gpid;
}
int getitemsize(int vt){
  if(vt<VT_MAX)return vtsize[vt];
  return 1;
}
/*-------------------------------------------------------------------------------------------------
  // init Part:
  // input:
  //  partmode: ??????????????!!!!!!!!!!!!!!!
  //  npe: number of PE
  //  pid: cur pid
  //  mpi_comm:mpi comm
  //    mask_:Mask for horizonal space/points in which 0 means land/useless points.????
  // halosize:
  // cycle_flag: Open boundary type.
  //   0 - None cycle boundary.
  //   1 - Cycled boundary in x-direction for (normal/displaced) polar grid.
  //   2 - Cycled boundary in x-direction and y-direction for triple polar grid.
  //   3 - Cycled boundaries are given manually.
  // scycle:
  //
  //-----------------------------------------------------------------------------------------------*/

void irrp_init(int partmode, int npe,int  pid,MPI_Comm mpi_comm,int*mask,int im,int jm,int  halosize,int  cycle_flag,int*scycle,int NS3){
  int*pemask=NULL;
  if(!halosize ) halosize =1;
  if(!cycle_flag){
    if(scycle) cycle_flag=3;
    else cycle_flag=0;
  }
  gsi->npe      = npe;                      // Number of PEs.
  gsi->pid      = pid;                      // Current ID of PE, start from 0.
  gsi->mpi_comm = mpi_comm;                 // Common world of MPI.
  NEWPN(gsi->npcs,gsi->npe);             // allocate npcs
  gsi->im       = im;                       // Size of grid matrix: first dimension.
  gsi->jm       = jm;                       // Size of grid matrix: second dimension.
  gsi->ijm      = gsi->im*gsi->jm;            // Maximum computing points.
  irrp_part_init(partmode,npe,pid,mpi_comm, mask,pemask,
                 im,jm, halosize, cycle_flag, scycle,NS3);
  gsi->parttype = 1;
  gsi->gnpc=0;
  for(int i=0;i<gsi->npe;i++){
    gsi->gnpc += gsi->npcs[i];                // Number of all computer points for all PEs.
  }
  //initgsd();
  //DBGB0("Initgsd end")
}
void initgsd(){
  int   ierr;
  if(gsd->gdispls)return;

  NEWPN(gsd->gdispls,gsi->npe);
  gsd->gdispls[1] = 0;
  DBGB0("Initgsd 1");
  for(int m=1;m<gsi->npe;m++){
    gsd->gdispls[m] = gsd->gdispls[m-1] + gsi->npcs[m-1];
  }
  DBGB0("Initgsd 2");

  NEWZN(gsd->pposg,1+gsi->gnpc);
  DBGB0("Initgsd 3");
  ierr=MPI_Gatherv(cpart->ppos+1, cpart->snpc, MPI_INT, gsd->pposg+1,
                   gsi->npcs, gsd->gdispls, MPI_INT, 0, gsi->mpi_comm );
  DBGB0("Initgsd 4");
  ierr=MPI_Bcast(gsd->pposg, gsi->gnpc+1, MPI_INT, 0, gsi->mpi_comm);
  DBGB0("Initgsd 5");
  gsd->gsvp=gsd->lvp=NULL;
}
void gsd_final(){
  if(gsi->npe==0)return;
  gsi->npe = 0;
  VFREE(gsd->gdispls);
  VFREE(gsd-> pposg) ;
  VFREE(gsd-> gsvp) ;
  VFREE(gsd->  lvp) ;
}
void irrp_init_data(){
  initgsd();
}
// Clear all exchange group
void irrp_exg_final_all(){
  for(int i=0; i<gsd->maxgroups;i++){
    irrp_exg_final(i);
  }
}
void irrp_final(int mode){
  irrp_part_final(mode);
  if(mode>0){
    irrp_exg_final_all();
    irrp_scatter_force_final_all();
    gsd_final();
  }
}

void test_(){

}
//-------------------------------------------------------------------------------------------------
// Init gsd Data
void irrp_output_pposg(){
  if(gsi->pid==0){
    FILE*fo=fopen("pposg.dat","wt");
    fprintf(fo,"PPOSG %d \n",gsi->gnpc);
    for(int i=1;i<=gsi->gnpc;i++){
      fprintf(fo,"%d ",gsd->pposg[i]);
    }
    fclose(fo);
  }
}
//-------------------------------------------------------------------------------------------------

#ifndef NOEXCHANGE__
// check var id of exchange group
void VerifyVarDef(exchange_group*exg,int*ivar_){
  vars_def*vardef; // nv
  int ivar=0;
  if(ivar_)ivar=*ivar_;                         
  int mm,i;
  if(ivar==0){
    if(exg->vardef)free(exg->vardef);
    exg->vardef=NULL; exg->mv=0;exg->nv=0;
    return;
  }
  vardef=exg->vardef;
  if(ivar<0){
    ivar=exg->mv+1;
    for(int i=0;i<exg->mv;i++){
      if(vardef[i].vartype<=0){ ivar=i+1; break; }
    }
  }
  if(ivar_)*ivar_=ivar;
  if(exg->mv<ivar){
    mm=exg->mv;if(mm<0)mm=0;
    exg->mv=ivar;
    RENEWZN(exg->vardef,exg->mv);
  }
}

void Init_exchange_group(exchange_group*exg,int nv){
  int i,iv;
  exg->mv=0;exg->nv=0;exg->state=0;
  exg->vardef=NULL;
  exg->sbuf  =NULL;
  exg->rbuf  =NULL;
  iv=nv;if(iv<0)iv=0;
  VerifyVarDef(exg,&iv);
}
// Init Exchange
// group_ind:index of Exchange group,<=0 ,Find Index
// nvar_: Number Vars in the group,optional ,default 4
void irrp_exg_init(int*group_ind_,int nvar) {
  // Init exchange_group,can
  exchange_group*exgt;
  if(nvar<=0)nvar=4;
  int group_ind=*group_ind_;
  // check group_ind
  if(group_ind<=0){
    for(int i=0;i<gsd->maxgroups;i++){
      if(gsd->exgroups[i].mv==0){ group_ind=i+1; break; }
    }
  }
  if(group_ind>gsd->maxgroups){
    int mm=gsd->maxgroups;gsd->maxgroups=group_ind+4;
    RENEWN(gsd->exgroups,gsd->maxgroups);
    for(int i=mm;i<gsd->maxgroups;i++){
      Init_exchange_group(gsd->exgroups+i,0);
    }
  }
  *group_ind_=group_ind;
  // Clear the group
  irrp_exg_final(group_ind);
  // Init the group
  Init_exchange_group(gsd->exgroups+group_ind,nvar);
}
// Clear exchange group
void irrp_exg_final(int group_ind){
  exchange_group*exg;
  exg=&gsd->exgroups[group_ind];
  for(int i=0;i<cpart->nids;i++){
    InitMpiPacket(exg->sbuf+i,0);
    InitMpiPacket(exg->rbuf+i,0);
  }
  VFREE(exg->sbuf);
  VFREE(exg->rbuf);
  VerifyVarDef(exg, 0);
  Init_exchange_group(exg,0);
}
// Set var info into exchange group
// group_ind:exchange group index,<=0 auto find it
// ivar index of var in the group
// vname:var name
// ivt:var type
// km :km
// var:var
void irrp_exg_setvar(int*group_ind_,int*ivar_,const char*vname,void*var,int ivt,int km){
  int ivar=1;
  if(ivar_)ivar=*ivar_;
  int group_ind=0;
  if(group_ind_)group_ind=*group_ind_;
  if(gsd->maxgroups==0 || (group_ind<=0|| group_ind>gsd->maxgroups) ){
    irrp_exg_init(&group_ind,ivar);
  }
  exchange_group*exg=&gsd->exgroups[group_ind];
  VerifyVarDef(exg,&ivar);
  if(ivar_)*ivar_=ivar;
  vars_def*vardef=exg->vardef;
  if(vardef[ivar].vartype==0)exg->nv=exg->nv+1;
  vardef[ivar].vartype=ivt;
  vardef[ivar].mpivt=mpidt[ivt];
  vardef[ivar].itemsize=vtsize[ivt];
  strncpy(vardef[ivar].vname,vname,16);
  vardef[ivar].km=km;
  exg->state=0;
  exg->vardef[ivar].vp=(BYTE*)var;
}
void irrp_exg_enddef(int group_ind){
  exchange_group*exg=&gsd->exgroups[group_ind];
  if(!exg->sbuf) NEWPN(exg->sbuf,cpart->nids);
  if(!exg->rbuf) NEWPN(exg->rbuf,cpart->nids);
  for(int inb=0;inb<cpart->nids;inb++){
    PI_Rsinfo*rsi=&cpart->rsinfo[inb];
    int lsize=0;
    for(int iv=0;iv<exg->mv;iv++){
      vars_def* var=&exg->vardef[iv];
      lsize+=var->itemsize*var->km;
    }
    InitMpiPacket(&exg->sbuf[inb],lsize*rsi->sn);
    InitMpiPacket(&exg->rbuf[inb],lsize*rsi->rn);
  }
  exg->state=1;
}

void exchange_group_action_prepare(int group_ind){
  if(group_ind<0|| group_ind>gsd->maxgroups ){
    printf("irrp_exg_action:group_ind is invalid %d %d\n",group_ind,gsd->maxgroups);
    irrp_abort(__FILE__, __LINE__,"");
  }
  if(!gsd->exgroups[group_ind].state){
    irrp_exg_enddef(group_ind);
  }
  if(cpart->requests==NULL){
    NEWPN(cpart->requests,cpart->nids);NEWPN(cpart->statuss,cpart->nids);
    NEWPN(cpart->requestr,cpart->nids);NEWPN(cpart->statusr,cpart->nids);
  }
  ZERON(cpart->requests,cpart->nids);ZERON(cpart->statuss,cpart->nids);
  ZERON(cpart->requestr,cpart->nids);ZERON(cpart->statusr,cpart->nids);
}
void irrp_exginf(int nst){
  for(int inb=0;inb<cpart->nids;inb++){
    PI_Rsinfo*rsi=&cpart->rsinfo[inb];
    printf("========%5d>%5d %5d\n",gsi->pid,rsi->id,rsi->sn*nst);
  }
}
void exchange_group_action_send(int group_ind){
  exchange_group*exg=&gsd->exgroups[group_ind];
  exg->state=2;
  for(int inb=0;inb<cpart->nids;inb++){
    int   ierr;
    PI_Rsinfo*rsi=&cpart->rsinfo[inb];
    int dstid =rsi->id;
    int tag=1000;
    mpipacket*pk=&exg->sbuf[inb];
    int*pnts=rsi->spnts;
    // pack Data
    InitMpiPacket(pk,10240);
    for(int iv=0;iv<exg->mv;iv++){
      vars_def*var=&exg->vardef[iv];
      int km=var->km,is=var->itemsize;
      for(int i=0;i<rsi->sn;i++){
        VMpiPacket(pk,var->itemsize*km);
        ierr=MPI_Pack(var->vp+pnts[i]*is,km,var->mpivt,pk->buf,pk->bsize,&pk->pos,gsi->mpi_comm);
      }
    }
    SetMpiPacketDSize(pk,0);
    ierr=MPI_Isend(pk->buf,pk->dsize,MPI_PACKED,dstid,tag,gsi->mpi_comm,&cpart->requests[inb]);
  }
}

void exchange_group_action_recv_unpack(int group_ind){
  exchange_group*exg=&gsd->exgroups[group_ind];
  exg->state=5;
  int   ierr;
  for(int inb=0;inb<cpart->nids;inb++){
    PI_Rsinfo*rsi=&cpart->rsinfo[inb];
    mpipacket*pk=&exg->rbuf[inb];
    //ierr=MPI_Get_count(cpart->statusr(1,inb),MPI_PACKED,lsize)
    SetMpiPacketDSize(pk,0);
    int*pnts=rsi->rpnts;
    // unpack
    for(int iv=0;iv<exg->mv;iv++){
      vars_def*var=&exg->vardef[iv];
      int km=var->km;
      for(int i=0;i<rsi->rn;i++){
        ierr=MPI_Unpack(pk->buf,pk->dsize,&pk->pos,var->vp+pnts[i]*var->itemsize,km   ,MPI_FLOAT  ,gsi->mpi_comm);
      }
    }
  }
}
void exchange_group_action_recvs(int group_ind){
  MPI_Status status;
  int   ierr;
  exchange_group*exg=&gsd->exgroups[group_ind];
  exg->state=3;
  for(int inb=0;inb<cpart->nids;inb++){
    int dstid =cpart->rsinfo[inb].id;
    int tag=1000;
    mpipacket*pk=&gsd->exgroups[group_ind].rbuf[inb];
    InitMpiPacket(pk,1);
    ierr=MPI_Recv(pk->buf,pk->bsize,MPI_PACKED,dstid,tag,gsi->mpi_comm,&status);
  }
  //ierr=MPI_Waitall(cpart->nids,cpart->requestr,cpart->statusr)
  exchange_group_action_recv_unpack(group_ind);
  ierr=MPI_Waitall(cpart->nids,cpart->requests,cpart->statuss);
}
void exchange_group_action_recv(int group_ind){
  exchange_group*exg=&gsd->exgroups[group_ind];
  exg->state=3;
  int   ierr;
  for(int inb=0;inb<cpart->nids;inb++){
    int dstid =cpart->rsinfo[inb].id;
    int tag=1000;
    mpipacket*pk=&gsd->exgroups[group_ind].rbuf[inb];
    InitMpiPacket(pk,1);
    ierr=MPI_Irecv(pk->buf,pk->bsize,MPI_PACKED,dstid,tag,gsi->mpi_comm,cpart->requestr+inb);
  }
}
void irrp_exg_start(int group_ind){
  exchange_group_action_prepare(group_ind);
  exchange_group_action_send(group_ind);
  exchange_group_action_recv(group_ind);
}
int irrp_exg_check(int group_ind){
  int  state, flag, ierr;
  state=1;
  if(gsd->exgroups[group_ind].state==3){
    ierr=MPI_Testall(cpart->nids,cpart->requestr,&flag,cpart->statusr);
    if(flag){
      state=1;
      gsd->exgroups[group_ind].state=4;
    }else{
      state=0;
    }
  }
  return state;
}
void irrp_exg_end(int group_ind){
  int ierr;
  ierr=MPI_Waitall(cpart->nids,cpart->requestr,cpart->statusr);
  ierr=MPI_Waitall(cpart->nids,cpart->requests,cpart->statuss);
  exchange_group_action_recv_unpack(group_ind);
}
void irrp_exg_action(int group_ind,int mode_){
  int mode=7;
  if(mode_>0)mode=mode_;
  if(mode==7){
    exchange_group_action_prepare(group_ind);
    exchange_group_action_send (group_ind);
    exchange_group_action_recvs(group_ind);
    exchange_group_action_recv_unpack(group_ind);
    return;
  }
  if(mode&1){
    exchange_group_action_prepare(group_ind);
    exchange_group_action_send(group_ind);
    exchange_group_action_recv(group_ind);
  }
  if(mode&2){
    int ierr;
    for(int n=0;n<cpart->nids;n++){
      ierr=MPI_Wait(cpart->requestr+n,cpart->statusr+n);
      ierr=MPI_Wait(cpart->requests+n,cpart->statuss+n);
    }
    //ierr=MPI_Waitall(cpart->nids,cpart->requestr,cpart->statusr)
    //ierr=MPI_Waitall(cpart->nids,cpart->requests,cpart->statuss)
  }
  if(mode&4){
    exchange_group_action_recv_unpack(group_ind);
  }
}
#endif

#ifndef NOGATHER
int Verifyroot(int root_){
  int i,ierr;
  int Verifyroot=0;
  if(root_>0){
    Verifyroot=root_;
    if(Verifyroot<0 || Verifyroot>=gsi->npe)Verifyroot=0;
  }
}
void vg_lfvar(int root,int vt,int flg){
  int isr= (gsi->pid==root);
  gsd->itemsize=getitemsize(vt);
  if(!gsd->gsvp&&isr)NEWPN(gsd->gsvp,gsd->itemsize*gsi->gnpc  );
  if(!gsd-> lvp&&flg)NEWPN(gsd-> lvp,gsd->itemsize*cpart->snpc);
}

#define VCPY(VT,dp_,sp_,n,doff,dpp,soff,spp) \
    VT*dp=((VT*)dp_)+doff,*sp=((VT*)sp_)+soff;\
    for(int i=0;i<n;i++) dp[dpp]=sp[spp];

#define VCPYS(dp_,sp_,n,doff,dpp,soff,spp) \
      switch(ivt){\
      case VT_FLOAT :{VCPY(float ,dp_,sp_,n,doff,dpp,soff,spp);break;}\
      case VT_DOUBLE:{VCPY(double,dp_,sp_,n,doff,dpp,soff,spp);break;}\
      case VT_INT   :{VCPY(int   ,dp_,sp_,n,doff,dpp,soff,spp);break;}\
      case VT_SHORT :{VCPY(short,dp_,sp_,n,doff,dpp,soff,spp);break;}\
      }
void irrp_gather(void*lvar,void*gvar,int ivt,int km,int k,int root_){
  //VT,intent(in)       lvar(km,cpart->iblv:cpart->npc);
  //VT,intent(out)      gvar(gsi->ijm);
  int   ierr,root,i;
  root=Verifyroot(root_);
  vg_lfvar(root,ivt,km!=1||gsi->parttype!=1);
  if(km==1&& gsi->parttype==1){
    ierr=MPI_Gatherv(lvar,cpart->snpc,mpidt[ivt],gsd->gsvp,gsi->npcs,gsd->gdispls,mpidt[ivt],root,gsi->mpi_comm);
  }else{
    int is=getitemsize(ivt);
    int npcp=cpart->npc-cpart->iblv+1;

    if(gsi->parttype==1){
      VCPYS(gsd->lvp,lvar,cpart->snpc,0,i,k,i*km);
    }else{
      int*icp=cpart->calpnts;
      VCPYS(gsd->lvp,lvar,cpart->snpc,0,i,k,icp[i]*km);
    }
    ierr=MPI_Gatherv(gsd->lvp,cpart->snpc,mpidt[ivt],gsd->gsvp,gsi->npcs,gsd->gdispls,mpidt[ivt],root,gsi->mpi_comm);
  }
  if(gsi->pid==root){
    int*ipp=gsd->pposg;
    VCPYS(gvar,gsd->gsvp,gsi->gnpc,0,ipp[i],0,i);
  }
}

void irrp_scatter(void*gvar,void*lvar,int ivt,int km,int k,int root_){
  int   ierr;
  int root=Verifyroot(root_);
  vg_lfvar(root,ivt,km!=1||gsi->parttype!=1);
  if(gsi->pid==root){
    int*ipp=gsd->pposg;
    VCPYS(gsd->gsvp,gvar,gsi->gnpc,0,i,0,ipp[i]);
  }
  if(km==1&& gsi->parttype==1){
    ierr=MPI_Scatterv(gsd->gsvp,gsi->npcs,gsd->gdispls,mpidt[ivt],lvar,cpart->snpc,mpidt[ivt],root,gsi->mpi_comm);
  }else{
    ierr=MPI_Scatterv(gsd->gsvp,gsi->npcs,gsd->gdispls,mpidt[ivt],gsd->lvp,cpart->snpc,mpidt[ivt],root,gsi->mpi_comm);
    if(gsi->parttype==1){
      VCPYS(lvar,gsd->lvp,cpart->snpc,k,i*km,0,i);
    }else{
      int*icp=cpart->calpnts;
      VCPYS(lvar,gsd->lvp,cpart->snpc,k,icp[i]*km,0,i);
    }
  }
}
#endif
//-------------------------------------------------------------------------------------------------
void scatter_force_final(force_info*fi){
  fi->nfp=0;fi->infp=0;fi->infpg=0;
  VFREE(fi->fpos   );
  VFREE(fi->fposg  );
  VFREE(fi->infps  );
  VFREE(fi->fdispls);
  VFREE(fi->gsvp   );
  VFREE(fi->lvp    );
  VFREE(fi->idx    );
  VFREE(fi->w      );
}

#ifndef NOSCATTER_EXT
void force_init(force_info*fi){
  int k,ierr;
  NEWZN(fi->infps  ,gsi->npe);
  NEWZN(fi->fdispls,gsi->npe+1);
  ierr=MPI_Gather(&fi->infp,1,MPI_INT,fi->infps,1,MPI_INT,0,gsi->mpi_comm);
  fi->infpg=0;
  for(int k=0;k<gsi->npe;k++)fi->infpg+=fi->infps[k];
  ierr=MPI_Bcast(&fi->infpg,1,MPI_INT,0,gsi->mpi_comm);
  ierr=MPI_Bcast(fi->infps,gsi->npe,MPI_INT,0,gsi->mpi_comm);
  fi->fdispls[0]=0;
  for(int k=0;k<gsi->npe;k++){
    fi->fdispls[k+1]=fi->fdispls[k]+fi->infps[k];
  }
  NEWPN(fi->fposg,fi->infpg+1);
  ierr=MPI_Gatherv(fi->fpos,fi->infp,MPI_INT ,fi->fposg,fi->infps,fi->fdispls,MPI_INT ,0,gsi->mpi_comm);
  ierr=MPI_Bcast(&fi->fposg,(fi->infpg+1),MPI_INT,0,gsi->mpi_comm);
}
void irrp_scatter_ext_init(){
  int   forceid;
  force_info*fi=gsd->scatterext;
  fi->infp=cpart->snp;
  fi->nfp=cpart->snpc;
  NEWPN(fi->fpos,fi->infp+1);
  CPYN (fi->fpos,cpart->ppos,fi->infp+1);
  force_init(fi);
}
void vf_lfvar(force_info*fi,int vt,int flg){
  fi->itemsize=getitemsize(vt);
  if(!fi->gsvp     )NEWPN(fi->gsvp,fi->itemsize*gsi->gnpc);
  if(!fi-> lvp&&flg)NEWPN(fi-> lvp,fi->itemsize*cpart->snpc);
}

void irrp_scatter_ext(void*gvar,void*lvar,int ivt,int km, int k, int root_){
  //VT,intent(in)       gvar(gsi->ijm);
  //VT,intent(out)      lvar(km, cpart->iblv:cpart->np);
  int   idx, i, j, n, st, ierr, root;
  root = Verifyroot(root_);
  force_info*fi =gsd->scatterext;
  if(fi->infp==0)irrp_scatter_ext_init();
  vf_lfvar(fi, ivt, km!=1|| gsi->parttype!=1 );
  if(gsi->pid==root){
    int*ipp=fi->fposg;
    VCPYS(fi->gsvp,gvar,fi->infpg,0,i,0,ipp[i]);
  }
  if(km==1 && gsi->parttype==1){
    ierr=MPI_Scatterv( fi->gsvp,fi->infps,fi->fdispls,mpidt[ivt],lvar,fi->infp,mpidt[ivt],root,gsi->mpi_comm);
  }else{
    ierr=MPI_Scatterv( fi->gsvp,fi->infps,fi->fdispls,mpidt[ivt],fi->lvp,fi->infp,mpidt[ivt],root,gsi->mpi_comm);
    if(gsi->parttype==1){
      VCPYS(lvar,fi->lvp,cpart->snp,k,i*km,0,i);
    }else{
      VCPYS(lvar,fi->lvp,cpart->snp,k,cpart->calpnts[i]*km,0,i);
    }
  }
}
void irrp_scatter_ext_final(){
  scatter_force_final(gsd->scatterext);
}
#endif
//-------------------------------------------------------------------------------------------------
#ifndef NOSCATTER_FORCE
#define VINTER(VT,dp_,sp_,n,doff,dpp) \
    VT*dp=((VT*)dp_)+doff,*sp=((VT*)sp_);\
    for(int i=0;i<n;i++) \
      dp[dpp]=sp[idx[i*4  ]]*w[i*4  ]+sp[idx[i*4+1 ]]*w[i*4+1]+\
              sp[idx[i*4+2]]*w[i*4+2]+ sp[idx[i*4+3]]*w[i*4+3];
#define VINTERS(dp_,sp_,n,doff,dpp) \
      switch(ivt){\
      case VT_FLOAT :{VINTER(float ,dp_,sp_,n,doff,dpp);break;}\
      case VT_DOUBLE:{VINTER(double,dp_,sp_,n,doff,dpp);break;}\
      case VT_INT   :{VINTER(int   ,dp_,sp_,n,doff,dpp);break;}\
      case VT_SHORT :{VINTER(short,dp_,sp_,n,doff,dpp);break;}\
      }
int sf_init(int forceid){
  force_info*fs;
  int mm,i;
  if(forceid<0){
    if(gsd->maxforces>0){
      for(int i=0;i<gsd->maxforces;i++){
        if(gsd->forces[i].nfp==0){
          forceid=i;
          break;
        }
      }
    }
    if(forceid<0)forceid=gsd->maxforces+1;
  }
  if(forceid>gsd->maxforces){
    mm=gsd->maxforces;
    gsd->maxforces=forceid+4;
    NEWPN(gsd->forces,gsd->maxforces);
    if(mm!=0){
      force_info*fs=gsd->forces;
      CPYN(gsd->forces,fs,mm);
      VFREE(fs);
    }
  }
  return forceid;
}

void irrp_scatter_force_final_all(){
  int   forceid;
  irrp_scatter_ext_final();
  for(int i=0;i<gsd->maxforces;i++){
    scatter_force_final(gsd->forces+i);
  }
}

void irrp_scatter_force_final(int forceid){
  if(forceid<0 || forceid>gsd->maxforces)return;
  scatter_force_final(gsd->forces+forceid);
}

void findi1i2(double*fp,int n,double x,int*i1_,int*i2_,double*p,double fcycle){
  int i1=n,i2=0;
  if(fp[1]>=fp[0]){
    for(int i=0;i<n;i++){ if(x<fp[i]){ i1 = i - 1; break; } }
  }else{
    for(int i=0;i<n;i++){ if(x>fp[i]){ i1 = i - 1; break; } }
  }
  i2=i1+1;
  if(i1<0 || i1==n){
    if(fcycle>1e-5){
      i1 = n; i2 = 1; 
      *p=(x-fp[i1])/(fp[i2]+fcycle-fp[i1]);
    }else{
      if(i1 < 0)i1 = 0;
      i2 = i1; p = 0;
    }
  }else{
    *p=(x-fp[i1])/(fp[i2]-fp[i1]);
  }
  *i1_=i1;*i2_=i2;
}

void InterInit1(force_info*fi,double*fix,double*fiy,double*fox,double*foy,int*fm,int nix,int  niy,int np,double xcycle){
  //double,intent(in)  xcycle,fix(nix), fiy(niy),fox(cpart->iblv:cpart->np), foy(cpart->iblv:cpart->np);
  //int ,intent(out)  fm(nix*niy);
  ZERON(fm,nix*niy);
  ZERON(fi->idx,4*fi->nfp);
  ZERON(fi->w  ,4*fi->nfp);
  for(int n=0;n<np;n++){
    double x,y,p,q;
    int i1,i2,j1,j2,ij;
    if(gsi->parttype==1){
      ij=n;
    }else{
      ij=cpart->calpnts[n];
    }
    x=fox[ij];y=foy[ij];
    findi1i2(fix,nix,x,&i1,&i2,&p,xcycle);
    findi1i2(fiy,niy,y,&j1,&j2,&q,0.);
    ij=i1+j1*nix;fm[ij]=1; fi->idx[n*4+0]=ij; fi->w[n*4+0]=(1-p)*(1-q);
    ij=i2+j1*nix;fm[ij]=1; fi->idx[n*4+1]=ij; fi->w[n*4+1]=   p *(1-q);
    ij=i1+j2*nix;fm[ij]=1; fi->idx[n*4+2]=ij; fi->w[n*4+2]=(1-p)*   q;
    ij=i2+j2*nix;fm[ij]=1; fi->idx[n*4+3]=ij; fi->w[n*4+3]=   p *   q;
  }
}

void InterInit2(force_info*fi,double*fix,double* fiy,double* fox,double* foy,int*fm,int nix,int  niy,int np,double xcycle){
  //double,intent(in)  xcycle, fix(nix, niy), fiy(niy, niy);
  //double,intent(in)  fox(cpart->iblv:cpart->np), foy(cpart->iblv:cpart->np);
  //int ,intent(out)  fm(nix*niy);
  //Fix Me
  printf("Orthogonal grid Force interpolation,Not surport Now\n");
  ZERON(fm,nix*niy);
  ZERON(fi->idx,4*fi->nfp);
  ZERON(fi->w  ,4*fi->nfp);
}


void irrp_scatter_force_init(int*forceid_,double*fix_,double*fiy_,double*fox,double*foy,
                              int nix_,int niy_,int nnx_,int nny_,int npflg,int xcycle,int root_){
  //int , intent(inout)   forceid;       // Index of forcing variable.
  //int , intent(in)      nix_, niy_;    // Size of data matrix.
  //int , intent(in)      npflg;         // flag for pnts need to prepare forcing.
  //                                     // 0 for computer pnts only, else for all pnts.
  //int , intent(in)      nnx_, nny_;    // size of input coordinate data of forcing.
  //double, intent(in)      fix_(nnx_);    // x-coordinate of input forcing. For curvlinear grid, nnx=nix*niy
  //double, intent(in)      fiy_(nny_);    // y-coordinate of input forcing. For curvlinear grid, nnx=nix*niy
  //double, intent(in)      xcycle;        // The value of cycled in x direction, i.e. 360.
  //double, intent(in)      fox(cpart->iblv:cpart->np);  // The x-coordinate of model.
  //double, intent(in)      foy(cpart->iblv:cpart->np);  // The y-coordinate of model.
  //int , intent(in), optional   root_;  // Root PE to scatter data.
  int ierr,nix,niy,nnx,nny,nxy;
  int*fm,*ij2s;
  double*fix,*fiy;
  int root = Verifyroot(root_);
  int forceid =-1;
  if(forceid_)forceid=*forceid_;
  *forceid_=forceid=sf_init(forceid);
  if(forceid_)*forceid_=forceid;
  force_info*fi = gsd->forces+forceid;
  if(npflg==0){
    fi->nfp = cpart->snpc;
  }else{
    fi->nfp = cpart->snp;
  }
  RENEWN(fi->idx,4*fi->nfp);
  RENEWN(fi->w  ,4*fi->nfp);
  mpipacket   pk;
  if(gsi->pid==root){
    InitMpiPacket(&pk,(4+nnx_+nny_+100)*8);  // Use a larger number for safe.
    ierr=MPI_Pack(&nix_,    1, MPI_INT   , pk.buf, pk.bsize, &pk.pos, gsi->mpi_comm);
    ierr=MPI_Pack(&niy_,    1, MPI_INT   , pk.buf, pk.bsize, &pk.pos, gsi->mpi_comm);
    ierr=MPI_Pack(&nnx_,    1, MPI_INT   , pk.buf, pk.bsize, &pk.pos, gsi->mpi_comm);
    ierr=MPI_Pack(&nny_,    1, MPI_INT   , pk.buf, pk.bsize, &pk.pos, gsi->mpi_comm);
    ierr=MPI_Pack(fix_,nnx_, MPI_DOUBLE, pk.buf, pk.bsize, &pk.pos, gsi->mpi_comm);
    ierr=MPI_Pack(fiy_,nny_, MPI_DOUBLE, pk.buf, pk.bsize, &pk.pos, gsi->mpi_comm);
  }
  bcast_packet(&pk, root, gsi->pid, gsi->mpi_comm);
  ierr=MPI_Unpack(pk.buf, pk.dsize, &pk.pos, &nix,   1, MPI_INT, gsi->mpi_comm);
  ierr=MPI_Unpack(pk.buf, pk.dsize, &pk.pos, &niy,   1, MPI_INT, gsi->mpi_comm);
  ierr=MPI_Unpack(pk.buf, pk.dsize, &pk.pos, &nnx,   1, MPI_INT, gsi->mpi_comm);
  ierr=MPI_Unpack(pk.buf, pk.dsize, &pk.pos, &nny,   1, MPI_INT, gsi->mpi_comm);
fm=NULL;
  nxy=nix*niy;
  NEWPN(fix,nnx);
  NEWPN(fiy,nny);
  NEWPN(ij2s,nxy);
  NEWZN(fm,nxy);
  ierr=MPI_Unpack(pk.buf, pk.dsize, &pk.pos, fix, nnx, MPI_DOUBLE   , gsi->mpi_comm);
  ierr=MPI_Unpack(pk.buf, pk.dsize, &pk.pos, fiy, nny, MPI_DOUBLE   , gsi->mpi_comm);
  InitMpiPacket(&pk,0);
  fi->nx = nix; fi->ny = niy; fi->nxy = nix*niy;
  if(nnx==nix){
    InterInit1(fi, fix, fiy, fox, foy, fm, nix, niy, fi->nfp, xcycle);
  }else{
    InterInit2(fi, fix, fiy, fox, foy, fm, nix, niy, fi->nfp, xcycle);
  }
  fi->infp=0;for(int k=0;k<nxy;k++)fi->infp+=fm[k];
  NEWZN(fi->fpos,fi->infp+1);
  int idx = 0;
  for(int ij = 0;ij< nxy;ij++){
    if(fm[ij] != 0){
      idx = idx + 1; ij2s[ij] = idx;
      fi->fpos[idx] = ij;
    }
  }
  for(int n = 0;n< fi->nfp*4;n++){
    int st = fi->idx[ n]; 
    fi->idx[n] = ij2s[st];
  }
  VFREE(fm ); 
  VFREE(ij2s);
  VFREE(fix); VFREE(fiy);
  force_init(fi);
}

force_info*GetForceInfo(int forceid){
  if(forceid<1||forceid>gsd->maxforces){
    printf("the gsd->forces Not Inited %d %d\n",forceid,gsd->maxforces);
    irrp_abort(__FILE__, __LINE__,"");
  }
  force_info*fi=gsd->forces+forceid;
  if(fi->nfp<=0){
    printf("the force Not Inited %d %d\n",forceid,gsd->maxforces);
    irrp_abort(__FILE__, __LINE__,"");
  }
  return fi;
}

void irrp_scatter_force(int fid,void*gvar,void*lvar,int ivt,int km,int k,int root_){
  //VT     , intent( in) gvar[fi->nxy];
  //VT     , intent(out) lvar(km,cpart->iblv:cpart->np);
  int i,n,st,ierr,ij,root;
  root=Verifyroot(root_);
  force_info*fi=&gsd->forces[fid];
  vf_lfvar(fi, ivt, 1 );
  if(gsi->pid==root){
    int*ipp=fi->fposg;
    VCPYS(fi->gsvp,gvar,fi->infpg,0,i,0,ipp[i]);
  }
  void*lv=fi->lvp;
  ierr=MPI_Scatterv(fi->gsvp,fi->infps,fi->fdispls,mpidt[ivt],lv,fi->infp,mpidt[ivt],root,gsi->mpi_comm);
  int*idx=fi->idx;
  double*w=fi->w;
  if(gsi->parttype==1){
    VINTERS(lvar,fi->lvp,fi->nfp,k,i*km);
  }else{
    int*ipp=cpart->calpnts;
    VINTERS(lvar,fi->lvp,fi->nfp,k,ipp[i]*km);
  }
}
#undef VINTER
#undef VINTERS
#endif
//Fortran Interface
#define _D(p,v) ((p)?*(p):v)
#define D1(p) ((p)?*(p):1)
#define D0(p) ((p)?*(p):0)
#define E0(p) ((p)?*(p)-1:0)
int setgpid_(int *gpid){
  return setgpid(_D(gpid,-1));
}
int getnpe_(){ return gsi->npe; }
int getnp_(){ return cpart->snp; }
void getnxny_(int*nx,int*ny){ *nx=cpart->nx; *ny=cpart->ny;  }
void irrp_output_pposg_(){irrp_output_pposg();}
void irrp_initt_(int*partmode,int*npe,int*pid,int*mpi_comm,int*mask,int*im,int*jm,int*halosize,int*cycle_flag,int*scycle,int*NS3){
  irrp_init(*partmode,*npe,*pid,MPICOMMF2c(*mpi_comm),mask,*im,*jm,*halosize,*cycle_flag,scycle,*NS3);
}
void irrp_init_data_(){ irrp_init_data();}
void irrp_exginf_(int*nst){irrp_exginf(*nst);}
void irrp_setpartmatrixt_(Rect*recto,int*ptype){
  irrp_SetPartMatrix(recto,ptype);
}
void  irrp_setpartserialt_(int*gnpc,int*npc,int*np,PI_Pos*plist,PI_NB8*nb8,int*iblv_,int*nwps_){
  irrp_SetPartSerial(gnpc,npc,np,plist,nb8,iblv_,nwps_);
}
void  irrp_getrectst_(Rect*recti,Rect*recto,Rect*rects_){
  irrp_getrects(recti,recto,rects_);
}
void irrp_final_(int*mode){ irrp_final(*mode);}
// Initialize exchage functions
void irrp_exg_init_(int*group_ind_,int*nvar) { irrp_exg_init(group_ind_,*nvar) ; }
// Act the exchange by group.
void irrp_exg_action_(int*group_ind,int*mode){ irrp_exg_action(*group_ind,*mode); }
// Finalize the exchange functions.
void irrp_exg_final_(int*group_ind){ irrp_exg_final(*group_ind);}
// start the exchange by group,equal irrp_exg_action(1).
void irrp_exg_start_(int*group_ind){ irrp_exg_start(*group_ind);}
// test mpi irecv  is end ,
// in some realizations of MPI standard ,big data packet
// completed at wait or test
int irrp_exg_check_(int*group_ind){ irrp_exg_check(*group_ind);}
// end the exchange by group,equal irrp_exg_action(6).
void irrp_exg_end_(int*group_ind){ irrp_exg_end(*group_ind);}
// Initialize the function of forcing scatter
void irrp_scatter_force_init_(int*forceid,double*fix,double*fiy,double*fox,double*foy,
                              int*nix,int*niy,int*nnx,int*nny,int*npflg,int*xcycle,int*root){
  irrp_scatter_force_init(forceid,fix,fiy,fox,foy,*nix,*niy,*nnx,*nny,*npflg,*xcycle,D0(root));
}
// Finalize the scatter_force Function.
void irrp_scatter_force_final_(int*forceid){irrp_scatter_force_final(*forceid);}
// Initialize scatter with extended pnts.
void irrp_scatter_ext_init_(){ irrp_scatter_ext_init();}
// Finalize scatter_ext Function
void irrp_scatter_ext_final_(){ irrp_scatter_ext_final();}
//::irrp_exg_setvar                 // Set/append var into a special exchange group.
void irrp_exg_setvart_(int*group_ind,int*ivar,char*var,int *ivt,int *km,const char*vname){
  irrp_exg_setvar(group_ind,ivar,vname,var,*ivt,D1(km));
}

// Gathering info/data from all PEs to root PE.
void irrp_gathert_(void*lvar,void*gvar,int*ivt,int*km,int*k,int*root){
  irrp_gather(lvar,gvar,*ivt,D1(km),E0(k),D0(root));
}
// Scatering info/data from root PE to all PEs.
void irrp_scattert_(void*gvar,void*lvar,int*ivt,int*km,int*k,int*root){
  irrp_scatter(gvar,lvar,*ivt,D1(km),E0(k),D0(root));
}
// Scattering the info with extended pnts.
void irrp_scatter_extt_(void*gvar,void*lvar,int*ivt,int*km,int*k,int*root){
  irrp_scatter_ext(gvar,lvar,*ivt,D1(km),E0(k),D0(root));
}
// Act the forcing scattering.
void irrp_scatter_forcet_(int*fid,void*gvar,void*lvar,int*ivt,int*km,int*k,int*root){
  irrp_scatter_force(*fid,gvar,lvar,*ivt,D1(km),E0(k),D0(root));
}
#undef _D
#undef E0
#undef D0
#undef D1
