#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <mpi.h>
#include "irrp_inner.h"
// /usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/openmpi/mpi.h


#define DBGB0(msg)    if(gsi->pid==0)printf("%4d %8.3fs %7d %s\n",gsi->pid,difftimer(0),iwalltime(),msg);
/*################################################################################################
 *         Copyright (C) 2013  Wei Zhao
 *         MODULE NAME : irrp_kernal_mod
 *         PRESENT VERSION : 2013-04-30
 * --- DEPEND: irrp_smpi_mod, irrp_pemask_mod
 * --- NOTE for describing of subroutine / function :
 *  A. The parameters bracketed with [], means optional parameter.
 *  B. The describe for the parameters of subroutine / function, started with:
 *   * It means input prameter;
 *   # It means output prameter;
 *   @ It means input and output prameter(it will be changed inside).
 *
 *-------------------------------------------------------------------------------------------------
 * ***                                 INTERFACE DESCRIBE                                       ***
 *-------------------------------------------------------------------------------------------------
 *
 *  1. sub. irrp_set_pemask : Do irregular partition using MASK and return the results in PEMASK.
 *
 *     irrp_part_init(partmode, npe, pid, mpi_comm, mask, halosize, cycle_flag, scycle)
 *
 *    * int  mask_   = Mask for horizonal space/points in which 0 means land/useless points.
 *    # int  pemask_ = Arrary to save partition reaults each point with the value equal to
 *                           the ID of the PE that this point is belonged. This ID of PEs is
 *                           started from 1 which is different with the way for MPI.
 *    * int  im_     = The first dimension size of mask_ and pemask_.
 *    * int  jm_     = The second dimension size of mask_ and pemask_.
 *    * int  npe_    = The total number of PEs used for this partition.
 *
 *-------------------------------------------------------------------------------------------------
 **************************************************************************************************
 *-----------------------------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------------------------------
 *
 *  public  irrp_part_init, irrp_part_final
 *  public  irrp_SetPartMatrix, irrp_SetPartSerial,irrp_getrects
 *  public  PI_Pos, PI_Rsinfo, PI_NB8
 *  public  PI_Ginfo,PI_Part,PI_Rsbuf,PI_Rplist
 *  public  gsi, cpart
 */
static int* scycle,NS3;
#ifdef SERIAL_TEST
static PI_Part*  parts=NULL;
#endif
static void  init_part();
void irrp_set_pemask(int *mask_,int*pemask_, int im_, int jm_, int npe_,int pid_);
static void  CheckBalance();
static void  set_part_ij2s(int *scycle,int NS3);
static void  get_nbp_recv(int m);
static void  SetPnts(int m);
static void  SetPartInf(int m);

static PI_Part*  part=NULL;
static void init_gsi(PI_Ginfo*pgsi){
  pgsi->mpi_comm=MPI_COMM_WORLD;
  pgsi->halosize=1;
}
// getrects
void  irrp_getrects(Rect*recti,Rect*recto,Rect*rects){
  if(recti)*recti= part->recti;
  if(recto)*recto= part->recto;
  if(rects){
    Rect*rts=part->pgsi->rects;
    for(int j=0;j<gsi->npe;j++){
      rects[j]=rts[j];
    }
  }
}
//points index convert:
//(GI,GJ) (i,j) in global
//(PI,PJ) (i,j) in Part
//GSIJ  serial num in global
//PSIJ  serial num in part
//PN   sorted serial num in part
//GSIJ TO (GI,GJ)
static void  gsij2gij(int gsij,int *gi,int *gj){
  *gi=gsij%gsi->im;
  *gj=gsij/gsi->im;
}
//(GI,GJ) TO gsij
static int sgij2gsij(int gi,int gj){
  return gi+gj*gsi->im;
}
static void adjgij(int *gi,int *gj){
  if(gsi->cycle_flag==3){
    //FIXME:
    return;
  }
  if(gsi->cycle_flag>0){
    if(*gi<0)*gi+=gsi->im;
    if(*gi>=gsi->im)*gi-=gsi->im;
  }else{
    if(*gi<0)*gi=0;
    if(*gi>=gsi->im)*gi=gsi->im-1;
  }
  if(gsi->cycle_flag>1){
    if(*gj<0)*gj+=gsi->jm;
    if(*gj>=gsi->jm)*gj-=gsi->jm;
  }else{
    if(*gj<0)*gj=0;
    if(*gj>=gsi->jm)*gj=gsi->jm-1;
  }
}
static int gij2gsij(int gi,int gj){
  adjgij(&gi,&gj);
  return gi+gj*gsi->im;
}
//(PI,PJ) to psij
static int pij2psij(int pi,int pj){
  int psij=pi+pj*part->nx;
  if(pi<0||pi>=part->nx||pj<0||pj>part->ny)return -1;
  return psij;
}
//(GI,GJ) to psij
static int gij2psij(int gi,int gj){
  return pij2psij(gi-part->ib,gj-part->jb);
}
//GISI to psij
static int gsij2psij(int gsij){
  int gi,gj;
  gsij2gij(gsij,&gi,&gj);
  return gij2psij(gi,gj);
}
//PN (gij) to psij
static int png2psij(int pn){
  int gij=part->ppos[pn];
  return gsij2psij(gij);
}
//(PI,PJ) to GIJ
static int pij2gsij(int i,int j){
  return gij2gsij(i+part->ib,j+part->jb);
}

/*------------------------------------------------------------------------------------------------
 * init Part:
 * input:
 *  partmode: 0:Normal part;other part method,Input pemask
 *  npe: number of PE
 *  pid: cur pid
 *  mpi_comm:mpi comm
 *    mask_:Mask for horizonal space/points in which 0 means land/useless points.????
 * halosize:
 * cycle_flag: Open boundary type.
 *   0 - None cycle boundary.
 *   1 - Cycled boundary in x-direction for (normal/displaced) polar grid.
 *   2 - Cycled boundary in x-direction and y-direction for triple polar grid.
 *   3 - Cycled boundaries are given manually.
 * scycle:
 *-----------------------------------------------------------------------------------------------*/
void starttimer(int iid);
void endtimer(int iid);
double difftimer(int iid);
void irrp_part_init(int partmode,int  npe,int  pid,MPI_Comm mpi_comm,
      int*  mask,int*  pemask,int im,int jm,int halosize,int cycle_flag,int* scycle_,int NS3_){
  int  i, j, ij, nn, m;
  gsi->npe      = npe;                      // Number of PEs.
  gsi->pid      = pid;                      // Current ID of PE, start from 0.
  gsi->mpi_comm = mpi_comm;                 // Common world of MPI.
  gsi->im       = im;                       // Size of grid matrix: first dimension.
  gsi->jm       = jm;                       // Size of grid matrix: second dimension.
  gsi->ijm      = gsi->im*gsi->jm;            // Maximum computing points.
  gsi->halosize = 1;                        // Halo size, default is 1.
  starttimer(0);
  if(halosize>0)gsi->halosize = halosize;
  // Record halo size given from outside.
  /*------------------------------------------------------------------------------------------
   * cycle flag is given from outside.
   *   0 - None cycle boundary.
   *   1 - Cycled boundary in x-direction for (normal/displaced) polar grid.
   *   2 - Cycled boundary in x-direction and y-direction for triple polar grid.
   *   3 - Cycled boundaries is given manually.
   *-----------------------------------------------------------------------------------------*/
  scycle=scycle_;
  NS3=NS3_;
  if(cycle_flag>=0){             // Record cycle flag: default is 0.
    gsi->cycle_flag = cycle_flag;
  }else if(scycle){
    gsi->cycle_flag = 3;         // Determined by scycle:
                                 // Cycled boundaries are given manually.
  }else{
    gsi->cycle_flag = 0;         // Default setting: None cycle boundary.
  }
  if(gsi->cycle_flag != 0 && gsi->halosize != 1){
    // Check the matchment of cycle flag and halo size.
    printf("The case of halosize !=1 not surport yet when cycle_flag !=0. \n");
    irrp_abort(__FILE__, __LINE__,"");
  }

  DBGB0("irrp_part_init begin");
  NEWZN(gsi->rects,4*gsi->npe);  // allocate rects
  NEWZN(gsi->npcs,gsi->npe);     // allocate npcs
  //DBGB0("irrp_part_init 1")
  nn = 0;
  for(int i = 0;i< gsi->ijm;i++){// sequentializing points from matrix form.
    if(mask[i] > 0)nn = nn +1;
  }
  //DBGB0("irrp_part_init 2 ")
  gsi->sumnp = nn;                          // Set number of sequentialized points.
  gsi->avenp = gsi->sumnp / (float)gsi->npe;  // Compute averaged number of points for each PE.
  NEWZN(gsi->mask,gsi->ijm);
  if(partmode == 0){                   // Do partion using irrp_pemask_mod
    DBGB0("irrp_set_pemask begin");
    irrp_set_pemask(mask,gsi->mask,  gsi->im, gsi->jm, gsi->npe,gsi->pid);
    gsi->balance = 1;
    if(pemask==NULL){
      memcpy(mask,gsi->mask,sizeof(int)*gsi->ijm);
    }else{
      memcpy(pemask,gsi->mask,sizeof(int)*gsi->ijm);
    }
    DBGB0("irrp_set_pemask end");
  }else{
    gsi->balance = 0;
    memcpy(mask,gsi->mask,sizeof(int)*gsi->ijm);
  }
#ifdef SERIAL_TEST
  NEWZN(parts,gsi->npe);                // allocate parts.
  DBGB0("SetPartInf");
  // Set rects and npcs.Set information for each partition.
  for(int m=0;m<gsi->npe;m++){
    part=parts+m;
    SetPartInf(m+1);
    gsi->rects[m]=part->recti;
    gsi->npcs[m]=part->snpc;
  }
  CheckBalance();
  DBGB0("set_part_ij2s");
  for(int m=0;m<gsi->npe;m++){
    part=parts+m;
    set_part_ij2s(scycle,NS3);
  }
  DBGB0("get_nbp_recv");
  for(int m=0;m<gsi->npe;m++){
    part=parts+m;
    get_nbp_recv(m+1);        // Set recv points for each partition.
  }
  DBGB0("SetPnts");
  for(int m=0;m<gsi->npe;m++){
    part=parts+m;
    SetPnts(m);             // Set all points for each partition
  }
#else
  part=cpart;
  DBGB0("SetPartInf");
  SetPartInf(gsi->pid+1);        // Set information for current partition.
  DBGB0("gather/bcast npc,recti");
  int ierr;
  ierr=MPI_Gather(&part->snpc,1,MPI_INTEGER,gsi->npcs ,1,MPI_INTEGER,0,gsi->mpi_comm);
  ierr=MPI_Gather(&part->recti,4,MPI_INTEGER,gsi->rects,4,MPI_INTEGER,0,gsi->mpi_comm);
  ierr=MPI_Bcast(gsi->npcs,gsi->npe,MPI_INTEGER,0,gsi->mpi_comm);
  ierr=MPI_Bcast(gsi->rects,4*gsi->npe,MPI_INTEGER,0,gsi->mpi_comm);
  CheckBalance();

  DBGB0("set_part_ij2s");
  set_part_ij2s(scycle,NS3);
  get_nbp_recv(gsi->pid+1);     // Set recv points for current partition.
  DBGB0("SetPnts");
  SetPnts(gsi->pid+1);          // Set all points for current partition.
  DBGB0("SetPnts End ");
#endif
  gsi->gnpc = 0;
  for(int i=0;i<gsi->npe;i++)gsi->gnpc+=gsi->npcs[i]; // Number of all computer points for all PEs.
  DBGB0("irrp_SetPartSerial");
  irrp_SetPartSerial(0,0,0,0,0,0,0);                 // Prepare point list in serial case as default.
                                           // It can be changed into matrix case later.
  DBGB0("End irrp_part_init");
}
// Init part struct
static void  init_part(){
  memset(part,0,sizeof(*part));
}
/*----------------------------------------------------------------------------------------------
 * sub. release_part: release buffers and point lists in part.
 *    release_part( mode)
 *   # int  mode = method of release. If not given, only release buffers.
 *---------------------------------------------------------------------------------------------*/
static void  release_part(int mode){
  int  i;
  if(part->rsbufs){
    for(int i = 0;i< part->nids;i++){
      VFREE(part->rsbufs[i].rbuf.pnts);
      VFREE(part->rsbufs[i].sbuf.pnts);
    }
    VFREE(part->rsbufs);
  }
  if(part->rsinfo){
    for(int i = 0;i< part->nids;i++){
      VFREE(part->rsinfo[i].mrpnts);
      VFREE(part->rsinfo[i].mspnts);
      VFREE(part->rsinfo[i].srpnts);
      VFREE(part->rsinfo[i].sspnts);
    }
    VFREE(part->rsinfo  );
  }
  VFREE(part->mpnttype);
  VFREE(part->calpnts  );
  VFREE(part->ijbe    );
  VFREE(part->ij2s     );
  VFREE(part->ppos    );
  VFREE(part->nb8     );
  VFREE(part->requests);
  VFREE(part->requestr);
  VFREE(part->statuss );
  VFREE(part->statusr );
}
/*-------------------------------------------------------------------------------------------------
 *End Part
 *mode: 0:part end ,leave vars for use
 *1:free all var,not use ????*/
void  irrp_part_final(int mode){
#ifdef SERIAL_TEST
  int  m;
  if(parts){
    for(int m=0;m<gsi->npe;m++){
      part=parts+m;
      release_part(mode);
    }
    if(mode )VFREE(parts);
  }
#else
  release_part( mode);          // Finalize cpart->
#endif
  VFREE(gsi->mask  );
  if(mode){
    VFREE(gsi->rects  );
    VFREE(gsi->npcs   );
  }
}
/*----------------------------------------------------------------------------------------------
 * sub. SetPartInf: Prepare rects and npcs for all PEs.
 *---------------------------------------------------------------------------------------------*/;
static void  SetPartInf(int m){
  Rect rect={-1,-1,-1,-1};
  int npc=0;
  init_part();
  for(int gj = 0;gj< gsi->jm;gj++){
    for(int gi = 0;gi< gsi->im;gi++){
      if(m == gsi->mask[sgij2gsij(gi,gj)]){      // Get ID for this point.
        if(rect.x0 < 0){
          rect.x0 = rect.x1 = gi;
          rect.y0 = rect.y1 = gj;
        }else{                                // Adjust for minimum and maximum gi/gj.
          if(rect.x0 > gi)rect.x0 = gi;
          if(rect.x1 < gi)rect.x1 = gi;
          if(rect.y0 > gj)rect.y0 = gj;
          if(rect.y1 < gj)rect.y1 = gj;
        }
        npc++;;         // Accumulation for npcs.
      }
    }
  }
  part->snpc =npc; part->recti=rect;
  rect.x0 -= gsi->halosize; rect.y0 -= gsi->halosize;
  rect.x1 += gsi->halosize; rect.y1 += gsi->halosize;
  if(gsi->cycle_flag<1){
    if(rect.x0 <0)rect.x0=0;
    if(rect.x1 >=gsi->im)rect.x1 =gsi->im-1;
  }
  if(gsi->cycle_flag<2){
    if(rect.y0 <0)rect.y0=0;
    if(rect.y1 >=gsi->jm)rect.y1 =gsi->jm-1;
  }
  part->recto=rect;
  part->ib = rect.x0; part->ie = rect.x1;
  part->jb = rect.y0; part->je = rect.y1;
  part->nx = rect.x1 - rect.x0 +1;
  part->ny = rect.y1 - rect.y0 +1;
  part->nxy = part->nx*part->ny;
}
static void  CheckBalance(){
  /*--------------------------------------------------------------------------------------------
  // Check balance and points omitted.
  // Only the first PE do the following checks.
  //-------------------------------------------------------------------------------------------*/
  if(gsi->pid == 0){
    if(gsi->balance == 1){               // This check is only for irregular partition.
      for(int n=0;n<gsi->npe;n++){       // Check absolute balance for each PE.
        if(abs(gsi->npcs[n] - gsi->avenp) > 1){
          printf("Error B Not Absolute balance %d %d %d %8.2f %dnwpcs:\n", n,gsi->npe, gsi->npcs[n], gsi->avenp,gsi->sumnp);
          int sn=0;
          for(int i=0;i<gsi->npe;i++){
            sn+=gsi->npcs[i];printf(" %d",gsi->npcs[i]);
          }
          printf(" %d\n",sn);
          irrp_abort(__FILE__, __LINE__,"");
        }
      }
    }
  }
}

static void  set_part_ij2s(int *scycle,int NS3){ //scycle;  [4*N]
  int  ibi,iei,jbi,jei;
  NEWPN(part->ij2s, part->nxy);
  memset(part->ij2s,-1,part->nxy);//default -1
  ibi=part->ib; iei=part->ie;
  jbi=part->jb; jei=part->je;
  for(int gj = jbi;gj<= jei;gj++){
    for(int gi = ibi;gi<=iei;gi++){
      int pij=gij2psij(gi,gj);
      int gij=gij2gsij(gi,gj);
      if(gsi->mask[gij]>0)part->ij2s[pij]=gij;
    }
  }
  if(gsi->cycle_flag == 3){
    if( scycle==NULL){
      printf( "The scycle is not given when cycle_flag == 3. \n");
      irrp_abort(__FILE__, __LINE__,"");
    }
    int *sc=scycle;
    for(int n=0;n<NS3;n++,sc+=4){
      int pij=gij2psij(sc[0],sc[1]);
      if(pij>0){
        int gij0=gij2gsij(sc[2],sc[3]);
        if(gsi->mask[gij0]>0)part->ij2s[pij]=gij0;
      }
    }
  }
}
/*----------------------------------------------------------------------------------------------
 * sub. set_recv_rplist: prepare the list of the receiving points for part.
 *    set_recv_rplist(m)
 *   * int  m    = ID of current PE, this number is begin from 1.
 *---------------------------------------------------------------------------------------------*/
static int set_recv_rplist(int m,PI_Rplist*rplist){
  double*  mt ;
  int*  flgs ;
  int *lpe,*rpe,*upe,*dpe;
  int ib,ie,jb,je;
  ib=part->recti.x0; ie=part->recti.x1;
  jb=part->recti.y0; je=part->recti.y1;
  NEWZN(upe,gsi->npe+1);
  NEWZN(dpe,gsi->npe+1);
  NEWZN(lpe,gsi->npe+1);
  NEWZN(rpe,gsi->npe+1);
  NEWZN(mt,part->nxy);
  NEWZN(flgs,part->nxy);
  for(int pj = 0;pj< part->ny;pj++){
    for(int pi = 0;pi< part->nx;pi++){
      int gsij=pij2gsij(pi,pj);
      if(gsi->mask[gsij] == m){
        int psij=pij2psij(pi,pj);
        mt[psij] = -1e30;        // The point belong to m-th PE is set as nearly infinite.
      }
    }
  }
  for(int pj = 0;pj< part->ny;pj++){
    for(int pi = 0;pi< part->nx;pi++){
      int psij=pij2psij(pi,pj);
      int gsij=pij2gsij(pi,pj);
      int is = part->ij2s[psij];//ZZ
      if(is >=0 && gsi->mask[gsij] == m){
        double  dm=1;
        for(int ih = gsi->halosize;ih>0;ih-- ){     // For each circle of surrounding points.
          for(int j1 = pj-ih;j1<= pj+ih;j1++){      // For the receiving points neighbored with inner points,
            for(int i1 = pi-ih;i1<= pi+ih;i1++){
              int gsij=pij2gsij(pi,pj);
              int mc=gsi->mask[gsij];
              if( mc== m)continue;
              int psij=pij2psij(i1,j1);
              if(psij<0)continue;
              mt[psij] +=  dm;
              int flg=0;
              /* Low                      High
               * 66EEEEEEEAA           * 01111111110
               * 6---------EEEEEEAA    * 2---------31111110
               * 7----------------A    * 2----------------2
               * 5----------------B    * 2----------------2
               * 57---------------B    * 03---------------2
               *  7---------------BA   *  2---------------32
               *  7----------------A   *  2----------------2
               *  7----------------B   *  2----------------2
               *  55DD-------------9   *  0113-------------2
               *     55DDDDDDDDDDD99   *     011111111111110
               *
               * */
              //1:down,2:up,4:left,8:right
              if(i1<pi)flg|=4; else if(i1>pi)flg|=8;else flg|=16;
              if(j1<pj)flg|=1; else if(j1>pj)flg|=2;else flg|=32;
              flgs[psij]|=flg;
            }
          }
        }
      }
    }
  }
  int npr = 0;
  for(int pj = 0;pj< part->ny;pj++){
    for(int pi = 0;pi< part->nx;pi++){
      int psij=pij2psij(pi,pj);
      int gsij=pij2gsij(pi,pj);
      int is = part->ij2s[psij];                   // part->ij2s(pi, pj) == 0 means land point.
      if(is ==0)continue;
      if(gsi->mask[gsij] == m)continue;
      if(gsi->mask[gsij] == 0)continue;            // Z must not land
      int flg=flgs[psij];
      if(flg > 0){
        PI_Rplist*rp=&rplist[npr];
        int it,jt,pe;
        gsij2gij(is,&it,&jt);
        rp->gis = is;
        rp->i   = it;
        rp->j   = jt;
        rp->pe  = pe=gsi->mask[is];
        if(jt<=jb)dpe[pe]++;
        else if(jt>=je)upe[pe]++;
        if(it<=ib)lpe[pe]++;
        else if(it>=ie)rpe[pe]++;
        //2:down,1:up,8:left,4:right
        /*
         *  >>>>>>>>>>>>>>>>>>
         *  99111111155            * 32222222223
         * ^9---------11111155     * 1---------02222223
         * ^8----------------5 ^   * 1----------------1
         * ^A----------------4 ^   * 1----------------1
         * ^A8---------------4 ^   * 30---------------1
         * ^ 8---------------45^   *  1---------------01
         * ^ 8----------------5^   *  1----------------1
         * ^ 8----------------4^   *  1----------------1
         *   AA22-------------6^   *  0220-------------1
         *      AA2222222222226    *     322222222222223
         *   >>>>>>>>>>>>>>>>>>
         */
        rp->flg = (~flg);
        npr = npr +1;
      }
    }
  }
  for(int ip=1;ip<=gsi->npe;ip++){
    if(dpe[ip]&&dpe[ip]<=2 ){
      if(ip>m)dpe[ip]=0;
      else if(lpe[ip]||rpe[ip])dpe[ip]=0;//fixme :if part width<5 may error
      else {
        if(lpe[ip])lpe[ip]=0;
        if(rpe[ip])rpe[ip]=0;
      }
    }
    if(upe[ip]&&upe[ip]<=2 ){
      if(ip<m)upe[ip]=0;
      else if(lpe[ip]||rpe[ip])upe[ip]=0;//fixme :if part width<5 may error
      else {
        if(lpe[ip])lpe[ip]=0;
        if(rpe[ip])rpe[ip]=0;
      }
    }
  }
  //sort_rplist(npr);
  {
    int *  seqs;
    NEWZN(seqs,npr);
    int  lf,nna,pe,pc,pu,pd,j;
    int  idxn, pidc;
    // zhaowei 2025-7-26
    // big enough for not 
    nna=(gsi->im+gsi->jm)*16;
    for(int ii=0;ii<npr;ii++){
      PI_Rplist*rp=&rplist[npr];
      int i=rp->i;
      int j=rp->j;
      int pe=rp->pe;
      int flg=rp->flg;
      //assume max 8 holo
      //2:down,1:up,8:left,4:right
      /*
       * +---+---+---+---+---+---+     *  >>>>>>>>>>>>>>>>>>
       * |D10|D12| 10| 12|A10|A12|     *  99111111155            * 32222222223
       * +---+---+---+---+---+---+     * ^9---------11111155     * 1---------02222223
       * |D9 |D11| 9 | 11|A9 |A11|     * ^8----------------5 ^   * 1----------------1
       * +---+---+---+---+---+---+     * ^A----------------4 ^   * 1----------------1
       * |D6 |D8 | 6 | 8 |A6 |A8 |     * ^A8---------------4 ^   * 30---------------1
       * +---+---+---+---+---+---+     * ^ 8---------------45^   *  1---------------01
       * |D5 |D7 | 5 | 7 |A5 |A7 |     * ^ 8----------------5^   *  1----------------1
       * +---+---+---+---+---+---+     * ^ 8----------------4^   *  1----------------1
       * |D2 |D4 | 2 | 4 |A2 |A4 |     *   AA22-------------6^   *  0220-------------1
       * +---+---+---+---+---+---+     *      AA2222222222226    *     322222222222223
       * |D1 |D3 | 1 | 3 |A1 |A3 |     *   >>>>>>>>>>>>>>>>>>
       * +---+---+---+---+---+---+
       * down:  &32; 32 12 x0 | 02 22 32
       */
      if     (dpe[pe]) lf=nna*0+i*8-j; // down boundart 1st ,left to right,if i is same,up to down
      else if(upe[pe]) lf=nna*3+i*8-j; // up boundart last ,left to right,if i is same,up to down
      else if(rpe[pe]) lf=nna*1+j*8-i; // right boundart 2nd ,down to up,if j is same,right to left
      else if(lpe[pe]) lf=nna*2+j*8-i; // left boundart 3rd ,down to up,if j is same,right to left
      seqs[ii]=lf;
    }

    for(int exchange=1;exchange;){
      exchange=0;
      for(int ii=1;ii<npr;ii++){
        if(seqs[ii-1]>seqs[ii]){
          int it=seqs[ii-1];        seqs[ii-1]=seqs[ii];              seqs[ii]=it;
          PI_Rplist rp=rplist[ii-1];rplist[ii-1]=rplist[ii];rplist[ii]=rp;
          exchange=1;
        }
      }
    }
    VFREE(seqs);
  }
  VFREE(upe); VFREE(dpe); VFREE(lpe); VFREE(rpe);
  VFREE(flgs);VFREE(mt ); 
  return npr;
}
static void  get_nbp_recv(int m){
  int * nbps, *nnbps, *rpind;
  PI_Rplist*rplist=NULL;
  NEWZN(rplist,part->nxy);
  int npr=set_recv_rplist(m,rplist);
  part->npr = npr;
  part->snp  = part->snpc + part->npr;
  NEWZN(nbps,gsi->npe);
  NEWZN(nnbps,gsi->npe);
  NEWZN(rpind,npr);
  int idxn = 0;
  //count neighbor part & rpnts
  for(int i = 0;i< npr;i++){
    int pidc = rplist[i].pe;
    int j;
    for(j = idxn-1 ;j>=0; j--){
      if(pidc==nbps[j])break;
    }
    if(j < 0)nbps[j= idxn++]=pidc;//add new neighbor
    nnbps[j] = nnbps[j] +1;//add pnts 
    rpind[i] = j;
  }
  part->nids = idxn;
  NEWZN(part->rsinfo,part->nids);
  for(int i = 0;i< part->nids;i++){
    part->rsinfo[i].id = nbps[i] - 1;//part index from 0
    part->rsinfo[i].rn = nnbps[i];
    NEWZN(part->rsinfo[i].srpnts,nnbps[i]);//alloc pnts buff
  }
  ZERON(nnbps,gsi->npe);
  for(int i = 0;i< npr;i++){  //set pnts buff
    int j= rpind[i],k= nnbps[j]++;
    part->rsinfo[j].srpnts[k] = rplist[i].gis;
  }
  VFREE(rpind);
  VFREE(nnbps);
  VFREE(nbps);
  VFREE(rplist);
}

#ifndef SERIAL_TEST
static void  exchange(int m ){
  int i, ierr;
  MPI_Request *requests, *requestr;
  MPI_Status*statusr,*statuss;
  NEWZN(requestr,part->nids); NEWZN(statusr, part->nids);
  NEWZN(requests,part->nids); NEWZN(statuss, part->nids);
  PI_Rsbuf *prsb=part->rsbufs;
  for(int i = 0;i< part->nids;i++){
    ierr=MPI_Isend(&prsb[i].sbuf, 2, MPI_INT,prsb[i].sbuf.id, 10000, gsi->mpi_comm, requests+i);
  }
  DBGB0("irrp exchange 1");
  for(int i = 0;i< part->nids;i++){
    ierr=MPI_Irecv(&prsb[i].rbuf, 2, MPI_INT , part->rsinfo[i].id, 10000, gsi->mpi_comm, requestr+i);
  }
  DBGB0("irrp exchange 2");

  ierr=MPI_Waitall(part->nids, requestr, statusr);
  ierr=MPI_Waitall(part->nids, requests, statuss);
  for(int i = 0;i< part->nids;i++){
    ierr=MPI_Isend(&prsb[i].sbuf.pnts,prsb[i].sbuf.rn, MPI_INT,prsb[i].sbuf.id, 10000, gsi->mpi_comm, requests+i);
  }
  DBGB0("irrp exchange 3");
  for(int i = 0;i< part->nids;i++){
    NEWZN(prsb[i].rbuf.pnts,prsb[i].rbuf.rn);
    ierr=MPI_Irecv(&prsb[i].rbuf.pnts,prsb[i].rbuf.rn, MPI_INT , prsb[i].rbuf.id, 10000, gsi->mpi_comm, requestr+i);
  }
  DBGB0("irrp exchange 4");

  ierr=MPI_Waitall(part->nids, requestr, statusr);
  ierr=MPI_Waitall(part->nids, requests, statuss);

  DBGB0("irrp exchange 5");
  VFREE(requestr); VFREE(statusr);
  VFREE(requests); VFREE(statuss);
}

static void  Getnbrecv(int m){
  PI_Rsbuf *prsb=part->rsbufs;
  for(int i = 0;i< part->nids;i++){
    prsb[i].sbuf.id = m;
    prsb[i].sbuf.rn = part->rsinfo[i].rn;
    NEWZN(prsb[i].sbuf.pnts,part->rsinfo[i].rn);
    for(int k = 0;k< part->rsinfo[i].rn;k++){
      prsb[i].sbuf.pnts[k] = part->rsinfo[i].srpnts[k];
    }
  }
  exchange(m);
}
#else
static void  Getnbrecv(int m){
  for(int i = 0;i< part->nids;i++){
    int m1 = part->rsinfo[i].id ;
    part->rsbufs[i].rbuf.id = m1;
    for(int j = 0;j< parts[m1].nids;j++){
      if(parts[m1].rsinfo[i].id == m){
        part->rsbufs[i].rbuf.id = m1;
        part->rsbufs[i].rbuf.rn = parts[m1].rsinfo[i].rn;
        NEWZN(prsb[i].rbuf.pnts,prsb[i].rbuf.rn);
        for(int k = 0;k< parts[m1].rsinfo[j].rn;k++){
          part->rsbufs[i].rbuf.pnts[k] = parts[m1].rsinfo[j].srpnts[k];
        }
        break;
      }
    }
  }
}
#endif

static void  SetPnts(int m){
  int *  iplist;     // (npr)
  int *  aplist;     // (npr)
  int *  g2l;        // (npr)
  int *  ij2s;
  int *ptf;
  //DBGB0("SetPnts1")
  NEWZN(part->rsbufs,part->nids);
  Getnbrecv(m);
  //DBGB0("SetPnts2")
  for(int i = 0;i< part->nids;i++){
    PI_Rsbuf*prsb=&part->rsbufs[i];
    part->rsinfo[i].sn    = prsb->rbuf.rn;
    NEWZN(part->rsinfo[i].sspnts,prsb->rbuf.rn);
    memcpy(part->rsinfo[i].sspnts,prsb->rbuf.pnts,prsb->rbuf.rn);
  }
  for(int i = 0;i< part->nids;i++){
    VFREE(part->rsbufs[i].rbuf.pnts);
    VFREE(part->rsbufs[i].sbuf.pnts);
  }
  VFREE(part->rsbufs);

  //DBGB0("SetPnts3")
  part->snp = part->snpc + part->npr;
  NEWZN(aplist,part->snp+1);
  NEWZN(part->ppos,part->snp+1);
  part->ppos[0]=0;
  NEWZN(iplist,part->snpc+1);
  NEWZN(ptf,part->nxy+1);
  NEWZN(g2l,part->nxy+1);
  int npidxm=0;
  DBGB0("SetPnts4")
  //cal pnt table ,iplist,g2l
  for(int j = 0;j< part->ny;j++){
    for(int i = 0;i< part->nx;i++){
      int gsij=pij2gsij(i,j);
      if(gsi->mask[gsij] == m){
        int psij=pij2psij(i,j);
        ptf[psij]=gsij;
        npidxm = npidxm +1;
        iplist[npidxm] = gsij;
        g2l[gsij] = npidxm;
      }
    }
  }
  //printf("%d %d %d %d \n",npidxm,part->snpc,part->ib,part->jb);
  DBGB0("SetPnts5")
  // Sort Send
  npidxm=0;
  // Sort Send pnts,add 
  for(int i = 0;i< part->nids;i++){
    for(int j=0;j<part->rsinfo[i].sn;j++){
      int is = part->rsinfo[i].sspnts[j];
      int psij=gsij2psij(is);
      int ls = g2l[is];
      if(iplist[ls]!=0){// send pnts maybe dump
        npidxm = npidxm +1;
        aplist[npidxm] = is;
        iplist[ls] = 0;//avoid dump sendpnts
      }
    }
  }
  //DBGB0("SetPnts6")
  part->nps = npidxm;
  // For inner pnts
  for(int ls = 1;ls<= part->snpc;ls++){
    if(iplist[ls] != 0){
      npidxm = npidxm +1;
      aplist[npidxm] = iplist[ls];
      iplist[ls] = 0;
    }
  }
  VFREE(iplist);
  // Sort recv pnts
  //DBGB0("SetPnts7")
  for(int i = 0;i< part->nids;i++){
    //part->rsinfo[i].re = npidxm +1
    for(int j=0;j<part->rsinfo[i].rn;j++){
      npidxm = npidxm +1;
      aplist[npidxm] = part->rsinfo[i].srpnts[j];
    }
    //part->rsinfo[i].re = npidxm
  }
  //DBGB0("SetPnts8")
  NEWZN(part->nb8,part->snp+1);
  part->nb8[0]=(PI_NB8){0,0,0,0,0,0,0,0};
  memset(g2l,0,part->nxy);
  ij2s = part->ij2s;
  //DBGB0("SetPnts9")
  //save ppos
  for(int ls = 1;ls<= part->snp;ls++){
    int is = aplist[ls];
    g2l[is] = ls; part->ppos[ls] = is;
  }
  //DBGB0("SetPnts10")
  int nx=gsi->im;
  for(int ls = 1;ls<= part->snpc;ls++){
    int is = aplist[ls]     ; 
    PI_NB8* nb8 = &part->nb8[ls];//ZZ
    nb8->ul = g2l[is-1+nx]; nb8->u = g2l[is-1+nx]; nb8->ur = g2l[is+1+nx];
    nb8-> l = g2l[is-1   ];                        nb8-> r = g2l[is+1   ];
    nb8->dl = g2l[is-1-nx]; nb8->d = g2l[is-1-nx]; nb8->dr = g2l[is+1-nx];
  }
  //DBGB0("SetPnts11")
  for(int i = 0;i< part->nids;i++){
    for(int j=0;j<part->rsinfo[i].sn;j++){
      int is = part->rsinfo[i].sspnts[j];
      part->rsinfo[i].sspnts[j] = g2l[is];
    }
    for(int j=0;j<part->rsinfo[i].rn;j++){
      int is = part->rsinfo[i].srpnts[j];
      part->rsinfo[i].srpnts[j] = g2l[is];
    }
  }
  //DBGB0("SetPnts12")
  for(int i = 0;i< part->nids;i++){
    VFREE(part->rsinfo[i].mrpnts);
    VFREE(part->rsinfo[i].mspnts);
  }
  //DBGB0("SetPnts13")
  VFREE(aplist);
  VFREE(g2l);
  //DBGB0("SetPnts15")
}

/*-----------------------------------------------------------------------------------------------
 *SetPartMatrix
 *return :
 *   recto:
 *   ptype: pointer Type */
void  irrp_SetPartMatrix(Rect*recto,int* ptype){
  int *  pnts;
  if(gsi->parttype != 3){
    gsi->parttype = 3;
    Rect ro=cpart->recto;
    Rect ri=cpart->recti;
    cpart->nx  = ro.x1 - ro.x0 +1;
    cpart->ny  = ro.y1 - ro.y0 +1;
    cpart->nxy = cpart->nx * cpart->ny;
    cpart->np  = cpart->nxy;
    cpart->npc = (ri.x1-ro.x0+1)*(ri.y1-ro.y0+1);//fixme
    cpart->iblv=0;
    RENEWZN(cpart->calpnts,cpart->snp+1);
    RENEWZN(cpart->mpnttype,cpart->nxy);
    RENEWZN(cpart->ijbe,cpart->nxy*2);

    for(int n = 1;n<=cpart->snp;n++){
      int psij = png2psij(n);
      cpart->calpnts[n] = psij;
      int pt=1;
      if(n > cpart->snpc)pt = 3;   
      else if(n > cpart->nps)pt = 2;
      cpart->mpnttype[psij] = pt; 
    }
    int *ibe=cpart->ijbe,*jbe=ibe+cpart->nxy;
    int *mt=cpart->mpnttype;
    for(int pj = 0;pj< part->ny;pj++){
      int it=pj*cpart->nx,mf=1;
      ibe[it++]=-1;
      for(int pi = 0;pi< part->nx;pi++){
        int psij=pij2psij(pi,pj);
        if(mf){
          if(mt[psij]){ ibe[it++]=pi; mf=0; }
        }else{
          if(!mt[psij]){ ibe[it++]=pi; mf=1; }
        }
      }
      ibe[it++]=-1;
    }
    for(int pi = 0;pi< part->nx;pi++){
      int jt=pi*cpart->ny,mf=1;
      jbe[jt++]=-1;
      for(int pj = 0;pj< part->ny;pj++){
        int psij=pij2psij(pi,pj);
        if(mf){
          if(mt[psij]){ jbe[jt++]=pi; mf=0; }
        }else{
          if(!mt[psij]){ jbe[jt++]=pi; mf=1; }
        }
      }
      jbe[jt++]=-1;
    }
    for(int i=0;i<cpart->nids;i++){
      PI_Rsinfo* rsi=&cpart->rsinfo[i];
      VFREE(rsi->mrpnts);
      VFREE(rsi->mspnts);
      NEWZN(rsi->mrpnts,rsi->rn);
      NEWZN(rsi->mspnts,rsi->sn);
      for(int n = 0;n< rsi->rn;n++){
        rsi->mrpnts[n]=png2psij(rsi->srpnts[n]);
      }
      for(int n=0;n<rsi->sn;n++){
        rsi->mspnts[n]=png2psij(rsi->sspnts[n]);
      }
      rsi->rpnts = rsi->mrpnts;
      rsi->spnts = rsi->mspnts;
    }
  }
  if(recto)*recto =cpart->recto;
  if(ptype){
    NEWZN(ptype,cpart->nxy);
    CPYN(ptype,cpart->mpnttype,cpart->nxy);
  }
}

/*------------------------------------------------------------------------------------------------
 *SetPartSerial
 *input:
 * iblv_:The start of sequentializing: from 1 or 0.
 *return :
 *   gnpc:number of Calpnts in all part
 *   npc:number of Calpnts in cur part
 *    np:number of pnts(Include recv pnts) in cur part
 * plist:attribs of pnts[ np]
 * nb8  : neibhor of pnts[ npc] */
void irrp_SetPartSerial(int *gnpc,int * npc,int * np,PI_Pos**plist,PI_NB8**nb8,int *iblv_,int *nwps_){
  //DBGB0("irrp_SetPartSerial 1")
  if(gsi->parttype != 1){
    gsi->parttype = 1;
    cpart->np = cpart->snp ;
    cpart->npc = cpart->snpc;
    VFREE(cpart->calpnts);
    for(int i = 0;i< cpart->nids;i++){
      PI_Rsinfo* rsi = &cpart->rsinfo[i];
      VFREE(rsi->mrpnts);
      VFREE(rsi->mspnts);
      rsi->rpnts = rsi->srpnts;
      rsi->spnts = rsi->sspnts;
      //     mpnttype ! for matrix part
    }
  }
  //DBGB0("irrp_SetPartSerial 2")
  if(npc ) *npc  = cpart->snpc;
  if(np  ) *np   = cpart->snp;
  if(gnpc) *gnpc = gsi->gnpc;
  if(plist){
    VFREE(*plist);
    NEWZN(*plist,cpart->np+1);
    for(int n = 1;n<= cpart->np;n++){
      gsij2gij(cpart->ppos[n],&(*plist)[n].i,&(*plist)[n].j);
    }
  }
  if(nb8){
    VFREE(*nb8);
    NEWZN(*nb8,cpart->np+1);
    memcpy(*nb8 , cpart->nb8,sizeof(**nb8)*(cpart->np+1));
  }
  cpart->iblv=1;
  if(iblv_&&*iblv_){
    cpart->iblv=1;
  }else{
    cpart->iblv=0;
  }
  if(nwps_){
    int *  itmp;
    NEWZN(itmp,cpart->nxy);
    for(int i=0;i<cpart->nids;i++){
      int *p=cpart->rsinfo[i].sspnts;
      for(int n=0;n<cpart->rsinfo[i].sn;n++){
        itmp[p[n]]=1;
      }
    }
    int nwps=0;
    for(int i=0;i<cpart->nxy;i++){
      nwps+=itmp[i];
    }
    VFREE(itmp);
    *nwps_=nwps;
  }
  DBGB0("irrp_SetPartSerial End");
}
