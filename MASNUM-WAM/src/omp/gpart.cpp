#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include <fcntl.h>
#include <mdk_sdma.h>
#include<sched.h>
#include<string.h>
#include<pthread.h>
#include "hmalloc.h"
#include "mthread.h"
//#include <sys/mman.h>
#include <arm_sve.h>
#include <malloc.h>
#include <math.h>

#define NEWN(T,N) (T*)malloc(sizeof(T)*(N))
#define NEWZN(p,N) p=(typeof(p))malloc(sizeof(*p)*(N));if(p)memset(p,0,sizeof(*p)*(N))
#define NEW(T,n) (T*) malloc(sizeof(T)*n)
//#define NEW(p,n) *(void*)malloc(sizeof(*p)*n)
#define RENEW(p,n) *(void**)&p= realloc(p,sizeof(*p)*n)
//#define RENEW(T,n) i(void*)p= realloc(p,sizeof(*p)*n)

#define SPM (1<<5)
#define NPM (SPM-1)
#define RPM (SPM*2)
#define RPN 64
#include "gpart.h"

extern int mpi_id;
__BEGIN_DECLS
int GetVInt(int volatile *volatile p);
int ntdelay(int n,int vt);
__END_DECLS
void prtieind(int*ieind,int nwps,int nwpc,int nwpa,int LXN,int LYN);
void pmap(int id,int LXN,int LYN,int *pemask);
void pgind(int ind,int *ieind,int*pemask,int NX,int NY);
void precvseg(int ind,recvcomm *recv_segs);
static void zunlock(int *p){ *p=0;}
static void zlock(int *p){ while(__sync_val_compare_and_swap(p,0,1))usleep(1);}
static int wlock=0;
/*
 *本函数基于ieind的值将数据分为NPART个组，并将分组信息记录在partition
 *ieind的值分为四个区间，ieind==0：该点不是数据点无需处理
 *0<ieind<=nwps：该点为计算点，计算时依赖从其他进程接收的数据
 *nwps<ieind<=nwpc：该点为不需要接收外部数据的纯计算点
 *nwpc<ieind<nwpa：该点为该进程计算结束时需要向外传输的数据点
 *分组时首先对（nwps，nwpc]的纯计算点进行分组，分成NPART组后再向外扩展
 *搜索每组的边界，最终将每组的全部数据点，以及与该组相邻其他组的计算点
 *分别存储在partition的NPART个数组中，recto存储每组数据在全局中的坐标信息
 *ixb返回该组x轴坐标最小值，ixn返回该组x轴长度
 *jyb返回该组y轴坐标最小值，jyn返回该组y轴长度
 */
int *gidg_=NULL;
void partitionmap(Part*part, int* ieind,int LXN,int LYN,int nwps,int nwpc,int nwpa,int NCX,int NCY,int adjust) {
  int NPART=NCX*NCY;
  int XYN=LXN*LYN;                                      // XYN代表要处理的矩阵的规模
  int *gidr,*gidg;
  gidr = NEWN(int,(XYN+1)*2 );         //存储分组id，前XYN+1存储先按y分组的id(1-NCY)
  memset(gidr, 0, (XYN+1)*2 * sizeof(int));
  gidg = gidr+(XYN+1);                                  //存储最终分到的全局id（1-NPART）
  gidg=gidg_+1;
  int total=nwpc-nwps;                                  //代表矩阵中有效计算点的数量，作为分组的总数
#define VMAX(v,n) if(v<n)v=n
#define VMIN(v,n) if(v>n)v=n
  int groupsize;
  if(adjust>=0) groupsize= total/NPART+1 + adjust;        //共NPART组，第2-NPART组每组多分adjust个计算点
  else groupsize=total/(NPART+0.01*(adjust))+1;           //第1组 减少计算能力 adjust%
  int groupsize00 = total- groupsize * NPART +groupsize;//由于其他组多分，第1组分到的计算点数量
  int count = 0, i , j , idx , iy , ix;                 //count用于统计当前组实际已被分到的数量
  int thresold,groupid = 1;                             //thresold代表当前组预先分配计算点的总数，分组id从1开始;
  int xoff[NCY+1];                                      //按列分NCY组时，记录每个列组最右侧数据所在y轴的坐标
  int colsize=groupsize*(NCY-1)+groupsize00;            //由于最终第1组少分计算点，第一次按列分组时第1组计算点数量
  thresold =colsize;
  xoff[0] = 0;
  /*
   * 只有满足nwps<ieind<=nwpc的点才是有效的需要处理的数据点，如果不是数据点就跳过处理下一点
   * 将数据按列分成NCY组，从y轴方向遍历，如果该列遍历结束该组数据仍未分够，继续从下一列开始遍历
   * 如果该组数据已达到分配数量，再对下一组进行分配，直到全部数据点都已分配。
   */
  for (ix = 0; ix < LXN; ix++) {                //x轴方向遍历
    for (iy = 0; iy < LYN&&groupid<5; iy++) {   //y轴方向遍历，如果不判断groupid<5会有部分数据点分组id被重复更新
      idx = iy * LXN + ix;                      //数据点坐标
      if (ieind[idx] == 0||ieind[idx]>nwpc||ieind[idx]<=nwps) continue;
      gidr[idx] = groupid;                      //存储当前分组id
      count++;
      if (count >= thresold) {
        /*
         *第1组的分配数量与其他组都不同，2-NCY组都相同
         *当前组分够后，thresold更新为2-NCY组的分配数量；
         */
        thresold = groupsize*NCY;
        xoff[groupid] = ix;  //记录当前组最右侧数据所在行序号
        groupid++;           //更新为下一组序号
        count = 0;
      }
    }
  }
  /*
   *将数据按列分成NCY组后，再把每个列组的数据也分为NCX组，
   *四个列组内部按照先y轴坐标后x轴坐标的顺序遍历数据，
   *如果该组数据已分够，结束当前组，再对下一组进行分配，
   *直到该列组已分为NCX个组，再遍历下一个列组，最终将数据分为NCX*NCY=NPART个组。
   */
  thresold = groupsize00;  //thresold初始化为第1组的分配数量
  for (i = 1 ; i <= NCY; i++) { //遍历列组，i代表第i个列组
    Rect *r=&part[groupid].recti;
    groupid=1;count = 0;
    *r={LXN,LYN,0,0,0,0};
    for (iy = 0; iy < LYN; iy++) {   //y轴坐标从0遍历
      for (ix = xoff[i-1]; ix <= xoff[i]&&groupid<5; ix++) {  //x轴坐标从当前列组记录下的起始坐标遍历到下一列组的起始坐标
        idx = iy * LXN + ix;
        if (gidr[idx]> i) break;     //遍历到下一列组的数据，跳出循环
        /*
         *这种情况只出现在每个列组的第一列，由于该列可能有部分数据被分配给上一列组
         *所以该列可能有数据来自上一列组，此时满足gidr[idx]<i，只跳过当前数据
         */
        if (gidr[idx]!= i) continue;
        if (ieind[idx] <= 0||ieind[idx]>nwpc||ieind[idx]<=nwps) continue;
        gidg[idx]=((i-1) *NCX) + groupid;
        VMAX(r->ixb,ix);VMIN(r->ixe,ix);
        VMAX(r->jyb,iy);VMIN(r->jye,iy);
        count++;
        if (count >= thresold) {
          /*
           *第1组的分配数量与其他组都不同，2-NPART组都相同
           *当前组分够后，thresold更新为2-NPART组的分配数量；
           */
          thresold = groupsize;
          groupid++; count = 0;
          r=&part[groupid].recti;
          *r={LXN,LYN,0,0,0,0};
        }
      }
    }
  }
  gidg[-1]=0;                       //设置越界点id为0
  for (iy = 0; iy < LYN; iy++) {    //y轴坐标遍历
    for (ix = 0; ix < LXN; ix++) {  //x轴坐标遍历
      int ig;
      idx = iy * LXN + ix;
      /*
       *将数据分为NPART个组后，对gidg内尚未分组的其他非计算点进行分类，分类基于ieind数组
       *如果该点的ieind值为0，代表该点为不处理的边界，不分入其他组，不做处理，gidg仍为0
       *如果该点的ieind值大于nwpc，代表该点为该进程与其他进程通信的边界，gidg设置为RPN
       *其他情况代表该点ieind值大于mwps并且小于等于nwpc，为某组数据与外部进行通信的边界
       *此时搜索该点周围8个方向的分组id，如果找到第ig（1-NPART）组，将该点分组设置为SPM+ig。
       */
      int lu=idx+LXN-1,u=idx+LXN,ru=idx+LXN+1;      //左上      上      右上
      int l =idx-1,              r =idx+1;          // 左      原点      右
      int ld=idx-LXN-1,d=idx-LXN,rd=idx-LXN+1;      //左下      下      右下
      if(ix==0    )lu=l=ld=-1;                  //左寻找越界
      if(ix>=LXN-1)ru=r=rd=-1;                  //右寻找越界
      if(iy==0    )ld=d=rd=-1;                  //下寻找越界
      if(iy>=LYN-1)lu=u=ru=-1;                  //上寻找越界
      if (ieind[idx]>nwpc) gidg[idx]=RPM;      //该点为当前进程对外发送数据的通信边界，设置为RPM
      if (ieind[idx]>nwps) continue;           //该点为当前进程对外发送数据的通信边界，设置为RPN
      if (gidg[idx]||ieind[idx]==0) continue;   //该点已被分组或该点为无需处理的边界点
      if(iy==1) {
        ig=gidg[u ];
        if(ig==0)ig=gidg[u+LXN ];
        if(ig!=0&&ig<RPM) {gidg[idx]=(ig&NPM)+SPM;continue;}
      }
      if(iy==LYN-2) {
        ig=gidg[d ];
        if(ig==0)ig=gidg[d-LXN ];
        if(ig!=0&&ig<RPM) {gidg[idx]=(ig&NPM)+SPM;continue;}
      }
      if(ix==1) {
        ig=gidg[r ];
        if(ig==0)ig=gidg[r+1 ];
        if(ig!=0&&ig<RPM) {gidg[idx]=(ig&NPM)+SPM;continue;}
      }
      if(ix==LXN-2) {
        ig=gidg[l ];
        if(ig==0)ig=gidg[l-1 ];
        if(ig!=0&&ig<RPM) {gidg[idx]=(ig&NPM)+SPM;continue;}
      }else
        if     (u >=0&&(ig=gidg[u ]) && ig<SPM) gidg[idx]=ig+SPM;
        else if(d >=0&&(ig=gidg[d ]) && ig<SPM) gidg[idx]=ig+SPM;
        else if(l >=0&&(ig=gidg[l ]) && ig<SPM) gidg[idx]=ig+SPM;
        else if(r >=0&&(ig=gidg[r ]) && ig<SPM) gidg[idx]=ig+SPM;
        else if(lu>=0&&(ig=gidg[lu]) && ig<SPM) gidg[idx]=ig+SPM;
        else if(ld>=0&&(ig=gidg[ld]) && ig<SPM) gidg[idx]=ig+SPM;
        else if(ru>=0&&(ig=gidg[ru]) && ig<SPM) gidg[idx]=ig+SPM;
        else if(rd>=0&&(ig=gidg[rd]) && ig<SPM) gidg[idx]=ig+SPM;
    }
  }
  //pmap(0,LXN,LYN,gidg);
  /*
   *将分组数据输出到文件，该段代码仅在测试时使用
   */
  for (int p = 1; p <=NPART; p++) {
    part[p].recto ={ LXN,LYN,0,0,0,0};  // 初始化坐标最小值
  }
  //遍历所有点的分组id，更新每个组的ixa，ixn，iya，iyn
  for (int iy = 0; iy < LYN; iy++) {
    for (int ix = 0; ix < LXN; ix++) {
      int id = gidg[iy * LXN + ix];
      int ind=id&NPM;    //&NPM：筛选出该组通信边界，也作为该组的一部分
      if (id >= 1 && id < RPM) {
        Rect *r=&part[ind].recto;
        VMAX(r->ixe,ix);VMIN(r->ixb,ix);
        VMAX(r->jye,iy);VMIN(r->jyb,iy);
      }
    }
  }
  /*
   *遍历结束后每组数据往外扩展一圈，以将该组周围一圈可能存在的通信边界也记录
   *同时ixn，jyn的含义也更新，ixn，jyn分别存储该组在x轴方向，y轴方向的长度
   */
  for (int p = 1; p <=NPART; p++) {
    Rect *r=&part[p].recto;
    // 左下角扩展一列 右上角扩展一列
    // 左下角扩展一行 右上角扩展一行
    r->ixb--;r->ixe++;
    r->jyb--;r->jye++;
    r->ixn=r->ixe-r->ixb+1;  // 存储x轴长度
    r->jyn=r->jye-r->jyb+1;  // 存储y轴长度
    r=&part[p].recti;
    r->ixn=r->ixe-r->ixb+1;
    r->jyn=r->jye-r->jyb+1;
  }

  /*
   *根据recto存储的每组数据规模以及坐标信息，将gidg存储的数据转存到NPART个数组中
   *转存之后再判断每组中的与该组分组id不同的点是否与该组数据相邻，
   *如果相邻就保留该点的分组id，如果不相邻就设置为0
   */
  for (int ip = 1; ip <=NPART; ip++) {                 //1-NPART个组
    int ixb=part[ip].recto.ixb,jyb=part[ip].recto.jyb;  //ixb：起始点x轴坐标，jyb：起始点y轴坐标
    int ixn=part[ip].recto.ixn,jyn=part[ip].recto.jyn;  //ixn：x轴方向长度  ，jyn：y轴方向长度
    int*pemask=part[ip].pemask=NEWN(int,ixn*jyn+1);
    for (iy = 0; iy < jyn; iy++) {    //组内y轴坐标
      for (ix = 0; ix < ixn; ix++) {  //组内x轴坐标
        int idxg=(iy+jyb)*LXN+ixb+ix; //组内坐标转化为全局坐标
        int idx=iy*ixn+ix;
        int idg=gidg[idxg];
        pemask[idx]=idg;                //转存数据
      }
    }
    for (iy = 0; iy < jyn; iy++) {
      for (ix = 0; ix < ixn; ix++) {
        int idx=iy*ixn+ix;            //局部坐标
        int ido=pemask[idx];            //分组id
        int id=ido&NPM;              //&NPM：筛选出该组通信边界，也作为该组的一部分
        if(id==ip)continue;           //该点位于该组内部，不处理
        /*
         *ip表示当前组的id，在搜索前先判断是否可以往上下左右扩展
         *由于之前的处理中每组的通信边界，0<ieind<nwps的数据为该组id+SPM
         *八个方向的pemask数据&NPM可以将每个节点的通信边界也作为该组数据的一部分
         *这样判断时就能将这部分数据信息也保留，同时对于nwpc<ieind的数据点以及
         *id为其他组的数据点，如果与该组相邻就仍保留id，如果不相邻就都设为0
         *该点周围8个方向坐标
         */
        int lu=idx+ixn-1,u=idx+ixn,ru=idx+ixn+1;      //左上      上      右上
        int l =idx    -1,          r =idx    +1;      // 左      原点      右
        int ld=idx-ixn-1,d=idx-ixn,rd=idx-ixn+1;      //左下      下      右下
        if(ix==0    )lu=l=ld=-1;                  //左寻找越界
        if(ix>=ixn-1)ru=r=rd=-1;                  //右寻找越界
        if(iy==0    )ld=d=rd=-1;                  //下寻找越界
        if(iy>=jyn-1)lu=u=ru=-1;                  //上寻找越界
        id=0;
        if     (u >=0&&ip==((pemask[u ])&NPM) ) id=ido;   //左侧与该点相邻
        else if(d >=0&&ip==((pemask[d ])&NPM) ) id=ido;   //右侧与该点相邻
        else if(l >=0&&ip==((pemask[l ])&NPM) ) id=ido;   //上侧与该点相邻
        else if(r >=0&&ip==((pemask[r ])&NPM) ) id=ido;   //下侧与该点相邻
        else if(lu>=0&&ip==((pemask[lu])&NPM) ) id=ido;   //左上与该点相邻
        else if(ld>=0&&ip==((pemask[ld])&NPM) ) id=ido;   //左下与该点相邻
        else if(ru>=0&&ip==((pemask[ru])&NPM) ) id=ido;   //右上与该点相邻
        else if(rd>=0&&ip==((pemask[rd])&NPM) ) id=ido;   //右下与该点相邻
        pemask[idx]=id;
      }
    }
    //pmap(ip,ixn,jyn,part[ip].pemask);
  }
  free(gidr);
}
void serial_comm(int ind,Part *part, Glb2Grp *glb2grp,int *ieind);
void sort_sendcomm(sendcomm* sc,int m);
void find_recv_block(Part *part,Glb2Grp *glb2grp,int src_id);
void find_send_block(Part *part,Glb2Grp *glb2grp,int src_id);

void serial_comm(int ind,Part *part, Glb2Grp *glb2grp,int *ieind){
  int GRXN=part->GRXN;
  int GRYN=part->GRYN;
  int block_id=part->block_id;
  int *pemask=part->pemask;
  int*visited=NEWN(int,GRXN*GRYN);
  memset(visited,0,sizeof(int)*GRXN*GRYN);
  int *grp_ieind=part->grp_ieind;
  xypos *iepos=part->iepos;
  int *grp2glb=part->grp2glb;
  int ixb=part->glb_ixb;
  int jyb=part->glb_jyb;
  int LYN=part->LYN;
  int LXN=part->LXN;
  int rb_ix0=GRXN;
  int rb_ix1=-1;
  int rb_iy0=GRYN;
  int rb_iy1=-1;

  //计算边界
  for (int jy=0;jy<GRYN;jy++){
    for (int ix=0;ix<GRXN;ix++){
      int id=pemask[jy*GRXN+ix]&NPM;
      if (id==block_id ){
        if(rb_iy0>jy) rb_iy0=jy;
        if(rb_iy1<jy) rb_iy1=jy;
        if(rb_ix0>ix) rb_ix0=ix;
        if(rb_ix1<ix) rb_ix1=ix;
      }
    }
  }
  int gidx=0;
#define RECPOS_L(idx,ix,jy)    \
  visited[idx]=1;       \
  grp_ieind[idx]=++gidx;\
  iepos[gidx].grp_iy=jy;\
  iepos[gidx].grp_ix=ix;\
  iepos[gidx].glb_iy=jy+jyb;\
  iepos[gidx].glb_ix=ix+ixb;\
  grp2glb[gidx]=ieind[(jy+jyb)*LXN+ix+ixb];
#define RECPOS_G(idx,ix,jy)    \
  glb2grp[(jy+jyb)*LXN+ix+ixb].id=block_id;\
  glb2grp[(jy+jyb)*LXN+ix+ixb].gidx=gidx;\


  //编号block_id+SPM层
  for (int jy=rb_iy0;jy<=rb_iy1;jy++){
    for (int ix=rb_ix0;ix<=rb_ix1;ix++){
      int idx=jy*GRXN+ix;
      if (!visited[idx] && pemask[idx]==block_id+SPM){
        RECPOS_L(idx,ix,jy);
        RECPOS_G(idx,ix,jy);
      }
    }
  }
  part->nwpt=gidx;
  //编号neighbor层
  int right_ptr=rb_ix1;
  for (int jy=rb_iy0;jy<=rb_iy1;jy++){
    for (int ix=right_ptr;ix>=rb_ix0;ix--){
      int idx=jy*GRXN+ix;
      int id=pemask[idx];
      if (!visited[idx] && id==block_id){
        right_ptr=ix;
        RECPOS_L(idx,ix,jy);
        RECPOS_G(idx,ix,jy);
      } 
      if (right_ptr==rb_ix0) break;
    }
  }
  int down_ptr=rb_iy0;
  for (int ix=rb_ix1;ix>=rb_ix0+1;ix--){
    for (int jy=down_ptr;jy<=rb_iy1;jy++){
      int idx=jy*GRXN+ix;
      int id=pemask[idx];
      if (!visited[idx] && id==block_id){
        down_ptr=jy;
        RECPOS_L(idx,ix,jy);
        RECPOS_G(idx,ix,jy);
      }
      if (down_ptr==rb_iy1) break;
    }
  }
  part->nwps=gidx;
  //编号inside层
  for (int jy=rb_iy0;jy<=rb_iy1;jy++){
    for (int ix=rb_ix0;ix<rb_ix1;ix++){
      int idx=jy*GRXN+ix;
      if (!visited[idx] && pemask[idx]==block_id){
        RECPOS_L(idx,ix,jy);
        RECPOS_G(idx,ix,jy);
      }
    }
  }
  part->nwpc=gidx;
  //编号outside层
  right_ptr=GRXN-1;
  for (int jy=0;jy<GRYN;jy++){
    for (int ix=right_ptr;ix>=0;ix--){
      int idx=jy*GRXN+ix;
      if (!visited[idx] && pemask[idx]!=0 && pemask[idx]!=RPN){
        RECPOS_L(idx,ix,jy);
        right_ptr=ix;
      }
      if (right_ptr==0) break;
    }
  }
  down_ptr=rb_iy0;
  for (int ix=GRXN-1;ix>=0;ix--){
    for (int jy=down_ptr;jy<GRYN;jy++){
      int idx=jy*GRXN+ix;
      if (!visited[idx] && pemask[idx]!=0 && pemask[idx]!=RPN){
        down_ptr=jy;
        RECPOS_L(idx,ix,jy);
      }
      if (down_ptr==GRYN-1) break;
    }
  }
  part->nwpo=gidx;

  //编号Z层
  for (int jy=0;jy<GRYN;jy++){
    for (int ix=0;ix<GRXN;ix++){
      int idx=jy*GRXN+ix;
      if (!visited[idx] && pemask[idx]==RPN){
        RECPOS_L(idx,ix,jy);
        glb2grp[(jy+jyb)*LXN+ix+ixb].id=0;
        glb2grp[(jy+jyb)*LXN+ix+ixb].gidx=ieind[(jy+jyb)*LXN+ix+ixb];
      }
    }
  }
  part->nwpa=gidx;
  free(visited);
}

void sort_sendcomm(sendcomm* sc,int m){
  int ij=0;
  int pm=0;
  if(sc->src_id>1)pm=1;
  for (int j=0;j<=m;j++){
    sendsegs *cur_segs=&sc->ssegs[j];
    if (sc->ssegs[j].dst_id!=-1){
      if(sc->src_id==0){
        if (sc->ssegs[j].dst_id==1)pm=1;
      }else if(sc->src_id==1){
        if (sc->ssegs[j].dst_id==0) pm=1;
      }
      if(ij!=j)sc->ssegs[ij]=sc->ssegs[j];
      ij++;
      //printf("sort_sendcomm %2.2d %2.2d %2.2d %2.2d =>%2.2d %4d %2.2d\n",
      //       mpi_id,ti->igrp,ti->indg,sc->src_id,cur_segs->dst_id,cur_segs->nsegs,j);fflush(stdout);
    }
  }
  if(pm==0){
    sendsegs *cur_segs=&sc->ssegs[ij];
    cur_segs->nsegs=-1;
    if(sc->src_id==0){
      cur_segs->dst_id=1;
    } else if(sc->src_id==1){
      cur_segs->dst_id=0;
    }
    ij++;
  }
  for (;ij<=m;ij++){
    sc->ssegs[ij].dst_id=-1;
  }
}

void find_recv_block(Part *part,Glb2Grp *glb2grp,int src_id){
	Part *pt=part+src_id;
  int GRXN=pt->GRXN;
  int GRYN=pt->GRYN;
  int LYN=pt->LYN;
  int LXN=pt->LXN;
  int dst_id=pt->block_id;
  int start,end;
  if (dst_id==0){
    start=1;
    end=pt->nwps;
  }
  else{
    start=pt->nwpc+1;
    end=pt->nwpa;
  }
  for(int i=start;i<=end;i++){
    int glb_ix=pt->iepos[i].glb_ix;
    int glb_iy=pt->iepos[i].glb_iy;
    int grp_ix=pt->iepos[i].grp_ix;
    int grp_iy=pt->iepos[i].grp_iy;
    int src_id=glb2grp[glb_iy*LXN+glb_ix].id;
    if(src_id==dst_id) continue;
    if(src_id<0){
      printf("%d %d %d %d %d %d\n",dst_id,i,glb_ix,glb_ix,glb_iy*LXN+glb_ix,src_id);
      continue;
    }
    int dst_number=pt->grp_ieind[grp_iy*GRXN+grp_ix];
    int src_number=glb2grp[glb_iy*LXN+glb_ix].gidx;
    recvsegs *cur_mes=&pt->recv_segs.rsegs[src_id];
    cur_mes->src_id=src_id;
    int block_n=cur_mes->nsegs;
    int flag=false;
    for (int j=0;j<=block_n;j++){
      if (src_number==cur_mes->segs[j].sib-1 
        && dst_number==cur_mes->segs[j].dib-1){
        cur_mes->segs[j].sib=src_number;
        cur_mes->segs[j].dib=dst_number;
        cur_mes->segs[j].sin++;
        cur_mes->segs[j].din++;
        flag=true;
        break;
      }else if (src_number==(cur_mes->segs[j].sib+cur_mes->segs[j].sin) && 
        dst_number==(cur_mes->segs[j].dib+cur_mes->segs[j].din)){
        cur_mes->segs[j].sin++;
        cur_mes->segs[j].din++;
        flag=true;
        break;
      }else{
        continue;
      }
    }
    if (!flag){
      block_n++; cur_mes->nsegs++;
      if(block_n>=MSEGS){
        printf("Error:%d recv nsegs %d >MSEGG %d, at dst_id %d from src_id %d\n",mpi_id,block_n,MSEGS,dst_id,src_id);
        for(int i=0;i<block_n;i++){
          printf("%d %d %d %d %d\n",i, 
                 cur_mes->segs[i].sib, cur_mes->segs[i].sin,
                 cur_mes->segs[i].dib, cur_mes->segs[i].din);
        }
        //fixme at 0.25*0.25 4nodes
        exit(0);
      }else{
				cur_mes->segs[block_n].sib=src_number;
				cur_mes->segs[block_n].sin=1;
				cur_mes->segs[block_n].dib=dst_number;
				cur_mes->segs[block_n].din=1;
      }
      //printf("find_recv_block %2.2d %2.2d => %2.2d %2.2d \n",
      //       mpi_id,src_id,dst_id,cur_mes->nsegs);fflush(stdout);
      
    }
  }
}

void find_send_block(Part *part,Glb2Grp *glb2grp,int src_id){
  int NPART=part->NPART;
	sendcomm *sc=&part[src_id].send_segs;
	sc->src_id=src_id;
  recvcomm *rc =&part[src_id].recv_segs;
  for (int j=0;j<=NPART;j++){
		recvsegs *cur_recv=&rc->rsegs[j];
		if(cur_recv->nsegs<0)continue;//not continued
		sc=&part[cur_recv->src_id].send_segs;
    if(sc->nlj<0)sc->nlj=0;
    sendsegs *cur_send=sc->ssegs+sc->nlj++;
    cur_send->dst_id=src_id;
    cur_send->nsegs=cur_recv->nsegs+1;
		for (int i=0;i<=cur_recv->nsegs;i++){
			cur_send->segs[i]=cur_recv->segs[i];
		}
	}
}
void prtieind(int*ieind,int nwps,int nwpc,int nwpa,int LXN,int LYN);
void pg2g(Part *part,Glb2Grp *gg);
static Part *dpart;
void GPart(Part *part,int *ieind,int NCX,int NCY,int nwps,int nwpc,int nwpa,int LXN,int LYN,int adjust){
  int NPART=NCX*NCY;
  NEWZN(gidg_,LXN*LYN+1);
  dpart=part;
  memset(part,0,sizeof(part[0])*(NPART+1));
  partitionmap(part,ieind,LXN,LYN,nwps,nwpc,nwpa,NCX,NCY,adjust);
  Glb2Grp *glb2grp;
  NEWZN(glb2grp,LXN*LYN);
  //memset(glb2grp,-1,LXN*LYN*sizeof(glb2grp[0]));
  {
    Part *pt=part+0;
    memset(&pt->recv_segs,-1,sizeof(recvcomm));
    memset(&pt->send_segs,-1,sizeof(sendcomm));
    pt->block_id=0;
    pt->LXN=LXN;
    pt->LYN=LYN;
    pt->GRXN=LXN;
    pt->GRYN=LYN;
    pt->NPART=NPART;

    pt->grp2glb=(int*)glb2grp;
    pt->grp_ieind=NEWN(int,LXN*LYN);
    memcpy(pt->grp_ieind,ieind,LXN*LYN*sizeof(int));
    pt->iepos=NEWN(xypos,LXN*LYN+1);
    for (int jy=0;jy<LYN;jy++){
      for (int ix=0;ix<LXN;ix++){
        int idx=jy*LXN+ix;
        int number=pt->grp_ieind[idx];
        pt->iepos[number].grp_ix=ix;
        pt->iepos[number].grp_iy=jy;
        pt->iepos[number].glb_ix=ix;
        pt->iepos[number].glb_iy=jy;
      }
    }
    pt->nwpt=0;
    pt->nwps=nwps;
    pt->nwpc=nwpc;
    pt->nwpo=nwpa;
    pt->nwpa=nwpa;
    pt->send_segs.src_id=0;
    pt->recv_segs.dst_id=0;
#ifdef USESDMA
    pt->send_segs.sdma.n=0;
    pt->send_segs.sdma.chn_handle=NULL;
    pt->send_segs.sdma.process_id=0;
    pt->send_segs.sdma.chn_id=1;
#endif
  }

  //编号排序
  for (int i=1;i<=NPART;i++){
    memset(&part[i].recv_segs,-1,sizeof(recvcomm));
    memset(&part[i].send_segs,-1,sizeof(sendcomm));
    part[i].send_segs.nlj=0;
    part[i].LXN=LXN;
    part[i].LYN=LYN;
    part[i].GRXN=part[i].recto.ixn;
    part[i].GRYN=part[i].recto.jyn;
    part[i].block_id=i;
    part[i].glb_ixb=part[i].recto.ixb;
    part[i].glb_jyb=part[i].recto.jyb;
    NEWZN(part[i].grp_ieind,part[i].recto.ixn*part[i].recto.jyn);
    NEWZN(part[i].grp2glb,(part[i].recto.ixn*part[i].recto.jyn+1));
    NEWZN(part[i].iepos,(part[i].recto.ixn*part[i].recto.jyn+1));
    part[i].recv_segs.dst_id=i;
    part[i].send_segs.src_id=i;
    part[i].NPART=NPART;
#ifdef USESDMA
    part[i].send_segs.sdma.n=0;
    part[i].send_segs.sdma.chn_handle=NULL;
    part[i].send_segs.sdma.process_id=0;
    part[i].send_segs.sdma.chn_id=i+1;
#endif
    serial_comm(i,&part[i],glb2grp,ieind);
  }
#if 0
  for(int i=0;i< part->LYN* part->LXN;i++){
    if( glb2grp[i].id<0&&ieind[i]>0){
      printf(" AA %d %d %d\n",i,ieind[i],gidg_[i+1]);
    }
  }
  prtieind(ieind,nwps,nwpc,nwpa,LXN,LYN);
  pg2g(part,glb2grp); 
#endif
  //构建接收段落体
  for (int i=0;i<=NPART;i++){
    find_recv_block(part,glb2grp,i);
  }
  //构建发送段落体
  for (int i=0;i<=NPART;i++){
    find_send_block(part,glb2grp,i);
  }
  for (int i=0;i<=NPART;i++){
    sort_sendcomm(&part[i].send_segs,NPART);
  }
}
#define VFREE(p) if(p) free(p)
void gp_freepart(Part *part){
  VFREE(part->pemask);
  VFREE(part->iepos);
  VFREE(part->grp2glb);
  VFREE(part->grp_ieind);
}
//==================================================
void gp_filldst(sendcomm **scs,int m){
  for (int i=0;i<m;i++){
    sendcomm* sc=scs[i];
    for (int j=0;j<m;j++){
      sendsegs *cur_segs=&sc->ssegs[j];
      int did=cur_segs->dst_id;
      if (did==-1)continue;//break;
      //printf("gp_filldst %2.2d %2.2d %2.2d %2.2d =>%2.2d %4d %2.2d\n",
      //       mpi_id,ti->igrp,ti->indg,sc->src_id,cur_segs->dst_id,cur_segs->nsegs,j);fflush(stdout);
      if (cur_segs->nsegs<=0)continue;//break;
      cur_segs->dst_ee=scs[cur_segs->dst_id]->src_ee;
      sendcomm *ss=scs[did];
      for (int k=0;k<m;k++){
        sendsegs *dst_segs=&ss->ssegs[k];
        if (dst_segs->dst_id==-1)continue;//break;
        if (dst_segs->dst_id==sc->src_id){
          cur_segs->pflag=dst_segs->rflag;
          break;
        }
      }
			//printf("WW 1240 %2d %2.2d %2d:%2d %2d:%2d %p\n",mpi_id,i,j ,did,sc->src_id,cur_segs->nsegs,cur_segs->pflag);
    }
  }
}
void gp_send_data(sendcomm* sc,int itemsize,int m,int RFS){
  void*src_copy,*dst_copy;
  int src_id=sc->src_id;
  src_copy=sc->src_ee;
  sc->sflag=RFS;
  //zlock(&wlock);
  //printf("senddata %2.2d %2.2d %2.2d %2.2d %2.2d %d S\n",mpi_id,ti->igrp,ti->indg,sc->src_id,m,RFS);fflush(stdout);
  for (int j=0;j<=m;j++){
    sendsegs *cur_segs=&sc->ssegs[j];
    //printf("senddata %2.2d %2.2d %2.2d %2.2d =>%2.2d %4d %2.2d %2.2d\n",
    //       mpi_id,ti->igrp,ti->indg,src_id,cur_segs->dst_id,cur_segs->nsegs,j,RFS);fflush(stdout);
    if (cur_segs->dst_id==-1)break;
    dst_copy=cur_segs->dst_ee;
    int *dst_flag=cur_segs->pflag;
    int nsegs=cur_segs->nsegs;
    if(nsegs>=0){
      struct seg *segs=cur_segs->segs;
      for (int k=0;k<nsegs;k++){
        //fixme sin sib dib should not little then 0
        if(segs[k].sin<=0||segs[k].sib<0||segs[k].dib<0){
          printf("CPYEE %d %2.2d => %2.2d %d %d %d %d %d \n",
                 mpi_id,src_id,cur_segs->dst_id,k,nsegs,segs[k].dib,segs[k].sib,segs[k].sin); fflush(stdout);
          printf("set FLG %d %2.2d %2.2d => %2.2d %d %p %p\n",mpi_id,ti->igrp,sc->src_id,cur_segs->dst_id,RFS,dst_flag,&sc->sflag);fflush(stdout);
          continue;
        }
        int seg_n    =itemsize*segs[k].sin;
        int src_start=itemsize*segs[k].sib;
        int dst_start=itemsize*segs[k].dib;
        //printf("CPY %d %2.2d => %2.2d %p %p %d %d %d %d %d %d\n",mpi_id,src_id,cur_segs->dst_id,(char*)dst_copy+dst_start,(char*)src_copy+src_start,dst_start,src_start,seg_n,segs[k].dib,segs[k].sib,segs[k].sin); fflush(stdout);
        //usleep(100);
        memcpy((char*)dst_copy+dst_start,(char*)src_copy+src_start,seg_n);
        //printf("CPY %d %2.2d => %2.2d %p %p %d %d %d End\n",mpi_id,src_id,cur_segs->dst_id,(char*)dst_copy+dst_start,(char*)src_copy+src_start,dst_start,src_start,seg_n); fflush(stdout);
      }
      //printf("set FLG %d %2.2d %2.2d => %2.2d %d\n"      ,mpi_id,ti->igrp,sc->src_id,cur_segs->dst_id,RFS);fflush(stdout);
      //printf("set FLG %d %2.2d %2.2d => %2.2d %d %p %p\n",mpi_id,ti->igrp,sc->src_id,cur_segs->dst_id,RFS,dst_flag,&sc->sflag);fflush(stdout);
			//usleep(5000);
      memcpy(dst_flag,&sc->sflag,16*sizeof(int));
      //printf("set FLG %d %2.2d %2.2d => %2.2d %d end\n"  ,mpi_id,ti->igrp,sc->src_id,cur_segs->dst_id,RFS);fflush(stdout);
    }
  }
  //zunlock(&wlock);usleep(1000);
}

#define MAXCHECK 5000000
extern "C" int Wait_LGE( volatile int *pv,int cv,int nt ,int line);
static inline uint64_t arm64_cntvct(void) ;
void gp_check_recv(sendcomm* sc,int m,int RFS){
  for (int j=0;j<=m;j++){
    sendsegs *cur_segs=&sc->ssegs[j];
    if (cur_segs->dst_id==-1)break;
    if (cur_segs->nsegs<=0)continue;
    if(sc->src_id==cur_segs->dst_id)continue; //Fixme 
    //zlock(&wlock);
    //printf("check FLG %d %d < %d %d\n",mpi_id,sc->src_id,cur_segs->dst_id,RFS);fflush(stdout);
    //zunlock(&wlock);usleep(1000);
    Wait_LGE(cur_segs->rflag,RFS,1,__LINE__);
  }
}

#ifdef USESDMA
//  sdma  =========
bool initsdma(zsdma *psdma){
  int ret=-1;
  char sdma_equip_name[128];
  psdma->eqid=sdma_nearest_id();
  sprintf(sdma_equip_name, "/dev/sdma%d", psdma->eqid);
  psdma->sdmafd = open(sdma_equip_name, O_RDWR);
  ret = sdma_get_process_id(psdma->sdmafd, &psdma->process_id);
  psdma->chn_handle = sdma_init_chn(psdma->sdmafd,psdma->chn_id); //alloc申请独占channel，对应free，init申请共享channel，对应deinit
  psdma->sqe_task=NULL; psdma->n= psdma->m=0;
  return 0;
}
void add_sdma_task(zsdma*aa, void*src, void*dest, int len){
  if (aa->n>=aa->m){
    int m=aa->m;
    aa->m+=10;
    RENEW(aa->sqe_task,aa->m);
    memset(aa->sqe_task+m,0,sizeof(sdma_sqe_task_t)*(aa->m-m));
  }
  sdma_sqe_task_t* sqe_task=aa->sqe_task+aa->n;
  sqe_task->src_addr = (unsigned long)src;
  sqe_task->dst_addr = (unsigned long)dest;
  sqe_task->src_process_id = aa->process_id;
  sqe_task->dst_process_id = aa->process_id;
  sqe_task->src_stride_len = 0;
  sqe_task->dst_stride_len = 0;
  sqe_task->stride_num = 0;
  sqe_task->length = len;
  sqe_task->opcode = 0x0;//Z???
  sqe_task->next_sqe=NULL;
  if(aa->n>0)aa->sqe_task[aa->n-1].next_sqe = sqe_task;//Z???
  aa->n++;
}
void setSdmaTask(sendcomm* sc,int itemsize,int m){
  zsdma*ps=&sc->sdma;
  void *src_copy,*dst_copy;
  uint64_t cookie_src;
  uint64_t cookie_dst;
  int src_id=sc->src_id;
  if(ps->chn_handle==NULL)initsdma(ps);
  src_copy=sc->src_ee;

  for (int j=0;j<=m;j++){
    sendsegs *cur_segs=&sc->ssegs[j];
    if (cur_segs->dst_id==-1)break;
    dst_copy=cur_segs->dst_ee;
    int nsegs=cur_segs->nsegs;
    if(nsegs>0){
      struct seg *segs=cur_segs->segs;
      for (int k=0;k<nsegs;k++){
        //fixme sin sib dib should not little then 0
        if(segs[k].sin<=0||segs[k].sib<0||segs[k].dib<0)continue;
        uint32_t seg_n    =itemsize*segs[k].sin;
        uint32_t src_start=itemsize*segs[k].sib;
        uint32_t dst_start=itemsize*segs[k].dib;
        sdma_pin_umem(ps->sdmafd,(void*)((char*)src_copy+src_start),(uint32_t)seg_n,&cookie_src);
        sdma_pin_umem(ps->sdmafd,(void*)((char*)dst_copy+dst_start),(uint32_t)seg_n,&cookie_dst);
        add_sdma_task(ps, (char*)src_copy+src_start, (char*)dst_copy+dst_start, seg_n);
      }
      sdma_pin_umem(ps->sdmafd,(char*)cur_segs->pflag,64,&cookie_src);
      sdma_pin_umem(ps->sdmafd,(char*)&sc->sflag,64,&cookie_dst);  
      add_sdma_task(ps, (char*)&sc->sflag, (char*)cur_segs->pflag, 16*sizeof(int));
    }
  }
}
void gp_sdma_send_data(sendcomm* sc,int itemsize,int m,int RFS){
  zsdma*ps=&sc->sdma;
  if(ps->n==0)setSdmaTask(sc,itemsize,m);
  sc->sflag=RFS;
  sdma_icopy_data(ps->chn_handle, ps->sqe_task, ps->n,&ps->request);
  sdma_iwait_chn(ps->chn_handle,&ps->request);
}
#endif

int gp_CheckSend(Part *part){
  //构建测试数据
  int NPART=part[0].NPART;
  Part *p0=part+0;
  int *se_data[NPART+1];
  int *gdata=se_data[0]=NEWN(int,(p0->nwpa+1));
  sendcomm*scs[NPART+1];
  for (int i=0;i<=NPART;i++){
    scs[i]=&part[i].send_segs;
  }
  scs[0]->src_ee=gdata;
  int LXN=p0->LXN,LYN=p0->LYN;
  int *glb_iepos;
  NEWZN(glb_iepos,(p0->nwpa+1));
  for (int ij=0;ij<LXN*LYN;ij++){
    glb_iepos[p0->grp_ieind[ij]]=ij;// ix=i%LXN;iy=i/LXN;
  }//nwpc+1...nwpa
  for(int iac=p0->nwpc+1;iac<=p0->nwpa;iac++){
    gdata[iac]=glb_iepos[iac];
  }
  for (int i=1;i<=NPART;i++){
    Part *pt=part+i;
    NEWZN(se_data[i],part[i].nwpa+1);
    int *p=se_data[i];
    for(int j=1;j<=pt->nwps;j++){
      p[j]=pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN;
    }
    for(int j=pt->nwps+1;j<=pt->nwpc;j++){
      p[j]=-(pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN);
    }
    part[i].send_segs.src_ee=p;

  }
  gp_filldst(scs,NPART);
  //开始发送
  int RFS=123;
  for (int i=0;i<=NPART;i++){
#ifdef USESDMA
    gp_sdma_send_data(&part[i].send_segs,sizeof(int),NPART,RFS);
#else
    gp_send_data(&part[i].send_segs,sizeof(int),NPART,RFS);
#endif
  }
  int err=0;
  for (int i=1;i<=NPART;i++){
    Part *pt=part+i;
    int *p=se_data[i];
    for(int j=1;j<=part[i].nwps;j++){
      if(p[j]!=pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN){
        printf("S: %d %d %d %d\n",i,j,p[j],pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN);
        err++;
      }
    }
    for(int j=pt->nwps+1;j<=pt->nwpc;j++){
      if(p[j]!=-(pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN)){
        printf("C: %d %d %d %d\n",i,j,p[j],pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN);
        err++;
      }
    }
    for(int j=pt->nwpc+1;j<=pt->nwpo;j++){
      if(p[j]!=pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN){
        printf("O: %d %d %d %d\n",i,j,p[j],pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN);
        err++;
      }
    }
    for(int j=pt->nwpo+1;j<=pt->nwpa;j++){
      if(p[j]!=pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN){
        printf("A: %d %d %d %d\n",i,j,p[j],pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN);
        err++;
      }
    }
  }
  int *p=gdata;
  for(int j=1;j<=part[0].nwps;j++){
    if(p[j]!=part[0].iepos[j].glb_ix+part[0].iepos[j].glb_iy*LXN){
      printf("S0: %d %d %d %d\n",0,j,p[j],part[0].iepos[j].glb_ix+part[0].iepos[j].glb_iy*LXN);
    }
  }
  free(glb_iepos);
  for (int i=0;i<=NPART;i++){
    free(se_data[i]);
  }
  return err;
}
static inline uint64_t arm64_cntvct(void) {
  uint64_t tsc;
  asm volatile("mrs %0, cntvct_el0" : "=r" (tsc));    //读取系统时间戳
  return tsc;
}
#ifdef GPART_TEST
int mpi_id=0;
#ifdef UTHREAD
int bindcpu(int id){
  cpu_set_t mask;  //CPU核的集合
  CPU_ZERO(&mask);    //置空
  CPU_SET(id,&mask);   //设置亲和力
  if (sched_setaffinity(0, sizeof(mask), &mask) == -1){//设置线程CPU亲和力
    printf("warning: could not set CPU affinity, continuing...\n");
    return -1;
  }
  return 0;
}
typedef struct Para{
  int cpu_id;
  sendcomm* sc;
  int itemsize;
  int NPART;
  int RFS;
  Part *part;
}Para;

void* threadFunction(void* arg){  
  Para *cur_para=(Para *)arg;
  int cpu_id=cur_para->cpu_id;
  int RFS=cur_para->RFS;
  bindcpu(cpu_id);
  sendcomm *sc=cur_para->sc;
  int NPART=cur_para->NPART;
  int itemsize=cur_para->itemsize;
  uint64_t send_start,send_end;
  send_start=arm64_cntvct();
#ifdef USESDMA
  gp_sdma_send_data(cur_para->sc,itemsize,NPART,RFS);
#else
  gp_send_data(cur_para->sc,itemsize,NPART,RFS);
#endif
  send_end=arm64_cntvct();
  printf("send time: %3.3d %2.2d %7.3fms\n",cpu_id,cpu_id/38,(send_end-send_start)/100000.);
  gp_check_recv(sc,NPART,RFS);
  return NULL;
}
int CheckSendThreads(Part *part){
  //构建测试数据
  int nitem=1;
  int NPART=part[0].NPART;
  Part *p0=part+0;
  int *se_data[NPART+1];
  int *gdata=se_data[0]=NEWN(int,(p0->nwpa+1)*nitem);
  sendcomm*scs[NPART+1];
  for (int i=0;i<=NPART;i++){
    scs[i]=&part[i].send_segs;
  }
  scs[0]->src_ee=gdata;
  for (int i=1;i<=NPART;i++){
    NEWZN(se_data[i],(part[i].nwpa+1)*nitem);
    scs[i]->src_ee=se_data[i];
  }
  gp_filldst(scs,NPART);
  int LXN=p0->LXN,LYN=p0->LYN;
  int *glb_iepos;
  NEWZN(glb_iepos,(p0->nwpa+1)*nitem);
  for (int ij=0;ij<LXN*LYN;ij++){
    glb_iepos[p0->grp_ieind[ij]]=ij;// ix=i%LXN;iy=i/LXN;
  }//nwpc+1...nwpa
  for(int iac=p0->nwpc+1;iac<=p0->nwpa;iac++){
    gdata[iac]=glb_iepos[iac];
  }
  for (int i=1;i<=NPART;i++){
    Part *pt=part+i;
    //    NEWZN(se_data[i],part[i].nwpa+1);
    int *p=se_data[i];
    for(int j=1;j<=pt->nwps;j++){
      p[j]=pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN;
    }
    for(int j=pt->nwps+1;j<=pt->nwpc;j++){
      p[j]=-(pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN);
    }
    part[i].send_segs.src_ee=p;
  }
  //开始发送
  pthread_t thread[NPART+1];
  Para mpy_param[NPART+1];
  for (int i=0;i<=NPART;i++){
    int id=(i-1)*38;
    if(i==0)id=0;
    else if(i==1)id=1;
    mpy_param[i].cpu_id=id;
    mpy_param[i].sc=&part[i].send_segs;
    mpy_param[i].itemsize=sizeof(int)*nitem;
    mpy_param[i].NPART=NPART;
    mpy_param[i].part=part;
    mpy_param[i].RFS=NPART;
  }
  for (int i=1;i<=NPART;i++){
    if (pthread_create(&thread[i],NULL,threadFunction,&mpy_param[i])!=0){
      printf("Failed to bind cpu%d\n",i);
    }
  }
  threadFunction(&mpy_param[0]);
  for (int i=1;i<=NPART;i++){
    pthread_join(thread[i],NULL);
  }

  int err=0;
  for (int i=1;i<=NPART;i++){
    Part *pt=part+i;
    int *p=se_data[i];
    for(int j=1;j<=part[i].nwps;j++){
      if(p[j]!=pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN){
        printf("S: %d %d %d %d\n",i,j,p[j],pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN);
        err++;
      }
    }
    for(int j=pt->nwps+1;j<=pt->nwpc;j++){
      if(p[j]!=-(pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN)){
        printf("C: %d %d %d %d\n",i,j,p[j],pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN);
        err++;
      }
    }
    for(int j=pt->nwpc+1;j<=pt->nwpo;j++){
      if(p[j]!=pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN){
        printf("O: %d %d %d %d\n",i,j,p[j],pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN);
        err++;
      }
    }
    for(int j=pt->nwpo+1;j<=pt->nwpa;j++){
      if(p[j]!=pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN){
        printf("A: %d %d %d %d\n",i,j,p[j],pt->iepos[j].glb_ix+pt->iepos[j].glb_iy*LXN);
        err++;
      }
    }
  }
  int *p=gdata;
  for(int j=1;j<=part[0].nwps;j++){
    if(p[j]!=part[0].iepos[j].glb_ix+part[0].iepos[j].glb_iy*LXN){
      printf("S0: %d %d %d %d\n",0,j,p[j],part[0].iepos[j].glb_ix+part[0].iepos[j].glb_iy*LXN);
    }
  }
  free(glb_iepos);
  for (int i=0;i<=NPART;i++){
    free(se_data[i]);
  }
  return 0;
}
#endif
void read_file(int **ieind,int *nwps,int *nwpc,int *nwpa,int *LXN,int *LYN){
  char fn[256];
  snprintf(fn, sizeof(fn), "ieind_in.txt");
  FILE *fo= fopen(fn,"r");
  fscanf(fo,"%d,%d,%d,%d,%d\n",nwps,nwpc,nwpa,LXN,LYN);
  int xn=*LXN,yn=*LYN;
  int *p=*ieind=NEWN(int,xn*yn);
  for(int iy=yn-1;iy>=0;iy--){
    for(int ix=0;ix<xn;ix++) fscanf(fo,"%8d", &p[iy*xn+ix]);
  }
  fclose(fo);
}
__attribute__ ((optnone)) int ntdelay(int n,int vt){
  for(int i=0;i<n;i++) vt=(vt*1357+2581);
  return vt;
}
__attribute__ ((optnone))int GetVInt(int volatile *volatile p){
  return *(int*)p;
}
Part*spart=NULL;
int main(){
  int nwps,nwpc,nwpa,*ieind;
  int LXN, LYN;
  Part part[MPART+1]={0};
  spart=part;
  read_file(&ieind,&nwps,&nwpc,&nwpa,&LXN,&LYN);
  GPart(part,ieind,4,4,nwps,nwpc,nwpa,LXN,LYN,2);
#ifdef UTHREAD
  int err=CheckSendThreads(part);
#else
  int err=gp_CheckSend(part);
#endif
  if (err!=0) printf("error!");
  for (int i=0;i<=MPART;i++){
    gp_freepart(&part[i]);
  }
  return 0;
}
#endif
// ================= output =============
/*
 *iepos[ 0:nwpa ]
 *ieind[0,LXN*LYN-1],ind=ix*LXN
 *nsp [ 0:nwpa ],flags,0:land,1:water,2:open boundary
 * */
void pmap(int id,int LXN,int LYN,int *pemask){
  char filename[128];
  snprintf(filename, sizeof(filename), "groupmap%d.dat",id);
  FILE* file = fopen(filename, "w");
  for (int row = LYN-1; row>=0; row--) {
    for(int col = 0; col < LXN; col++){
      int id=pemask[row*LXN+col];
      char c;
      if(id==0)c='_';
      else if(id<16)c='a'+id-1;
      else if(id<32)c='A'+id-1;
      else c='*';
      fprintf(file,"%c",c);
    }
    fprintf(file, "\n");
  }
  fclose(file);
}
void prtieind(int*ieind,int nwps,int nwpc,int nwpa,int LXN,int LYN){
    char filename[1024];
    snprintf(filename, sizeof(filename), "ieind_%2.2d.txt", mpi_id);
    FILE* file = fopen(filename, "w");
    fprintf(file,"%d,%d,%d,%d,%d\n",nwps,nwpc,nwpa,LXN,LYN);
    char cm[32]="0+1*";
    for (int iy = LYN-1; iy>=0; iy--) {
      for(int ix = 0; ix < LXN; ix++){
        int ia=ieind[iy*LXN+ix];
        fprintf(file,"%6d",ia);
      }
      fprintf(file,"\n");
    }
    fclose(file);
#if 0
    snprintf(filename, sizeof(filename), "iepos_%2.2d.txt", gppar.mpi_id);
    file = fopen(filename, "w");
    for (int iy = LYN-1; iy>=0; iy--) {
      for(int ix = 0; ix < LXN; ix++){
        int x=iepos[iy*LXN+ix].ix;
        int y=iepos[iy*LXN+ix].iy;
        fprintf(file,"%d,%d",x,y);
      }
      fprintf(file,"\n");
    }
    fclose(file);
#endif
}
void pg2g(Part *part,Glb2Grp *gg){
  int LYN=part->LYN;
  int LXN=part->LXN;
  char filename[1024];
  snprintf(filename, sizeof(filename), "g2g_%2.2d.txt", mpi_id);
  FILE* fo = fopen(filename, "w");
  fprintf(fo,"%d %d\n",LXN,LYN);
  for (int jy = LYN-1; jy>=0; jy--) {
    for (int ix=0;ix<LXN;ix++){
      int idx=jy*LXN+ix;
      fprintf(fo,"%5d:%d ",gg[idx].gidx,gg[idx].id);
    }
    fprintf(fo,"\n");
  }
  fprintf(fo,"\n");
}
#if 0
void pgind(int ind,int *ieind,int*pemask,int NX,int NY){
  printf("PART %2.2d %4d %4d\n",ind,NX,NY);
  for (int jy=0;jy<NY;jy++){
    for (int ix=0;ix<NX;ix++){
      printf("%5d(%5d) ",grp_ieind[jy*NX+ix],part[i].pemask[jy*NX+ix]);
    }
    printf("\n");
  }
  printf("\n");
}

void precvseg(int ind,recvcomm *recv_segs){
  for (int j=0;j<17;j++){
    if (rsegs->rsegs[j].src_id!=-1){
      printf("dst id:%d src id:%d nsegs:%d",rsegs->dst_id,rsegs->rsegs[j].src_id,rsegs->rsegs[j].nsegs);
      printf("\n");
      for (int k=0;k<10;k++){
        printf("sib[%d]:%d sin[%d]:%d ",k,rsegs->rsegs[j].segs[k].sib,k,rsegs->rsegs[j].segs[k].sin);
      }
      printf("\n");
      for (int k=0;k<10;k++){
        printf("dib[%d]:%d din[%d]:%d ",k,rsegs->rsegs[j].segs[k].dib,k,rsegs->rsegs[j].segs[k].din);
      }
      printf("\n\n");
    }
  }
}
#endif
