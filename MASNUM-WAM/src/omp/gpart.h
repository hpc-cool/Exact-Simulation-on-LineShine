#include <mdk_sdma.h>
#define MPART 16
#define MSEGS 20
typedef struct xypos{
  int grp_ix, grp_iy;
  int glb_ix,glb_iy;
  int glb_iac;
}xypos;

typedef struct Glb2Grp{
  int id, gidx;
}Glb2Grp;

struct seg{
  int sib,sin,dib,din;//16
};
struct recvsegs{
  int src_id;
  int nsegs;
  struct seg segs[MSEGS];
};
typedef struct recvcomm{
  int dst_id;
  struct recvsegs rsegs[MPART+1];
}recvcomm;
/* sendcomm sendsegs
 *  void*  src_ee
 *  sendsegs ssegs
 *    dst_id
 *    nsegs
 *    dst_ee
 *    segs
 *      sin,sib,dib
 * */
#ifdef USESDMA
typedef struct zsdma{
  int sdmafd,eqid;  //
  int chn_id;
  void *chn_handle;//8+8
  uint32_t process_id ;//16+4
  int n,m; //20+8
  sdma_request_t request;// 28+12 =40
  sdma_sqe_task_t *sqe_task;// 40+8=48
}zsdma;
#endif
struct sendsegs{
  int dst_id;
  int nsegs;
  void *dst_ee;//8+8 16
  int *pflag;//16+8
  long fill;//24 +8
  struct seg segs[MSEGS];//32+16*10
  int rflag[16];//16*12=64*3 +64 =64*4
};
typedef struct sendcomm{
  int sflag,src_id,nlj ;//8
  void*src_ee;//8+8
#ifdef USESDMA
  zsdma sdma;
#endif
  struct sendsegs ssegs[MPART+1];
}sendcomm;
typedef struct Rect{
  int ixb,jyb,ixe,jye,ixn,jyn;
}Rect;
typedef struct Part{
  int block_id;
  int LXN;
  int LYN;
  int NPART;
  int GRXN;
  int GRYN;
  int glb_ixb;
  int glb_jyb;

  Rect recti,recto;
  recvcomm recv_segs;
  sendcomm send_segs;
  int *pemask;
  xypos *iepos;
  int *grp_ieind;
  int *grp2glb;
  int nwpt;//capital
  int nwps;//neighbor
  int nwpc;//inside
  int nwpo;//outside
  int nwpa;//Z
}Part;
__BEGIN_DECLS
void GPart(Part *part,int *ieind,int NCX,int NCY,int nwps,int nwpc,int nwpa,int LXN,int LYN,int adjust);
void gp_filldst(sendcomm **scs,int m);
void gp_check_recv(sendcomm* sc,int m,int RFS);
void gp_send_data(sendcomm* sc,int itemsize,int m,int RFS);
void gp_sdma_send_data(sendcomm* sc,int itemsize,int m,int RFS);
int  gp_CheckSend(Part *part);
void gp_pmap(int id,int LXN,int LYN,int *pemask);
void gp_freepart(Part *part);
__END_DECLS
