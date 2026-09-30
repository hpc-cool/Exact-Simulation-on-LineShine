//typedef struct MPI_Fint MPI_FComm;
#include "irrp.h"
#define MPICOMMF2c(comm) MPI_Comm_f2c(comm)
#define NEWN(T,N)    (T*)malloc(sizeof(T)*(N))
#define NEWPN(p,N)   p=((typeof(p)) malloc(sizeof(*(p))*(N)))
#define NEWZN(p,N)   memset( NEWPN(p,N),0,sizeof(*(p))*(N))
#define RENEWN(p,N)  p=((typeof(p)) realloc(p,sizeof(*(p))*(N)))
#define RENEWZN(p,N) memset(RENEWN(p,N),0,sizeof(*(p))*(N))
#define ZERON(p,N)   memset(p,0,sizeof(*(p))*(N))
#define CPYN(p,s,N)  memcpy(p,s,sizeof(*(p))*(N))
#define VFREE(p)     if(p)free(p);p=NULL
typedef unsigned char BYTE ;
typedef struct Rect Rect;
typedef struct mpipacket{
  int bsize ,pos ,dsize ,*buf;
}mpipacket;
typedef struct PI_Pos{
  int  i,j;
}PI_Pos;
typedef struct PI_NB8{
  int  r;         // Right neighbor point.
  int  ur;        // Up-right neighbor point.
  int  u;         // Up neighbor point.
  int  ul;        // Up-left neighbor point.
  int  l;         // Left neighbor point.
  int  dl;        // Down-left neighbor point.
  int  d;         // Down neighbor point.
  int  dr;        // Down-right neighbor point.
}PI_NB8;

typedef struct SRBuf{
  int id;
  int rn;
  int *pnts;
} SRBuf;
typedef struct PI_Rsbuf{
  int  nrb;       // Size of receive buffer.
  int  nsb;       // Size of send buffer.
  SRBuf rbuf;      // (nrb), buffer for receiving.
  SRBuf sbuf;      // (nsb), buffer for sending.
}PI_Rsbuf;

typedef struct PI_Rsinfo{
    int  id;                           // ID of the destination PE.
    int  sn;                           // Number of sending points.
    int  rn;                           // Number of receiving points.
    int *  sspnts;           // List of sending   points for serial part.
    int *  srpnts;           // List of receiving points for serial part.
                             // on part End :spnts==sspnts rpnts==srpnts
    int *  mspnts;           // List of sending   points for matrix part.
    int *  mrpnts;           // List of receiving points for matrix part
    int *  spnts;            // List of sending   points for current part.
    int *  rpnts;            // List of receiving points for current part
}PI_Rsinfo;

typedef struct PI_Rplist{
  int gis,i,j,pe,flg;//flg,1:down,2:right,4:left,8,up
}PI_Rplist;
typedef struct PI_Ginfo PI_Ginfo;
typedef struct PI_Part{
  PI_Ginfo*pgsi;
  int  snpc;            // computing points of current partion for serial mode.
  int  snp;             // points of current patition for serial mode.
  int  nps;             // sending pnts of current partition.
  int  npr;             // receiving points 
  int  npc;             // conputing points := snpc for serial mode,=nxy for rectangle mode
  int  np;              // points for  := snp for serial mode,=for rectangle mode
  int  iblv ;           // The start of sequentializing: from 1 or 0.
  int  nx, ny;          // Size of matrix: first and second dimension.
  int  nxy;             // nxy = nx * ny.
  /*--------------------------------------------------------------------------------------------
   * pnttype:
   * serial ,nodump :(0:np),0:space/land,1:open bound             3:calc,4:recv1, 4+:recv1+
   * serial ,dump   :(0:np),0:space/land,1:open bound 2:dump pnts,3:calc,4:recv1, 4+:recv1+
   * matrix : (0:im*jm)
   *-------------------------------------------------------------------------------------------*/
  int *ijbe;
  int *mpnttype;         // Point type for matrix part
  int *calpnts;          // (snp), ! for matrix part
                         // calpnts SN;
                         //      on parttype=0,noneeded,calpnts[1:snpc]=1:snpc
                         //      othertype Need to Set for gather/scatter
                         //   sort:nps,npi=snpc-nps,npr1,npr2 ...
                         //  needed for gather/scatter
  Rect  recti;         // rectangle of the partion only for computing points.
  Rect  recto;         // rectangle of the partion for all points.
  int  ib ,ie ,jb ,je;
  //ZNOUSE int *  ij2s
  int *  ij2s;           // (1-halosize:im+halosize, 1-halosize:jm+halosize)
  int  nids;             // Number of PEs for sending or receiving.
  int *ppos;             // (0:np)
  PI_NB8*nb8;            //(0:np), to record neighbor points for each points.
                         // neibhor part info needed for exchange
  PI_Rsinfo*rsinfo;      // (nids)
  PI_Rsbuf *rsbufs;      // (nids), buffer for recv and send.

  MPI_Request *requests; // (nids)
  MPI_Request *requestr; // (nids)
  MPI_Status *statuss;   // (nids)
  MPI_Status *statusr;   // (nids)
}PI_Part;
typedef struct PI_Ginfo{  // Self-defined type to record partition information.
  MPI_Comm  mpi_comm ;
  int  parttype   ;     // Default is for serial, 3 for matrix.
                        // model part mem shape,
                        //       1: serial without dump points;
                        //       2: serial with dump points,for nopack sending;
                        //       3: rectangle shape
  int  balance    ;     // Type of balance. Default is 0. 1: absolute balance.
  int  npe        ;     // Numbers of PEs.
  int  pid        ;     // Current PE id.
  int  halosize   ;     // size of halo (outer boundary)
  int  cycle_flag ;     // Cycle type for 4 boundaries of grid matrix.
  int  im, jm     ;     // Size of grid matrix, first & second dimension.
  int  ijm        ;     // Maximum computing points.
  int  sumnp      ;     // Sum of all points.
  double  avenp   ;     // Averaged points number for each PE.
  int  gnpc       ;     // Total computing points. It may include dump points.
  int *npcs       ;     // (npe)
  int *mask       ;     // (im, jm)
  Rect*rects      ;     // (npe)
                        // Temp Array no needed when partion end.
}PI_Ginfo ;
extern PI_Ginfo *gsi;//={0};
extern PI_Part  *cpart;
