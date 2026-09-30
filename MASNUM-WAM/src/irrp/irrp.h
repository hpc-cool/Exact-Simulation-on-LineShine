//typedef struct MPI_Fint MPI_FComm;
typedef unsigned char BYTE ;
typedef struct PI_Part PI_Part;
typedef struct PI_Ginfo PI_Ginfo;
typedef struct PI_Pos PI_Pos;
typedef struct PI_NB8 PI_NB8;
typedef struct mpipacket mpipacket;
typedef struct Rect{
  int x0,y0,x1,y1;
}Rect;

void irrp_set_pemask(int *mask_,int*pemask_, int im_, int jm_, int npe_,int pid_);
void irrp_part_init(int partmode,int  npe,int  pid,MPI_Comm mpi_comm,int*  mask,int*  pemask,int im,int jm,int halosize,int cycle_flag,int* scycle,int NS3);
void irrp_part_final(int mode);
// Set the pointers arranging method as Matrix way.
void irrp_SetPartMatrix(Rect*recto,int* ptype);
// Set the pointers arranging method as Serial way.
void irrp_SetPartSerial(int *gnpc,int * npc,int * np,PI_Pos**plist,PI_NB8**nb8,int *iblv_,int *nwps_);
void irrp_getrects(Rect*recti,Rect*recto,Rect*rects);
void irrp_output_pposg();
// Initialize the irrp package.
void irrp_init(int partmode, int npe,int  pid,MPI_Comm mpi_comm,int * mask,int im,int jm,int  halosize,int  cycle_flag,int  *scycle,int NS3);
void irrp_init_data();
// Finalize/Deinit the irrp package.
void irrp_final(int mode);
void irrp_exginf(int nst);
//public   irrp_gather;                     // Gathering info/data from all PEs to root PE.
//public   irrp_scatter;                    // Scatering info/data from root PE to all PEs.
// Initialize exchage functions
void irrp_exg_init(int *group_ind,int nvar) ;
// Set/append var into a special exchange group.
void irrp_exg_setvar(int*group_ind_,int*ivar_,const char*vname,void*var,int ivt,int km);
// Act the exchange by group.
void irrp_exg_action(int group_ind,int mode_);
// Finalize the exchange functions.
void irrp_exg_final(int group_ind);
// start the exchange by group,equal irrp_exg_action(1).
void irrp_exg_start(int group_ind);
// test mpi irecv  is end ,
// in some realizations of MPI standard ,big data packet
// completed at wait or test
int irrp_exg_check(int group_ind);
// end the exchange by group,equal irrp_exg_action(6).
void irrp_exg_end(int group_ind);
// Initialize the function of forcing scatter
void irrp_scatter_force_init(int *forceid_,double*fix_,double*fiy_,double*fox,double* foy,
                              int nix_,int niy_,int nnx_,int nny_,int npflg,int xcycle,int root_);
//irrp_scatter_force;              // Act the forcing scattering.
// Finalize the scatter_force Function.
void irrp_scatter_force_final(int forceid);
// Initialize scatter with extended pnts.
void irrp_scatter_ext_init();
  //irrp_scatter_ext;                // Scattering the info with extended pnts.
// Finalize scatter_ext Function
void irrp_scatter_ext_final();
void irrp_exg_final_all();
void irrp_scatter_force_final_all();
// Gathering info/data from all PEs to root PE.
void irrp_gather(void*lvar,void*gvar,int ivt,int km,int k,int root_);
// Scatering info/data from root PE to all PEs.
void irrp_scatter(void*gvar,void*lvar,int ivt,int km,int k,int root_);
// Scattering the info with extended pnts.
void irrp_scatter_ext(void*gvar,void*lvar,int ivt,int km, int k, int root_);
// Act the forcing scattering.
void irrp_scatter_force(int fid,void*gvar,void*lvar,int ivt,int km,int k,int root_);

void starttimer();
double Difftimer();
double wav_irtc();
int  iwalltime() ;
void irrp_abort(const char *file, int line,const char*msg);
void InitMpiPacket(mpipacket*pk,int lsize);
void VMpiPacket(mpipacket*pk,int lsize);
void SetMpiPacketDSize(mpipacket*pk, int dsize);
void bcast_packet(mpipacket*pk,int root,int pid,MPI_Comm mpicomm);
void next(int flg,int pid,int npe,MPI_Comm mpicomm);

void prev(int flg,int pid,int npe,MPI_Comm mpicomm);
void prevnext (int flg,int pid,int npe,MPI_Comm mpicomm);
void barr(int pid,int npe,MPI_Comm mpicomm);
