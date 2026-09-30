#pragma once
#include <mpi.h>
#include <stdlib.h>

__BEGIN_DECLS
extern int mpmd_global_rank;
extern int mpmd_global_size;
extern int mpmd_local_rank;
extern int mpmd_local_size;
extern MPI_Comm MPI_COMM_LOCAL;
struct ZComm{
  int tag;
  int rootid;
  int nrank;
};
extern struct ZComm *tags;
int mpmd_get_global_rank_();
int mpmd_get_global_size_();
int mpmd_get_global_rank() ;
int mpmd_get_global_size() ;

int mpmd_get_local_rank_() ;
int mpmd_get_local_size_() ;
int mpmd_get_local_rank()  ;
int mpmd_get_local_size()  ;

int mpmd_get_comm_pid_(int *prog_tag);
int mpmd_get_comm_pid(int prog_tag);

MPI_Comm  mpmd_start(int prog_tag);
MPI_Fint  mpmd_start_(int *prog_tag);

MPI_Comm mpmd_get_local_comm();
MPI_Fint mpmd_get_local_comm_();

MPI_Comm mpmd_get_global_comm();
MPI_Fint mpmd_get_global_comm_();

void mpmd_end(); 
void mpmd_end_();
__END_DECLS
