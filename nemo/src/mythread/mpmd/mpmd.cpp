#include "mpmd.h"
#include <string.h>
#include <stdlib.h>
int mpmd_global_rank = -1;
int mpmd_global_size = -1;
int mpmd_local_rank = -1;
int mpmd_local_size = -1;
int  comm_type_size = 0;
MPI_Comm MPI_COMM_LOCAL = MPI_COMM_NULL; // MPI通信子初始化为空
struct ZComm *tags;

int mpmd_get_mpmd_global_rank_() { return mpmd_global_rank; }
int mpmd_get_mpmd_global_size_() { return mpmd_global_size; }
int mpmd_get_mpmd_global_rank()  { return mpmd_global_rank; }
int mpmd_get_mpmd_global_size()  { return mpmd_global_size; }

int mpmd_get_mpmd_local_rank_() { return mpmd_local_rank; }
int mpmd_get_mpmd_local_size_() { return mpmd_local_size; }
int mpmd_get_mpmd_local_rank()  { return mpmd_local_rank; }
int mpmd_get_mpmd_local_size()  { return mpmd_local_size; }


MPI_Fint mpmd_get_global_comm_(){ return MPI_Comm_c2f(MPI_COMM_WORLD); }
MPI_Comm  mpmd_get_global_comm(){
    return MPI_COMM_WORLD;
}
// prog_tag 要求从0开始连续递增
MPI_Fint mpmd_start_(int *prog_tag) { return  MPI_Comm_c2f(mpmd_start(*prog_tag) );}

MPI_Comm mpmd_start(int prog_tag)
{
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpmd_global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpmd_global_size);
    MPI_Comm_split(MPI_COMM_WORLD, prog_tag, mpmd_global_rank, &MPI_COMM_LOCAL);
    MPI_Comm_rank(MPI_COMM_LOCAL, &mpmd_local_rank);
    MPI_Comm_size(MPI_COMM_LOCAL, &mpmd_local_size);
    
    int *comm_types = (int*)malloc(mpmd_global_size * sizeof(int));
    int *all_mpmd_local_ranks=  (int*)malloc(mpmd_global_size * sizeof(int));

    MPI_Gather(&prog_tag, 1, MPI_INT, comm_types, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&mpmd_local_rank, 1, MPI_INT, all_mpmd_local_ranks, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(comm_types, mpmd_global_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(all_mpmd_local_ranks, mpmd_global_size, MPI_INT, 0, MPI_COMM_WORLD);

    tags = (struct ZComm*)malloc(mpmd_global_size * sizeof(struct ZComm));
    memset(tags,0,sizeof(int)*mpmd_global_size);
    for(int i = 0; i < mpmd_global_size; i++){
      int tag = comm_types[i];
      int j = 0;
      for( ; j < comm_type_size; j++) {
        if(tags[j].tag == tag) {
            if(all_mpmd_local_ranks[i] == 0) {
                tags[j].rootid = i;
            }
            tags[j].nrank++;
            break;
        }
      }
      if(j == comm_type_size){
        tags[comm_type_size].tag = tag;
         tags[j].nrank = 1;
        comm_type_size++;
      }
    }
    tags = (struct ZComm*) realloc(tags, comm_type_size  * sizeof(struct ZComm));  // 合法：同属 C 内存体系
    MPI_Barrier(MPI_COMM_WORLD);
    free(comm_types);
    free(all_mpmd_local_ranks);
    return MPI_COMM_LOCAL;
}

int mpmd_get_comm_pid_(int *prog_tag) { return mpmd_get_comm_pid(*prog_tag); }

int mpmd_get_comm_pid(int prog_tag) {
    
    for(int i = 0; i < comm_type_size; i++)
    {
        if( prog_tag == tags[i].tag){
            return  tags[i].rootid;
        }
    }
    MPI_Abort(MPI_COMM_WORLD, -555);
    return 0;
  
}
MPI_Fint  mpmd_get_local_comm_() { return MPI_Comm_c2f(MPI_COMM_LOCAL);; }
MPI_Comm  mpmd_get_local_comm(){
    return MPI_COMM_LOCAL;
}
void mpmd_end_() { mpmd_end(); }
void mpmd_end()
{
    if (tags != nullptr) {
       free(tags);
        tags = nullptr;  
    }
    MPI_Comm_free(&MPI_COMM_LOCAL);
    MPI_Finalize();
}
