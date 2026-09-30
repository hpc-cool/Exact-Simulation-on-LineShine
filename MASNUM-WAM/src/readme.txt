1度分辨率单节点大约需要1780MB DDR内存
1/20 1.7*400
MemTotal:       565G
30G 32G


void c_allocate_vee_(int *res) :std.c
void c_setpointers_(double*eec,double*eet,ImplschP*ip,spc_interg*ps_vs,geo_interg*pg_vs,int*ipos12,propinf*ipos8,windvs*pwvs);

void C_SetGPar(FC_PARA *fcp,struct Iepos*iepos,int*ieind,int*nsp){
  set pointers to gppar
    InitThreadsGPar(); :svars.c
    set to pd;

void c_init_distinct_() :std.c
  initmtpar();:svars.c
  InitCPropgats();:propagat.cpp
svars.c:
void InitThreadsGPar();
  set pd;gpart
void  initmtpar():
  fillsendcomm();
void c_initimplsch_(ImplschPar *pisp_,ImplschP*pip,float*pwvs){
void c_setwind(){ // ZZZ should opt to sdma 
void setflag_(int *hist_eot,int *stop_now,int *rest_eot){
void caliacbe(){
void c_sendboundary_(int *RFG){
void c_checkboundary_(int *RFG){

void _threadmain_(HTHREADINFO ti){
  
==========================
1:link files in srcdir to here:
	exec:
		bash src/lnr
	or exec:
		ln src/Makefile .
		ln src/nwamctl.ini .
2:check var in Makefile is correct
  MACH  : HWARM:in huawei  
          NORM :other
  NO_MPI :1 : test for no mpi  
          0 : other
  HBM     :1 :use HBM in huawei
           0 : for other
	NETCDF_BASE :The directory where the netcdf library is located
               /usr: default
  LIBNETCDF  :-lnetcdf -lnetcdff:default
              -lnetcdf :for c&F lib in 1 file
  minimake:0 :default
           1 :speed for debugging
  NJ      :-j 10 :default
           -j 1  :for many compile error
  defmake :
    nml:  origin code(fortran), 
    nmlc: Change the core computing code to C code
    iai4s: nmlc + output for ai4s
    omp:  nmlc+ 1 level threads
    arm1: nmlc +vecterlize
    arm2: arm1+ 1 level threads
    arm3: arm1+ 2 level threads
    irrp: develop/testting change irrp code to c code
3:check files exists
  etopo5.nc: in dir init/
4:make
  see var defmake
  make
  make nml
  make nmlc
  make omp
  make ai4s
  or in huawei mach 
  make arm1
  make arm2
  make arm3
  for sample:
  HBM=0 NO_MPI=0 make nml
5:check parameter file 
  nwamctl.ini :in current dir
  nwam_nml :excute file
6:run 
  mpi_run -np 8 ./nwam_nml





