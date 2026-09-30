#ifndef WAVEDEF_H_INCLUDED
# define WAVEDEF_H_INCLUDED
# define MBVDEP   50
# define SIGMALVL 0
# define CKL 32
# define CJNTHET 24
# define CMKJ (CJNTHET*CKL)
# define NKP 8
# define KLPAT NKP
# define NTSPLIT 10
# define MAXAVHIST 16
# define IMDPSO  4
# define BVIMDPS 1
# define MAXRFB (2000*1000*1000)
# if 0 
   === Include File for Fortran and C Code
   ==== #include "wavedef.h" =======
#  define TEST_PART
#  define CP_SWF90
#  define NO_MPI
# endif
# ifdef MMTHREAD
#  define MTHREAD
# endif
# ifdef C_CALCULATE
#  define C_PROPAGAT 
#  define C_IMPLSCH 
# endif
# ifndef LOGSCURR
#  define LOGSCURR 0
# else
#  define LOGSCURR 1
# endif
# ifdef USEHBM
#  define USE_C_ALLOC 
# endif
#endif
