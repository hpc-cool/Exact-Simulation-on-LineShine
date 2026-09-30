#Mach Variable defines,can defined in environment variables
# surport: default,JNCSSW2,SW4CC,SW5CC,THGZ,THGZMIC
#for Make mode ;surport nml,nmlc,ai4s,omp,sw,mic

MACH    ?= HWARM   # in huawei machine
MACH    ?= NORM      # normal linux
NO_MPI  ?= 0				# 1:test for no mpi  
DEBUG		?= 0        # 1:for debug build
minimake?= 0				# 1:notcheck dependce,for speedup make in mini modify
defmake ?= arm3     # default tagets
HBM     ?= 1				# 1:use HBM in huawei machine
NJ      ?= 10       # make tasks, -j 1 for many warning or error
XNL     ?= 1
EXE     ?= nwam
ALOG    ?= 0
LIBNETCDF ?=-lnetcdf -lnetcdff  #netcdf lib 

#for var may have space,cannot compare
NJ      ?= 1 
minimake:= $(strip $(minimake))
MACH    := $(strip $(MACH))
NO_MPI  := $(strip $(NO_MPI))
DEBUG		:= $(strip $(DEBUG))
minimake:= $(strip $(minimake))
defmake := $(strip $(defmake))
HBM     := $(strip $(HBM))
LPID    := $(strip $(LPID))
XNL     := $(strip $(XNL))
ALOG    := $(strip $(ALOG))
NJ 			:=-j $(strip $(NJ))
ifeq ($(MACH),NORM)
	HBM=0
endif

ifeq ($(DEBUG),1)
	CDEBUG=-g -DDEBUG -O0 
	FDEBUG=-g -DDEBUG -O0 -ggdb -fcheck=bounds -fopenmp
	LDEBUG=-g 
else
	CDEBUG= -O3 -g -DNODBGINFO
	FDEBUG= -O3 -g -DNODBGINFO #-fcheck=bounds
	LDEBUG= -g
endif
ifeq ($(ALOG),1)
	DFLAG := $(DFLAG) -DALLLOG
	PR:=l$(PR)
endif
ifeq ($(HBM),1)
	DFLAG := $(DFLAG) -DUSEHBM
	PR:=h$(PR)
endif
NO_MPI := $(strip $(NO_MPI))
ifeq ($(NO_MPI),1)
	DFLAG := $(DFLAG) -DNO_MPI
	PR    := n$(PR)
endif
ifeq ($(XNL),1)
	MODEL_DEF:= $(MODEL_DEF) -DUSEXNL
	SRC      :=$(SRC) mod_mkgrid.F90  mod_xnl4v5.F90
	CPPSRC   :=$(CPPSRC) xnlcal.cpp
	PR       :=x$(PR)
endif
.PHONY:defmake all zall nml nmlc ai4s omp arm0 arm1 arm2 arm3

defmake:$(defmake)


irrp: src #develop/testting change irrp code to c code
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
irrpt:src #origin code(fortran), 
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
omp:src  #nmlc+ 1 level threads
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
#ai4s:src  #  nmlc + output for ai4s
#	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
nml:src  #origin code(fortran), 
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
nmlc:src  #Change the core computing code to C code
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
arm0:src  #nmlc 
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
arm1:src  #nmlc +vecterlize
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
arm2:src  #arm1+ 1 level threads
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
arm3:src  #arm1+ 2 level threads
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all

clean_irrp:
	rm -rf objs_$(PR)irrp nwam_$(PR)irrp
clean_nml:
	rm -rf objs_$(PR)nml nwam_$(PR)nml
clean_nmlc:
	rm -rf objs_$(PR)nmlc nwam_$(PR)nmlc
clean_ai4s:
	rm -rf objs_$(PR)ai4s nwam_$(PR)ai4s
clean_omp:
	rm -rf objs_$(PR)omp nwam_$(PR)omp
clean_arm:
	rm -rf objs_$(PR)arm nwam_$(PR)arm
clean_arm1:
	rm -rf objs_$(PR)arm1 nwam_$(PR)arm1
clean_arm2:
	rm -rf objs_$(PR)arm2 nwam_$(PR)arm2
clean_arm3:
	rm -rf objs_$(PR)arm3 nwam_$(PR)arm3
clean_all:    clean_nml clean_nmlc clean_omp clean_arm clean_arm1 clean_arm2 clean_arm2 clean_irrp
clean : clean_$(defmake)

#MODEL_DEF:=  -DC_CALCULATE -DNPROPAGAT #  -DOCHECK

#default :
FC = mpif90
CC = mpicc
CXX =g++
CPP =cpp
LD = mpif90
DFLAG := $(DFLAG) 
#-DFORTRANUNDERSCORE #FORTRANUNDERSCORE is default
FFLG = -cpp -ffree-line-length-0 #-fallow-argument-mismatch
CFLG =
FMATH=
#OMPFLG=-openmp

MACH_State=default

ifeq ($(MACH),NORM)
	MACH_State=$(MACH)
	#The directory where the netcdf library is located
	NETCDF_BASE ?=/usr
	# fortran is suport type(*) 
	#GF13 := $(shell echo 'subroutine ta;integer ia(10) ; real(4) fa(10);call abf(ia);call abf(fa);end' >s.F90;if ! $(FC) -c s.F90 2>/dev/null;then echo '-DGF13';fi;rm -f s.F90 )
	#GF13 := $(shell echo 'subroutine t;interface;subroutine a(p);type(*)::p(..);end;end interface;end'>s.F90;if $(FC) -c s.F90 2>/dev/null;then echo '-DGF13';fi;rm -f s.F90 )
	GF13:=-DGF13
	FFLG := $(FFLG) -g $(GF13) -fallow-argument-mismatch #-fallow-invalid-boz,-std=legac
	dsts=nml nmlc omp #ai4s irrp #irrpt
endif
ifeq ($(MACH),HWARM)
	hpctooldir ?= /pacific_fs/HPCKit/latest
	hpcthirddir?= /pacific_fs/HPCKit/third
	NETCDF_BASE ?=$(HOME)/local#$(NETCDF_CLANG_PATH)
	MACH_State=$(MACH)
	CPU=-mcpu=hip11
	FC=mpif90 $(CPU)
	#export NAN_REPORT_CONFIG=error+detail
	#-mcpu=linxicore6100
	CFLG :=$(CFLG) $(CPU) #-march=armv9.2-a+sve2
	CFLG :=$(CFLG) -I$(HOME)/local/include 
	#-finstrument-functions
	CC =mpicc
	CXX =$(CC)
	CPP=clang-cpp
	LD = mpif90  -mcpu=hip11 -lfpe --rtlib=compiler-rt
	FFLG := $(FFLG) -g #-fp-detect-exceptions
	SLIBS= -lpthread -lksvml -lm -lsdma_dk  -L$(HOME)/local/lib 
	#dsts=nml nmlc omp arm1 arm2 arm3 #ai4s 
	dsts=arm0 arm1 arm2 arm3 #ai4s 
endif
NETCDF_BASE ?=/usr


ifdef MACH
	ifeq ($(MACH_State),default)
	  $(error MACH SET error ,$(MACH) is not surport ,you can change it or add it)
	endif
endif
zall: $(dsts)

LPID := $(strip $(LPID))

OBJD=objs_$(LPID)

CINC:= $(CINC) svars.h std.h
#CSRC : *.c source file 
#CPPSRC : *.cpp source file 
#CINC   : depends include file for c/cpp,
#BASEP  : based on LPID,add to VPATH
ifeq ($(LPID),nml)
	MODEL_DEF:=  $(MODEL_DEF) 
	CSRC  :=$(CSRC) 
	CPPSRC:= $(CPPSRC) 
	BASEP  =omp
endif
ifeq ($(LPID),nmlc)
	MODEL_DEF:=  $(MODEL_DEF) -DC_CALCULATE
	CINC   := $(CINC) zarm.h
	CSRC  :=$(CSRC) svars.c 
	CPPSRC:= $(CPPSRC) implsch.cpp  propagat.cpp output_cal.cpp boundary.cpp
	BASEP  =omp
endif
ifeq ($(LPID),arm0)
	MODEL_DEF :=  $(MODEL_DEF) -DC_CALCULATE 
	CINC   := $(CINC) zarm.h
	CSRC   :=$(CSRC) svars.c
	CPPSRC := $(CPPSRC) propagat.cpp implsch.cpp output_cal.cpp boundary.cpp
	BASEP=omp
endif
ifeq ($(LPID),ai4s)
	MODEL_DEF :=  $(MODEL_DEF) -DC_CALCULATE
	CSRC   :=$(CSRC) svars.c
	CPPSRC := $(CPPSRC) implsch.cpp  propagat.cpp output_cal.cpp boundary.cpp
	BASEP   =omp
endif
ifeq ($(LPID),arm1)
	MODEL_DEF :=  $(MODEL_DEF) -DC_CALCULATE -DARM_VEC
	CINC   := $(CINC) zarm.h
	CSRC   :=$(CSRC) svars.c
	CPPSRC := $(CPPSRC) propagat.cpp implsch.cpp output_cal.cpp boundary.cpp
	BASEP=omp
endif
ifeq ($(LPID),arm2)
	MODEL_DEF :=  $(MODEL_DEF) -DC_CALCULATE -DMTHREAD -DARM_VEC
	CINC   := $(CINC) zarm.h
	CINC   := $(CINC) mthread.h
	CSRC   :=$(CSRC) svars.c mthread.c
	CPPSRC := $(CPPSRC) propagat.cpp implsch.cpp output_cal.cpp boundary.cpp
	BASEP=omp
endif
ifeq ($(LPID),arm3)
	MODEL_DEF :=  $(MODEL_DEF) -DC_CALCULATE -DMMTHREAD -DARM_VEC
	CINC:= $(CINC) zarm.h
	CINC:= $(CINC) mthread.h
	CINC:= $(CINC) gpart.h
	#FSRC = split_mod.F90
	CSRC:=$(CSRC) svars.c mthread.c
	CPPSRC:= $(CPPSRC) propagat.cpp implsch.cpp output_cal.cpp boundary.cpp gpart.cpp
	BASEP=omp
endif
ifeq ($(LPID),omp)
	#UOMPFLG=$(OMPFLG)
	MODEL_DEF :=  $(MODEL_DEF) -DC_CALCULATE -DMTHREAD
	CINC  := $(CINC) mthread.h
	CSRC  :=$(CSRC) svars.c mthread.c
	CPPSRC:= $(CPPSRC) implsch.cpp  propagat.cpp output_cal.cpp boundary.cpp
endif

ifeq ($(LPID),irrp)
	#UOMPFLG=$(OMPFLG)
	MODEL_DEF :=  $(MODEL_DEF) 
	CINC  := $(CINC) 
	CSRC  :=$(CSRC) irrp_smpi.c irrp_split.c irrp_kernal.c irrp_package.c
	CPPSRC:= $(CPPSRC)
	BASEP=omp
endif
ifeq ($(LPID),irrpt)
	#UOMPFLG=$(OMPFLG)
	MODEL_DEF :=  $(MODEL_DEF) 
	CINC   := $(CINC) 
	CPPSRC := $(CPPSRC)
	BASEP=omp
endif

DFLAG := $(DFLAG) $(MODEL_DEF)

FFLAGS = $(DFLAG) $(FDEBUG) $(FFLG) $(INCLDIR) $(UOMPFLG)
CFLAGS = $(DFLAG) $(CDEBUG) $(CFLG) $(INCLDIR) $(UOMPFLG)
LFOPTS = $(SLIBS) $(LDEBUG) -L${NETCDF_BASE}/lib $(LIBNETCDF) $(UOMPFLG)

ifneq ($(NO_MPI),1)
	SRCT = wav_mpi_mod.F90 irrp_smpi_mod.F90 irrp_split_mod.F90 irrp_kernal_mod.F90  irrp_package_mod.F90
	CSRCT = mympi.c
endif
	CSRCT := hmalloc.c zfile.c $(CSRCT) ctools.c
	SRCT :=$(SRCT) wav_cal_mod.F90 netcdf_mod.F90
	SRC :=$(SRC) debughlp_mod.F90  varcommon_mod.F90 \
	      partition_mod.F90 dataio_mod.F90 windin_mod.F90 \
	      output_cal_mod.F90 output_loc_mod.F90 output_nc_mod.F90 output_mod.F90 \
	      boundary_mod.F90 implsch_mod.F90  propagat_mod.F90 \
	      restart_mod.F90  $(FSRC) \
	      platform_init_mod.F90 waveinit_mod.F90 wavemdl_mod.F90 main.F90
	#
# msg_mod.F90
CSRC := std.c $(CSRC)

FOBJS = ${SRC:.F90=.o}
COBJS = $(CPPSRC:.cpp=.o) $(CSRC:.c=.o)
FOBJST = ${SRCT:.F90=.o}
COBJST = $(CSRCT:.c=.o)

OBJST =$(FOBJST) $(COBJST) #$(SSOBJS) $(SCOBJS)
OBJSM =$(FOBJS) $(COBJS) #$(SSOBJS) $(SCOBJS)
OBJS =$(OBJST) $(OBJSM) #$(SSOBJS) $(SCOBJS)

ifdef BASEP
	vbp=$(subst  :, ../src/,:$(BASEP))
	INCDB=$(subst  :, -I../src/,:$(BASEP))
endif
INCLDIR := -I. -I../src/$(LPID) $(INCDB) -I../src/tools -I../src $(INCLDIR)
#=Set src file Path 
VPATH = ../src/$(LPID) $(vbp) ../src/io ../src/tools ../src ../src/xnl

#default/dump target:
#check  directory is exist,and make real/dump target in objs directory,
#and copy EXE file to current directory

.SUFFIXES:
.SUFFIXES:.F90 .f90 .o .c .oo .ss .s

%.o:%.F90
	$(FC) $(FFLAGS) -c $<

%.o:%.c
	$(CC) $(CFLAGS) -c $<

%.o:%.cpp
	$(CXX) $(CFLAGS) -c $<

$(FOBJS) : ../src/wavedef.h ../src/wavedeff.h
$(COBJS) : ../src/wavedef.h ../src/wavedefc.h

#CPP1 = $(CPP) -P -w  -traditional-cpp $(GF13) # -U__FILE__ -U__LINE__
#CPPE = $(CPP) -P -w  -traditional-cpp
CPP1 = $(CPP) -P -w   $(GF13) # -U__FILE__ -U__LINE__
CPPE = $(CPP) -P -w  
# Function: -P : Inhibit generation of linemarkers.
#           -w : Avoid warning due to undefine __FILE__ and __LINE__.
#           -U__FILE__ -U__LINE__ : undefine __FILE__ and __LINE__ ,for use in future.

SEDP = sed -e 's/\\$$/~\\/g;s/\/\//~~/g'
# Function: Add char '~' before  '\' at end of line
#    where: \\ used as \ at command line.
#           $$ treated as $ in makefile.
#           $ means end of line at command of sed.

SEDN = sed -e "s/~~/\/\//g;s/ *;~/\n/g ; s/~/\&\n/g ; s/^ *&//g;s/\#\#\s*//g"
SEDE = sed -e "/^\s*$$/{ N;  /^\s*$$/D }"
# Function: (1) Replace ";~" and any spaces into a newline
#           (2) Replace "~" into "&\n"
#           (3) Delete spaces before & for the case that there
#               is none valid characters before this &.

irrp_package_mod.o:irrp_package_mod.F90
	$(SEDP) $< | $(CPP1) -D__FILE__=\"$*.F90\" | $(SEDN) | $(SEDE)  >$*.f90
	$(FC) $(FFLAGS) -c $*.f90

wav_mpi_mod.o:wav_mpi_mod.F90
	$(SEDP) $< | $(CPP1) -D__FILE__=\"$*.F90\" | $(SEDN) | $(SEDE)  >$*.f90
	$(FC) $(FFLAGS) -c $*.f90

netcdf_mod.o:netcdf_mod.F90
	$(SEDP) $< | $(CPP1) -D__FILE__=\"$*.F90\" | $(SEDN) | $(SEDE) |$(CPPE) -D__FILE__=\"$*.f90\" >$*.f90
	$(FC) $(FFLAGS) -I$(NETCDF_BASE)/include -c $*.f90

ifneq ($(minimake),1)

ifneq ($(NO_MPI),1)
irrp_split_mod.o :irrp_smpi_mod.o
irrp_kernal_mod.o:irrp_smpi_mod.o irrp_split_mod.o
irrp_package_mod.o:irrp_kernal_mod.o irrp_smpi_mod.o
partition_mod.o  :wav_mpi_mod.o irrp_split_mod.o irrp_smpi_mod.o irrp_kernal_mod.o  irrp_package_mod.o
debughlp_mod.o   : wav_mpi_mod.o
wav_mpi_mod.o    :
endif

gpart.o 				 :$(CINC) 
svars.o 				 :$(CINC)
mthread.o				 :$(CINC)
propagat.o 			 :$(CINC)
implsch.o 			 :$(CINC)
output_cal.o 		 :$(CINC)
boundary.o			 :$(CINC)
xnlcal.o  			 :$(CINC)
wav_cal_mod.o    :
debughlp_mod.o   :
varcommon_mod.o  : defcheck debughlp_mod.o  netcdf_mod.o wav_cal_mod.o
boundary_mod.o   : varcommon_mod.o
implsch_mod.o    : varcommon_mod.o partition_mod.o
propagat_mod.o   : varcommon_mod.o partition_mod.o

output_cal_mod.o : varcommon_mod.o
output_nc_mod.o  : varcommon_mod.o partition_mod.o wav_cal_mod.o netcdf_mod.o
output_loc_mod.o : varcommon_mod.o partition_mod.o wav_cal_mod.o netcdf_mod.o
output_mod.o     : varcommon_mod.o output_cal_mod.o partition_mod.o output_nc_mod.o output_loc_mod.o   netcdf_mod.o

partition_mod.o  : varcommon_mod.o
dataio_mod.o     : varcommon_mod.o partition_mod.o netcdf_mod.o wav_cal_mod.o
windin_mod.o     : varcommon_mod.o dataio_mod.o wav_cal_mod.o

restart_mod.o   : varcommon_mod.o partition_mod.o windin_mod.o
platform_init_mod.o:varcommon_mod.o implsch_mod.o propagat_mod.o windin_mod.o output_mod.o partition_mod.o  restart_mod.o

waveinit_mod.o  : varcommon_mod.o platform_init_mod.o boundary_mod.o implsch_mod.o propagat_mod.o windin_mod.o output_mod.o partition_mod.o  restart_mod.o
wavemdl_mod.o   : varcommon_mod.o waveinit_mod.o platform_init_mod.o boundary_mod.o implsch_mod.o propagat_mod.o windin_mod.o output_mod.o partition_mod.o  restart_mod.o
main.o          : debughlp_mod.o wavemdl_mod.o
mod_xnl4v5.o    : mod_mkgrid.o
else
#propagat_mod.o   : varcommon_mod.o
endif

$(OBJSM) :defcheck

#real/dump targets
all : $(EXE)

#check if defines are changed?
defcheck: DUMP
	@echo  $(CFLAGS) $(FFLAGS) >defcheck.new;cmp -s defcheck defcheck.new || { mv defcheck.new defcheck ; echo LPID= $(LPID) MACH=$(MACH) `cat defcheck`; }
DUMP:

#real target,depend defcheck & files $(OBJS)
$(EXE) : defcheck $(OBJS)  ../src/Makefile
	$(LD) -o  $@ $(OBJS)  $(LFOPTS)
	cp $(EXE) ..
bak:
	sdate=`date +"%y%m%d_%H%M"`; \
	mkdir -p bak/src$${sdate} ;\
	cp -rp src/* bak/src$${sdate};)\
	(cd ../bak;zip -r src$${sdate} src$${sdate} )
