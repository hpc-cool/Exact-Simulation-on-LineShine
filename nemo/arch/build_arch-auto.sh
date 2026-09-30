#!/bin/bash
#
set -u
#
#- Choice of the options ---
#
while [ $# -gt 0 ]
do
    case $( echo $1 | tr '[:upper:]' '[:lower:]' ) in
	-h|--help) cat <<'EOF'

Usage:
------
./build_arch-auto.sh


build_arch-auto.sh should automatically find all necessary information by itself...
However the following environment variables or their corresponding optional arguments can be used
to force the default and automatic definition of the variables used by build_arch-auto.sh

Environment variables that can be defined
-----------------------------------------

  NETCDF_C_prefix : prefix of the NetCDF-C library ($NCDF_F_PREFIX/lib) and include ($NCDF_C_PREFIX/include)
	     - set it to "no" if you don't want to use NetCDF library
	     - if not defined : we use "nc-config --prefix"
  	     Can also be specified with the optional argument --NETCDF_C_prefix

  NETCDF_F_prefix : prefix of the NetCDF-Fortran library ($NCDF_F_PREFIX/lib) and include ($NCDF_F_PREFIX/include)
	     - if not defined : we use "nf-config --prefix"
	     - if nf-config not found : we use "nc-config  --flibs" to find it
	     - not used if NETCDF_C_prefix="no"
  	     Can also be specified with the optional argument --NETCDF_F_prefix

  HDF5_prefix : prefix of the HFD5 library ($HDF5_PREFIX/lib)
	     - if not defined : we use "nc-config  --cflags" and "h5pcc -showconfig" to find it
	     - not used if NETCDF_C_prefix="no"
  	     Can also be specified with the optional argument --HDF5_prefix

  XIOS_prefix : prefix of the XIOS library ($XIOS_PREFIX/lib) and include ($XIOS_PREFIX/inc)
  	     Can also be specified with the optional argument --XIOS_prefix

  OASIS_prefix : prefix of the OASIS library ($OASIS_PREFIX/lib) and
	     include ($OASIS_PREFIX/build/lib/mct and $OASIS_PREFIX/build/lib/psmile.MPI1)
  	     Can also be specified with the optional argument --OASIS_prefix

  LIBMpath : path of the m library (standard C library of basic mathematical functions)
  	     Can also be specified with the optional argument --LIBMpath

  CURLpath : path of the curl library
  	     Can also be specified with the optional argument --CURLpath

  ZLIBpath : path of the z library
  	     Can also be specified with the optional argument --ZLIBpath

  FCnemo : fortran compiler
  	   Can also be specified with the optional argument --FCnemo

  CCnemo : C compiler (used only for AGRIF conv)
  	   Can also be specified with the optional argument --CCnemo

  CPPnemo : cpp compiler
  	    Can also be specified with the optional argument --CPPnemo

  MKnemo : make 
  	   Can also be specified with the optional argument --MKnemo

  ARnemo : archiver
           Can also be specified with the optional argument --ARnemo

EOF
	   exit 0 ;;

	--netcdf_c_prefix) NETCDF_C_prefix=${2} ; shift ;;
	--netcdf_f_prefix) NETCDF_F_prefix=${2} ; shift ;;
	--hdf5_prefix)     HDF5_prefix=${2}     ; shift ;;
	--xios_prefix)     XIOS_prefix=${2}     ; shift ;;
	--oasis_prefix)    OASIS_prefix=${2}    ; shift ;;
	--libmpath)        LIBMpath=${2}        ; shift ;;
	--curlpath)        CURLpath=${2}        ; shift ;;
	--zlibpath)        ZLIBpath=${2}        ; shift ;;
	--fcnemo)          FCnemo=${2}          ; shift ;;
	--ccnemo)          CCnemo=${2}          ; shift ;;
	--cppnemo)         CPPnemo=${2}         ; shift ;;
	--mknemo)          MKnemo=${2}          ; shift ;;
	--arnemo)          ARnemo=${2}          ; shift ;;
	
	*) echo -e "\033[0;31m\nERROR: \"$1\" BAD OPTION\033[0m\n"
	   exit 2 ;;
    esac
    
    shift
done
#
#-----------------------------------------------------
#
echo_red () {
    echo
    while [ $# -gt 0 ] ; do echo -e "\033[0;31m$1\033[0m" ; shift ; done
}

echo_green () {
    echo
    while [ $# -gt 0 ] ; do echo -e "\033[0;32m$1\033[0m" ; shift ; done
}
#
echo_orange () {
    echo
    while [ $# -gt 0 ] ; do echo -e "\033[0;33m$1\033[0m" ; shift ; done
}
#
# check the exit status of which 
err_which () { which ${1} > /dev/null 2>&1 ; echo $? ; }
#
# call err_which and print and error message
chk_which () {
    if [ $( err_which ${!1} ) -ne 0 ]
    then
	echo_red "ERROR: ${2^} ${!1} not found." \
		 "       please define your ${2,} path with 'export $1=...'"
	exit 1
    fi
}
#
# find fortran wrapper
find_fortran_wrapper () {
    if   [ $( err_which mpiifort ) -eq 0 ] ; then FCnemo=mpiifort
    elif [ $( err_which mpif90   ) -eq 0 ] ; then FCnemo=mpif90
    elif [ $( err_which ftn      ) -eq 0 ] ; then FCnemo=ftn
    else
	echo_red "ERROR: we found neither \"mpiifort\" nor \"mpif90\" nor \"ftn\"." \
		 "       please define your fortran compiler path with 'export FCnemo=...'"
	exit 1
    fi
}
# find the fortran compiler associated with the fortran wrapper $FCnemo
find_fortran_compiler () {
    if   [ $( $FCnemo --version | head -n 1 | grep -ci "\(gcc\|gnu\)" ) -eq 1 ] ; then ftncomp="gnu"
    elif [ $( $FCnemo --version | head -n 1 | grep -ci        "ifort" ) -eq 1 ] ; then ftncomp="intel"
    elif [ $( $FCnemo --version | head -n 1 | grep -ci "Cray Fortran" ) -eq 1 ] ; then ftncomp="cray"
    else
	echo_red "ERROR: the fortran wrapper $FCnemo does not correspond to the gnu, the intel or the cray compiler" \
		 "       please use one of these compilers or add your compiler in $0"
	exit 1
    fi
}
#
#-----------------------------------------------------
# NetCDF
#-----------------------------------------------------
#
NETCDF_C_prefix=${NETCDF_C_prefix:-notdef}
NETCDF_F_prefix=${NETCDF_F_prefix:-notdef}
if [ "$NETCDF_C_prefix" == "no" ]
then
    echo_orange "WARNING: You chose to compile without any NetCDF Library" \
		"     1)  You must use --nonetcdf when calling makenemo" \
		"     2)  You won't be able to read/write any input/output file when running nemo" \
		"         -> You can run only the BENCH test" \
		"     3)  You can neither use XIOS nor OASIS" \
    
    NETCDF_F_prefix=no
    HDF5_prefix=no
    XIOS_prefix=no
    OASIS_prefix=no
    NCDF_INC=""
    NCDF_LIB=""
    XIOS_INC=""
    XIOS_LIB=""    
    OASIS_INC=""
    OASIS_LIB=""
else
    # find NETCDF_C_prefix and NC_CONFIG
    if [ "$NETCDF_C_prefix" == "notdef" ]              # use nc-config to define NETCDF_C_prefix
    then
	if [ $( err_which nc-config ) -ne 0 ]          # nc-config not found in $PATH
	then
	    echo_red "ERROR: nc-config not found." \
		     "       please define either your path to the NetCDF-C library with 'export NETCDF_C_prefix=...'" \
		     "       or specify you don't want to use NetCDF with 'export NETCDF_C_prefix=no'" \
	    exit 2
	fi    
       	NC_CONFIG=nc-config
	NETCDF_C_prefix=$( $NC_CONFIG --prefix  )
    else                                               # use NETCDF_C_prefix to define NC_CONFIG
	NC_CONFIG=$NETCDF_C_prefix/bin/nc-config       # assume that nc-config is in $NETCDF_C_prefix/bin
	nbok=$( ls $NC_CONFIG 2>/dev/null | wc -l )    # check if we have nc-config
	if [ $nbok -eq 0 ]
	then
	    echo_red "ERROR: nc-config not found in $NETCDF_C_prefix/bin" ; exit 2
	fi
    fi

    # find NETCDF_F_prefix and NF_CONFIG
    if [ "$NETCDF_F_prefix" == "notdef" ]              # use nf-config to define NETCDF_F_prefix
    then
	if [ $( err_which nf-config ) -ne 0 ]          # nf-config not found in $PATH
	then
	    NF_CONFIG=${NC_CONFIG/nc-config/nf-config} # try to find nf-config at the same place than nc-config
	    if [ $( err_which $NF_CONFIG ) -ne 0 ]     # nf-config not located at the same place than nc-config
	    then
		if [ "$( $NC_CONFIG --all | grep -c "has-fortran" )" -gt 0 ]  # we can use nc-config instead of nf-config?
		then
		    # do we have NetCDF-Fortran?
		    if [ "$( $NC_CONFIG --has-fortran )" != "yes" ]
		    then
			echo_red "ERROR: no netcdf-fortran " ; exit 2
		    fi
		    NF_CONFIG=$NC_CONFIG   # use nc-config instead of nf-config
		    #                        define NETCDF_F_prefix with nc-config
		    NETCDF_F_prefix=$( $NF_CONFIG --flibs | sed -e "s/.*\(-L\|-rpath,\)\([^ ]*\)\/lib  *-lnetcdff.*/\2/" )
		    [ -z $NETCDF_F_prefix ] && NETCDF_F_prefix=$NETCDF_C_prefix   # empty -> we try NETCDF_C_prefix
		else
		    echo_red "ERROR: nf-config not found." \
			     "       please define either your path to the NetCDF-Fortran library with 'export NETCDF_F_prefix=...'" \
			     "       or specify you don't want to use NetCDF with 'export NETCDF_C_prefix=no'" \
			     exit 2
		fi
	    else
		NETCDF_F_prefix=$( $NF_CONFIG --prefix  )   # use nf-config to define NETCDF_F_prefix
	    fi
	else
       	    NF_CONFIG=nf-config
	    NETCDF_F_prefix=$( $NF_CONFIG --prefix  )   # use nf-config to define NETCDF_F_prefix
	fi
    else
	NF_CONFIG=$NETCDF_F_prefix/bin/nf-config   # assume that nf-config is in $NETCDF_F_prefix/bin
	nbok=$( ls $NF_CONFIG 2>/dev/null | wc -l )   # check if we have nf-config
	if [ $nbok -eq 0 ]
	then
	    echo_red "ERROR: nf-config not found in $NETCDF_F_prefix/bin" ; exit 2
	fi
    fi

    # do we have the proper path to the netcdf library?
    nbok=$( ls $NETCDF_C_prefix/lib/libnetcdf* 2>/dev/null | wc -l )
    if [ $nbok -eq 0 ]
    then
	echo_red "ERROR: netcdf library not found in $NETCDF_C_prefix/lib" ; exit 2
    else
	echo_green "NETCDF_C_prefix=$NETCDF_C_prefix"
    fi

    # do we have NetCDF-F90 interface?
    if [ "$( $NF_CONFIG --has-f03 )" != "yes" ]
    then
	if [ "$( $NF_CONFIG --has-f90 )" != "yes" ]
	then
	    echo_red "ERROR: no netcdf-fortran F90 interface" ; exit 2
	fi
    fi
  
    # do we have the proper path to the netcdff library?
    nbok=$( ls $NETCDF_F_prefix/lib/libnetcdff* 2>/dev/null | wc -l )
    if [ $nbok -eq 0 ]
    then
	echo_red "ERROR: netcdff library not found in $NETCDF_F_prefix/lib" ; exit 2
    else
        echo_green "NETCDF_F_prefix=$NETCDF_F_prefix"
    fi

    # do we have nc4?
    if [ "$( $NC_CONFIG --has-nc4 )" != "yes" ]
    then
	echo_red "ERROR: no nc4 interface in your netcdf library" ; exit 2
    fi
    
    # do we have hdf5?
    if [ "$( $NC_CONFIG --has-hdf5 )" != "yes" ]
    then
	echo_red "ERROR: no hdf5 interface in your netcdf library" ; exit 2
    fi
    
    HDF5_prefix=${HDF5_prefix:-notdef}
    if [ "$HDF5_prefix" == "notdef" ]
    then
	# look for libhdf5 in the paths found in nc-config --cflags
	for dd in $( $NC_CONFIG --cflags | sed -e "s/-I//g" )
	do
	    nbok=$( ls $dd/../lib/libhdf5* 2>/dev/null | wc -l )
	    [ $nbok -gt 0 ] && HDF5_prefix=$( dirname $dd )
	done
    fi
    if [ "$HDF5_prefix" == "notdef" ]   # not found any libhdf5 file... try "h5pcc -showconfig"
    then
	if [ $( err_which h5pcc ) -eq 0 ]   # do we have h5pcc?
	then
	    HDF5_prefix=$( h5pcc -showconfig | grep "Installation point" | sed -e "s/.*: //" )
	    [ -z $HDF5_prefix ] && HDF5_prefix=notdef # empty -> back to "notdef"
	fi
    fi
    if [ "$HDF5_prefix" == "notdef" ]   # still not found any libhdf5 file... try $NETCDF_C_prefix
    then
	nbok=$( ls $NETCDF_C_prefix/lib/libhdf5* 2>/dev/null | wc -l )
	[ $nbok -gt 0 ] && HDF5_prefix=$NETCDF_C_prefix
    fi
    if [ "$HDF5_prefix" == "notdef" ]   # still not found any libhdf5 file...
    then
	echo_red "ERROR: HDF5_prefix not found." \
		 "       please define it with 'export HDF5_prefix=...'" ; exit 2
    fi
    
    # do we have the proper path to the hdf5 library?
    nbok=$( ls $HDF5_prefix/lib/libhdf5* 2>/dev/null | wc -l )
    if [ $nbok -eq 0 ]
    then
	echo_red "ERROR: hdf5 library not found in $HDF5_prefix/lib" ; exit 2
    else
	echo_green "HDF5_prefix=$HDF5_prefix"
    fi
    
    # do we have netcdf parallel?
    if [ "$( $NC_CONFIG --has-parallel )" != "yes" ]
    then
	echo_orange "WARNING: your netcdf library cannot write in parallel" \
		    "         you can use only 'the one_file' mode in XOS"
    fi
    if [ $( err_which h5pcc ) -eq 0 ]
    then
	if [ $( h5pcc -showconfig | grep "Parallel HDF5" | grep -c yes ) -ne 1 ]
	then
	    echo_orange "WARNING: your hdf5 library was configured without --enable-parallel" \
			"         you cannot use only the 'one_file' mode in XOS"
	fi
    fi
    
    # curl
    CURLpath=${CURLpath:-notdef}
    if [ "$( $NC_CONFIG --dap )" == "yes" ]
    then
	if [ "$CURLpath" == "notdef" ]
	then
	    if [ $( err_which curl ) -eq 0 ]
	    then
		CURLpath="-L$( dirname $( ls -1 $( dirname $( dirname $( which curl ) ) )/*/libcurl* | head -n 1 ) )"
		echo_green "CURLpath=$CURLpath"
	    else
		CURLpath=""
	    fi    
	fi
	CURLlib="${CURLpath:-""} -lcurl"
    else
	CURLlib=""
    fi
    
    # Zlib: compression library
    Zlib="${ZLIBpath:-""} -lz"

    # libM: standard C library of basic mathematical functions
    libM="${LIBMpath:-""} -lm"
    
    # NCDF_INC and NCDF_LIB
    NCDF_INC="-I%NCDF_F_PREFIX/include -I%NCDF_C_PREFIX/include"
    #NCDF_LIB="-L%NCDF_F_PREFIX/lib -lnetcdff -L%NCDF_C_PREFIX/lib -lnetcdf -Wl,-rpath=%NCDF_F_PREFIX/lib -L%HDF5_PREFIX/lib -lhdf5_hl -lhdf5 $CURLlib $Zlib $libM -Wl,-rpath=%HDF5_PREFIX/lib"
    NCDF_LIB="-L%NCDF_F_PREFIX/lib -lnetcdff -L%NCDF_C_PREFIX/lib -lnetcdf -L%HDF5_PREFIX/lib -lhdf5_hl -lhdf5 $CURLlib $Zlib $libM"

fi
#
#-----------------------------------------------------
# fortran Compileur
#-----------------------------------------------------
#
nompi=0 # default
FCnemo=${FCnemo:-notdef}
if [ "$FCnemo" == "notdef" ]
then
    if [ "$NETCDF_C_prefix" == "no" ]
    then
	find_fortran_wrapper
    else
	FCnemo=$( $NC_CONFIG --fc ) # use nc-config to find the fortran compiler
	[ -z $FCnemo ] && find_fortran_wrapper # if $FCnemo is empty try to find it...
	chk_which FCnemo "fortran compiler"

	# does the fortran compiler start with the 3 letters mpi or ftn
	ismpi=$( basename $FCnemo | cut -c 1-3 | grep -Ec 'mpi|ftn' )
	if [ $ismpi -eq 0 ]
	then
	    echo_orange "WARNING: the fortran compiler provided by $( basename $NF_CONFIG ), \"$FCnemo\", is not starting with \"mpi\" or \"ftn\"." \
			"         we look for mpiifort, mpif90 or ftn..."
	    FCnemo_org=$FCnemo
	    find_fortran_compiler
	    ftncomp_org=$ftncomp
	    find_fortran_wrapper
	    find_fortran_compiler
	    if [ "$ftncomp" != "$ftncomp_org" ]
	    then
		echo_red "ERROR: \"$FCnemo\" and \"$FCnemo_org\" are refering to different compilers: \"$ftncomp\" and \"$ftncomp_org\"" \
			 "       please define your fortran compiler path with 'export FCnemo=...'"
		exit 1
	    fi
	fi
    fi
else
    ismpi=$( basename $FCnemo | cut -c 1-3 | grep -Ec 'mpi|ftn' ) # the fortran compiler does not start with mpi or ftn...
    if [ $ismpi -eq 0 ]
    then
	echo_orange "WARNING: the fortran compiler you provided does not start with mpi or ftn." \
		    "         you must compile nemo with key_mpi_off and you cannot run nemo in parallel"
	nompi=1
    fi
fi
echo_green "FCnemo=$FCnemo"

find_fortran_compiler
#
#-----------------------------------------------------
# C Compileur (not critical, needed only for the conv of AGRIF)
#-----------------------------------------------------
#
# first guess
case "$ftncomp" in
    intel) guess=${CCnemo:-icc}	;;
    gnu)   guess=${CCnemo:-gcc}	;;
    cray)  guess=${CCnemo:-cc}	;;
    *)     guess=${CCnemo:-""}	;;
esac
listcc="$guess gcc cc icc"  # try $guess and other usuals C compilers... 
CCnemo="notfound"
for ii in $listcc
do
    [ $( err_which $ii ) -eq 0 ] && [ "$CCnemo" == "notfound" ] && CCnemo=$ii
done
if [ "$CCnemo" == "notfound" ]
then
    echo_red "ERROR: we found neither \"gcc\" nor \"icc\" nor \"cc\"." \
	     "       please define your C compiler with 'export CCnemo=...'"
    exit 1
fi
#
#-----------------------------------------------------
# CPP
#-----------------------------------------------------
#
CPPnemo=${CPPnemo:-cpp}
chk_which CPPnemo "pre-compiler"
echo_green "CPPnemo=$CPPnemo"

[ "$ftncomp" == "gnu" ] || [ "$ftncomp" == "cray" ] && CPPnemo="$CPPnemo -Dkey_nosignedzero"
[ $nompi     -eq  1   ] && CPPnemo="$CPPnemo -Dkey_mpi_off"
#
#-----------------------------------------------------
# FCFLAGS
#-----------------------------------------------------
#
case "$ftncomp" in
    intel)
	PROD_FCFLAGS="-i4 -r8 -O3 -fp-model strict -xHost -fno-alias"
	DEBUG_FCFLAGS="-i4 -r8 -g -O0 -debug all -traceback -fp-model strict -ftrapuv -check all,noarg_temp_created -fpe-all0 -ftz -init=arrays,snan,huge"
	echo_orange "WARNING: We assume you will execute NEMO on the same machine you compiled it. " \
		    "         If it is not the case, replace -xHost by the appropiate option in PROD_FCFLAGS"
	;;
    gnu)
	PROD_FCFLAGS="-fdefault-real-8 -O3 -march=native -funroll-all-loops -fcray-pointer -ffree-line-length-none" 
	DEBUG_FCFLAGS="-fdefault-real-8 -Og -g -fbacktrace -funroll-all-loops -fcray-pointer -ffree-line-length-none -fcheck=all -finit-real=nan -ffpe-trap=invalid,zero,overflow -ffpe-summary=invalid,zero,overflow"
	rev=$( $FCnemo --version | head -n 1 | sed -e "s/.* \([0-9]*\).*/\1/" )
	# if gfortran version >= 10 : add -fallow-argument-mismatch 
	if [ $rev -ge 10 ]
	then
	    PROD_FCFLAGS="$PROD_FCFLAGS -fallow-argument-mismatch" 
	    DEBUG_FCFLAGS="$DEBUG_FCFLAGS -fallow-argument-mismatch"
	fi
	echo_orange "WARNING: We assume you will execute NEMO on the same machine you compiled it. " \
		    "         If it is not the case, replace -march=native by the appropiate option in PROD_FCFLAGS"
	;;
    cray)
	PROD_FCFLAGS="-s real64 -s integer32 -O2 -hflex_mp=intolerant"
	DEBUG_FCFLAGS="-s real64 -s integer32 -Ovector0 -hfp0 -O0 -hflex_mp=intolerant -e mCI -G0 -m2 -rl -Rcdsp -N1023"
	;;
    *)
	echo_red "ERROR: compilation options for $ftncomp are not defined..."
	exit 1
	;;
esac
#
#-----------------------------------------------------
# CFLAGS
#-----------------------------------------------------
#
case "${ftncomp}" in
    gnu)
        PROD_CFLAGS="-O0 -fcommon"
        ;;
    *)
        PROD_CFLAGS="-O0"
        ;;
esac
#
#-----------------------------------------------------
# XIOS and OASIS
#-----------------------------------------------------
#
if [ "$NETCDF_C_prefix" != "no" ]
then

    # XIOS
    XIOS_prefix=${XIOS_prefix:-notdef}
    if [ "$XIOS_prefix" == "notdef" ]
    then
	echo_orange "WARNING: XIOS_prefix not specified" \
		    "         either define it with 'export XIOS_prefix=...'" \
 		    "         or don't use key_xios when compiling nemo (-> you will not be able to use XIOS)."
	XIOS_INC=""
	XIOS_LIB=""
    else
	XIOS_INC="-I%XIOS_PREFIX/inc"
	#XIOS_LIB="-L%XIOS_PREFIX/lib -lxios -lstdc++ -Wl,-rpath=%XIOS_PREFIX/lib"
	XIOS_LIB="-L%XIOS_PREFIX/lib -lxios -lstdc++"
    fi

    # OASIS
    OASIS_prefix=${OASIS_prefix:-notdef}
    if [ "$OASIS_prefix" == "notdef" ]
    then
	echo_orange "WARNING: OASIS_prefix not specified" \
		    "         either define it with 'export OASIS_prefix=...'" \
		    "         or don't use key_oasis3 when compiling nemo (-> you will not be able to run coupled simulations)."
	OASIS_INC=""
	OASIS_LIB=""
    else
	OASIS_INC="-I%OASIS_PREFIX/build/lib/mct -I%OASIS_PREFIX/build/lib/psmile.MPI1"
	OASIS_LIB="-L%OASIS_PREFIX/lib -lpsmile.MPI1 -lmct -lmpeu -lscrip"    
    fi
fi
#
#-----------------------------------------------------
#   make
#-----------------------------------------------------
#
MKnemo=${MKnemo:-make}
mkerr=$( err_which make )
gmkerr=$( err_which gmake )

if [ $( err_which $MKnemo ) -ne 0 ]   # no make found
then
    if [ $mkerr -eq 0 ]   # make found
    then
	echo_orange "WARNING: $MKnemo not found, we switch to make"
	MKnemo=make
    else
	if [ $gmkerr -ne 0 ]   # gmake found
	then
	    echo_orange "WARNING: $MKnemo not found, we switch to gmake"
	    MKnemo=gmake
	else
	    echo_red "ERROR: gnu make not found, please define its path with 'export MKnemo=...'"
	    exit 1
	fi
    fi
fi
#
isgnu=$( $MKnemo --version | head -n 1 | grep -ic gnu )   # make is gmake
if [ $isgnu -ne 1 ]
then
    if [ $gmkerr -ne 0 ]   # gmake found
    then
	echo_orange "WARNING: $MKnemo is not gnu make, we switch to gmake"
	MKnemo=gmake
    else
	echo_red "ERROR: $MKnemo is not gnu make and gmake not found." \
		 "       please define gnu make path with 'export MKnemo=...'"
	exit 1
    fi
fi
echo_green "MKnemo=$MKnemo"
#
#-----------------------------------------------------
# ar command
#-----------------------------------------------------
#
ARnemo=${ARnemo:-ar}
chk_which ARnemo "archiver"
echo_green "ARnemo=$ARnemo"
#
#
#-----------------------------------------------------
#-----------------------------------------------------
# write arch file
#-----------------------------------------------------
#-----------------------------------------------------
#
#
archname=arch-auto.fcm
cat > $archname << EOF
#
# This arch file was automatically created by $0
# $( date ) 
#
%NCDF_C_PREFIX       $NETCDF_C_prefix
%NCDF_F_PREFIX       $NETCDF_F_prefix
%HDF5_PREFIX         $HDF5_prefix
%XIOS_PREFIX         $XIOS_prefix
%OASIS_PREFIX        $OASIS_prefix

%NCDF_INC            $NCDF_INC
%NCDF_LIB            $NCDF_LIB
%XIOS_INC            $XIOS_INC
%XIOS_LIB            $XIOS_LIB
%OASIS_INC           $OASIS_INC
%OASIS_LIB           $OASIS_LIB

%CPP	             $CPPnemo
%FC                  $FCnemo 
%FCFLAGS             $PROD_FCFLAGS
### comment out the following line if you want to use the debugging compilation options
#%FCFLAGS            $DEBUG_FCFLAGS
%FFLAGS              %FCFLAGS
%LD                  %FC
%LDFLAGS             -Wl,-rpath=%HDF5_PREFIX/lib -Wl,-rpath=%NCDF_F_PREFIX/lib -Wl,-rpath=%XIOS_PREFIX/lib
%FPPFLAGS            -P -traditional
%AR                  $ARnemo
%ARFLAGS             rs
%MK                  $MKnemo
%USER_INC            %XIOS_INC %OASIS_INC %NCDF_INC
%USER_LIB            %XIOS_LIB %OASIS_LIB %NCDF_LIB

%CC                  $CCnemo
%CFLAGS              $PROD_CFLAGS
EOF
#
# Additional module search command for Cray Fortran to enable successful parallel builds
if [ $ftncomp == "cray" ] ; then
    echo "bld::tool::fc_modsearch -J" >> $archname
fi
#
echo
echo "Content of the created $archname:"
echo "-------------------------------------"
echo
cat $archname
echo
