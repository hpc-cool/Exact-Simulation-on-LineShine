# NEMO coupling with the Biogeochemical Flux Model (BFM)

## What is the BFM?

The Biogeochemical Flux Model (BFM) is a numerical model for the simulation of the dynamics of major biogeochemical properties in marine ecosystems (see www.bfm-community.eu). BFM is open source software freely available under the GNU Public License. 
The model can be used in standalone mode to simulate a 0-D system or coupled with other OGCM.
The coupling with NEMO is maintained by CMCC as part of the NEMO System Team activity.

## How to get the BFM code

Access to the code is provided trough the BFM website http://www.bfm-community.eu along with instructions on how to install and use it within the Documentation `Quick Guide`. 
It is recommended to run the STANDALONE test cases before using the NEMO-BFM coupled system.

## Compile NEMO with BFM 

NEMO-BFM is compiled from the BFM configuration script exploiting the NEMO FCM compilation environment. This is done to allow BFM users to create new configurations in NEMO that are not part of the NEMO standard distribution code.
The BFM configuration shipped with NEMO is `GYRE_BFM` (described in next section).

Make sure to define in your shell enviroment the following variables with code root path:
- `BFMDIR`, pointing to the root of BFM source code
- `NEMODIR`, pointing to the root of BFM source code

Check that the appropriate ARCHFILE used for the NEMO compilation with makenemo is associated to the ARCH field within the configuration file of the selected BFM preset, e.g. `$BFMDIR/build/configurations/GYRE_BFM/configuration`

Here below an example of the commands sequence for `GYRE_BFM` preset (-p):
```
$> export BFMDIR=/home/user/bfm
$> export NEMODIR=/home/user/nemo
$> cd $BFMDIR/build
$> ./bfm_configure.sh -gcd -p GYRE_BFM
```

The script will generate (-g) the BFM code, then launch makenemo for compilation (-c) and create the run directory (-d) in $BFMDIR/run.

To get information on how to use the BFM configuration script execute the following:
```$> ./bfm_configure.sh -h```

## GYRE_BFM standard configuration

The distributed standard test case is GYRE_BFM, a version of GYRE with a full-blown BFM. 
It is a demonstration simulation and it is not meant to produce any published result. GYRE_BFM runs with analytical input data only.
The namelists for the BFM are not distributed with NEMO but are generated directly by the BFM, in directory `$BFMDIR/run/gyre_bfm`. 
The generation of the BFM namelist also copy the required NEMO namelist and namelist_top files to this directory. 
This is why there are no namelist files found in the standard run directory `$NEMODIR/cfgs/GYRE_BFM/EXPREF` 

Please refer to the README file in the preset directory for more information.

## Contacts

Please visit www.bfm-community.eu for further informations and address any technical query to the BFM System Team `bfm_st@lists.cmcc.it`
