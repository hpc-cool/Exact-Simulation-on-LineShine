!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~C
!----                LAGFD-WAM source program                          ----C
!----                       Version 2000                               ----C
!----   Developed by wave modelling group in LAGFD,FIO(QINGDAO),SOA    ----C
!----             Modified by Yang Yongzeng in 2000                    ----C
!----             Modified by Zhao Wei in 2007                         ----C
!----                                                                  ----C
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~C
#include "wavedeff.h"
program Masnum_WAM
  use wavemdl_mod
  use debughlp_mod
  IMPLICIT NONE
  !call system("/bin/sh tenv");
  !call initcpu
  !call c_checkstructs ;
  !call wav_stdio
  call startinit
  call RunWaveMdl
end program Masnum_WAM !}
