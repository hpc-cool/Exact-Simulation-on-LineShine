MODULE wet_dry

   !! includes updates to namelist namwad for diagnostic outputs of ROMS wetting and drying

   !!==============================================================================
   !!                       ***  MODULE  wet_dry  ***
   !! Wetting and drying includes initialisation routine and routines to
   !! compute and apply flux limiters and preserve water depth positivity
   !! only effects if wetting/drying is on ( ln_wd_dl==.true. )
   !!==============================================================================
   !! History :  3.6  ! 2014-09  ((H.Liu)  Original code
   !!                 ! will add the runoff and periodic BC case later
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   wad_init      : initialisation of wetting and drying
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce  , ONLY: ln_rnf   ! surface boundary condition: ocean
   USE sbcrnf         ! river runoff 
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! MPP library
   USE timing         ! timing of the main modules

   IMPLICIT NONE
   PRIVATE

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! critical depths,filters, limiters,and masks for  Wetting and Drying
   !! ---------------------------------------------------------------------

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   wdmask   !: u- and v- limiter 
   !                                                           !  (can include negative depths)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   wdramp, wdrampu, wdrampv !: for hpg limiting

   LOGICAL,  PUBLIC  ::   ln_wd_dl    !: Wetting/drying dl activation switch (T:on,F:off)
   REAL(wp), PUBLIC  ::   rn_wdmin0   !: depth at which wetting/drying starts
   REAL(wp), PUBLIC  ::   rn_wdmin1   !: minimum water depth on dried cells
   REAL(wp), PUBLIC  ::   r_rn_wdmin1 !: 1/minimum water depth on dried cells 
   REAL(wp), PUBLIC  ::   rn_wdmin2   !: tolerance of minimum water depth on dried cells
   REAL(wp), PUBLIC  ::   rn_wd_sbcdep   !: Depth at which to taper sbc fluxes
   REAL(wp), PUBLIC  ::   rn_wd_sbcfra   !: Fraction of SBC at taper depth
   REAL(wp), PUBLIC  ::   rn_wdld     !: land elevation below which wetting/drying will be considered
   LOGICAL,  PUBLIC  ::   ln_wd_dl_bc !: DL scheme: True implies 3D velocities are set to the barotropic values at points 
                                      !: where the flow is from wet points on less than half the barotropic sub-steps  
   LOGICAL,  PUBLIC  ::  ln_wd_dl_rmp !: use a ramp for the dl flux limiter between 2 rn_wdmin1 and rn_wdmin1 (rather than a cut-off at rn_wdmin1)      
   REAL(wp), PUBLIC  ::   ssh_ref     !: height of z=0 with respect to the geoid; 

   LOGICAL,  PUBLIC  ::   ll_wd = .FALSE. !: Wetting/drying activation switch (ln_wd_dl) <- default def if wad_init not called

   PUBLIC   wad_init                  ! initialisation routine called by step.F90

   !! * Substitutions
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE wad_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE wad_init  ***
      !!                    
      !! ** Purpose :   read wetting and drying namelist and print the variables.
      !!
      !! ** input   : - namwad namelist
      !!----------------------------------------------------------------------
      INTEGER  ::   ios, ierr   ! Local integer
      !!
      NAMELIST/namwad/ ln_wd_dl, rn_wdmin0,   rn_wdmin1,    rn_wdmin2,   rn_wdld,   &
         &             ln_wd_dl_bc, ln_wd_dl_rmp, rn_wd_sbcdep,rn_wd_sbcfra
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_ref(INDEX(numnam_ref,"namwad")-1:), namwad, IOSTAT = ios, ERR = 905)
905   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namwad in reference namelist' ) 
      READ  ( numnam_cfg(INDEX(numnam_cfg,"namwad")-1:), namwad, IOSTAT = ios, ERR = 906)
906   IF( ios >  0 )   CALL ctl_nam ( ios , 'namwad in configuration namelist' )
      IF(lwm) WRITE ( numond, namwad )
      !
      IF( rn_wd_sbcfra>=1 )   CALL ctl_stop( 'STOP', 'rn_wd_sbcfra >=1 : must be < 1' )
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'wad_init : Wetting and drying initialization through namelist read'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '   Namelist namwad'
         WRITE(numout,*) '      Logical for Dir. Lim wd option   ln_wd_dl     = ', ln_wd_dl
         WRITE(numout,*) '      Depth at which wet/drying starts rn_wdmin0    = ', rn_wdmin0
         WRITE(numout,*) '      Minimum wet depth on dried cells rn_wdmin1    = ', rn_wdmin1
         WRITE(numout,*) '      Tolerance of min wet depth       rn_wdmin2    = ', rn_wdmin2
         WRITE(numout,*) '      land elevation threshold         rn_wdld      = ', rn_wdld
         WRITE(numout,*) '      T => baroclinic u,v=0 at dry pts: ln_wd_dl_bc = ', ln_wd_dl_bc     
         WRITE(numout,*) '      use a ramp for rwd limiter:  ln_wd_dl_rwd_rmp = ', ln_wd_dl_rmp
         WRITE(numout,*) '      cut off depth sbc for wd   rn_wd_sbcdep       = ', rn_wd_sbcdep
         WRITE(numout,*) '      fraction to start sbc wgt rn_wd_sbcfra        = ', rn_wd_sbcfra
      ENDIF
      IF( .NOT. ln_read_cfg ) THEN
         IF(lwp) WRITE(numout,*) '      No configuration file so seting ssh_ref to zero  '
         ssh_ref=0._wp
      ENDIF

      r_rn_wdmin1 = 1 / rn_wdmin1
      IF( ln_wd_dl ) THEN
         ll_wd = .TRUE.
         ALLOCATE( wdmask(jpi,jpj),   STAT=ierr )
         ALLOCATE( wdramp(jpi,jpj), wdrampu(jpi,jpj), wdrampv(jpi,jpj), STAT=ierr ) 
         IF( ierr /= 0 ) CALL ctl_stop('STOP', 'wad_init : Array allocation error')
      ENDIF

      !
   END SUBROUTINE wad_init

   !!==============================================================================
END MODULE wet_dry
