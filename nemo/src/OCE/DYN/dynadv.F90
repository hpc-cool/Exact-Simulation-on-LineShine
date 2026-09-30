MODULE dynadv
   !!==============================================================================
   !!                       ***  MODULE  dynadv  ***
   !! Ocean active tracers:  advection scheme control
   !!==============================================================================
   !! History :  1.0  !  2006-11  (G. Madec)  Original code
   !!            3.3  !  2010-10  (C. Ethe, G. Madec)  reorganisation of initialisation phase
   !!            3.6  !  2015-05  (N. Ducousso, G. Madec)  add Hollingsworth scheme as an option 
   !!            4.0  !  2017-07  (G. Madec)  add a linear dynamics option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_adv      : compute the momentum advection trend 
   !!   dyn_adv_init : control the different options of advection scheme
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE dynadv_cen2     ! centred flux form advection      (dyn_adv_cen2 routine)
   USE dynadv_ubs      ! UBS flux form advection          (dyn_adv_ubs  routine)
   USE dynadv_up3      ! UBS flux form advection (NEW)    (dyn_adv_up3  routine)
   USE dynkeg          ! kinetic energy gradient          (dyn_keg      routine)
   USE dynzad          ! vertical advection               (dyn_zad      routine)
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_adv       ! routine called by step module
   
   PUBLIC dyn_adv_lll       ! routine called by step module
   PUBLIC dyn_adv_init  ! routine called by opa  module
 
   !                                   !!* namdyn_adv namelist *
   LOGICAL, PUBLIC ::   ln_dynadv_OFF   !: linear dynamics (no momentum advection)
   LOGICAL, PUBLIC ::   ln_dynadv_vec   !: vector form
   INTEGER, PUBLIC ::      nn_dynkeg       !: scheme of grad(KE): =0 C2 ; =1 Hollingsworth
   LOGICAL, PUBLIC ::   ln_dynadv_cen2  !: flux form - 2nd order centered scheme flag
   LOGICAL, PUBLIC ::   ln_dynadv_ubs   !: flux form - 3rd order UBS scheme flag (OLD)
   LOGICAL, PUBLIC ::   ln_dynadv_up3   !: flux form - 3rd order UBS scheme flag (NEW)

   INTEGER, PUBLIC ::   n_dynadv   !: choice of the formulation and scheme for momentum advection
   !                               !  associated indices:
   INTEGER, PUBLIC, PARAMETER ::   np_LIN_dyn = 0   ! no advection: linear dynamics
   INTEGER, PUBLIC, PARAMETER ::   np_VEC_c2  = 1   ! vector form : 2nd order centered scheme
   INTEGER, PUBLIC, PARAMETER ::   np_FLX_c2  = 2   ! flux   form : 2nd order centered scheme
   INTEGER, PUBLIC, PARAMETER ::   np_FLX_ubs = 3   ! flux   form : 3rd order Upstream Biased Scheme (OLD)
   INTEGER, PUBLIC, PARAMETER ::   np_FLX_up3 = 4   ! flux   form : 3rd order Upstream Biased Scheme (NEW)
   
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynadv.F90 14053 2020-12-03 13:48:38Z techene $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_adv( kt, Kbb, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv  ***
      !!                
      !! ** Purpose :   compute the ocean momentum advection trend.
      !!
      !! ** Method  : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the advection term following n_dynadv
      !!
      !!      NB: in flux form advection (ln_dynadv_cen2 or ln_dynadv_ubs=T) 
      !!      a metric term is add to the coriolis term while in vector form 
      !!      it is the relative vorticity which is added to coriolis term
      !!      (see dynvor module).
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt               ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kbb, Kmm, Krhs   ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv         ! ocean velocities and RHS of momentum equation
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start( 'dyn_adv' )
      !
      SELECT CASE( n_dynadv )    !==  compute advection trend and add it to general trend  ==!
      CASE( np_VEC_c2  )     
         CALL dyn_keg     ( kt, nn_dynkeg,      Kmm, puu, pvv, Krhs )    ! vector form : horizontal gradient of kinetic energy
         CALL dyn_zad     ( kt,                 Kmm, puu, pvv, Krhs )    ! vector form : vertical advection
      CASE( np_FLX_c2  ) 
         CALL dyn_adv_cen2( kt,                 Kmm, puu, pvv, Krhs )    ! 2nd order centered scheme
      CASE( np_FLX_ubs )   
         CALL dyn_adv_ubs ( kt,            Kbb, Kmm, puu, pvv, Krhs )    ! 3rd order OLD UBS  scheme (UP3)
      CASE( np_FLX_up3 )
         CALL dyn_adv_up3 ( kt,            Kbb, Kmm, puu, pvv, Krhs )    ! 3rd order NEW UBS  scheme (UP3)
      END SELECT
      !
      IF( ln_timing )   CALL timing_stop( 'dyn_adv' )
      !
   END SUBROUTINE dyn_adv

   SUBROUTINE dyn_adv_lll( kt, Kbb, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv  ***
      !!                
      !! ** Purpose :   compute the ocean momentum advection trend.
      !!
      !! ** Method  : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the advection term following n_dynadv
      !!
      !!      NB: in flux form advection (ln_dynadv_cen2 or ln_dynadv_ubs=T) 
      !!      a metric term is add to the coriolis term while in vector form 
      !!      it is the relative vorticity which is added to coriolis term
      !!      (see dynvor module).
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt               ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kbb, Kmm, Krhs   ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv         ! ocean velocities and RHS of momentum equation
      !!----------------------------------------------------------------------
      
      TD_INIT
      
      IF( ln_timing ) then 
         MThRun 
           CALL timing_start( 'dyn_adv' )
         EndMThRun
      end if
      
      SELECT CASE( n_dynadv )    !==  compute advection trend and add it to general trend  ==!
      CASE( np_VEC_c2  )     
         CALL dyn_keg_lll     ( kt, nn_dynkeg,      Kmm, puu, pvv, Krhs )    ! vector form : horizontal gradient of kinetic energy
     
         CALL dyn_zad_lll     ( kt,                 Kmm, puu, pvv, Krhs )    ! vector form : vertical advection
        
      CASE( np_FLX_c2  ) 
         MThRun
            CALL dyn_adv_cen2( kt,                 Kmm, puu, pvv, Krhs )    ! 2nd order centered scheme
         EndMThRun
      CASE( np_FLX_ubs )   
         MThRun
            CALL dyn_adv_ubs ( kt,            Kbb, Kmm, puu, pvv, Krhs )    ! 3rd order OLD UBS  scheme (UP3)
         EndMThRun
      CASE( np_FLX_up3 )
         MThRun
            CALL dyn_adv_up3 ( kt,            Kbb, Kmm, puu, pvv, Krhs )    ! 3rd order NEW UBS  scheme (UP3)
         EndMThRun
      END SELECT
      MThRun
         IF( ln_timing )   CALL timing_stop( 'dyn_adv' )
      EndMThRun
   END SUBROUTINE dyn_adv_lll

   
   
   SUBROUTINE dyn_adv_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_init  ***
      !!                
      !! ** Purpose :   Control the consistency between namelist options for 
      !!              momentum advection formulation & scheme and set n_dynadv
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ios   ! Local integer
      !
      NAMELIST/namdyn_adv/ ln_dynadv_OFF, ln_dynadv_vec, nn_dynkeg, ln_dynadv_cen2, ln_dynadv_ubs, ln_dynadv_up3
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_adv_init : choice/control of the momentum advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      READ  ( numnam_ref(INDEX(numnam_ref,"namdyn_adv")-1:), namdyn_adv, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn_adv in reference namelist' )
      READ  ( numnam_cfg(INDEX(numnam_cfg,"namdyn_adv")-1:), namdyn_adv, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdyn_adv in configuration namelist' )
      IF(lwm) WRITE ( numond, namdyn_adv )

      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*) '   Namelist namdyn_adv : chose a advection formulation & scheme for momentum'
         WRITE(numout,*) '      linear dynamics : no momentum advection          ln_dynadv_OFF  = ', ln_dynadv_OFF
         WRITE(numout,*) '      Vector form: 2nd order centered scheme           ln_dynadv_vec  = ', ln_dynadv_vec
         WRITE(numout,*) '         with Hollingsworth scheme (=1) or not (=0)       nn_dynkeg   = ', nn_dynkeg
         WRITE(numout,*) '      flux form: 2nd order centred scheme              ln_dynadv_cen2 = ', ln_dynadv_cen2
         WRITE(numout,*) '                 3rd order UBS scheme (OLD)            ln_dynadv_ubs  = ', ln_dynadv_ubs
         WRITE(numout,*) '                 3rd order UBS scheme (NEW)            ln_dynadv_up3  = ', ln_dynadv_up3
      ENDIF

      ioptio = 0                      ! parameter control and set n_dynadv
      IF( ln_dynadv_OFF  ) THEN   ;   ioptio = ioptio + 1   ;   n_dynadv = np_LIN_dyn   ;   ENDIF
      IF( ln_dynadv_vec  ) THEN   ;   ioptio = ioptio + 1   ;   n_dynadv = np_VEC_c2    ;   ENDIF
      IF( ln_dynadv_cen2 ) THEN   ;   ioptio = ioptio + 1   ;   n_dynadv = np_FLX_c2    ;   ENDIF
      IF( ln_dynadv_ubs  ) THEN   ;   ioptio = ioptio + 1   ;   n_dynadv = np_FLX_ubs   ;   ENDIF
      IF( ln_dynadv_up3  ) THEN   ;   ioptio = ioptio + 1   ;   n_dynadv = np_FLX_up3   ;   ENDIF

      IF( ioptio /= 1 )   CALL ctl_stop( 'choose ONE and only ONE advection scheme' )
      IF( nn_dynkeg /= nkeg_C2 .AND. nn_dynkeg /= nkeg_HW )   CALL ctl_stop( 'KEG scheme wrong value of nn_dynkeg' )
#if defined key_qcoTest_FluxForm
      IF( ln_dynadv_vec  ) THEN CALL ctl_stop( 'STOP', 'key_qcoTest_FluxForm requires flux form advection' )
#endif

      IF(lwp) THEN                    ! Print the choice
         WRITE(numout,*)
         SELECT CASE( n_dynadv )
         CASE( np_LIN_dyn )   ;   WRITE(numout,*) '   ==>>>   linear dynamics : no momentum advection used'
         CASE( np_VEC_c2  )   ;   WRITE(numout,*) '   ==>>>   vector form : keg + zad + vor is used' 
            IF( nn_dynkeg == nkeg_C2  )   WRITE(numout,*) '              with Centered standard keg scheme'
            IF( nn_dynkeg == nkeg_HW  )   WRITE(numout,*) '              with Hollingsworth keg scheme'
         CASE( np_FLX_c2  )   ;   WRITE(numout,*) '   ==>>>   flux form   : 2nd order scheme is used'
         CASE( np_FLX_ubs )   ;   WRITE(numout,*) '   ==>>>   flux form   : OLD   UBS scheme is used'
         CASE( np_FLX_up3 )   ;   WRITE(numout,*) '   ==>>>   flux form   : NEW   UBS scheme is used'
         END SELECT
      ENDIF
      !
   END SUBROUTINE dyn_adv_init

  !!======================================================================
END MODULE dynadv
