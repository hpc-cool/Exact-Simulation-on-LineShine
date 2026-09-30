MODULE dynadv_up3
   !!======================================================================
   !!                       ***  MODULE  dynadv_up3  ***
   !! Ocean dynamics: Update the momentum trend with the flux form advection
   !!                 trend using a 3rd order upstream biased scheme
   !!======================================================================
   !! History :  2.0  ! 2006-08  (R. Benshila, L. Debreu)  Original code
   !!            3.2  ! 2009-07  (R. Benshila)  Suppression of rigid-lid option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_adv_up3   : flux form momentum advection using    (ln_dynadv=T)
   !!                   an 3rd order Upstream Biased Scheme or Quick scheme
   !!                   combined with 2nd or 4th order finite differences 
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   REAL(wp), PARAMETER :: gamma1 = 1._wp/3._wp  ! =1/4 quick      ; =1/3  3rd order UBS

   PUBLIC   dyn_adv_up3   ! routine called by step.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynadv_up3.F90 14834 2021-05-11 09:24:44Z hadcv $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_adv_up3( kt, Kbb, Kmm, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_up3  ***
      !!
      !! ** Purpose :   Compute the now momentum advection trend in flux form
      !!              and the general trend of the momentum equation.
      !!
      !! ** Method  :   The scheme is the one implemeted in ROMS. It depends 
      !!      on two parameter gamma1 and gamma2. The former control the 
      !!      upstream baised part of the scheme and the later the centred 
      !!      part:     gamma1 = 0    pure centered  (no diffusive part)
      !!                       = 1/4  Quick scheme
      !!                       = 1/3  3rd order Upstream biased scheme
      !!      For stability reasons, the first term of the fluxes which cor-
      !!      responds to a second order centered scheme is evaluated using  
      !!      the now velocity (centered in time) while the second term which  
      !!      is the diffusive part of the scheme, is evaluated using the 
      !!      before velocity (forward in time). 
      !!      Default value (hard coded in the begining of the module) are 
      !!      gamma1=1/3 and gamma2=1/32.
      !!
      !! ** Action : - (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) updated with the 3D advective momentum trends
      !!
      !! Reference : Shchepetkin & McWilliams, 2005, Ocean Modelling. 
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt              ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kbb, Kmm, Krhs  ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv        ! ocean velocities and RHS of momentum equation
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zui, zvj, zfuj, zfvi, zl_u, zl_v   ! local scalars
      REAL(wp) ::   zfwi, zzfu_uw_kp1, zlu_uw_kp1   ! local scalars
      REAL(wp) ::   zfwj, zzfv_vw_kp1, zlv_vw_kp1   ! local scalars
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) :: zfu_f, zfu
      REAL(dp), DIMENSION(A2D(nn_hls),jpk) :: zfu_t
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) :: zfv_f, zfv
      REAL(dp), DIMENSION(A2D(nn_hls),jpk) :: zfv_t
      REAL(dp), DIMENSION(A2D(nn_hls),jpk) :: ztrdu, ztrdv
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) :: zlu_uu, zlu_uv
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) :: zlv_vv, zlv_vu
      REAL(dp), DIMENSION(A2D(nn_hls)) :: zfw, zfu_uw, zfv_vw, zlu_uw, zlv_vw
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn_adv_ubs : UBS flux form momentum advection'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         ENDIF
      ENDIF
      !
      zfu_t(:,:,:) = 0._wp
      zfv_t(:,:,:) = 0._wp
      zfu_f(:,:,:) = 0._wp
      zfv_f(:,:,:) = 0._wp
      !
      zlu_uu(:,:,:) = 0._wp
      zlv_vv(:,:,:) = 0._wp 
      zlu_uv(:,:,:) = 0._wp 
      zlv_vu(:,:,:) = 0._wp
      !
      IF( l_trddyn ) THEN           ! trends: store the input trends
         ztrdu(:,:,:) = puu(:,:,:,Krhs)
         ztrdv(:,:,:) = pvv(:,:,:,Krhs)
      ENDIF
      !                                      ! =========================== !
      DO jk = 1, jpkm1                       !  Laplacian of the velocity  !
         !                                   ! =========================== !
         !                                         ! horizontal volume fluxes
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            zfu(ji,jj,jk) = e2u(ji,jj) * e3u(ji,jj,jk,Kmm) * puu(ji,jj,jk,Kmm)
            zfv(ji,jj,jk) = e1v(ji,jj) * e3v(ji,jj,jk,Kmm) * pvv(ji,jj,jk,Kmm)
         END_2D
         !            
         DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )                       ! laplacian
            ! round brackets added to fix the order of floating point operations
            ! needed to ensure halo 1 - halo 2 compatibility
            zlu_uu(ji,jj,jk) = ( ( puu (ji+1,jj  ,jk,Kbb) - puu (ji  ,jj  ,jk,Kbb) &
               &                 )                                                 & ! bracket for halo 1 - halo 2 compatibility
               &               + ( puu (ji-1,jj  ,jk,Kbb) - puu (ji  ,jj  ,jk,Kbb) &
               &                 )                                                 & ! bracket for halo 1 - halo 2 compatibility
               &               ) * umask(ji  ,jj  ,jk)
            zlv_vv(ji,jj,jk) = ( ( pvv (ji  ,jj+1,jk,Kbb) - pvv (ji  ,jj  ,jk,Kbb) &
               &                 )                                                 & ! bracket for halo 1 - halo 2 compatibility
               &               + ( pvv (ji  ,jj-1,jk,Kbb) - pvv (ji  ,jj  ,jk,Kbb) &
               &                 )                                                 & ! bracket for halo 1 - halo 2 compatibility
               &               ) * vmask(ji  ,jj  ,jk)
            zlu_uv(ji,jj,jk) = (  puu (ji  ,jj+1,jk,Kbb) - puu (ji  ,jj  ,jk,Kbb)  ) * fmask(ji  ,jj  ,jk)   &
               &             - (  puu (ji  ,jj  ,jk,Kbb) - puu (ji  ,jj-1,jk,Kbb)  ) * fmask(ji  ,jj-1,jk)
            zlv_vu(ji,jj,jk) = (  pvv (ji+1,jj  ,jk,Kbb) - pvv (ji  ,jj  ,jk,Kbb)  ) * fmask(ji  ,jj  ,jk)   &
               &             - (  pvv (ji  ,jj  ,jk,Kbb) - pvv (ji-1,jj  ,jk,Kbb)  ) * fmask(ji-1,jj  ,jk)
         END_2D
      END DO
      !
      IF( nn_hls == 1 ) CALL lbc_lnk( 'dynadv_ubs', zlu_uu(:,:,:), 'U', -1.0_wp , zlu_uv(:,:,:), 'U', -1.0_wp,  &
                                                &   zlv_vv(:,:,:), 'V', -1.0_wp , zlv_vu(:,:,:), 'V', -1.0_wp   )
      !
      !                                      ! ====================== !
      !                                      !  Horizontal advection  !
      DO jk = 1, jpkm1                       ! ====================== !
         !                                         ! horizontal volume fluxes
         DO_2D( 1, 1, 1, 1 )
            zfu(ji,jj,jk) = 0.25_wp * e2u(ji,jj) * e3u(ji,jj,jk,Kmm) * puu(ji,jj,jk,Kmm)
            zfv(ji,jj,jk) = 0.25_wp * e1v(ji,jj) * e3v(ji,jj,jk,Kmm) * pvv(ji,jj,jk,Kmm)
         END_2D
         !
         DO_2D( 1, 0, 1, 0 )                       ! horizontal momentum fluxes at T- and F-point
            zui = ( puu(ji,jj,jk,Kmm) + puu(ji+1,jj  ,jk,Kmm) )
            zvj = ( pvv(ji,jj,jk,Kmm) + pvv(ji  ,jj+1,jk,Kmm) )
            !
            IF( zui > 0 ) THEN   ;   zl_u = zlu_uu(ji  ,jj,jk)
            ELSE                 ;   zl_u = zlu_uu(ji+1,jj,jk)
            ENDIF
            IF( zvj > 0 ) THEN   ;   zl_v = zlv_vv(ji,jj  ,jk)
            ELSE                 ;   zl_v = zlv_vv(ji,jj+1,jk)
            ENDIF
            !
            zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj  ,jk) )                              &
               &                * ( zui - gamma1 * zl_u)
            zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji  ,jj+1,jk) )                             &
               &                * ( zvj - gamma1 * zl_v)
            !
            zfuj = ( zfu(ji,jj,jk) + zfu(ji  ,jj+1,jk) )
            zfvi = ( zfv(ji,jj,jk) + zfv(ji+1,jj  ,jk) )
            IF( zfuj > 0 ) THEN   ;    zl_v = zlv_vu( ji  ,jj  ,jk )
            ELSE                  ;    zl_v = zlv_vu( ji+1,jj  ,jk )
            ENDIF
            IF( zfvi > 0 ) THEN   ;    zl_u = zlu_uv( ji,jj  ,jk )
            ELSE                  ;    zl_u = zlu_uv( ji,jj+1,jk )
            ENDIF
            !
            zfv_f(ji  ,jj  ,jk) = ( zfvi )   &
               &                * ( puu(ji,jj,jk,Kmm) + puu(ji  ,jj+1,jk,Kmm) - gamma1 * zl_u )
            zfu_f(ji  ,jj  ,jk) = ( zfuj )   &
               &                * ( pvv(ji,jj,jk,Kmm) + pvv(ji+1,jj  ,jk,Kmm) - gamma1 * zl_v )
         END_2D
         DO_2D( 0, 0, 0, 0 )                       ! divergence of horizontal momentum fluxes
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)    &
               &                           + zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  ) * r1_e1e2u(ji,jj)   &
               &                           / e3u(ji,jj,jk,Kmm)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)    &
               &                           + zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  ) * r1_e1e2v(ji,jj)   &
               &                           / e3v(ji,jj,jk,Kmm)
         END_2D
      END DO
      IF( l_trddyn ) THEN                          ! trends: send trends to trddyn for diagnostic
         ztrdu(:,:,:) = puu(:,:,:,Krhs) - ztrdu(:,:,:)
         ztrdv(:,:,:) = pvv(:,:,:,Krhs) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_keg, kt, Kmm )
         zfu_t(:,:,:) = puu(:,:,:,Krhs)
         zfv_t(:,:,:) = pvv(:,:,:,Krhs)
      ENDIF
      !                                ! ==================== !
      !                                !  Vertical advection  !
      !                                ! ==================== !
      !
      !                                   != surface vertical fluxes =! (jk = 1)
      !
      DO_2D( 0, 1, 0, 1 )
         zfw(ji,jj) = e1e2t(ji,jj) * ww(ji,jj,1)
      END_2D
      !
      IF( ln_linssh ) THEN                       ! linear free surface: advection through the surface z=0
         DO_2D( 0, 0, 0, 0 )
            zfu_uw(ji,jj) = 0.5_wp * ( zfw(ji,jj) + zfw(ji+1,jj  ) ) * puu(ji,jj,1,Kmm)
            zfv_vw(ji,jj) = 0.5_wp * ( zfw(ji,jj) + zfw(ji  ,jj+1) ) * pvv(ji,jj,1,Kmm)
         END_2D
      ELSE                                       ! non linear free: surface advective fluxes set to zero
         DO_2D( 0, 0, 0, 0 )
            zfu_uw(ji,jj) = 0._wp
            zfv_vw(ji,jj) = 0._wp
         END_2D
      ENDIF
      !
      DO_2D( 0, 0, 0, 0 )                        ! uniform shear in the 1st layer : dz(u(k=1)) = dz(u(k=2)) ==> zlu = 0
         zlu_uw(ji,jj) = 0._wp
         zlv_vw(ji,jj) = 0._wp
      END_2D
      !
      DO jk = 1, jpk-2                    != divergence of advective fluxes =! (jk = 1 to jpk-2)
         DO_2D( 0, 1, 0, 1 )                  ! vertical transport at level k+1
            zfw(ji,jj) = e1e2t(ji,jj) * ww(ji,jj,jk+1)
         END_2D
         !
         DO_2D( 0, 0, 0, 0 )
            zlu_uw_kp1 = ( puu(ji,jj,jk  ,Kbb) - puu(ji,jj,jk+1,Kbb) ) * wumask(ji,jj,jk+1)   &
               &       - ( puu(ji,jj,jk+1,Kbb) - puu(ji,jj,jk+2,Kbb) ) * wumask(ji,jj,jk+2)
            zlv_vw_kp1 = ( pvv(ji,jj,jk  ,Kbb) - pvv(ji,jj,jk+1,Kbb) ) * wvmask(ji,jj,jk+1)   &
               &       - ( pvv(ji,jj,jk+1,Kbb) - pvv(ji,jj,jk+2,Kbb) ) * wvmask(ji,jj,jk+2)
            !
            zfwi = zfw(ji,jj) + zfw(ji+1,jj)
            IF( zfwi > 0 ) THEN   ;   zl_u = zlu_uw_kp1
            ELSE                  ;   zl_u = zlu_uw(ji,jj)
            ENDIF
            zfwj = zfw(ji,jj) + zfw(ji,jj+1)
            IF( zfwj > 0 ) THEN   ;   zl_v = zlv_vw_kp1
            ELSE                  ;   zl_v = zlv_vw(ji,jj)
            ENDIF
            !                                    ! vertical flux at level k+1
            zzfu_uw_kp1 = 0.25_wp * ( zfw(ji,jj) + zfw(ji+1,jj  ) ) * ( puu(ji,jj,jk+1,Kmm) + puu(ji,jj,jk,Kmm) - gamma1 * zl_u )
            zzfv_vw_kp1 = 0.25_wp * ( zfw(ji,jj) + zfw(ji  ,jj+1) ) * ( pvv(ji,jj,jk+1,Kmm) + pvv(ji,jj,jk,Kmm) - gamma1 * zl_v )
            !                                    ! divergence of vertical momentum flux
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - ( zfu_uw(ji,jj) - zzfu_uw_kp1 ) * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - ( zfv_vw(ji,jj) - zzfv_vw_kp1 ) * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)
            !                                    ! store vertical flux for next level calculation
            zfu_uw(ji,jj) = zzfu_uw_kp1
            zfv_vw(ji,jj) = zzfv_vw_kp1
            !
            zlu_uw(ji,jj) = zlu_uw_kp1
            zlv_vw(ji,jj) = zlv_vw_kp1
         END_2D
      END DO
      !
      jk = jpkm1                          != compute last level =! (zzfu_uw_kp1 = zzfu_uw_kp1 = 0)
      DO_2D( 0, 0, 0, 0 )
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - zfu_uw(ji,jj) * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - zfv_vw(ji,jj) * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)
      END_2D
      !
      IF( l_trddyn ) THEN                       ! save the vertical advection trend for diagnostic
         zfu_t(:,:,:) = puu(:,:,:,Krhs) - zfu_t(:,:,:)
         zfv_t(:,:,:) = pvv(:,:,:,Krhs) - zfv_t(:,:,:)
         CALL trd_dyn( zfu_t, zfv_t, jpdyn_zad, kt, Kmm )
      ENDIF
      !                                         ! Control print
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' ubs2 adv - Ua: ', mask1=umask,   &
         &                                  tab3d_2=pvv(:,:,:,Krhs), clinfo2=           ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
   END SUBROUTINE dyn_adv_up3

   !!==============================================================================
END MODULE dynadv_up3
