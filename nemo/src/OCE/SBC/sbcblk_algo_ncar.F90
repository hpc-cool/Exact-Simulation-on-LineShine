MODULE sbcblk_algo_ncar
   !!======================================================================
   !!                   ***  MODULE  sbcblk_algo_ncar  ***
   !! Computes:
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m Ubzu
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of Large & Yeager 2008
   !!
   !!       Routine turb_ncar maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!                         L. Brodeau, 2015
   !!=====================================================================
   !! History :  3.6  !  2016-02  (L.Brodeau) successor of old turb_ncar of former sbcblk_core.F90
   !!            4.2  !  2020-12  (L. Brodeau) Introduction of various air-ice bulk parameterizations + improvements
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   turb_ncar  : computes the bulk turbulent transfer coefficients
   !!                   adjusts t_air and q_air from zt to zu m
   !!                   returns the effective bulk wind speed at 10m
   !!----------------------------------------------------------------------
   USE oce 
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce, ONLY: ln_cdgw
   USE sbcwave, ONLY: cdn_wave ! wave module
   USE phycst          ! physical constants
   USE sbc_phy         ! Catalog of functions for physical/meteorological parameters in the marine boundary layer

   IMPLICIT NONE
   PRIVATE
   REAL(wp), DIMENSION(:,:), allocatable ::   zCdN, zCeN, zChN        ! 10m neutral latent/sensible coefficient
   REAL(wp), DIMENSION(:,:), allocatable ::   zsqrt_Cd, zsqrt_CdN   ! root square of Cd and Cd_neutral
   REAL(wp), DIMENSION(:,:), allocatable ::   zeta_u        ! stability parameter at height zu
   REAL(wp), DIMENSION(:,:), allocatable ::   ztmp0, ztmp1, ztmp2
   PUBLIC :: TURB_NCAR, turb_ncar_lll   ! called by sbcblk.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"

   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE turb_ncar_lll(kt,  zt, zu, sst, t_zt, ssq, q_zt, U_zu, &
      &                     Cd, Ch, Ce, t_zu, q_zu, Ubzu,       &
      &                     nb_iter, CdN, ChN, CeN               )
      INTEGER , INTENT(in   )                     ::   kt
      REAL(wp), INTENT(in   )                     ::    zt       ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                             [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature                [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   ssq      ! sea surface specific humidity           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity at zt             [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ubzu    ! bulk wind speed at zu                     [m/s]
      !
      INTEGER , INTENT(in   ), OPTIONAL                     :: nb_iter  ! number of iterations
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   CdN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   ChN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   CeN
      !
      INTEGER :: nbit, jit,ji,jj                   ! iterations...
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      TD_INIT
      !
      !!----------------------------------------------------------------------------------
      nbit = nb_iter0
      
      if(kt == 1) then 
         MThRun
            allocate(zCdN(jpi,jpj), zCeN(jpi,jpj), zChN(jpi,jpj), zsqrt_Cd(jpi,jpj), zsqrt_CdN(jpi,jpj), &
               &  zeta_u(jpi,jpj), ztmp0(jpi,jpj), ztmp1(jpi,jpj), ztmp2(jpi,jpj))
         EndMThRun
         SYNCA
      end if
      IF( PRESENT(nb_iter) ) nbit = nb_iter
      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp ) ! testing "zu == zt" is risky with double precision

      
      Ubzu(YJPI,YJPJ) = MAX( 0.5_wp , U_zu(YJPI,YJPJ) ) 
     
      YDO_2D( 1, 1, 1, 1)
          ztmp0(ji,jj) = virt_temp(t_zt(ji,jj), q_zt(ji,jj)) - virt_temp(sst(ji,jj), ssq(ji,jj))                 
      YEND_2D

      ztmp1(YJPI,YJPJ) = 0.5_wp + SIGN(0.5_wp, ztmp0(YJPI,YJPJ))   
      IF( ln_cdgw ) THEN      
        
         cdn_wave(YJPI,YJPJ) = cdn_wave(YJPI,YJPJ) + rsmall * ( 1._wp - tmask(ji,jj,1) )
         zCdN(YJPI,YJPJ)   = cdn_wave(YJPI,YJPJ)
        
      ELSE
         YDO_2D( 1, 1, 1, 1)
            zCdN(ji,jj) = cd_n10_ncar_lll( Ubzu(ji,jj) )
         YEND_2D
      ENDIF
 
      zsqrt_CdN(YJPI,YJPJ) = SQRT( zCdN(YJPI,YJPJ) )
      Cd(YJPI,YJPJ) = zCdN(YJPI,YJPJ)
      
      zsqrt_Cd(YJPI,YJPJ) = zsqrt_CdN(YJPI,YJPJ)
   
      YDO_2D( 1, 1, 1, 1)
         Ce(ji,jj) = ce_n10_ncar_lll( zsqrt_CdN(ji,jj) )
         Ch(ji,jj) = ch_n10_ncar_lll( zsqrt_CdN(ji,jj) , ztmp1(ji,jj) )   
      YEND_2D

      IF( ln_cdgw ) THEN
        
         zCeN(YJPI,YJPJ) = Ce(YJPI,YJPJ)
         zChN(YJPI,YJPJ) = Ch(YJPI,YJPJ)
        
      ENDIF
      
          t_zu(YJPI,YJPJ) = MAX( t_zt(YJPI,YJPJ) ,  180._wp )   
          q_zu(YJPI,YJPJ) = MAX( q_zt(YJPI,YJPJ) , 1.e-6_wp )   

      DO jit = 1, nbit
                 
         ztmp1(YJPI,YJPJ) = t_zu(YJPI,YJPJ) - sst(YJPI,YJPJ)   
         ztmp2(YJPI,YJPJ) = q_zu(YJPI,YJPJ) - ssq(YJPI,YJPJ)

         ztmp0(YJPI,YJPJ) = zsqrt_Cd(YJPI,YJPJ) * Ubzu(YJPI,YJPJ)    
         ztmp1(YJPI,YJPJ) = Ch(YJPI,YJPJ)/zsqrt_Cd(YJPI,YJPJ) * ztmp1(YJPI,YJPJ)   
         ztmp2(YJPI,YJPJ) = Ce(YJPI,YJPJ)/zsqrt_Cd(YJPI,YJPJ) * ztmp2(YJPI,YJPJ)   
         YDO_2D( 1, 1, 1, 1)
            ztmp0(ji,jj) = One_on_L_lll( t_zu(ji,jj), q_zu(ji,jj), ztmp0(ji,jj), ztmp1(ji,jj), ztmp2(ji,jj) )
         YEND_2D 
  
         zeta_u(YJPI,YJPJ) = zu * ztmp0(YJPI,YJPJ)
         zeta_u(YJPI,YJPJ) = SIGN( MIN(ABS(zeta_u(YJPI,YJPJ)), 10._wp), zeta_u(YJPI,YJPJ) )

    
         IF( .NOT. l_zt_equal_zu ) THEN
           
            ztmp0(YJPI,YJPJ) = zt * ztmp0(YJPI,YJPJ) 
            ztmp0(YJPI,YJPJ) = SIGN( MIN(ABS(ztmp0(YJPI,YJPJ)), 10._wp), ztmp0(YJPI,YJPJ) )  
           
            YDO_2D( 1, 1, 1, 1)
               ztmp0(ji,jj) = LOG(zt/zu) + psi_h_ncar_lll(zeta_u(ji,jj)) - psi_h_ncar_lll(ztmp0(ji,jj))   
            YEND_2D                 
           
            t_zu(YJPI,YJPJ) = t_zt(YJPI,YJPJ) - ztmp1(YJPI,YJPJ)/vkarmn * ztmp0(YJPI,YJPJ)  
            q_zu(YJPI,YJPJ) = q_zt(YJPI,YJPJ) - ztmp2(YJPI,YJPJ)/vkarmn * ztmp0(YJPI,YJPJ)    ! L&Y 2004 Eq. (9c)
            q_zu(YJPI,YJPJ) = MAX(0._wp, q_zu(YJPI,YJPJ))  ! 湿度非负保护
         END IF

         YDO_2D( 1, 1, 1, 1)
            ztmp2(ji,jj) = psi_m_ncar_lll(zeta_u(ji,jj))  ! 动量通量的稳定度函数
         YEND_2D           

         IF( ln_cdgw ) THEN      ! 表面波情况
            ! 更新Cd的平方根
            zsqrt_Cd(YJPI,YJPJ) = vkarmn / ( vkarmn / zsqrt_CdN(YJPI,YJPJ) - ztmp2(YJPI,YJPJ) )
            Cd(YJPI,YJPJ)   = zsqrt_Cd(YJPI,YJPJ) * zsqrt_Cd(YJPI,YJPJ)
            ! 计算Ch/Ce调整项：调用psi_h_ncar_lll
            YDO_2D( 1, 1, 1, 1)
               ztmp0(ji,jj) = (LOG(zu/10._wp) - psi_h_ncar_lll(zeta_u(ji,jj))) / vkarmn / zsqrt_CdN(ji,jj)
            YEND_2D 
            ztmp2(YJPI,YJPJ) = zsqrt_Cd(YJPI,YJPJ) / zsqrt_CdN(YJPI,YJPJ)
            ztmp1(YJPI,YJPJ) = 1._wp + zChN(YJPI,YJPJ) * ztmp0(YJPI,YJPJ)
            Ch(YJPI,YJPJ)    = zChN(YJPI,YJPJ) * ztmp2(YJPI,YJPJ) / ztmp1(YJPI,YJPJ)  ! L&Y 2004 eq. (10b)
            ztmp1(YJPI,YJPJ) = 1._wp + zCeN(YJPI,YJPJ) * ztmp0(YJPI,YJPJ)
            Ce(YJPI,YJPJ)    = zCeN(YJPI,YJPJ) * ztmp2(YJPI,YJPJ) / ztmp1(YJPI,YJPJ)  ! L&Y 2004 eq. (10c)

         ELSE  
            YDO_2D( 1, 1, 1, 1)
               ztmp0(ji,jj) = MAX( 0.25_wp , UN10_from_CD_lll(zu, Ubzu(ji,jj), Cd(ji,jj), ppsi=ztmp2(ji,jj)) )
               zCdN(ji,jj) = cd_n10_ncar_lll(ztmp0(ji,jj))
            YEND_2D 

            zsqrt_CdN(YJPI,YJPJ) = SQRT(zCdN(YJPI,YJPJ))

            !! 更新拖曳系数C_D (L&Y 2004 Eq. (10a))
            ztmp1(YJPI,YJPJ)  = 1._wp + zsqrt_CdN(YJPI,YJPJ)/vkarmn * (LOG(zu/10._wp) - ztmp2(YJPI,YJPJ))   
            Cd(YJPI,YJPJ)     = MAX( zCdN(YJPI,YJPJ) / ( ztmp1(YJPI,YJPJ)*ztmp1(YJPI,YJPJ) ), Cx_min )

            !! 更新C_H和C_E的中间变量：调用psi_h_ncar_lll
            zsqrt_Cd(YJPI,YJPJ) = SQRT( Cd(YJPI,YJPJ) )
           
            YDO_2D( 1, 1, 1, 1)
               ztmp0(ji,jj) = ( LOG(zu/10._wp) - psi_h_ncar_lll(zeta_u(ji,jj)) ) / vkarmn / zsqrt_CdN(ji,jj)
            YEND_2D 
            
            ztmp2(YJPI,YJPJ) = zsqrt_Cd(YJPI,YJPJ) / zsqrt_CdN(YJPI,YJPJ)

            ! 更新稳定度标识（1=稳定，0=不稳定）
            ztmp1(YJPI,YJPJ) = 0.5_wp + SIGN(0.5_wp, zeta_u(YJPI,YJPJ))                                
            ! 计算中性Ch (L&Y 2004 eq. (6c-6d))
            zChN(YJPI,YJPJ)  = 1.e-3_wp * zsqrt_CdN(YJPI,YJPJ) * (18._wp*ztmp1(YJPI,YJPJ) + 32.7_wp*(1._wp - ztmp1(YJPI,YJPJ)))  
            ! 计算中性Ce (L&Y 2004 eq. (6b))
            zCeN(YJPI,YJPJ)  = 1.e-3_wp * (34.6_wp * zsqrt_CdN(YJPI,YJPJ))                             

            ! 更新Ch和Ce (L&Y 2004 eq. (10b-10c))
            Ch(YJPI,YJPJ)    = MAX( zChN(YJPI,YJPJ)*ztmp2(YJPI,YJPJ) / ( 1._wp + zChN(YJPI,YJPJ)*ztmp0(YJPI,YJPJ) ) , Cx_min )
            Ce(YJPI,YJPJ)    = MAX( zCeN(YJPI,YJPJ)*ztmp2(YJPI,YJPJ) / ( 1._wp + zCeN(YJPI,YJPJ)*ztmp0(YJPI,YJPJ) ) , Cx_min )

         ENDIF
      END DO !DO jit = 1, nbit
      IF(PRESENT(CdN)) CdN(YJPI,YJPJ) = zCdN(YJPI,YJPJ)
      IF(PRESENT(CeN)) CeN(YJPI,YJPJ) = zCeN(YJPI,YJPJ)
      IF(PRESENT(ChN)) ChN(YJPI,YJPJ) = zChN(YJPI,YJPJ)

   END SUBROUTINE turb_ncar_lll

   SUBROUTINE turb_ncar( kt,   zt, zu, sst, t_zt, ssq, q_zt, U_zu, &
      &                     Cd, Ch, Ce, t_zu, q_zu, Ubzu,       &
      &                     nb_iter, CdN, ChN, CeN               )
      INTEGER , INTENT(in   )                     ::   kt
      REAL(wp), INTENT(in   )                     ::   zt       ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                             [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature                [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   ssq      ! sea surface specific humidity           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity at zt             [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ubzu    ! bulk wind speed at zu                     [m/s]
      !
      INTEGER , INTENT(in   ), OPTIONAL                     :: nb_iter  ! number of iterations
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   CdN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   ChN
      REAL(wp), INTENT(  out), OPTIONAL, DIMENSION(jpi,jpj) ::   CeN
      !
      INTEGER :: nbit, jit                    ! iterations...
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   zCdN, zCeN, zChN        ! 10m neutral latent/sensible coefficient
      REAL(wp), DIMENSION(jpi,jpj) ::   zsqrt_Cd, zsqrt_CdN   ! root square of Cd and Cd_neutral
      REAL(wp), DIMENSION(jpi,jpj) ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(jpi,jpj) ::   ztmp0, ztmp1, ztmp2
      !!----------------------------------------------------------------------------------
      nbit = nb_iter0
      IF( PRESENT(nb_iter) ) nbit = nb_iter

      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp ) ! testing "zu == zt" is risky with double precision

      Ubzu = MAX( 0.5_wp , U_zu )   !  relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s

      !! First guess of stability:
      ztmp0 = virt_temp(t_zt, q_zt) - virt_temp(sst, ssq) ! air-sea difference of virtual pot. temp. at zt
      ztmp1 = 0.5_wp + SIGN(0.5_wp,ztmp0)                 ! ztmp1 = 1 if dTv > 0  => STABLE, 0 if unstable

      !! Neutral coefficients at 10m:
      IF( ln_cdgw ) THEN      ! wave drag case
         cdn_wave(:,:) = cdn_wave(:,:) + rsmall * ( 1._wp - tmask(:,:,1) )
         zCdN   (:,:) = cdn_wave(:,:)
      ELSE
      zCdN = cd_n10_ncar( Ubzu )
      ENDIF

      zsqrt_CdN = SQRT( zCdN )

      !! Initializing transf. coeff. with their first guess neutral equivalents :
      Cd = zCdN
      Ce = ce_n10_ncar( zsqrt_CdN )
      Ch = ch_n10_ncar( zsqrt_CdN , ztmp1 )   ! ztmp1 is stability (1/0)
      zsqrt_Cd = zsqrt_CdN

      IF( ln_cdgw ) THEN
         zCeN = Ce
         zChN = Ch
      ENDIF

      !! Initializing values at z_u with z_t values:
      t_zu = MAX( t_zt ,  180._wp )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 1.e-6_wp )   !               "


      DO jit = 1, nbit
         !
         ztmp1 = t_zu - sst   ! Updating air/sea differences
         ztmp2 = q_zu - ssq

         ! Updating turbulent scales :   (L&Y 2004 Eq. (7))
         ztmp0 = zsqrt_Cd*Ubzu       ! u*
         ztmp1 = Ch/zsqrt_Cd*ztmp1    ! theta*
         ztmp2 = Ce/zsqrt_Cd*ztmp2    ! q*

         ! Estimate the inverse of Obukov length (1/L) at height zu:
         ztmp0 = One_on_L( t_zu, q_zu, ztmp0, ztmp1, ztmp2 )

         !! Stability parameters :
         zeta_u   = zu*ztmp0
         zeta_u   = sign( min(abs(zeta_u),10._wp), zeta_u )

         !! Shifting temperature and humidity at zu (L&Y 2004 Eq. (9b-9c))
         IF( .NOT. l_zt_equal_zu ) THEN
            ztmp0 = zt*ztmp0 ! zeta_t !
            ztmp0 = SIGN( MIN(ABS(ztmp0),10._wp), ztmp0 )  ! Temporaty array ztmp0 == zeta_t !!!
            ztmp0 = LOG(zt/zu) + psi_h_ncar(zeta_u) - psi_h_ncar(ztmp0)                   ! ztmp0 just used as temp array again!
            t_zu = t_zt - ztmp1/vkarmn*ztmp0    ! ztmp1 is still theta*  L&Y 2004 Eq. (9b)
            !!
            q_zu = q_zt - ztmp2/vkarmn*ztmp0    ! ztmp2 is still q*      L&Y 2004 Eq. (9c)
            q_zu = MAX(0._wp, q_zu)
         END IF

         ! Update neutral wind speed at 10m and neutral Cd at 10m (L&Y 2004 Eq. 9a)...
         !   In very rare low-wind conditions, the old way of estimating the
         !   neutral wind speed at 10m leads to a negative value that causes the code
         !   to crash. To prevent this a threshold of 0.25m/s is imposed.
         ztmp2 = psi_m_ncar(zeta_u)
         IF( ln_cdgw ) THEN      ! surface wave case
            zsqrt_Cd = vkarmn / ( vkarmn / zsqrt_CdN - ztmp2 )
            Cd   = zsqrt_Cd * zsqrt_Cd
            ztmp0 = (LOG(zu/10._wp) - psi_h_ncar(zeta_u)) / vkarmn / zsqrt_CdN
            ztmp2 = zsqrt_Cd / zsqrt_CdN
            ztmp1 = 1._wp + zChN * ztmp0
            Ch    = zChN * ztmp2 / ztmp1  ! L&Y 2004 eq. (10b)
            ztmp1 = 1._wp + zCeN * ztmp0
            Ce    = zCeN * ztmp2 / ztmp1  ! L&Y 2004 eq. (10c)

         ELSE
         ztmp0 = MAX( 0.25_wp , UN10_from_CD(zu, Ubzu, Cd, ppsi=ztmp2) ) ! U_n10 (ztmp2 == psi_m_ncar(zeta_u))

         zCdN = cd_n10_ncar(ztmp0)
         zsqrt_CdN = sqrt(zCdN)

         !! Update of transfer coefficients:

         !! C_D
         ztmp1  = 1._wp + zsqrt_CdN/vkarmn*(LOG(zu/10._wp) - ztmp2)   ! L&Y 2004 Eq. (10a) (ztmp2 == psi_m(zeta_u))
         Cd     = MAX( zCdN / ( ztmp1*ztmp1 ), Cx_min )

         !! C_H and C_E
         zsqrt_Cd = SQRT( Cd )
         ztmp0 = ( LOG(zu/10._wp) - psi_h_ncar(zeta_u) ) / vkarmn / zsqrt_CdN
         ztmp2 = zsqrt_Cd / zsqrt_CdN

         ztmp1 = 0.5_wp + SIGN(0.5_wp,zeta_u)                                ! update stability
         zChN  = 1.e-3_wp * zsqrt_CdN*(18._wp*ztmp1 + 32.7_wp*(1._wp - ztmp1))  ! L&Y 2004 eq. (6c-6d)
         zCeN  = 1.e-3_wp * (34.6_wp * zsqrt_CdN)                             ! L&Y 2004 eq. (6b)

         Ch    = MAX( zChN*ztmp2 / ( 1._wp + zChN*ztmp0 ) , Cx_min ) ! L&Y 2004 eq. (10b)
         Ce    = MAX( zCeN*ztmp2 / ( 1._wp + zCeN*ztmp0 ) , Cx_min ) ! L&Y 2004 eq. (10c)

         ENDIF

      END DO !DO jit = 1, nbit

      IF(PRESENT(CdN)) CdN(:,:) = zCdN(:,:)
      IF(PRESENT(CeN)) CeN(:,:) = zCeN(:,:)
      IF(PRESENT(ChN)) ChN(:,:) = zChN(:,:)

   END SUBROUTINE turb_ncar


   FUNCTION cd_n10_ncar( pw10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral drag coefficient at 10m as a function
      !! of neutral wind  speed at 10m
      !!
      !! Origin: Large & Yeager 2008, Eq. (11)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pw10           ! scalar wind speed at 10m (m/s)
      REAL(wp), DIMENSION(jpi,jpj)             :: cd_n10_ncar
      !
      INTEGER  ::     ji, jj     ! dummy loop indices
      REAL(wp) :: zgt33, zw, zw6 ! local scalars
      !!----------------------------------------------------------------------------------
      !
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            !
            zw  = pw10(ji,jj)
            zw6 = zw*zw*zw
            zw6 = zw6*zw6
            !
            ! When wind speed > 33 m/s => Cyclone conditions => special treatment
            zgt33 = 0.5_wp + SIGN( 0.5_wp, (zw - 33._wp) )   ! If pw10 < 33. => 0, else => 1
            !
            cd_n10_ncar(ji,jj) = 1.e-3_wp * ( &
               &       (1._wp - zgt33)*( 2.7_wp/zw + 0.142_wp + zw/13.09_wp - 3.14807E-10_wp*zw6) & ! wind <  33 m/s
               &      +    zgt33   *      2.34_wp )                                                 ! wind >= 33 m/s
            !
            cd_n10_ncar(ji,jj) = MAX( cd_n10_ncar(ji,jj), Cx_min )
            !
      END_2D
      !
   END FUNCTION cd_n10_ncar


   FUNCTION ch_n10_ncar( psqrtcdn10 , pstab )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral heat transfer coefficient at 10m      !!
      !! Origin: Large & Yeager 2008, Eq. (9) and (12)

      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: ch_n10_ncar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: psqrtcdn10 ! sqrt( CdN10 )
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pstab      ! stable ABL => 1 / unstable ABL => 0
      !!----------------------------------------------------------------------------------
      IF( ANY(pstab < -0.00001) .OR. ANY(pstab >  1.00001) ) THEN
         PRINT *, 'ERROR: ch_n10_ncar@mod_blk_ncar.f90: pstab ='
         PRINT *, pstab
         STOP
      END IF
      !
      ch_n10_ncar = MAX( 1.e-3_wp * psqrtcdn10*( 18._wp*pstab + 32.7_wp*(1._wp - pstab) )  , Cx_min )   ! Eq. (9) & (12) Large & Yeager, 2008
      !
   END FUNCTION ch_n10_ncar

   FUNCTION ce_n10_ncar( psqrtcdn10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral heat transfer coefficient at 10m      !!
      !! Origin: Large & Yeager 2008, Eq. (9) and (13)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: ce_n10_ncar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: psqrtcdn10 ! sqrt( CdN10 )
      !!----------------------------------------------------------------------------------
      ce_n10_ncar = MAX( 1.e-3_wp * ( 34.6_wp * psqrtcdn10 ) , Cx_min )
      !
   END FUNCTION ce_n10_ncar

   FUNCTION cd_n10_ncar_lll( pw10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral drag coefficient at 10m as a function
      !! of neutral wind  speed at 10m
      !!
      !! Origin: Large & Yeager 2008, Eq. (11)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pw10           ! scalar wind speed at 10m (m/s)
      REAL(wp)             :: cd_n10_ncar_lll
      !
      INTEGER  ::     ji, jj     ! dummy loop indices
      REAL(wp) :: zgt33, zw, zw6 ! local scalars
      !!----------------------------------------------------------------------------------
      !
     
            !
      zw  = pw10
      zw6 = zw*zw*zw
      zw6 = zw6*zw6
      !
      ! When wind speed > 33 m/s => Cyclone conditions => special treatment
      zgt33 = 0.5_wp + SIGN( 0.5_wp, (zw - 33._wp) )   ! If pw10 < 33. => 0, else => 1
      !
      cd_n10_ncar_lll = 1.e-3_wp * ( &
         &       (1._wp - zgt33)*( 2.7_wp/zw + 0.142_wp + zw/13.09_wp - 3.14807E-10_wp*zw6) & ! wind <  33 m/s
         &      +    zgt33   *      2.34_wp )                                                 ! wind >= 33 m/s
      !
      cd_n10_ncar_lll = MAX( cd_n10_ncar_lll, Cx_min )

   END FUNCTION cd_n10_ncar_lll


   FUNCTION ch_n10_ncar_lll( psqrtcdn10 , pstab )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral heat transfer coefficient at 10m      !!
      !! Origin: Large & Yeager 2008, Eq. (9) and (12)

      !!----------------------------------------------------------------------------------
      REAL(wp)             :: ch_n10_ncar_lll
      REAL(wp), INTENT(in) :: psqrtcdn10 ! sqrt( CdN10 )
      REAL(wp), INTENT(in) :: pstab      ! stable ABL => 1 / unstable ABL => 0
      !!----------------------------------------------------------------------------------
      IF( pstab < -0.00001 .OR. pstab >  1.00001 ) THEN
         PRINT *, 'ERROR: ch_n10_ncar@mod_blk_ncar.f90: pstab ='
         PRINT *, pstab
         STOP
      END IF
      !
      ch_n10_ncar_lll = MAX( 1.e-3_wp * psqrtcdn10*( 18._wp*pstab + 32.7_wp*(1._wp - pstab) )  , Cx_min )   ! Eq. (9) & (12) Large & Yeager, 2008
      !
   END FUNCTION ch_n10_ncar_lll

   FUNCTION ce_n10_ncar_lll( psqrtcdn10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral heat transfer coefficient at 10m      !!
      !! Origin: Large & Yeager 2008, Eq. (9) and (13)
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: ce_n10_ncar_lll
      REAL(wp), INTENT(in) :: psqrtcdn10 ! sqrt( CdN10 )
      !!----------------------------------------------------------------------------------
      ce_n10_ncar_lll = MAX( 1.e-3_wp * ( 34.6_wp * psqrtcdn10 ) , Cx_min )
      !
   END FUNCTION ce_n10_ncar_lll


   FUNCTION psi_m_ncar( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for momentum
      !!    !! Psis, L&Y 2004, Eq. (8c), (8d), (8e)
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m_ncar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zx2, zx, zpsi_unst, zpsi_stab,  zstab   ! local scalars
      !!----------------------------------------------------------------------------------
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            zta = pzeta(ji,jj)
            !
            zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 - 16z)^0.5
            zx2 = MAX( zx2 , 1._wp )
            zx  = SQRT(zx2)                          ! (1 - 16z)^0.25
            zpsi_unst = 2._wp*LOG( (1._wp + zx )*0.5_wp )   &
               &            + LOG( (1._wp + zx2)*0.5_wp )   &
               &          - 2._wp*ATAN(zx) + rpi*0.5_wp
            !
            zpsi_stab = -5._wp*zta
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
            !
            psi_m_ncar(ji,jj) =          zstab  * zpsi_stab &  ! (zta > 0) Stable
               &              + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable
            !
            !
      END_2D
   END FUNCTION psi_m_ncar


   FUNCTION psi_h_ncar( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !!    !! Psis, L&Y 2004, Eq. (8c), (8d), (8e)
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h_ncar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zta, zx2, zpsi_unst, zpsi_stab, zstab  ! local scalars
      !!----------------------------------------------------------------------------------
      !
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            !
            zta = pzeta(ji,jj)
            !
            zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 -16z)^0.5
            zx2 = MAX( zx2 , 1._wp )
            zpsi_unst = 2._wp*LOG( 0.5_wp*(1._wp + zx2) )
            !
            zpsi_stab = -5._wp*zta
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
            !
            psi_h_ncar(ji,jj) =          zstab  * zpsi_stab &  ! (zta > 0) Stable
               &              + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable
            !
      END_2D
   END FUNCTION psi_h_ncar

   FUNCTION psi_m_ncar_lll( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for momentum
      !!    !! Psis, L&Y 2004, Eq. (8c), (8d), (8e)
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp) :: psi_m_ncar_lll
      REAL(wp), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zx2, zx, zpsi_unst, zpsi_stab,  zstab   ! local scalars
      !!----------------------------------------------------------------------------------
      
      zta = pzeta
      !
      zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 - 16z)^0.5
      zx2 = MAX( zx2 , 1._wp )
      zx  = SQRT(zx2)                          ! (1 - 16z)^0.25
      zpsi_unst = 2._wp*LOG( (1._wp + zx )*0.5_wp )   &
         &            + LOG( (1._wp + zx2)*0.5_wp )   &
         &          - 2._wp*ATAN(zx) + rpi*0.5_wp
      !
      zpsi_stab = -5._wp*zta
      !
      zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
      !
      psi_m_ncar_lll =          zstab  * zpsi_stab &  ! (zta > 0) Stable
         &              + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable

   END FUNCTION psi_m_ncar_lll


   FUNCTION psi_h_ncar_lll( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !!    !! Psis, L&Y 2004, Eq. (8c), (8d), (8e)
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp) :: psi_h_ncar_lll
      REAL(wp), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zta, zx2, zpsi_unst, zpsi_stab, zstab  ! local scalars
      !!----------------------------------------------------------------------------------
      !
     
      !
      zta = pzeta
      !
      zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 -16z)^0.5
      zx2 = MAX( zx2 , 1._wp )
      zpsi_unst = 2._wp*LOG( 0.5_wp*(1._wp + zx2) )
      !
      zpsi_stab = -5._wp*zta
      !
      zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
      !
      psi_h_ncar_lll =          zstab  * zpsi_stab &  ! (zta > 0) Stable
         &              + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable
            !
   END FUNCTION psi_h_ncar_lll

   !!======================================================================
END MODULE sbcblk_algo_ncar
