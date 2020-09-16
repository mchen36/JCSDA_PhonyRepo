!
! Helper module containing the foam-related utility routines for the
! CRTM implementation of FASTEM4 and FASTEM5
!
!
! CREATION HISTORY:
!       Written by:     Original FASTEM1-5 authors
!
!       Refactored by:  Paul van Delst, November 2011
!                       paul.vandelst@noaa.gov
!

MODULE MW_Ocean_Foam_Module

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds     , ONLY: fp => CSEM_fp
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: MW_ANN_Foam_Coverage
  PUBLIC :: MW_ANN_Foam_Coverage_TL
  PUBLIC :: MW_ANN_Foam_Coverage_AD
  PUBLIC :: MW_ANN_Foam_Reflectivity


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: Foam_Utility_Module.f90 18621 2012-04-09 01:17:04Z paul.vandelst@noaa.gov $'

  ! Literal constants
  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp
  REAL(fp), PARAMETER :: TWO  = 2.0_fp
  

CONTAINS


  ! ===================================================================
  ! Foam coverage.
  !
  !   Monahan, E.C., and O'Muircheartaigh, I.G., (1986)
  !     Whitecaps and the passive remote sensing of the ocean surface,
  !     International Journal of Remote Sensing, 7, pp627-642.
  !
  ! The neutral stability condition is used here (i.e. the difference
  ! between the skin and air temperature is assumed to be zero) so
  ! that the form of the foam coverage equation is the same as in
  ! Tang (1974) and Liu et al. (1998)..
  !
  !   Liu, Q. et al. (1998) Monte Carlo simulations of the
  !     microwave emissivity of the sea surface.
  !     JGR, 103(C11), pp24983-24989
  !
  !   Tang, C. (1974) The effect of droplets in the air-sea
  !     transition zone on the sea brightness temperature.
  !     J. Phys. Oceanography, 4, pp579-593.
  !
  ! ===================================================================
  ! Forward model
  SUBROUTINE MW_ANN_Foam_Coverage(wind_speed, coverage)
    REAL(fp)              , INTENT(IN)  :: wind_speed
    REAL(fp)              , INTENT(OUT) :: coverage
    REAL(fp) :: FCCoeff(2) = (/1.95E-5, 2.55/)  !CRTM
    IF ( wind_speed < ZERO ) THEN
      coverage = ZERO
      RETURN
    END IF
    coverage = FCCoeff(1) * (wind_speed**FCCoeff(2))
   END SUBROUTINE MW_ANN_Foam_Coverage
  
  ! Tangent-linear model
  SUBROUTINE MW_ANN_Foam_Coverage_TL(wind_speed, wind_speed_TL, coverage_TL)
    REAL(fp)              , INTENT(IN)  :: wind_speed
    REAL(fp)              , INTENT(IN)  :: wind_speed_TL
    REAL(fp)              , INTENT(OUT) :: coverage_TL
    REAL(fp) :: FCCoeff(2) = (/1.95E-5, 2.55/)  !CRTM
    IF ( wind_speed < ZERO ) THEN
      coverage_TL = ZERO
      RETURN
    END IF
    coverage_TL = FCCoeff(1)*FCCoeff(2) * (wind_speed**(FCCoeff(2)-ONE)) * wind_speed_TL
  END SUBROUTINE MW_ANN_Foam_Coverage_TL

  ! Adjoint model
  SUBROUTINE MW_ANN_Foam_Coverage_AD(wind_speed, coverage_AD, wind_speed_AD)
    REAL(fp)              , INTENT(IN)     :: wind_speed     ! Input
    REAL(fp)              , INTENT(IN OUT) :: coverage_AD    ! Input
    REAL(fp)              , INTENT(IN OUT) :: wind_speed_AD  ! Output
    REAL(fp) :: FCCoeff(2) = (/1.95E-5, 2.55/)  !CRTM
    IF ( wind_speed < ZERO ) THEN
      coverage_AD = ZERO
      RETURN
    END IF
    wind_speed_AD = wind_speed_AD + &
                    FCCoeff(1)*FCCoeff(2) * (wind_speed**(FCCoeff(2)-ONE)) * coverage_AD
    coverage_AD = ZERO
  END SUBROUTINE MW_ANN_Foam_Coverage_AD


  ! =============================================================
  ! Foam reflectivity
  !
  ! See section d in
  !
  !   Kazumori, M. et al. (2008) Impact Study of AMSR-E Radiances
  !     in the NCEP Global Data Assimilation System,
  !     Monthly Weather Review, 136, pp541-559
  !
  ! Function dependence is on zenith angle only so no TL
  ! or AD routine.
  ! =============================================================
  SUBROUTINE MW_ANN_Foam_Reflectivity( &
    Zenith_Angle, &
    Frequency   , &
    Rv          , &
    Rh            )
    ! Arguments
    REAL(fp)              , INTENT(IN)  :: Zenith_Angle
    REAL(fp)              , INTENT(IN)  :: Frequency
    REAL(fp)              , INTENT(OUT) :: Rv, Rh
   
    REAL(fp) :: FRCoeff(6) = (/0.93, -1.748E-3, -7.336E-5, -1.004E-7, 0.40, -5.0E-2/)
    
    ! Local variables
    REAL(fp) :: factor
    
    ! The vertical component is a fixed value
    Rv = ONE - FRCoeff(1)  ! Fixed nadir emissivity
    
    ! The horizontal component uses a regression equation
    ! to compute a factor modifying the nadir emissivity
    factor = ONE + Zenith_Angle*(FRCoeff(2) + &
                     Zenith_Angle*(FRCoeff(3) + &
                       Zenith_Angle*FRCoeff(4)  )  )
    Rh = ONE - factor*FRCoeff(1)
    
    ! Frequency correction
    factor = FRCoeff(5) * EXP(FRCoeff(6)*Frequency)
    Rv = Rv * factor
    Rh = Rh * factor
   END SUBROUTINE MW_ANN_Foam_Reflectivity  

END MODULE MW_Ocean_Foam_Module
