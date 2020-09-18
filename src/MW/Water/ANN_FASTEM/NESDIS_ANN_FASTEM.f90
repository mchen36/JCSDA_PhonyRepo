!
! NESDIS_ANN_FASTEM
!
! Module implementing two-scale ocean surface emissivity model. 
! Machine learning regressions are employed to approximate those 
! calculations which are computationally expensive in the original
! physical space.
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 30-08-2017
!                       ming.chen@noaa.gov

MODULE NESDIS_ANN_FASTEM

  USE CSEM_Type_Kinds,           ONLY: fp => CSEM_fp 
  USE CSEM_Exception_Handler,    ONLY: SUCCESS, FAILURE 
  USE ANN_MWwaterCoeff_Reader,   ONLY: ANN_MWwaterCoeff_Load, &
                                       ANN_MWwaterCoeff_CleanUp
       
  USE ANN_MLPCoeff_Module,       ONLY: ANN_MLPCoeff_Type,   &
                                       ANN_MLP_Recon,       &
                                       ANN_MLP_Recon_TL,    &
                                       ANN_MLP_Recon_AD,    &
                                       nn_iVar_type => iVar_type

  USE MW_Ocean_Foam_Module,      ONLY: MW_ANN_Foam_Coverage,    &
                                       MW_ANN_Foam_Coverage_TL, &
                                       MW_ANN_Foam_Coverage_AD, &
                                       MW_ANN_Foam_Reflectivity
  USE ANN_Reflection_Correction_Module, ONLY: rcVar_type => iVar_type, &
                                       ANN_Reflection_Correction   , &
                                       ANN_Reflection_Correction_TL, &
                                       ANN_Reflection_Correction_AD

  !USE Ocean_Permittivity

  USE FASTEM_Fresnel, &
    ONLY: fVar_type => iVar_type , &
          FASTEM_Fresnel_Reflectivity   , &
          FASTEM_Fresnel_Reflectivity_TL, &
          FASTEM_Fresnel_Reflectivity_AD

  USE Liu, &
    ONLY: pVar_type => iVar_type, &
          Ocean_Permittivity    => Liu_Ocean_Permittivity   , &
          Ocean_Permittivity_TL => Liu_Ocean_Permittivity_TL, &
          Ocean_Permittivity_AD => Liu_Ocean_Permittivity_AD
 
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Compute_ANN_Fastem
  PUBLIC :: Compute_ANN_Fastem_TL
  PUBLIC :: Compute_ANN_Fastem_AD
  PUBLIC :: ANN_Fastem_INIT
  PUBLIC :: ANN_Fastem_Destroy
  PUBLIC :: iVar_type
  PUBLIC :: ANN_MWwaterCoeff_INIT
  
  REAL(fp), PARAMETER :: PI  = 3.141592653589793238462643_fp
  REAL(fp), PARAMETER :: ONE = 1.0_fp, ZERO =0.0_fp 
  ! Stokes component information
  ! ...The number of Stokes components
  INTEGER, PARAMETER :: N_STOKES = 4
  ! ...The vector indices
  INTEGER, PARAMETER :: Iv_IDX = 1 ! Describes vertical polarization
  INTEGER, PARAMETER :: Ih_IDX = 2 ! Describes horizontal polarization
  INTEGER, PARAMETER :: U_IDX  = 3 ! Describes plane of polarization
  INTEGER, PARAMETER :: V_IDX  = 4 ! Describes ellipticity of polarization

  ! Group name of the coefficient set 
  INTEGER,      PARAMETER :: N_SUBGRP = 14
  CHARACTER(LEN=15) :: SUB_GROUP_NAME(N_SUBGRP) = &
                   [CHARACTER(LEN=15):: 'TV_MEAN', 'TH_MEAN', &
                   'TV_Azimu_Cos1', 'TV_Azimu_Cos2','TV_Azimu_Cos3', &
                   'TH_Azimu_Cos1', 'TH_Azimu_Cos2','TH_Azimu_Cos3', &
                   'U_Azimu_Sin1 ', 'U_Azimu_Sin2 ','U_Azimu_Sin3 ', &
                   'V_Azimu_Sin1 ', 'V_Azimu_Sin2 ','V_Azimu_Sin3 ' ]
  

  ! --------------------------------------
  ! Structure definition to hold internal
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    ! Validity indicator
    LOGICAL :: Is_Valid = .FALSE.
    ! Forward model input values
    REAL(fp) :: Frequency    = ZERO
    REAL(fp) :: Zenith_Angle = ZERO
    REAL(fp) :: Temperature  = ZERO
    REAL(fp) :: Salinity     = ZERO
    REAL(fp) :: Wind_Speed   = ZERO
 
    ! The azimuth angle term
    REAL(fp) :: Azimuth = ZERO
    ! The cos azimuth angle term
    REAL(fp) :: C_Azimuth(3) = ZERO
    ! The sin azimuth angle term
    REAL(fp) :: S_Azimuth(3) = ZERO
    ! The zenith angle term
    REAL(fp) :: cos_z   = ONE
    ! The permittivity term
    COMPLEX(fp) :: Permittivity = ZERO
    ! The Fresnel reflectivity terms
    REAL(fp) :: Rv_Fresnel = ZERO
    REAL(fp) :: Rh_Fresnel = ZERO
    ! Foam Terms
    REAL(fp) :: Rv_Foam = ZERO
    REAL(fp) :: Rh_Foam = ZERO
    REAL(fp) :: Foam_Cover = ZERO
    ! Coherent correction reflectivities
    REAL(fp) :: Rv_Coherent = ZERO
    REAL(fp) :: Rh_Coherent = ZERO
    ! Incoherent correction coefficients
    REAL(fp) :: Y(3,N_STOKES)

    ! Final reflectivities
    REAL(fp) :: Rv = ZERO
    REAL(fp) :: Rh = ZERO
    ! Azimuthal emissivity
    REAL(fp) :: e_Azimuth(N_STOKES) = ZERO
    ! Anisotropic downward radiation correction
    REAL(fp) :: Rv_Mod = ZERO
    REAL(fp) :: Rh_Mod = ZERO
    ! The final emissivity
    REAL(fp) :: e(N_STOKES) = ZERO
    ! ...Optional
    LOGICAL  :: Transmittance_Valid = .FALSE.
    REAL(fp) :: Transmittance       = ZERO

    ! Internal variables for subcomponents
    TYPE(pVar_type)    :: pVar
    TYPE(fVar_type)    :: fVar
    TYPE(rcVar_type)   :: rcVar
  
    TYPE(nn_iVar_type) :: nn_iVar(N_SUBGRP)
  END TYPE iVar_type

  ! --------------------------------------------------
  ! The shared microwave water surface emissivity coefficient data
  ! --------------------------------------------------
  TYPE(ANN_MLPCoeff_Type), TARGET, SAVE :: NN_MWwaterC(N_SUBGRP)

  LOGICAL :: ANN_MWwaterCoeff_INIT

CONTAINS


  SUBROUTINE Compute_ANN_Fastem( &
    Frequency    , &  ! Input
    Zenith_Angle , &  ! Input
    Azimuth_Angle, &  ! Input
    Temperature  , &  ! Input
    Salinity     , &  ! Input
    Wind_Speed   , &  ! Input
    iVar         , &  ! Internal variable output
    Emissivity   , &  ! Output
    Reflectivity , &  ! Output
    Transmittance  )  ! Optional Input

    ! Arguments
    REAL(fp),                INTENT(IN)  :: Frequency
    REAL(fp),                INTENT(IN)  :: Zenith_Angle
    REAL(fp),                INTENT(IN)  :: Temperature
    REAL(fp),                INTENT(IN)  :: Salinity
    REAL(fp),                INTENT(IN)  :: Wind_Speed
    REAL(fp),                INTENT(IN)  :: Azimuth_Angle
    TYPE(iVar_type),         INTENT(INOUT) :: iVar
    REAL(fp),                INTENT(OUT) :: Emissivity(:)
    REAL(fp),                INTENT(OUT) :: Reflectivity(:)
    REAL(fp), OPTIONAL,      INTENT(IN)  :: Transmittance
 
    REAL(fp) :: X(3), Y1,Y2
    REAL(fp) :: I0(4) = (/ONE, ONE, ZERO, ZERO/)
    INTEGER  :: k, ipol


    ! Setup
    ! ...Save forward input variables for TL and AD calculations
    iVar%Frequency    = Frequency
    iVar%Zenith_Angle = Zenith_Angle
    iVar%Temperature  = Temperature
    iVar%Salinity     = Salinity
    iVar%Wind_Speed   = Wind_Speed

    ! ...Save derived variables
    iVar%cos_z = COS(Zenith_Angle/180.0*PI)
    iVar%Azimuth = Azimuth_Angle/180.0*PI

    ! Permittivity calculation
    CALL Ocean_Permittivity( Temperature, Salinity, Frequency, &
                             iVar%Permittivity, &
                             iVar%pVar )

    ! Fresnel reflectivity calculation
    CALL FASTEM_Fresnel_Reflectivity( iVar%Permittivity, iVar%cos_z, &
                               iVar%Rv_Fresnel, iVar%Rh_Fresnel, &
                               iVar%fVar )


    ! Foam reflectivity calculation
    CALL MW_ANN_Foam_Reflectivity( &
           Zenith_Angle, &
           Frequency   , &
           iVar%Rv_Foam, &
           iVar%Rh_Foam  )

    ! Foam coverage calculation
    CALL MW_ANN_Foam_Coverage( Wind_Speed,  iVar%Foam_Cover )
    X = (/ Frequency, Wind_Speed, Zenith_Angle/)
    !X = (/ Zenith_Angle, Wind_Speed, Frequency/)
    DO k = 1, 3
       X(k) = (X(k) - NN_MWwaterC(1)%vmin(k))  /  &
           (NN_MWwaterC(1)%vmax(k) - NN_MWwaterC(1)%vmin(k))
    END DO   
    CALL ANN_MLP_Recon(NN_MWwaterC(1), X, Y1, iVar%nn_iVar(1))
    CALL ANN_MLP_Recon(NN_MWwaterC(2), X, Y2, iVar%nn_iVar(2))

    iVar%Rv_Coherent = Y1/NN_MWwaterC(1)%fscale 
    iVar%Rh_Coherent = Y2/NN_MWwaterC(2)%fscale 
    iVar%Rv = exp(iVar%Rv_Coherent) * iVar%Rv_Fresnel
    iVar%Rh = exp(iVar%Rh_Coherent) * iVar%Rh_Fresnel

    iVar%C_Azimuth = (/cos(iVar%Azimuth), cos(2.00*iVar%Azimuth), cos(3.00*iVar%Azimuth)/)
    iVar%S_Azimuth = (/sin(iVar%Azimuth), sin(2.00*iVar%Azimuth), sin(3.00*iVar%Azimuth)/)

    DO k = 3, 12, 3
       ipol = (k-3)/3+1
       CALL ANN_MLP_Recon(NN_MWwaterC(k+0), X, iVar%Y(1,ipol), iVar%nn_iVar(k+0))
       CALL ANN_MLP_Recon(NN_MWwaterC(k+1), X, iVar%Y(2,ipol), iVar%nn_iVar(k+1))
       CALL ANN_MLP_Recon(NN_MWwaterC(k+2), X, iVar%Y(3,ipol), iVar%nn_iVar(k+2))
       iVar%Y(1,ipol) = iVar%Y(1,ipol)/NN_MWwaterC(k+0)%fscale
       iVar%Y(2,ipol) = iVar%Y(2,ipol)/NN_MWwaterC(k+1)%fscale
       iVar%Y(3,ipol) = iVar%Y(3,ipol)/NN_MWwaterC(k+2)%fscale
      
        IF( ipol < 3)  THEN 
          Reflectivity(ipol) = DOT_PRODUCT(iVar%Y(:,ipol),iVar%C_Azimuth)
        ELSE
          Reflectivity(ipol) = DOT_PRODUCT(iVar%Y(:,ipol),iVar%S_Azimuth)
        ENDIF
    END DO

    ! Compute the first two Stokes components of the reflectivity

    Reflectivity(Iv_IDX)  =  Reflectivity(Iv_IDX) + &
                           (ONE-iVar%Foam_Cover) * iVar%Rv + iVar%Foam_Cover*iVar%Rv_Foam
    Reflectivity(Ih_IDX)  =  Reflectivity(Ih_IDX) + &
                           (ONE-iVar%Foam_Cover) * iVar%Rh + iVar%Foam_Cover*iVar%Rh_Foam

    ! Compute the full Stokes components of the Emissivity
    Emissivity       =  I0 - Reflectivity
    iVar%e           =  Emissivity

    ! Anisotropic downward radiation correction calculation
    iVar%Transmittance_Valid = .FALSE.
    iVar%Rv_Mod = ONE
    iVar%Rh_Mod = ONE
    IF ( PRESENT(Transmittance) ) THEN
      IF ( Transmittance > ZERO .AND. Transmittance < ONE ) THEN
        CALL ANN_Reflection_Correction( &
               iVar%Frequency , &
               iVar%cos_z     , &
               iVar%Wind_Speed, &
               Transmittance  , &
               iVar%Rv_Mod    , &
               iVar%Rh_Mod    , &
               iVar%rcVar       )
        iVar%Transmittance_Valid = .TRUE.
        iVar%Transmittance       = Transmittance
        Reflectivity(Iv_IDX)     = iVar%Rv_Mod *  (ONE - Emissivity(Iv_IDX))
        Reflectivity(Ih_IDX)     = iVar%Rh_Mod *  (ONE - Emissivity(Ih_IDX))
      END IF
    END IF

    iVar%Is_Valid    = .TRUE.
  END SUBROUTINE Compute_ANN_Fastem

  SUBROUTINE Compute_ANN_Fastem_TL( &
    Azimuth_Angle_TL, &  ! Input
    Temperature_TL  , &  ! Input
    Salinity_TL     , &  ! Input
    Wind_Speed_TL   , &  ! Input
    iVar            , &  ! Internal variable output
    Emissivity_TL   , &  ! Output
    Reflectivity_TL , &  ! Output
    Transmittance_TL  )  ! Optional Input

    ! Arguments
    REAL(fp),                INTENT(IN)  :: Temperature_TL
    REAL(fp),                INTENT(IN)  :: Salinity_TL
    REAL(fp),                INTENT(IN)  :: Wind_Speed_TL
    REAL(fp),                INTENT(IN)  :: Azimuth_Angle_TL
    TYPE(iVar_type),         INTENT(IN)  :: iVar
    REAL(fp),                INTENT(OUT) :: Emissivity_TL(:)
    REAL(fp),                INTENT(OUT) :: Reflectivity_TL(:)
    REAL(fp), OPTIONAL,      INTENT(IN)  :: Transmittance_TL
 
    INTEGER  :: k, ipol

    REAL(fp) :: Rv_Fresnel_TL, Rh_Fresnel_TL
    REAL(fp) :: Rv_Foam_TL,  Rh_Foam_TL, Foam_Cover_TL 
    REAL(fp) :: Rv_Coherent_TL, Rh_Coherent_TL
    REAL(fp) :: Rv_TL, Rh_TL
    REAL(fp) :: Rv_Mod_TL,  Rh_Mod_TL
    REAL(fp) :: X_TL(3), Y_TL(3)

    COMPLEX(fp) :: Permittivity_TL

    ! Check internal structure
    IF ( .NOT. iVar%Is_Valid ) THEN
      Emissivity_TL   = ZERO
      Reflectivity_TL = ZERO
      RETURN
    END IF

 
    ! Permittivity calculation
    CALL Ocean_Permittivity_TL( Temperature_TL, Salinity_TL, iVar%Frequency, &
                                Permittivity_TL, &
                                iVar%pVar)

    ! Fresnel reflectivity calculation
    CALL FASTEM_Fresnel_Reflectivity_TL( permittivity_TL, iVar%cos_z, &
                                  Rv_Fresnel_TL, Rh_Fresnel_TL, &
                                  iVar%fVar )

    ! Foam reflectivity "calculation"
    Rv_Foam_TL = ZERO
    Rh_Foam_TL = ZERO

    ! Foam coverage calculation
    CALL MW_ANN_Foam_Coverage_TL( &
           iVar%Wind_Speed, &
           Wind_Speed_TL, &
           Foam_Cover_TL )

    X_TL(1) = ZERO
    X_TL(2) = Wind_Speed_TL / (NN_MWwaterC(1)%vmax(2) - NN_MWwaterC(1)%vmin(2))
    X_TL(3) = ZERO
   
    CALL ANN_MLP_Recon_TL(NN_MWwaterC(1), X_TL, Y_TL(1), iVar%nn_iVar(1))
    CALL ANN_MLP_Recon_TL(NN_MWwaterC(2), X_TL, Y_TL(2), iVar%nn_iVar(2))

    Rv_Coherent_TL = Y_TL(1)/NN_MWwaterC(1)%fscale 
    Rh_Coherent_TL = Y_TL(2)/NN_MWwaterC(2)%fscale 
    Rv_TL = (Rv_Coherent_TL * iVar%Rv_Fresnel + Rv_Fresnel_TL)*exp(iVar%Rv_Coherent)
    Rh_TL = (Rh_Coherent_TL * iVar%Rh_Fresnel + Rh_Fresnel_TL)*exp(iVar%Rh_Coherent)
    

    Y_TL = ZERO
    DO k = 3, 12, 3
       ipol = (k-3)/3+1
       CALL ANN_MLP_Recon_TL(NN_MWwaterC(k+0), X_TL, Y_TL(1), iVar%nn_iVar(k+0))
       CALL ANN_MLP_Recon_TL(NN_MWwaterC(k+1), X_TL, Y_TL(2), iVar%nn_iVar(k+1))
       CALL ANN_MLP_Recon_TL(NN_MWwaterC(k+2), X_TL, Y_TL(3), iVar%nn_iVar(k+2))
       Y_TL(1) = Y_TL(1)/NN_MWwaterC(k+0)%fscale
       Y_TL(2) = Y_TL(2)/NN_MWwaterC(k+1)%fscale
       Y_TL(3) = Y_TL(3)/NN_MWwaterC(k+2)%fscale
      
        IF(ipol < 3)  THEN 
           Reflectivity_TL(ipol) = DOT_PRODUCT(Y_TL,iVar%C_Azimuth) + &
                           DOT_PRODUCT(iVar%Y(:,ipol),iVar%S_Azimuth*(/-1.0,-2.0,-3.0/)*Azimuth_Angle_TL)
        ELSE
           Reflectivity_TL(ipol) = DOT_PRODUCT(Y_TL,iVar%S_Azimuth) + &
                           DOT_PRODUCT(iVar%Y(:,ipol),iVar%C_Azimuth*(/1.0, 2.0, 3.0/)*Azimuth_Angle_TL)

        ENDIF
    END DO


    ! Compute the first two Stokes components of the reflectivity
                           
    Reflectivity_TL(Iv_IDX)  = Reflectivity_TL(Iv_IDX) + (ONE-iVar%Foam_Cover) * Rv_TL + &
                               Foam_Cover_TL * (iVar%Rv_Foam - iVar%Rv)
    Reflectivity_TL(Ih_IDX)  = Reflectivity_TL(Ih_IDX) + (ONE-iVar%Foam_Cover) * Rh_TL  + &
                               Foam_Cover_TL * (iVar%Rh_Foam - iVar%Rh)

    ! Compute the full Stokes components of the Emissivity
    Emissivity_TL       =  - Reflectivity_TL

    ! Anisotropic downward radiation correction calculation
    IF ( PRESENT(Transmittance_TL) .AND. iVar%Transmittance_Valid ) THEN
      CALL ANN_Reflection_Correction_TL( &
             Wind_Speed_TL   , &
             Transmittance_TL, &
             Rv_Mod_TL       , &
             Rh_Mod_TL       , &
             iVar%rcVar        )

      ! ...reflectivities
      Reflectivity_TL(Iv_IDX)      = (ONE-iVar%e(Iv_IDX))*Rv_Mod_TL - iVar%Rv_Mod*Emissivity_TL(Iv_IDX)
      Reflectivity_TL(Ih_IDX)      = (ONE-iVar%e(Ih_IDX))*Rh_Mod_TL - iVar%Rh_Mod*Emissivity_TL(Ih_IDX)

    ELSE
      Rv_Mod_TL = ZERO
      Rh_Mod_TL = ZERO
    END IF
   
  END SUBROUTINE Compute_ANN_Fastem_TL


  SUBROUTINE Compute_ANN_Fastem_AD( &
    Azimuth_Angle_AD, &  ! Input
    Temperature_AD  , &  ! Input
    Salinity_AD     , &  ! Input
    Wind_Speed_AD   , &  ! Input
    iVar            , &  ! Internal variable output
    Emissivity_AD   , &  ! Output
    Reflectivity_AD , &  ! Output
    Transmittance_AD  )  ! Optional Input

    ! Arguments
    REAL(fp),                INTENT(INOUT)  :: Temperature_AD
    REAL(fp),                INTENT(INOUT)  :: Salinity_AD
    REAL(fp),                INTENT(INOUT)  :: Wind_Speed_AD
    REAL(fp),                INTENT(INOUT)  :: Azimuth_Angle_AD
    TYPE(iVar_type),         INTENT(IN)     :: iVar
    REAL(fp),                INTENT(INOUT)  :: Emissivity_AD(:)
    REAL(fp),                INTENT(INOUT)  :: Reflectivity_AD(:)
    REAL(fp), OPTIONAL,      INTENT(INOUT)  :: Transmittance_AD
 
    INTEGER  :: k, ipol

    REAL(fp) :: Rv_Fresnel_AD, Rh_Fresnel_AD
    REAL(fp) :: Rv_Foam_AD,  Rh_Foam_AD, Foam_Cover_AD 
    REAL(fp) :: Rv_Coherent_AD, Rh_Coherent_AD
    REAL(fp) :: Rv_AD, Rh_AD
    REAL(fp) :: Rv_Mod_AD,  Rh_Mod_AD
    REAL(fp) :: X_AD(3), Y_AD(3)

    COMPLEX(fp) :: Permittivity_AD

    ! Check internal structure
    IF ( .NOT. iVar%Is_Valid ) THEN
      Emissivity_AD   = ZERO
      Reflectivity_AD = ZERO
      RETURN
    END IF
    ! Anisotropic downward radiation correction calculation
    IF ( PRESENT(Transmittance_AD) .AND. iVar%Transmittance_Valid ) THEN
      ! Compute the full Stokes components of the Emissivity
      Rh_Mod_AD               = (ONE-iVar%e(Ih_IDX))*Reflectivity_AD(Ih_IDX) 
      Rv_Mod_AD               = (ONE-iVar%e(Iv_IDX))*Reflectivity_AD(Iv_IDX)
      Emissivity_AD(Iv_IDX)   =  Emissivity_AD(Iv_IDX) - iVar%Rv_Mod*Reflectivity_AD(Iv_IDX)
      Emissivity_AD(Ih_IDX)   =  Emissivity_AD(Ih_IDX) - iVar%Rh_Mod*Reflectivity_AD(Ih_IDX)
      CALL ANN_Reflection_Correction_AD( &
             Rv_Mod_AD       , &
             Rh_Mod_AD       , &
             Wind_Speed_AD   , &
             Transmittance_AD, &
             iVar%rcVar        ) 
      Reflectivity_AD(Iv_IDX)  = ZERO
      Reflectivity_AD(Ih_IDX)  = ZERO
    ELSE
      Rv_Mod_AD = ZERO
      Rh_Mod_AD = ZERO
    END IF
 
    Reflectivity_AD       = Reflectivity_AD - Emissivity_AD
  
    Foam_Cover_AD = ZERO
    ! Compute the first two Stokes components of the reflectivity
    Foam_Cover_AD =  Foam_Cover_AD + Reflectivity_AD(Iv_IDX) * (iVar%Rv_Foam - iVar%Rv)
    Foam_Cover_AD =  Foam_Cover_AD + Reflectivity_AD(Ih_IDX) * (iVar%Rh_Foam - iVar%Rh)
    Rv_AD  = ZERO ;  Rh_AD = ZERO
    Rv_AD = Rv_AD + (ONE-iVar%Foam_Cover) * Reflectivity_AD(Iv_IDX)
    Rh_AD = Rh_AD + (ONE-iVar%Foam_Cover) * Reflectivity_AD(Ih_IDX)
    
    Y_AD = ZERO ; X_AD = ZERO
    DO k = 3, 12, 3
       ipol = (k-3)/3+1
       IF( ipol < 3)  THEN 
           Y_AD =  Reflectivity_AD(ipol) * iVar%C_Azimuth
           Azimuth_Angle_AD = Azimuth_Angle_AD +  &
                             DOT_PRODUCT(iVar%Y(:,ipol),iVar%S_Azimuth*(/-1.0, -2.0, -3.0/)*Reflectivity_AD(ipol))
       ELSE
           Y_AD = Reflectivity_AD(ipol) * iVar%S_Azimuth
           Azimuth_Angle_AD = Azimuth_Angle_AD +  &
                             DOT_PRODUCT(iVar%Y(:,ipol),iVar%C_Azimuth*(/1.0, 2.0, 3.0/)*Reflectivity_AD(ipol))
       ENDIF
       Y_AD(1) = Y_AD(1)/NN_MWwaterC(k+0)%fscale
       Y_AD(2) = Y_AD(2)/NN_MWwaterC(k+1)%fscale
       Y_AD(3) = Y_AD(3)/NN_MWwaterC(k+2)%fscale

       CALL ANN_MLP_Recon_AD(NN_MWwaterC(k+0), X_AD, Y_AD(1), iVar%nn_iVar(k+0))
       CALL ANN_MLP_Recon_AD(NN_MWwaterC(k+1), X_AD, Y_AD(2), iVar%nn_iVar(k+1))
       CALL ANN_MLP_Recon_AD(NN_MWwaterC(k+2), X_AD, Y_AD(3), iVar%nn_iVar(k+2))

    END DO
   
    Rv_Coherent_AD =  Rv_AD * iVar%Rv_Fresnel* exp(iVar%Rv_Coherent)
    Rh_Coherent_AD =  Rh_AD * iVar%Rh_Fresnel* exp(iVar%Rh_Coherent)
    Rv_Fresnel_AD  =  Rv_AD * exp(iVar%Rv_Coherent)
    Rh_Fresnel_AD  =  Rh_AD * exp(iVar%Rh_Coherent)
    
    Y_AD(1)  = Rv_Coherent_AD/NN_MWwaterC(1)%fscale 
    Y_AD(2)  = Rh_Coherent_AD/NN_MWwaterC(2)%fscale 
   
    CALL ANN_MLP_Recon_AD(NN_MWwaterC(1), X_AD, Y_AD(1), iVar%nn_iVar(1))
    CALL ANN_MLP_Recon_AD(NN_MWwaterC(2), X_AD, Y_AD(2), iVar%nn_iVar(2))

    Wind_Speed_AD =Wind_Speed_AD + X_AD(2)/ (NN_MWwaterC(1)%vmax(2) - NN_MWwaterC(1)%vmin(2))
    X_AD  = ZERO
   

    ! Foam coverage calculation
    CALL MW_ANN_Foam_Coverage_AD( &
           iVar%Wind_Speed, &
           Foam_Cover_AD, &
           Wind_Speed_AD )

    ! Foam reflectivity "calculation"
    Rv_Foam_AD = ZERO
    Rh_Foam_AD = ZERO

    ! Fresnel reflectivity calculation
    Permittivity_AD = ZERO
    CALL FASTEM_Fresnel_Reflectivity_AD( Rv_Fresnel_AD, Rh_Fresnel_AD, iVar%cos_z, &
                                  Permittivity_AD, &
                                  iVar%fVar )
    ! Permittivity calculation
    CALL Ocean_Permittivity_AD( Permittivity_AD, iVar%Frequency, &
                                Temperature_AD, Salinity_AD, &
                                iVar%pVar )

    Emissivity_AD   = ZERO
    Reflectivity_AD = ZERO
  END SUBROUTINE Compute_ANN_Fastem_AD


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!      ANN_FASTEM_Init
!
! PURPOSE:
!       Function to load FASTEM coefficient NETDCF files
!
!       This function must be called before calling other FASTEM functions
!
! CALLING SEQUENCE:
!       Error_Status = ANN_FASTEM_Init(    
!                              MWwaterCoeff_File                      , &  ! Output
!                              Version)                                 &  ! Input
! INPUTS:
!
!      MWwaterCoeff_File:FASTEM coefficient file (full-path)
!                        UNITS:      N/A
!                        TYPE:       CHARACTER
!                        DIMENSION:  SCALAR
!                        ATTRIBUTES: INTENT(INOUT)
!
!
!      Version:          version index
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!
!      ANN_MWwaterCoeff_INIT: Global-scope variable
!                        UNITS:      N/A
!                        TYPE:       LOGICAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: Global 
!
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the Message_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!
!:sdoc-:
!----------------------------------------------------------------------------------
  
  FUNCTION ANN_FASTEM_Init( &
    MWwaterCoeff_File)      &
  RESULT( Error_Status  )
    CHARACTER(LEN=*)  :: MWwaterCoeff_File
    INTEGER  :: Error_Status

    !PRINT*, 'loading MWwaterCoeff data from '//TRIM(MWwaterCoeff_File)
  
    IF(ANN_MWwaterCoeff_INIT)  THEN
       PRINT*, 'ANN_MWwaterCoeff already loaded, reloading ...'
       RETURN
    ENDIF
    
    Error_Status = ANN_MWwaterCoeff_Load( &
                   TRIM(MWwaterCoeff_File), &
                   NN_MWwaterC      , &  ! InOUT
                   SUB_GROUP_NAME   , &  ! Input
                   N_SUBGRP)             ! Input

    IF ( Error_Status /= SUCCESS ) THEN
        PRINT*, 'Error loading MWwaterCoeff data from '//TRIM(MWwaterCoeff_File)
        STOP
    END IF
    ANN_MWwaterCoeff_INIT = .TRUE.
  END FUNCTION ANN_FASTEM_Init 


  FUNCTION ANN_FASTEM_Destroy () RESULT( Error_Status  )
    INTEGER  :: Error_Status
    Error_Status = SUCCESS
    !PRINT*, ' Closeing MWwaterCoeff data ....'
    IF(.NOT. ANN_MWwaterCoeff_INIT) THEN 
        PRINT*, 'Error Closeing MWwaterCoeff data ....'
        Error_Status= FAILURE 
        RETURN
    ENDIF
    ! Destroy the structure
    Error_Status = ANN_MWwaterCoeff_CleanUp(NN_MWwaterC)
    
    ANN_MWwaterCoeff_INIT = .FALSE.
  END FUNCTION ANN_FASTEM_Destroy
 
END MODULE NESDIS_ANN_FASTEM
