!
! ANN_MLPCoeff_Module
!
! Module defining the ANN_MLPCoeff object.
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 30-08-2017
!                       ming.chen@noaa.gov

MODULE ANN_MLPCoeff_Module

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds         , ONLY: fp=>CSEM_fp 
  USE CSEM_Exception_Handler  , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message

  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Datatypes
  PUBLIC :: ANN_MLPCoeff_Type

  ! Procedures
  PUBLIC :: ANN_MLPCoeff_Destroy
  PUBLIC :: ANN_MLPCoeff_Create
  PUBLIC :: ANN_MLP_Recon
  PUBLIC :: ANN_MLP_Recon_TL
  PUBLIC :: ANN_MLP_Recon_AD

  PUBLIC :: iVar_type


  ! -----------------
  ! Module parameters
  ! -----------------
  ! String lengths
  INTEGER,  PARAMETER :: ML = 256 ! Message length
  INTEGER,  PARAMETER :: SL =  80 ! String length


  ! ---------------------------------
  ! NN MLPCoeff data type definition
  ! ---------------------------------
  !:tdoc+:
  TYPE :: ANN_MLPCoeff_Type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! type components
    INTEGER  :: n_features = 1              ! NN inputs
    INTEGER  :: n_neurons  = 1              ! NN neurons
    REAL(fp) :: fscale     = 1.0            ! scale factor
    
    REAL(fp), ALLOCATABLE :: w1(:,:)        ! NN hidden Layer 1 weights
    REAL(fp), ALLOCATABLE :: w2(:)          ! NN output Layer   weights
    REAL(fp), ALLOCATABLE :: b1(:)          ! NN hidden Layer 1 intercepts
    REAL(fp) :: b2                          ! NN output Layer   intercepts
    REAL(fp), ALLOCATABLE :: vmax(:)        ! Upper bounds of input features
    REAL(fp), ALLOCATABLE :: vmin(:)        ! Lower bounds of input features
  END TYPE ANN_MLPCoeff_Type
  !:tdoc-:

  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.

    ! type components
    INTEGER :: n_neurons = -1
    REAL(fp), ALLOCATABLE :: neurons(:)
  END TYPE iVar_type


CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                           ## PUBLIC PROCEDURES ##                          ##
!##                                                                            ##
!################################################################################
!################################################################################


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       ANN_MLPCoeff_Destroy
!
! PURPOSE:
!       Subroutine to re-initialize NN MLPCoeff objects.
!
! CALLING SEQUENCE:
!       CALL ANN_MLPCoeff_Destroy( ANN_MLPCoeff )
!
! OBJECTS:
!       ANN_MLPCoeff: Re-initialized MLPCoeff structure.
!                     UNITS:      N/A
!                     TYPE:       ANN_MLPCoeff_type
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE ANN_MLPCoeff_Destroy( self )
    TYPE(ANN_MLPCoeff_type), INTENT(IN OUT) :: self
    INTEGER  :: Error_Status

    DEALLOCATE(self%w1,             &
               self%w2,             &
               self%b1,       &
               self%vmax,              &
               self%vmin,              &
               STAT=Error_Status)
    self%Is_Allocated = .FALSE.
  END SUBROUTINE ANN_MLPCoeff_Destroy


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       ANN_MLPCoeff_Create
!
! PURPOSE:
!       Subroutine to create a valid instance of an ANN_MLPCoeff object.
!
! CALLING SEQUENCE:
!       CALL ANN_MLPCoeff_Create( ANN_MLPCoeff )
!
! OBJECTS:
!       ANN_MLPCoeff:     ANN_MLPCoeff_Type object structure.
!                           UNITS:      N/A
!                           TYPE:       ANN_MLPCoeff_Type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE ANN_MLPCoeff_Create( &
    self,           &
    n_neurons,      &
    n_features)

    ! Arguments

    INTEGER,   INTENT(IN)  :: n_features
    INTEGER,   INTENT(IN)  :: n_neurons
    TYPE(ANN_MLPCoeff_type), INTENT(IN OUT) :: self

    INTEGER  :: Error_Status

    self%Is_Allocated = .FALSE.
    ALLOCATE(self%w1(n_neurons, n_features),   &
             self%w2(n_neurons),               &
             self%b1(n_neurons),               &
             self%vmax(n_features),            &
             self%vmin(n_features),            &
             STAT=Error_Status)
    IF(Error_Status /= 0) THEN
      WRITE(*,*)"ANN_MLPCoeff_Create: Memory allocation failed ...."
      RETURN
    ENDIF
    
    self%n_features = n_features
    self%n_neurons  = n_neurons

    ! Set allocation indicator
    self%Is_Allocated = .TRUE.
  END SUBROUTINE ANN_MLPCoeff_Create

  SUBROUTINE ANN_MLP_Recon(clf, X, Y, ivar)
    TYPE(ANN_MLPCoeff_type), INTENT(IN) :: clf
    TYPE(iVar_Type),      INTENT(INOUT) :: ivar
    REAL(fp), INTENT(IN)  :: X(:)
    REAL(fp), INTENT(OUT) :: Y 

    INTEGER  :: i, j, Error_Status

    IF(ivar%Is_Allocated) THEN
      IF(ivar%n_neurons /= clf%n_neurons) THEN
        WRITE(*,*)"Warning, deallocating ANN_MLP_Recon ivar%neurons ..."
        DEALLOCATE(ivar%neurons, STAT = Error_Status)
        ivar%Is_Allocated = .FALSE.
      ENDIF
    ENDIF

    IF(.NOT. ivar%Is_Allocated) THEN
      ALLOCATE(ivar%neurons(clf%n_neurons), STAT = Error_Status)
      IF(Error_Status /= 0)THEN
        WRITE(*,*)"ANN_MLP_Recon memory allocation failed ..."
        RETURN
      ENDIF
      ivar%n_neurons = clf%n_neurons
      ivar%Is_Allocated = .TRUE. 
    ENDIF

    ! summation over neuraons
    Y = 0.0
    DO i = 1, clf%n_neurons
       ivar%neurons(i) = 0.0
       DO j = 1, clf%n_features
          ivar%neurons(i) = ivar%neurons(i) + clf%w1(i,j)*X(j)
       ENDDO
       ivar%neurons(i) = ivar%neurons(i) + clf%b1(i)
       Y = Y + tanh(ivar%neurons(i))*clf%w2(i)
    ENDDO
    Y = Y + clf%b2

  END SUBROUTINE ANN_MLP_Recon

  SUBROUTINE ANN_MLP_Recon_TL(clf, X_TL, Y_TL, iVar)
    TYPE(ANN_MLPCoeff_type), INTENT(IN) :: clf
    TYPE(iVar_Type),         INTENT(IN) :: iVar
    REAL(fp), INTENT(IN)  :: X_TL(:)
    REAL(fp), INTENT(OUT) :: Y_TL 
    REAL(fp) :: neuron_TL
    INTEGER  :: i, j

    IF((.NOT. ivar%Is_Allocated) .OR. &
      (ivar%n_neurons /= clf%n_neurons)) THEN
      WRITE(*,*)"ANN_MLP_Recon_TL failed ..."
      RETURN  
    ENDIF
    ! summation over neuraons
    Y_TL = 0.0
    DO i = 1, clf%n_neurons
       neuron_TL = 0.0
       DO j = 1, clf%n_features
          neuron_TL = neuron_TL + clf%w1(i,j)*X_TL(j)
       ENDDO
       Y_TL = Y_TL + (1.0-(tanh(ivar%neurons(i)))**2)*clf%w2(i)*neuron_TL
    ENDDO


  END SUBROUTINE ANN_MLP_Recon_TL
 
  SUBROUTINE ANN_MLP_Recon_AD(clf, X_AD, Y_AD, iVar)
    TYPE(ANN_MLPCoeff_type), INTENT(IN) :: clf
    TYPE(iVar_Type),         INTENT(IN) :: iVar
    REAL(fp), INTENT(INOUT)  :: X_AD(:)
    REAL(fp), INTENT(INOUT)  :: Y_AD
    REAL(fp) :: neuron_AD
    INTEGER  :: i, j

    IF((.NOT. ivar%Is_Allocated) .OR. &
      (ivar%n_neurons /= clf%n_neurons)) THEN
      WRITE(*,*)"ANN_MLP_Recon_AD failed ..."
      RETURN  
    ENDIF

    ! summation over neuraons
    DO i = clf%n_neurons, 1, -1
      neuron_AD = (1.0-(tanh(ivar%neurons(i)))**2)*clf%w2(i)*Y_AD
       DO j = clf%n_features, 1, -1
          X_AD(j) = X_AD(j) + clf%w1(i,j)*neuron_AD
       ENDDO
    ENDDO
    Y_AD = 0.0
   
  END SUBROUTINE ANN_MLP_Recon_AD



END MODULE ANN_MLPCoeff_Module
