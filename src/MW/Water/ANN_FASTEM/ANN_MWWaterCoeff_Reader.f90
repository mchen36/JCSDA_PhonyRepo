!
! ANN_MWwaterCoeff_Reader
!
! Module containing the shared CSEM microwave water surface emissivity
! data and their load/destruction routines. 
!
! PUBLIC DATA:
!   MWwaterC:  Data structure containing the microwave water surface
!              emissivity data.
!
! SIDE EFFECTS:
!       Routines in this module modify the contents of the public
!       data structure MWwaterC.
!
! RESTRICTIONS:
!       Routines in this module should only be called during the
!       CSEM initialisation.
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 08-28-2017
!                       ming.chen@noaa.gov
!

MODULE ANN_MWwaterCoeff_Reader
  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use
  USE CSEM_Type_Kinds,          ONLY: fp=>CSEM_fp
  USE CSEM_Exception_Handler
  USE ANN_MLPCoeff_Module,      ONLY: ANN_MLPCoeff_Type,  &
                                ANN_MLPCoeff_Create,      &
                                ANN_MLPCoeff_Destroy

  USE netcdf
  ! Disable implicit typing
  IMPLICIT NONE
 
  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE

  ! Procedures
  PUBLIC :: ANN_MWwaterCoeff_Load
  PUBLIC :: ANN_MWwaterCoeff_CleanUp
  
  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: ANN_MWwaterCoeff.f90 21141 2017-09-18 17:40:43Z ming.chen@noaa.gov $'
  ! Message string length
  INTEGER, PARAMETER :: ML = 512
  
  
  ! maximum variable dimensions
  INTEGER,      PARAMETER :: MAX_VAR_DIMS = 4
  
  ! Global attribute names. Case sensitive
  CHARACTER(*), PARAMETER :: SCALING_NAME       = 'scaling' 
  CHARACTER(*), PARAMETER :: DESCRIPTION_NAME   = 'description'

  ! Dimension names
  CHARACTER(*), PARAMETER :: FEATURE_DIMNAME    = 'n_features'
  CHARACTER(*), PARAMETER :: NEURONS_DIMNAME    = 'n_neurons'

  ! Variable names. Case sensitive.
  CHARACTER(*), PARAMETER :: WEIGHT0_VARNAME    = 'coef0'
  CHARACTER(*), PARAMETER :: WEIGHT1_VARNAME    = 'coef1'
  CHARACTER(*), PARAMETER :: INTERCEPT0_VARNAME = 'intercept0'
  CHARACTER(*), PARAMETER :: INTERCEPT1_VARNAME = 'intercept1'

  ! Variable description attribute.


  ! Variable _FillValue attribute.
  CHARACTER(*), PARAMETER :: FILLVALUE_ATTNAME = '_FillValue'

  

CONTAINS


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       ANN_MWwaterCoeff_Load
!
! PURPOSE:
!       Function to load the microwave water surface emissivity data into
!       the public data structure MWwaterC
!
! CALLING SEQUENCE:
!       Error_Status = ANN_MWwaterCoeff_Load( &
!                        Filename,                              &
!                        File_Path         = File_Path        , &
!                        Quiet             = Quiet            , &
!                        Process_ID        = Process_ID       , &
!                        Output_Process_ID = Output_Process_ID  )
!
! INPUT ARGUMENTS:
!       Filename:           Name of the MWwaterCoeff file.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!
! OPTIONAL INPUT ARGUMENTS:
!       File_Path:          Character string specifying a file path for the
!                           input data file. If not specified, the current
!                           directory is the default.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Quiet:              Set this logical argument to suppress INFORMATION
!                           messages being printed to stdout
!                           If == .FALSE., INFORMATION messages are OUTPUT [DEFAULT].
!                              == .TRUE.,  INFORMATION messages are SUPPRESSED.
!                           If not specified, default is .FALSE.
!                           UNITS:      N/A
!                           TYPE:       LOGICAL
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Process_ID:         Set this argument to the MPI process ID that this
!                           function call is running under. This value is used
!                           solely for controlling INFORMATIOn message output.
!                           If MPI is not being used, ignore this argument.
!                           This argument is ignored if the Quiet argument is set.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Output_Process_ID:  Set this argument to the MPI process ID in which
!                           all INFORMATION messages are to be output. If
!                           the passed Process_ID value agrees with this value
!                           the INFORMATION messages are output. 
!                           This argument is ignored if the Quiet argument
!                           is set.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:       The return value is an integer defining the error
!                           status. The error codes are defined in the
!                           Message_Handler module.
!                           If == SUCCESS the data load was successful
!                              == FAILURE an unrecoverable error occurred.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data
!       structure MWwaterC.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION ANN_MWwaterCoeff_Load( &
    Filename         , &  ! Input
    NN_MWwaterC      , &  ! Input
    SUB_GROUP_NAME   , &  ! Input
    N_SUBGRP         , &  ! Input
    File_Path        , &  ! Optional input
    Quiet            , &  ! Optional input
    Process_ID       , &  ! Optional input
    Output_Process_ID) &  ! Optional input
  RESULT( Error_Status )
    ! Arguments
    CHARACTER(*),           INTENT(IN) :: Filename
    INTEGER,                INTENT(IN) :: N_SUBGRP
    CHARACTER(LEN=15),      INTENT(IN) :: SUB_GROUP_NAME(N_SUBGRP) 
    
    CHARACTER(*), OPTIONAL, INTENT(IN) :: File_Path
    LOGICAL     , OPTIONAL, INTENT(IN) :: Quiet             
    INTEGER     , OPTIONAL, INTENT(IN) :: Process_ID
    INTEGER     , OPTIONAL, INTENT(IN) :: Output_Process_ID

    TYPE(ANN_MLPCoeff_Type),INTENT(INOUT) :: NN_MWwaterC(N_SUBGRP)
    
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ANN_MWwaterCoeff_Load'
    ! Local variables
    CHARACTER(ML) :: Message, pid_msg
    CHARACTER(ML) :: MWwaterCoeff_File
    LOGICAL :: noisy

    CHARACTER (LEN = 8) :: scale
   
    INTEGER :: GIDS(N_SUBGRP+1)
    INTEGER :: ndim_subgrp(N_SUBGRP)
    INTEGER :: dims_subgrp(MAX_VAR_DIMS,N_SUBGRP)
    INTEGER :: varid, i, dim1(1),dim2(2)
   
    ! Setup 
    Error_Status = SUCCESS
    ! ...Assign the filename to local variable
    MWwaterCoeff_File = ADJUSTL(Filename)
    ! ...Add the file path
    IF ( PRESENT(File_Path) ) MWwaterCoeff_File = TRIM(ADJUSTL(File_Path))//TRIM(MWwaterCoeff_File)
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Check the MPI Process Ids
    IF ( noisy .AND. PRESENT(Process_ID) .AND. PRESENT(Output_Process_ID) ) THEN
      IF ( Process_Id /= Output_Process_Id ) noisy = .FALSE.
    END IF
    ! ...Create a process ID message tag for error messages
    IF ( PRESENT(Process_Id) ) THEN
      WRITE( pid_msg,'("; Process ID: ",i0)' ) Process_ID
    ELSE
      pid_msg = ''
    END IF
  
  
    Error_Status = INQ_ANNCoeff_File( &
          MWwaterCoeff_File,     &
          SUB_GROUP_NAME,        &  ! Input
          N_SUBGRP,              &  ! Input
          ndim_subgrp,           &
          dims_subgrp,           &
          GIDS)  

    ! Allocate the output structure
    DO i = 1, N_SUBGRP
       CALL  ANN_MLPCoeff_Create(NN_MWwaterC(i),dims_subgrp(1,i),dims_subgrp(2,i))
    ENDDO

    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error occurred allocating MWwaterC structure.'
      PRINT*,Message
      Error_Status = check( nf90_close(GIDS(1)) )
      RETURN
    END IF

    ! Read the MWwaterCoeff data file
    DO i = 1, N_SUBGRP
       dim1 = dims_subgrp(1:1,i) ; dim2 = dims_subgrp(1:2,i)
       Error_Status = check(nf90_inq_varid(GIDS(i+1), 'coef0', varid))
       Error_Status = check(nf90_get_var(GIDS(i+1), varid,  &
              NN_MWwaterC(i)%w1,  count=dim2))
       Error_Status = check(nf90_inq_varid(GIDS(i+1), 'coef1', varid))
       Error_Status = check(nf90_get_var(GIDS(i+1), varid,  &
              NN_MWwaterC(i)%w2,  count=dim1))

       Error_Status = check(nf90_inq_varid(GIDS(i+1), 'intercepts0', varid))
       Error_Status = check(nf90_get_var(GIDS(i+1), varid,  &
              NN_MWwaterC(i)%b1,  count=dim1))
       Error_Status = check(nf90_inq_varid(GIDS(i+1), 'intercepts1', varid))
       Error_Status = check(nf90_get_var(GIDS(i+1), varid,  &
              NN_MWwaterC(i)%b2))

       Error_Status = check(nf90_inq_varid(GIDS(i+1), 'vmin', varid))
       Error_Status = check(nf90_get_var(GIDS(i+1), varid,  &
              NN_MWwaterC(i)%vmin))
       Error_Status = check(nf90_inq_varid(GIDS(i+1), 'vmax', varid))
       Error_Status = check(nf90_get_var(GIDS(i+1), varid,  &
              NN_MWwaterC(i)%vmax))

       Error_Status = check(nf90_get_att( GIDS(i+1),nf90_global, "scaling", scale))
       READ(scale,"(f8.0)")NN_MWwaterC(i)%fscale
        !NN_MWwaterC(i)%w1(:,:) =     NN_MWwaterC(i)%w1(:,(/3,2,1/))
        !NN_MWwaterC(i)%vmax = NN_MWwaterC(i)%vmax((/3,2,1/))
        !NN_MWwaterC(i)%vmin  = NN_MWwaterC(i)%vmin((/3,2,1/))
    ENDDO

    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error occurred allocating MWwaterC structure.'
      PRINT*,Message
      Error_Status = check(nf90_close(GIDS(1)) )
      RETURN
    END IF
   
 
 
  END FUNCTION ANN_MWwaterCoeff_Load

 
!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       ANN_MWwaterCoeff_Destroy
!
! PURPOSE:
!       Function to deallocate the public data structure MWwaterC containing
!       the CSEM microwave water surface emissivity data.
!
! CALLING SEQUENCE:
!       Error_Status = ANN_MWwaterCoeff_Destroy( Process_ID = Process_ID )
!
! OPTIONAL INPUTS:
!       Process_ID:       Set this argument to the MPI process ID that this
!                         function call is running under. This value is used
!                         solely for controlling message output. If MPI is not
!                         being used, ignore this argument.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:     The return value is an integer defining the error
!                         status. The error codes are defined in the
!                         Message_Handler module.
!                         If == SUCCESS the deallocation of the public data
!                                       structure was successful
!                            == FAILURE an unrecoverable error occurred.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data
!       structure NN_MWwaterC.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION ANN_MWwaterCoeff_CleanUp( NN_MWwaterC, Process_ID ) RESULT( err_stat )
    ! Arguments
    TYPE(ANN_MLPCoeff_Type),INTENT(INOUT) :: NN_MWwaterC(:)
    INTEGER, OPTIONAL, INTENT(IN) :: Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ANN_MWwaterCoeff_Destroy'
    ! Local variables
    CHARACTER(ML) :: pid_msg
    INTEGER :: i

    ! Setup
    err_stat = SUCCESS
    ! ...Create a process ID message tag for error messages
    IF ( PRESENT(Process_Id) ) THEN
      WRITE( pid_msg,'("; Process ID: ",i0)' ) Process_ID
    ELSE
      pid_msg = ''
    END IF
    
 
    ! Destroy the structure
    DO i = 1, size(NN_MWwaterC)
       CALL ANN_MLPCoeff_Destroy(NN_MWwaterC(i))
    END DO

  END FUNCTION ANN_MWwaterCoeff_CleanUp

  
!##################################################################################
!##################################################################################
!##                                                                              ##
!##                          ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################


  FUNCTION INQ_ANNCoeff_File( &
       FILE_NAME,     &
       SUB_GROUP_NAME,&  ! Input
       N_SUBGRP,      &
       ndim_subgrp,   &
       dims_subgrp,   &
       GIDS)          & 
    RESULT ( Error_status)
  
    CHARACTER(*), INTENT(IN)  :: FILE_NAME
    INTEGER,      INTENT(IN)  :: N_SUBGRP
    INTEGER,      INTENT(OUT) :: ndim_subgrp(N_SUBGRP)
    INTEGER,      INTENT(OUT) :: dims_subgrp(MAX_VAR_DIMS,N_SUBGRP)
    INTEGER,      INTENT(OUT),   OPTIONAL  :: GIDS(N_SUBGRP+1)
 
    CHARACTER(LEN=15),      INTENT(IN) :: SUB_GROUP_NAME(N_SUBGRP) 

    INTEGER :: ncid,  varid
    LOGICAL :: Existence
    INTEGER :: i, j, len, ndims, numgrps
    INTEGER :: ncids(N_SUBGRP+1)
    INTEGER :: grpids(N_SUBGRP+1)
    INTEGER :: dimids(MAX_VAR_DIMS)
    INTEGER :: Error_status

  
    Error_Status = SUCCESS
    INQUIRE( FILE = TRIM( FILE_NAME), EXIST = Existence )
    IF ( .NOT. Existence ) THEN
      PRINT*,'File '//TRIM( FILE_NAME )//' not found.'
      Error_Status = FAILURE
      RETURN
    END IF

    !Open the NetCDF file:
    Error_Status = check(nf90_open( TRIM( FILE_NAME ), nf90_nowrite, ncid))
    Error_Status = check(nf90_inq_grps(ncid, numgrps, ncids))
    grpids(1) = ncid

    IF(numgrps /= N_SUBGRP) THEN
       PRINT*, 'Number of sub-groups is not equal to the defined '
       Error_Status = check( nf90_close(ncid) ) 
       PRINT*, 'NC file does not have the target dataset group.....' 
       RETURN    
    END IF
    DO i = 1, N_SUBGRP
        Error_Status = check(nf90_inq_grp_ncid(grpids(1), TRIM(SUB_GROUP_NAME(i)), grpids(i+1)))
        IF(Error_Status /= SUCCESS) THEN 
          Error_Status = check( nf90_close(ncid) ) 
          print*, TRIM(SUB_GROUP_NAME(i)) // '  is not the sub-group name...'
          RETURN
        END IF 
        !Error_Status = check(nf90_inq_dimid(grpids(i+1), FEATURE_DIMNAME,  dimid))
        !Error_Status = check(nf90_inquire_dimension(grpids(i+1), dimid, len=len))
        Error_Status = check(nf90_inq_varid(grpids(i+1), WEIGHT0_VARNAME, varid))
        Error_Status = check(nf90_inquire_variable(grpids(i+1), varid,  ndims=ndims, dimids=dimids))
        ndim_subgrp(i) = ndims
        DO j = 1, ndims
          Error_Status = check(nf90_inquire_dimension(grpids(i+1), dimids(j), len=len))
          dims_subgrp(j,i) = len
        END DO

    END DO
 
    IF(PRESENT(GIDS)) THEN
      GIDS =  grpids 
      RETURN
    ENDIF
    
    Error_Status = check( nf90_close(ncid) )

  END FUNCTION  INQ_ANNCoeff_File


  
  FUNCTION check(status) RESULT ( Err_Status )

    INTEGER, INTENT ( IN) :: status
    INTEGER :: err_status
    err_status = SUCCESS
    IF (status /= nf90_noerr) THEN
      PRINT *, trim(nf90_strerror(status))
      err_status = FAILURE
    END IF
  END FUNCTION check


END MODULE ANN_MWwaterCoeff_Reader
