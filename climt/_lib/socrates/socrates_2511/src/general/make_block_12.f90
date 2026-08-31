! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 12.
!
! Method:
!       The type of ice crystal is set and a file containing a 
!       parametrization of the single scattering properties is
!       read.
!- ---------------------------------------------------------------------
SUBROUTINE make_block_12(Spectrum, ierr)

  USE realtype_rd
  USE def_spectrum
  USE rad_pcf
  USE file_type_pcf
  USE def_std_io_icf

  IMPLICIT NONE

  TYPE (StrSpecData), Intent(INOUT), TARGET :: Spectrum
!   Spectral file to be assigned
  INTEGER, Intent(INOUT) :: ierr
!   Error flag

! Local variables.
  CHARACTER (LEN=1) :: char_yn
!   Character response variable
  CHARACTER (LEN=36) :: char_36
!   Character response variable
  LOGICAL :: l_range_set
!   Flag to indicate validity range has been set
  INTEGER :: n_parameter
!   Number of parameters
  INTEGER :: ios
!   I/O error flag
  INTEGER :: i_input_type
!   Type of input file
  INTEGER :: iu_cloud_fit
!   Unit number for reading the fit to cloud data
  INTEGER :: i_ice
!   Identification number of ice crystal
  INTEGER :: i, k
!   Loop variables
  REAL (RealK) :: ice_dim_temp
!   Temporary dimension of ice crystal used in ordering sizes


  IF (.NOT.Spectrum%Basic%l_present(12)) THEN
    IF (ALLOCATED(Spectrum%Ice%l_ice_type)) &
        DEALLOCATE(Spectrum%Ice%l_ice_type)
    ALLOCATE(Spectrum%Ice%l_ice_type(Spectrum%Dim%nd_ice_type))
    Spectrum%Ice%l_ice_type = .FALSE.
    IF (ALLOCATED(Spectrum%Ice%i_ice_parm)) &
        DEALLOCATE(Spectrum%Ice%i_ice_parm)
    ALLOCATE(Spectrum%Ice%i_ice_parm(Spectrum%Dim%nd_ice_type))
    IF (ALLOCATED(Spectrum%Ice%parm_min_dim)) &
        DEALLOCATE(Spectrum%Ice%parm_min_dim)
    ALLOCATE(Spectrum%Ice%parm_min_dim(Spectrum%Dim%nd_ice_type))
    IF (ALLOCATED(Spectrum%Ice%parm_max_dim)) &
        DEALLOCATE(Spectrum%Ice%parm_max_dim)
    ALLOCATE(Spectrum%Ice%parm_max_dim(Spectrum%Dim%nd_ice_type))
    IF (ALLOCATED(Spectrum%Ice%n_phf)) &
        DEALLOCATE(Spectrum%Ice%n_phf)
    ALLOCATE(Spectrum%Ice%n_phf(Spectrum%Dim%nd_ice_type))
    IF (ALLOCATED(Spectrum%Ice%parm_list)) &
        DEALLOCATE(Spectrum%Ice%parm_list)
    ALLOCATE(Spectrum%Ice%parm_list(Spectrum%Dim%nd_cloud_parameter, &
      Spectrum%Dim%nd_band, Spectrum%Dim%nd_ice_type))
  END IF

  CALL get_free_unit(ierr, iu_cloud_fit)

  WRITE(iu_stdout, '(/a)') 'Enter the number for the type of ' &
    //'ice crystals for which data are to be provided.'
1 READ(iu_stdin, *, iostat=ios) i_ice
  IF (ios /= 0) THEN
    write (iu_err, '(a)') 'unrecognized response: '
    write (iu_stdout, '(a)') 'please re-enter.'
    goto 1
  ENDIF
  IF ( (i_ice <= 0).OR.(i_ice > Spectrum%Dim%nd_ice_type) ) THEN
    WRITE(iu_err, '(a)') '+++ response is out of range: '
    WRITE(iu_stdout, '(a)') 'please re-enter.'
    goto 1
  ENDIF

  IF (Spectrum%Ice%l_ice_type(i_ice)) THEN
    WRITE(iu_stdout, '(/a)') 'Information for this type ' &
      //'of ice crystal is already present.'
    WRITE(iu_stdout, '(a)') 'Do you wish to overwrite (y/n), or remove? (r)'
2   read(iu_stdin, '(a)') char_yn
    IF ( (char_yn == 'N').OR.(char_yn == 'n') ) THEN
      WRITE(iu_stdout, '(a)') 'Enter a new number.'
      RETURN
    ELSE IF ( (char_yn == 'R').OR.(char_yn == 'r') ) THEN
      Spectrum%Ice%l_ice_type(i_ice)=.FALSE.
      ! Set the presence flag for the data.
      Spectrum%Basic%l_present(12)=.FALSE.
      DO i=1, Spectrum%Dim%nd_ice_type
        Spectrum%Basic%l_present(12) = Spectrum%Basic%l_present(12) &
          .OR. Spectrum%Ice%l_ice_type(i)
      END DO
      RETURN
    ELSE IF ( (char_yn /= 'Y').AND.(char_yn /= 'y') ) THEN
      WRITE(iu_err, '(a)') '+++ Unrecognized response: '
      WRITE(iu_stdout, '(a)') 'Please re-enter.'
      goto 2
    ENDIF
  ENDIF


! Open the file of input data.
8 CALL open_file_in(ierr, iu_cloud_fit, &
    'Enter the name of the file containing the fitted ' &
    //'parameters for ice crystals.')
  IF (ierr /= i_normal) RETURN

  READ(iu_cloud_fit, '(13x, i5)') i_input_type

  IF ( (i_input_type /= IT_file_cloud_fit).AND. &
       (i_input_type /= IT_file_cloud_fit_phf) ) THEN
    WRITE(iu_err, '(a)') '+++ This file is of an invalid type.'
    goto 8
  ENDIF

  READ(iu_cloud_fit, '(//, a)', iostat=ios) char_36
  l_range_set = .FALSE.
  IF (char_36 == 'Min and max characteristic dimension') THEN
    BACKSPACE(iu_cloud_fit)
    READ(iu_cloud_fit, '(38x, 2(f12.5))', iostat=ios) &
      Spectrum%Ice%parm_min_dim(i_ice), Spectrum%Ice%parm_max_dim(i_ice)
    IF (ios == 0) THEN
      l_range_set = .TRUE.
    ELSE
      ! Recover from error and ask for range to be set explicitly
      ios = 0
    END IF
  END IF
  READ(iu_cloud_fit, '(33x, i5, /, 31x, i5)', iostat=ios) &
    Spectrum%Ice%i_ice_parm(i_ice), n_parameter

  IF (i_input_type == IT_file_cloud_fit_phf) THEN
    READ(iu_cloud_fit, '(31x, i3)') Spectrum%Ice%n_phf(i_ice)
  ELSE
    Spectrum%Ice%n_phf(i_ice)=1
  ENDIF

  IF (ios /= 0) THEN
    WRITE(iu_err, '(/a)') '*** Error: The format of the data is invalid.'
  ENDIF

  DO i=1, Spectrum%Basic%n_band
    READ(iu_cloud_fit, '(/, 4(4x, 1pe12.5))', iostat=ios) &
      (Spectrum%Ice%parm_list(k, i, i_ice), k=1, n_parameter)
    IF (ios /= 0) THEN
      WRITE(iu_err, '(/a)') &
        '*** Error: The format of the data is invalid.'
    ENDIF
  ENDDO

  CLOSE(iu_cloud_fit)

  IF (.NOT. l_range_set) THEN
    WRITE(iu_stdout, '(/a)') &
      'Enter the range of validity of the parametrization.'
    READ(iu_stdin, *) Spectrum%Ice%parm_min_dim(i_ice), &
                      Spectrum%Ice%parm_max_dim(i_ice)
  END IF
  IF (Spectrum%Ice%parm_min_dim(i_ice) > &
      Spectrum%Ice%parm_max_dim(i_ice)) THEN
    ice_dim_temp=Spectrum%Ice%parm_max_dim(i_ice)
    Spectrum%Ice%parm_max_dim(i_ice)=Spectrum%Ice%parm_min_dim(i_ice)
    Spectrum%Ice%parm_min_dim(i_ice)=ice_dim_temp
  ENDIF

  Spectrum%Ice%l_ice_type(i_ice)=.TRUE.

! Set the presence flag for the data.
  Spectrum%Basic%l_present(12)=.FALSE.
  DO i=1, Spectrum%Dim%nd_ice_type
    Spectrum%Basic%l_present(12) = Spectrum%Basic%l_present(12) &
      .OR. Spectrum%Ice%l_ice_type(i)
  ENDDO

END SUBROUTINE make_block_12
