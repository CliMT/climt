! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 10.
!
! Method:
!       The type of droplet is set. A fit is read from a file.
!
!- ---------------------------------------------------------------------
SUBROUTINE make_block_10(Spectrum, ierr)

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
  INTEGER :: n_parameter
!   Number of parameters
  INTEGER :: ios
!   I/O error flag
  INTEGER :: i_input_type
!   Type of input file
  INTEGER :: iu_cloud_fit
!   Unit number for reading the fit to cloud data
  INTEGER :: i_drop
!   Identification number of drops
  INTEGER :: i, k
!   Loop variables
  REAL (RealK) :: drop_dim_temp
!   Temporary dimension of droplet used in ordering sizes


  IF (.NOT.Spectrum%Basic%l_present(10)) THEN
    IF (ALLOCATED(Spectrum%Drop%l_drop_type)) &
        DEALLOCATE(Spectrum%Drop%l_drop_type)
    ALLOCATE(Spectrum%Drop%l_drop_type(Spectrum%Dim%nd_drop_type))
    Spectrum%Drop%l_drop_type = .FALSE.
    IF (ALLOCATED(Spectrum%Drop%i_drop_parm)) &
        DEALLOCATE(Spectrum%Drop%i_drop_parm)
    ALLOCATE(Spectrum%Drop%i_drop_parm(Spectrum%Dim%nd_drop_type))
    IF (ALLOCATED(Spectrum%Drop%parm_min_dim)) &
        DEALLOCATE(Spectrum%Drop%parm_min_dim)
    ALLOCATE(Spectrum%Drop%parm_min_dim(Spectrum%Dim%nd_drop_type))
    IF (ALLOCATED(Spectrum%Drop%parm_max_dim)) &
        DEALLOCATE(Spectrum%Drop%parm_max_dim)
    ALLOCATE(Spectrum%Drop%parm_max_dim(Spectrum%Dim%nd_drop_type))
    IF (ALLOCATED(Spectrum%Drop%n_phf)) &
        DEALLOCATE(Spectrum%Drop%n_phf)
    ALLOCATE(Spectrum%Drop%n_phf(Spectrum%Dim%nd_drop_type))
    IF (ALLOCATED(Spectrum%Drop%parm_list)) &
        DEALLOCATE(Spectrum%Drop%parm_list)
    ALLOCATE(Spectrum%Drop%parm_list(Spectrum%Dim%nd_cloud_parameter, &
      Spectrum%Dim%nd_band, Spectrum%Dim%nd_drop_type))
  END IF

  CALL get_free_unit(ierr, iu_cloud_fit)

  WRITE(iu_stdout, '(/a)') 'enter the number for the type of ' &
    //'droplets for which data are to be provided.'
1 read(iu_stdin, *, iostat=ios) i_drop
  IF (ios /= 0) THEN
    write (iu_err, '(a)') 'unrecognized response: '
    write (iu_stdout, '(a)') 'please re-enter.'
    goto 1
  ENDIF
  IF ( (i_drop <= 0).OR.(i_drop > Spectrum%Dim%nd_drop_type) ) THEN
    WRITE(iu_err, '(a)') '+++ response is out of range: '
    WRITE(iu_stdout, '(a)') 'please re-enter.'
    goto 1
  ENDIF

  IF (Spectrum%Drop%l_drop_type(i_drop)) THEN
    WRITE(iu_stdout, '(/a)') 'Information for this type ' &
      //'of droplet is already present.'
    WRITE(iu_stdout, '(a)') 'Do you wish to overwrite (y/n), or remove? (r)'
2   read(iu_stdin, '(a)') char_yn
    IF ( (char_yn == 'N').OR.(char_yn == 'n') ) THEN
      WRITE(iu_stdout, '(a)') 'Enter a new number.'
      RETURN
    ELSE IF ( (char_yn == 'R').OR.(char_yn == 'r') ) THEN
      Spectrum%Drop%l_drop_type(i_drop)=.FALSE.
      ! Set the presence flag for the data.
      Spectrum%Basic%l_present(10)=.FALSE.
      DO i=1, Spectrum%Dim%nd_drop_type
        Spectrum%Basic%l_present(10) = Spectrum%Basic%l_present(10) &
          .OR. Spectrum%Drop%l_drop_type(i)
      ENDDO
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
    //'parameters for the droplets.')
  IF (ierr /= i_normal) RETURN

  READ(iu_cloud_fit, '(13x, i5)') i_input_type

  IF ( (i_input_type /= IT_file_cloud_fit).AND. &
       (i_input_type /= IT_file_cloud_fit_phf) ) THEN
    WRITE(iu_err, '(a)') '+++ This file is of an invalid type.'
    goto 8
  ENDIF

  READ(iu_cloud_fit, '(///, 33x, i5, /, 31x, i5)', iostat=ios) &
    Spectrum%Drop%i_drop_parm(i_drop), n_parameter

  IF (i_input_type == IT_file_cloud_fit_phf) THEN
    READ(iu_cloud_fit, '(31x, i3)') Spectrum%Drop%n_phf(i_drop)
  ELSE
    Spectrum%Drop%n_phf(i_drop)=1
  ENDIF

  IF (ios /= 0) THEN
    WRITE(iu_err, '(/a)') '*** Error: Invalid format in the input data.'
    ierr=i_err_fatal
    RETURN
  ENDIF

  DO i=1, Spectrum%Basic%n_band
    READ(iu_cloud_fit, '(/, 4(4x, 1pe12.5))', iostat=ios) &
      (Spectrum%Drop%parm_list(k, i, i_drop), k=1, n_parameter)
    IF (ios /= 0) THEN
      WRITE(iu_err, '(/a)') &
        '*** Error: Invalid format in the input data.'
    ENDIF
  ENDDO

  CLOSE(iu_cloud_fit)

  WRITE(iu_stdout, '(/a)') &
    'Enter the range of validity of the parametrization.'
  READ(iu_stdin, *) Spectrum%Drop%parm_min_dim(i_drop), &
                    Spectrum%Drop%parm_max_dim(i_drop)
  IF (Spectrum%Drop%parm_min_dim(i_drop) > &
      Spectrum%Drop%parm_max_dim(i_drop)) THEN
    drop_dim_temp=Spectrum%Drop%parm_max_dim(i_drop)
    Spectrum%Drop%parm_max_dim(i_drop)=Spectrum%Drop%parm_min_dim(i_drop)
    Spectrum%Drop%parm_min_dim(i_drop)=drop_dim_temp
  ENDIF

  Spectrum%Drop%l_drop_type(i_drop)=.TRUE.

! Set the presence flag for the data.
  Spectrum%Basic%l_present(10)=.FALSE.
  DO i=1, Spectrum%Dim%nd_drop_type
    Spectrum%Basic%l_present(10) = Spectrum%Basic%l_present(10) &
      .OR. Spectrum%Drop%l_drop_type(i)
  ENDDO

END SUBROUTINE make_block_10
