! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 1.
!
! Description:
!   This routine defines the wavelengths of the spectral bands
!
!- ----------------------------------------------------------------------------
SUBROUTINE make_block_1(Spectrum, l_interactive, ierr)

  USE realtype_rd
  USE def_spectrum
  USE def_std_io_icf
  USE error_pcf

  IMPLICIT none


! Dummy arguments.
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive operation
  TYPE(StrSpecData), Intent(INOUT), Target :: Spectrum

! Local variables.
  CHARACTER (LEN=1) :: char_unit
!   Unit type
  INTEGER :: ios
!   Reading error flag
  INTEGER :: i, j
!   Loop variable
  INTEGER :: range_bands
!   Number of bands in range
  LOGICAL :: l_metre = .FALSE.
!   Logical for units of metres
  LOGICAL :: l_inverse_cm = .FALSE.
!   Logical for units of inverse cm
  LOGICAL :: l_micron = .FALSE.
!   Logical for units of microns
  LOGICAL :: l_range = .FALSE.
!   Logical for specifying bands as a range
  REAL (RealK) :: range_low, range_high, range_inc
!   Range lower and upper limits and increment
  REAL (RealK) :: wavelength_temp
!   Temporary for interchange
  TYPE (StrSpecBasic), Pointer :: SpBasic
!   Pointer to basic components of the spectrum
  TYPE (StrSpecDim), Pointer :: SpDim
!   Pointer to dimensions within the spectrum


  SpDim   => Spectrum%Dim
  ALLOCATE(Spectrum%Basic%wavelength_long(SpDim%nd_band))
  ALLOCATE(Spectrum%Basic%wavelength_short(SpDim%nd_band))
  SpBasic => Spectrum%Basic

! Obtain the limits for each band.
  WRITE(iu_stdout, '(4(/A))') &
    'Specify the band limits, type either:', &
    '  r) to define a range of bands of equal size', &
    '  m, c, or u) to set individual band limits in units of', &
    '     metres, inverse centimetres or microns respectively'
  DO
    READ(iu_stdin, '(A)') char_unit
    IF ( (char_unit.eq.'m').OR.(char_unit.eq.'M') ) THEN
      l_metre=.TRUE.
      EXIT
    ELSE IF ( (char_unit.eq.'c').OR.(char_unit.eq.'C') ) THEN
      l_inverse_cm=.TRUE.
      EXIT
    ELSE IF ( (char_unit.eq.'u').OR.(char_unit.eq.'U') ) THEN
      l_micron=.TRUE.
      EXIT
    ELSE IF ( (char_unit.eq.'r').OR.(char_unit.eq.'R') ) THEN
      l_range=.TRUE.
      EXIT
    ELSE
      WRITE(iu_err, '(a)') 'Unknown reponse:'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(a)') 'Please re-enter.'
      ELSE
        ierr=i_err_fatal
        RETURN
      ENDIF
    ENDIF
  ENDDO

  IF (l_range) THEN
    i = 0
    DO
      WRITE(iu_stdout, '(A)') &
        'Enter units followed by lower and upper limits and increment:'
      DO
        READ(iu_stdin, *, IOSTAT=ios) &
          char_unit, range_low, range_high, range_inc
        IF ( (char_unit.eq.'m').OR.(char_unit.eq.'M') ) THEN
          l_metre=.TRUE.
          l_inverse_cm=.FALSE.
          l_micron=.FALSE.
          EXIT
        ELSE IF ( (char_unit.eq.'c').OR.(char_unit.eq.'C') ) THEN
          l_inverse_cm=.TRUE.
          l_metre=.FALSE.
          l_micron=.FALSE.
          EXIT
        ELSE IF ( (char_unit.eq.'u').OR.(char_unit.eq.'U') ) THEN
          l_micron=.TRUE.
          l_metre=.FALSE.
          l_inverse_cm=.FALSE.
          EXIT
        ELSE
          WRITE(iu_err, '(a)') 'Unknown reponse:'
          IF (l_interactive) THEN
            WRITE(iu_stdout, '(a)') 'Please re-enter.'
          ELSE
            ierr=i_err_fatal
            RETURN
          END IF
        END IF
      END DO
      range_bands = NINT( (range_high - range_low) / range_inc )
      IF (i + range_bands > SpBasic%n_band) THEN
        WRITE(iu_err, '(a)') 'Specified range exceeds number of bands.'
        CYCLE
      END IF
      DO j = 1, range_bands
        i = i + 1
        SpBasic%wavelength_short(i) = range_low + range_inc*REAL(j-1, RealK)
        SpBasic%wavelength_long(i)  = range_low + range_inc*REAL(j,   RealK)
!       Convert to metres if in inverse cm or microns.
        IF (l_inverse_cm) THEN
          SpBasic%wavelength_short(i) = &
            1.0e-02_RealK / SpBasic%wavelength_short(i)
          SpBasic%wavelength_long(i) = &
            1.0e-02_RealK / SpBasic%wavelength_long(i)
        ELSE IF (l_micron) THEN
          SpBasic%wavelength_short(i) = &
            1.0e-06_RealK * SpBasic%wavelength_short(i)
          SpBasic%wavelength_long(i) = &
            1.0e-06_RealK * SpBasic%wavelength_long(i)
        ENDIF
!       Interchange the limits if they are in the wrong order.
        IF (SpBasic%wavelength_short(i) > SpBasic%wavelength_long(i)) THEN
          wavelength_temp = SpBasic%wavelength_short(i)
          SpBasic%wavelength_short(i) = SpBasic%wavelength_long(i)
          SpBasic%wavelength_long(i) = wavelength_temp
        ENDIF
      END DO
      IF (i == SpBasic%n_band) EXIT
    END DO
  ELSE
    DO i=1, SpBasic%n_band
      WRITE(iu_stdout, '(A21, 1X, I5)') 'Enter limits for band ', i
      DO
        READ(iu_stdin, *, IOSTAT=ios) &
          SpBasic%wavelength_short(i), SpBasic%wavelength_long(i)
        IF (ios == 0) THEN
          EXIT
        ELSE IF (l_interactive) THEN
          WRITE(iu_stdout, '(A)') 'Please re-enter.'
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDDO

!     Convert to metres if in inverse cm or microns.
      IF (l_inverse_cm) THEN
        SpBasic%wavelength_short(i) = &
          1.0e-02_RealK / SpBasic%wavelength_short(i)
        SpBasic%wavelength_long(i) = &
          1.0e-02_RealK / SpBasic%wavelength_long(i)
      ELSE IF (l_micron) THEN
        SpBasic%wavelength_short(i) = &
          1.0e-06_RealK * SpBasic%wavelength_short(i)
        SpBasic%wavelength_long(i) = &
          1.0e-06_RealK * SpBasic%wavelength_long(i)
      ENDIF
!     Interchange the limits if they are in the wrong order.
      IF (SpBasic%wavelength_short(i) > SpBasic%wavelength_long(i)) THEN
        wavelength_temp = SpBasic%wavelength_short(i)
        SpBasic%wavelength_short(i) = SpBasic%wavelength_long(i)
        SpBasic%wavelength_long(i) = wavelength_temp
      ENDIF
    ENDDO
  END IF

  WRITE(iu_stdout, '(A//)') 'All bands specified.'

END SUBROUTINE make_block_1
