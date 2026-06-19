! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a file of refractive indices.
!
Subroutine get_refract_index &
!
(l_interactive, nd_refract, n_refract, wavelength_refract, &
 re_refract, im_refract, &
 l_spline_interp, d2_re_refract, d2_im_refract, &
 ierr)
!
! Description:
!   This subroutine reads a file containing the 
!   complex refractive indices for the material
!   of the scatterer. A flag for interpolation
!   od the refractive indices is set.
!
! Modules used:
  USE def_std_io_icf
  USE error_pcf
  USE realtype_rd
  USE spline_fit_mod, ONLY: spline_fit
!
  IMPLICIT NONE
!
!
  LOGICAL, Intent(IN)    :: l_interactive
!                             Flag for interactive operation
  INTEGER, Intent(IN)    :: nd_refract
!                             Size allocated for refractive indices
  INTEGER, Intent(INOUT) :: ierr
!                             Error flag
  INTEGER, Intent(OUT)   :: n_refract
!                             Number of wavlengths for refractive indices
!
  REAL (RealK), Intent(OUT) :: wavelength_refract(nd_refract)
!                                Wavelengths for refractive indices 
  REAL (RealK), Intent(OUT) :: re_refract(nd_refract)
!                                Real parts of the refractive indices 
  REAL (RealK), Intent(OUT) :: im_refract(nd_refract)
!                                Imaginary parts of the refractive 
!                                indices 
!
  LOGICAL, Intent(OUT)   :: l_spline_interp
!                             Flag for interpolation using cubic splines
  REAL (RealK), Intent(OUT) :: d2_re_refract(nd_refract)
!                                Second derivatives of the real parts 
!                                of the refractive indices 
  REAL (RealK), Intent(OUT) :: d2_im_refract(nd_refract)
!                                Second derivatives of the imaginary
!                                parts of the refractive indices 
!
!
! Local variables:
  CHARACTER (LEN=80) :: line
!                         Line of input text
  CHARACTER (LEN=1)  :: char_yn
!                         Single character input
!
  INTEGER :: iu_input
!              Unit number for input file
  INTEGER :: ios
!              I/O status of input file
  LOGICAL :: l_data_region
!              Flag for region of data in the file
  LOGICAL :: l_exit
!              Flag for exit from loop determining interpolation
!              of the refractive indices.
!
!
!- End of header
!
!
!
  CALL get_free_unit(ierr, iu_input)
  IF ( ierr /= i_normal ) RETURN
  CALL open_file_in(ierr, iu_input, &
         'Enter the name of the file of refractive indices.')
  IF ( ierr /= i_normal ) RETURN
!
  l_data_region = .FALSE.
  n_refract     = 0
  ios           = 0
!
  DO ; IF ( ios /= 0 ) EXIT
!
    read(iu_input, '(A)', iostat=ios) line
!
    IF ( l_data_region ) THEN
!
      IF ( line(1:4) .NE. '*END' ) THEN
        BACKSPACE(iu_input)
        n_refract = n_refract + 1
        IF ( n_refract > nd_refract ) THEN
          WRITE(iu_err, '(/A)') '*** Error: Refractive ' &
                //'indices are specified at too many wavelengths.'
          ierr = i_err_fatal
          RETURN
        ENDIF
        read(iu_input, *) wavelength_refract(n_refract), &
          re_refract(n_refract), im_refract(n_refract)
      ELSE
        l_data_region = .FALSE.
      ENDIF
!
    ELSE
!
      IF ( line(1:11) == '*BEGIN_DATA' ) THEN
        l_data_region = .TRUE.
      ENDIF
!
    ENDIF
!
  ENDDO
!
  CLOSE(iu_input)
!
!
! Choose between linear interpolation of refractive indices
! to the scattering wavelengths or cubic splines.
  WRITE(iu_stdout, '(/A)') &
    'How should the refractive index be interpolated?'
  WRITE(iu_stdout, '(/A)') &
    'Enter C for cubic splines or L for linear interpolation.'
  l_exit=.FALSE.
  DO ; IF (l_exit) EXIT
!
    READ(iu_stdin, '(A)') char_yn
    IF ( (char_yn.EQ.'C').OR.(char_yn.EQ.'c') ) THEN
      l_spline_interp = .TRUE.
      l_exit          = .TRUE.
!     Calculate second derivatives of the refractive indices
!     for interpolation by splines
      CALL spline_fit(n_refract, wavelength_refract, &
        re_refract, d2_re_refract)
      CALL spline_fit(n_refract, wavelength_refract, &
        im_refract, d2_im_refract)
    ELSE IF ( (char_yn.EQ.'L').OR.(char_yn.EQ.'l') ) THEN
      l_spline_interp = .FALSE.
      l_exit          = .TRUE.
!     Zero the second derivatives of the refractive indices
!     to obtain linear fits from generic splining routines.
      d2_re_refract(:) = 0.0E+00_RealK
      d2_im_refract(:) = 0.0E+00_RealK
    ELSE IF (l_interactive) THEN
      WRITE(IU_ERR, '(/A/)') &
        '+++ Unrecognized response: please re-enter.'
    ELSE
      l_exit = .TRUE.
    ENDIF
!
  ENDDO
!
!
!
  RETURN
!
END
