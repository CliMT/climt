! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a file of wavelengths for scattering.
!
Subroutine get_wavelengths &
!
(nd_wavelength_scat, n_wavelength, wavelength, ierr)
!
! Description:
!   This subroutine reads a file of wavelengths for
!   calculations of single scattering calculations.
!
! Modules used:
  USE def_std_io_icf
  USE error_pcf
  USE realtype_rd
!
  IMPLICIT NONE
!
!
  INTEGER, Intent(IN)    :: nd_wavelength_scat
!                             Size allocated for wavelengths
  INTEGER, Intent(INOUT) :: ierr
!                             Error flag
  INTEGER, Intent(OUT)   :: n_wavelength
!                             Number of wavlengths read
!
  REAL (RealK), Intent(OUT) :: wavelength(nd_wavelength_scat)
!                                Wavelengths for scattering calculations
!
!
! Local variables:
  CHARACTER (LEN=80) :: line
!                         Line of input text
!
  INTEGER :: iu_input
!              Unit number for input file
  INTEGER :: ios
!              I/O status of input file
  LOGICAL :: l_data_region
!              Flag for region of data in the file
!
!
!- End of header
!
!
!
  CALL get_free_unit(ierr, iu_input)
  IF ( ierr /= i_normal ) RETURN
  CALL open_file_in(ierr, iu_input, &
                      'Enter the name of the file of wavelengths.')
  IF ( ierr /= i_normal ) RETURN
!
  l_data_region = .FALSE.
  n_wavelength  = 0
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
        n_wavelength = n_wavelength + 1
        IF ( n_wavelength > nd_wavelength_scat ) THEN
          WRITE(iu_err, '(/A)') '*** Error: There are ' &
                //'too many wavelengths in the input file.'
          ierr = i_err_fatal
          RETURN
        ENDIF
        read(iu_input, *) wavelength(n_wavelength)
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
!
  RETURN
!
END
