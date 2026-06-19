! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to open a file of averaged scattering data for output.
!
SUBROUTINE open_average_90(l_write, iu_average, l_interactive, ierr)
!
! Description:
!   This subroutine opens a file for the output of averaged 
!   scattering data.
!
! Modules used.
  USE realtype_rd
  USE def_std_io_icf
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
!
! Dummy variables.
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive operation
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  INTEGER, Intent(OUT) :: iu_average
!   Unit number for output of averaged data
  LOGICAL, Intent(OUT) :: l_write
!   Writing flag
!
!
! Local variables.
!
  CHARACTER (LEN=1) :: char_yn
!   Response variable
!
! Subroutines called:
  EXTERNAL open_file_out
!
!
!
! Determine whether outpur for each block is required.
  WRITE(iu_stdout, '(/a)') 'Do you wish to print the averaged ' // &
    'scattering properties for each block? (y/n)'
  WRITE(iu_stdout, '(3x, a)') 'This must be selected for aerosols.'
  DO
    READ(iu_stdin, "(a)") char_yn
    IF ( (char_yn == 'n') .OR. (char_yn == 'n') ) THEN
      l_write= .FALSE. 
      RETURN
    ELSE IF ( (char_yn == 'y') .OR. (char_yn == 'y') ) THEN
      EXIT
    ELSE 
      WRITE(iu_err, '(a)') '+++ illegal reponse:'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(a)') 'Please re-enter.'
      ELSE
        ierr=i_err_fatal
        RETURN
      ENDIF
    ENDIF
  ENDDO
!
! Set the writing flag to .TRUE.: this will be checked as each block
! is averaged.
  l_write = .TRUE. 
!
  CALL get_free_unit(ierr, iu_average)
  CALL open_file_out(ierr, iu_average, &
    'Give the name of the file to contain the averaged scattering data.')
  IF (ierr /= i_normal) RETURN
!
!
!
  RETURN
END SUBROUTINE open_average_90
