! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to open anon-specific file for output.
!
SUBROUTINE open_file_out_90 &
!
(iunit, l_interactive, request, file_out, ierr)
!
!
! Description:
!   This subroutine constructs a formatted file containing a k-fit.
!
! Method:
!   Straightforward.
!
! Modules used.
  USE def_std_io_icf
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
! Dummy arguments.
  INTEGER, Intent(IN) :: iunit
!   Unit number
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive operation
  CHARACTER  (LEN=*), Intent(IN) :: request
!   Request prompting for file
  CHARACTER  (LEN=*), Intent(OUT) :: file_out
!   Name of output file
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
!
!
!
! Local variables.
  CHARACTER  (LEN=1) :: char_yn
!   Terminal response
  INTEGER :: ios
!   I/O error flag
  LOGICAL :: l_open
!   Open flag
  LOGICAL :: l_exist
!   Existence flag
!
!
!
! Check whether the unit for the output file is already open.
  INQUIRE(UNIT=iunit, OPENED=l_open)
  IF (l_open) THEN
    WRITE(iu_err, '(/, a16, i5, a)') &
      '*** Error: The unit ', iunit, ' is already open elsewhere.'
    ierr=i_err_io
    RETURN
  ENDIF
!
! Request the name of the output file.
  WRITE(iu_stdout, '(/a)') request
  READ(iu_stdin, '(a)') file_out
!
! Check whether the file exists.
  DO
    INQUIRE(FILE=file_out, EXIST=l_exist)
    IF (l_exist) THEN
      WRITE(iu_err, '(/a)') 'This file already exists.'
      INQUIRE(FILE=file_out, OPENED=l_open)
      IF (l_open) THEN
        WRITE(iu_stdout, '(a)') 'This file is already open.'
      ENDIF
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(/a/)') 'Do you wish to overwrite? (y/n).'
        READ(iu_stdin, '(a)') char_yn
        IF ( (char_yn /= 'y') .AND. (char_yn /= 'Y') ) THEN
          WRITE(iu_stdout, '(/a/)') &
            'Please specify another name or quit (?q).'
          READ(iu_stdin, '(a)') file_out
          IF ( (file_out(1:2) == '?Q') .OR. &
               (file_out(1:2) == '?q') ) THEN
            WRITE(iu_stdout, '(/a/)') '*** Routine terminated.'
            ierr=i_err_fatal
            RETURN
          ENDIF
        ELSE
          EXIT
        ENDIF
      ELSE
        ierr=i_err_io
        RETURN
      ENDIF
    ELSE
      EXIT
    ENDIF
  ENDDO
!
! The file may now safely be opened.
  OPEN(UNIT=iunit, FILE=file_out, IOSTAT=ios, STATUS='unknown')
  IF (ios /= 0) THEN
    WRITE(iu_err, '(3(/a))') &
      '*** Error: The file', file_out, ' could not be opened.'
    ierr=i_err_io
    RETURN
  ENDIF
!
!
!
  RETURN
END SUBROUTINE open_file_out_90
