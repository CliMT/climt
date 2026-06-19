! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to open a non-specific file for output.
!
! Method:
!       Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE open_file_out(ierr, iunit, request)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!     Dummy arguments.
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
      INTEGER, Intent(IN) ::
     &    iunit
!           Unit number
      CHARACTER !, Intent(IN)
     &    request*(*)
!           Request prompting for file
!
!     Local varaibles.
!
      LOGICAL
     &    lock_code
!           Logical to forbid interactive looping
      EXTERNAL
     &    lock_code
!
      CHARACTER
     &    char_yn*1
!           Terminal response
     &  , file_out*132
!           Name of output file
      INTEGER
     &    ios
!           I/O error flag
      LOGICAL
     &    l_open
!           Open flag
     &  , l_exist
!           Existence flag
!
!
!     Check whether the unit for the output file is already open.
      INQUIRE(unit=iunit, opened=l_open)
      IF (l_open) THEN
        WRITE(iu_err, '(/, a16, i5, a)') '*** Error: The unit ', iunit
     &    , ' is already open elsewhere.'
        ierr=i_err_io
        RETURN
      ENDIF
!
!     Request the name of the output file.
      WRITE(*, '(/a)') request
      READ(iu_stdin, '(a)') file_out
!
!     Check whether the file exists.
100   inquire(file=file_out, exist=l_exist)
      IF (l_exist) THEN
        INQUIRE(file=file_out, opened=l_open)
        WRITE(iu_err, '(/a)') 'This file already exists.'
        IF (l_open) THEN
          WRITE(*, '(a)') 'This file is already open.'
        ENDIF
        IF (lock_code(.true.)) THEN
          ierr=i_err_io
          RETURN
        ELSE
          WRITE(*, '(/a/)') 'Do you wish to overwrite? (Y/N).'
          READ(iu_stdin, '(a)') char_yn
          IF ( (char_yn /= 'Y').AND.(char_yn /= 'y') ) THEN
            WRITE(*, '(/a/)')
     &        'Please specify another name or quit (?q).'
            READ(iu_stdin, '(a)') file_out
            IF ( (file_out(1:2) == '?Q').OR.
     &           (file_out(1:2) == '?q') ) THEN
              WRITE(*, '(/a/)') '*** Program terminated.'
              ierr=i_err_fatal
              RETURN
            ELSE
              goto 100
            ENDIF
          ENDIF
        ENDIF
      ENDIF
!
!     The file may now safely be opened.
      OPEN(UNIT=iunit, FILE=file_out, iostat=ios
     &  , status='unknown'
     &  )
      IF (ios /= 0) THEN
        WRITE(iu_err, '(3(/a))') '*** Error: The file', file_out
     &  , ' could not be opened.'
        ierr=i_err_io
        RETURN
      ENDIF
!
!
!
      RETURN
      END
