! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to open a non-specific file for input.
!
! Method:
!       Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE open_file_in(ierr, iunit, request)
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
!     Dummy variables
      CHARACTER !, Intent(IN)
     &    request*(*)
!           Prompt to user
      INTEGER, Intent(IN) ::
     &    iunit
!           Number of unit
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
!
!     Local variables
!
      LOGICAL
     &    lock_code
!           Logical to forbid interactive looping
      EXTERNAL
     &    lock_code
!
      CHARACTER
     &    file_in*132
!           Name of input file
      LOGICAL
     &    l_exist
!           Existence flag for file
     &  , l_open
!           Open flag for file
      INTEGER
     &    ios
!           I/O error flag
     &  , j
!           Loop variable
!
!
!     Obtain the name of the file of data to be read in.
      WRITE(*, '(/a)') request
!
1     do j=1, 132
!        Clear the input filename.
         file_in(j:j)=' '
      ENDDO
      READ(iu_stdin, '(a)', iostat=ios) file_in
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ Unexpected response:'
        IF (lock_code(.true.)) THEN
          ierr=i_err_io
          RETURN
        ELSE
          WRITE(*, '(a)') 'Please retype.'
          goto 1
        ENDIF
      ENDIF
!
!     "?Q" is the response to quit.
      IF ( (file_in(1:2) == '?Q').OR.(file_in(1:2).eq.'?q') ) THEN
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Check that the unit for input is free.
      INQUIRE(UNIT=iunit, OPENED=l_open)
      IF (l_open) THEN
        WRITE(iu_err, '(/a)') '*** Error: The unit for reading the '
     &    //'file is already open elsewhere.'
        ierr=i_err_io
        RETURN
      ENDIF
!
!     Check for the existence of the file.
      INQUIRE(FILE=file_in, EXIST=l_exist)
      IF (.NOT.l_exist) THEN
        WRITE(iu_err, '(a)') 'This file does not exist:'
        IF (lock_code(.true.)) THEN
          ierr=i_err_exist
          RETURN
        ELSE
          WRITE(iu_err, '(a)')
     &      'Please give another name or type "?q" to quit.'
          goto 1
        ENDIF
      ENDIF
!
!     Now open the file.
      OPEN(UNIT=iunit, FILE=file_in, iostat=ios
     &  , status='old'
     &  )
      IF (ios /= 0) THEN
        WRITE(iu_err, '(3(/a))')
     &    '*** Error: The file ', file_in, 'could not be opened.'
        ierr=i_err_io
        RETURN
      ENDIF
!
!
!
      RETURN
      END
