! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to check for the existence of the locking file.
!
! Method:
!   The inquire statement is used. an argument must be given
!   to the function, otherwise it will not be called.
!
!- ---------------------------------------------------------------------
      FUNCTION lock_code(l_apply)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
!
!
      IMPLICIT NONE
!
!
!
!     Dummy arguments.
      LOGICAL, Intent(IN) ::
     &    l_apply
!           Argument to apply test: an argument is required for the
!           function actually to be called, even though it will
!           always be true.
      LOGICAL ::
     &    lock_code
!           Logical for existence of lock-file.
!
!     Local variables.
      CHARACTER
     &    lock_file*30
!           Name of lock file
!
      DATA lock_file/'radiation_code.lock           '/
!
!
      IF (l_apply) INQUIRE(FILE=lock_file, EXIST=lock_code)
!
      IF (lock_code) THEN
        WRITE(iu_err, '(/a, /a)') '*** Error in input.'
     &    , 'This is treated as fatal since a locking file exists.'
      ENDIF
!
!
!
      RETURN
      END
