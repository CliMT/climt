! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to return interactive characteristics fo operation.
LOGICAL FUNCTION set_interactive &
!
()
!
! Description:
!   This function returns the flag indicating whether the program
!   runs interactively.
!
! Method:
!   Straightforward.
!
!
!
!
!
  IMPLICIT NONE
!
!
! Local variables.
  CHARACTER (LEN=20), Parameter :: lock_file = 'radiation_code.lock '
!           Name of locking file
!
!
!- End of Header
!
   INQUIRE(FILE=lock_file, EXIST=set_interactive)
   set_interactive = .NOT.set_interactive
!
!
!
  RETURN
END FUNCTION set_interactive
