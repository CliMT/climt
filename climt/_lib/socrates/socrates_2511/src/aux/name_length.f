! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to determine the length of non-blank substring.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      FUNCTION name_length(name)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
!
!
      IMPLICIT NONE
!
!
!     Function
      INTEGER ::
     &    name_length
!           Length of name
!
!     Dummy argument.
      CHARACTER !, Intent(IN)
     &    name*(*)
!           Name to be striped
!
!     Local variable
      INTEGER
     &    j
!           Loop variable
!
!
      j=1
1     if (name(j:j) == ' ') then
        name_length=j
        RETURN
      ELSE IF (j < len(name)) THEN
        j=j+1
        goto 1
      ELSE
        name_length=len(name)
        RETURN
      ENDIF
!
!
!
      END
