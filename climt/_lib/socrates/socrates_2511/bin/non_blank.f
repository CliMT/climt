! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to determine whether a line is blank.
!
! Method:
!   Straightforward
!
!- ---------------------------------------------------------------------
      function non_blank(line)
!
!
!
      IMPLICIT NONE
!
!
!     Dummy arguments
      CHARACTER !, Intent(IN)
     &    line*(*)
!           Input line
!
      LOGICAL ::
     &    non_blank
!           Functional logical
!
!     Local variables
      INTEGER
     &    j
!           Loop index
!
!
      j=len(line)
      non_blank=.false.
1     if (j.gt.0) then
        IF (line(j:j) /= ' ') THEN
          non_blank=.true.
          RETURN
        ELSE
          j=j-1
          goto 1
        ENDIF
      ELSE
        RETURN
      ENDIF
!
!
!
      END
