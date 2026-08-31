! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to remove blanks at the ends of a line.
!
! Purpose:
!   This subroutine removes blank at the beginning and end of a line.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE remove_blank(in_string, out_string, out_length)
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
!     Dummy arguments:
      CHARACTER !, Intent(IN)
     &    in_string*(*)
!           Input string
      CHARACTER !, Intent(OUT)
     &    out_string*(*)
!           Output string
      INTEGER, Intent(OUT) ::
     &    out_length
!           Length of output string
!
!     Local variables:
      INTEGER
     &    i_begin
!           Beginning of non-blank substring
     &  , i_end
!           End of non-blank substring
!
!
!
      i_begin=1
      DO WHILE ( (i_begin < len(in_string)).AND.
     &  (in_string(i_begin: i_begin) == ' ') )
        i_begin=i_begin+1
      ENDDO
!
      i_end=len(in_string)
      DO WHILE ( (i_end > 1).AND.
     &  (i_end > i_begin).AND.
     &  (in_string(i_end: i_end) == ' ') )
        i_end=i_end-1
      ENDDO
!
      out_length=i_end+1-i_begin
      out_string(1: out_length)=in_string(i_begin: i_end) 
!
!
!
      RETURN
      END
