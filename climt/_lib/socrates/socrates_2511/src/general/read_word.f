! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a blank-delimited substring from a string.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE read_word(ierr
     &  , string, length
     &  , i_begin, i_end, word
     &  )
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
!     Dummy arguments
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
      CHARACTER
     &    string*(*)
!           Main string
      INTEGER, Intent(IN) ::
     &    length
!           Length of string
      INTEGER, Intent(INOUT) ::
     &    i_begin
!           Beginning of word
      CHARACTER !, Intent(OUT)
     &    word*16
!           Word read (right justified in field)
      INTEGER, Intent(OUT) ::
     &    i_end
!           End of word
!
!     Local variables
      INTEGER
     &    k
!           Loop variable
!
!
!     Advance the beginning of the word to a non-blank character.
      IF (i_begin > length) THEN
        WRITE(iu_err, '(/a, a)')
     &    '*** error: attempt to read beyond END of string'
     &    , string
        ierr=i_err_fatal
        RETURN
      ENDIF
!
      DO WHILE ( (string(i_begin:i_begin) == ' ').AND.
     &           (i_begin < length) )
        i_begin=i_begin+1
      ENDDO
      i_end=i_begin
      DO WHILE ( (string(i_end+1:i_end+1) /= ' ').AND.
     &           (i_end < length) )
        i_end=i_end+1
      ENDDO
!
!     Check limits on length of word.
      IF ((i_end-i_begin+1) > 16) THEN
        WRITE(iu_err, '(/a, a)')
     &    '*** error: a substring is too long in the string:'
     &    , string
        ierr=i_err_fatal
      ENDIF
      DO k=0, i_end-i_begin
        word(16-k:16-k)=string(i_end-k:i_end-k)
      ENDDO
      DO k=i_end-i_begin+1, 15
        word(16-k:16-k)=' '
      ENDDO
!
!
!
      RETURN
      END
