! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to split a CDL line into numbers
!
! Purpose:
!   This subroutine extracts numbers from a line of CDL data.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE split_cdl_line(ierr, in_string, type
     &  , n, x_int, x_fl, i_start, l_end, n_read)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!
!     Dummy arguments:
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error Flag
!
      CHARACTER !, Intent(IN)
     &    in_string*(*)
!           Input string
     &  , type*(*)
!           Type of field
      INTEGER, Intent(IN) ::
     &    n
!           Size of data array
     &  , i_start
!           Starting point in array
      LOGICAL, Intent(INOUT) ::
     &    l_end
!           Flag for end of the data
      INTEGER, Intent(INOUT) ::
     &    x_int(n)
!           Integral data
      REAL  (RealK), Intent(INOUT) ::
     &    x_fl(n)
!           Floating point data
      INTEGER, Intent(OUT) ::
     &    n_read
!           Number of elements read
!
!     Local variables:
      INTEGER
     &    i_begin
!           Beginning the substring
     &  , j
!           Loop variable
     &  , ios
!           I/O error flag
!
!
!
      n_read=0
!
      j=1
      i_begin=1
      DO WHILE( (j <= len(in_string)).AND.(.NOT.l_end) )
        IF ( (in_string(j:j) == ',').OR.(in_string(j:j) == ';') ) THEN
          IF (type(1:min(5, len(type))) == 'float') THEN
            READ(in_string(i_begin:j-1), *, iostat=ios) 
     &        x_fl(i_start+n_read)
          ELSE IF (type(1:min(3, len(type))) == 'int') THEN
            READ(in_string(i_begin:j-1), *, iostat=ios) 
     &        x_int(i_start+n_read)
          ENDIF
          IF (ios /= 0) THEN
            ierr=i_err_io
            RETURN
          ENDIF
          n_read=n_read+1
          IF (in_string(j:j) == ';') THEN
            l_end=.true.
            RETURN
          ELSE
            i_begin=j+1
          ENDIF
        ENDIF
        j=j+1
      ENDDO
!
!
!
      RETURN
      END
