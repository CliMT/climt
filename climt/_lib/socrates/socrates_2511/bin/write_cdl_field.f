! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to write a field of data to a CDL file.
!
! Purpose:
!   This subroutine writes a single field of numerical data to
!   a CDL file.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE write_cdl_field(iunit 
     &  , name_data, type_data
     &  , n_data, data_int, data_fl)
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
      INTEGER, Intent(IN) ::
     &    iunit
!           Unit for output
     &  , n_data
!           Number of data in field
      INTEGER, Intent(IN) ::
     &    data_int(n_data)
!           Values of integer data
      REAL  (RealK), Intent(IN) ::
     &    data_fl(n_data)
!           Values of floating-point data
      CHARACTER !, Intent(IN)
     &    name_data*(*)
!           Name of the data variable
     &  , type_data*(*)
!           Type of the data variable
!
!     Local variables:
      INTEGER
     &    n_not_last
!           Number of non-final lines
     &  , n_left
!           Number of data left for final line
     &  , n_start
!           Starting datum for current line of output
     &  , len_num
!           Length of the string holding the number
     &  , n_per_line
!           Number of data to write on each line
     &  , i
!           Loop variable
     &  , k
!           Loop variable
      CHARACTER
     &    ch_num*5
!           Character string for current number
     &  , fmt_elem*5
!           String holding format for output of a single variable
     &  , fmt_str*80
!           Full format string for output
!
!
!
!
!     Set the format of the written data.
      IF (type_data(1:min(5, len(type_data))) == 'float') THEN
        n_per_line=4
        fmt_elem='e13.6'
      ELSE IF (type_data(1:min(3, len(type_data))) == 'int') THEN
        n_per_line=8
        fmt_elem='i5'
      ENDIF
!
!     Determine the number of lines of each type required. The
!     last line will be exceptional
      n_not_last=(n_data-1)/n_per_line
      n_left=n_data-n_per_line*n_not_last
      n_start=1
!
!     Write the field out in three groups. The first line will
!     contain the name of the field, other full lines follow.
!     The last line will have the terminating semicolon and 
!     may be short. 
      IF (n_not_last > 0) THEN
!
!       Convert the number of data per line into a string.
        WRITE(ch_num, '(i5)') n_per_line
        CALL remove_blank(ch_num, ch_num, len_num)
!
!       The first line is not the last.
        fmt_str='(t4, a14, t18, " = ",'//ch_num(1: len_num)
     &    //'('//fmt_elem//', ",") )'
        IF (type_data(1:min(5, len(type_data))) == 'float') THEN
          WRITE(iunit, fmt=fmt_str) name_data
     &      , (data_fl(k), k=n_start, n_start+n_per_line-1)
        ELSE IF (type_data(1:min(3, len(type_data))) == 'int') THEN
          WRITE(iunit, fmt=fmt_str) name_data
     &      , (data_int(k), k=n_start, n_start+n_per_line-1)
        ENDIF
!       Middle lines.
        DO i=2, n_not_last
          n_start=n_start+n_per_line
          fmt_str='(t21, '//ch_num(1: len_num)
     &      //'('//fmt_elem//', ",") )'
          IF (type_data(1:min(5, len(type_data))) == 'float') THEN
            WRITE(iunit, fmt=fmt_str) 
     &        (data_fl(k), k=n_start, n_start+n_per_line-1)
          ELSE IF (type_data(1:min(3, len(type_data))) == 'int') THEN
            WRITE(iunit, fmt=fmt_str) 
     &        (data_int(k), k=n_start, n_start+n_per_line-1)
          ENDIF
        ENDDO
!
!       If N_LEFT exceeds 1 we need a repeated group in the format
!       string terminated with a comma: the last part of the string
!       must always end with a semi-colon.
        IF (n_left > 1) THEN
          WRITE(ch_num, '(i5)') n_left-1
          CALL remove_blank(ch_num, ch_num, len_num)
          fmt_str='(t21, '//ch_num(1: len_num)
     &      //'('//fmt_elem//', ","), '//fmt_elem//', ";" )'
        ELSE
          fmt_str='(t21, '//fmt_elem//', ";" )'
        ENDIF
        n_start=n_start+n_per_line
        IF (type_data(1:min(5, len(type_data))) == 'float') THEN
          WRITE(iunit, fmt=fmt_str) 
     &      (data_fl(k), k=n_start, n_start+n_left-1)
        ELSE IF (type_data(1:min(3, len(type_data))) == 'int') THEN
          WRITE(iunit, fmt=fmt_str) 
     &      (data_int(k), k=n_start, n_start+n_left-1)
        ENDIF
!
      ELSE
!
!       There is only one line.
!       If N_DATA exceeds 1 we need a repeated group in the format
!       string terminated with a comma: the last part of the string
!       must always end with a semi-colon.
!
!       Convert the number of data per line into a string.
        WRITE(ch_num, '(i5)') n_per_line
        CALL remove_blank(ch_num, ch_num, len_num)
!
        IF (n_data > 1) THEN
          WRITE(ch_num, '(i5)') n_data-1
          fmt_str='(t4, a14, t18, " = ",'//ch_num(5-len_num:5)
     &      //'('//fmt_elem//', ","), '//fmt_elem//',";" )'
        ELSE
          fmt_str='(t4, a14, t18, " = ",'//fmt_elem//',";" )'
        ENDIF
        IF (type_data(1:min(5, len(type_data))) == 'float') THEN
          WRITE(iunit, fmt=fmt_str) name_data
     &      , (data_fl(k), k=n_start, n_start+n_data-1)
        ELSE IF (type_data(1:min(3, len(type_data))) == 'int') THEN
          WRITE(iunit, fmt=fmt_str) name_data
     &      , (data_int(k), k=n_start, n_start+n_data-1)
        ENDIF
!
      ENDIF
!
!
!
      RETURN
      END
