! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to write a file in CDL format.
!
! Purpose:
!   This subroutine selects a unit, opens a file and writes
!   out a field in CDL format.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE write_cdl(ierr
     &  , filename
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  , n_dimension, dimension_name, dimension_type, dimension_unit
     &  , dimension_long, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , n_var, var_name, var_type, var_unit, var_long
     &  , n_dimension_var, list_dimension_var
     &  , n_data, data_int, data_fl
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE error_pcf
      USE def_std_io_icf
      USE filenamelength_mod, ONLY: filenamelength
!
      IMPLICIT NONE
!
!
!
!
!     Dummy arguments:
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
!
      INCLUDE 'cdl_struc.finc'
!
      CHARACTER !, Intent(IN)
     &    filename*(*)
!           Name of the output file
!
!
!     Local variables:
      CHARACTER (LEN=filenamelength) :: file_out
!           Compressed name of output
      CHARACTER
     &    text*80
!           Text for output file
     &  , string(2)*80
!           Temporary strings
      INTEGER
     &    iunit
!           Unit number for output
     &  , length_out
!           Length of name of output
     &  , ios
!           I/O error flag
     &  , length(3)
!           Length of character strings
     &  , length_text
!           Length of text string
     &  , i
!           Loop variable
     &  , iv
!           Loop variable
      LOGICAL
     &    l_exist
!           Existence flag for file
!
!     Subroutines called
      EXTERNAL
     &    get_free_unit, remove_blank, write_cdl_field
!
!
!
!     Find a free unit and open the unit for output.
      CALL get_free_unit(ierr, iunit)
      IF (ierr /= i_normal) RETURN
!
!     Check whether the file exists.
      CALL remove_blank(filename, file_out, length_out)
      INQUIRE(file=file_out(1: length_out), exist=l_exist)
      IF (l_exist) THEN
        WRITE(iu_err, '(3(/a))') '*** Error: The file '
     &    , file_out(1: length_out)
     &    , 'already exists: it will not be overwritten.'
        ierr=i_err_io
        RETURN
      ENDIF
!
!     Open the file for output
      OPEN(unit=iunit, file=file_out(1: length_out), iostat=ios
     &   )
      IF (ios /= 0) THEN
        WRITE(iu_err, '(/a, /a, /a, i3, a, i3)')
     &    '*** Error: The file'
     &    , file_out(1: length_out)
     &    , 'could not be opened on unit ', iunit, ': iostat = ', ios
        ierr=i_err_io
        RETURN
      ENDIF
!
!
      text='netcdf '//file_out(1: length_out)//'{'
      WRITE(iunit, '(a80)') text
!
!     Dimensions:
      WRITE(iunit, '(/, a11)') 'dimensions:'
      DO i=1, n_dimension
        WRITE(iunit, '(t5, a14, t19, a3, i3, a1)')
     &    dimension_name(i), ' = ', dimension_size(i), ';'
      ENDDO
!
!     Variables:
      WRITE(iunit, '(//, a10)') 'variables:'
      DO i=1, n_dimension
!
        CALL remove_blank(dimension_type(i), string(1), length(1))
        CALL remove_blank(dimension_name(i), string(2), length(2))
        text=string(1)(1: length(1))//' '//
     &    string(2)(1: length(2))//'('//
     &    string(2)(1: length(2))//');'
        WRITE(iunit, '(t5, a74)') text
!
!       Attributes:
        CALL remove_blank(dimension_name(i), string(1), length(1))
        CALL remove_blank(dimension_unit(i), string(2), length(2))
        text=string(1)(1: length(1))//':units = "'//
     &    string(2)(1: length(2))//'";'
        WRITE(iunit, '(t14, a65)') text
        CALL remove_blank(dimension_long(i), string(2), length(2))
        text=string(1)(1: length(1))//':title = "'//
     &    string(2)(1: length(2))//'";'
        WRITE(iunit, '(t14, a64)') text
!
      ENDDO
!
!     Declaration of the Fields:
!     The dimenions are written in reverse order since CDL stores
!     arrays like C rather than FORTRAN: when the array is written
!     out normally from the program this will give the right effect.
      DO iv=1, n_var
!
        CALL remove_blank(var_type(iv), string(1), length(1))
        CALL remove_blank(var_name(iv), string(2), length(2))
        text=string(1)(1: length(1))//' '//string(2)(1: length(2))
     &    //'('
        length_text=length(1)+1+length(2)+1
        DO i=n_dimension_var(iv), 2, -1
          CALL remove_blank(dimension_name(list_dimension_var(i, iv))
     &      , string(2), length(2))
          text=text(1: length_text)//string(2)(1: length(2))//','
          length_text=length_text+length(2)+1
        ENDDO
        CALL remove_blank(dimension_name(list_dimension_var(1, iv))
     &    , string(2), length(2))
        text=text(1: length_text)//string(2)(1: length(2))//');'
        length_text=length_text+length(2)+2
        WRITE(iunit, '(/, t5, a75)') text
!
!       Attributes of the fields:
        CALL remove_blank(var_name(iv), string(1), length(1))
        CALL remove_blank(var_unit(iv), string(2), length(2))
        text=string(1)(1: length(1))//':units = "'//
     &    string(2)(1: length(2))//'";'
        WRITE(iunit, '(t14, a65)') text
        CALL remove_blank(var_long(iv), string(2), length(2))
        text=string(1)(1: length(1))//':title = "'//
     &    string(2)(1: length(2))//'";'
        WRITE(iunit, '(t14, a)') text
!
      ENDDO
!
!
!     Data:
      WRITE(iunit, '(/, "data:")')
      DO i=1, n_dimension
        CALL remove_blank(dimension_name(i), string(1), length(1))
        CALL write_cdl_field(iunit, string(1)(1: length(1))
     &    , dimension_type(i), dimension_size(i)
     &    , dimension_array_int(1, i), dimension_array_fl(1, i))
      ENDDO
      DO iv=1, n_var
        CALL remove_blank(var_name(iv), string(1), length(1))
        CALL write_cdl_field(iunit, string(1)(1: length(1))
     &    , var_type(iv), n_data(iv), data_int(1, iv), data_fl(1, iv))
      ENDDO
!
!
      WRITE(iunit, '(/, "}")')
!
      CLOSE(iunit)
!
!
!
      RETURN
      END
