! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a file in CDL format.
!
! Purpose:
!   This subroutine selects a unit, opens the file and reads
!   in the field in CDL format. NOTE that for ease of programming
!   in FORTRAN CDL-syntax is actually restricted here. In the section 
!   declaring variables only one variable is permitted on each line.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE read_cdl(ierr
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
      CHARACTER !, Intent(IN)
     &    filename*(*)
!           Name of input file
!
      INCLUDE 'cdl_struc.finc'
!
!
!     Local variables:
      CHARACTER (LEN=filenamelength) :: file_in
!           Compressed name of output
      CHARACTER
     &    text*132
!           Text for output file
     &  , name_temp*80
!           Temporary string
      INTEGER
     &    iunit
!           Unit number for input
     &  , length_in
!           Length of name of input
     &  , ios1
!           I/O error flag
     &  , ios2
!           I/O error flag
     &  , i
!           Loop variable
     &  , iv
!           Loop variable
     &  , j
!           Loop variable
     &  , i_end
!           Parsing integer
     &  , i_equal
!           Parsing integer
     &  , i_begin_attribute
!           Parsing integer
     &  , i_end_attribute
!           Parsing integer
     &  , i_open
!           Parsing integer
     &  , i_close
!           Parsing integer
     &  , i_start
!           Parsing integer
     &  , i_colon
!           Parsing integer
     &  , i_comma(nd_cdl_dimen)
!           Parsing integer
     &  , n_comma
!           Parsing integer
     &  , i_first
!           Parsing integer
     &  , i_last
!           Parsing integer
     &  , length_temp
!           Length of temporary variable
     &  , i_dimension
!           Current dimension
     &  , size_temp
!           Temporary size variable
     &  , n_read
!           Number of elements read
      LOGICAL
     &    l_exist
!           Existence flag for file
     &  , l_non_comment
!           Flag for non-commented part of line
     &  , l_declaration
!           Flag for declarative line
     &  , l_attribution
!           Flag for attributive line
     &  , l_parse
!           Flag to parse line
     &  , l_end
!           Logical detecting end of input data
     &  , l_loop
!           Flag to control the execution of WHILE-loops with
!           split tests: this is required on some machines
!           to prevent out of bounds messages
!
!     Subroutines called
      EXTERNAL
     &    get_free_unit, remove_blank
!
!
!
!     Initial setting of dimensions and variables for error checking.
      n_dimension=0
      n_var=0
!
!     Initialize the character strings
      DO i=1, nd_cdl_dimen
        DO j=1, len(dimension_name(i))
          dimension_name(i)(j:j)=' '
        ENDDO
        DO j=1, len(dimension_unit(i))
          dimension_unit(i)(j:j)=' '
        ENDDO
        DO j=1, len(dimension_long(i))
          dimension_long(i)(j:j)=' '
        ENDDO
      ENDDO
      DO i=1, nd_cdl_var
        DO j=1, len(var_name(i))
          var_name(i)(j:j)=' '
        ENDDO
        DO j=1, len(var_unit(i))
          var_unit(i)(j:j)=' '
        ENDDO
        DO j=1, len(var_long(i))
          var_long(i)(j:j)=' '
        ENDDO
      ENDDO
!
!     Find a free unit and open the unit for input.
      CALL get_free_unit(ierr, iunit)
      IF (ierr /= i_normal) RETURN
!
!     Check whether the file exists.
      CALL remove_blank(filename, file_in, length_in)
      INQUIRE(file=file_in(1: length_in), exist=l_exist)
      IF (.NOT.l_exist) THEN
        WRITE(iu_stdout, '(3(/a))') 'Warning: The file '
     &    , file_in(1: length_in)
     &    , 'does not exist.'
        ierr=i_err_exist
        RETURN
      ENDIF
!
!     Open the file for input.
      OPEN(unit=iunit, file=file_in(1: length_in), iostat=ios1
     &   )
      IF (ios1 /= 0) THEN
        WRITE(iu_err, '(/a, /a, /a, i3, a, i3)')
     &    '*** Error: The file'
     &    , file_in(1: length_in)
     &    , 'could not be opened on unit ', iunit, ': iostat = ', ios1
        ierr=i_err_io
        RETURN
      ENDIF
!
!
!     Search the file for the block of dimensions.
      READ(iunit, '(a)', END=999) text
      DO WHILE (text(1: 11) /= 'dimensions:')
        READ(iunit, '(a)', END=999) text
      ENDDO
!     Now read successive lines of dimensions until we reach the next
!     block.
      n_dimension=0
      READ(iunit, '(a)', END=999) text
      DO WHILE (text(1:10) /= 'variables:')
!       Parse the line.
        j=1
        l_non_comment=.true.
        i_end=0
        i_equal=0
        DO WHILE (l_non_comment.AND.(j <= len(text)))
          IF (text(j:j) == ';') THEN
            i_end=j
          ELSE IF (text(j:j) == '=') THEN
            i_equal=j
          ELSE IF (text(j:min(len(text),j+1)) == '//') THEN
            l_non_comment=.false.
            i_end=j
          ENDIF
          j=j+1
        ENDDO
!       Now check for valid syntax.
        IF ( (i_equal > 1).AND.(i_end > i_equal+1) ) THEN
!         This is a line of data.
          name_temp=' '
          name_temp(1: i_equal-1)=text(1: i_equal-1)
          READ(text(i_equal+1: i_end-1), *, iostat=ios2) size_temp
          IF ( (ios1 == 0).AND.(ios2 == 0) ) THEN
            IF ( (n_dimension >= 7).OR.
     &           (n_dimension >= nd_cdl_dimen) ) THEN
              WRITE(iu_err, '(4(/a), 1x, i5, a1)')
     &          '*** Error reading file'
     &          , file_in
     &          , 'no more than 7 dimensions are '
     &          //'permitted in a cdl-file.'
     &          , 'this code was compiled with space for '
     &          , nd_cdl_dimen, '.'
              ierr=i_err_io
              goto 999
            ENDIF
            n_dimension=n_dimension+1
            CALL remove_blank(name_temp, dimension_name(n_dimension)
     &        , length_temp)
            dimension_size(n_dimension)=size_temp
          ENDIF
        ENDIF
!
        READ(iunit, '(a)', END=999) text
!
      ENDDO
!
!     Now move on to the group of variables: no variables have yet been 
!     declared.
      n_var=0
      READ(iunit, '(a)', END=999) text
      DO WHILE (text(1:5) /= 'data:')
!
        j=1
        l_declaration=.false.
        l_attribution=.false.
        l_non_comment=.true.
        i_begin_attribute=0
        i_end_attribute=0
        i_end=0
        i_open=0
        i_close=0
        i_colon=0
        i_equal=0
        i_start=0
        n_comma=0
!
        DO WHILE (l_non_comment.AND.(j <= len(text)))
          IF (text(j:j) == ';') THEN
            i_end=j
          ELSE IF (text(j:min(len(text), j+5)) == 'float ') THEN
            IF (.NOT.l_attribution) l_declaration=.true.
            i_start=j+5
          ELSE IF (text(j:min(len(text), j+3)) == 'int ') THEN
!           If the string "int" appears in an attribution we could
!           have confusion.
            IF (.NOT.l_attribution) l_declaration=.true.
            i_start=j+3
          ELSE IF (text(j:min(len(text), j+1)) == '//') THEN
            l_non_comment=.false.
            i_end=j
          ELSE IF (text(j:j) == '=') THEN
            i_equal=j
          ELSE IF (text(j:j) == ':') THEN
            l_attribution=.true.
            i_colon=j
          ELSE IF (text(j:j) == '"') THEN
            IF (i_begin_attribute > 0) THEN
              i_end_attribute=j-1
            ELSE
              i_begin_attribute=j+1
            ENDIF
          ELSE IF (text(j:j) == ',') THEN
             n_comma=n_comma+1
             i_comma(n_comma)=j
          ELSE IF (text(j:j) == '(') THEN
             i_open=j
          ELSE IF (text(j:j) == ')') THEN
             i_close=j
          ENDIF
!
          j=j+1
!
        ENDDO
!
        IF (l_declaration) THEN
!
!         Check that the line seems sensible.
          IF ( .NOT.( (i_open-1 > i_start).AND.
     &                (i_close > i_open+1).AND.
     &                (i_end > i_close) ) ) THEN
            WRITE(iu_err, '(3(/a))')
     &        '*** Error reading file'
     &        , file_in
     &        , 'there is an illegal variable declaration.'
            ierr=i_err_io
            goto 999
          ENDIF
!
!         Extract the name of the variable.
          CALL remove_blank(text(i_start:i_open-1)
     &      , name_temp, length_temp)
!         Check against the existing dimensions to see if
!         we have a coordinate variable.
          i=1
          l_loop=(i <= n_dimension)
          IF (l_loop) l_loop=( l_loop.AND.
     &      (name_temp(1:length_temp) /= 
     &      dimension_name(i)(1:length_temp)) )
          DO WHILE (l_loop)
            i=i+1
            l_loop=(i <= n_dimension)
            IF (l_loop) l_loop=( l_loop.AND.
     &        (name_temp(1:length_temp) /= 
     &        dimension_name(i)(1:length_temp)) )
          ENDDO
          IF (i > n_dimension) THEN
!           This is a data variable
            IF (n_var < nd_cdl_var) THEN
              n_var=n_var+1
            ELSE
              WRITE(iu_err, '(/a, a, /a, i3, a1)')
     &          '*** Error: The file ', file_in(1:length_in)
     &          , 'contains too many variables.'
     &          //'space is allocated for only ', nd_cdl_var, '.'
              ierr=i_err_fatal
              RETURN
            ENDIF
            var_name(n_var)(1:length_temp)=name_temp(1:length_temp)
            DO j=length_temp+1, len(var_name(n_var))
              var_name(n_var)(j:j)=' '
            ENDDO
!           Construct the list of dimensions
            n_dimension_var(n_var)=n_comma+1
            i_dimension=0
            l_parse=.true.
            i_first=i_open+1
            IF (n_comma == 0) THEN
              i_last=i_close-1
            ELSE
              i_last=i_comma(1)-1
            ENDIF
            DO WHILE (l_parse)
              CALL remove_blank(text(i_first:i_last)
     &          , name_temp, length_temp)
              i=1
              l_loop=(i <= n_dimension)
              IF (l_loop) l_loop=( l_loop.AND.
     &          (name_temp(1:length_temp) /= 
     &          dimension_name(i)(1:length_temp)) )
              DO WHILE (l_loop)
                i=i+1
                l_loop=(i <= n_dimension)
                IF (l_loop) l_loop=( l_loop.AND.
     &            (name_temp(1:length_temp) /= 
     &            dimension_name(i)(1:length_temp)) )
              ENDDO
              IF (i > n_dimension) THEN
                WRITE(iu_err, '(3(/a))')
     &            '*** Error reading file'
     &            , file_in(1:length_in)
     &            , 'an undeclared dimension is used.'
                ierr=i_err_io
                goto 999
              ELSE
                i_dimension=i_dimension+1
                list_dimension_var(i_dimension, n_var)=i
              ENDIF
!             Check for termination or advance to the next group.
              i_first=i_last+2
              IF (i_dimension == n_dimension_var(n_var)) THEN
                l_parse=.false.
              ELSE IF (i_dimension == n_comma) THEN
                i_last=i_close-1
              ELSE
                i_last=i_comma(i_dimension+1)-1
              ENDIF
            ENDDO
!           Assign the type of the variable.
            CALL remove_blank(text(1:i_start-1)
     &        , name_temp, length_temp)
            var_type(n_var)(1:length_temp)=name_temp(1:length_temp)
          ELSE
!           Check that the CDL convention for coordinate variables
!           is obeyed.
            CALL remove_blank(text(i_open+1:i_close-1)
     &        , name_temp, length_temp)
            IF (name_temp(1:length_temp) /= 
     &        dimension_name(i)(1:length_temp)) THEN
              WRITE(iu_err, '(3(/a))')
     &          '*** Error reading file'
     &          , file_in(1:length_in)
     &          , 'a coordinate variable does not '
     &          //'match its dimension.'
                ierr=i_err_io
                goto 999
            ENDIF
!           Assign the type of the dimension.
            CALL remove_blank(text(1:i_start-1)
     &        , name_temp, length_temp)
            dimension_type(i)(1:length_temp)=name_temp(1:length_temp)
          ENDIF
!
        ELSE IF (l_attribution) THEN
!
!         Check that the line seems sensible.
          IF ( .NOT.( (i_colon > 1).AND.
     &                (i_equal > i_colon+1).AND.
     &                (i_begin_attribute > i_equal+1).AND.
     &                (i_end_attribute >= i_begin_attribute) ) ) THEN
            WRITE(iu_err, '(3(/a))')
     &        '*** Error reading file'
     &        , file_in
     &        , 'there is an illegal variable attribution.'
          ENDIF
!
!
!         Get the name of the variable.
          CALL remove_blank(text(1:i_colon-1)
     &      , name_temp, length_temp)
          i=1
          l_loop=(i <= n_dimension)
          IF (l_loop) l_loop=( l_loop.AND.
     &      (name_temp(1:length_temp) /= 
     &      dimension_name(i)(1:length_temp)) )
          DO WHILE (l_loop)
            i=i+1
            l_loop=(i <= n_dimension)
            IF (l_loop) l_loop=( l_loop.AND.
     &        (name_temp(1:length_temp) /= 
     &        dimension_name(i)(1:length_temp)) )
          ENDDO
          IF (i > n_dimension) THEN
            iv=1
            l_loop=(iv <= n_var)
            IF (l_loop) l_loop=( l_loop.AND.
     &        (name_temp(1:length_temp) /= 
     &        var_name(iv)(1:length_temp)) )
            DO WHILE (l_loop)
              iv=iv+1
              l_loop=(iv <= n_var)
              IF (l_loop) l_loop=( l_loop.AND.
     &          (name_temp(1:length_temp) /= 
     &          var_name(iv)(1:length_temp)) )
            ENDDO
            IF (iv > n_var) THEN
              WRITE(iu_err, '(3(/a))')
     &          '*** Error reading file'
     &          , file_in(1:length_in)
     &          , 'an attribution to an undeclared variable is made.'
              ierr=i_err_io
              goto 999
            ENDIF
!           Check for valid attributes.
            CALL remove_blank(text(i_colon+1:i_equal-1)
     &        , name_temp, length_temp)
            IF (name_temp(1:5) == 'units') THEN
              var_unit(iv)(1:i_end_attribute+1-i_begin_attribute)
     &          =text(i_begin_attribute:i_end_attribute)
            ELSE IF (name_temp(1:5) == 'title') THEN
              var_long(iv)(1:i_end_attribute+1-i_begin_attribute)
     &          =text(i_begin_attribute:i_end_attribute)
            ENDIF
          ELSE
!           Check for valid attributes.
            CALL remove_blank(text(i_colon+1:i_equal-1)
     &        , name_temp, length_temp)
            IF (name_temp(1:5) == 'units') THEN
              dimension_unit(i)(1:i_end_attribute+1-i_begin_attribute)
     &          =text(i_begin_attribute:i_end_attribute)
            ELSE IF (name_temp(1:5) == 'title') THEN
              dimension_long(i)(1:i_end_attribute+1-i_begin_attribute)
     &          =text(i_begin_attribute:i_end_attribute)
            ENDIF
          ENDIF
!
        ENDIF
!
        READ(iunit, '(a)', END=999) text
!
      ENDDO
!
!
!     Now read the data until an error condition is reached, which
!     should be EOF.
      ios1=0
      DO WHILE (ios1 == 0)
!
!       Process lines until an uncommented = sign is found. A 
!       backspace is required as a previous reading of a field
!       may have brought us to the beginning of the new field.
        backspace(iunit)
        i_equal=0
        DO WHILE (i_equal == 0)
          READ(iunit, '(a)', END=999) text
          j=1
          l_end=.false.
          DO WHILE ( (j <= len(text)).AND.
     &      (i_equal == 0).AND.(.NOT.l_end) )
            IF (text(j:j) == '=') THEN
              i_equal=j
            ELSE IF (text(j:min(len(text), j+1)) == '//') THEN
              l_end=.true.
            ENDIF
            j=j+1
          ENDDO
        ENDDO
!
        CALL remove_blank(text(1:i_equal-1), name_temp, length_temp)
        i=1
        l_loop=(i <= n_dimension)
        IF (l_loop) l_loop=( l_loop.AND.(name_temp(1:length_temp) /= 
     &    dimension_name(i)(1:length_temp)) )
        DO WHILE (l_loop)
          i=i+1
          l_loop=(i <= n_dimension)
          IF (l_loop) l_loop=( l_loop.AND.(name_temp(1:length_temp) /= 
     &      dimension_name(i)(1:length_temp)) )
        ENDDO
        IF (i > n_dimension) THEN
          iv=1
          l_loop=(iv <= n_var)
          IF (l_loop) l_loop=( l_loop.AND.(name_temp(1:length_temp) /= 
     &      var_name(iv)(1:length_temp)) )
          DO WHILE (l_loop)
            iv=iv+1
            l_loop=(iv <= n_var)
            IF (l_loop) l_loop=( l_loop.AND.
     &        (name_temp(1:length_temp) /= 
     &        var_name(iv)(1:length_temp)) )
          ENDDO
          IF (iv > n_var) THEN
            WRITE(iu_err, '(/a)') 
     &        '*** Error: the variable'
     &        , name_temp(1:length_temp)
     &        , 'is not recognized in the file'
     &        , file_in(1:length_in)
            ierr=i_err_io
            RETURN
          ENDIF
        ENDIF
!
!       Shift the text to start at the numeric data and overwrite the
!       end with spaces.
        j=1
        DO WHILE (j <= len(text)-i_equal) 
          text(j:j)=text(j+i_equal:j+i_equal)
          j=j+1
        ENDDO
        DO WHILE (j <= len(text))
          text(j:j)=' '
          j=j+1
        ENDDO
!
!
        IF (i > n_dimension) THEN
!
!         This will be a variable indicated by the number IV set above.
!
          n_data(iv)=1
          DO j=1, n_dimension_var(iv)
            n_data(iv)
     &        =n_data(iv)*dimension_size(list_dimension_var(j, iv))
          ENDDO
!
          l_end=.false.
          i_start=1
          DO WHILE (.NOT.l_end)
            CALL split_cdl_line(ierr, text, var_type(iv), nd_cdl_data
     &        , data_int(1, iv), data_fl(1, iv), i_start, l_end, n_read)
            IF (ierr /= i_normal) THEN
              WRITE(iu_err, '(2(/a))')
     &          '*** Error: invalid data in file'
     &          , file_in(1:length_in)
              ierr=i_err_io
              goto 999
            ENDIF
            i_start=i_start+n_read
            READ(iunit, '(a)', iostat=ios1) text
          ENDDO
!
        ELSE
!
          l_end=.false.
          i_start=1
          DO WHILE (.NOT.l_end)
            CALL split_cdl_line(ierr, text, dimension_type(i)
     &        , dimension_size(i), dimension_array_int(1, i)
     &        , dimension_array_fl(1, i), i_start, l_end, n_read)
            IF (ierr /= i_normal) THEN
              WRITE(iu_err, '(2(/a))')
     &          '*** Error: invalid data in file'
     &          , file_in(1:length_in)
              ierr=i_err_io
              goto 999
            ENDIF
            i_start=i_start+n_read
            READ(iunit, '(a)', iostat=ios1) text
          ENDDO
!
        ENDIF
!
      ENDDO
!
!
!
999   continue
      CLOSE(iunit)
!
!
!
      RETURN
      END
