! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find and assign a dimeniosn in a CDL-file.
!
! Purpose:
!   This subroutine searches for a dimension with a given name
!   in a CDL-file.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE find_dimen_cdl(ierr
     &  , file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , dim_find, len_dim, type_find, len_type
     &  , l_assigned, n_assigned, l_assign, l_required
     &  , id, array_int, array_fl
     &  , nd_array
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
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
!     Dimensions of arrays:
!
      INTEGER, Intent(IN) ::
     &    nd_array
!           Allowed size retrieved arrays
!
!
!     INCLUDE HEADER FILES.
!
!
!
!     Dummy arguments:
!
      INTEGER   !,Intent(OUT)
     &    ierr
!            Error flag
!
      INCLUDE 'cdl_struc.finc'
!
      CHARACTER !, Intent(IN)
     &    file_name*(*)
!           Name of input file
     &  , dim_find*(*)
!           Name of dimension to find
     &  , type_find*(*)
!           Type fo dimension to find
      INTEGER, Intent(IN) ::
     &    len_dim
!           Length of name of dimension to find
     &  , len_type
!           Length of name of type of dimension
!
      LOGICAL, Intent(IN) ::
     &    l_assigned
!           Flag for pre-assigned variable
     &  , l_assign
!           Flag to assigned variable
     &  , l_required
!           Flag to invoke test for a mandatory variable
!
!
      INTEGER, Intent(INOUT) ::
     &    n_assigned
!           Size assigned to the dimension
      INTEGER, Intent(OUT) ::
     &    id
!           CDL `index' of dimension
!
      INTEGER, Intent(INOUT) ::
     &    array_int(nd_array)
!           Integer array assigned
      REAL  (RealK), Intent(INOUT) ::
     &    array_fl(nd_array)
!           Floating-point array assigned
!
!
!     Local Variables
!
      INTEGER
     &    l
!           Loop variable
      LOGICAL
     &    l_found
!           Flag for dimension found 
!
!
!
!
      l_found=.false.
      id=1
      DO WHILE ( (id <= n_dimension).AND.(.NOT.l_found) )
        IF (dimension_name(id)(1:len_dim) == dim_find(1:len_dim)) THEN
!
          l_found=.true.
          IF (dimension_type(id)(1:len_type) /= 
     &      type_find(1:len_type)) THEN
            WRITE(iu_err, '(2(/a))')
     &        '*** Error: Illegal type for dimension '
     &        //dim_find(1:len_dim)//' while reading the file '
     &        , file_name
            ierr=i_err_fatal
            RETURN
          ENDIF
!         Check the size of the fields.
          IF (l_assigned) THEN
            IF (n_assigned /= dimension_size(id)) THEN
              WRITE(iu_err, '(2(/a))')
     &          '*** Error: The size of dimension '
     &          //dim_find(1:len_dim)
     &          //' is inconsistent with the existing value in '
     &          , file_name
              ierr=i_err_fatal
              RETURN
            ENDIF
          ENDIF
          IF ((.NOT.l_assigned).AND.l_assign) THEN
            IF (dimension_size(id) <= nd_array) THEN
              n_assigned=dimension_size(id)
            ELSE
              WRITE(iu_err, '(2(/a))')
     &          '*** Error: The size of dimension '//dim_find(1:len_dim)
     &          //' is too large in ', file_name
              ierr=i_err_fatal
              RETURN
            ENDIF
            IF (type_find(1:min(5, len_type)) == 'float') THEN
              DO l=1, n_assigned
                array_fl(l)=dimension_array_fl(l, id)
              ENDDO
            ELSE IF (type_find(1:min(3, len_type)) == 'int') THEN
              DO l=1, n_assigned
                array_int(l)=dimension_array_int(l, id)
              ENDDO
            ENDIF
          ENDIF
!
        ELSE
!
          id=id+1
!
        ENDIF
      ENDDO
!
!
      IF (l_required) THEN
        IF (.NOT.l_found) THEN
          WRITE(iu_err, '(2(/a))')
     &      '*** Error: The dimension '//dim_find(1:len_dim)
     &      //' was not found in the file ', file_name
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
!
!
!
      RETURN
      END
