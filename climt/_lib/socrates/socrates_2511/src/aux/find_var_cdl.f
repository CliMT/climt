! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find a variable in a CDL-file.
!
! Purpose:
!   This subroutine searches for a variable with a given name
!   in a CDL-file.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE find_var_cdl(ierr
     &  , file_name
     &  , n_var, var_name, var_type
     &  , var_find, len_var, type_find, len_type
     &  , l_required
     &  , iv
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
     &  , var_find*(*)
!           Name of variable to find
     &  , type_find*(*)
!           Type of variable to find
      INTEGER, Intent(IN) ::
     &    len_var
!           Length of name of variable to find
     &  , len_type
!           Length of name of type of variable
!
      LOGICAL, Intent(IN) ::
     &    l_required
!           Flag to invoke test for a mandatory variable
!
!
      INTEGER, Intent(OUT) ::
     &    iv
!           CDL `index' of variable
!
!
!     Local Variables
!
      LOGICAL
     &    l_found
!           Flag for dimension found 
!
!
!
!
      l_found=.false.
      iv=1
      DO WHILE ( (iv <= n_var).AND.(.NOT.l_found) )
        IF (var_name(iv)(1:len_var) == var_find(1:len_var)) THEN
!
          l_found=.true.
          IF (var_type(iv)(1:len_type) /= 
     &      type_find(1:len_type)) THEN
            WRITE(iu_err, '(2(/a))')
     &        '*** Error: Illegal type for variable '
     &        //var_find(1:len_var)//' while reading the file '
     &        , file_name
            ierr=i_err_fatal
            RETURN
          ENDIF
!
        ELSE
!
          iv=iv+1
!
        ENDIF
      ENDDO
!
!
      IF (l_required) THEN
        IF (.NOT.l_found) THEN
          WRITE(iu_err, '(2(/a))')
     &      '*** Error: The variable '//var_find(1:len_var)
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
