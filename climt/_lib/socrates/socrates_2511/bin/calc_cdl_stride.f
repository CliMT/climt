! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the stride in a CDL-array.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      FUNCTION calc_cdl_stride(n_dimension, list, id
     &  , dimension_size, nd_cdl_dimen
     &  )
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
!     Declaration of variables:
!
!     Sizes of arrays
      INTEGER, Intent(IN) ::
     &    nd_cdl_dimen
!           Allowed size for CDL-dimensions
!
      INTEGER, Intent(IN) ::
     &    n_dimension
!           Number of dimensions in the field
     &  , list(nd_cdl_dimen)
!           List of CDL-dimensions
     &  , id
!           Index of the dimension
     &  , dimension_size(nd_cdl_dimen)
!           SIzes of the CDL-dimensions
!
      INTEGER   !,Intent(OUT)
     &    calc_cdl_stride
!            Calculated stride
!
!
!     Local Variables
!
      INTEGER
     &    i
!           Loop variable
      LOGICAL
     &    l_loop
!           Logical flag controlling the WHILE-loop
!
!
!
      i=n_dimension
      calc_cdl_stride=1
      l_loop=(list(i) /= id)
      DO WHILE (l_loop)
        calc_cdl_stride=calc_cdl_stride*dimension_size(list(i))
        i=i-1
        l_loop=(i > 0)
        IF (l_loop) l_loop=(list(i) /= id)
      ENDDO
!
!
!
      RETURN
      END
