! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module setting maximum dimensions for netCDF-files.

MODULE dimensions_cdf_ucf

  IMPLICIT NONE

  INTEGER, Parameter :: npd_cdl_dimen       = 5
!   Maximum number of netCDF dimensions
  INTEGER            :: npd_cdl_dimen_size  = 502
!   Maximum size of a dimension
  INTEGER            :: npd_cdl_data        = 10000
!   Maximum size of netCDF data
  INTEGER, Parameter :: npd_cdl_var         = 5
!   Maximum number of netCDF variables in a file

END MODULE dimensions_cdf_ucf
