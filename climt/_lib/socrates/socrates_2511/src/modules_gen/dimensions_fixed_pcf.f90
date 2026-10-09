! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set internal dimensions tied to algorithms, mostly for clouds.

MODULE dimensions_fixed_pcf

  IMPLICIT NONE

  INTEGER, Parameter :: npd_cloud_component        =  4
!   Number of components of clouds.
  INTEGER, Parameter :: npd_cloud_type             =  4
!   Number of permitted types of clouds.
  INTEGER, Parameter :: npd_overlap_coeff          = 18
!   Number of overlap coefficients for cloud
  INTEGER, Parameter :: npd_source_coeff           =  2
!   Number of coefficients for two-stream sources
  INTEGER, Parameter :: npd_region                 =  3
!   Number of regions in a layer

END MODULE dimensions_fixed_pcf
