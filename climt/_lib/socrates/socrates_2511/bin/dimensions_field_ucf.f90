! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module setting the dimensions of physical fields

MODULE dimensions_field_ucf

  IMPLICIT NONE

  INTEGER, Parameter :: npd_latitude = 3
!       Number of latitudes
  INTEGER, Parameter :: npd_longitude = 1
!       Number of longitudes
  INTEGER, Parameter :: npd_profile = 3
!       Number of atmospheric profiles
  INTEGER, Parameter :: npd_layer = 170
!       Number of atmospheric layers
  INTEGER, PARAMETER :: npd_channel = 2
!       Number of spectral channels permitted for output
  INTEGER, Parameter :: npd_column = 6
!       Maximum number of cloudy subcolumns
  INTEGER, Parameter :: npd_direction = 63
!       Maximum number of directions for radiances
  INTEGER, Parameter :: npd_max_order = 101
!       Maximum order of spherical harmonics used
  INTEGER, PARAMETER :: npd_brdf_basis_fnc = 2
!       Number of BRDF basis functions
  INTEGER, PARAMETER :: npd_brdf_trunc = 5
!       Order of BRDF truncation
  INTEGER :: npd_profile_aerosol_prsc = 9
!       Size allocated for profiles of prescribed
!       cloudy optical properties
  INTEGER :: npd_profile_cloud_prsc = 9
!       Size allocated for profiles of prescribed
!       aerosol optical properties
  INTEGER :: npd_opt_level_aerosol_prsc = 170
!       Size allocated for levels of prescribed
!       cloudy optical properties
  INTEGER :: npd_opt_level_cloud_prsc = 170
!       Size allocated for levels of prescribed
!       aerosol optical properties

END MODULE dimensions_field_ucf
