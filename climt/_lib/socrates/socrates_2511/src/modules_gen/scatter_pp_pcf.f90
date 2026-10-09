! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set types of scattering object.

MODULE scatter_pp_pcf

! Description:
!
! This module defines identifiers for the types of atmospheric
! objects which scatter radiation.

  IMPLICIT NONE

  INTEGER, Parameter :: npd_scatter_type = 3
!           Number of tyeps of scatterer

  INTEGER, Parameter :: ip_type_unassigned = 0
!           Index of unassigned scatterer
  INTEGER, Parameter :: ip_type_aerosol = 1
!           Index of aerosol
  INTEGER, Parameter :: ip_type_droplet = 2
!           Index of droplets
  INTEGER, Parameter :: ip_type_ice = 3
!           Index of ice crystals

  CHARACTER(LEN=20), Parameter, Dimension(0: npd_scatter_type) :: &
    name_scatter_type = (/ &
      'Unassigned          ', &
      'Aerosol             ', &
      'Droplets            ', &
      'Ice crystals        ' /)

END MODULE scatter_pp_pcf
