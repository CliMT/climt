! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Definitions of the scattering database
!
MODULE def_sct_db
!
! Description:
!
! This module defines the record format used by the databases of
! scattering propeties of ice crystals.
!
!- End of header
!
!
!
  USE realtypefx_rd
!
!
  IMPLICIT NONE
!
!
!
  INTEGER, PARAMETER  ::  sct_db_recl = 512
!   Record length used in the scattering database
!
  INTEGER, PARAMETER  ::  npd_sct_db_angle = 1801
!   Size allocated for the angles of the phase function
!
! All units are to be SI.
  TYPE str_sctdb_entry
    REAL (Real4)  ::  dm
!     Mean maximum dimension of the crystal
    REAL (Real4)  ::  wavelength
!     Wavelength of incident radiation
    REAL (Real4)  ::  csca
!     Scattering cross-section
    REAL (Real4)  ::  cext
!     Extinction cross-section
    REAL (Real4)  ::  omega
!     Albedo of single scattering
    REAL (Real4)  ::  asymm
!     Asymmetry of the phase function
    INTEGER       ::  n_angle
!     Number of angles at which the phase function is specified
    LOGICAL       ::  l_uniform
!     Flag for uniform spacing of the data
    REAL (Real4)  ::  d_angle
!     Angular increment with even spacing
    REAL (Real4)  ::  phf_angle(npd_sct_db_angle)
!     Angles at which the phase function is specified
    REAL (Real4)  ::  phf(npd_sct_db_angle)
!     The phase function at the specified angles
!
  END TYPE str_sctdb_entry
!
!
END MODULE def_sct_db
