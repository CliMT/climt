! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to declare a structure for single scattering properties.
!
MODULE def_s_scat_prop
!
! Description:
!
! This module contains the declaration of the structure
! for single scattering properties.
!
!- End of header
!
! Modules used:
  USE realtype_rd
!
!
!
  TYPE StrSingleScatSzD
!
    INTEGER :: i_scatter_type
!     Type of scattering particle
    CHARACTER  (LEN=20) :: name_type
!     Name of the type of scatterer
    INTEGER :: i_component
!     Index of component, if an aerosol
    CHARACTER  (LEN=20) :: name_component
!     Name of the component, if an aerosol

    INTEGER :: code
!     Code indicating which data are available
!
!   Elements from the size distribution:
    REAL  (RealK) :: mass_mixing_ratio
!     Mass mixing ratio of scatterers
    REAL  (RealK) :: air_density
!     Air density used for unit conversions
    REAL  (RealK) :: particle_density
!     Particle density used for unit conversions
    REAL  (RealK) :: number_density
!     Number density of scatterers
    REAL  (RealK) :: vol_frac
!     Volume fraction of scattering species
    REAL  (RealK) :: proj_area
!     Projected area of scattering species
    REAL  (RealK) :: dim_char
!     Characteristic dimension: this is a generic variable, which
!     may hold the effective radius, mean maximum dimension or some
!     other measure of size as appropriate
    REAL  (RealK) :: re
!     Effective radius
    REAL  (RealK) :: de
!     Effective dimension
    REAL  (RealK) :: dmm
!     Mean maximum dimension
!
!   Information for aerosols
    REAL  (RealK) :: humidity
!     Humidity at which moist aerosol properties
!     have been calculated
    REAL  (RealK) :: vol_frac_dry
!     Volume fraction of dry aerosol
    REAL  (RealK) :: re_dry
!     Effective radius of dry aerosol
!
!   Generic single scattering properties:
    INTEGER :: n_phf_term
!     Number of terms in the phase function
!
!   Spectrally dependent properties:
    INTEGER :: n_wavenumber
!     Number of wavelengths
    REAL  (RealK), Pointer :: wavenum(:)
!     Wavelengths of data (Unit: m^-1)
    REAL  (RealK), Pointer :: abs(:)
!     Absorption coefficients (Unit: m^-1)
    REAL  (RealK), Pointer :: scat(:)
!     Scattering coefficient (Unit: m^-1)
    REAL  (RealK), Pointer :: phf(:, :)
!     Moments of the phase function
!
  END TYPE StrSingleScatSzD
!
!
END MODULE def_s_scat_prop
