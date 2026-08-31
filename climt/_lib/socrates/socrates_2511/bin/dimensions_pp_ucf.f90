! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set dimensions for preprocessing routines.

MODULE dimensions_pp_ucf

  IMPLICIT NONE

  INTEGER, Parameter ::  npd_pt = 2000
!   Size allocated for p and T pairs
  INTEGER, Parameter ::  npd_amount = 330
!   Size allocated for amounts
  INTEGER, Parameter ::  npd_partial_p = 11
!   Size allocated for partial pressures
  INTEGER, Parameter ::  npd_amount_contin =  npd_amount / npd_partial_p
!   Size allocated for continuum amounts
  INTEGER, Parameter ::  npd_wavenumber_trans = 10000
!   Size allocated for wavenumbers in a band
  INTEGER, Parameter ::  npd_scatt_angle = 1001
!   Size allocated for scattering angles
  INTEGER, Parameter ::  npd_size_scat   = 500
!   Size allocated for snumber of different sizes in the database
  INTEGER, Parameter ::  npd_log_deriv = 250000
!   Size allocated for size of array of logarithmic 
!   derivatives in the Mie scattering code
  INTEGER, Parameter ::  npd_refract = 600
!   Size allocated for (complex) refractive indices 
  INTEGER, Parameter ::  npd_wavelength_scat = 500
!   Size allocated for scattering wavelengths
  INTEGER, Parameter ::  npd_mie_block = 1500
!   Size allocated for blocks of Mie scattering data
  INTEGER, Parameter ::  npd_thermal_abscissa = 3000000
!   Size allocated for number of points for integration with
!   respect to temperature
  INTEGER, Parameter ::  npd_general_data = 100
!   Size allocated for array allowed for general data
  INTEGER, Parameter ::  npd_gauss_point = 51
!   Size allocated for points of Gaussian integration
  INTEGER, Parameter ::  npd_phi_point = 51
!   Size allocated for integration in the phi-direction
  INTEGER, Parameter ::  npd_brdf_order = 51
!   Size allocated for orders of BRDFs
  INTEGER, Parameter ::  npd_size=5000
!   Size allocated for elements of a size distribution

END MODULE dimensions_pp_ucf
