! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module setting sizes of spectral arrays.
!
! Description:
! This module contains the default sizes for spectral arrays. In the
! course of time it should become redundant.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------
MODULE dimensions_spec_ucf

IMPLICIT NONE

INTEGER, PARAMETER :: npd_band = 301
! Number of spectral bands
INTEGER, PARAMETER :: npd_exclude = 1
! Numer of excluded bands
INTEGER, PARAMETER :: npd_k_term = 150
! Number of esft terms
INTEGER, PARAMETER :: npd_type = 20
! Number of data types
INTEGER, PARAMETER :: npd_species = 11
! Number of gaseous species
INTEGER, PARAMETER :: npd_scale_variable = 10
! Number of scaling variables
INTEGER, PARAMETER :: npd_continuum = 4
! Number of continua
INTEGER, PARAMETER :: npd_drop_type = 6
! Number of drop types
INTEGER, PARAMETER :: npd_ice_type = 16
! Number of ice crystal types
INTEGER, PARAMETER :: npd_aerosol_species = 50
! Number of aerosol species
INTEGER, PARAMETER :: npd_thermal_coeff = 5001
! Number of thermal coefficients
INTEGER, PARAMETER :: npd_fit_temp = 5001
! Number of temperature datapoints
INTEGER, PARAMETER :: npd_cloud_parameter = 504
! Number of cloud parameters
INTEGER, PARAMETER :: npd_humidities = 21
! Number of humidities
INTEGER, PARAMETER :: npd_gas_frac = 21
! Number of gas fractions
INTEGER, PARAMETER :: npd_phase_term = 501
! Number of terms in the phase function
INTEGER, PARAMETER :: npd_aod_wavel = 6
! Number of wavelengths for aerosol optical depth

END MODULE dimensions_spec_ucf
