! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to define identifers for SRA aerosol models.
!
MODULE aerosol_model_pcf
!
! Description:
!
!   This module defines identifiers for the aerosol models used
!   in the SRA climatological models set out in WCRP 55. The
!   models are conbinations of aerosols, individual components.
!
!   IMPORTANT: The initialization of fraction_component depends
!   on the ordering of the components. I don't know how to do this
!   better in F90.
!
!- End of header
!
!
! Modules used
  USE realtype_rd
  USE dimensions_spec_ucf
  USE rad_pcf
!
!
!
  IMPLICIT NONE
!
!
  INTEGER, Parameter :: NPD_aerosol_model   = 6
!   Size allocated for models
!
  INTEGER, Parameter :: I_urban_industrial  = 1
!   Urban-industrial model
  INTEGER, Parameter :: I_continental       = 2
!   Continental model
  INTEGER, Parameter :: I_maritime          = 3
!   Maritime model
  INTEGER, Parameter :: I_stratospheric     = 4
!   Stratospheric model
  INTEGER, Parameter :: I_volcanic_ash      = 5
!   Volcanic ash model
  INTEGER, Parameter :: I_volcanic_acid     = 6
!   Volcanic acid model
!
  INTEGER :: im
!   Variable for implicit DO-loop
!
  REAL  (RealK), Parameter, &
    Dimension(NPD_aerosol_component, NPD_aerosol_model) :: &
    fraction_component = RESHAPE( (/ &
    0.61_RealK, 0.17_RealK, 0.0_RealK, 0.22_RealK, 0.0_RealK, 0.0_RealK, &
      (0.0_RealK, im=1,NPD_aerosol_component-6), &
    0.29_RealK, 0.70_RealK, 0.0_RealK, 0.01_RealK, 0.0_RealK, 0.0_RealK, &
      (0.0_RealK, im=1,NPD_aerosol_component-6), &
    0.05_RealK, 0.0_RealK, 0.95_RealK, 0.0_RealK, 0.0_RealK, 0.0_RealK, &
      (0.0_RealK, im=1,NPD_aerosol_component-6), &
    0.0_RealK, 0.0_RealK, 0.0_RealK, 0.0_RealK, 0.0_RealK, 1.0_RealK, &
      (0.0_RealK, im=1,NPD_aerosol_component-6), &
    0.0_RealK, 0.0_RealK, 0.0_RealK, 0.0_RealK, 1.0_RealK, 0.0_RealK, &
      (0.0_RealK, im=1,NPD_aerosol_component-6), &
    0.0_RealK, 0.0_RealK, 0.0_RealK, 0.0_RealK, 0.0_RealK, 1.0_RealK, &
      (0.0_RealK, im=1,NPD_aerosol_component-6) /), &
    (/ NPD_aerosol_component, NPD_aerosol_model /) )
!
!
!
END MODULE aerosol_model_pcf
