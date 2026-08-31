! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Diffusivity factors used in two-stream schemes

MODULE diffusivity_factor

  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE

! Diffusivity factor for elsasser's scheme
  REAL (RealK), PARAMETER :: elsasser_factor = 1.66e+00_RealK

! Diffusivity factor for use with equivalent extinction
  REAL (RealK), PARAMETER :: diffusivity_factor_minor = 1.66e+00_RealK

! Diffusivity factor for use with the contribution function
  REAL (RealK), PARAMETER :: diffusivity_factor_cf = 2.0_RealK

END MODULE diffusivity_factor
