! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to define possible weighting schemes.

MODULE weighting_pcf

  IMPLICIT NONE

  INTEGER, Parameter :: IP_weight_planck = 1
!    Planckian weighting
  INTEGER, Parameter :: IP_weight_d_planck = 2
!    Differential Planckian weighting
  INTEGER, Parameter :: IP_weight_solar = 3
!    Solar weighting
  INTEGER, Parameter :: IP_weight_uniform = 4
!    Uniform weighting

END MODULE weighting_pcf
