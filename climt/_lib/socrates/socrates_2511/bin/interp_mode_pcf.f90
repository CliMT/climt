! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to define modes of interpolation.

MODULE interp_mode_pcf

! Description:
!
!   This module defines modes of interpolation. Broadly speaking,
!   these may be classified according to the use of logarithms
!   for the dependent or independent variables and whether
!   cubic splines are preferred to linear interpolation.
!
!   Note by way of example that interpolation of temperature in the
!   logarithm of the pressure and of logarithms of mixing ratios
!   in logarithms of the pressure generally work best.

  IMPLICIT NONE

  INTEGER, Parameter :: IP_1_lin_lin  = 1
!   Linear Direct-Direct interpolation
  INTEGER, Parameter :: IP_1_log_lin  = 2
!   Linear Logarithmic-Direct interpolation
  INTEGER, Parameter :: IP_1_lin_log  = 3
!   Linear Direct-Logarithmic interpolation
  INTEGER, Parameter :: IP_1_log_log  = 4
!   Linear Logarithmic-Logarithmic interpolation

  INTEGER, Parameter :: IP_3_lin_lin  = 5
!   Direct-Direct interpolation with cubic splines
  INTEGER, Parameter :: IP_3_log_lin  = 6
!   Logarithmic-Direct interpolation with cubic splines
  INTEGER, Parameter :: IP_3_lin_log  = 7
!   Direct-Logarithmic interpolation with cubic splines
  INTEGER, Parameter :: IP_3_log_log  = 8
!   Logarithmic-Logarithmic interpolation with cubic splines

END MODULE interp_mode_pcf
