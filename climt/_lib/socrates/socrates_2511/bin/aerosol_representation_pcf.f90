! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to representations for aerosols.
!
MODULE aerosol_representation_pcf
!
! Description:
!
!   This module defines representations of the variation of aerosol
!   mixing ratios with height as used in the SRA atmospheres.
!
!- End of header
!
!
  IMPLICIT NONE
!
!
  INTEGER, Parameter :: I_null         = 0
!   Null initializer
  INTEGER, Parameter :: I_constant     = 1
!   No variation with height
  INTEGER, Parameter :: I_linear       = 2
!   Linrear variation with height
  INTEGER, Parameter :: I_exponential  = 3
!   Exponential decay with height
!
!
!
END MODULE aerosol_representation_pcf
