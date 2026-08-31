! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set types for units used in input files.

MODULE unit_list_pcf

! Description:
!   This module defines the units for input files such as
! instrumental response functions.

  IMPLICIT NONE

  INTEGER, Parameter :: IP_unit_metre = 1
!   Identifier for metres.
  INTEGER, Parameter :: IP_unit_inverse_cm = 2
!   Identifier for inverses centimeters
  INTEGER, Parameter :: IP_unit_micron = 3
!   Identifier for microns
  INTEGER, Parameter :: IP_unit_nanometre = 4
!   Identifier for nanometres

END MODULE unit_list_pcf
