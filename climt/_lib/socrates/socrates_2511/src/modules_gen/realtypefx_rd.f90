! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set the precision of fixed types of REALs.
! This module is used for I/O in specific formats and the kinds
! should not be changed.

MODULE realtypefx_rd 

  IMPLICIT NONE

  INTEGER, PARAMETER :: Real4=SELECTED_REAL_KIND(6, 25)
  INTEGER, PARAMETER :: Real8=SELECTED_REAL_KIND(15, 307)

END MODULE realtypefx_rd
