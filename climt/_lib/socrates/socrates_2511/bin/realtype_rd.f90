! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set the precision of real variables.

MODULE realtype_rd

  IMPLICIT NONE

  ! Internal Real precision within Socrates
  INTEGER, PARAMETER :: RealK=SELECTED_REAL_KIND(15, 307)

  ! External Real precision for variables passed through the Runes interface
  INTEGER, PARAMETER :: RealExt=SELECTED_REAL_KIND(15, 307)

END MODULE realtype_rd
