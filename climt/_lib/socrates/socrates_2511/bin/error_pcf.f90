! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set error flags in the radiation code.

MODULE error_pcf

  IMPLICIT NONE

  INTEGER, PARAMETER  ::  i_normal = 0
!                           Error free condition
  INTEGER, PARAMETER  ::  i_err_fatal = 1
!                           Fatal error: immediate return
  INTEGER, PARAMETER  ::  i_abort_calculation = 2
!                           Calculation aborted
  INTEGER, PARAMETER  ::  i_missing_data = 3
!                           Missing data error: conditional
  INTEGER, PARAMETER  ::  i_err_io = 4
!                           I/O error
  INTEGER, PARAMETER  ::  i_err_range = 5
!                           Interpolation range error
  INTEGER, PARAMETER  ::  i_err_exist = 6
!                           Existence error
  INTEGER, PARAMETER  ::  i_warning = -1
!                           Non-fatal warning

END MODULE error_pcf
