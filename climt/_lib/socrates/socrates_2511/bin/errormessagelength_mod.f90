! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Module for allocation of error message length.

MODULE errormessagelength_mod

USE filenamelength_mod, ONLY: filenamelength

IMPLICIT NONE

! Description:
!   The module sets the max length of several filename-related strings.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

INTEGER, PARAMETER      :: errormessagelength = 2*filenamelength

!-------------------------------------------------------------------

END MODULE errormessagelength_mod
