! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE parkind1

!
! Description:
!   Dummy module to replace the DrHook library. Defines data types
!   which would otherwise be declared by the DrHook library.
!
! Owner: UM System team
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.


  IMPLICIT NONE

  INTEGER, PARAMETER :: jpim = SELECTED_INT_KIND(9) 
  INTEGER, PARAMETER :: jprb = SELECTED_REAL_KIND(13,300)

END MODULE PARKIND1
