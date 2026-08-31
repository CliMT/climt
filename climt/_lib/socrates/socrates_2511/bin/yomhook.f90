! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE yomhook

!
! Description:
!   Dummy module to replace the DrHook library.
!
! Owner: UM System team
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.

  USE parkind1, ONLY: jpim, jprb
  IMPLICIT NONE

  LOGICAL, PARAMETER :: lhook = .FALSE.

CONTAINS

  SUBROUTINE dr_hook(name,code,handle)
    IMPLICIT NONE

    !Arguments
    CHARACTER(len=*),   INTENT(in)    :: name
    INTEGER(kind=jpim), INTENT(in)    :: code
    REAL(kind=jprb),    INTENT(inout) :: handle

    !Nothing to do

    RETURN
  END SUBROUTINE dr_hook

END MODULE yomhook
