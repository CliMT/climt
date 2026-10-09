! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure for interpolated quantum yields
MODULE def_qy

USE realtype_rd, ONLY: RealK

IMPLICIT NONE

TYPE StrQy
  REAL (RealK), ALLOCATABLE :: qy(:, :, :)
!   Quantum yield interpolated to model grid
END TYPE StrQy

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE allocate_qy(photol, n_profile, n_layer, n_wl)

IMPLICIT NONE

TYPE (StrQy), INTENT(INOUT) :: photol
INTEGER, INTENT(IN) :: n_profile, n_layer, n_wl

IF (.NOT. ALLOCATED(photol%qy)) ALLOCATE(photol%qy(n_profile, n_layer, n_wl))

END SUBROUTINE allocate_qy
!------------------------------------------------------------------------------
SUBROUTINE deallocate_qy(photol)

IMPLICIT NONE

TYPE (StrQy), INTENT(INOUT) :: photol

IF (ALLOCATED(photol%qy)) DEALLOCATE(photol%qy)

END SUBROUTINE deallocate_qy

END MODULE def_qy
