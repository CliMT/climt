! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to declare a structure for refractive index data.

MODULE def_refract

  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE


  TYPE StrRefract

    INTEGER :: n_points = 0
!     Number of wavelength points
    REAL (RealK), ALLOCATABLE :: wavelength(:)
!     Wavelength at which the refractive index is specified
    REAL (RealK), ALLOCATABLE :: re_part(:)
!     Real part of refractive index
    REAL (RealK), ALLOCATABLE :: im_part(:)
!     Imaginary part of refractive index

  END TYPE StrRefract

CONTAINS


SUBROUTINE allocate_refract(refract)

  IMPLICIT NONE

  TYPE(StrRefract), INTENT(INOUT) :: refract

  ALLOCATE(refract%wavelength(refract%n_points))
  ALLOCATE(refract%re_part(refract%n_points))
  ALLOCATE(refract%im_part(refract%n_points))

END SUBROUTINE allocate_refract


SUBROUTINE deallocate_refract(refract)

  IMPLICIT NONE

  TYPE(StrRefract), INTENT(INOUT) :: refract

  DEALLOCATE(refract%wavelength)
  DEALLOCATE(refract%re_part)
  DEALLOCATE(refract%im_part)

END SUBROUTINE deallocate_refract

END MODULE def_refract
