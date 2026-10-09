! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to generate a Gaussian random variable
!
FUNCTION rand_gauss &

(seed) &

RESULT (gauss)

! Description:
!   This routine generates a ransom variable with zero mean and unit
!   variance using the Box Muller method.

  USE realtype_rd

  IMPLICIT NONE


! Dummy arguments
  INTEGER, Intent(IN) :: seed
!   Seed for random number generator
  REAL  (RealK) :: gauss
!   Generated number

! Local variables
  INTEGER ::  start = 0
!   Starting value for generator
  REAL  (RealK) :: gset
  REAL  (RealK) :: fac
  REAL  (RealK) :: rsq
  REAL  (RealK) :: v1
  REAL  (RealK) :: v2
  REAL  (RealK) :: harvest

  SAVE start, gset


  SELECT CASE(start)
    CASE(0)
      DO
        call random_seed
        call random_number(harvest)
        v1=2.0*harvest-1.0
        call random_seed
        call random_number(harvest)
        v2=2.0*harvest-1.0
        rsq=v1**2+v2**2
        IF ( (rsq > 0) .AND. (rsq < 1) ) EXIT
      ENDDO
      fac=SQRT(-2.0*LOG(rsq)/rsq)
      gset=v1*fac
      gauss=v2*fac
      start=1
    CASE DEFAULT
      gauss=gset
      start=0
  END SELECT

  RETURN
END FUNCTION rand_gauss
