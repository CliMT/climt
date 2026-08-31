! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates Voigt lineshape
!
SUBROUTINE voigt_profile ( &
     nu, line_centre, &                      ! Input arguments
     S_adj, alpha_lorentz, alpha_doppler, &  ! Input arguments
     kres)                                   ! Output argument

! Description:
!
! Calculates Voigt lineshape
!
! Method:
!
! Uses algorithm of 
! Schreier (1992), JQSRT, 48, 743
! based on formulation of
! Humlicek (1982), JQSRT, 27, 437
!

  USE realtype_rd
  
  IMPLICIT NONE

! Subroutine arguments

! Scalar arguments with intent(in):
  REAL  (RealK), Intent(IN) :: nu
  REAL  (RealK), Intent(IN) :: line_centre
  REAL  (RealK), Intent(IN) :: S_adj
  REAL  (RealK), Intent(IN) :: alpha_lorentz
  REAL  (RealK), Intent(IN) :: alpha_doppler

! Scalar arguments with intent(out):
  REAL  (RealK), Intent(OUT) :: kres

! Local parameters:

  REAL  (RealK), PARAMETER :: recip_sqrt_pi = 0.5641895
  REAL  (RealK), PARAMETER :: sqrt_ln2      = 0.83255461

! Local scalars:
  COMPLEX  (RealK) :: t, u
  COMPLEX  (RealK) :: prob_func

  REAL  (RealK) :: S0
  REAL  (RealK) :: x, y
  REAL  (RealK) :: s1, s2

!- End of header

  S0 = (sqrt_ln2/alpha_doppler) * (recip_sqrt_pi) * S_adj
  y  = alpha_lorentz * (sqrt_ln2/alpha_doppler)

  x  = (nu - line_centre) * (sqrt_ln2/alpha_doppler)

  s1 = ABS(x) + y
  s2 = 0.195*ABS(x) - 0.176
  t  = CMPLX(y,-x)

  IF (s1 >=  15.0) THEN
     ! Region I
     prob_func  = approx1(t)

  ELSE IF (s1 < 15.0 .AND. s1 >=  5.5) THEN
     ! Region II
     u = t*t
     prob_func  = approx2(t,u)

  ELSE IF (s1 < 5.5 .AND. y >= s2) THEN
     ! Region III
     prob_func  = approx3(t)

  ELSE
     ! Region IV
     u = t*t
     prob_func  = exp(u) - approx4(t,u)
  END IF

  kres = S0 * REAL(prob_func)

CONTAINS
  ! Humlicek formulations for four regions of Voigt profile

  COMPLEX (RealK) FUNCTION approx1(t) RESULT (z)
    ! Region I

    IMPLICIT NONE

    COMPLEX  (RealK) :: t

    z = (t*0.5641896) / (0.5 + (t*t))

    RETURN

  END FUNCTION approx1

  COMPLEX (RealK) FUNCTION approx2(t, u) RESULT (z)
    ! Region II

    IMPLICIT NONE

    COMPLEX  (RealK) :: t, u

    z = (t * (1.410474 + u*0.5641896)) / (0.75 + (u * (3.0 + u)))

    RETURN

  END FUNCTION approx2

  COMPLEX (RealK) FUNCTION approx3(t) RESULT (z)
    ! Region III

    IMPLICIT NONE

    COMPLEX  (RealK) :: t

    z = (16.4955 + t*(20.20933 + t*(11.96482 + t*(3.778987 + t*0.5642236))))/ &
        (16.4955 + t*(38.82363 + t*(39.27121 + t*(21.69274 + t*(6.699398 + t)))))
  
    RETURN

  END FUNCTION approx3

  COMPLEX (RealK) FUNCTION approx4(t, u) RESULT (z)
    ! Region IV

    IMPLICIT NONE

    COMPLEX  (RealK) :: t, u

    z = (t*(36183.31 - u*(3321.9905 - u*(1540.787 - u &
         *(219.0313 - u*(35.76683  - u*(1.320522 - u*0.56419))))))) / &
         (32066.6 - u*(24322.84 - u*(9022.228 - u*(2186.181 &
         - u*(364.2191 - u*(61.57037 - u*(1.841439 - u)))))))

    RETURN

  END FUNCTION approx4

END SUBROUTINE voigt_profile
