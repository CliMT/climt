! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to perform a parabolic fit.
!
SUBROUTINE fit_parabola_90 &
!
(s1, s2, s3, &
 f1, f2, f3, &
 s0, l_fit)
!
! Description:
!   Three values, F1, F2 and F3 are supplied at points S1, S2 and
!   S3. A parabola, Q0+Q1*([S]-S0)**2 is fitted.
!
!
! Modules used:
  USE realtype_rd
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments.
  REAL  (RealK), Intent(IN) :: s1
!   Data point
  REAL  (RealK), Intent(IN) :: s2
!   Data point
  REAL  (RealK), Intent(IN) :: s3
!   Data point
  REAL  (RealK), Intent(IN) :: f1
!   Data value
  REAL  (RealK), Intent(IN) :: f2
!   Data value
  REAL  (RealK), Intent(IN) :: f3
!   Data value
  REAL  (RealK), Intent(OUT) :: s0
!   Position of minimum
  LOGICAL, Intent(OUT) :: l_fit
!   Flag for successful fit
!
! Local variables.
  REAL  (RealK) :: a
!   Algebraic term
  REAL  (RealK) :: b
!   Algebraic term
  REAL  (RealK) :: c
!   Algebraic term
  REAL  (RealK) :: d
!   Algebraic term
!
!
!
! Locate the turning point.
  a=(f1-f2)*(s2-s3)
  b=a*(s2+s3)
  c=(f3-f2)*(s2-s1)
  d=c*(s2+s1)
  IF (abs(a-c) > 0.0_RealK) THEN
    s0=0.5_RealK*(b-d)/(a-c)
!   Check that it is a minimum.
    IF (abs(f1-f2) > 0.0_RealK) THEN
      IF (((f1-f2)*(s1-s2)*(s1+s2-2.0_RealK*s0)) > 0.0_RealK) THEN
        l_fit = .TRUE. 
      ELSE
        l_fit = .FALSE. 
      ENDIF
    ELSE
      IF (((f2-f3)*(s2-s3)*(s2+s3-2.0_RealK*s0)) > 0.0_RealK) THEN
        l_fit = .TRUE. 
      ELSE
        l_fit = .FALSE. 
      ENDIF
    ENDIF
  ELSE
    l_fit = .FALSE. 
!   Logical tests elsewhere may involve both l_fit and s0. On some
!   machines it is important to ensure that s0 has a legal value
!   even if l_fit is .FALSE.
    s0 = 0.0_RealK
  ENDIF
!
!
!
  RETURN
END
