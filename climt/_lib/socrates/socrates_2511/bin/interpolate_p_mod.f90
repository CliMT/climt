! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE interpolate_p_mod

!   This module defines modes of interpolation. Broadly speaking,
!   these may be classified according to the use of logarithms
!   for the dependent or independent variables and whether
!   cubic splines are preferred to linear interpolation.
!
!   Note by way of example that interpolation of temperature in the
!   logarithm of the pressure and of logarithms of mixing ratios
!   in logarithms of the pressure generally work best.

IMPLICIT NONE

! Interpolation modes:
  INTEGER, Parameter :: IP_1_lin_lin  = 1
!   Linear Direct-Direct interpolation
  INTEGER, Parameter :: IP_1_log_lin  = 2
!   Linear Logarithmic-Direct interpolation
  INTEGER, Parameter :: IP_1_lin_log  = 3
!   Linear Direct-Logarithmic interpolation
  INTEGER, Parameter :: IP_1_log_log  = 4
!   Linear Logarithmic-Logarithmic interpolation

  INTEGER, Parameter :: IP_3_lin_lin  = 5
!   Direct-Direct interpolation with cubic splines
  INTEGER, Parameter :: IP_3_log_lin  = 6
!   Logarithmic-Direct interpolation with cubic splines
  INTEGER, Parameter :: IP_3_lin_log  = 7
!   Direct-Logarithmic interpolation with cubic splines
  INTEGER, Parameter :: IP_3_log_log  = 8
!   Logarithmic-Logarithmic interpolation with cubic splines

CONTAINS

SUBROUTINE interpolate_p(n, p, a, x, y, y2, pp, aa, i_mode, l_splined)

! Subroutine to interpolate to the value of a field at a point.
!
! Method:
!   The field to be evaluated is passed to the routine, together
!   with its second derivatives if a splined fit is used. Its
!   value at the pressure supplied is evaluated according to the
!   type of interpolation prescribed. The result is returned.

! Modules to set types of variables:
  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_normal, i_err_fatal, i_err_range
  USE ereport_mod, ONLY: ereport
  USE spline_evaluate_mod, ONLY: spline_evaluate
  USE spline_fit_mod, ONLY: spline_fit

  IMPLICIT NONE

! Dummy arguments.
  INTEGER, INTENT(IN)         :: n
!   Number of levels
  INTEGER, INTENT(IN)         :: i_mode
!   Mode of interpolation
  LOGICAL, INTENT(INOUT)      :: l_splined
!   True if field already splined
  REAL (RealK), INTENT(IN)    :: p(n)
!   Pressure levels
  REAL (RealK), INTENT(IN)    :: a(n)
!   Field to be interpolated
  REAL (RealK), INTENT(IN)    :: pp
!   Interpolating pressure
  REAL (RealK), INTENT(INOUT) :: x(n)
!   Converted abscissa
  REAL (RealK), INTENT(INOUT) :: y(n)
!   Converted ordinate
  REAL (RealK), INTENT(INOUT) :: y2(n)
!   Second derivative of ordinate
  REAL (RealK), INTENT(OUT)   :: aa
!   Interpolated value

! Local variables.
  INTEGER                     :: i
!   Loop variable
  REAL (RealK)                :: xx
!   Converted evalaution point
  REAL (RealK)                :: yy
!   Raw interpolate

  INTEGER                       :: ierr = i_normal
  CHARACTER (LEN = 80)          :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'INTERPOLATE_P'


! Perform the initial spline fit on the first call.
  IF (.NOT. l_splined) THEN
    IF (     (i_mode == IP_1_lin_lin)                                         &
        .OR. (i_mode == IP_3_lin_lin) ) THEN
!     Linear-linear fit.
      DO i = 1, n
        x(i) = p(i)
        y(i) = a(i)
      END DO
    ELSE IF (     (i_mode == IP_1_log_lin)                                    &
             .OR. (i_mode == IP_3_log_lin) ) THEN
!     Logarithmic-linear fit.
      DO i = 1, n
        x(i) = log(p(i))
        y(i) = a(i)
      END DO
    ELSE IF (     (i_mode == IP_1_lin_log)                                    &
             .OR. (i_mode == IP_3_lin_log) ) THEN
!     Linear-logarithmic fit.
      DO i = 1, n
        x(i) = p(i)
        y(i) = log(a(i))
      END DO
    ELSE IF (     (i_mode == IP_1_log_log)                                    &
             .OR. (i_mode == IP_3_log_log) ) THEN
!     Logarithmic-logarithmic fit.
      DO i = 1, n
        x(i) = log(p(i))
        y(i) = log(a(i))
      END DO
    ELSE
      cmessage = '*** Error: Unrecognized fitting mode.'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    END IF

    IF (     (i_mode == IP_1_lin_lin)                                         &
        .OR. (i_mode == IP_1_log_lin)                                         &
        .OR. (i_mode == IP_1_lin_log)                                         &
        .OR. (i_mode == IP_1_log_log) ) THEN
!     This is a spline fit without second derivatives.
      DO i = 1, n
        y2(i) = 0.0_RealK
      END DO
    ELSE IF (     (i_mode == IP_3_lin_lin)                                    &
             .OR. (i_mode == IP_3_log_lin)                                    &
             .OR. (i_mode == IP_3_lin_log)                                    &
             .OR. (i_mode == IP_3_log_log) ) THEN
!     Cubic splines.
      CALL spline_fit(n, x, y, y2)
    END IF
    l_splined = .true.

  END IF

! Calculate the splining point from PP.
  IF (     (i_mode == IP_1_lin_lin)                                           &
      .OR. (i_mode == IP_3_lin_lin) ) THEN
    xx = pp
  ELSE IF (     (i_mode == IP_1_log_lin)                                      &
           .OR. (i_mode == IP_3_log_lin) ) THEN
    xx = log(pp)
  ELSE IF (     (i_mode == IP_1_lin_log)                                      &
           .OR. (i_mode == IP_3_lin_log) ) THEN
    xx = pp
  ELSE IF (     (i_mode == IP_1_log_log)                                      &
           .OR. (i_mode == IP_3_log_log) ) THEN
    xx = log(pp)
  ELSE
    cmessage = '*** Error: Unrecognized fitting mode.'
    ierr = i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)
  END IF

  CALL spline_evaluate(ierr, n, x, y, y2, xx, yy)
  IF (ierr /= i_normal) THEN
    IF (ierr == i_err_range) THEN
      IF (xx < x(1)) THEN
        yy = y(1)
      ELSE IF (xx > x(n)) THEN
        yy = y(n)
      END IF
!     Recover form the error
      ierr=i_normal
    ELSE
      cmessage = '*** Error: Unable to recover from spline error.'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    END IF
  END IF

! Convert back to the actual data field.
  IF (     (i_mode == IP_1_lin_lin)                                           &
      .OR. (i_mode == IP_3_lin_lin) ) THEN
    aa = yy
  ELSE IF (     (i_mode == IP_1_log_lin)                                      &
           .OR. (i_mode == IP_3_log_lin) ) THEN
    aa = yy
  ELSE IF (     (i_mode == IP_1_lin_log)                                      &
           .OR. (i_mode == IP_3_lin_log) ) THEN
    aa = exp(yy)
  ELSE IF (     (i_mode == IP_1_log_log)                                      &
           .OR. (i_mode == IP_3_log_log) ) THEN
    aa = exp(yy)
  END IF

END SUBROUTINE interpolate_p

END MODULE interpolate_p_mod
