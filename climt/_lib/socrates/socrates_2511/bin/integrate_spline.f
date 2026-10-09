! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to evaluate the integral of a cubic spline.
!
! Method:
!   The range of integration is checked against the range of the
!   spline. If the range of integration lies within the range of 
!   the spline the integral is evaluated noting that between each
!   pair of knots the spline is a cubic which can be integrated
!   directly.
!
!- ---------------------------------------------------------------------
      SUBROUTINE integrate_spline(ierr, x1, x2
     &  , n_data, x, y, y2, integral)
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE rad_pcf
      USE def_std_io_icf

      IMPLICIT NONE


!     Dummy arguments.
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
      INTEGER, Intent(IN) ::
     &    n_data
!           Number of data points
      REAL  (RealK), Intent(IN) ::
     &    x1
!           Lower limit of integration
     &  , x2
!           Upper limit of integration
     &  , x(n_data)
!           Spline knotting points
     &  , y(n_data)
!           Values of spline
     &  , y2(n_data)
!           Second derivative
      REAL  (RealK), Intent(OUT) ::
     &    integral
!           Evaluated integral
!
!     Local variables.
      INTEGER
     &    i_low
!           Lower knot point
     &  , i_high
!           Upper knot point
     &  , i_mid
!           Bisecting point
      REAL  (RealK) ::
     &    xl
!           Lower limit of integral
     &  , xh
!           Upper limit of integral
     &  , x_bar
!           Mid-point of range of integration
     &  , delta
!           Width of interval
     &  , delta_1
!           Width of integration
     &  , ia
!           Integral of coefficient
     &  , ib
!           Integral of coefficient
     &  , i3a
!           Integral of coefficient
     &  , i3b
!           Integral of coefficient
!
!
      integral=0.0_RealK
!
!     An integral with one data point cannot be performed.
      IF (n_data == 1) THEN
        WRITE(iu_err, '(/a)') 
     &    '*** error: spline has only one knot point.'
        ierr=i_abort_calculation
        RETURN
      ENDIF
!
!     Check that the integrand lies within the domain of the spline
!     and find the first interval.
      IF ( (x1 > x(n_data)).OR.(x2 < x(1)) ) THEN
!       Domain of integrand lies outside range of spline.
        ierr=i_err_range
        RETURN
      ELSE IF (x1 <= x(1)) THEN
        xl=x(1)
        xh=min(x2, x(2))
        i_low=1
      ELSE
!       Find the interval in which x1 lies.
        i_low=1
        i_high=n_data
1         i_mid=(i_low+i_high)/2
          IF (x1 > x(i_mid)) THEN
            i_low=i_mid
          ELSE
            i_high=i_mid
          ENDIF
          IF ( (i_high-i_low) > 1) goto 1
!
        xl=x1
        xh=min(x(i_high), x2)
      ENDIF
      i_high=i_low+1
!
2       delta=(x(i_high)-x(i_low))
        delta_1=xh-xl
        x_bar=0.5_RealK*(xl+xh)
        ia=(delta_1/delta)*(x(i_high)-x_bar)
        ib=(delta_1/delta)*(x_bar-x(i_low))
        i3a=0.25_RealK*((x(i_high)-xl)**4
     &    -(x(i_high)-xh)**4)/delta**3
        i3b=0.25_RealK*((x(i_low)-xh)**4-(x(i_low)-xl)**4)/delta**3
        integral=integral+y(i_low)*ia+y(i_high)*ib
     &    +y2(i_low)*delta**2*(i3a-ia)/6.0_RealK
     &    +y2(i_high)*delta**2*(i3b-ib)/6.0_RealK
!       Determine next interval.
        IF (xh < x2) THEN
          i_low=i_low+1
          i_high=i_low+1
          IF (i_high > n_data) THEN
!           Range of spline is exhausted, so stop adding terms.
            RETURN
          ENDIF
          xl=x(i_low)
          xh=min(x2, x(i_high))
          goto 2
        ENDIF
!
!
!
      RETURN
      END
