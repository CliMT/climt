! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate weightings across a contiguous band.
!
SUBROUTINE weightings_single_90(ierr, &
  i_weight, &
  SolarSpec, t, &
  include_instrument_response, filter, &
  nu_low, nu_high, &
  n_nu, nu, &
  i_begin, n_int_weight, weight_coeff, &
  integral_weight)
!
! Description:
!   The weighting coefficients are calculated and the integral
!   is formed.
!
!
!
!
! Modules to set types of variables:
  USE realtype_rd
  USE def_std_io_icf
  USE def_inst_flt
  USE def_solarspec
  USE dimensions_pp_ucf
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
!
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
!
  INTEGER, Intent(IN) :: i_weight
!   Weighting option
  INTEGER, Intent(IN) :: n_nu
!   Number of frequency points
!
  INTEGER, Intent(OUT) :: i_begin
!   Beginning of data range
  INTEGER, Intent(OUT) :: n_int_weight
!   Number of points of integral
!
!
  REAL  (RealK), Intent(IN) :: t
!   Temperature of black body
  REAL  (RealK), Intent(IN) :: nu_low
!   Lower frequency of band
  REAL  (RealK), Intent(IN) :: nu_high
!   Upper frequency of band
  REAL  (RealK), Intent(IN) :: nu(n_nu)
!   Frequencies of input data
!
  TYPE (StrSolarSpec), Intent(IN) :: SolarSpec
!   Solar spectrum
!
  LOGICAL, Intent(IN) :: include_instrument_response
!   Flag to include the instrumental response function
  TYPE  (StrFiltResp), Intent(IN) :: filter
!   Instrumental response function
!
  REAL  (RealK), Intent(OUT) :: weight_coeff(0: n_nu)
!   Weighting coefficients
  REAL  (RealK), Intent(OUT) :: integral_weight
!   Integral of weighting function
!
!
! Local variables
  INTEGER :: i
!   Loop variable
  INTEGER :: i_end
!   End of data range
!
  REAL  (RealK), Allocatable, Dimension(:) :: nu_wgt
!   Frequencies for evaluating weighting function
  REAL  (RealK), Allocatable, Dimension(:) :: wgt
!   Weighting function: note that this is the actual weighting
!   function as appears in a continuous integral, not the weights
!   applied to a set of quadrature points
!
  REAL  (RealK) :: fraction_lower
!   Fraction of lower extreme interval included
  REAL  (RealK) :: fraction_upper
!   Fraction of upper extreme interval included
  REAL  (RealK) :: width_lower
!   Width of lower interval
  REAL  (RealK) :: width_upper
!   Width of upper interval
  REAL  (RealK) :: weight_lower
!   Weighting for lower point
  REAL  (RealK) :: weight_upper
!   Weighting for upper point
!
!
! Functions called:
!
  INTERFACE
!
    SUBROUTINE rad_weight_90(i_weight, nu, SolarSpec, T, weight)
!     Function to calculate the array of radiant weightings
!
      USE def_solarspec
!
      INTEGER, Intent(IN)              :: i_weight
      REAL  (RealK), Intent(IN)        :: nu(:)
      TYPE  (StrSolarSpec), Intent(IN) :: SolarSpec
      REAL  (RealK), Intent(IN)        :: T
!
      REAL  (RealK), Intent(OUT)       :: weight(:)
!
    END SUBROUTINE rad_weight_90
!
!
  END INTERFACE
!
!
!
! Find the points in the array of frequencies which lie just
! outside the band. begin by finding points just within the band.
  i_begin=1
  DO
!   This loop must terminate because of the check in the calling
!   routine.
    IF (nu(i_begin) > nu_low) EXIT
    i_begin=i_begin+1
  ENDDO
  i_begin=i_begin-1
!
  i_end=n_nu
  DO
!   This loop must terminate because of the check in the calling
!   routine.
    IF (nu(i_end) < nu_high) EXIT
    i_end=i_end-1
  ENDDO
  i_end=i_end+1
!
!
!
! Any scheme for numerical quadrature will involve assigning
! weightings to be applied at the points where data are available.
! there will therefore be n_int_weight = i_end-i_begin intervals 
! to be weighted, so the weighting coefficients may be numbered 
! from 0 to n_int_weight.
!
  n_int_weight = i_end-i_begin
! The weighting scheme involves trapezoidal integration over the
! available intervals. Linear interpolation of the values to be
! weighted is assumed.
!
! Define the array of frequencies where weights are to be calculated.
  ALLOCATE(nu_wgt(0:n_int_weight))
  ALLOCATE(wgt(0:n_int_weight))
  nu_wgt(0)            = nu_low
  nu_wgt(n_int_weight) = nu_high
  nu_wgt(1:n_int_weight-1) = nu(i_begin+1:i_end-1)
  CALL rad_weight_90(i_weight, nu_wgt, SolarSpec, t, wgt)
  IF (include_instrument_response) CALL apply_response_int
  
  DEALLOCATE(nu_wgt)
!
! Now define the weighting coefficients themselves.
!
!
  SELECT CASE(n_int_weight)
    CASE (3:)
!
!     Normal case: 3 or more weighting intervals.
!
!     Two lower extreme points:
      width_lower=nu(i_begin+1)-nu(i_begin)
      fraction_lower=(nu(i_begin+1)-nu_low) /  width_lower
!     The squaring of fraction_lower is correct as one factor
!     comes from interpolation and one from the width of the
!     interval.
      weight_coeff(0) = &
        0.5_RealK * wgt(0) * width_lower * fraction_lower**2
      width_upper=nu(i_begin+2)-nu(i_begin+1)
      weight_coeff(1) = 0.5_RealK * &
        (wgt(1)*(width_lower*fraction_lower+width_upper) + &
        wgt(0)*width_lower*fraction_lower*(1.0_RealK-fraction_lower))
!
!     Middle points:
      DO i=2, n_int_weight-2
        width_lower=width_upper
        width_upper=nu(i_begin+i+1)-nu(i_begin+i)
        weight_coeff(i) = 0.5_RealK * &
          (width_lower+width_upper) * wgt(i)
      ENDDO
!
!     Upper extreme points:
      width_lower=width_upper
      width_upper=nu(i_end)-nu(i_end-1)
      fraction_upper=(nu_high-nu(i_end-1)) / width_upper
      weight_coeff(n_int_weight-1) = 0.5_RealK * &
        (wgt(n_int_weight-1) * &
        (width_lower+width_upper*fraction_upper) + &
        wgt(n_int_weight)*width_upper*fraction_upper * &
        (1.0_RealK-fraction_upper))
      weight_coeff(n_int_weight) = 0.5_RealK * &
        wgt(n_int_weight)*fraction_upper**2*width_upper
!
    CASE(2)
!
!     With two intervals the central weight is exceptional
!
      width_lower=nu(i_begin+1)-nu(i_begin)
      fraction_lower=(nu(i_begin+1)-nu_low) / width_lower
      width_upper=nu(i_end)-nu(i_end-1)
      fraction_upper=(nu_high-nu(i_end-1)) / width_upper
!
!     Contributions over lower interval.
!
      weight_coeff(0) = &
        0.5_RealK*wgt(0)*width_lower*fraction_lower**2
      weight_coeff(1) = &
        0.5_RealK*width_lower*fraction_lower * &
        (wgt(0)*(1.0_RealK-fraction_lower)+wgt(1))
!
!     Contributions from upper interval.
!
      weight_coeff(1)=weight_coeff(1) + &
        0.5_RealK*width_upper*fraction_upper * &
        (wgt(1)+wgt(2)*(1.0_RealK-fraction_upper))
      weight_coeff(2)=0.5_RealK*width_upper*wgt(2)*fraction_upper**2
!
    CASE(1)
!
!     Only the edge points contribute
!
!     Calculate fractions of the interval delimited by the
!     edges of the band.
      width_upper=nu(i_end)-nu(i_begin)
      fraction_lower=(nu_low-nu(i_begin)) / width_upper
      fraction_upper=(nu_high-nu(i_begin)) / width_upper
      weight_coeff(0)=0.5_RealK*(nu_high-nu_low) * &
        (wgt(0)*(1.0_RealK-fraction_lower) + &
        wgt(1)*(1.0_RealK-fraction_upper))
      weight_coeff(1)=0.5_RealK*(nu_high-nu_low) * &
        (wgt(0)*fraction_lower+wgt(1)*fraction_upper)
!
  END SELECT
!
  DEALLOCATE(wgt)
!
!
! Evaluate integral of weighting function.
  integral_weight = SUM(weight_coeff(0:n_int_weight))
!
!
!
  RETURN
!
!
!
CONTAINS
!
!
!
  SUBROUTINE apply_response_int
    USE spline_evaluate_mod, ONLY: spline_evaluate
    IMPLICIT NONE
!
!
!
    REAL  (RealK) :: response_0
!     Value of the response function at the point of evaluation
!
!
!
    DO i=0, n_int_weight
      CALL spline_evaluate(ierr, filter%n_pts, &
        filter%wavenumber, filter%response, filter%d2_response, &
        nu_wgt(i), response_0)
      IF (ierr == i_err_range) THEN
!       Recover harmlessly form errors of interpolation.
        response_0 = 0.0_RealK
        ierr = i_normal
      ENDIF
      IF (ierr /= i_normal) RETURN
      wgt(i) = wgt(i) * response_0
    ENDDO
!
!
!
  END SUBROUTINE apply_response_int

!
!
!
END SUBROUTINE weightings_single_90
