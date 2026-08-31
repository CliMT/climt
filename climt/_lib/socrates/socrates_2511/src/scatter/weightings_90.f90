! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate weightings across a spectral band.
!
SUBROUTINE weightings_90(ierr, &
  i_weight, &
  SolarSpec, t, &
  include_instrument_response, filter, &
  nu_low, nu_high, &
  n_exclude, nu_low_exclude, nu_high_exclude, &
  n_nu, nu, &
  i_begin, n_int_weight, weight_coeff)
!
! Description:
!   The array nu holds the points at which weightings
!   are required. This routine deals with excluded regions
!   in the band. The actual weights are calculated by 
!   WEIGHTINGS_SINGLE for a contiguous band.
!
!
!
! Modules used:
  USE realtype_rd
  USE def_std_io_icf
  USE def_inst_flt
  USE def_solarspec
  USE dimensions_spec_ucf
  USE dimensions_pp_ucf
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
!
  INTEGER, Intent(IN) :: i_weight
!   Method of weightng in frequency
  TYPE (StrSolarSpec), Intent(IN) :: SolarSpec
!   Solar spectrum
  INTEGER, Intent(IN) :: n_nu
!   Number of frequencies
  INTEGER, Intent(IN) :: n_exclude
!   Number of regions to exclude
!
  INTEGER, Intent(OUT) :: i_begin
!   Beginning of data range
  INTEGER, Intent(OUT) :: n_int_weight
!   Number of weighting intervals
!
  LOGICAL, Intent(IN) :: include_instrument_response
!   Flag to include the instrumental response function
  TYPE  (StrFiltResp), Intent(IN) :: filter
!   Instrumental response function
!
  REAL  (RealK), Intent(IN) :: t
!   Temperature of black body
  REAL  (RealK), Intent(IN) :: nu_low
!   Lowest frequency for weighting
  REAL  (RealK), Intent(IN) :: nu_high
!   Highest frequency for weighting
  REAL  (RealK), Intent(IN), Dimension(n_nu) :: nu
!   Abscissae of input data
  REAL  (RealK), Intent(IN), Dimension(n_exclude) :: nu_low_exclude
!   Shorter limits of excluded regions.
  REAL  (RealK), Intent(IN), Dimension(n_exclude) :: nu_high_exclude
!   Longer limits of excluded regions.
!
  REAL  (RealK), Intent(OUT), Dimension(0: n_nu) :: weight_coeff
!   Coefficients at weighting frequencies
!
!
! Local variables
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: i_begin_exclude
!   Initial point of excluded region
  INTEGER :: n_int_weight_exclude
!   Number of intervals in excluded region
  INTEGER :: offset
!   Offset of beginning of excluded region
!
  REAL  (RealK) :: integral_weight
!   Integral of weighting function
  REAL  (RealK) :: integral_weight_exclude
!   Integral of weighting function across excluded region.
  REAL  (RealK) :: weight_coeff_exclude(0: n_nu)
!   Weighting coefficeints for excluded regions.
!
! Subroutines called:
  EXTERNAL &
    weightings_single_90
!
!
!
! Check that the range covered by the data is greater than the
! band to be treated.
  IF ( (nu(n_nu) < nu_high) .OR. (nu(1) > nu_low) ) THEN
    WRITE(iu_err, '(/a)') &
      '*** Error: data do not cover the range to be weighted.'
    ierr=i_err_fatal
    RETURN
  ENDIF
!
  CALL weightings_single_90(ierr, &
    i_weight, &
    SolarSpec, t, &
    include_instrument_response, filter, &
    nu_low, nu_high, &
    n_nu, nu, &
    i_begin, n_int_weight, weight_coeff, &
    integral_weight)
  IF (ierr /= i_normal) return
!
! Remove contributions from the excluded regions.
  DO i=1, n_exclude
    CALL weightings_single_90(ierr, &
      i_weight, &
      SolarSpec, t, &
      include_instrument_response, filter, &
      nu_low_exclude(i), nu_high_exclude(i), &
      n_nu, nu, &
      i_begin_exclude, n_int_weight_exclude, &
      weight_coeff_exclude, &
      integral_weight_exclude)
    IF (ierr /= i_normal) return
!
!   Adjust the weighting coefficients and the integral.
    integral_weight=integral_weight-integral_weight_exclude
    offset=i_begin_exclude-i_begin
    DO j=0, n_int_weight_exclude
      weight_coeff(offset+j) = weight_coeff(offset+j) - &
        weight_coeff_exclude(j)
    ENDDO
  ENDDO
!
! Normalize by the integral of the weighting coefficients.
  weight_coeff(0:n_int_weight) = weight_coeff(0:n_int_weight) / &
    integral_weight
!
!
!
  RETURN
END SUBROUTINE weightings_90
