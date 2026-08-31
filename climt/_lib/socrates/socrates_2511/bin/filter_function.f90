! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to caclulate weights for a filter function.
!
! Method:
!   Appropriate weights for each band in the spectral file
!   are calculted from a filter function which is read in.
!   Splines are used to represent the function, so many
!   points should be used for functions which are not smooth.
!
!- ----------------------------------------------------------------------------
SUBROUTINE filter_function(ierr, n_band, wave_length_short, wave_length_long,  &
                           weight_band, nd_band)

  USE realtype_rd, ONLY: RealK
  USE error_pcf, ONLY: i_normal, i_err_fatal, i_err_range
  USE def_inst_flt, ONLY: StrFiltResp
  USE spline_fit_mod, ONLY: spline_fit

  IMPLICIT NONE


! Dummy arguments.
  INTEGER, Intent(INOUT) :: ierr
!       Error flag
  INTEGER, Intent(IN) :: n_band
!       Number of bands
  INTEGER, Intent(IN) :: nd_band
!       Size allocated for spectral bands
  REAL (RealK), Intent(IN) :: wave_length_short(nd_band)
!       Shortwave limits of each band
  REAL (RealK), Intent(IN) :: wave_length_long(nd_band)
!       Longwave limits of each band
  REAL (RealK), Intent(OUT) :: weight_band(nd_band)
!       Weights for each band

! Local variables
  INTEGER :: i
!       Loop variable
  REAL (RealK) :: wavenumber_low
!       Low wavenumber of band
  REAL (RealK) :: wavenumber_high
!       High wavenumber of band
  TYPE (StrFiltResp) :: filter
!       Instrumental response function


! Read the filter function
  CALL read_instrument_response_90(filter, ierr)

! Establish a spline fit to the filter function.
  CALL spline_fit(filter%n_pts, filter%wavenumber, filter%response,            &
    filter%d2_response)

! Calculate the weighting function for each band by integrating
! the spline fit
  DO i=1, n_band
    wavenumber_low=1.0_RealK/wave_length_long(i)
    wavenumber_high=1.0_RealK/wave_length_short(i)
    CALL integrate_spline(ierr, wavenumber_low, wavenumber_high,               &
      filter%n_pts, filter%wavenumber, filter%response,                        &
      filter%d2_response, weight_band(i))
    IF (ierr /= i_normal) THEN
      IF (ierr == i_err_range) THEN
!       If the band lies outside the spline range the weighting is 0.
!       We recover from this error.
        weight_band(i)=0.0_RealK
        ierr=i_normal
      ELSE
        ierr=i_err_fatal
      END IF
    END IF
    weight_band(i)=weight_band(i)/(wavenumber_high-wavenumber_low)
  END DO

END SUBROUTINE filter_function
