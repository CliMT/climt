! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 2.
!
! Method:
!   A solar spectrum, giving the irradiance against wavelength,
!   is read in. The total flux in the spectrum is determined,
!   using a Rayleigh-Jeans distribution for the far infra-red.
!   The fraction of this flux in each band is calculated. The
!   band will probably not cover the whole spectrum. If 
!   simulating observations this is as required, since the
!   spectral cut-offs of the observing instruments must be 
!   matched. In atmospheric models the whole flux should be
!   calculated and there is therefore an option to enhance
!   the fluxes in the outside bands to include the full flux.
!
!- ---------------------------------------------------------------------
SUBROUTINE make_block_2(Spectrum, SolarSpec, ierr)

  USE realtype_rd
  USE def_spectrum
  USE def_solarspec
  USE rad_pcf
  USE dimensions_pp_ucf
  USE def_inst_flt, ONLY: StrFiltResp

  IMPLICIT NONE

  TYPE (StrSpecData), INTENT(INOUT), TARGET :: Spectrum
!   Spectral file to be assigned
  TYPE (StrSolarSpec), INTENT(INOUT) :: SolarSpec
!   Solar spectrum
  INTEGER, INTENT(INOUT) :: ierr
!   Error flag

! Local Variables
  INTEGER :: ios
!   IO status
  CHARACTER (LEN=1) :: char_yn
!   Character response variable
  LOGICAL :: l_enhance
!   Enhance outer bands
  LOGICAL :: l_filter
!   Flag for filter function
  TYPE (StrFiltResp) :: filter
!   Instrumental response function

  IF (ALLOCATED(Spectrum%Solar%solar_flux_band)) &
      DEALLOCATE(Spectrum%Solar%solar_flux_band)
  ALLOCATE(Spectrum%Solar%solar_flux_band(Spectrum%Dim%nd_band))

  WRITE(*, '(/A)') 'Is a filter function required (Y/N)?'
  DO
    READ(*, '(a)') char_yn
    IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
      CALL read_instrument_response_90(filter, ierr)
      l_filter=.TRUE.
      EXIT
    ELSE IF ( (char_yn == 'N').OR.(char_yn == 'n') ) THEN
      l_filter=.FALSE.
      EXIT
    ELSE
      WRITE(*, '(a)') '+++ Unrecognised response: '
      WRITE(*, '(a)') 'Please re-type.'
    END IF
  END DO

! Obtain the solar spectrum if data are not already present.
  IF (SolarSpec%n_points > 0) THEN
    WRITE(*, '(/a/)')'Previous solar spectrum will be used.'
  ELSE
    CALL read_solar_spectrum(SolarSpec, ierr)
    IF (ierr /= i_normal) RETURN
  END IF

! The radiance outside the nominal limits of the spectrum can be
! assigned to the edging bands.
  WRITE(*, '(/a)') 'Assign solar flux outside given bands ' &
    //'to outside bands? (y/n)'
  DO
    READ(*, '(a)') char_yn
    IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
      l_enhance=.TRUE.
      EXIT
    ELSE IF ( (char_yn == 'N').OR.(char_yn == 'n') ) THEN
      l_enhance=.FALSE.
      EXIT
    ELSE
      WRITE(*, '(a)') '+++ Unrecognised response: '
      WRITE(*, '(a)') 'Please re-type.'
    END IF
  END DO

  CALL make_block_2_1(Spectrum, SolarSpec, filter, &
    l_filter, l_enhance, .TRUE., ierr)

END SUBROUTINE make_block_2
