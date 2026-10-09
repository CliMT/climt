! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

!+ Subroutine to read the schema describing a spectral file.
!
SUBROUTINE read_schema_spectrum(ierr
    , n_band, n_absorb, type_absorb
    , n_aerosol, type_aerosol
    , wave_length_short, wave_length_long
    , n_band_absorb, index_absorb
    , n_band_continuum, index_continuum
    , l_exclude, n_band_exclude, index_exclude
    )
!
! Method:
!      The schema is read to find directives which are interpreted.
!      A separate subroutine is called to process the directive
!      *BAND.
!
!
!
! Modules to set types of variables:
  USE realtype_rd
  USE dimensions_spec_ucf
  USE dimensions_io_ucf
  USE error_pcf
  USE gas_list_pcf
  USE def_std_io_icf
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments.
!
  INTEGER :: , Intent(INOUT) :: ierr
!   Error flag
!
! Elements of spectrum
  INTEGER :: , Intent(OUT) ::
      n_band
!       Number of spectral bands
! List of excluded bands
  logical, Intent(OUT) ::
      l_exclude
!       Flag to exclude bands
  INTEGER :: , Intent(OUT) ::
      n_band_exclude(npd_band)
!       Number bands excluded from this band
    , index_exclude(npd_exclude, npd_band)
!       List of excluded bands
! Gaseous absorption
  INTEGER :: , Intent(OUT) ::
      n_absorb
!       Number of absorbers
    , type_absorb(npd_species)
!       Type numbers of absorbing species
    , n_band_absorb(npd_band)
!       Number of absorbers in this band
    , index_absorb(npd_species, npd_band)
!       List of absorbers active in specific bands
! Continuum processes
  INTEGER :: , Intent(OUT) ::
      n_band_continuum(npd_band)
!       Number of continua in this band
    , index_continuum(npd_band, npd_continuum)
!       List of continua active in specific bands
! Aerosol processes
  INTEGER :: , Intent(OUT) ::
      n_aerosol
!       Number of aerosols
    , type_aerosol(npd_aerosol_species)
!       Types of aerosols
  real  (RealK), Intent(OUT) ::
      wave_length_short(npd_band)
!       Short end of bands
    , wave_length_long(npd_band)
!       Long end of bands
!
!
! Local variables.
  CHARACTER
      line*80
!       Line read from file
    , ch_unit*6
!       Character string identifying unit for limits of bands
  logical
      l_set_unit
!       Flag set to true when units have been defined
  INTEGER :: 
      ios
!       Flag recording status of i/o
    , index_type(npd_gases)
!       Indexing numbers of types of gases
    , i
!       Loop variable
    , j
!       Loop variable
!
! Subroutines called:
  EXTERNAL
      open_file_in, set_band, read_list
!
!
  data index_type/npd_gases*0/
!
!
!
! Open the file containing the schema for the spectrum.
  CALL get_free_unit(ierr, iu_schema)
  IF (ierr /= i_normal) RETURN
  CALL open_file_in(ierr, iu_schema, &
    'Give the name of the file containing the schema for the spectrum.')
!
! Initialize the properties of the spectrum.
  n_band=0
  n_absorb=0
  n_aerosol=0
  l_exclude= .FALSE. 
!
! Read each line looking for directives. Before any bands are set
! we need one setting the unit for wavelength.
  l_set_unit= .FALSE. 
  ios=0
  DO while (ios == 0)
    READ(iu_file_in, '(a)', iostat=ios) line
    IF (ios /= 0) goto 1
!
!   Check for directives setting wavelengths.
    IF (line(1:6) == '*METRE') THEN
      l_set_unit= .TRUE. 
      ch_unit='metre '
    ELSE IF (line(1:7) == '*MICRON') THEN
      l_set_unit= .TRUE. 
      ch_unit='micron'
    ELSE IF (line(1:7) == '*INV_CM') THEN
      l_set_unit= .TRUE. 
      ch_unit='inv_cm'
    ELSE IF (line(1:5) == '*BAND') THEN
      IF (l_set_unit) THEN
        n_band=n_band+1
        CALL set_band(ierr
          , line(7:80), ch_unit
          , n_band, n_absorb
          , n_band_absorb, index_absorb
          , n_band_continuum, index_continuum
          , wave_length_short, wave_length_long
          , type_absorb, index_type
          , l_exclude, n_band_exclude, index_exclude
          )
        IF (ierr /= i_normal) return
      ELSE
        WRITE(iu_err, '(/a)') 'A unit for wavelength must be '
          //'set before bands may be declared.'
        ierr=i_err_fatal
        RETURN
      ENDIF
    ELSE IF (line(1:8) == '*AEROSOL') THEN
      CALL read_list(ierr, line(10:80), 70
        , n_aerosol, type_aerosol)
      IF (ierr == i_err_io) THEN
        WRITE(iu_err, '(/a)') '*** Error: incorrect aerosol data'
        ierr=i_err_fatal
        RETURN
      ENDIF
      IF (ierr /= i_normal) return
    ENDIF
  ENDDO
!
! Check that excluded bands are in range.
  DO i=1, n_band
    DO j=1, n_band_exclude(i)
      IF ( (index_exclude(j, i) < 1) .OR. 
           (index_exclude(j, i) > n_band) ) THEN
        WRITE(iu_err, '(/a, i5, a)') 'An excluded band in band '
          , i, ' is out of range.'
        ierr=i_err_fatal
        RETURN
      ENDIF
    ENDDO
  ENDDO
!
!
!
  RETURN
END SUBROUTINE read_schema_spectrum
