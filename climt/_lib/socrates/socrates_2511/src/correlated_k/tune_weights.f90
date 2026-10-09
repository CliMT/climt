! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to tune weights for a particular atmosphere.
!
PROGRAM tune_weights
!
! Description:
!    This code tunes the weights of the terms in a c-k fit to match
! the heating rates calculated using a more highly resolved spectrum.
!
! Method:
!   A file of atmospheres is read in and heating rates are calculated
! for each k-term: the weights for the k-terms are then determined using
! a least-squares fitting procedure.
!
! NOTE: This program antedates the convert of the input system from CDL
! to F90. File suffixes are hard-wired. A final version of this
! program requires tidying up. For similar reasons solar tuning is done
! with an overhead sun. 
!
! Modules used:
  USE realtype_rd
  USE def_spectrum
  USE def_atm
  USE def_std_io_icf
  USE gas_list_pcf
  USE rad_pcf
  USE rad_ccf, ONLY: grav_acc, cp_air_dry, seconds_per_day
!
!
  IMPLICIT NONE
!
!
  TYPE (StrSpecData) :: RefSpectrum
!   Reference spectral file
  CHARACTER (LEN=80) :: file_spectral_tune
!   Name of the spectral file to be tuned
  TYPE (StrSpecData), Target :: TuneSpectrum
!   Reference spectral file
  LOGICAL :: l_interactive
!   Flag for interactive operstion
  INTEGER :: ierr = i_normal
!   Error flag
  INTEGER :: ios
!   I/O error flag
!
  INTEGER :: isolir
!   Spectral region
  INTEGER :: i_gas
!   Idenitifier for gas to be tuned
  LOGICAL :: l_continuum
!   Flag to include the self-broadened continuum in the tuning process
  INTEGER :: tune_band
!   Band of low-resolution spectral file to be tuned
  INTEGER :: first_band
!   First band in high-resolution spectral file marking the range
!   of bands correponding to the band to be tuned
  INTEGER :: last_band
!   Last band in high-resolution spectral file marking the range
!   of bands correponding to the band to be tuned
  REAL  (RealK), Dimension(:), Pointer :: weight_band_ref
!   Weighting of the reference bands
!
  INTEGER, Pointer :: n_k
!   Number of terms in the band being looked at
  INTEGER :: i_index_ref
!   Index of gas to be tuned in the reference file
  INTEGER :: i_index_tune
!   Index of gas to be tuned in the tuning file
!
  REAL :: offset
!   Offset to heating rates
  REAL :: gamma
!   Weighting applied to surface flux differences
!
  TYPE  (StrAtm) :: Atm
!   Atmospheric profile
  REAL  (RealK), Allocatable :: broad(:, :)
!   Continuum broadening density
  INTEGER :: first_layer
!   First layer in the range to be used for tuning
  INTEGER :: last_layer
!   Last layer in the range to be used for tuning
!
  REAL  (RealK), Allocatable :: flux_up(:, :)
!   Upward fluxes
  REAL  (RealK), Allocatable :: flux_down(:, :)
!   Downward fluxes
!
  REAL  (RealK), Allocatable :: hr_ref(:, :)
!   Heating rates calculated for the reference spectral file
  REAL  (RealK), Allocatable :: hr_tune(:, :, :)
!   Heating rates calculated for the reference spectral file
  REAL  (RealK), Allocatable :: flux_surf_ref(:)
!   The net surface flux from the reference spectrum 
  REAL  (RealK), Allocatable :: flux_surf_tune(:, :)
!   The net surface fluxes from each k-term of the spectrum 
!   to be tuned
!
  REAL  (RealK), Allocatable :: a(:, :)
!   Array of coefficients for SVD
  REAL  (RealK), Allocatable :: w(:)
!   Singular values of SVD decoposition
  REAL  (RealK), Allocatable :: v(:, :)
!   Upper triangular matrix of the SVD decomposition
  REAL  (RealK), Allocatable :: b(:)
!   RHSs for SVD equations
  REAL  (RealK), Allocatable :: x(:)
!   Solution of SVD equations
  REAL  (RealK), Allocatable :: wrk(:)
!   Workspace for SVD decomposition
  REAL  (RealK) :: threshold
!   Tolerance for zeroing singular values in SVD
!
! External functions:
  LOGICAL :: set_interactive
!   Function to set the flag for interactive operation
  EXTERNAL &
    set_interactive
!
!- End of Header
!
!
!
! Set the flag for interactive operation
  l_interactive=set_interactive()
!
! Obtain the spectral files
  CALL get_spectra
  IF (ierr == i_err_fatal) STOP
!
!
! Enter the reference details for the calculation
  CALL get_details
  IF (ierr == i_err_fatal) STOP
!
! Define the set of atmospheres in which tuning to be done.
  CALL get_atm
  IF (ierr == i_err_fatal) STOP
!
  WRITE(iu_stdout, "(a)") "Enter the spectral region."
  READ(iu_stdin, *) isolir
!
! Allocate the flux arrays.
  ALLOCATE(flux_up(Atm%n_profile, 0:Atm%n_layer))
  ALLOCATE(flux_down(Atm%n_profile, 0:Atm%n_layer))
  ALLOCATE(hr_ref(Atm%n_profile, Atm%n_layer))
  ALLOCATE(hr_tune(Atm%n_profile, Atm%n_layer, n_k))
  ALLOCATE(flux_surf_ref(Atm%n_profile))
  ALLOCATE(flux_surf_tune(Atm%n_profile, n_k))
!
! Calculate the reference heating rates
  CALL ref_hrate
!
! Calculate the heating rates for tuning
  CALL tuning_hrate
!
! Define the relative weighting of surface fluxes relative 
! to heating rates. Gamma appears to overwhelm the contributions from the
! heating rates and is set to 0 for now.
  gamma = 0.0 * LOG(ABS(Atm%p(1, last_layer)/Atm%p(1, first_layer)))
!
! Form the matrix for least squares calculation.
  CALL form_svd_matrix
!
! Solve by SVD decomposition.
  ALLOCATE(v(n_k+1, n_k+1))
  ALLOCATE(w(n_k+1))
  ALLOCATE(wrk(n_k+1))
  ALLOCATE(x(n_k+1))
  CALL svd_decompose(ierr, &
    a, n_k+1, n_k+1, &
    n_k+1, n_k+1, w, v, wrk)
  threshold = 1.0E+03_RealK * EPSILON(threshold) * MAXVAL(w)
  WHERE (w < threshold)
    w = 0.0_RealK
  ENDWHERE
  CALL back_substitute( &
    a, w, v, n_k+1, n_k+1, &
    n_k+1, n_k+1, b, &
    x, wrk)
    TuneSpectrum%Gas%w(1:n_k, tune_band, i_index_tune) = x(1:n_k)
!
! Write out the adjusted spectrum
  CALL out_spectrum(file_spectral_tune, TuneSpectrum, ierr)
!
!
!
  STOP
!
!
!
CONTAINS
!
!
!
  SUBROUTINE get_spectra

!   Local variables
    CHARACTER (LEN=80) :: file_spectral_ref
!     Name of the spectral file

!   Read in the reference spectral file.
    WRITE(*, "(a)") "Enter the name of the reference spectral file."
    READ(*, "(a)") file_spectral_ref
    CALL read_spectrum(file_spectral_ref, RefSpectrum)

!   Read in the spectral file to be tuned.
    WRITE(*, "(a)") &
      "Enter the name of the spectral file to be tuned."
    READ(*, "(a)") file_spectral_tune
    CALL read_spectrum(file_spectral_tune, TuneSpectrum)

  END SUBROUTINE get_spectra
!
!
!
  SUBROUTINE get_details
!
!
!
!   Local variables:
    LOGICAL :: l_found
!     Flag used in locating a gas
    INTEGER :: i_swap
!     Temporary swapping variable
    INTEGER :: i
!     Loop variable
    INTEGER :: j
!     Loop variable
!
    REAL  (RealK) :: f_min
!     Minimum frequency (measured in wavenumbers) covered by both
!     the current reference band and the band to be tuned
    REAL  (RealK) :: f_max
!     Maximum frequency (measured in wavenumbers) covered by both
!     the current reference band and the band to be tuned
    REAL  (RealK) :: f_min_e
!     Minimum frequency (measured in wavenumbers) of the excluded
    REAL  (RealK) :: f_max_e
!     Maximum frequency (measured in wavenumbers) of the excluded
!     band
!
!
!
!   Choose the gas and range of bands  
    WRITE(*, "(a)") "Enter the identifier for the gas to be considered."
    DO
      READ(*, *, IOSTAT=ios) i_gas
      IF (ios /= 0) THEN
        WRITE(iu_err, "(/a)") "Erroneous response"
        IF (l_interactive) THEN
          WRITE(*, "(a)") "Please re-enter."
        ELSE
          ierr=i_err_fatal
          STOP
        ENDIF
      ELSE
!       Check that the gas is in both spectral files.
        CALL find_index(RefSpectrum%Gas%n_absorb, &
                        RefSpectrum%Gas%type_absorb, &
                        i_index_ref)
        IF (ierr /= i_normal) THEN
          WRITE(iu_err, "(a)") &
            "*** Error: This gas is not in the reference spectrum."
          STOP
        ENDIF
        CALL find_index(TuneSpectrum%Gas%n_absorb, &
                        TuneSpectrum%Gas%type_absorb, &
                        i_index_tune)
        IF (ierr /= i_normal) THEN
          WRITE(iu_err, "(a)") &
            "*** Error: This gas is not in the spectrum to be tuned."
          RETURN
        ENDIF
        EXIT
      ENDIF
    ENDDO
!
!   Select the band to be tuned.
    WRITE(iu_stdout, '(/a)') 'Enter the number of the band to be tuned.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) tune_band
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ Erroneous response.'
        IF(l_interactive) THEN
          WRITE(iu_stdout, '(a)') 'Please re-enter.'
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ELSE
        l_found =  .FALSE.
        DO j=1, TuneSpectrum%Gas%n_band_absorb(tune_band)
          IF (i_index_tune == &
              TuneSpectrum%Gas%index_absorb(j, tune_band)) THEN
            l_found = .TRUE.
          ENDIF
        ENDDO
        IF (.NOT.l_found) THEN
          WRITE(iu_err, "(/A)") &
            "*** Error: This band does not contain the requested gas."
          RETURN
         ENDIF
        EXIT
      ENDIF
    ENDDO
!
!   For convenience of abbreviation, define a pointer 
!   to the number of terms.
    n_k => TuneSpectrum%Gas%i_band_k(tune_band, i_index_tune)  
!
!   Select the range of bands to be covered in the reference file.
    WRITE(iu_stdout, '(/a)') &
      'Enter the corresponding range of bands in the reference file.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) first_band, last_band
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ Erroneous response.'
        IF(l_interactive) THEN
          WRITE(iu_stdout, '(a)') 'Please re-enter.'
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ELSE
        IF (first_band > last_band) THEN
          i_swap=first_band
          first_band=last_band
          last_band=i_swap
        ENDIF
!       Check limits.
        IF ( (first_band < 1) .OR. &
             (last_band > RefSpectrum%Basic%n_band) ) THEN
          WRITE(iu_err, '(a)') '+++ Response out of range.'
          IF (l_interactive) THEN
            WRITE(iu_stdout, '(a)') 'Please re-enter.'
          ELSE
            ierr=i_err_fatal
            RETURN
          ENDIF
        ENDIF
        EXIT
      ENDIF
    ENDDO
!
!   Because the reference bands may not exactly match the bands in
!   the file to be tuned we apply a weighting based on the coverage
!   in frequency space.
    ALLOCATE(weight_band_ref(first_band:last_band))
    DO i = first_band, last_band
      f_max = 1.0_RealK / &
        MAX( RefSpectrum%Basic%wavelength_short(i), &
             TuneSpectrum%Basic%wavelength_short(tune_band) )                   
      f_min = 1.0_RealK / &
        MIN( RefSpectrum%Basic%wavelength_long(i), &
             TuneSpectrum%Basic%wavelength_long(tune_band) )                   
!
!     Check for excluded regions. We assume that the excluded region
!     is wider than any band in the reference file.
      DO j = 1, TuneSpectrum%Basic%n_band_exclude(tune_band)
!
        f_max_e = 1.0_RealK / &
          TuneSpectrum%Basic%wavelength_short( &
            TuneSpectrum%Basic%index_exclude(j, tune_band) )
        f_min_e = 1.0_RealK / &
          TuneSpectrum%Basic%wavelength_long( &
            TuneSpectrum%Basic%index_exclude(j, tune_band) )
!
        IF ( (f_min < f_min_e) .AND. (f_max > f_min_e) ) THEN
          f_max = f_min_e
        ELSE IF ( (f_min >= f_min_e) .AND. (f_max < f_max_e) ) THEN
          f_max = f_min
        ELSE IF ( (f_min < f_max_e) .AND. (f_max >= f_max_e) ) THEN
          f_min = f_max_e
        ENDIF
      ENDDO
!
      weight_band_ref(i) = ( f_max - f_min ) / &
        ( 1.0_RealK / RefSpectrum%Basic%wavelength_short(i) - &
          1.0_RealK / RefSpectrum%Basic%wavelength_long(i) )
!
    ENDDO
!
    WRITE(iu_stdout, "(a)") &
      "Enter offset to heating rates (K/day)."
    DO
      READ(iu_stdin, *, IOSTAT=ios) offset
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ Erroneous response.'
        IF(l_interactive) THEN
          WRITE(iu_stdout, '(a)') 'Please re-enter.'
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO
!   Convert to Wm-2/kg:
    offset = offset * cp_air_dry / seconds_per_day
    
!
!
!
  END SUBROUTINE get_details
!
!
!
  SUBROUTINE find_index(n_absorb, type_absorb, i_index)
!
!
!
    INTEGER, Intent(IN) :: n_absorb
    INTEGER, Intent(IN), Dimension(n_absorb) :: type_absorb
    INTEGER, Intent(OUT) :: i_index
!
    INTEGER :: i
!
!
!
    i_index=0
    i=1
    DO
      IF (type_absorb(i) == i_gas) THEN
        i_index=i
        EXIT
      ELSE IF (i == n_absorb) THEN
        EXIT
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF (i_index == 0) THEN
      WRITE(iu_err, "(a)") "This gas is not in the spectral file."
      ierr = i_err_fatal
      RETURN
    ENDIF
!
!
!
  END SUBROUTINE find_index
!
!
!
  SUBROUTINE get_atm
!
!
  USE dimensions_field_ucf
  USE dimensions_cdl_ucf
!
! Local variables
  CHARACTER (LEN=80) :: base_name
!   Base name of CDL input files
  CHARACTER (LEN=80) :: cdl_file_name
!   Base name of CDL input files
  CHARACTER (LEN=12) :: suffix
!   Suffix appended to CDL files indicating the type
  CHARACTER (LEN=24) :: name_vert_coord
!   Name of the vertical coordinate
  LOGICAL :: l_vert_coord = .FALSE.
!   Flag indicating that the vertical coordinate has been set
  LOGICAL :: l_vert_coord_level = .FALSE.
!   Flag indicating that the vertical coordinate for fields at the
!   edges of layers has been set
  INTEGER :: n_lat = 0
!   Number of latitudes
  INTEGER :: n_lon = 0
!   Number of longitudes
  INTEGER :: n_level
!   Dummy number of levels
  INTEGER :: i
!   Loop variable
!
!
!
! Provisional code: once the CDL routines have been converted to F90,
! this can be removed, as allocation will be done in the CDL routines.
  ALLOCATE(Atm%lat(npd_profile))
  ALLOCATE(Atm%lon(npd_profile))
  ALLOCATE(Atm%p(npd_profile, npd_layer))
  ALLOCATE(Atm%t(npd_profile, npd_layer))
  ALLOCATE(Atm%p_level(npd_profile, 0:npd_layer))
  ALLOCATE(Atm%t_level(npd_profile, 0:npd_layer))
  ALLOCATE(Atm%gas_mix_ratio(npd_profile, npd_layer, 1))
  ALLOCATE(Atm%mass(npd_profile, npd_layer))
!
! Read in the atmosphere to be used for tuning.
  WRITE(iu_stdout, "(a)") "Enter the basename for the CDL files."
  READ(iu_stdin, "(a)") base_name
!
! Read the temperatures at the centres of layers
  cdl_file_name = TRIM(base_name)//'.t'
  CALL assign_input_vert_cdl(ierr, &
    cdl_file_name, 'temperatures', &
    l_vert_coord, name_vert_coord, &
    .TRUE., Atm%n_layer, .NOT.l_vert_coord, &
    n_lat, Atm%lat, n_lon, Atm%lon, 1, &
    Atm%n_profile, Atm%n_layer, &
    Atm%p, Atm%t, &
    npd_profile, npd_latitude, npd_longitude, 1, npd_layer, &
    npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var)
  IF (ierr /= i_normal) STOP
!
! Read the temperatures at the edges of layers.
  cdl_file_name = TRIM(base_name)//'.tl'
  CALL assign_input_vert_cdl(ierr, &
    cdl_file_name, 'edge temperatures', &
    l_vert_coord_level, name_vert_coord, &
    .TRUE., Atm%n_layer+1, .TRUE., &
    n_lat, Atm%lat, n_lon, Atm%lon, 0, &
    Atm%n_profile, n_level, &
    Atm%p_level, Atm%t_level, &
    npd_profile, npd_latitude, npd_longitude, 0, npd_layer, &
    npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var)
  IF (ierr /= i_normal) STOP
!
! Calculate the mass using hydrostatic balance.
  DO i=1, Atm%n_layer
    Atm%mass(1:Atm%n_profile, i) = &
      (Atm%p_level(1:Atm%n_profile, i) - &
       Atm%p_level(1:Atm%n_profile, i-1)) / grav_acc
  ENDDO
!
! Get the gaseous mixing ratios: temporarily, we use hardwirded
! suffixes.
  SELECT CASE(i_gas)
    CASE(IP_H2O)
      suffix='.q'
    CASE(IP_CO2)
      suffix='.co2'
    CASE(IP_O3)
      suffix='.o3'
    CASE(IP_CH4)
      suffix='.ch4'
    CASE(IP_N2O)
      suffix='.n2o'
    CASE(IP_O2)
      suffix='.o2'
  END SELECT
!
  cdl_file_name = TRIM(base_name)//TRIM(suffix)
  CALL assign_input_vert_cdl(ierr, &
    cdl_file_name, 'gaseous mixing ratios', &
    l_vert_coord, name_vert_coord, &
    .TRUE., Atm%n_layer, .FALSE., &
    n_lat, Atm%lat, n_lon, Atm%lon, 1, &
    Atm%n_profile, Atm%n_layer, &
    Atm%p, Atm%gas_mix_ratio, &
    npd_profile, npd_latitude, npd_longitude, 1, npd_layer, &
    npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var)
  IF (ierr /= i_normal) STOP
!
  WRITE(iu_stdout, '(a)') &
    'Enter range of layers for tuning.'
  DO
    READ(iu_stdin, *, IOSTAT=ios) first_layer, last_layer
    IF (ios /= 0) THEN
      WRITE(iu_err, '(a)') '+++ Erroneous response.'
      IF(l_interactive) THEN
        WRITE(iu_stdout, '(a)') 'Please re-enter.'
      ELSE
        ierr=i_err_fatal
        RETURN
      ENDIF
    ELSE
      EXIT
    ENDIF
  ENDDO
!
! Keep data in range.
  first_layer = MAX(first_layer, 1)
  last_layer  = MIN(last_layer, Atm%n_layer)
!
  WRITE(iu_stdout, '(a)') &
    'Should the continuum be used in tuning? (T/F)'
  READ(iu_stdin, *) l_continuum
  IF (l_continuum) THEN
    ALLOCATE(broad(npd_profile, npd_layer))
    broad = Atm%gas_mix_ratio(:, :, 1) * (28.966/18.0153 - 1.0)
    broad = 1.0 + broad
    broad = Atm%gas_mix_ratio(:, :, 1) * Atm%p * (28.966/18.0153) / (&
          28.966e-03 * 287.04 * Atm%t * broad )
  ENDIF
!
!
!
  END SUBROUTINE get_atm
!
!
!
  SUBROUTINE ref_hrate
!
!
!
!   Local variables:
    INTEGER :: i_band
!     Loop variable
    INTEGER :: i_k
!     Loop variable
    INTEGER :: i
!     Loop variable
    REAL :: normalize
!     Variable to normalize the fluxes
!
!
    hr_ref(:, :) = 0.0_RealK
    flux_surf_ref(:) = 0.0_RealK
!
    DO i_band=first_band, last_band
      DO i_k=1, RefSpectrum%Gas%i_band_k(i_band, i_index_ref)
!
!       Calculate fluxes
        SELECT CASE(isolir)
          CASE(IP_solar)
            CALL monochromatic_solar &
            (RefSpectrum%Gas%k(i_k, i_band, i_index_ref), &
             RefSpectrum%Gas%i_scale_fnc(i_band, i_index_ref), &
             RefSpectrum%Gas%scale(:, i_k, i_band, i_index_ref), &
             RefSpectrum%Gas%p_ref(i_index_ref, i_band), &
             RefSpectrum%Gas%t_ref(i_index_ref, i_band))
!           The monochromatic calculation is based on unit incident flux
!           and renormalization for the solar flux is required.
            normalize = RefSpectrum%Solar%solar_flux_band(i_band) / &
                    TuneSpectrum%Solar%solar_flux_band(tune_band)
          CASE(IP_infra_red)
            CALL monochromatic_IR &
              (RefSpectrum%Planck%thermal_coeff(:, i_band), &
               RefSpectrum%Planck%t_ref_planck, &
               RefSpectrum%Gas%k(i_k, i_band, i_index_ref), &
               RefSpectrum%Gas%i_scale_fnc(i_band, i_index_ref), &
               RefSpectrum%Gas%scale(:, i_k, i_band, i_index_ref), &
               RefSpectrum%Gas%p_ref(i_index_ref, i_band), &
               RefSpectrum%Gas%t_ref(i_index_ref, i_band), &
               RefSpectrum%Cont%k_cont(i_band, i_index_ref), &
               RefSpectrum%Cont%i_scale_fnc_cont(i_band, i_index_ref), &
               RefSpectrum%Cont%scale_cont(:, i_band, i_index_ref), &
               RefSpectrum%Cont%p_ref_cont(i_index_ref, i_band), &
               RefSpectrum%Cont%t_ref_cont(i_index_ref, i_band) &
              )
!           No normalization is required here, because the Planckian fit
!           is used directly.
            normalize = 1.0_RealK
        END SELECT
!
!       Convert to heating rates and store the net surface flux.
        DO i=1, Atm%n_layer
          hr_ref(:, i) = hr_ref(:, i) + weight_band_ref(i_band) * &
                         normalize * &
                         RefSpectrum%Gas%w(i_k, i_band, i_index_ref) * &
                         (flux_up(:, i) - flux_up(:, i-1) + &
                          flux_down(:, i-1) - flux_down(:, i)) / &
                         Atm%mass(:, i)
        ENDDO
        flux_surf_ref(:) = flux_surf_ref(:) + weight_band_ref(i_band) * &
                           normalize * &
                           RefSpectrum%Gas%w(i_k, i_band, i_index_ref) * &
                            (flux_down(:, Atm%n_layer) - &
                             flux_up(:, Atm%n_layer))
!
      ENDDO
    ENDDO
!
!
!
  END SUBROUTINE ref_hrate
!
!
!
  SUBROUTINE tuning_hrate
!
!
!
!   Local variables:
    INTEGER :: i_band
!     Loop variable
    INTEGER :: i_k
!     Loop variable
    INTEGER :: i
!     Loop variable
!
!
!
    DO i_k=1, TuneSpectrum%Gas%i_band_k(tune_band, i_index_tune)
!
!     Calculate fluxes
      SELECT CASE(isolir)
        CASE(IP_solar)
          CALL monochromatic_solar &
            (TuneSpectrum%Gas%k(i_k, tune_band, i_index_tune), &
             TuneSpectrum%Gas%i_scale_fnc(tune_band, i_index_tune), &
             TuneSpectrum%Gas%scale(:, i_k, tune_band, i_index_tune), &
             TuneSpectrum%Gas%p_ref(i_index_tune, tune_band), &
             TuneSpectrum%Gas%t_ref(i_index_tune, tune_band))
        CASE(IP_infra_red)
          CALL monochromatic_IR &
            (TuneSpectrum%Planck%thermal_coeff(:, tune_band), &
             TuneSpectrum%Planck%t_ref_planck, &
             TuneSpectrum%Gas%k(i_k, tune_band, i_index_tune), &
             TuneSpectrum%Gas%i_scale_fnc(tune_band, i_index_tune), &
             TuneSpectrum%Gas%scale(:, i_k, tune_band, i_index_tune), &
             TuneSpectrum%Gas%p_ref(i_index_tune, tune_band), &
             TuneSpectrum%Gas%t_ref(i_index_tune, tune_band), &
             TuneSpectrum%Cont%k_cont(tune_band, i_index_tune), &
             TuneSpectrum%Cont%i_scale_fnc_cont(tune_band, i_index_tune), &
             TuneSpectrum%Cont%scale_cont(:, tune_band, i_index_tune), &
             TuneSpectrum%Cont%p_ref_cont(i_index_tune, tune_band), &
             TuneSpectrum%Cont%t_ref_cont(i_index_tune, tune_band) &
            )
      END SELECT
!
!     Convert to heating rates
      DO i=1, Atm%n_layer
        hr_tune(:, i, i_k) = (flux_up(:, i) - flux_up(:, i-1) + &
                             flux_down(:, i-1) - flux_down(:, i)) / &
                             Atm%mass(:, i)
      ENDDO
      flux_surf_tune(:, i_k) = flux_down(:, Atm%n_layer) - &
                               flux_up(:, Atm%n_layer)
!
    ENDDO
!
!
!
  END SUBROUTINE tuning_hrate
!
!
!
  SUBROUTINE monochromatic_solar(kabs, i_scale_fnc, &
    scale_cf, p_ref, t_ref)
!
!
!
    REAL  (RealK) :: kabs
!     Absorption coefficient
    REAL  (RealK), Dimension(:) :: scale_cf
!     Coefficients of the scaling function
    REAL  (RealK) :: p_ref
!     Reference pressure
    REAL  (RealK) :: t_ref
!     Reference temperature
    REAL  (RealK) :: mu_0
!     Cosine of the solar zenith angle
    INTEGER :: i_scale_fnc
!     Type of scaling function
    INTEGER :: i
!     Loop variable
!
!
!
!   Calculate fluxes assuming an overhead sun and a 
!   nonreflecting surface.
!
!   Upward fluxes are 0.
    flux_up(1:Atm%n_profile, 0:Atm%n_layer) = 0.0_RealK
!
!   Unit incident flux and a zenith angle of 60 degrees sun are assumed.
    flux_down(1:Atm%n_profile, 0) = 1.0_RealK
    mu_0 = 0.5_RealK
    DO i=1, Atm%n_layer
      flux_down(1:Atm%n_profile, i) = flux_down(1:Atm%n_profile, i-1) * &
        EXP( - kabs * scaling_fnc(i_scale_fnc, scale_cf, p_ref, t_ref, &
          Atm%p(1:Atm%n_profile, i), Atm%t(1:Atm%n_profile, i)) * &
          ( 1.0_RealK / mu_0 ) * &
          Atm%gas_mix_ratio(1:Atm%n_profile, i, 1) * &
          Atm%mass(1:Atm%n_profile, i) )
    ENDDO
!
!
!
  END SUBROUTINE monochromatic_solar
!
!
!
  SUBROUTINE monochromatic_IR(therm_coeff, t_ref_planck, &
    kline, i_scale_fnc, scale_cf, p_ref, t_ref, &
    kcont, i_scale_fnc_c, scale_cf_c, p_ref_c, t_ref_c)
!
!
!
    INTEGER :: i_scale_fnc
!     Type of scaling function
    REAL  (RealK), Intent(IN) :: kline
!     Line absorption coefficient
    REAL  (RealK), Dimension(:) :: scale_cf
!     Coefficients of the scaling function
    REAL  (RealK) :: p_ref
!     Reference pressure
    REAL  (RealK) :: t_ref
!     Reference temperature
!
    INTEGER :: i_scale_fnc_c
!     Type of continuum scaling function
    REAL  (RealK), Intent(IN) :: kcont
!     Continuum absorption coefficient
    REAL  (RealK), Dimension(:) :: scale_cf_c
!     Coefficients of the scaling function
    REAL  (RealK) :: p_ref_c
!     Reference pressure
    REAL  (RealK) :: t_ref_c
!     Reference temperature
!
    REAL  (RealK), Intent(IN) :: therm_coeff(:)
!     Coefficients fitting the Planckian
    REAL  (RealK), Intent(IN) :: t_ref_planck
!     Reference temperature for the Planckian
    REAL  (RealK), Parameter :: diff = 1.66_RealK
!     Diffusivity factor (fixed at Elsasser's value)
    REAL  (RealK), Dimension(Atm%n_profile, 0: Atm%n_layer) :: bf
!     Flux Planckian functions at edges
    REAL  (RealK), Dimension(Atm%n_profile) :: kabs
!     Overall absorption coefficient
    REAL  (RealK), Dimension(Atm%n_profile) :: d_tau
!     Optical depths of the layers
    REAL  (RealK), Dimension(Atm%n_profile) :: trans
!     Transmissions through the current layer
    INTEGER :: i
!     Loop variable
!
!
!
!   Calculate IR fluxes.
!
!   Calculate Planckian fluxes:
    bf = planckian(Atm%t_level(1:Atm%n_profile, 0:Atm%n_layer), &
                t_ref_planck, &
                therm_coeff)
!
!   Upward fluxes:
    flux_up(1:Atm%n_profile, Atm%n_layer) = bf(1:Atm%n_profile, Atm%n_layer)
    DO i=Atm%n_layer-1, 0, -1
      kabs = kline * &
        scaling_fnc(i_scale_fnc, scale_cf, p_ref, t_ref, &
          Atm%p(1:Atm%n_profile, i+1), Atm%t(1:Atm%n_profile, i+1))
      IF (l_continuum) THEN
        kabs = kabs + kcont * &
          scaling_fnc(i_scale_fnc_c, scale_cf_c, p_ref_c, t_ref_c, &
          Atm%p(1:Atm%n_profile, i+1), Atm%t(1:Atm%n_profile, i+1)) * &
          broad(1:Atm%n_profile, i+1)
      ENDIF
      d_tau(1:Atm%n_profile) = kabs * &
        Atm%gas_mix_ratio(1:Atm%n_profile, i+1, 1) * &
        Atm%mass(1:Atm%n_profile, i+1) 
      trans(1:Atm%n_profile) = EXP( -diff * d_tau(1:Atm%n_profile) )
      flux_up(1:Atm%n_profile, i) = flux_up(1:Atm%n_profile, i+1) * &
        trans(1:Atm%n_profile) + bf(1:Atm%n_profile, i) - &
        trans(1:Atm%n_profile) * bf(1:Atm%n_profile, i+1) + &
        ( (bf(1:Atm%n_profile, i+1) - bf(1:Atm%n_profile, i)) * &
          (1.0_RealK + SQRT(EPSILON(bf)) - trans(1:Atm%n_profile)) ) / &
          ( d_tau(1:Atm%n_profile) * diff + SQRT(EPSILON(bf))  )
    ENDDO
!
!   Downward fluxes:
    flux_down(1:Atm%n_profile, 0) = 0.0_RealK
    DO i=1, Atm%n_layer
      kabs = kline * &
        scaling_fnc(i_scale_fnc, scale_cf, p_ref, t_ref, &
          Atm%p(1:Atm%n_profile, i), Atm%t(1:Atm%n_profile, i))
      IF (l_continuum) THEN
        kabs = kabs + kcont * &
          scaling_fnc(i_scale_fnc_c, scale_cf_c, p_ref_c, t_ref_c, &
          Atm%p(1:Atm%n_profile, i), Atm%t(1:Atm%n_profile, i)) * &
          broad(1:Atm%n_profile, i)
      ENDIF
      d_tau(1:Atm%n_profile) = kabs * &
        Atm%gas_mix_ratio(1:Atm%n_profile, i, 1) * &
        Atm%mass(1:Atm%n_profile, i) 
      trans(1:Atm%n_profile) = EXP( -diff * d_tau(1:Atm%n_profile) )
      flux_down(1:Atm%n_profile, i) = flux_down(1:Atm%n_profile, i-1) * &
        trans(1:Atm%n_profile) + bf(1:Atm%n_profile, i) - &
        trans(1:Atm%n_profile) * bf(1:Atm%n_profile, i-1) - &
        ( (bf(1:Atm%n_profile, i) - bf(1:Atm%n_profile, i-1)) * &
          (1.0_RealK + SQRT(EPSILON(bf)) - trans(1:Atm%n_profile)) ) / &
          ( d_tau(1:Atm%n_profile) * diff + SQRT(EPSILON(bf))  )
    ENDDO
!
!
!
  END SUBROUTINE monochromatic_IR
!
!
!
  FUNCTION scaling_fnc(i_scale_fnc, cf, p_ref, t_ref, p, t) &
!
  RESULT (scaling)
!
!
!
    INTEGER :: i_scale_fnc
    REAL  (RealK) :: cf(:)
    REAL  (RealK) :: p_ref
!     Reference pressure
    REAL  (RealK) :: t_ref
!     Reference temperature
    REAL  (RealK) :: p(:)
    REAL  (RealK) :: t(:)
    REAL  (RealK), Pointer :: scaling(:)
    INTEGER :: k
!
!
    ALLOCATE(scaling(SIZE(p)))
    scaling(:) = 1.0_RealK
    SELECT CASE(i_scale_fnc)
      CASE (IP_scale_power_law)
        scaling(:) = EXP( cf(1)*LOG(p(:)/p_ref) + &
                          cf(2)*LOG(t(:)/t_ref) )
      CASE (IP_scale_power_quad)
        scaling(:) = EXP(cf(1)*LOG(p(:)/p_ref)) * &
          (1.0_RealK + cf(2)*(t(:)/t_ref-1.0_RealK) + &
            cf(3) * (t(:)/t_ref-1.0_RealK) * (t(:)/t_ref-1.0_RealK))
      CASE (IP_scale_doppler_quad)
        scaling(:) = EXP( cf(1)*LOG( (p(:)+cf(2)) / &
                                     (p_ref+cf(2)) ) ) * &
          (1.0_RealK + cf(3)*(t(:)/t_ref-1.0_RealK) + &
            cf(4) * (t(:)/t_ref-1.0_RealK) * (t(:)/t_ref-1.0_RealK))
  END SELECT
!
!
!
  END FUNCTION scaling_fnc
!
!
!
  FUNCTION planckian(t, t_ref_planck, cf) RESULT (bf)
!
    REAL  (RealK) :: t(:, :)
    REAL  (RealK) :: t_ref_planck
    REAL  (RealK), Pointer :: bf(:, :)
    REAL  (RealK) :: cf(:)
    INTEGER :: k
!
!
    ALLOCATE(bf(size(t, 1), LBOUND(t, 2):UBOUND(t, 2)))
    bf(:, :) = cf(SIZE(cf))
    DO k=SIZE(cf)-1, 1, -1
      bf = cf(k) + bf * ( t / t_ref_planck)
    ENDDO
!
!
!
  END FUNCTION planckian
!
!
!
  SUBROUTINE form_svd_matrix
!
!
!   Local variables
    INTEGER :: i_k
    INTEGER :: j_k
    INTEGER :: i
    INTEGER :: l
!
    ALLOCATE(a(n_k+1, n_k+1))
    ALLOCATE(b(n_k+1))
    a(:, :) = 0.0_RealK
    b(:)    = 0.0_RealK
    DO i_k=1, n_k
!
!     Begin with the surface terms
      DO l=1, Atm%n_profile
        b(i_k) = b(i_k) + &
                 gamma * flux_surf_tune(l, i_k) / flux_surf_ref(l) 
      ENDDO
      DO i=first_layer, last_layer
        DO l=1, Atm%n_profile
          b(i_k) = b(i_k) + hr_ref(l, i) * hr_tune(l, i, i_k) * &
            (1.0_RealK/(offset+ABS(hr_ref(l, i)))**2) * &
            LOG(Atm%p_level(l, i)/Atm%p_level(l, i-1))
        ENDDO
      ENDDO
!
      DO j_k=1, n_k 
!
!       Begin with the surface terms
        DO l=1, Atm%n_profile
          a(i_k, j_k) = a(i_k, j_k) + &
                                gamma * flux_surf_tune(l, i_k) * &
                                flux_surf_tune(l, j_k) / &
                                flux_surf_ref(l)**2
        ENDDO
!
        DO i=first_layer, last_layer
          DO l=1, Atm%n_profile
            a(i_k, j_k) = a(i_k, j_k) + hr_tune(l, i, i_k) * &
                                        hr_tune(l, i, j_k) * &
              (1.0_RealK/(offset+ABS(hr_ref(l, i)))**2) * &
              LOG(Atm%p_level(l, i)/Atm%p_level(l, i-1))
          ENDDO
        ENDDO 
!     
      ENDDO
!
      a(i_k, n_k+1) = 1.0_RealK
!
    ENDDO
!
!   Impose the constraint on the sum of the terms
    a(n_k+1, 1:n_k) = 1.0_RealK
    b(n_k+1)        = 1.0_RealK
!
!
!
  END SUBROUTINE form_svd_matrix
!
!
!
END PROGRAM tune_weights

