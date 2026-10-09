! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Program to generate correlated-k data for a spectral file.
!
PROGRAM corr_k

  USE realtype_rd, ONLY: RealK
  USE rad_pcf
  USE ck_fit_pcf, ONLY: ip_ck_none
  USE dimensions_pp_ucf
  USE dimensions_spec_ucf
  USE def_spectrum
  USE def_solarspec
  USE def_inst_flt
  USE def_std_io_icf
  USE def_hitran_record
  USE hitran_cnst, ONLY: read_parsum_dat

  IMPLICIT NONE
!
!
! Local scalars:
  INTEGER :: start_time(8) 
!   Start/finish of program
  INTEGER :: end_time(8)   
!   End of program
  INTEGER :: iu_lbl
!   Unit number for input of the LbL database in HITRAN format
  INTEGER :: iu_cia
!   Unit number for input of the CIA database in HITRAN format
  INTEGER :: ierr = i_normal
!   Error flag
  INTEGER :: ios
!   I/O error status
  INTEGER :: n_pt_pair
!   Number of pairs for pressure and temperature
  INTEGER :: n_p
!   Number of unique pressures
  INTEGER :: i_gas
!   Actual type of the gas to be considered
  INTEGER :: i_index
!   Index of the gas in the spectral file
  INTEGER :: n_selected_band
!   Number of bands of the spectral file selected for deriving
!   k-distributions
  INTEGER, Allocatable :: list_band(:)
!   List of selected bands
  INTEGER :: i_type_residual
!   Identifier of the type of residual to be minimized
  INTEGER :: i_scale_fnc
!   Identifier of the scaling function used
  INTEGER :: n_omp_threads
!   Number of OpenMP threads to use
!
  LOGICAL :: l_fit_line_data
!   Flag requiring the fitting of line data
  LOGICAL :: l_fit_frn_continuum
!   Flag requiring the fitting of foreign-broadened continuum data
  LOGICAL :: l_fit_self_continuum
!   Flag requiring the fitting of self-broadened continuum data
  LOGICAL :: l_fit_cont_data
!   Flag requiring the fitting of generalised continuum data
  LOGICAL :: l_access_HITRAN
!   Flag for HITRAN database file
  LOGICAL :: l_access_xsc
!   Flag for HITRAN cross-section file
  LOGICAL :: l_access_cia
!   Flag for HITRAN CIA file
  LOGICAL :: l_lbl_exist
!   Flag for lbl absorption coefficient file
!
  LOGICAL :: l_interactive
!   Flag for interactive operation
  LOGICAL :: include_instrument_response
!   Flag for including the resposnse of an observing instrument in
!   the k-distribution
  LOGICAL :: include_h2o_foreign_continuum
!   Flag to generate continuum data
  LOGICAL :: l_use_h2o_frn_param
!   Flag to use foreign broadened H2O continuum parametrisation
  LOGICAL :: l_use_h2o_self_param
!   Flag to use self broadened H2O continuum parametrisation
  LOGICAL :: l_scale_pT
!   Flag for deriving a scaling function in pressure and temperature
!
  LOGICAL :: l_cont_line_abs_weight
!   Flag to use line absorption data for weighting in continuum transmissions
  INTEGER :: i_gas_1
!   Actual type of first continuum gas to be considered
  INTEGER :: i_index_1
!   Index of first continuum gas in the spectral file
  INTEGER :: i_gas_2
!   Actual type of second continuum gas to be considered
  INTEGER :: i_index_2
!   Index of second continuum gas in the spectral file
!
  CHARACTER (LEN=256) :: file_spectral
!   Name of the spectral file
  CHARACTER (LEN=1) :: char_if
!   Character response variable
!
  CHARACTER  (LEN=132) :: file_k
!   File for output of the k-values
  CHARACTER  (LEN=132) :: file_monitor
!   File for output of monitoring information
  CHARACTER  (LEN=132) :: file_lbl
!   File for output of LbL absorption coefficients
!
  REAL  (RealK) :: p(npd_pt)
!   Pressure at which correlated-k data will be generated
  REAL  (RealK) :: t(npd_pt)
!   Temperature at which correlated-k data will be generated
  REAL  (RealK), Allocatable :: p_ref(:)
!   Reference pressures for the scaling in each band
  REAL  (RealK), Allocatable :: t_ref(:)
!   Reference temperature for the scaling in each band
!
  INTEGER :: i_line_prof_corr
!   Line profile correction type
!
  LOGICAL :: l_self_broadening
!   Flags to include effects of self-broadening
  INTEGER :: n_gas_frac
!   Number of gas fractions at which to tabulate ESFT terms
  REAL  (RealK) :: gas_frac(npd_gas_frac)
!   List of gas fractions at which to tabulate ESFT terms
!
  TYPE  (StrSolarSpec) :: SolarSpec
!   Solar Spectrum
!
!
  REAL  (RealK) :: nu_inc_0
!   Default spacing for line integration
  REAL  (RealK) :: line_cutoff
!   Cutoff for choosing lines
  LOGICAL :: l_ckd_cutoff
!   Line absorption is to be adjusted for consistency with the CKD
!   continuum
  INTEGER :: i_weight
!   Type of weighting
  INTEGER :: i_ck_fit
!   Type of correlated-k fit
  REAL  (RealK) :: tol
!   Tolerance for the c-k fit
  REAL  (RealK) :: max_path
!   Maximum pathlength to be considered
  REAL  (RealK) :: max_path_wgt
!   Maximum pathlength to be considered for the absorber used for weighting
!   in continuum transmissions
!
! Continuum data:
  REAL  (RealK) :: umin_c
!   Minimum pathlength for continuum absorption
  REAL  (RealK) :: umax_c
!   Maximum pathlength for continuum absorption
  INTEGER :: n_path_c
!   Number of pathlengths for continuum absorption
  INTEGER :: n_pp
!   Number of partial pressures for continuum absorption
!
! The k-fit
  INTEGER, Allocatable :: n_k(:)
!   Number of k-terms in the fit
  REAL  (RealK), Allocatable :: w_k(:, :)
!   Weights for each k-term
  REAL  (RealK), Allocatable :: k_ave(:, :)
!   Mean k-value across the band
  REAL  (RealK), Allocatable :: k_opt(:, :)
!   Optimal k-value across the band
  REAL  (RealK), Allocatable :: scale(:, :, :)
!   Parameters of the scaling functions
!
  REAL  (RealK), Allocatable :: k_opt_self(:)
!   Optimal k-value for self-broadened continuum
  REAL  (RealK), Allocatable :: k_opt_frn(:)
!   Optimal k-value for foreign-broadened continuum
  REAL  (RealK), Allocatable :: scale_cont(:, :, :)
!   Parameters of the scaling function for the continuum
!
! Mapping
  LOGICAL :: l_load_map
!   Use pre-defined mapping of wavenumbers to g-space
  LOGICAL :: l_load_wgt
!   Use pre-defined k-term weights
  LOGICAL :: l_save_map
!   Save mapping of wavenumbers to g-space
  CHARACTER(LEN=132) :: file_map
!   Name of file with mapping
!
  INTEGER :: iu_k_out
!   Unit number for output of the k-fit
  INTEGER :: iu_monitor
!   Unit number for output of detailed monitoring information
!
  TYPE  (StrSpecData) :: Spectrum
!   Spectral data
  TYPE  (StrFiltResp) :: filter
!   Instrumental response function
!
! External functions:
  LOGICAL :: set_interactive
!   Function to set the flag for interactive operation
  EXTERNAL &
    set_interactive
!
! External subroutines:
  INTERFACE
!
    SUBROUTINE corr_k_single &
!
      (i_gas, i_index, i_gas_1, i_index_1, i_gas_2, i_index_2, &
       n_band, n_selected_band, list_band, &
       n_pt_pair, n_p, p_calc, t_calc, p_ref, t_ref, l_scale_pt, &
       iu_lbl, iu_cia, nu_inc_0, line_cutoff, l_ckd_cutoff, &
       wavelength_long, wavelength_short, &
       n_band_exclude, index_exclude, &
       i_weight, solarspec, &
       include_h2o_foreign_continuum, &
       l_use_h2o_frn_param, l_use_h2o_self_param, l_cont_line_abs_weight, &
       l_access_HITRAN, l_access_xsc, l_access_cia, l_lbl_exist, &
       l_fit_line_data, l_fit_self_continuum, l_fit_frn_continuum, &
       l_fit_cont_data, n_path_c, umin_c, umax_c, n_pp, &
       include_instrument_response, filter, &
       i_line_prof_corr, l_self_broadening, n_gas_frac, gas_frac, &
       i_ck_fit, tol, max_path, max_path_wgt, &
       nd_k_term, n_k, w_k, k_ave, k_opt, &
       k_opt_self, k_opt_frn, &
       i_type_residual, i_scale_function, scale_vector, &
       scale_cont, &
       iu_k_out, file_k, iu_monitor, file_monitor, file_lbl, &
       l_load_map, l_load_wgt, l_save_map, file_map, &
       n_omp_threads, ierr &
      )
!
      USE realtype_rd
      USE def_solarspec
      USE def_inst_flt
!
      INTEGER, Intent(IN) :: n_band
      INTEGER, Intent(IN) :: n_selected_band
      INTEGER, Intent(IN), Dimension(:) :: list_band
      REAL  (RealK), Intent(IN), Dimension(:) :: wavelength_long
      REAL  (RealK), Intent(IN), Dimension(:) :: wavelength_short
      INTEGER, Intent(IN), Dimension(:) :: n_band_exclude
      INTEGER, Intent(IN), Dimension(:, :) :: index_exclude
      INTEGER, Intent(IN) :: i_gas
      INTEGER, Intent(IN) :: i_index
      INTEGER, Intent(IN) :: i_gas_1
      INTEGER, Intent(IN) :: i_index_1
      INTEGER, Intent(IN) :: i_gas_2
      INTEGER, Intent(IN) :: i_index_2
      INTEGER, Intent(IN) :: iu_lbl
      INTEGER, Intent(IN) :: iu_cia
      REAL  (RealK), Intent(IN) :: nu_inc_0
      INTEGER, Intent(IN) :: i_weight
      INTEGER, Intent(IN) :: i_ck_fit
      LOGICAL, Intent(IN) :: l_scale_pt
      REAL  (RealK), Intent(IN) :: tol
      REAL  (RealK), Intent(IN) :: max_path
      REAL  (RealK), Intent(IN) :: max_path_wgt
      REAL  (RealK), Intent(IN) :: line_cutoff
      LOGICAL, Intent(IN) :: l_ckd_cutoff
      INTEGER, Intent(IN) :: n_pt_pair
      INTEGER, Intent(IN) :: n_p
      REAL  (RealK), Intent(IN), Dimension(:) :: p_calc
      REAL  (RealK), Intent(IN), Dimension(:) :: t_calc
      REAL  (RealK), Intent(IN), Dimension(:) :: p_ref
      REAL  (RealK), Intent(IN), Dimension(:) :: t_ref
      TYPE  (StrSolarSpec), Intent(IN) :: SolarSpec
      LOGICAL, Intent(IN) :: include_h2o_foreign_continuum
      LOGICAL, Intent(IN) :: l_use_h2o_frn_param
      LOGICAL, Intent(IN) :: l_use_h2o_self_param
      LOGICAL, Intent(IN) :: l_cont_line_abs_weight
      LOGICAL, Intent(IN) :: l_access_HITRAN
      LOGICAL, Intent(IN) :: l_access_xsc
      LOGICAL, Intent(IN) :: l_access_cia
      LOGICAL, Intent(IN) :: l_lbl_exist
      LOGICAL, Intent(IN) :: l_fit_line_data
      LOGICAL, Intent(IN) :: l_fit_frn_continuum
      LOGICAL, Intent(IN) :: l_fit_self_continuum
      LOGICAL, Intent(IN) :: l_fit_cont_data
      INTEGER, Intent(IN) :: n_path_c
      REAL  (RealK), Intent(IN) :: umin_c
      REAL  (RealK), Intent(IN) :: umax_c
      INTEGER, Intent(IN) :: n_pp
      LOGICAL, Intent(IN) :: include_instrument_response
      TYPE  (StrFiltResp), Intent(IN) :: filter
      INTEGER, Intent(IN) :: i_line_prof_corr
      LOGICAL, Intent(IN) :: l_self_broadening
      INTEGER, Intent(IN) :: n_gas_frac
      REAL  (RealK), Intent(IN), Dimension(:) :: gas_frac
      INTEGER, Intent(IN) :: n_omp_threads
      INTEGER, Intent(INOUT) :: ierr
      INTEGER, Intent(IN) :: nd_k_term
      INTEGER, Intent(INOUT), Dimension(:) :: n_k
      REAL  (RealK), Intent(OUT), Dimension(:, :) :: w_k
      REAL  (RealK), Intent(OUT), Dimension(:, :) :: k_ave
      REAL  (RealK), Intent(OUT), Dimension(:, :) :: k_opt
      REAL  (RealK), Intent(OUT), Dimension(:) :: k_opt_self
      REAL  (RealK), Intent(OUT), Dimension(:) :: k_opt_frn
      INTEGER, Intent(IN) :: i_type_residual
      INTEGER, Intent(IN) :: i_scale_function
      REAL  (REALK), Intent(OUT), Dimension(:, :, :) :: scale_vector
      REAL  (REALK), Intent(OUT), Dimension(:, :, :) :: scale_cont
      INTEGER, Intent(IN) :: iu_k_out
      CHARACTER  (LEN=*), Intent(IN) :: file_k
      INTEGER, Intent(IN) :: iu_monitor
      CHARACTER  (LEN=*), Intent(IN) :: file_monitor
      CHARACTER  (LEN=*), Intent(IN) :: file_lbl
      LOGICAL, Intent(IN) :: l_load_map
      LOGICAL, Intent(IN) :: l_load_wgt
      LOGICAL, Intent(IN) :: l_save_map
      CHARACTER  (LEN=*), Intent(IN) :: file_map
!
    END SUBROUTINE corr_k_single
!
  END INTERFACE
!
!- End of header
!
!
!
! Set the flag for interactive operation
  l_interactive=set_interactive()
!
! Request the name of the LbL file.
  WRITE(iu_stdout, '(/a)') &
    'Give the name of the LbL absorption coefficient file.'
  READ(iu_stdin, '(a)') file_lbl
  INQUIRE(FILE=file_lbl, EXIST=l_lbl_exist)
!
! Acquire the line data-base.
  WRITE(iu_stdout, "(/a,/a,/a,/a,/a,/a)") &
    "Will a HITRAN database be provided? Type:", &
    "  L for line database (.par)", &
    "  B for bespoke line database (.bpar)", &
    "  X for cross-section database (.xsc)", &
    "  U for bespoke UV cross-section database (.uvxsc)", &
    "  N for none."
  DO
!
    READ(iu_stdin, '(A)') char_if
!
    IF ( (char_if == 'L') .OR. (char_if == 'l') .OR. &
         (char_if == 'Y') .OR. (char_if == 'y') ) THEN
!
      l_access_HITRAN=.TRUE.
      l_access_xsc=.FALSE.
      CALL get_free_unit(ierr, iu_lbl)
      IF (ierr /= i_normal) STOP
      CALL open_file_in(ierr, iu_lbl, &
        "Give the name of the HITRAN .par database.")
      IF (ierr /= i_normal) STOP
      CALL read_parsum_dat
      EXIT
!
    ELSE IF ( (char_if == 'B') .OR. (char_if == 'b') ) THEN
!
      l_access_HITRAN=.TRUE.
      l_access_xsc=.FALSE.
      ! A bespoke input format for data selected from hitran.org
      hitran_record_format = hitran_bespoke_frmt
      CALL get_free_unit(ierr, iu_lbl)
      IF (ierr /= i_normal) STOP
      CALL open_file_in(ierr, iu_lbl, &
        "Give the name of the bespoke HITRAN .bpar database.")
      IF (ierr /= i_normal) THEN
        WRITE(iu_err, '(A, i5)') 'Error in open_file_in: ', ierr        
        STOP
      END IF
      CALL read_parsum_dat
      EXIT
!
    ELSE IF ( (char_if == 'X') .OR. (char_if == 'x') ) THEN
!
      l_access_HITRAN=.FALSE.
      l_access_xsc=.TRUE.
      CALL get_free_unit(ierr, iu_lbl)
      IF (ierr /= i_normal) STOP
      CALL open_file_in(ierr, iu_lbl, &
        "Give the name of the HITRAN .xsc database.")
      IF (ierr /= i_normal) STOP
      EXIT
!
    ELSE IF ( (char_if == 'U') .OR. (char_if == 'u') ) THEN
!
      l_access_HITRAN=.FALSE.
      l_access_xsc=.TRUE.
      ! A bespoke input format is used for UV cross-section data
      xsc_header_format = uvxsc_header_frmt
      xsc_data_format = uvxsc_data_frmt
      CALL get_free_unit(ierr, iu_lbl)
      IF (ierr /= i_normal) THEN
        WRITE(iu_err, '(A, i5)') 'Error in get_free_unit: ', ierr
        STOP
      END IF
      CALL open_file_in(ierr, iu_lbl, &
        "Give the name of the .uvxsc database.")
      IF (ierr /= i_normal) THEN
        WRITE(iu_err, '(A, i5)') 'Error in open_file_in: ', ierr        
        STOP
      END IF
      EXIT
!
    ELSE IF ( (char_if == 'N') .OR. (char_if == 'n') ) THEN
!
      l_access_HITRAN=.FALSE.
      l_access_xsc=.FALSE.
      EXIT
!
    ELSE
!
      WRITE(iu_err, '(/A)') '*** Error: Unrecognized input'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(/A)') 'Please re-enter.'
      ELSE
        STOP
      ENDIF
!
    ENDIF
!
  ENDDO

! Acquire the line data-base.
  WRITE(iu_stdout, "(/a)") &
    "Will a HITRAN CIA (.cia) database be provided? (Y/N)"
  DO
!
    READ(iu_stdin, '(A)') char_if
!
    IF ( (char_if == 'Y') .OR. (char_if == 'y') ) THEN
!
      l_access_cia=.TRUE.
      CALL get_free_unit(ierr, iu_cia)
      IF (ierr /= i_normal) STOP
      CALL open_file_in(ierr, iu_cia, &
        "Give the name of the HITRAN .cia database.")
      IF (ierr /= i_normal) STOP
      EXIT
      EXIT
!
    ELSE IF ( (char_if == 'N') .OR. (char_if == 'n') ) THEN
!
      l_access_cia=.FALSE.
      EXIT
!
    ELSE
!
      WRITE(iu_err, '(/A)') '*** Error: Unrecognized input'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(/A)') 'Please re-enter.'
      ELSE
        STOP
      ENDIF
!
    ENDIF
!
  ENDDO

! Read in the spectral file.
  WRITE(*, "(/ a)") "Enter the name of the spectral file."
  READ(*, "(a)") file_spectral
  CALL read_spectrum(file_spectral, Spectrum)


! Read in instrument response if required
!
  WRITE(iu_stdout, "(/A)") &
    "Is an instrumental response required? (Y/N)"
  Inst: DO
!
    READ(iu_stdin, '(A)') char_if
!
    IF ( (char_if == 'Y') .OR. (char_if == 'y') ) THEN
!
      include_instrument_response=.TRUE.
      CALL read_instrument_response_90(filter, ierr)
      IF (ierr /= i_normal) STOP
      EXIT
!
    ELSE IF ( (char_if == 'N') .OR. (char_if == 'n') ) THEN
!
      include_instrument_response=.FALSE.
      EXIT
!
    ELSE
!
      WRITE(iu_err, '(/A)') '*** Error: Unrecognized input'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(/A)') 'Please re-enter.'
      ELSE
        STOP
      ENDIF
!
    ENDIF
!
  ENDDO Inst
!
! Allocate arrays with spectrally dependent sizes.
  ALLOCATE(list_band(Spectrum%Dim%nd_band))
  ALLOCATE(p_ref(Spectrum%Dim%nd_band))
  ALLOCATE(t_ref(Spectrum%Dim%nd_band))
  ALLOCATE(n_k(Spectrum%Dim%nd_band))
!
! Select the gas to be considered, the range of bands and
! the pressures and temperatures.
  CALL set_condition_ck_90(l_interactive, &
    npd_pt, n_pt_pair, n_p, p, t, &
    Spectrum%Dim%nd_band, Spectrum%Dim%nd_species, &
    Spectrum%Gas%n_absorb, Spectrum%Gas%type_absorb, &
    i_gas, i_index, i_gas_1, i_index_1, i_gas_2, i_index_2, &
    Spectrum%Basic%n_band, Spectrum%Gas%n_band_absorb, &
    Spectrum%Gas%index_absorb, &
    Spectrum%Dim%nd_continuum, Spectrum%Cont%n_band_continuum, &
    Spectrum%Cont%index_continuum, &
    l_fit_line_data, l_fit_frn_continuum, l_fit_self_continuum, &
    l_fit_cont_data, umin_c, umax_c, n_path_c, n_pp, &
    l_access_HITRAN, l_access_xsc, l_access_cia, &
    include_h2o_foreign_continuum, &
    l_use_h2o_frn_param, l_use_h2o_self_param, l_cont_line_abs_weight, &
    n_selected_band, list_band, &
    i_ck_fit, tol, max_path, max_path_wgt, n_k, nu_inc_0, line_cutoff, &
    l_ckd_cutoff, l_scale_pT, i_type_residual, i_scale_fnc, p_ref, t_ref, &
    l_load_map, l_load_wgt, l_save_map, file_map, &
    i_line_prof_corr, l_self_broadening, n_gas_frac, gas_frac, npd_gas_frac, &
    ierr)
!
! Allocate arrays for the k-fit, now that the size of the scaling 
! vector is known.
  ALLOCATE(w_k(npd_k_term, Spectrum%Dim%nd_band))
  ALLOCATE(k_opt(npd_k_term, Spectrum%Dim%nd_band))
  ALLOCATE(k_ave(npd_k_term, Spectrum%Dim%nd_band))
  ALLOCATE(scale(n_scale_variable(i_scale_fnc), &
    npd_k_term, Spectrum%Dim%nd_band))
  ALLOCATE(k_opt_self(Spectrum%Dim%nd_band))
  ALLOCATE(k_opt_frn(Spectrum%Dim%nd_band))
  ALLOCATE(scale_cont(MAX(n_scale_variable(i_scale_fnc), 1), 1, &
    Spectrum%Dim%nd_band))

  IF ((i_ck_fit /= ip_ck_none) .OR. &
      l_fit_self_continuum .OR. l_fit_frn_continuum) THEN
!   Select the weighting to be applied.
    CALL select_weight_ck_90(i_weight, SolarSpec, l_interactive, ierr)
!   
!   Set the output file.
    CALL get_free_unit(ierr, iu_k_out)
    IF (ierr /= i_normal) STOP
    CALL open_file_out_90(iu_k_out, l_interactive, &
      "Give the name of the output file.", &
      file_k, ierr)
    IF (ierr /= i_normal) STOP
  END IF
    
! Define the output file of detailed monitoring information.
  CALL get_free_unit(ierr, iu_monitor)
  IF (ierr /= i_normal) STOP
  CALL open_file_out_90(iu_monitor, l_interactive, &
    "Give the name of the monitoring file.", &
    file_monitor, ierr)
  IF (ierr /= i_normal) STOP

! The opening routine opens the file generically; however, it is more
! useful here to open the file afresh for each band, so it is closed
! first. (But the call was required to get the name of the file.)
  IF ((i_ck_fit /= ip_ck_none) .OR. &
      l_fit_self_continuum .OR. l_fit_frn_continuum) THEN
    CLOSE(iu_k_out)
  END IF
  CLOSE(iu_monitor)

! Aquire number of OpenMP threads to use
  WRITE(iu_stdout, "(/A)") &
    "Specify the number of OpenMP threads to use."
  DO
!
    READ(iu_stdin, *, IOSTAT=ios) n_omp_threads
!
    IF ( ios == 0 .AND. n_omp_threads >= 1 ) THEN
      EXIT
    ELSE
      WRITE(iu_err, '(/A)') '*** Error: Invalid input'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(/A)') 'Please re-enter.'
      ELSE
        STOP
      ENDIF
    ENDIF
!
  ENDDO
!
!
!
  WRITE(*,"(a)")"==============================="
  WRITE(*,"(a)")"corr_k : Execution starts "
  CALL DATE_AND_TIME(VALUES=start_time)
  WRITE(*,"(a4,i2.2,a1,i2.2,a1,i2.2,a4,i2.2,a1,i2.2,a1,i4)")     &
       " at ",start_time(5),":",start_time(6),":",start_time(7), &
       " on ",start_time(3),"/",start_time(2),"/",start_time(1)
!
  CALL corr_k_single(i_gas, i_index, i_gas_1, i_index_1, i_gas_2, i_index_2, &
    Spectrum%Basic%n_band, n_selected_band, list_band, &
    n_pt_pair, n_p, p, t, p_ref, t_ref, l_scale_pt, &
    iu_lbl, iu_cia, nu_inc_0, line_cutoff, l_ckd_cutoff, &
    Spectrum%Basic%wavelength_long, Spectrum%Basic%wavelength_short, &
    Spectrum%Basic%n_band_exclude, Spectrum%Basic%index_exclude, &
    i_weight, solarspec, &
    include_h2o_foreign_continuum, &
    l_use_h2o_frn_param, l_use_h2o_self_param, l_cont_line_abs_weight, &
    l_access_HITRAN, l_access_xsc, l_access_cia, l_lbl_exist, &
    l_fit_line_data, l_fit_self_continuum, l_fit_frn_continuum, &
    l_fit_cont_data, n_path_c, umin_c, umax_c, n_pp, &
    include_instrument_response, filter, &
    i_line_prof_corr, l_self_broadening, n_gas_frac, gas_frac, &
    i_ck_fit, tol, max_path, max_path_wgt, &
    npd_k_term, n_k, w_k, k_ave, k_opt, &
    k_opt_self, k_opt_frn, &
    i_type_residual, i_scale_fnc, scale, scale_cont, &
    iu_k_out, file_k, iu_monitor, file_monitor, file_lbl, &
    l_load_map, l_load_wgt, l_save_map, file_map, &
    n_omp_threads, ierr )

  IF (include_instrument_response) THEN
    DEALLOCATE(filter%wavenumber)
    DEALLOCATE(filter%response)
    DEALLOCATE(filter%d2_response)
  ENDIF

  CLOSE(iu_k_out)
  IF (l_access_HITRAN .OR. l_access_xsc) CLOSE(iu_lbl)
  IF (l_access_cia) CLOSE(iu_cia)

  WRITE(*,"(a)")"==============================="
  WRITE(*,"(a)")"corr_k : Execution ends   "
  CALL DATE_AND_TIME(VALUES=end_time)
  WRITE(*,"(a4,i2.2,a1,i2.2,a1,i2.2,a4,i2.2,a1,i2.2,a1,i4)") &
       " at ",end_time(5),":",end_time(6),":",end_time(7),   &
       " on ",end_time(3),"/",end_time(2),"/",end_time(1)
  WRITE(*,"(a)")"==============================="

END PROGRAM corr_k
