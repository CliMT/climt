! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to generate c-k data for one gas and band.
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
 i_type_residual, i_scale_function, scale_vector, scale_cont, &
 iu_k_out, file_k, iu_monitor, file_monitor, file_lbl, &
 l_load_map, l_load_wgt, l_save_map, file_map, &
 n_omp_threads, ierr &
)

! Modules used:
  USE realtype_rd
  USE rad_pcf
  USE rad_ccf
  USE def_solarspec
  USE def_std_io_icf
  USE def_inst_flt
  USE def_hitran_record
  USE ck_parm_acf
  USE ck_fit_pcf
  USE hitran_cnst, ONLY: atomic_mass_unit, molar_gas_constant
  USE gas_list_pcf
  USE weighting_pcf, ONLY: ip_weight_planck, ip_weight_d_planck
  USE h2o_continuum, ONLY: self_continuum, foreign_continuum, sat_vap_press
  USE line_prof_corr_mod, ONLY: line_prof_corr, set_line_prof_corr_cnst
  USE errormessagelength_mod, ONLY: errormessagelength
  USE ereport_mod, ONLY: ereport

  IMPLICIT NONE 


  TYPE StrLineParam
    REAL  (RealK) :: line_centre
    REAL  (RealK) :: S_adj
    REAL  (RealK) :: alpha_lorentz
    REAL  (RealK) :: alpha_lorentz_air
    REAL  (RealK) :: alpha_lorentz_self
    REAL  (RealK) :: alpha_doppler
  END TYPE StrLineParam
!
  TYPE (StrLineParam), Allocatable, Dimension(:) :: adj_line_parm
!   Adjusted line parameters
!
  REAL  (RealK) :: upper_cutoff
  REAL  (RealK) :: lower_cutoff
!
!
!
! Dummy arguments:
  INTEGER, Intent(IN) :: n_band
!   Number of bands in the spectral file
  INTEGER, Intent(IN) :: n_selected_band
!   Number of bands selected for integration
  INTEGER, Intent(IN), Dimension(:) :: list_band
!   List of bands to be considered
  REAL  (RealK), Intent(IN), Dimension(:) :: wavelength_long
!   Longer wavelength limits of spectral bands
  REAL  (RealK), Intent(IN), Dimension(:) :: wavelength_short
!   Shorter wavelength limits of spectral bands
  INTEGER, Intent(IN), Dimension(:)    :: n_band_exclude
!   Number of exclusions from each band
  INTEGER, Intent(IN), Dimension(:, :) :: index_exclude
!   List of excluded bands within each region
!
  INTEGER, Intent(IN) :: i_gas
!   Identifier for gas
  INTEGER, Intent(IN) :: i_index
!   Index of the gas in the spectral file
  INTEGER, Intent(IN) :: i_gas_1
!   Identifier for first continuum gas
  INTEGER , Intent(IN) :: i_index_1
!   Index of first continuum gas in the spectral file
  INTEGER, Intent(IN) :: i_gas_2
!   Identifier for second continuum gas
  INTEGER , Intent(IN) :: i_index_2
!   Index of second continuum gas in the spectral file
!
  INTEGER, Intent(IN) :: iu_lbl
!   Unit number for input from the database
  INTEGER, Intent(IN) :: iu_cia
!   Unit number for input from the HITRAN CIA database
!
! Options for generating data
  LOGICAL, Intent(IN) :: l_fit_line_data
!   Controlling flag to fit line data
  LOGICAL, Intent(IN) :: l_fit_self_continuum
!   Controlling flag to fit self-broadened data
  LOGICAL, Intent(IN) :: l_fit_frn_continuum
!   Controlling flag to fit foreign-broadened continuum data
  LOGICAL, Intent(IN) :: l_fit_cont_data
!   Controlling flag to fit generalised continuum data
  LOGICAL, Intent(IN) :: l_access_HITRAN
!   Logical set to true if program has access to HITRAN LbL database
  LOGICAL, Intent(IN) :: l_access_xsc
!   Logical set to true if program has access to HITRAN XSC database
  LOGICAL, Intent(IN) :: l_access_cia
!   Logical set to true if program has access to HITRAN CIA database
  LOGICAL, Intent(IN) :: l_lbl_exist
!   Flag for lbl absorption coefficient file
!
  REAL  (RealK), Intent(IN) :: nu_inc_0
!   Default spacing for line integration
  INTEGER, Intent(IN) :: i_weight
!   Type of weighting
  INTEGER, Intent(IN) :: i_ck_fit
!   Type of correlated-k fit
  LOGICAL, Intent(IN) :: l_scale_pt
!   Flag to carry out pressure and temperature scaling
  REAL  (RealK), Intent(IN) :: tol
!   Tolerance for fitting k-terms
  REAL  (RealK), Intent(IN) :: max_path
!   Maximum pathlength to be considered for the absorber
  REAL  (RealK), Intent(IN) :: max_path_wgt
!   Maximum pathlength to be considered for the absorber used for weighting
!   in continuum transmissions
  REAL  (RealK), Intent(IN) :: line_cutoff
!   Cutoff for choosing lines
  LOGICAL, Intent(IN) :: l_ckd_cutoff
!   Cutoff is to be consistent with the CKD continuum model
!
! Variables required for continuum data
  INTEGER, Intent(IN) :: n_path_c
!   Number of paths required for continuum data
  REAL  (RealK), Intent(IN) :: umin_c
!   Minimum path length fopr continuum data
  REAL  (RealK), Intent(IN) :: umax_c
!   Maximum path length fopr continuum data
  INTEGER, Intent(IN) :: n_pp
!   Number of partial pressures
!
! Scaling conditions
  INTEGER, Intent(IN) :: n_pt_pair
!   Number of pairs of pressure and temperature
  INTEGER, Intent(IN) :: n_p
!   Number of unique pressures
  REAL  (RealK), Intent(INOUT), Dimension(:) :: p_calc
!   Pressures for c-k calculation
  REAL  (RealK), Intent(INOUT), Dimension(:) :: t_calc
!   Temperatures for c-k calculation
!
  REAL  (RealK), Intent(IN), Dimension(:) :: p_ref
!   Reference pressures for scaling in each band
  REAL  (RealK), Intent(IN), Dimension(:) :: t_ref
!   Reference temperatures for scaling in each band
!
  TYPE  (StrSolarSpec), Intent(IN) :: SolarSpec
!   Solar Spectrum
!
  LOGICAL, Intent(IN) :: include_h2o_foreign_continuum
!   Flag to include the foreign-broadened H2O continuum (with a 
!   partial pressure of 0) with the line data
  LOGICAL, Intent(IN) :: l_use_h2o_frn_param
!   Flag to use foreign broadened H2O continuum parametrisation
  LOGICAL, Intent(IN) :: l_use_h2o_self_param
!   Flag to use self broadened H2O continuum parametrisation
  LOGICAL, Intent(IN) :: l_cont_line_abs_weight
!   Flag to use line absorption data for weighting in continuum transmissions
!
  LOGICAL, Intent(IN) :: include_instrument_response
!   Flag to include the instrumental response function
  TYPE  (StrFiltResp), Intent(IN) :: filter
!   Instrumental response function
!
  INTEGER, Intent(IN) :: i_line_prof_corr
!   Line profile correction type
!
! Self-broadening options
  LOGICAL, Intent(IN) :: l_self_broadening
!   Flags to include effects of self-broadening
  INTEGER, Intent(IN) :: n_gas_frac
!   Number of gas fractions at which to tabulate ESFT terms
  REAL  (RealK), Intent(IN), Dimension(:) :: gas_frac
!   List of gas fractions at which to tabulate ESFT terms
!
  INTEGER, Intent(IN) :: n_omp_threads
!   Number of OpenMP threads to use
!
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
!
! The k-fit
  INTEGER, Intent(IN) :: nd_k_term
!   Size allocated for k-terms
  INTEGER, Intent(INOUT), Dimension(:) :: n_k
!   Number of k-terms in the fit
  REAL  (RealK), Intent(OUT), Dimension(:, :) :: w_k
!   Weights for each k-term
  REAL  (RealK), Intent(OUT), Dimension(:, :) :: k_ave
!   Mean k-value across the band
  REAL  (RealK), Intent(OUT), Dimension(:, :) :: k_opt
!   Optimal k-value across the band
  REAL  (RealK), Intent(OUT), Dimension(:), Target :: k_opt_self
!   Optimal self-broadened continuum absorption coefficient
  REAL  (RealK), Intent(OUT), Dimension(:), Target :: k_opt_frn
!   Optimal foreign-broadened continuum absorption coefficient
!
! The scaling function:
  INTEGER, Intent(IN) :: i_type_residual
!   Type of residual used
  INTEGER, Intent(IN) :: i_scale_function
!   Type of scaling function
  REAL  (REALK), Intent(OUT), Dimension(:, :, :) :: scale_vector
!   Coefficients of the scaling function
  REAL  (REALK), Intent(OUT), Dimension(:, :, :) :: scale_cont
!   Coefficients of the scaling function for the continuum
!
  INTEGER, Intent(IN) :: iu_k_out
!   Unit number for output of the k-fit
  CHARACTER  (LEN=*), Intent(IN) :: file_k
!   File for output of the k-fit
  INTEGER, Intent(IN) :: iu_monitor
!   Unit number for output of detailed monitoring information
  CHARACTER  (LEN=*), Intent(IN) :: file_monitor
!   File for output of the monitoring information
  CHARACTER  (LEN=*), Intent(IN) :: file_lbl
!   File for output of the lbl absorption coefficients
!
  LOGICAL, Intent(IN) :: l_load_map
!   Use pre-defined mapping of wavenumbers to g-space
  LOGICAL, Intent(IN) :: l_load_wgt
!   Use pre-defined k-term weights
  LOGICAL, Intent(IN) :: l_save_map
!   Save mapping of wavenumbers to g-space
  CHARACTER(LEN=*), Intent(IN) :: file_map
!   Name of file with mapping
!
!
! Local variables
  INTEGER :: alloc_status
!   Status flag for allocation of HITRAN arrays
  INTEGER :: i, j, k, ipb
!   Loop variable
  INTEGER :: ibb
!   Loop variable for bands
  INTEGER :: ib
!   Loop variable for bands
  INTEGER :: isb
!   Loop variable for sub-bands
  INTEGER :: ip, it, ipt
!   Loop variable for pressure and temperature
  INTEGER :: igf
!   Loop variable for gas fraction
  INTEGER :: jx
!   Loop variable for k-terms
  INTEGER :: ik
!   Loop variable for excluded regions
  INTEGER :: h_gas
!   Identifier for the gas in HITRAN
  INTEGER :: h_isotopes(npd_isotopes)
!   List of isotopes to be read
  INTEGER :: num_lines_in_band
!   Number of lines extracted from the database for the gas
  INTEGER :: num_cia_lines_in_band
!   Number of lines extracted from the CIA database for the gas
  REAL  (RealK), Allocatable, Dimension(:) :: k_cutoff
!   Absorption by each line at its cutoff
!
  INTEGER :: i0
!   Pointer to position in weighting array
  INTEGER :: i1
!   Pointer to position in weighting array
  INTEGER :: n_nu_k
!   Number of frequency points in the k-interval
  INTEGER :: n_path
!   Number of absorber paths
  INTEGER, Target :: ig(0: nd_k_term)
  INTEGER :: ig_tmp(0: nd_k_term)
!   Pointers to g-quadrature points
  INTEGER :: n_k_last, n_pre
  INTEGER :: pstart(n_pt_pair), pend(n_pt_pair), pcount(n_pt_pair)
  REAL (RealK) :: pmax, pmin, pstep
!
  TYPE(StrHitranRec), DIMENSION(:), ALLOCATABLE :: hitran_data
!   Information on lines read from the Hitran line database
  TYPE(StrXscRec), DIMENSION(:), ALLOCATABLE :: xsc
!   Information on cross-sections read from the Hitran database
  TYPE(StrCIARec), DIMENSION(:), ALLOCATABLE :: cia
!   Information on CIA read from the Hitran database
!
! Weighting:
  INTEGER :: n_nu, n_nu_written, hitran_lines(n_selected_band)
!   Number of frequencies for weighting
  REAL  (RealK), Allocatable :: nu_wgt(:), nu_wgt_all(:)
!   Frequencies of weighting points
  REAL  (RealK), Allocatable :: kabs(:), kabs_all(:,:), kabs_all_sb(:,:,:)
!   Absorption coefficients at weighting frequencies
  REAL  (RealK), Allocatable :: kabs_lines(:), kabs_all_lines(:,:)
!   Line absorption coefficients at weighting frequencies
  REAL  (RealK), Allocatable :: kabs_rank1(:), kabs_rank2(:)
!   Rankings for absorption coeffs at weighting freqs
  REAL  (RealK) :: rank_work
  INTEGER, Allocatable :: kabs_mask(:)
!   Mask for binning absorption coefficients by their scaling behaviour
  REAL  (RealK), Allocatable :: k_for(:)
!   Foreign continuum absorption coefficients at weighting frequencies
  INTEGER, Allocatable :: map(:), gmap(:), pmap(:), kmap(:)
!   Mapping array in frequency space used to order the absorption
!   coefficients
  REAL  (RealK), Allocatable :: wgt(:), wgt_sum(:), wgt_ref(:), wgt_sv(:)
!   Weightings at different frequencies
  REAL  (RealK), Allocatable :: wgt_unsorted(:)
!   Copy of the unsorted weights to calculate sub-band weights
  REAL  (RealK), Allocatable :: wgt_k(:)
!   Product of weighting and absorption coefficients
  REAL  (RealK) :: integ_wgt, integ_wgt_tmp
!   Integral of the weighting function across the whole band
  REAL  (RealK) :: integ_k
!   Integral of the weighting function across the part of the band
!   represented by the k-term
  REAL  (RealK) :: response_0
!   Value of the response function at the point of evaluation
!
  REAL  (RealK), Allocatable :: band_min(:)
!   Low frequency limit of bands in m-1
  REAL  (RealK), Allocatable :: band_max(:)
!   High frequency limit of bands in m-1
  REAL  (RealK) :: band_width
!   Width of bands in wavenumbers, after removing excluded regions
  REAL  (RealK) :: nu_band_adjust
!   Wavenumber used for adjusting/checking band limits
  REAL  (RealK) :: lower_band_limit
!   Lower limit of region to be read
  REAL  (RealK) :: upper_band_limit
!   Upper limit of region to be read
  REAL  (RealK) :: nu_start
!   Lowest frequency of the current segment of the band
  REAL  (RealK) :: nu_end
!   Highest frequency of the current segment of the band
  REAL  (RealK), Pointer :: trans_ref(:)
!   Transmissions at the reference conditions across the whole band
!   including all k-terms
  REAL  (RealK), Pointer :: trans_pt_k(:, :, :)
!   Actual transmissions for a single k-term only (used when fitting
!   a scaling function) at the supplied range of pressures and 
!   temperatures
  REAL  (RealK), Pointer :: trans_calc(:)
!   Calculated transmissions
  REAL  (RealK) :: nu_inc
!   Actual increment applied across the band
  REAL  (RealK) :: kabs_line
!   Contribution to absorption from a single line

! Sub-bands
  INTEGER :: n_sub_band
  INTEGER, Allocatable :: sub_band_k(:)
  REAL (RealK) :: sub_band_w_sum
  REAL (RealK), Allocatable :: sub_band_w(:)
  REAL (RealK), Allocatable :: sub_band_nu_short(:)
  REAL (RealK), Allocatable :: sub_band_nu_long(:)

! Continuum variables:
  REAL  (RealK), Dimension(:), Allocatable :: u_c
!   Pathlengths for continuum absorption
  REAL  (RealK), Pointer, Dimension(:, :) :: u_fit_c
!   Products of the mass of absorber and the partial pressure in 
!   continua
  REAL  (RealK), Dimension(:), Allocatable :: trans_c
!   Continuum transmissions
  REAL  (RealK), Pointer, Dimension(:, :) :: trans_fit_c
!   Continuum transmissions to be fitted
  REAL  (RealK), Pointer, Dimension(:) :: trans_app_c
!   Approximated continuum transmission
!
  REAL  (RealK), Allocatable :: rdiff(:)
!   Array of squared differences from reference conditions
  INTEGER :: ipt_ref, index_pt_ref(n_p), n_t(n_p), ipoint
  LOGICAL :: mask_pt_ref(n_pt_pair)
!   Pointer for reference pressure/temperature

  REAL  (RealK) :: umin
!   Minimum pathlength of absorber
  REAL  (RealK) :: umax
!   Maximum pathlength of absorber
  REAL  (RealK), Dimension(:, :), Allocatable :: u_l
  REAL  (RealK), Dimension(:), Allocatable :: u_l_ref
!   Pathlengths of absorber
  LOGICAL :: l_wgt_scale_sqrt
  REAL  (RealK) :: u_wgt_scale
!   In order to apply a line transmission weighting to continuum absorption
!   transmissions the column mass of the weighting gas is needed. This is
!   calculated as u_gas = max_path_wgt*sqrt(u_cont/max_path) for  self-
!   broadened continua (l_wgt_scale_sqrt == .TRUE.) and 
!   as u_gas = max_path_wgt*u_cont/max_path if a foreign-broadened continuum
!   (l_wgt_scale_sqrt == .FALSE.). u_wgt_scale is used to store the quantity
!   max_path_wgt/sqrt(max_path) or max_path_wgt/max_path.

  INTEGER :: i_index_c
!   Variable indicating the type of the continuum
  REAL  (RealK) :: k_ave_tmp(nd_k_term)
!   Mean k-value across the band: temporary value at conditions 
!   other than the reference
  REAL  (RealK) :: k_opt_tmp(nd_k_term)
!   Optimal k-value across the band: temporary value at conditions
!   other than the reference
  REAL  (RealK) :: kopt_all(nd_k_term, n_pt_pair, n_selected_band)
  REAL  (RealK) :: kopt_all_sb(nd_k_term, n_pt_pair, &
                               n_gas_frac, n_selected_band)
  REAL  (RealK), ALLOCATABLE :: k_lookup(:,:,:,:)
  REAL  (RealK), ALLOCATABLE :: p_lookup(:), t_lookup(:,:)
  REAL  (RealK) :: ln_p_calc(n_pt_pair), max_p_calc
  REAL  (RealK), Pointer :: k_cont
!   Continuum absorption coefficient

  LOGICAL :: l_transparent_fit
!   Logical is true if a transparent fit should be applied

  LOGICAL :: l_calc_cont
!   Logical is true if calculating continuum data is required

  INTEGER :: ncidin_lbl,ncidout_lbl
!   NetCDF file identifiers for input and output lbl files

  INTEGER :: ncidin_map,ncidout_map
!   NetCDF file identifiers for input and output mapping files

  REAL  (RealK) :: err_norm
!   Error norm for fitting
  REAL  (RealK) :: rms_residual
!   R.m.s. residual error in transmission
  REAL  (RealK) :: grid_offset
!   Offset of externally supplied wavenumber grid

  REAL :: start_band, finish_band, timer1, timer2
!   Timers
  
  LOGICAL :: l_debug = .FALSE.
!  LOGICAL :: l_debug = .TRUE.
  LOGICAL :: l_output_reference_weight = .FALSE.

  CHARACTER (LEN=errormessagelength) :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'CORR_K_SINGLE'

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
    SUBROUTINE map_heap_func(a, map)
!
      USE realtype_rd
!
      REAL  (RealK), Intent(IN), Dimension(:) :: a
!
      INTEGER, Intent(OUT), Dimension(:) :: map
!
    END SUBROUTINE map_heap_func
!
!
    FUNCTION trans_k_dist(n_nu, k, nu_inc, wgt, integ_wgt, n_path, u) &
      RESULT (trans)
!
      USE realtype_rd
!
      INTEGER, Intent(IN) :: n_nu
      INTEGER, Intent(IN) :: n_path
      REAL  (RealK), Intent(IN), Dimension(n_nu) :: k
      REAL  (RealK), Intent(IN) :: nu_inc
      REAL  (RealK), Intent(IN), Dimension(n_nu) :: wgt
      REAL  (RealK), Intent(IN) :: integ_wgt
      REAL  (RealK), Intent(IN), Dimension(n_path) :: u
!
      REAL  (RealK), Dimension(n_path) :: trans

!
    END FUNCTION trans_k_dist
!
!
    SUBROUTINE set_g_point_90(n_nu, nu_inc, kabs, wgt, integ_wgt, &
      i_ck_fit, tol, max_path, l_kabs_wgt, kabs_wgt, &
      l_wgt_scale_sqrt, u_wgt_scale, nd_k_term, iu_monitor, &
      n_k, w_k, k_opt, k_ave, ig, ierr)
!
      USE realtype_rd
!
      INTEGER, Intent(IN) :: n_nu
      INTEGER, Intent(IN) :: nd_k_term
      INTEGER, Intent(IN) :: iu_monitor
      INTEGER, Intent(IN) :: i_ck_fit
      INTEGER, Intent(INOUT) :: ig(0:nd_k_term)
      INTEGER, Intent(INOUT) :: n_k
      INTEGER, Intent(INOUT) :: ierr
      REAL  (RealK), Intent(IN) :: nu_inc
      REAL  (RealK), Intent(IN) :: kabs(:)
      REAL  (RealK), Intent(IN) :: wgt(:)
      REAL  (RealK), Intent(IN) :: integ_wgt
      REAL  (RealK), Intent(IN) :: tol
      REAL  (RealK), Intent(IN) :: max_path
      REAL  (RealK), Intent(IN) :: kabs_wgt(:)
      REAL  (RealK), Intent(IN) :: u_wgt_scale
      REAL  (RealK), Intent(OUT) :: w_k(:)
      REAL  (RealK), Intent(OUT) :: k_opt(:)
      REAL  (RealK), Intent(OUT) :: k_ave(:)
      LOGICAL, Intent(IN) :: l_kabs_wgt
      LOGICAL, Intent(IN) :: l_wgt_scale_sqrt
!
    END SUBROUTINE set_g_point_90
!
    SUBROUTINE optimal_k &
      (n_nu, nu_inc, k, wgt, integ_k, k_mean, tol, &
       l_k_wgt, k_wgt, l_wgt_scale_sqrt, u_wgt_scale, &
       k_opt, error, ierr)
!
      USE realtype_rd
!
      INTEGER, Intent(IN) :: n_nu
      INTEGER, Intent(INOUT) :: ierr
      REAL  (RealK), Intent(IN) :: nu_inc
      REAL  (RealK), Intent(IN) :: k(:)
      REAL  (RealK), Intent(IN) :: wgt(:)
      REAL  (RealK), Intent(IN) :: integ_k
      REAL  (RealK), Intent(IN) :: k_mean
      REAL  (RealK), Intent(IN) :: tol
      REAL  (RealK), Intent(IN) :: k_wgt(:)
      REAL  (RealK), Intent(IN) :: u_wgt_scale
      REAL  (RealK), Intent(OUT) :: k_opt
      REAL  (RealK), Intent(OUT) :: error
      LOGICAL, Intent(IN) :: l_k_wgt
      LOGICAL, Intent(IN) :: l_wgt_scale_sqrt
!
    END SUBROUTINE optimal_k
!
!
    SUBROUTINE scale_ck_fit_90 &
!
      (ierr, iu_monitor, &
       n_other, &
       n_path, u_set, trans_set, p_set, t_set, &
       p_0, t_0, k_ref, &
       i_type_residual, i_scale_function, scale_vector, &
       rms_residual)
!
      USE realtype_rd
!
      INTEGER, Intent(INOUT) :: ierr
      INTEGER, Intent(IN) :: iu_monitor
      INTEGER, Intent(IN) :: i_type_residual
      INTEGER, Intent(IN) :: i_scale_function
      INTEGER, Intent(IN) :: n_other
      INTEGER, Intent(IN) :: n_path
!
      REAL  (RealK), Intent(IN) :: p_0
      REAL  (RealK), Intent(IN) :: t_0
      REAL  (RealK), Intent(IN) :: k_ref
      REAL  (RealK), Intent(IN), Dimension(:, :) :: u_set
      REAL  (RealK), Intent(IN), Dimension(:, :) :: trans_set
      REAL  (RealK), Intent(IN), Dimension(:) :: p_set
      REAL  (RealK), Intent(IN), Dimension(:) :: t_set
!
      REAL  (RealK), Intent(OUT), Dimension(:) :: scale_vector
      REAL  (RealK), Intent(OUT) :: rms_residual
!
    END SUBROUTINE scale_ck_fit_90
!
!
    SUBROUTINE write_fit_90 &
      (iu_k_out, l_continuum, l_cont_gen, l_self_broadening, i_band, &
       i_gas, i_index, i_index_1, i_index_2, &
       p, t, n_points, amount, transmittance, trans_calc, &
       n_k, k, w_k, i_scale, i_scale_function, scale_vector)
!
      USE realtype_rd
!
      INTEGER, Intent(IN) :: iu_k_out
!
      INTEGER, Intent(IN) :: i_band
      INTEGER, Intent(IN) :: i_gas
      INTEGER, Intent(IN) :: i_index
      INTEGER, Intent(IN) :: i_index_1
      INTEGER, Intent(IN) :: i_index_2
      INTEGER, Intent(IN) :: n_points
      INTEGER, Intent(IN) :: n_k
      INTEGER, Intent(IN) :: i_scale
      INTEGER, Intent(IN) :: i_scale_function
      LOGICAL, Intent(IN) :: l_continuum
      LOGICAL, Intent(IN) :: l_cont_gen
      LOGICAL, Intent(IN) :: l_self_broadening
      REAL  (RealK), Intent(IN) :: p
      REAL  (RealK), Intent(IN) :: t
      REAL  (RealK), Intent(IN), Dimension(:) :: amount
      REAL  (RealK), Intent(IN), Dimension(:) :: transmittance
      REAL  (RealK), Intent(IN), Dimension(:) :: trans_calc
      REAL  (RealK), Intent(IN), Dimension(:) :: k
      REAL  (RealK), Intent(IN), Dimension(:) :: w_k
      REAL  (RealK), Intent(IN), Dimension(:, :) :: scale_vector
!
    END SUBROUTINE write_fit_90
!
!
  END INTERFACE


! Check that either a lbl file has been provided or the code has access to
! HITRAN
  IF (l_fit_line_data.OR.l_fit_frn_continuum.OR.l_fit_self_continuum.OR. &
      (l_fit_cont_data.AND.l_cont_line_abs_weight)) THEN
    IF (.NOT.l_lbl_exist.AND..NOT.l_access_hitran.AND..NOT.l_access_xsc) THEN
      WRITE(cmessage,'(a)') &
        'Please provide either a HITRAN or line-by-line file.'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    ELSE IF (l_lbl_exist.AND.(l_access_hitran.OR.l_access_xsc)) THEN
      WRITE(*,'(a)') 'Line-by-line file exists and HITRAN file ' // &
          'provided. Using line-by-line file.'
    ENDIF
  ENDIF
  IF (l_fit_cont_data) THEN
    IF (l_cont_line_abs_weight) THEN
      IF (.NOT.l_lbl_exist.AND..NOT.l_access_hitran) THEN
        WRITE(cmessage,'(a)') &
          'Please provide either a HITRAN line database, ' // &
          'or line-by-line file.'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      ELSE IF (.NOT.l_access_cia.AND. &
               .NOT.l_use_h2o_frn_param.AND..NOT.l_use_h2o_self_param) THEN
        WRITE(cmessage,'(a)') 'Please provide a HITRAN CIA file or use ' // &
            'built-in water vapour continuum.'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    ELSE IF (.NOT.l_lbl_exist.AND..NOT.l_access_cia.AND. &
             .NOT.l_use_h2o_frn_param.AND..NOT.l_use_h2o_self_param) THEN
      WRITE(cmessage,'(a)') 'Please provide a HITRAN CIA file use ' // &
          'built-in water continuum.'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    ELSE IF (l_lbl_exist.AND.l_access_cia) THEN
      WRITE(*,'(a)') 'Line-by-line file exists and HITRAN file ' // &
          'provided. Using line-by-line file.'
    END IF
  END IF

! Check that scaling is look-up table if using generalised continuum formulation
  IF (l_fit_cont_data .AND. i_scale_function /= ip_scale_t_lookup .AND. &
      i_ck_fit /= ip_ck_none) THEN
    WRITE(cmessage,'(a)') 'A look-up table must be used with the ' // &
        'generalised continuum formulation.'
    ierr = i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)
  END IF

! Check that only a single set of temperatures has been provided
  IF (l_fit_cont_data .AND. n_p > 1) THEN
    WRITE(cmessage,'(a)') 'Only a single pressure should be specified ' // &
        'with the generalised continuum formulation.'
    ierr = i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)
  END IF

! Check that the scaling is a look-up table if self-broadening is included
  IF (l_fit_line_data .AND. l_self_broadening .AND. &
      i_scale_function /= IP_scale_lookup .AND. &
      i_ck_fit /= ip_ck_none) THEN
    WRITE(cmessage,'(a)') &
        'Scaling must be a look-up table when self-broadening ' // &
        'is included.'
    ierr = i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)
  ENDIF

! Perform preliminary allocation of arrays.
  IF (l_fit_line_data .OR. l_fit_cont_data) THEN
    ALLOCATE(u_l(2*nd_k_term+1, n_pt_pair))
    ALLOCATE(u_l_ref(2*nd_k_term+1))
    ln_p_calc=LOG(p_calc(1:n_pt_pair))
  ENDIF
  w_k=0.0_RealK
  k_opt_tmp=0.0_RealK
  kopt_all=0.0_RealK
  kopt_all_sb=0.0_RealK

  IF ( l_fit_self_continuum .OR. l_fit_frn_continuum ) THEN
!
!   Define the set of paths we will use for continuum data
    ALLOCATE(u_c(n_path_c))
    ALLOCATE(trans_c(n_path_c))
    DO i = 1, n_path_c
      u_c(i) = umin_c * EXP( REAL(i-1, RealK) * &
        LOG( umax_c / umin_c ) / REAL((n_path_c-1), RealK) )
    ENDDO
!
!   The arrays of points to be fitted are sized for the amounts
!   of absorber and the number of partial pressures.
    ALLOCATE(u_fit_c(n_path_c * n_pp, n_pt_pair))
    ALLOCATE(trans_fit_c(n_path_c * n_pp, n_pt_pair))
!
  ENDIF

! Find the median profile of P/Ts:
  i = 1
  n_t = 1
  DO ipt = 2, n_pt_pair
    IF (p_calc(ipt) == p_calc(ipt-1)) THEN
      n_t(i) = n_t(i)+1
    ELSE
      i=i+1
    END IF
  END DO
  mask_pt_ref = .FALSE.
  DO i = 1,n_p
    index_pt_ref(i) = SUM(n_t(1:i-1)) + 1 + n_t(i)/2
    mask_pt_ref(index_pt_ref(i)) = .TRUE.
  END DO

! Convert band limits from wavelength (m) to wavenumber (m-1)
  ALLOCATE(band_min(n_band))
  ALLOCATE(band_max(n_band))
  DO ib=1,n_band
!   Round to 8 significant figures due to accuracy of input format
    band_min(ib) = 1.0E+08 / &
      10.0_RealK**AINT(LOG10( 1.0_RealK / wavelength_long(ib) ))
    band_max(ib) = 1.0E+08 / &
      10.0_RealK**AINT(LOG10( 1.0_RealK / wavelength_short(ib)))
    band_min(ib) = ANINT(band_min(ib)/wavelength_long(ib))/band_min(ib)
    band_max(ib) = ANINT(band_max(ib)/wavelength_short(ib))/band_max(ib)
  END DO

! Initialise lbl file
  IF (l_lbl_exist) THEN
    CALL input_lbl_band_cdf_init
    IF (ierr /= i_normal) THEN
      cmessage='Error in input_lbl_band_cdf_init'
      CALL ereport(RoutineName, ierr, cmessage)
    END IF
    nu_band_adjust=nu_wgt_all(1)-nu_inc/2.0_RealK
  ELSE
    nu_inc=nu_inc_0
    nu_band_adjust=band_min(1)
  END IF

! Check that band limits agree with wavenumber grid in file. Note that the
! use of a tolerance of 1e-9 due to use of the format e16.9 for wavelengths
! in the spectral file
  DO ibb=1, n_selected_band
    ib=list_band(ibb)
    IF (abs(band_min(ib)-(nint((band_min(ib)-nu_band_adjust)/nu_inc)*nu_inc + &
      nu_band_adjust)) >= 1.0E-09_RealK*band_min(ib) .OR. &
      abs(band_max(ib)-(nint((band_max(ib)-nu_band_adjust)/nu_inc)*nu_inc + &
      nu_band_adjust)) >= 1.0E-09_RealK*band_max(ib)) THEN

!     If necessary adjust the band limits to enclose an integer number
!     of increments
      band_min(ib) = ANINT(1.0_RealK/(nu_inc*wavelength_long(ib)))*nu_inc
      band_max(ib) = ANINT(1.0_RealK/(nu_inc*wavelength_short(ib)))*nu_inc

      IF (l_lbl_exist) THEN
        ! Ensure externally supplied grid aligns with band increments
        grid_offset = MINVAL(nu_wgt_all - band_min(ib), &
          mask=(nu_wgt_all - band_min(ib) >= 1.0E-09_RealK))
        IF (grid_offset > nu_inc + 1.0E-09_RealK) THEN
          WRITE(cmessage,'(a,i0)') &
            'The line-by-line file does not contain sufficient data'//&
            ' for band ', ib
          ierr = i_err_fatal
          CALL ereport(RoutineName, ierr, cmessage)
        END IF
        band_min(ib) = band_min(ib) + grid_offset - nu_inc/2.0_RealK
        band_max(ib) = band_max(ib) + grid_offset - nu_inc/2.0_RealK
      END IF

      WRITE(*,'(a,i4,a,2f16.3,a)') 'Band ', ib, ' limits adjusted to:', &
        band_min(ib), band_max(ib), ' m-1'
    END IF
  END DO

! Calculate scaling of continuum column mass to gas column mass
  IF ((i_ck_fit /= ip_ck_none) .AND. &
      l_fit_cont_data .AND. l_cont_line_abs_weight) THEN
    l_wgt_scale_sqrt = i_gas_1 == i_gas_2
    IF (l_wgt_scale_sqrt) THEN
      u_wgt_scale = max_path_wgt/sqrt(max_path)
    ELSE
      u_wgt_scale = max_path_wgt/max_path
    END IF
  END IF

! Check if calculating continuum data is required (if provided access)
  l_calc_cont = .NOT. l_lbl_exist .OR. &
    (l_fit_cont_data .AND. l_cont_line_abs_weight)

! Initialise lbl output file
  IF (.NOT.l_lbl_exist) CALL output_lbl_band_cdf_init

! Initialise map file
  IF (l_load_map .OR. l_load_wgt) THEN
    CALL input_map_band_cdf_init
  ELSE IF (l_save_map) THEN
    CALL output_map_band_cdf_init
  END IF

  Bands: DO ibb=1, n_selected_band

!   The output file is opened afresh for each band to enable
!   recovery from crashes.
    IF (ibb == 1) THEN
      IF ((i_ck_fit /= ip_ck_none) .OR. &
          l_fit_self_continuum .OR. l_fit_frn_continuum) THEN
        OPEN(UNIT=iu_k_out, FILE=file_k, POSITION='REWIND')
      END IF
      OPEN(UNIT=iu_monitor, FILE=file_monitor, POSITION='REWIND')
    ELSE
      IF ((i_ck_fit /= ip_ck_none) .OR. &
          l_fit_self_continuum .OR. l_fit_frn_continuum) THEN
        OPEN(UNIT=iu_k_out, FILE=file_k, POSITION='APPEND')
      END IF
      OPEN(UNIT=iu_monitor, FILE=file_monitor, POSITION='APPEND')
    ENDIF
!
    ib=list_band(ibb)
!
    WRITE(*,"(a)") "==============================="
    WRITE(*,'(a, i5)') "Processing band ", ib
    WRITE(iu_monitor,"(a)") "==============================="
    WRITE(iu_monitor,'(a, i5)') "Fitting in band: ", ib
    CALL cpu_time(start_band)

!   Set the wavenumbers of the weighting points in the band (in nu_wgt).
    CALL set_wgt_int
    IF (ierr /= i_normal) THEN
      cmessage='Error in set_wgt_int'
      CALL ereport(RoutineName, ierr, cmessage)
    END IF

    ! Allocate sub-band arrays
    ALLOCATE(sub_band_k(n_nu))
    ALLOCATE(sub_band_w(n_nu))
    ALLOCATE(sub_band_nu_short(n_nu))
    ALLOCATE(sub_band_nu_long(n_nu))

!   Find P/T closest to the reference conditions.
    IF (l_scale_pT) THEN
      ALLOCATE(rdiff(1:n_pt_pair))
      rdiff = &
            ( (p_calc(1:n_pt_pair)-p_ref(ib))/p_ref(ib) )**2 + &
            ( (t_calc(1:n_pt_pair)-t_ref(ib))/t_ref(ib) )**2
      ipt_ref=MINLOC(rdiff,1)
      DEALLOCATE(rdiff)
    ELSE
!     Use first P/T pair for output if no scaling is done.
      ipt_ref = 1
    END IF

!   Processing depends on the purpose of the call.
    IF (l_fit_line_data .OR. l_fit_frn_continuum .OR. &
        l_fit_self_continuum .OR. l_fit_cont_data) THEN
!
      l_transparent_fit=.FALSE.
      IF (l_lbl_exist) THEN
        CALL cpu_time(timer1)
        CALL input_lbl_band_cdf ! Read lbl file for current band
        CALL cpu_time(timer2)
        WRITE(iu_monitor,*) 'CPU time reading LBL file: ',timer2-timer1
        IF (l_self_broadening) THEN
          l_transparent_fit=SUM(kabs_all_sb) < EPSILON(kabs_all_sb)
        ELSE
          l_transparent_fit=SUM(kabs_all) < EPSILON(kabs_all)
        END IF
      ELSE IF (l_access_hitran) THEN
        CALL access_hitran_int
        IF (ierr /= i_normal) THEN
          cmessage='Error in access_hitran_int'
          CALL ereport(RoutineName, ierr, cmessage)
        END IF
        hitran_lines(ibb)=num_lines_in_band
        l_transparent_fit=num_lines_in_band.EQ.0
      ELSE IF (l_access_xsc) THEN
        CALL access_xsc_int
        IF (ierr /= i_normal) THEN
          cmessage='Error in access_xsc_int'
          CALL ereport(RoutineName, ierr, cmessage)
        END IF
        hitran_lines(ibb)=num_lines_in_band
        l_transparent_fit=num_lines_in_band.EQ.0
      END IF
      IF (l_calc_cont .AND. l_access_cia) THEN
        CALL access_cia_int
        IF (ierr /= i_normal) THEN
          cmessage='Error in access_cia_int'
          CALL ereport(RoutineName, ierr, cmessage)
        END IF
        l_transparent_fit=num_cia_lines_in_band.EQ.0
      ENDIF

      IF (l_transparent_fit) THEN
!
!       A null transparent fit can be used.
        IF (l_fit_line_data .OR. l_fit_cont_data) CALL fit_transparent_int
        IF (.NOT.l_lbl_exist) kabs_all=0.0
        IF (ierr /= i_normal) THEN
          cmessage='Error in fit_transparent_int'
          CALL ereport(RoutineName, ierr, cmessage)
        END IF

        IF ( (l_fit_self_continuum .OR. l_fit_frn_continuum) .AND. &
             (i_weight /= ip_weight_planck) .AND. &
             (i_weight /= ip_weight_d_planck) ) THEN
          CALL rad_weight_90(i_weight, nu_wgt, SolarSpec, t_calc(1), wgt)
          CALL apply_response_int
          IF (ierr /= i_normal) THEN
            cmessage='Error in apply_response_int'
            CALL ereport(RoutineName, ierr, cmessage)
          END IF
        END IF
!       Continuum absorption must be fitted anyway.
        IF (l_fit_self_continuum) THEN
          DO ipt=1, n_pt_pair
            IF ( (i_weight == ip_weight_planck) .OR. &
                 (i_weight == ip_weight_d_planck) ) THEN
              CALL rad_weight_90(i_weight, nu_wgt, SolarSpec, t_calc(ipt), wgt)
              CALL apply_response_int
              IF (ierr /= i_normal) THEN
                cmessage='Error in apply_response_int'
                CALL ereport(RoutineName, ierr, cmessage)
              END IF
            END IF
            CALL calc_self_trans_int
            IF (ierr /= i_normal) THEN
              cmessage='Error in calc_self_trans_int'
              CALL ereport(RoutineName, ierr, cmessage)
            END IF
          ENDDO
        ENDIF
        IF (l_fit_frn_continuum) THEN
          DO ipt=1, n_pt_pair
            IF ( (i_weight == ip_weight_planck) .OR. &
                 (i_weight == ip_weight_d_planck) ) THEN
              CALL rad_weight_90(i_weight, nu_wgt, SolarSpec, t_calc(ipt), wgt)
              CALL apply_response_int
              IF (ierr /= i_normal) THEN
                cmessage='Error in apply_response_int'
                CALL ereport(RoutineName, ierr, cmessage)
              END IF
            END IF
            CALL calc_frn_trans_int
            IF (ierr /= i_normal) THEN
              cmessage='Error in calc_frn_trans_int'
              CALL ereport(RoutineName, ierr, cmessage)
            END IF
          ENDDO
        ENDIF

!       There will be no sub-bands in this case.
        n_sub_band=0

      ELSE
!
        IF (.NOT.l_lbl_exist.AND.l_access_hitran) THEN
          WRITE(*,"(a,f12.6,2x,f12.6)") &
            "Wavenumbers of min and max lines are: ", &
            hitran_data(1) % frequency, &
            hitran_data(num_lines_in_band) % frequency
          WRITE(iu_monitor,"(a,f12.6,2x,f12.6)") &
            "Wavenumbers of min and max lines are: ", &
            hitran_data(1) % frequency, &
            hitran_data(num_lines_in_band) % frequency
!         
!         Set aside space for the adjusted line parameters:
!         this will be reused each time T and p change.
          ALLOCATE(adj_line_parm(num_lines_in_band))

!         Loop over gas fractions
          DO igf=1, n_gas_frac

!           Loop over each pair of pressure and temperature
            P_and_T: DO ipt=1, n_pt_pair

!             Adjust the line parameters of each line for the particular
!             ambient conditions.
!$OMP PARALLEL DO PRIVATE(i) NUM_THREADS(n_omp_threads)
              DO i = 1, num_lines_in_band
                CALL adjust_path ( &
                  hitran_data(i) % mol_num, &
                  hitran_data(i) % iso_num, &
                  hitran_data(i) % frequency, &
                  hitran_data(i) % intensity, &
                  hitran_data(i) % air_broadhw, &
                  hitran_data(i) % lstate_energy, &
                  hitran_data(i) % air_broad_coeff, &
                  hitran_data(i) % self_broadhw, &
                  t_calc(ipt), p_calc(ipt), &
                  gas_frac(igf), &
                  adj_line_parm(i) % line_centre, &
                  adj_line_parm(i) % S_adj, &
                  adj_line_parm(i) % alpha_lorentz, &
                  adj_line_parm(i) % alpha_lorentz_air, &
                  adj_line_parm(i) % alpha_lorentz_self, &
                  adj_line_parm(i) % alpha_doppler)
              ENDDO
!$OMP END PARALLEL DO
!
!             If using the CKD continuum we now calculate the absorption
!             of each line at its cutoff: the check on i_gas is used for
!             safety.
              IF ( l_ckd_cutoff .AND. (i_gas == IP_H2O) ) THEN
                ALLOCATE(k_cutoff(num_lines_in_band))
!$OMP PARALLEL DO PRIVATE(i) NUM_THREADS(n_omp_threads)
                DO i = 1, num_lines_in_band
                  CALL voigt_profile ( &
                    adj_line_parm(i) % line_centre+line_cutoff,  &
                    adj_line_parm(i) % line_centre, &
                    adj_line_parm(i) % S_adj, &
                    adj_line_parm(i) % alpha_lorentz, &
                    adj_line_parm(i) % alpha_doppler, &
                    k_cutoff(i))
                ENDDO
!$OMP END PARALLEL DO
              ENDIF

!             Set line profile correction parameters for current condition
              CALL set_line_prof_corr_cnst(t_calc(ipt), i_line_prof_corr)

              CALL calc_line_abs_int
              IF (ierr /= i_normal) THEN
                cmessage='Error in calc_line_abs_int'
                CALL ereport(RoutineName, ierr, cmessage)
              END IF
              kabs_all(:,ipt)=kabs
              IF (l_self_broadening) &
                kabs_all_sb(:,ipt,igf)=kabs

              IF ( l_ckd_cutoff .AND. (i_gas == IP_H2O) ) THEN
                DEALLOCATE(k_cutoff)
              ENDIF

            END DO p_and_T

          END DO

          DEALLOCATE(adj_line_parm)

        ELSE IF (.NOT.l_lbl_exist.AND.l_access_xsc) THEN

!         Loop over each pair of pressure and temperature
          DO ipt=1, n_pt_pair
            CALL calc_xsc_abs_int
            kabs_all(:,ipt)=kabs
          END DO
!         Deallocate structure holding raw cross-section data
          DO i = 1, num_lines_in_band
            DEALLOCATE(xsc(i)%data)
          END DO
          DEALLOCATE(xsc)
        
        END IF

        IF (l_fit_cont_data .AND. l_cont_line_abs_weight) THEN
!         Save line absorption for use in weighting
          kabs_all_lines = kabs_all
        END IF

        IF (l_calc_cont .AND. l_access_cia) THEN
!         Loop over temperatures
          DO ipt=1, n_pt_pair
            CALL calc_cia_abs_int
            kabs_all(:,ipt)=kabs
          END DO
!         Deallocate structure holding raw cross-section data
          DO i = 1, num_cia_lines_in_band
            DEALLOCATE(cia(i)%wavenumber)
            DEALLOCATE(cia(i)%data)
          END DO
          DEALLOCATE(cia)

        ELSE IF (l_calc_cont .AND. l_use_h2o_frn_param) THEN

!         Loop over temperatures
!$OMP PARALLEL DO PRIVATE(ipt, kabs) NUM_THREADS(n_omp_threads)
          DO ipt=1, n_pt_pair
            CALL foreign_continuum(t_calc(ipt), 0.0_RealK, 0.0_RealK, &
              n_nu, nu_wgt, .TRUE., kabs)
            kabs_all(:,ipt)=kabs
          END DO
!$OMP END PARALLEL DO

        ELSE IF (l_calc_cont .AND. l_use_h2o_self_param) THEN

!         Loop over temperatures
!$OMP PARALLEL DO PRIVATE(ipt, kabs) NUM_THREADS(n_omp_threads)
          DO ipt=1, n_pt_pair
            CALL self_continuum(t_calc(ipt), 0.0_RealK, 0.0_RealK, &
              n_nu, nu_wgt, .TRUE., kabs)
            kabs_all(:,ipt)=kabs
          END DO
!$OMP END PARALLEL DO

        END IF

        IF (l_fit_self_continuum .OR. l_fit_frn_continuum .OR. &
            (i_ck_fit /= ip_ck_none)) THEN
          IF ( i_weight /= ip_weight_planck .AND. &
               i_weight /= ip_weight_d_planck ) THEN
            ! The weighting does not depend on temperature and can be
            ! calculated once for each wavenumber point in the band.
            CALL rad_weight_90(i_weight, nu_wgt, SolarSpec, t_calc(1), wgt)
            CALL apply_response_int
            wgt_sv=wgt
          END IF
        END IF

        IF (i_ck_fit /= ip_ck_none) THEN
!         Define the mapping for correlated-k.

!         If self-broadening is included, use middle gas fraction.
          IF (l_self_broadening) &
            kabs_all=kabs_all_sb(:,:,n_gas_frac/2+1)

          ALLOCATE(kabs_rank1(n_nu))
          ALLOCATE(kabs_rank2(n_nu))

          max_p_calc=MAXVAL(p_calc(1:n_pt_pair))
          DO i=1,n_nu
!           Try to group together absorption coefficients
!           that have a similar scaling behaviour.

!           Rank absorption coefficients by the pressure at which
!           they peak:
            kabs_rank1(i)=ln_p_calc(MAXLOC(kabs_all(i,:),1,mask_pt_ref))

!           Rank by the pressure where optical depth = 1.

            ipt=index_pt_ref(1)
            IF (l_fit_cont_data) THEN
!             The continuum absorption coefficient is pressure independent,
!             i.e. we can calculate the pressure at which the optical depth = 1
!             directly
              IF (kabs_all(i,ipt) < TINY(kabs_all(i,ipt))) THEN
                kabs_rank2(i) = -1.0_RealK
              ELSE
                kabs_rank2(i) = max_p_calc/SQRT(max_path*kabs_all(i,ipt))
              END IF

            ELSE
!             Calculate pressure at which the optical depth = 1
              kabs_rank2(i)=kabs_all(i,ipt)*p_calc(ipt)
              IF (kabs_rank2(i) > max_p_calc/max_path .OR. &
                  n_p == 1) THEN
                IF (kabs_all(i,ipt) < TINY(kabs_all(i,ipt))) THEN
                  kabs_rank2(i) = -1.0_RealK
                ELSE
                  kabs_rank2(i) = max_p_calc / (max_path*kabs_all(i,ipt))
                END IF
              ELSE
                DO ip=2,n_p
                  ipt=index_pt_ref(ip)
                  rank_work = kabs_rank2(i) + kabs_all(i,ipt) &
                    *(p_calc(ipt)-p_calc(index_pt_ref(ip-1)))
                  IF (rank_work > max_p_calc/max_path .OR. &
                      ip == n_p) THEN
                    IF (kabs_all(i,ipt) < TINY(kabs_all(i,ipt))) THEN
                      kabs_rank2(i) = -1.0_RealK
                    ELSE
                      kabs_rank2(i) = (max_p_calc/max_path &
                        - kabs_rank2(i)) / kabs_all(i,ipt) &
                        + p_calc(index_pt_ref(ip-1))
                    END IF
                    EXIT
                  ELSE
                    kabs_rank2(i) = rank_work
                  END IF
                END DO
              END IF
            END IF
          END DO

          IF (l_load_map) THEN
!           Use pre-defined mapping

!           Set absorption coefficients for this bin using an
!           effective value for optical depth = 1
            IF (l_fit_cont_data) THEN
              WHERE (kabs_rank2 < TINY(kabs_rank2))
                kabs = 0.0_RealK
              ELSEWHERE
                kabs = max_p_calc**2/(max_path*kabs_rank2**2)
              END WHERE
            ELSE
              WHERE (kabs_rank2 < TINY(kabs_rank2))
                kabs = 0.0_RealK
              ELSEWHERE
                kabs = max_p_calc/(max_path*kabs_rank2)
              END WHERE
            END IF
          
!           Read mapping, g-points and reference k-term weights
            CALL input_map_band_cdf

          ELSE
!           Derive new mapping

            IF (l_load_wgt) THEN
!             Use pre-defined pressure mapping and reference k-term weights
              CALL input_wgt_band_cdf

            ELSE
!             Derive new pressure mapping

!             Max/min in the pressure ranking
              pmax=MAXVAL(kabs_rank1)
              pmin=MINVAL(kabs_rank1)

              ALLOCATE(kabs_mask(n_nu))
              kabs_mask=0

              IF (i_ck_fit == ip_ck_bin) THEN
!               Split mapping into a number of bins in the pressure ranking
                n_pre=MAX(INT(pmax-pmin)/3,1)
              ELSE
                n_pre=1
              END IF
              DO
                pstep=(pmax-pmin)/REAL(n_pre, RealK)
                DO ipb=1,n_pre
                  WHERE (kabs_rank1 >= pmin+pstep*REAL(ipb-1, RealK))
                    kabs_mask=ipb
                  END WHERE
                END DO
!               If absorption coefficients are optically thin place in last bin
                WHERE (kabs_rank2 >= max_p_calc)
                  kabs_mask=n_pre
                END WHERE
                DO ipb=1,n_pre
                  pcount(ipb)=COUNT(kabs_mask==ipb)
                END DO
                IF (MINVAL(pcount(1:n_pre)) > 0) EXIT
                IF (n_pre == 1) EXIT
                n_pre=n_pre-1
              END DO

              pstart(1)=1
              DO ipb=1,n_pre
!               Find start and end points of current pressure bin
                pcount(ipb)=COUNT(kabs_mask==ipb)
                pend(ipb)=pstart(ipb)+pcount(ipb)-1
                IF (ipb < n_pre) pstart(ipb+1)=pend(ipb)+1
              END DO

              CALL map_heap_func(REAL(kabs_mask, RealK), pmap)

              DEALLOCATE(kabs_mask)

            END IF

!           Set absorption coefficients for this bin using an
!           effective value for optical depth = 1
            IF (l_fit_cont_data) THEN
              WHERE (kabs_rank2 < TINY(kabs_rank2))
                kabs = 0.0_RealK
              ELSEWHERE
                kabs = max_p_calc**2/(max_path*kabs_rank2**2)
              END WHERE
            ELSE
              WHERE (kabs_rank2 < TINY(kabs_rank2))
                kabs = 0.0_RealK
              ELSEWHERE
                kabs = max_p_calc/(max_path*kabs_rank2)
              END WHERE
            END IF

            DO ipb=1,n_pre
              WRITE(iu_monitor,*) 'pstart, pend: ', pstart(ipb),pend(ipb)

!             Reorder absorption coefficients for this bin
              CALL map_heap_func( kabs(pmap(pstart(ipb):pend(ipb))), &
                gmap(pstart(ipb):pend(ipb)) )
              gmap(pstart(ipb):pend(ipb)) = pstart(ipb)-1 + &
                gmap(pstart(ipb):pend(ipb))
            END DO

!           Define the global mapping from the original wavenumber
!           spectrum to the binned and reordered g-space.
            map = pmap(gmap)

          END IF

          DEALLOCATE(kabs_rank1)

          IF ( (i_weight == ip_weight_planck) .OR. &
               (i_weight == ip_weight_d_planck) ) THEN
            ! Set the weights to use the fraction of the Plankian
            ! at the reference temperature where optical depth = 1
            ! (inefficiently coded)
            ALLOCATE(wgt_sum(n_p))
            DO ipt=1,n_p
              CALL rad_weight_90(i_weight, nu_wgt(:), &
                                SolarSpec, t_calc(index_pt_ref(ipt)), wgt)
              wgt_sum(ipt)=SUM(wgt)
            END DO
            DO i=1,n_nu
              IF (kabs_rank2(i) < TINY(kabs_rank2(i))) THEN
                ipoint = n_p
              ELSE
                ipoint=MINLOC( ABS( &
                  ln_p_calc(index_pt_ref) - LOG(kabs_rank2(i)) ), 1 )
              END IF
              CALL rad_weight_90(i_weight, nu_wgt(i:i), SolarSpec, &
                t_calc(index_pt_ref(ipoint)), wgt(i:i))
              wgt(i) = wgt(i)/wgt_sum(ipoint)
            END DO
            DEALLOCATE(wgt_sum)
            CALL apply_response_int
          ELSE
            wgt=wgt_sv
          END IF

          DEALLOCATE(kabs_rank2)

!         Save the unsorted weights to calculate the sub-band weights
          ALLOCATE(wgt_unsorted(n_nu))
          wgt_unsorted=wgt
!         Adjust the boundary weights using the correct bounds
          wgt_unsorted(1)=wgt_unsorted(1) &
            *(nu_wgt(1)+nu_inc/2.0_RealK-1.0/wavelength_long(ib))/nu_inc
          wgt_unsorted(n_nu)=wgt_unsorted(n_nu) &
            *(1.0/wavelength_short(ib)-nu_wgt(n_nu)+nu_inc/2.0_RealK)/nu_inc

!         Sort the weighting function.
          kabs=kabs(map)
          wgt=wgt(map)
          IF (l_fit_cont_data .AND. l_cont_line_abs_weight) &
            kabs_lines=kabs_all_lines(map,index_pt_ref(1))

!         Don't consider path lengths greater than max_path*10
!         for calculating the error on the k-term fit:
          umax_kopt=max_path*10.0_RealK
          umin_kopt=umin_kopt_default

          IF (l_load_map .OR. l_load_wgt) THEN
!           The g-points are given in the mapping file, use them to calculate
!           k_opt
            CALL calc_k_opt_ref

          ELSE
!           Set the g-points within each pressure bin
            n_k_last=0
            ig(0)=0
            integ_wgt=nu_inc * SUM(wgt(1:n_nu))
            DO ipb=1,n_pre
!             Integrate the sorted weightings across the band.
              integ_wgt_tmp=nu_inc * SUM(wgt(pstart(ipb):pend(ipb)))

!             Set the number of g-points based on the error in the
!             transmission.
              CALL set_g_point_90(pcount(ipb), nu_inc, &
                kabs(pstart(ipb):pend(ipb)), &
                wgt(pstart(ipb):pend(ipb)), integ_wgt, &
                i_ck_fit, tol*REAL(ipb,RealK), max_path, &
                l_fit_cont_data .AND. l_cont_line_abs_weight, &
                kabs_lines(pstart(ipb):pend(ipb)), &
                l_wgt_scale_sqrt, u_wgt_scale, &
                nd_k_term, iu_monitor, &
                n_k(ib), w_k(n_k_last+1:nd_k_term,ib), &
                k_opt(n_k_last+1:nd_k_term,ib), &
                k_ave(n_k_last+1:nd_k_term,ib), &
                ig_tmp(0:nd_k_term), ierr)

                ig(n_k_last+1:n_k_last+n_k(ib))= &
                  ig_tmp(1:n_k(ib))+ig(n_k_last)
                n_k_last=n_k_last+n_k(ib)
            END DO
            n_k(ib)=n_k_last
          END IF
          WRITE(*,*) 'Number of k-terms in band: ',n_k(ib)
          WRITE(iu_monitor,*) 'Number of k-terms in band: ',n_k(ib)
          CALL ck_trans_fit

!         Save mapping if requested
          IF (l_save_map) CALL output_map_band_cdf

          ! Calculate sub-band mapping
          ALLOCATE( kmap(n_nu) )
          ! kmap holds the k-term number for each nu
          DO ik=1, n_k(ib)
            kmap(map( ig(ik-1)+1 : ig(ik) )) = ik
          END DO
          ! Sub-bands are found where k-terms are mapped to contiguous
          ! regions in frequency
          k=0
          n_sub_band=0
          DO i=1, n_nu
            IF (kmap(i) /= k) THEN
              k=kmap(i)
              n_sub_band=n_sub_band+1
              sub_band_k(n_sub_band)=k
              sub_band_w(n_sub_band)=wgt_unsorted(i)
              sub_band_nu_short(n_sub_band)=nu_wgt(i)-nu_inc/2.0_RealK
            ELSE
              sub_band_w(n_sub_band)=sub_band_w(n_sub_band)+wgt_unsorted(i)
            END IF
            sub_band_nu_long(n_sub_band)=nu_wgt(i)+nu_inc/2.0_RealK
          END DO
          sub_band_w_sum=SUM(sub_band_w(1:n_sub_band))
          DO isb=1, n_sub_band
            sub_band_w(isb)=sub_band_w(isb)/sub_band_w_sum
          END DO
          DEALLOCATE(kmap)
          DEALLOCATE(wgt_unsorted)

        END IF

        IF (l_fit_self_continuum .OR. l_fit_frn_continuum .OR. &
            (i_ck_fit /= ip_ck_none)) THEN

          DO ipt=1, n_pt_pair
!           Calculate the weighting function across the band.
            IF ( i_weight == ip_weight_planck .OR. &
                 i_weight == ip_weight_d_planck ) THEN
              CALL rad_weight_90(i_weight, nu_wgt, SolarSpec, t_calc(ipt), wgt)
              CALL apply_response_int
            ELSE
              wgt=wgt_sv
            END IF
          
            IF (l_fit_cont_data .AND. l_cont_line_abs_weight) &
              kabs_lines=kabs_all_lines(1:n_nu,ipt)
          
            IF (ipt == ipt_ref) THEN
              wgt_ref(1:n_nu) = nu_inc * wgt(1:n_nu)
            END IF
          
!           Perform the appropriate fits.
            IF (l_fit_self_continuum) THEN
              kabs=kabs_all(1:n_nu,ipt)
              integ_wgt=nu_inc * SUM(wgt(1:n_nu))
              CALL calc_self_trans_int
              IF (ierr /= i_normal) THEN
                cmessage='Error in calc_self_trans_int'
                CALL ereport(RoutineName, ierr, cmessage)
              END IF
            ELSE IF (l_fit_frn_continuum) THEN
              kabs=kabs_all(1:n_nu,ipt)
              integ_wgt=nu_inc * SUM(wgt(1:n_nu))
              CALL calc_frn_trans_int
              IF (ierr /= i_normal) THEN
                cmessage='Error in calc_frn_trans_int'
                CALL ereport(RoutineName, ierr, cmessage)
              END IF
            ELSE IF (l_fit_line_data .OR. l_fit_cont_data) THEN
              IF (l_self_broadening) THEN
                wgt_sv=wgt
                DO igf=1, n_gas_frac
                  wgt=wgt_sv
                  kabs=kabs_all_sb(:,ipt,igf)
                  CALL ck_fit_k
                  IF (ierr /= i_normal) THEN
                    write(cmessage,'(a,2i3)') &
                      'Error in ck_fit_k,: ipt, igf =', ipt, igf
                    CALL ereport(RoutineName, ierr, cmessage)
                  END IF
                  kopt_all_sb(:,ipt,igf,ibb)=kopt_all(:,ipt,ibb)
                END DO
              ELSE
                kabs=kabs_all(1:n_nu,ipt)
                CALL ck_fit_k
                IF (ierr /= i_normal) THEN
                  write(cmessage,'(a,i3)') 'Error in ck_fit_k: ipt =', ipt
                  CALL ereport(RoutineName, ierr, cmessage)
                END IF
              END IF
            END IF
          ENDDO
          
          IF (l_fit_frn_continuum) k_cont => k_opt_frn(ib)
          IF (l_fit_self_continuum) k_cont => k_opt_self(ib)
          
          IF (l_scale_pT) THEN
!         
            IF (l_fit_line_data) THEN
              SELECT CASE(i_scale_function)
              CASE (IP_scale_power_law, IP_scale_power_quad,  &
                    IP_scale_doppler_quad )
                CALL fit_scale_line_int
                IF (ierr /= i_normal) THEN
                  cmessage='Error after fit_scale_line_int'
                  CALL ereport(RoutineName, ierr, cmessage)
                END IF
              CASE (IP_scale_dbl_pow_law, IP_scale_dbl_pow_quad,  &
                    IP_scale_dbl_dop_quad )
                CALL fit_scale_line_int2
                IF (ierr /= i_normal) THEN
                  cmessage='Error after fit_scale_line_int2'
                  CALL ereport(RoutineName, ierr, cmessage)
                END IF
              END SELECT
            END IF
            IF (l_fit_frn_continuum .OR. l_fit_self_continuum) THEN
              CALL fit_scale_cont_int
              IF (ierr /= i_normal) THEN
                cmessage='Error in fit_scale_cont_int'
                CALL ereport(RoutineName, ierr, cmessage)
              END IF
            END IF
          END IF
          
          IF (l_fit_line_data .AND. l_scale_pT) DEALLOCATE(trans_pt_k)
        END IF

      ENDIF ! End else block for transparent fit
!
    ENDIF
!
    DEALLOCATE(wgt)
    DEALLOCATE(wgt_sv)
    DEALLOCATE(kabs)
    DEALLOCATE(map)
    DEALLOCATE(gmap)
    DEALLOCATE(pmap)
    DEALLOCATE(kabs_lines)
!
!   Write out the calculated fit.
    IF (i_ck_fit /= ip_ck_none) THEN

      CALL write_fit_90(iu_k_out, .FALSE., l_fit_cont_data, &
        l_self_broadening, ib, i_gas, i_index, i_index_1, i_index_2, &
        p_calc(ipt_ref), t_calc(ipt_ref), &
        n_path, u_l_ref, trans_ref, trans_calc, &
        n_k(ib), k_opt(:, ib), w_k(:, ib), IP_scale_term, &
        i_scale_function, scale_vector(:, :, ib))

      IF (i_scale_function == IP_scale_lookup) THEN
!       Write look-up table to file:
        ALLOCATE(p_lookup( n_p ))
        ALLOCATE(t_lookup( MAXVAL(n_t(1:n_p)), n_p ))
        ALLOCATE(k_lookup( MAXVAL(n_t(1:n_p)), n_p, n_gas_frac, n_k(ib) ))
        DO ik=1, n_k(ib)
          ipt=0
          DO ip=1, n_p
            p_lookup(ip)=p_calc(index_pt_ref(ip))
            DO it=1, n_t(ip)
              ipt=ipt+1
              t_lookup(it,ip)=t_calc(ipt)
              IF (l_self_broadening) THEN
                DO igf=1, n_gas_frac
                  k_lookup(it,ip,igf,ik)=kopt_all_sb(ik,ipt,igf,ibb)
                END DO
              ELSE
                k_lookup(it,ip,1,ik)=kopt_all(ik,ipt,ibb)
              END IF
            END DO
          END DO
        END DO
        WRITE(iu_k_out,'(2(a,i4),a)') 'Lookup table: ', &
          n_p, ' pressures, ', MAXVAL(n_t(:)), ' temperatures.'
        DO ip=1, n_p
          WRITE(iu_k_out,'(6(1PE13.6))') p_lookup(ip), &
            t_lookup(1:n_t(ip),ip)
        END DO
        IF (l_self_broadening) THEN
          WRITE(iu_k_out,'(/,a,i4,a)') 'Lookup table: ', &
            n_gas_frac, ' gas fractions.'
          WRITE(iu_k_out,'(6(1PE13.6))') gas_frac(1:n_gas_frac)
        END IF
        WRITE(iu_k_out,'(/,3(a,i4))') 'Band: ',ib,', gas: ',i_gas, &
          ', k-terms: ',n_k(ib)
        DO ik=1, n_k(ib)
          DO igf=1, n_gas_frac
            DO ip=1, n_p
              WRITE(iu_k_out,'(6(1PE13.6))') k_lookup(1:n_t(ip),ip,igf,ik)
            END DO
          END DO
        END DO
        DEALLOCATE(k_lookup)
        DEALLOCATE(t_lookup)
        DEALLOCATE(p_lookup)
      ELSE IF (i_scale_function == IP_scale_t_lookup) THEN
!       Write temperature look-up table to file (n_p = 1):
        WRITE(iu_k_out,'(a,i4,a)') 'Lookup table: ', &
          n_t(1), ' temperatures.'
        WRITE(iu_k_out,'(6(1PE13.6))') t_calc(1:n_t(1))
        WRITE(iu_k_out,'(/,4(a,i4))') 'Band: ', ib, &
          ', gases: ', i_gas_1, ' and ', i_gas_2, &
          ', k-terms: ', n_k(ib)
        DO ik=1, n_k(ib)
          WRITE(iu_k_out,'(6(1PE13.6))') kopt_all(ik,1:n_t(1),ibb)
        END DO
      END IF

      IF (n_sub_band > 1) THEN
        ! Write out sub-band mappings ordered in increasing wavelength
        WRITE(iu_k_out, '(/a,i6,a)') 'Sub-band mapping: ', &
          n_sub_band, ' sub-bands.'
        WRITE(iu_k_out, '(a)') &
        'Sub-band  k-term       weight       wavelength_short   wavelength_long'
        WRITE(iu_k_out, '(2i8, 3(2x,1PE16.9))') 1, &
          sub_band_k(n_sub_band), sub_band_w(n_sub_band), &
          wavelength_short(ib), 1.0_RealK/sub_band_nu_short(n_sub_band)
        DO isb=n_sub_band-1, 2, -1
          WRITE(iu_k_out, '(2i8, 3(2x,1PE16.9))') n_sub_band-isb+1, &
            sub_band_k(isb), sub_band_w(isb), &
            1.0_RealK/sub_band_nu_long(isb), 1.0_RealK/sub_band_nu_short(isb)
        END DO
        WRITE(iu_k_out, '(2i8, 3(2x,1PE16.9))') n_sub_band, &
          sub_band_k(1), sub_band_w(1), &
          1.0_RealK/sub_band_nu_long(1), wavelength_long(ib)
      END IF

    ENDIF

    IF (l_fit_frn_continuum) i_index_c = IP_frn_continuum
    IF (l_fit_self_continuum) i_index_c = IP_self_continuum
    IF (l_fit_frn_continuum .OR. l_fit_self_continuum) THEN
      ALLOCATE(trans_app_c(n_path_c*n_pp))
      trans_app_c(:) = EXP( -k_cont * u_fit_c(:, ipt_ref) )
      CALL write_fit_90(iu_k_out, .TRUE., .FALSE., &
        .FALSE., ib, i_gas, i_index_c, 0, 0, &
        p_calc(ipt_ref), t_calc(ipt_ref), &
        n_path_c*n_pp, u_fit_c(:, ipt_ref), &
        trans_fit_c(:, ipt_ref), trans_app_c(:), &
        1, (/ k_cont /), (/ 1.0_RealK /), IP_scale_term, &
        i_scale_function, scale_cont(:, :, ib))
      DEALLOCATE(trans_app_c)
    ENDIF
    
    IF (l_fit_cont_data .AND. l_cont_line_abs_weight) &
      kabs_all=kabs_all_lines
    IF (.NOT.l_lbl_exist) CALL output_lbl_band_cdf
    DEALLOCATE(wgt_ref)
    DEALLOCATE(nu_wgt)
    DEALLOCATE(kabs_all)
    IF (l_fit_cont_data .AND. l_cont_line_abs_weight) DEALLOCATE(kabs_all_lines)
    IF (l_self_broadening) DEALLOCATE(kabs_all_sb)
!
    IF ((i_ck_fit /= ip_ck_none) .OR. &
        l_fit_self_continuum .OR. l_fit_frn_continuum) THEN
      CLOSE(iu_k_out)
    END IF
    CALL cpu_time(finish_band)
    WRITE(iu_monitor,*) 'CPU time for band: ', finish_band-start_band
    CLOSE(iu_monitor)
!
    IF (l_access_HITRAN .AND. .NOT.l_lbl_exist) THEN
      DEALLOCATE(hitran_data, STAT = alloc_status)
      IF (alloc_status /= 0) THEN
        WRITE(*,"(a)") "Error deallocating array for HITRAN line data"
        EXIT
      ENDIF
    ENDIF

    IF (i_ck_fit /= ip_ck_none) DEALLOCATE(trans_ref)
    IF (i_ck_fit /= ip_ck_none) DEALLOCATE(trans_calc)

    DEALLOCATE(sub_band_nu_long)
    DEALLOCATE(sub_band_nu_short)
    DEALLOCATE(sub_band_w)
    DEALLOCATE(sub_band_k)

  ENDDO Bands

  IF (i_ck_fit /= ip_ck_none) CALL output_ck_cdf
  IF (ierr /= i_normal) THEN
    cmessage='Error after output_ck_cdf'
    CALL ereport(RoutineName, ierr, cmessage)
  END IF
  IF (l_lbl_exist) DEALLOCATE(nu_wgt_all)
  DEALLOCATE(band_min)
  DEALLOCATE(band_max)
  
  CALL close_lbl_files
  CALL close_map_files


CONTAINS


  SUBROUTINE access_hitran_int
!
!
    IMPLICIT NONE

    INTEGER, SAVE :: number_lines = 0
    INTEGER :: io_status    ! Error code for file I/O

    alloc_status = 0
!
!   Convert our identifier for the gas to that used by HITRAN.
    h_gas=hitran_number(i_gas)
    h_isotopes=hitran_isotopes(:, i_gas)
    IF (h_gas <= 0) THEN
      WRITE(iu_err, "(/a)") &
        "*** Error: This gas does not exist in the database."
      ierr=i_err_fatal
      RETURN
    ENDIF

!   Determine total number of lines in HITRAN file
    IF (ibb == 1) THEN
      DO
        READ(iu_lbl, *, IOSTAT = io_status)
        IF (io_status < 0) EXIT
        number_lines = number_lines + 1
      END DO
    END IF

    ALLOCATE  ( &
         hitran_data(number_lines), &
         STAT = alloc_status)
!
    IF (alloc_status /= 0) THEN
      WRITE(*,"(a)")"Error allocating array for HITRAN line data."
      ierr=i_err_fatal
      RETURN
    END IF
!
!   Extend the region of reading to capture lines outside the band
!   which will contribute.
    lower_band_limit = (band_min(ib) - line_cutoff) * 0.01_RealK
    upper_band_limit = (band_max(ib) + line_cutoff) * 0.01_RealK

!   To cover the whole range, the database must be rewound.
    REWIND(iu_lbl)
    CALL read_hitran ( &
          iu_lbl,            &
          iu_monitor,        &
          h_gas, h_isotopes, &
          number_lines,      &
          lower_band_limit,  &
          upper_band_limit,  &
          num_lines_in_band, &
          hitran_data )
!
    WRITE(*,"(a,i6)") "Number of HITRAN lines in band: ", &
      num_lines_in_band
    WRITE(iu_monitor,"(a,i6)") "Number of HITRAN lines in band: ", &
      num_lines_in_band
!
!
!
  END SUBROUTINE access_hitran_int



  SUBROUTINE access_xsc_int

    IMPLICIT NONE

    TYPE(StrXscHead) :: single_header
    REAL (RealK), ALLOCATABLE :: xsc_data(:)
    INTEGER :: io_status    ! Error code for file I/O

    lower_band_limit = band_min(ib) * 0.01_RealK
    upper_band_limit = band_max(ib) * 0.01_RealK

!   To cover the whole range, the database must be rewound for bands
!   after the first.
    IF (ibb > 1) REWIND(iu_lbl)

!   First find the number of records that contain cross-sections for
!   this band
    num_lines_in_band = 0
    DO
      READ(iu_lbl, xsc_header_format, IOSTAT = io_status) single_header
      IF (io_status < 0) EXIT ! EOF
      IF (single_header%wavenumber_min <= upper_band_limit .AND. &
          single_header%wavenumber_max >= lower_band_limit) THEN
        num_lines_in_band = num_lines_in_band + 1
      END IF
      ALLOCATE(xsc_data(single_header%no_pts))
      READ(iu_lbl, xsc_data_format, IOSTAT = io_status) xsc_data
      DEALLOCATE(xsc_data)
    END DO

!   Now allocate structure and read in required data
    IF (num_lines_in_band > 0) THEN
      ALLOCATE(xsc(num_lines_in_band))
      REWIND(iu_lbl)
      num_lines_in_band = 0
      DO
        READ(iu_lbl, xsc_header_format, IOSTAT = io_status) single_header
        IF (io_status < 0) EXIT ! EOF
        ALLOCATE(xsc_data(single_header%no_pts))
        READ(iu_lbl, xsc_data_format, IOSTAT = io_status) xsc_data
        IF (single_header%wavenumber_min <= upper_band_limit .AND. &
            single_header%wavenumber_max >= lower_band_limit) THEN
          num_lines_in_band = num_lines_in_band + 1
          xsc(num_lines_in_band)%head = single_header
          ALLOCATE(xsc(num_lines_in_band)%data(single_header%no_pts))
          WHERE (xsc_data < 0.0_RealK)
            xsc(num_lines_in_band)%data = 0.0_RealK
          ELSEWHERE
            xsc(num_lines_in_band)%data = xsc_data
          END WHERE
        END IF
        DEALLOCATE(xsc_data)
      END DO
    END IF

  END SUBROUTINE access_xsc_int



  SUBROUTINE access_cia_int

    TYPE(StrCIAHead) :: single_header
    REAL (RealK), ALLOCATABLE :: cia_wavenumber(:)
    REAL (RealK), ALLOCATABLE :: cia_data(:)
    INTEGER :: io_status    ! Error code for file I/O

    lower_band_limit = band_min(ib) * 0.01_RealK
    upper_band_limit = band_max(ib) * 0.01_RealK

!   To cover the whole range, the database must be rewound for bands
!   after the first.
    IF (ibb > 1) REWIND(iu_cia)

!   First find the number of records that contain cross-sections for
!   this band
    num_cia_lines_in_band = 0
    DO
      READ(iu_cia, cia_header_frmt, IOSTAT = io_status) single_header
      IF (io_status < 0) EXIT ! EOF
      IF (single_header%wavenumber_min <= upper_band_limit .AND. &
          single_header%wavenumber_max >= lower_band_limit) THEN
        num_cia_lines_in_band = num_cia_lines_in_band + 1
      END IF
      ALLOCATE(cia_wavenumber(single_header%no_pts))
      ALLOCATE(cia_data(single_header%no_pts))
      READ(iu_cia, cia_data_frmt, IOSTAT = io_status) &
        (cia_wavenumber(i), cia_data(i), i = 1, single_header%no_pts)
      DEALLOCATE(cia_wavenumber)
      DEALLOCATE(cia_data)
    END DO

!   Now allocate structure and read in required data
    IF (num_cia_lines_in_band > 0) THEN
      ALLOCATE(cia(num_cia_lines_in_band))
      REWIND(iu_cia)
      num_cia_lines_in_band = 0
      DO
        READ(iu_cia, cia_header_frmt, IOSTAT = io_status) single_header
        IF (io_status < 0) EXIT ! EOF
        ALLOCATE(cia_wavenumber(single_header%no_pts))
        ALLOCATE(cia_data(single_header%no_pts))
        READ(iu_cia, cia_data_frmt, IOSTAT = io_status)  &
          (cia_wavenumber(i), cia_data(i), i = 1, single_header%no_pts)
        IF (single_header%wavenumber_min <= upper_band_limit .AND. &
            single_header%wavenumber_max >= lower_band_limit) THEN
          num_cia_lines_in_band = num_cia_lines_in_band + 1
          cia(num_cia_lines_in_band)%head = single_header
          ALLOCATE(cia(num_cia_lines_in_band)%wavenumber(single_header%no_pts))
          ALLOCATE(cia(num_cia_lines_in_band)%data(single_header%no_pts))
          cia(num_cia_lines_in_band)%wavenumber = cia_wavenumber
          WHERE (cia_data < 0.0_RealK)
            cia(num_cia_lines_in_band)%data = 0.0_RealK
          ELSEWHERE
            cia(num_cia_lines_in_band)%data = cia_data
          END WHERE
        END IF
        DEALLOCATE(cia_wavenumber)
        DEALLOCATE(cia_data)
      END DO
    END IF

  END SUBROUTINE access_cia_int



  SUBROUTINE fit_transparent_int
!
!
    IMPLICIT NONE 
!
    n_k(ib)=1
    w_k(1, ib)=1.0_RealK
    k_opt(1, ib)=0.0_RealK
!
    SELECT CASE(i_scale_function)
      CASE (IP_scale_power_law)
        scale_vector(1:2, 1, ib) = &
          (/ 1.0_RealK, 0.5_RealK /)
      CASE (IP_scale_power_quad)
        scale_vector(1:3, 1, ib) = &
          (/ 1.0_RealK, 0.0_RealK, 0.0_RealK /)
      CASE (IP_scale_doppler_quad)
        scale_vector(1:4, 1, ib) = &
          (/ 1.0_RealK, -2.0_RealK, 0.0_RealK, 0.0_RealK /)
      CASE (IP_scale_dbl_pow_law)
        scale_vector(1:6, 1, ib) = &
          (/ 1.0_RealK, 0.5_RealK, &
             1.0_RealK, 0.5_RealK, &
             p_ref(ib), t_ref(ib) /)
      CASE (IP_scale_dbl_pow_quad)
        scale_vector(1:8, 1, ib) = &
          (/ 1.0_RealK, 0.0_RealK, 0.0_RealK, &
             1.0_RealK, 0.0_RealK, 0.0_RealK, &
             p_ref(ib), t_ref(ib) /)
      CASE (IP_scale_dbl_dop_quad)
        scale_vector(1:10, 1, ib) = &
          (/ 1.0_RealK, -2.0_RealK, 0.0_RealK, 0.0_RealK, &
             1.0_RealK, -2.0_RealK, 0.0_RealK, 0.0_RealK, &
             p_ref(ib), t_ref(ib) /)
    END SELECT
!
!   Transmission data of at least the most basic kind are required
!   for the output file.
    n_path     = 1
    u_l_ref(1)      = 0.0_RealK
    ALLOCATE(trans_ref(n_path))
    trans_ref  = (/ 1.0_RealK /)
    ALLOCATE(trans_calc(n_path))
    trans_calc = (/ 1.0_RealK /)
!
!
!
  END SUBROUTINE fit_transparent_int
!
!
!
  SUBROUTINE set_wgt_int

    IMPLICIT NONE 

!   Set up the array of wavenumbers at the weighting points within the
!   band. The band should contain a whole number of intervals (band limits
!   have been adjusted in the calling routine if required). Care is taken
!   to ensure excluded bands are properly dealt with.
    
!   Calculate band width
    band_width = band_max(ib) - band_min(ib)
    DO jx = 1, n_band_exclude(ib)
      band_width = band_width + &
        band_min(index_exclude(jx, ib)) - band_max(index_exclude(jx, ib))
    ENDDO

!   Calculate the number of wavenumber intervals in the band and allocate
!   the array to hold the wavenumbers at the mid-point of each interval.
    n_nu=NINT(band_width/nu_inc)
    ALLOCATE(nu_wgt(n_nu))
    i  = 1
    jx = 1
    nu_wgt(1)=band_min(ib)+nu_inc/2.0_RealK

    IF (n_band_exclude(ib) == 0) THEN
      nu_end = band_max(ib)
    ELSE
      nu_end   = band_min(index_exclude(jx, ib))
      nu_start = band_max(index_exclude(jx, ib))
    ENDIF
    DO
      i=i+1
      nu_wgt(i)=NINT((nu_wgt(i-1)+nu_inc-nu_wgt(1))/nu_inc)*nu_inc+nu_wgt(1)
      IF (i == n_nu) EXIT
      IF (nu_wgt(i) > nu_end) THEN
        IF (n_band_exclude(ib) == 0) THEN
!         We should have exited anyway, but this covers rounding errors.
          EXIT
        ELSE
          nu_wgt(i)=nu_start+nu_inc/2.0_RealK
          IF (jx < n_band_exclude(ib)) THEN
            jx = jx+1
            nu_end   = band_min(index_exclude(jx, ib))
            nu_start = band_max(index_exclude(jx, ib))
          ELSE
            nu_end = band_max(ib)
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    ALLOCATE(kabs_all(n_nu,n_pt_pair))
    IF (l_self_broadening) ALLOCATE(kabs_all_sb(n_nu,n_pt_pair,n_gas_frac))
    ALLOCATE(wgt_ref(n_nu))
    ALLOCATE(wgt(n_nu))
    ALLOCATE(wgt_sv(n_nu))
    ALLOCATE(kabs(n_nu))
    ALLOCATE(map(n_nu))
    ALLOCATE(gmap(n_nu))
    ALLOCATE(pmap(n_nu))
    ALLOCATE(kabs_lines(n_nu))
    IF (l_fit_cont_data .AND. l_cont_line_abs_weight) &
      ALLOCATE(kabs_all_lines(n_nu,n_pt_pair))

  END SUBROUTINE set_wgt_int
!
!
!
  SUBROUTINE apply_response_int
    USE spline_evaluate_mod, ONLY: spline_evaluate
!
!
    IMPLICIT NONE 
!
    IF (include_instrument_response) THEN
      DO i=1, n_nu
        CALL spline_evaluate(ierr, filter%n_pts, &
          filter%wavenumber, filter%response, filter%d2_response, &
          nu_wgt(i), response_0)
        IF (ierr == i_err_range) THEN
!         The filter function is taken to be 0 outside 
!         the explicit range. We therefore zero the response and
!         recover from the error.
          response_0 = 0.0
          ierr = i_normal
        ENDIF
        IF (ierr /= i_normal) RETURN
        wgt(i) = wgt(i) * response_0
      ENDDO
    ENDIF
!
!
!
  END SUBROUTINE apply_response_int
!
!
!
  SUBROUTINE calc_line_abs_int
!
    IMPLICIT NONE
!
!   Local variables.
    INTEGER :: jx
!     Loop index for excluded bands
    INTEGER :: j_nu_first, j_nu_last
!     First and last nu index at which to include line.
    INTEGER :: n_nu_excluded_upper, n_nu_excluded_lower
!     Number of excluded wavenumbers left of upper and lower wavenumber
!
!   Calculate the monochromatic absorption coefficients at each
!   wavenumber.
    IF (l_self_broadening) THEN
      WRITE(*, '(A, 1PE10.3, A, 1PE10.3, A, 0PF4.2, A)') &
        'Calculation of absorption coefficients at ', p_calc(ipt), &
        ' Pa, ', t_calc(ipt), ' K and ',  gas_frac(igf), &
        ' gas fraction.'
      WRITE(iu_monitor, '(A, 1PE10.3, A, 1PE10.3, A, 0PF4.2, A)') &
        'Calculation of absorption coefficients at ', p_calc(ipt), &
        ' Pa, ', t_calc(ipt), ' K and ', gas_frac(igf), &
        ' gas fraction.'
    ELSE
      WRITE(*, '(A, 1PE10.3, A8, 1PE10.3, A2)') &
        'Calculation of absorption coefficients at ', p_calc(ipt), &
        ' Pa and ', t_calc(ipt), ' K.'
      WRITE(iu_monitor, '(A, 1PE10.3, A8, 1PE10.3, A2)') &
        'Calculation of absorption coefficients at ', p_calc(ipt), &
        ' Pa and ', t_calc(ipt), ' K.'
    END IF
    kabs(1:n_nu)=0.0_RealK
!
    DO i = 1, num_lines_in_band
!
!     Find the range of wavenumbers the current line contributes to.
      upper_cutoff  = adj_line_parm(i) % line_centre + line_cutoff
      lower_cutoff  = adj_line_parm(i) % line_centre - line_cutoff
!
!     Find the index of the first and last wavenumber at which to include line.
!
      n_nu_excluded_upper = 0
      n_nu_excluded_lower = 0
      DO jx = 1, n_band_exclude(ib)
!
        IF (upper_cutoff > band_min(index_exclude(jx, ib)) .AND. &
            upper_cutoff < band_max(index_exclude(jx, ib))) THEN
!         upper_cutoff is inside excluded band
          n_nu_excluded_upper = n_nu_excluded_upper + &
            NINT((upper_cutoff - band_min(index_exclude(jx, ib)))/nu_inc)
        ELSE IF (upper_cutoff >= band_max(index_exclude(jx, ib))) THEN
!         upper_cutoff is above maximum in wavenumber excluded band
          n_nu_excluded_upper = n_nu_excluded_upper + &
            NINT((band_max(index_exclude(jx, ib)) - &
            band_min(index_exclude(jx, ib)))/nu_inc)
        END IF
!
        IF (lower_cutoff > band_min(index_exclude(jx, ib)) .AND. &
            lower_cutoff < band_max(index_exclude(jx, ib))) THEN
!         lower_cutoff is inside excluded band
          n_nu_excluded_lower = n_nu_excluded_lower + &
            NINT((lower_cutoff - band_min(index_exclude(jx, ib)))/nu_inc)
        ELSE IF (lower_cutoff >= band_max(index_exclude(jx, ib))) THEN
!         lower_cutoff is above maximum in wavenumber excluded band
          n_nu_excluded_lower = n_nu_excluded_lower + &
            NINT((band_max(index_exclude(jx, ib)) - &
            band_min(index_exclude(jx, ib)))/nu_inc)
        END IF
!
      END DO
!
      j_nu_first = MAX(CEILING((lower_cutoff - nu_wgt(1))/ &
        nu_inc + 1.0_RealK) - n_nu_excluded_lower, 1)
      j_nu_last = MIN(FLOOR((upper_cutoff - nu_wgt(1))/ &
        nu_inc + 1.0_RealK) - n_nu_excluded_upper, n_nu)
!
!$OMP PARALLEL DO PRIVATE(j, kabs_line) NUM_THREADS(n_omp_threads)
      DO j = j_nu_first, j_nu_last
!
        CALL voigt_profile ( &
          nu_wgt(j), &
          adj_line_parm(i) % line_centre,   &
          adj_line_parm(i) % S_adj,         &
          adj_line_parm(i) % alpha_lorentz, &
          adj_line_parm(i) % alpha_doppler, &
          kabs_line)
!
!       If using the CKD continuum model, we reduce the line
!       absorption by its value at the cut-off: we should never
!       be in a situation where this is applied outside the cut-off.
        IF (l_ckd_cutoff) kabs_line = kabs_line-k_cutoff(i)
!
        kabs(j) = kabs(j) + kabs_line*line_prof_corr( &
          nu_wgt(j), &
          adj_line_parm(i) % line_centre, &
          adj_line_parm(i) % alpha_lorentz_air, &
          adj_line_parm(i) % alpha_lorentz_self, &
          i_line_prof_corr)
!
      END DO
!$OMP END PARALLEL DO
!
    END DO
!
!   Add the foreign continuum to the line data if necessary.
    IF (include_h2o_foreign_continuum) THEN
      ALLOCATE(k_for(n_nu))
      CALL foreign_continuum ( &
        t_calc(ipt), p_calc(ipt), &
        0.0_RealK, &
        n_nu, nu_wgt, &
        .FALSE., &
        k_for)
      kabs = kabs + k_for
      DEALLOCATE(k_for)
    END IF
!
!
!
  END SUBROUTINE calc_line_abs_int
!
!
!
  SUBROUTINE calc_xsc_abs_int

    IMPLICIT NONE

    INTEGER :: l, ll
    LOGICAL :: p_test, p_test2, t_test, t_test2
    REAL (RealK) :: p_lk(4), t_lk(4), ptorr, waveno
    REAL (RealK) :: p_int(4), t_int(4), k_int(4)
    REAL (RealK), PARAMETER :: rmdi = -1.0_RealK

!   Calculate the monochromatic absorption coefficients at each
!   wavenumber.
    WRITE(*, '(A, 1PE10.3, A8, 1PE10.3, A2)') &
      'Calculation of absorption coefficients at ', p_calc(ipt), &
      ' Pa and ', t_calc(ipt), ' K.'
    WRITE(iu_monitor, '(A, 1PE10.3, A8, 1PE10.3, A2)') &
      'Calculation of absorption coefficients at ', p_calc(ipt), &
      ' Pa and ', t_calc(ipt), ' K.'
    kabs(1:n_nu)=0.0_RealK
    ptorr=p_calc(ipt)*760.0_RealK/101325.0_RealK

!$OMP PARALLEL DO                                            &
!$OMP PRIVATE(j, i, l, ll, waveno, p_lk, t_lk, k_int, p_int, &
!$OMP         t_int, p_test, p_test2, t_test, t_test2)       &
!$OMP NUM_THREADS(n_omp_threads)
!   Look-up absorption cross-section for P/T/wavenumber
    DO j = 1, n_nu
!     Lookup wavenumer in cm-1
      waveno = nu_wgt(j) * 0.01_RealK
      p_lk = HUGE(p_lk)
      t_lk = HUGE(t_lk)
      k_int = rmdi
      DO i = 1, num_lines_in_band
        IF ( xsc(i)%head%wavenumber_min <= waveno .AND. &
             xsc(i)%head%wavenumber_max >= waveno ) THEN
          DO ll = 1, 4
            SELECT CASE (ll)
            CASE (1)
              p_test = xsc(i)%head%pressure >= ptorr .AND. &
                   ABS(xsc(i)%head%pressure - ptorr) < p_lk(ll)
              p_test2 = xsc(i)%head%pressure >= ptorr .AND. &
                   ABS(xsc(i)%head%pressure - ptorr) == p_lk(ll)
              t_test = xsc(i)%head%temperature >= t_calc(ipt)
              t_test2 = t_test .AND. &
                   ABS(xsc(i)%head%temperature - t_calc(ipt)) <= t_lk(ll)
            CASE (2)
              p_test = xsc(i)%head%pressure < ptorr .AND. &
                   ABS(xsc(i)%head%pressure - ptorr) < p_lk(ll)
              p_test2 = xsc(i)%head%pressure < ptorr .AND. &
                   ABS(xsc(i)%head%pressure - ptorr) == p_lk(ll)
              t_test = xsc(i)%head%temperature >= t_calc(ipt)
              t_test2 = t_test .AND. &
                   ABS(xsc(i)%head%temperature - t_calc(ipt)) <= t_lk(ll)
            CASE (3)
              p_test = xsc(i)%head%pressure >= ptorr .AND. &
                   ABS(xsc(i)%head%pressure - ptorr) < p_lk(ll)
              p_test2 = xsc(i)%head%pressure >= ptorr .AND. &
                   ABS(xsc(i)%head%pressure - ptorr) == p_lk(ll)
              t_test = xsc(i)%head%temperature < t_calc(ipt)
              t_test2 = t_test .AND. &
                   ABS(xsc(i)%head%temperature - t_calc(ipt)) <= t_lk(ll)
            CASE (4)
              p_test = xsc(i)%head%pressure < ptorr .AND. &
                   ABS(xsc(i)%head%pressure - ptorr) < p_lk(ll)
              p_test2 = xsc(i)%head%pressure < ptorr .AND. &
                   ABS(xsc(i)%head%pressure - ptorr) == p_lk(ll)
              t_test = xsc(i)%head%temperature < t_calc(ipt)
              t_test2 = t_test .AND. &
                   ABS(xsc(i)%head%temperature - t_calc(ipt)) <= t_lk(ll)
            END SELECT
            IF ((p_test.AND.t_test) .OR. (p_test2.AND.t_test2)) THEN
              p_lk(ll) = ABS(xsc(i)%head%pressure-ptorr)
              t_lk(ll) = ABS(xsc(i)%head%temperature-t_calc(ipt))
              l = NINT( REAL(xsc(i)%head%no_pts-1, RealK) &
                * (waveno-xsc(i)%head%wavenumber_min) &
                / (xsc(i)%head%wavenumber_max-xsc(i)%head%wavenumber_min) &
                ) + 1
              k_int(ll) = xsc(i)%data(l)
              p_int(ll) = xsc(i)%head%pressure
              t_int(ll) = xsc(i)%head%temperature
            END IF
          END DO
        END IF
      END DO
      IF (ANY(k_int /= rmdi)) THEN
        CALL bi_interp( t_int, LOG(p_int+EPSILON(p_int)), k_int, &
            rmdi, t_calc(ipt), LOG(ptorr+EPSILON(ptorr)), kabs(j) )
!       Convert from  cm2 molecule-1  to  m2 kg-1
        kabs(j)=kabs(j) &
          / (molar_weight(i_gas)*atomic_mass_unit*1.0E+04_RealK)
      END IF
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE calc_xsc_abs_int



  SUBROUTINE calc_cia_abs_int

    USE interp1d_mod, ONLY: interp1d

    IMPLICIT NONE

    REAL(RealK) :: waveno
!     Wavenumber in cm-1
    REAL(RealK) :: t_cia(num_cia_lines_in_band)
!     CIA temperatures with data at current wavenumber
    INTEGER :: n_t_cia
!     Number of CIA temperatures with data at current wavenumber
    REAL(RealK) :: kabs_t1, kabs_t2
!     Temporary storage for interpolation of absorption coefficient
    REAL(RealK) :: wgt_t1, wgt_t2
!     Temperature interpolation weights
    INTEGER :: map_cia(num_cia_lines_in_band)
!     Mapping of CIA data entries to strictly increasing tempeatures


!   Calculate the monochromatic absorption coefficients at each
!   wavenumber.
    WRITE(*, '(A, 1PE10.3, A2)') &
      'Calculation of CIA at ', t_calc(ipt), ' K.'
    WRITE(iu_monitor, '(A, 1PE10.3, A2)') &
      'Calculation of CIA at ', t_calc(ipt), ' K.'
    kabs(1:n_nu)=0.0_RealK

!   Look-up absorption cross-section for P/T/wavenumber
    DO j = 1, n_nu
!     Lookup wavenumer in cm-1
      waveno = nu_wgt(j) * 0.01_RealK
      
!     Reset interpolation quantities
      t_cia = 0.0_RealK
      n_t_cia = 0

      DO i = 1, num_cia_lines_in_band
        IF ( cia(i)%head%wavenumber_min <= waveno .AND. &
             cia(i)%head%wavenumber_max >= waveno ) THEN
!         Data for this wavenumber exists in this entry
          t_cia(i) = cia(i)%head%temperature
          n_t_cia = n_t_cia + 1
        END IF
      END DO

      IF (n_t_cia > 0) THEN

!       Sort the data so that the temperature is increasing
        CALL map_heap_func(t_cia(1:n_t_cia), map_cia(1:n_t_cia))

!       Perform interpolation in wavenumber and temperature
        IF (t_calc(ipt) <= t_cia(map_cia(1))) THEN
          kabs(j) = interp1d(cia(map_cia(1))%wavenumber, &
            cia(map_cia(1))%data, waveno, cia(map_cia(1))%head%no_pts)
        ELSE IF (t_calc(ipt) >= t_cia(map_cia(n_t_cia))) THEN
          kabs(j) = interp1d(cia(map_cia(n_t_cia))%wavenumber, &
            cia(map_cia(n_t_cia))%data, waveno, &
            cia(map_cia(n_t_cia))%head%no_pts)
        ELSE
          DO i = 1, n_t_cia - 1
            IF (t_calc(ipt) >= t_cia(map_cia(i)) .AND. &
                t_calc(ipt) <= t_cia(map_cia(i+1))) THEN
              kabs_t1 = interp1d(cia(map_cia(i))%wavenumber, &
                cia(map_cia(i))%data, waveno, cia(map_cia(i))%head%no_pts)
              kabs_t2 = interp1d(cia(map_cia(i+1))%wavenumber, &
                cia(map_cia(i+1))%data, waveno, cia(map_cia(i+1))%head%no_pts)
              wgt_t1 = (t_cia(map_cia(i+1)) - t_calc(ipt))/ &
                (t_cia(map_cia(i+1)) - t_cia(map_cia(i)))
              wgt_t2 = 1.0_RealK - wgt_t1
              kabs(j) = wgt_t1*kabs_t1 + wgt_t2*kabs_t2
            END IF
          END DO
        END IF

!       Convert from  cm5 molecule-2  to  m5 kg-2
        kabs(j) = 1.0E-10_RealK*kabs(j)/ &
          (molar_weight(i_gas_1)*molar_weight(i_gas_2)*atomic_mass_unit**2)

      END IF

    END DO

  END SUBROUTINE calc_cia_abs_int



  SUBROUTINE calc_k_opt_ref

    IMPLICIT NONE

!   Calculate for a fixed set of quadrature points.
    DO ik=1, n_k(ib)
      i0 = ig(ik-1)+1
      i1 = ig(ik)
      n_nu_k=i1-i0+1
!     Calculate the simple mean k-value across this interval.
      integ_k =  nu_inc * SUM(wgt(i0:i1))
      ALLOCATE(wgt_k(i0:i1))
      wgt_k(i0:i1)=kabs(i0:i1)*wgt(i0:i1)
      k_ave(ik,ib) = nu_inc * SUM(wgt_k(i0:i1))/integ_k
      DEALLOCATE(wgt_k)
      CALL optimal_k(n_nu_k, &
        nu_inc, kabs(i0:i1), wgt(i0:i1), &
        integ_k, k_ave(ik,ib), tol, &
        l_fit_cont_data .AND. l_cont_line_abs_weight, kabs_lines(i0:i1), &
        l_wgt_scale_sqrt, u_wgt_scale, k_opt(ik,ib), err_norm, ierr)
      IF (ierr /= i_normal) RETURN

    ENDDO

  END SUBROUTINE calc_k_opt_ref



  SUBROUTINE ck_trans_fit

    IMPLICIT NONE 

!   Integrate the sorted weightings across the band.
    integ_wgt=nu_inc * SUM(wgt(1:n_nu))

!   Scaling will be carried out on the transmissions over a
!   range of path-lengths: these are chosen with regard to
!   the range of k-terms, looking at transmissions between
!   tol and 1-tol.
    umax=umax_kopt
    IF (k_opt(1,ib) > 0.0_RealK) umax = &
      MIN(umax_kopt, -LOG(tol)/k_opt(1,ib))
    umin=umin_kopt
    IF (k_opt(n_k(ib),ib) > 0.0_RealK) umin = &
      MAX(umin_kopt, -LOG(1.0_RealK-tol)/k_opt(n_k(ib),ib))
!   Optimization over paths is likely to be simpler with 
!   fewer k-terms.
    n_path=2*n_k(ib)+1
    DO i=1, n_path
      u_l_ref(i)=umin*EXP( REAL(i-1, RealK) * &
      LOG(umax/umin) / REAL((n_path-1), RealK))
    ENDDO
    ALLOCATE(trans_ref(n_path))
    ALLOCATE(trans_calc(n_path))
    IF (l_scale_pT) ALLOCATE(trans_pt_k(n_path, n_pt_pair, n_k(ib)))
!
!   Calculate the reference fitted transmissions
    trans_ref(:) = &
      trans_k_dist(n_nu, kabs, nu_inc, wgt, integ_wgt, n_path, u_l_ref)
    DO i=1, n_path
      trans_calc(i)=0.0_RealK
      DO ik=1, n_k(ib)
        trans_calc(i) = trans_calc(i) + &
          w_k(ik, ib) * EXP(-k_opt(ik, ib) * u_l_ref(i))
      ENDDO
    ENDDO

  END SUBROUTINE ck_trans_fit

  SUBROUTINE ck_fit_k

    IMPLICIT NONE

    IF (.NOT.l_load_map) THEN
!     Sort the absorption and weighting functions
!     with reordering at each pressure level
      DO ipb=1,n_pre
        CALL map_heap_func( kabs(pmap(pstart(ipb):pend(ipb))), &
          gmap(pstart(ipb):pend(ipb)))
        gmap(pstart(ipb):pend(ipb)) = pstart(ipb)-1 + &
          gmap(pstart(ipb):pend(ipb))
      END DO
      map = pmap(gmap)
    END IF

    kabs=kabs(map)
    wgt=wgt(map)
    IF (l_fit_cont_data .AND. l_cont_line_abs_weight) &
      kabs_lines=kabs_lines(map)

!   Integrate the sorted weightings across the band.
    integ_wgt=nu_inc * SUM(wgt(1:n_nu))

!   The quadrature points are determined for the first pair of
!   p and T. We may now determine fits for all pairs of p and T
!   or find a scaling function to fit the transmissions.
!   Calculate for a fixed set of quadrature points.
    DO ik=1, n_k(ib)
      i0 = ig(ik-1)+1
      i1 = ig(ik)
      n_nu_k=i1-i0+1
!     Calculate the simple mean k-value across this interval.
      integ_k =  nu_inc * SUM(wgt(i0:i1))
      ALLOCATE(wgt_k(i0:i1))
      wgt_k(i0:i1)=kabs(i0:i1)*wgt(i0:i1)
      k_ave_tmp(ik) = nu_inc * SUM(wgt_k(i0:i1))/integ_k
      DEALLOCATE(wgt_k)
      CALL optimal_k(n_nu_k, &
        nu_inc, kabs(i0:i1), wgt(i0:i1), &
        integ_k, k_ave_tmp(ik), tol, &
        l_fit_cont_data .AND. l_cont_line_abs_weight, kabs_lines(i0:i1), &
        l_wgt_scale_sqrt, u_wgt_scale, k_opt_tmp(ik), err_norm, ierr)
      IF (ierr /= i_normal) RETURN

      IF (l_scale_pT) THEN
!       Fitting is done over the same amounts of absorber as at the
!       reference conditions.
        u_l(:, ipt) = u_l_ref(:)
!       Store the transmissions for use in determining a fit
        trans_pt_k(:, ipt, ik) = &
          trans_k_dist(n_nu_k, &
            kabs(i0:i1), nu_inc, wgt(i0:i1), &
            integ_k, n_path, u_l)
      ENDIF
    ENDDO

    kopt_all(:,ipt,ibb)=k_opt_tmp(:)

  END SUBROUTINE ck_fit_k
!
!
  SUBROUTINE calc_self_trans_int
!
!
    IMPLICIT NONE 
!
!   Local variables
    INTEGER :: i_pp
!     Loop variable
    INTEGER :: il
!     Starting position of data in the long array
    INTEGER :: ih
!     Finishing position of data in the long array
    INTEGER :: ju
!     Indexing variable
    INTEGER :: jv
!     Indexing variable
    REAL  (RealK) :: e_sat
!     Saturation vapour pressure of water vapour
    REAL  (RealK) :: pp
!     Partial pressure of water vapour
    REAL  (RealK), Dimension(:), Allocatable :: trans_line
!     Line transmissions
    REAL  (RealK), Dimension(:), Allocatable :: k_self
!     Self-broadened continuum absorption at individual frequencies
    REAL  (RealK), Dimension(:), Allocatable :: ktot
!     Total (lines + continuum) absorption at individual frequencies
!
    INTERFACE
!
      SUBROUTINE exponent_fit_90(n, u, trans, k, ierr)
!
        USE realtype_rd
!
        INTEGER, Intent(INOUT) :: ierr
        INTEGER, Intent(IN) :: n
        REAL  (RealK), Intent(IN), Dimension(:) :: u
        REAL  (RealK), Intent(IN), Dimension(:) :: trans
        REAL  (RealK), Intent(OUT) :: k
!
      END SUBROUTINE exponent_fit_90
!
!
      FUNCTION trans_k_dist(n_nu, k, nu_inc, wgt, integ_wgt, n_path, u) &
        RESULT (trans)
!  
        USE realtype_rd
!  
        INTEGER, Intent(IN) :: n_nu
        INTEGER, Intent(IN) :: n_path
        REAL  (RealK), Intent(IN), Dimension(n_nu) :: k
        REAL  (RealK), Intent(IN) :: nu_inc
        REAL  (RealK), Intent(IN), Dimension(n_nu) :: wgt
        REAL  (RealK), Intent(IN) :: integ_wgt
        REAL  (RealK), Intent(IN), Dimension(n_path) :: u
!  
        REAL  (RealK), Dimension(n_path) :: trans
   
!  
      END FUNCTION trans_k_dist
!
!
    END INTERFACE
!
!
!
!   Firstly, calculate the line transmission along these paths.
    ALLOCATE(trans_line(n_path_c))
    trans_line = &
      trans_k_dist(n_nu, kabs, nu_inc, wgt, integ_wgt, n_path_c, u_c)
!
!   The continuum is calculated at a range of partial pressures up to 
!   the saturation value.
    e_sat = sat_vap_press(t_calc(ipt), p_calc(ipt))
    ALLOCATE(k_self(n_nu))
    ALLOCATE(ktot(n_nu))

!$OMP PARALLEL DO                                              &
!$OMP PRIVATE(i_pp, pp, k_self, ktot, trans_c, il, ih, ju, jv) &
!$OMP NUM_THREADS(n_omp_threads)
    DO i_pp = 1, n_pp
!
      pp = e_sat * REAL(i_pp, RealK) / REAL(n_pp, RealK) 
!
!     Calculate the self-broadened continuum coefficients at this
!     partial pressure.
      CALL self_continuum ( &
        t_calc(ipt), p_calc(ipt), &
        pp, &
        n_nu, nu_wgt, &
        .FALSE., &
        k_self)
!
!     Now calculate the total transmission
      ktot = kabs + k_self
      trans_c = &
        trans_k_dist(n_nu, ktot, nu_inc, wgt, integ_wgt, n_path_c, u_c)
!
!     Feed these data into the full array of transmissions.
      il = 1 + (i_pp - 1) * n_path_c
      ih = i_pp * n_path_c
      trans_fit_c(il:ih, ipt) = trans_c / (trans_line + TINY(trans_line) )
!     Eventually, the continuum coefficient will be in units of 
!     m5/(mol.kg).
      u_fit_c(il:ih, ipt) = u_c * &
        ( pp / (molar_gas_constant * t_calc(ipt)) )
!
!     Where the line transmission is very small, the calculated continuum
!     transmission (trans_c/trans_line) may be unreliable. In such cases
!     we set the transmission to 1 and the amount to 0. This is not the
!     most efficient course of action, but allows all sets of data passed
!     to the fitting routine to be of the same length.
      DO ju = 1, n_path_c
        IF (trans_line(ju) < SQRT(EPSILON(trans_line))) THEN
          jv = ju + (i_pp - 1) * n_path_c
          u_fit_c(jv, ipt)     = 0.0_RealK
          trans_fit_c(jv, ipt) = 1.0_RealK
        ENDIF
      ENDDO

    ENDDO
!$OMP END PARALLEL DO

    DEALLOCATE(trans_line)
    DEALLOCATE(k_self)
    DEALLOCATE(ktot)

    IF (ipt == ipt_ref) THEN
      CALL exponent_fit_90(n_path_c * n_pp, u_fit_c(:, ipt), &
                        trans_fit_c(:, ipt), &
                        k_opt_self(ib), ierr)
      IF (ierr /= i_normal) RETURN
    ENDIF
!
!
!
  END SUBROUTINE calc_self_trans_int
!
!
!
  SUBROUTINE calc_frn_trans_int
!
!
    IMPLICIT NONE 
!
!   Local variables
    INTEGER :: i_pp
!     Loop variable
    INTEGER :: il
!     Starting position of data in the long array
    INTEGER :: ih
!     Finishing position of data in the long array
    INTEGER :: ju
!     Indexing variable
    INTEGER :: jv
!     Indexing variable
    REAL  (RealK) :: e_sat
!     Saturation vapour pressure of water vapour
    REAL  (RealK) :: pp
!     Partial pressure of water vapour
    REAL  (RealK), Dimension(:), Allocatable :: trans_line
!     Line transmissions
    REAL  (RealK), Dimension(:), Allocatable :: k_frn
!     Foreign-broadened continuum absorption at individual frequencies
    REAL  (RealK), Dimension(:), Allocatable :: ktot
!     Total (lines + continuum) absorption at individual frequencies
!
    INTERFACE
!
      SUBROUTINE exponent_fit_90(n, u, trans, k, ierr)
!
        USE realtype_rd
!
        INTEGER, Intent(INOUT) :: ierr
        INTEGER, Intent(IN) :: n
        REAL  (RealK), Intent(IN), Dimension(:) :: u
        REAL  (RealK), Intent(IN), Dimension(:) :: trans
        REAL  (RealK), Intent(OUT) :: k
!
      END SUBROUTINE exponent_fit_90
!
!
      FUNCTION trans_k_dist(n_nu, k, nu_inc, wgt, integ_wgt, n_path, u) &
        RESULT (trans)
!  
        USE realtype_rd
!  
        INTEGER, Intent(IN) :: n_nu
        INTEGER, Intent(IN) :: n_path
        REAL  (RealK), Intent(IN), Dimension(n_nu) :: k
        REAL  (RealK), Intent(IN) :: nu_inc
        REAL  (RealK), Intent(IN), Dimension(n_nu) :: wgt
        REAL  (RealK), Intent(IN) :: integ_wgt
        REAL  (RealK), Intent(IN), Dimension(n_path) :: u
!  
        REAL  (RealK), Dimension(n_path) :: trans
   
!  
      END FUNCTION trans_k_dist
!
!
    END INTERFACE
!
!
!
!   Firstly, calculate the line transmission along these paths.
    ALLOCATE(trans_line(n_path_c))
    trans_line = &
      trans_k_dist(n_nu, kabs, nu_inc, wgt, integ_wgt, n_path_c, u_c)
!
!   The continuum is calculated at a range of partial pressures up to 
!   the saturation value.
    e_sat = sat_vap_press(t_calc(ipt), p_calc(ipt))
    ALLOCATE(k_frn(n_nu))
    ALLOCATE(ktot(n_nu))

!$OMP PARALLEL DO                                             &
!$OMP PRIVATE(i_pp, pp, k_frn, ktot, trans_c, il, ih, ju, jv) &
!$OMP NUM_THREADS(n_omp_threads)
    DO i_pp = 1, n_pp
!
      pp = e_sat * REAL(i_pp, RealK) / REAL(n_pp, RealK) 
!
!     Calculate the foreign-broadened continuum coefficients at this
!     partial pressure.
      CALL foreign_continuum ( &
        t_calc(ipt), p_calc(ipt), &
        pp, &
        n_nu, nu_wgt, &
        .FALSE., &
        k_frn)
!
!     Now calculate the total transmission
      ktot = kabs + k_frn
      trans_c = &
        trans_k_dist(n_nu, ktot, nu_inc, wgt, integ_wgt, n_path_c, u_c)
!
!     Feed these data into the full array of transmissions.
      il = 1 + (i_pp - 1) * n_path_c
      ih = i_pp * n_path_c
      trans_fit_c(il:ih, ipt) = trans_c / (trans_line + TINY(trans_line) )
!     Eventually, the continuum coefficient will be in units of 
!     m5/(mol.kg).
      u_fit_c(il:ih, ipt) = u_c * &
        ( (p_calc(ipt) - pp) / (molar_gas_constant * t_calc(ipt)) )
!
!     Where the line transmission is very small, the calculated continuum
!     transmission (trans_c/trans_line) may be unreliable. In such cases
!     we set the transmission to 1 and the amount to 0. This is not the
!     most efficient course of action, but allows all sets of data passed
!     to the fitting routine to be of the same length.
      DO ju = 1, n_path_c
        IF (trans_line(ju) < SQRT(EPSILON(trans_line))) THEN
          jv = ju + (i_pp - 1) * n_path_c
          u_fit_c(jv, ipt)     = 0.0_RealK
          trans_fit_c(jv, ipt) = 1.0_RealK
        ENDIF
      ENDDO

    ENDDO
!$OMP END PARALLEL DO

    DEALLOCATE(trans_line)
    DEALLOCATE(k_frn)
    DEALLOCATE(ktot)

    IF (ipt == ipt_ref) THEN
      CALL exponent_fit_90(n_path_c * n_pp, u_fit_c(:, ipt), &
                        trans_fit_c(:, ipt), &
                        k_opt_frn(ib), ierr)
      IF (ierr /= i_normal) RETURN
    ENDIF
!
!
!
  END SUBROUTINE calc_frn_trans_int
!
!
!
  SUBROUTINE fit_scale_line_int
!
!
    IMPLICIT NONE 
!
    DO ik=1, n_k(ib)
!
      WRITE(iu_monitor, '(/a, /a, i3)') &
       "===================", &
       "Calculating scaling for term ", ik
!
!     Ininitialize the parameters of the scaling function.
      SELECT CASE(i_scale_function)
        CASE (IP_scale_power_law)
          scale_vector(1:2, ik, ib) = &
            (/ 1.0_RealK, 0.5_RealK /)
        CASE (IP_scale_power_quad)
          scale_vector(1:3, ik, ib) = &
          (/ 1.0_RealK, 0.0_RealK, 0.0_RealK /)
        CASE (IP_scale_doppler_quad)
          scale_vector(1:4, ik, ib) = &
            (/ 1.0_RealK, -2.0_RealK, 0.0_RealK, 0.0_RealK /)
      END SELECT
!
!     No if the absorption coefficient for the term is 0, there is
!     no optimal scaling function, so no scaling function can be
!     determined. We relate this to the machine's precision.
      IF (k_opt(ik, ib) > SQRT(EPSILON(k_opt))) THEN
        CALL scale_ck_fit_90(ierr, iu_monitor, &
          n_pt_pair, &
          n_path, u_l(:, 1: n_pt_pair), &
          trans_pt_k(:, 1: n_pt_pair, ik), &
          p_calc(1: n_pt_pair), t_calc(1: n_pt_pair), &
          p_calc(ipt_ref), t_calc(ipt_ref), kopt_all(ik, ipt_ref, ibb), &
          i_type_residual, i_scale_function, scale_vector(:, ik, ib), &
          rms_residual)
        IF (ierr == i_abort_calculation) THEN
!         Recover from mild errors.
          ierr=i_normal
        ELSE IF (ierr /= i_normal) THEN
          RETURN
        ENDIF
      ENDIF

!     If the fit has been done using a single pressure then the pressure
!     scaling will not have changed from the initial value and must now
!     be set to zero (no pressure scaling).
      IF (n_p == 1) THEN
        SELECT CASE(i_scale_function)
          CASE (IP_scale_power_law, IP_scale_power_quad)
            scale_vector(1, ik, ib) = 0.0_RealK
        END SELECT
      END IF

    ENDDO
!
!
!
  END SUBROUTINE fit_scale_line_int

  SUBROUTINE fit_scale_line_int2

    IMPLICIT NONE

    INTEGER :: index_k_ref, i_scale_function2
    REAL  (RealK) :: err_norm_old

    IF (i_scale_function == ip_scale_dbl_pow_law)  &
        i_scale_function2 = ip_scale_power_law
    IF (i_scale_function == ip_scale_dbl_pow_quad) &
        i_scale_function2 = ip_scale_power_quad
    IF (i_scale_function == ip_scale_dbl_dop_quad) &
        i_scale_function2 = ip_scale_doppler_quad

    DO ik=1, n_k(ib)

      WRITE(iu_monitor, '(/a, /a, i3)') &
       "===================", &
       "Calculating scaling for term ", ik

      err_norm=HUGE(err_norm)*0.5

!     Loop over pressures
      DO i=2,n_p-1
        index_k_ref=index_pt_ref(i)

        IF (l_debug) &
          print*,'Fitting 2 scaling functions split at pressure: ', &
          p_calc(index_k_ref)
        
        err_norm_old=err_norm
        
!       Fit 2 scaling functions, one either side of the maximum k.
        
!       Ininitialize the parameters of the first scaling function.
        SELECT CASE(i_scale_function2)
          CASE (IP_scale_power_law)
            scale_vector(1:2, ik, ib) = &
              (/ 1.0_RealK, 0.5_RealK /)
          CASE (IP_scale_power_quad)
            scale_vector(1:3, ik, ib) = &
            (/ 1.0_RealK, 0.0_RealK, 0.0_RealK /)
          CASE (IP_scale_doppler_quad)
            scale_vector(1:4, ik, ib) = &
              (/ 1.0_RealK, -2.0_RealK, 0.0_RealK, 0.0_RealK /)
        END SELECT
        
!       Now if the absorption coefficient for the term is 0, there is
!       no optimal scaling function, so no scaling function can be
!       determined. We relate this to the machine's precision.
        IF (kopt_all(ik, index_k_ref, ibb) > SQRT(EPSILON(k_opt)) &
          .AND. p_calc(index_k_ref) > p_calc(1) ) THEN
          CALL scale_ck_fit_90(ierr, iu_monitor, &
            index_k_ref-1, &
            n_path, u_l(:, 1:index_k_ref-1), &
            trans_pt_k(:, 1:index_k_ref-1, ik), &
            p_calc(1:index_k_ref-1), t_calc(1:index_k_ref-1), &
            p_calc(index_k_ref), t_calc(index_k_ref), &
            kopt_all(ik, index_k_ref, ibb), &
            i_type_residual, i_scale_function2, scale_vector(1: &
              n_scale_variable(i_scale_function2), ik, ib), &
            rms_residual)
          IF (ierr == i_abort_calculation) THEN
!           Recover from mild errors.
            ierr=i_normal
          ELSE IF (ierr /= i_normal) THEN
            RETURN
          ENDIF
          err_norm=rms_residual*(index_k_ref-1)
        ENDIF
        
!       Ininitialize the parameters of the second scaling function.
        SELECT CASE(i_scale_function2)
          CASE (IP_scale_power_law)
            scale_vector(3:4, ik, ib) = &
              (/ 1.0_RealK, 0.5_RealK /)
          CASE (IP_scale_power_quad)
            scale_vector(4:6, ik, ib) = &
            (/ 1.0_RealK, 0.0_RealK, 0.0_RealK /)
          CASE (IP_scale_doppler_quad)
            scale_vector(5:8, ik, ib) = &
              (/ 1.0_RealK, -2.0_RealK, 0.0_RealK, 0.0_RealK /)
        END SELECT
        IF (kopt_all(ik, index_k_ref, ibb) > SQRT(EPSILON(k_opt)) &
          .AND. p_calc(index_k_ref) < p_calc(n_pt_pair) ) THEN
          CALL scale_ck_fit_90(ierr, iu_monitor, &
            n_pt_pair-index_k_ref, &
            n_path, u_l(:, index_k_ref+1: n_pt_pair), &
            trans_pt_k(:, index_k_ref+1: n_pt_pair, ik), &
            p_calc(index_k_ref+1: n_pt_pair), &
            t_calc(index_k_ref+1: n_pt_pair), &
            p_calc(index_k_ref), t_calc(index_k_ref), &
            kopt_all(ik, index_k_ref, ibb), &
            i_type_residual, i_scale_function2, &
            scale_vector(n_scale_variable(i_scale_function2)+1: &
              n_scale_variable(i_scale_function2)*2, ik, ib), &
            rms_residual)
          IF (ierr == i_abort_calculation) THEN
            ierr=i_normal
          ELSE IF (ierr /= i_normal) THEN
            RETURN
          ENDIF
          err_norm=err_norm + rms_residual*(n_pt_pair-index_k_ref)
        ENDIF
        
        IF (err_norm < err_norm_old) THEN
!         Set reference P, T, and absorption for each k-term:
          scale_vector(n_scale_variable(i_scale_function2)*2+1,ik,ib) = &
            p_calc(index_k_ref)
          scale_vector(n_scale_variable(i_scale_function2)*2+2,ik,ib) = &
            t_calc(index_k_ref)
          k_opt(ik,ib) = kopt_all(ik, index_k_ref, ibb)
        END IF

      END DO

    ENDDO

  END SUBROUTINE fit_scale_line_int2



  SUBROUTINE fit_scale_cont_int
!
!
    IMPLICIT NONE 
!
    WRITE(iu_monitor, '(/a, /a)') &
     "===================", &
     "Calculating scaling for continuum "
!
!   Ininitialize the parameters of the scaling function.
    SELECT CASE(i_scale_function)
      CASE (IP_scale_power_law)
        scale_cont(1:2, 1, ib) = &
          (/ 1.0_RealK, 0.5_RealK /)
      CASE (IP_scale_power_quad)
        scale_cont(1:3, 1, ib) = &
        (/ 1.0_RealK, 0.0_RealK, 0.0_RealK /)
      CASE (IP_scale_doppler_quad)
        scale_cont(1:4, 1, ib) = &
          (/ 1.0_RealK, -2.0_RealK, 0.0_RealK, 0.0_RealK /)
    END SELECT
!
!   Note that if the absorption coefficient for the term is 0, there is
!   no optimal scaling function, so no scaling function can be
!   determined. We relate this to the machine's precision.
    IF (k_cont > SQRT(EPSILON(k_cont))) THEN
      CALL scale_ck_fit_90(ierr, iu_monitor, &
        n_pt_pair, &
        n_path_c * n_pp, u_fit_c(:, 1: n_pt_pair), &
        trans_fit_c(:, 1: n_pt_pair), &
        p_calc(1: n_pt_pair), t_calc(1: n_pt_pair), &
        p_calc(ipt_ref), t_calc(ipt_ref), k_cont, &
        i_type_residual, i_scale_function, scale_cont(:, 1, ib), &
        rms_residual)
      IF (ierr == i_abort_calculation) THEN
!       Recover from mild errors.
        ierr=i_normal
      ELSE IF (ierr /= i_normal) THEN
        RETURN
      ENDIF
    ENDIF
!
!
!
  END SUBROUTINE fit_scale_cont_int


  SUBROUTINE input_lbl_band_cdf_init

    USE netcdf
    IMPLICIT NONE
    INTEGER :: dimid1                  ! dimension ID
    INTEGER :: varid                   ! variable ID
    CHARACTER(LEN=2) :: dim_name='nu'  ! Name of wavenumber dimension
    INTEGER :: dim_len                 ! Length of wavenumber dimension
    REAL (RealK) :: p_calc_in(n_pt_pair), t_calc_in(n_pt_pair)
    REAL (RealK) :: max_nu_wgt, min_nu_wgt
    REAL (RealK) :: gas_frac_in(n_gas_frac)

!   Open the file for reading
    CALL nf(nf90_open(TRIM(file_lbl),NF90_NOWRITE,ncidin_lbl))

!   Get length of wavenumber dimension
    CALL nf(nf90_inq_dimid(ncidin_lbl, dim_name, dimid1))
    CALL nf(nf90_inquire_dimension(ncidin_lbl, dimid1, dim_name, dim_len))

!   Allocate array with all wavenumbers
    ALLOCATE(nu_wgt_all(dim_len))

!   Read and get step in wavenumber array 
    CALL nf(nf90_inq_varid(ncidin_lbl,'nu',varid))
    CALL nf(nf90_get_att(ncidin_lbl,varid,'step',nu_inc))
    CALL nf(nf90_get_var(ncidin_lbl,varid,nu_wgt_all))

    min_nu_wgt = MINVAL(nu_wgt_all)
    max_nu_wgt = MAXVAL(nu_wgt_all)
!   Check bands against wavenumbers available
    DO ibb=1,n_selected_band
      ib = list_band(ibb)
      IF ((band_min(ib) - min_nu_wgt) <= &
        (-nu_inc - min_nu_wgt*EPSILON(nu_inc)).OR. &
        (band_max(ib) - max_nu_wgt) >= &
        (nu_inc + max_nu_wgt*EPSILON(nu_inc))) THEN

        WRITE(*,'(a)') 'Band limits out of wavenumber range: '
        WRITE(*,*) 'Band Min, Band Max: ',band_min(ib),band_max(ib)
        WRITE(*,*) 'Wavenumber range: ', &
          min_nu_wgt - nu_inc/2.0_RealK, &
          max_nu_wgt + nu_inc/2.0_RealK
        WRITE(*,'(a)') 'Please adjust band limits.'
        ierr = i_err_fatal
        RETURN
      END IF
    END DO

!   Read pressures and temperatures
    CALL nf(nf90_inq_varid(ncidin_lbl,'p_calc',varid)) 
    CALL nf(nf90_get_var(ncidin_lbl,varid,p_calc_in))
    CALL nf(nf90_inq_varid(ncidin_lbl,'t_calc',varid))
    CALL nf(nf90_get_var(ncidin_lbl,varid,t_calc_in))

!   Check pressures and temperatures in lbl file
    DO ipt=1,n_pt_pair
!     Correct in case single precision is used in file
      IF (abs(p_calc_in(ipt)-p_calc(ipt)) >= &
        p_calc(ipt)*EPSILON(REAL(1e+0)) .OR. &
        abs(t_calc_in(ipt)-t_calc(ipt)) >= t_calc(ipt)*EPSILON(REAL(1e+0))) THEN
        WRITE(*,'(a)') 'P,T points in lbl file do not match input values.'
        ierr = i_err_fatal
        RETURN
      ENDIF
    END DO

    IF (l_self_broadening) THEN
!     Read gas fractions
      CALL nf(nf90_inq_varid(ncidin_lbl,'gas_frac',varid)) 
      CALL nf(nf90_get_var(ncidin_lbl,varid,gas_frac_in))

!     Check gas fractions in lbl file
      DO igf=1,n_gas_frac
!       Correct in case single precision is used in file
        IF (ABS(gas_frac_in(igf)-gas_frac(igf)) > &
          gas_frac(igf)*EPSILON(REAL(1e+0))) THEN
          WRITE(*,'(a)') 'Gas fractions in lbl file do not match input values.'
          ierr = i_err_fatal
          RETURN
        ENDIF
      END DO
    END IF

  END SUBROUTINE input_lbl_band_cdf_init


  SUBROUTINE input_lbl_band_cdf

    USE netcdf
    IMPLICIT NONE
    INTEGER :: varid                   ! variable ID
    INTEGER :: i_nu, nu_index, nu_count
    INTEGER :: nu_wgt_map(n_nu), nu_wgt_all_map(SIZE(nu_wgt_all))

    nu_wgt_map=NINT((nu_wgt-nu_wgt(1))/nu_inc)
    nu_wgt_all_map=NINT((nu_wgt_all-nu_wgt(1))/nu_inc)

!   Get absorption coefficients for current band
    CALL nf(nf90_inq_varid(ncidin_lbl,'kabs',varid))

    i_nu=1
    DO
      IF (i_nu > n_nu) EXIT
      nu_index = MINLOC(nu_wgt_all_map, 1, &
        nu_wgt_all_map == nu_wgt_map(i_nu))
      IF (nu_index == 0) THEN
        WRITE(*,'(a)') &
          '***Error: required wavenumber not found in LbL file'
        STOP
      END IF
      nu_count=1
      DO
!       Read contiguous blocks of data from LbL file
        IF (i_nu     + nu_count <= n_nu .AND. &
            nu_index + nu_count <= SIZE(nu_wgt_all)) THEN
          IF (nu_wgt_map(i_nu+nu_count) == &
            nu_wgt_all_map(nu_index+nu_count)) THEN
            nu_count=nu_count+1
          ELSE
            EXIT
          END IF
        ELSE
          EXIT
        END IF
      END DO
      WRITE(*,'(2(a,f12.2))') 'Reading wavenumbers:', &
        nu_wgt(i_nu),' to', nu_wgt(i_nu+nu_count-1)
      IF (l_self_broadening) THEN
        CALL nf(nf90_get_var(ncidin_lbl,varid, &
          kabs_all_sb(i_nu:i_nu+nu_count-1,:,:), &
          start=(/nu_index, 1, 1/), &
          count=(/nu_count, n_pt_pair, n_gas_frac/)))
      ELSE
        CALL nf(nf90_get_var(ncidin_lbl,varid, &
          kabs_all(i_nu:i_nu+nu_count-1,:), &
          start=(/nu_index, 1/), &
          count=(/nu_count, n_pt_pair/)))
      END IF
      i_nu=i_nu+nu_count
    END DO

  END SUBROUTINE input_lbl_band_cdf


  SUBROUTINE input_map_band_cdf_init

    USE netcdf
    IMPLICIT NONE
    INTEGER :: varid

!   Open the file for reading
    CALL nf(nf90_open(TRIM(file_map),NF90_NOWRITE,ncidin_map))

!   Read number of k-terms
    CALL nf(nf90_inq_varid(ncidin_map,'n_k',varid))
    CALL nf(nf90_get_var(ncidin_map,varid,n_k,start=(/1/),count=(/n_band/)))

  END SUBROUTINE input_map_band_cdf_init


  SUBROUTINE input_map_band_cdf

    USE netcdf
    IMPLICIT NONE
    INTEGER :: varid                   ! variable ID

!   Read g-points
    CALL nf(nf90_inq_varid(ncidin_map,'ig',varid))
    CALL nf(nf90_get_var(ncidin_map,varid,ig, &
      start=(/1,ib/),count=(/n_k(ib)+1,1/)))

!   Read mapping of frequencies
    CALL nf(nf90_inq_varid(ncidin_map,'map',varid))
    CALL nf(nf90_get_var(ncidin_map,varid,map, &
      start=(/1,ib/),count=(/n_nu,1/)))

!   Read reference k-term weights
    CALL nf(nf90_inq_varid(ncidin_map,'w_k',varid))
    CALL nf(nf90_get_var(ncidin_map,varid,w_k(1:n_k(ib),ib), &
      start=(/1,ib/),count=(/n_k(ib),1/)))

  END SUBROUTINE input_map_band_cdf


  SUBROUTINE input_wgt_band_cdf

    USE netcdf
    IMPLICIT NONE
    INTEGER :: varid
    INTEGER :: n_pre_arr(1)

!   Read number of pressure bins
    CALL nf(nf90_inq_varid(ncidin_map,'n_pre',varid))
    CALL nf(nf90_get_var(ncidin_map,varid,n_pre_arr, &
      start=(/ib/),count=(/1/)))
    n_pre = n_pre_arr(1)

!   Read g-points
    CALL nf(nf90_inq_varid(ncidin_map,'ig',varid))
    CALL nf(nf90_get_var(ncidin_map,varid,ig, &
      start=(/1,ib/),count=(/n_k(ib)+1,1/)))

!   Read reference k-term weights
    CALL nf(nf90_inq_varid(ncidin_map,'w_k',varid))
    CALL nf(nf90_get_var(ncidin_map,varid,w_k(1:n_k(ib),ib), &
      start=(/1,ib/),count=(/n_k(ib),1/)))

!   Read pressure bin data
    CALL nf(nf90_inq_varid(ncidin_map,'pstart',varid))
    CALL nf(nf90_get_var(ncidin_map,varid,pstart, &
      start=(/1,ib/),count=(/n_pre,1/)))
    CALL nf(nf90_inq_varid(ncidin_map,'pend',varid))
    CALL nf(nf90_get_var(ncidin_map,varid,pend, &
      start=(/1,ib/),count=(/n_pre,1/)))
    CALL nf(nf90_inq_varid(ncidin_map,'pcount',varid))
    CALL nf(nf90_get_var(ncidin_map,varid,pcount, &
      start=(/1,ib/),count=(/n_pre,1/)))

!   Read pressure bin mapping
    CALL nf(nf90_inq_varid(ncidin_map,'pmap',varid))
    CALL nf(nf90_get_var(ncidin_map,varid,pmap, &
      start=(/1,ib/),count=(/n_nu,1/)))

  END SUBROUTINE input_wgt_band_cdf


  SUBROUTINE output_lbl_band_cdf_init

    USE netcdf
    IMPLICIT NONE
    INTEGER :: dimid1, dimid2, dimid3  ! dimension ID
    INTEGER :: varid                   ! variable ID
    INTEGER :: n_nu_tot                ! Total number of frequency points

    ! Calculate total number of frequency points
    n_nu_tot=0
    DO ibb=1,n_selected_band
      ib=list_band(ibb)
      band_width = band_max(ib) - band_min(ib)
      DO jx = 1, n_band_exclude(ib)
        band_width = band_width + &
          band_min(index_exclude(jx, ib)) - band_max(index_exclude(jx, ib))
      ENDDO
      n_nu=NINT(band_width/nu_inc)
      n_nu_tot=n_nu_tot+n_nu
    END DO

!   Create the file and open for writing
    CALL nf(nf90_create(TRIM(file_lbl),NF90_NOCLOBBER,ncidout_lbl))

!   Create dimensions
    CALL nf(nf90_def_dim(ncidout_lbl, 'nu', n_nu_tot, dimid1))
    CALL nf(nf90_def_dim(ncidout_lbl, 'pt_pair', n_pt_pair, dimid2))
    IF (l_self_broadening) &
      CALL nf(nf90_def_dim(ncidout_lbl, 'gas_frac', n_gas_frac, dimid3))

!   Create variables and write p_calc and t_calc
    CALL nf(nf90_def_var(ncidout_lbl, 'p_calc', NF90_FLOAT, dimid2, varid))
    CALL nf(nf90_put_att(ncidout_lbl, varid, 'title', 'pressure' ))
    CALL nf(nf90_put_att(ncidout_lbl, varid, 'long_name', 'pressure' ))
    CALL nf(nf90_put_att(ncidout_lbl, varid, 'units', 'Pa' ))
    CALL nf(nf90_enddef(ncidout_lbl))
    CALL nf(nf90_put_var(ncidout_lbl, varid, p_calc(1:n_pt_pair) ))
    CALL nf(nf90_redef(ncidout_lbl))

    CALL nf(nf90_def_var(ncidout_lbl, 't_calc', NF90_FLOAT, dimid2, varid))
    CALL nf(nf90_put_att(ncidout_lbl, varid, 'title', 'temperature' ))
    CALL nf(nf90_put_att(ncidout_lbl, varid, 'long_name', 'temperature' ))
    CALL nf(nf90_put_att(ncidout_lbl, varid, 'units', 'K' ))
    CALL nf(nf90_enddef(ncidout_lbl))
    CALL nf(nf90_put_var(ncidout_lbl, varid, t_calc(1:n_pt_pair) ))
    CALL nf(nf90_redef(ncidout_lbl))

    CALL nf(nf90_def_var(ncidout_lbl, 'nu', NF90_DOUBLE, dimid1, varid))
    CALL nf(nf90_put_att(ncidout_lbl, varid, 'title', 'wavenumber'))
    CALL nf(nf90_put_att(ncidout_lbl, varid, 'long_name', 'wavenumber'))
    CALL nf(nf90_put_att(ncidout_lbl, varid, 'units', 'm-1'))
    CALL nf(nf90_put_att(ncidout_lbl, varid, 'step', nu_inc))

    IF (l_self_broadening) THEN
      CALL nf(nf90_def_var(ncidout_lbl, 'gas_frac', NF90_FLOAT, dimid3, varid))
      CALL nf(nf90_put_att(ncidout_lbl, varid, 'title', 'gas fraction'))
      CALL nf(nf90_put_att(ncidout_lbl, varid, 'long_name', 'gas fraction'))
      CALL nf(nf90_enddef(ncidout_lbl))
      CALL nf(nf90_put_var(ncidout_lbl, varid, &
        gas_frac(1:n_gas_frac) ))
      CALL nf(nf90_redef(ncidout_lbl))
    END IF

    IF (l_output_reference_weight) THEN
      CALL nf(nf90_def_var(ncidout_lbl, 'wgt', NF90_FLOAT, dimid1, varid))
      CALL nf(nf90_put_att(ncidout_lbl, varid, 'title', 'weight'))
      CALL nf(nf90_put_att(ncidout_lbl, varid, 'long_name', 'reference weight'))
    END IF

    IF (l_self_broadening) THEN
      CALL nf(nf90_def_var(ncidout_lbl, 'kabs', NF90_FLOAT, &
        (/dimid1,dimid2,dimid3/), varid))
    ELSE
      CALL nf(nf90_def_var(ncidout_lbl, 'kabs', NF90_FLOAT, &
        (/dimid1,dimid2/), varid))
    END IF
    CALL nf(nf90_put_att(ncidout_lbl, varid, 'title', 'absorption' ))
    CALL nf(nf90_put_att(ncidout_lbl, varid, 'long_name', 'absorption' ))
    IF (l_fit_cont_data) THEN
      CALL nf(nf90_put_att(ncidout_lbl, varid, 'units', 'm5 kg-2' ))
    ELSE
      CALL nf(nf90_put_att(ncidout_lbl, varid, 'units', 'm2 kg-1' ))
    END IF

    CALL nf(nf90_enddef(ncidout_lbl))

!   Initialise number of frequencies written to file
    n_nu_written=0

  END SUBROUTINE output_lbl_band_cdf_init


  SUBROUTINE output_lbl_band_cdf

    USE netcdf
    IMPLICIT NONE
    INTEGER :: varid                   ! variable ID

!   Write frequency array and absorption coefficient
    CALL nf(nf90_inq_varid(ncidout_lbl,'nu',varid))
    CALL nf(nf90_put_var(ncidout_lbl,varid,nu_wgt, &
      start=(/n_nu_written+1/), count=(/n_nu/)))
    IF (l_output_reference_weight) THEN
      CALL nf(nf90_inq_varid(ncidout_lbl,'wgt',varid))
      CALL nf(nf90_put_var(ncidout_lbl,varid,wgt_ref, &
        start=(/n_nu_written+1/), count=(/n_nu/)))
    END IF
    CALL nf(nf90_inq_varid(ncidout_lbl,'kabs',varid))
    IF (l_self_broadening) THEN
      CALL nf(nf90_put_var(ncidout_lbl,varid,kabs_all_sb, &
        start=(/n_nu_written+1, 1, 1/), &
        count=(/n_nu ,n_pt_pair, n_gas_frac/)))
    ELSE
      CALL nf(nf90_put_var(ncidout_lbl,varid,kabs_all, &
        start=(/n_nu_written+1, 1/), count=(/n_nu,n_pt_pair/)))
    END IF

    n_nu_written=n_nu_written+n_nu

  END SUBROUTINE output_lbl_band_cdf


  SUBROUTINE output_map_band_cdf_init

    USE netcdf
    IMPLICIT NONE
    INTEGER :: dimid1, dimid2, dimid3, dimid4 ! dimension ID
    INTEGER :: varid                          ! variable ID
    INTEGER :: n_nu_band(n_selected_band)     ! number of frequency points
    LOGICAL :: l_map_exist                    ! flag for mapping file existing 

!   Calculate total number of frequency points
    DO ibb=1, n_selected_band
      ib=list_band(ibb)
      band_width = band_max(ib) - band_min(ib)
      DO jx = 1, n_band_exclude(ib)
        band_width = band_width + &
          band_min(index_exclude(jx, ib)) - band_max(index_exclude(jx, ib))
      ENDDO
      n_nu_band(ibb)=NINT(band_width/nu_inc)
    END DO

    INQUIRE(FILE=file_map, EXIST=l_map_exist)

!   Create or open for writing
    IF (l_map_exist) THEN
      CALL nf(nf90_open(TRIM(file_map),NF90_WRITE,ncidout_map))

    ELSE
      CALL nf(nf90_create(TRIM(file_map), &
                          IOR(NF90_NOCLOBBER, NF90_64BIT_OFFSET), &
                          ncidout_map))

!     Create dimensions
      CALL nf(nf90_def_dim(ncidout_map, 'nu', MAXVAL(n_nu_band), dimid1))
      CALL nf(nf90_def_dim(ncidout_map, 'nd_k_term', nd_k_term+1, dimid2))
      CALL nf(nf90_def_dim(ncidout_map, 'pre', n_pt_pair, dimid3))
      CALL nf(nf90_def_dim(ncidout_map, 'band', n_band, dimid4))

      CALL nf(nf90_def_var(ncidout_map, 'nu', NF90_DOUBLE, &
        (/ dimid1,dimid4 /), varid))
      CALL nf(nf90_put_att(ncidout_map, varid, 'title', 'wavenumber'))
      CALL nf(nf90_put_att(ncidout_map, varid, 'long_name', 'wavenumber'))
      CALL nf(nf90_put_att(ncidout_map, varid, 'units', 'm-1'))
      CALL nf(nf90_put_att(ncidout_map, varid, 'step', nu_inc))

      CALL nf(nf90_def_var(ncidout_map, 'n_k', NF90_INT, (/ dimid4 /), varid))
      CALL nf(nf90_put_att(ncidout_map, varid, 'title', 'number of k-terms'))
      CALL nf(nf90_put_att(ncidout_map, varid, 'long_name', &
        'number of k-terms'))

      CALL nf(nf90_def_var(ncidout_map, 'n_pre', NF90_INT, (/ dimid4 /), varid))
      CALL nf(nf90_put_att(ncidout_map, varid, 'title', 'number of pbins'))
      CALL nf(nf90_put_att(ncidout_map, varid, 'long_name', &
        'number of pressure bins'))

      CALL nf(nf90_def_var(ncidout_map, 'ig', NF90_INT, &
        (/ dimid2,dimid4 /), varid))
      CALL nf(nf90_put_att(ncidout_map, varid, 'title', 'g-points'))
      CALL nf(nf90_put_att(ncidout_map, varid, 'long_name', &
        'g-quadrature points'))

      CALL nf(nf90_def_var(ncidout_map, 'w_k', NF90_DOUBLE, &
        (/dimid2,dimid4/), varid))
      CALL nf(nf90_put_att(ncidout_map, varid, 'title', 'weights' ))
      CALL nf(nf90_put_att(ncidout_map, varid, 'long_name', 'weights' ))

      CALL nf(nf90_def_var(ncidout_map, 'pstart', NF90_INT, &
        (/ dimid3,dimid4 /), varid))
      CALL nf(nf90_put_att(ncidout_map, varid, 'title', 'start of pbin'))
      CALL nf(nf90_put_att(ncidout_map, varid, 'long_name', &
        'start of pressure bin'))

      CALL nf(nf90_def_var(ncidout_map, 'pend', NF90_INT, &
        (/ dimid3,dimid4 /), varid))
      CALL nf(nf90_put_att(ncidout_map, varid, 'title', 'end of pbin'))
      CALL nf(nf90_put_att(ncidout_map, varid, 'long_name', &
        'end of pressure bin'))

      CALL nf(nf90_def_var(ncidout_map, 'pcount', NF90_INT, &
        (/ dimid3,dimid4 /), varid))
      CALL nf(nf90_put_att(ncidout_map, varid, 'title', 'size of pbin'))
      CALL nf(nf90_put_att(ncidout_map, varid, 'long_name', &
        'number of points in pressure bin'))

      CALL nf(nf90_def_var(ncidout_map, 'pmap', NF90_INT, &
        (/dimid1,dimid4/), varid))
      CALL nf(nf90_put_att(ncidout_map, varid, 'title', &
        'pressure bin mapping' ))
      CALL nf(nf90_put_att(ncidout_map, varid, 'long_name', &
        'pressure bin mapping' ))

      CALL nf(nf90_def_var(ncidout_map, 'map', NF90_INT, &
        (/dimid1,dimid4/), varid))
      CALL nf(nf90_put_att(ncidout_map, varid, 'title', 'mapping' ))
      CALL nf(nf90_put_att(ncidout_map, varid, 'long_name', 'mapping' ))

      CALL nf(nf90_enddef(ncidout_map))
    END IF

  END SUBROUTINE output_map_band_cdf_init


  SUBROUTINE output_map_band_cdf

    USE netcdf
    IMPLICIT NONE
    INTEGER :: varid                   ! variable ID

!   Write frequency array mapping arrays
    CALL nf(nf90_inq_varid(ncidout_map,'nu',varid))
    CALL nf(nf90_put_var(ncidout_map,varid,nu_wgt, &
      start=(/1,ib/), count=(/n_nu,1/)))
    CALL nf(nf90_inq_varid(ncidout_map,'n_k',varid))
    CALL nf(nf90_put_var(ncidout_map,varid,(/n_k(ib)/), &
      start=(/ib/), count=(/1/)))
    CALL nf(nf90_inq_varid(ncidout_map,'n_pre',varid))
    CALL nf(nf90_put_var(ncidout_map,varid,(/n_pre/), &
      start=(/ib/), count=(/1/)))
    CALL nf(nf90_inq_varid(ncidout_map,'ig',varid))
    CALL nf(nf90_put_var(ncidout_map,varid,ig, &
      start=(/1,ib/), count=(/n_k(ib)+1,1/)))
    CALL nf(nf90_inq_varid(ncidout_map,'w_k',varid))
    CALL nf(nf90_put_var(ncidout_map,varid,w_k(1:n_k(ib),ib), &
      start=(/1,ib/), count=(/n_k(ib),1/)))
    CALL nf(nf90_inq_varid(ncidout_map,'pstart',varid))
    CALL nf(nf90_put_var(ncidout_map,varid,pstart(1:n_pre), &
      start=(/1,ib/), count=(/n_pre,1/)))
    CALL nf(nf90_inq_varid(ncidout_map,'pend',varid))
    CALL nf(nf90_put_var(ncidout_map,varid,pend(1:n_pre), &
      start=(/1,ib/), count=(/n_pre,1/)))
    CALL nf(nf90_inq_varid(ncidout_map,'pcount',varid))
    CALL nf(nf90_put_var(ncidout_map,varid,pcount(1:n_pre), &
      start=(/1,ib/), count=(/n_pre,1/)))
    CALL nf(nf90_inq_varid(ncidout_map,'pmap',varid))
    CALL nf(nf90_put_var(ncidout_map,varid,pmap, &
      start=(/1,ib/), count=(/n_nu,1/)))
    CALL nf(nf90_inq_varid(ncidout_map,'map',varid))
    CALL nf(nf90_put_var(ncidout_map,varid,map, &
      start=(/1,ib/), count=(/n_nu,1/)))

  END SUBROUTINE output_map_band_cdf


  SUBROUTINE output_ck_cdf

    USE netcdf
    IMPLICIT NONE
    Integer :: ncid                            ! netCDF file ID
    Integer :: dimid1, dimid2, dimid3, dimid4  ! dimension ID
    Integer :: varid                           ! variable ID

!   Create the file and open for writing
    Call nf(nf90_create(Trim(file_k)//'.nc',NF90_NOCLOBBER,ncid))

!   Write dimensions
    Call nf(nf90_def_dim(ncid, 'nd_k_term', nd_k_term, dimid1))
    Call nf(nf90_def_dim(ncid, 'pt_pair', n_pt_pair, dimid2))
    IF (l_self_broadening) &
      CALL nf(nf90_def_dim(ncid, 'gas_frac', n_gas_frac, dimid3))
    Call nf(nf90_def_dim(ncid, 'band', n_selected_band, dimid4))

    Call nf(nf90_def_var(ncid, 'band', NF90_INT, dimid4, varid))
    Call nf(nf90_put_att(ncid, varid, 'title', 'band number' ))
    Call nf(nf90_enddef(ncid))
    Call nf(nf90_put_var(ncid, varid, list_band(1:n_selected_band) ))
    Call nf(nf90_redef(ncid))

!   Write variables
    Call nf(nf90_def_var(ncid, 'n_k', NF90_INT, dimid4, varid))
    Call nf(nf90_put_att(ncid, varid, 'title', 'number of k-terms' ))
    Call nf(nf90_enddef(ncid))
    Call nf(nf90_put_var(ncid, varid, n_k(list_band(1:n_selected_band)) ))
    Call nf(nf90_redef(ncid))

    IF (l_self_broadening) THEN
      Call nf(nf90_def_var(ncid, 'gas_frac', NF90_FLOAT, dimid3, varid))
      Call nf(nf90_put_att(ncid, varid, 'title', 'gas fraction'))
      Call nf(nf90_put_att(ncid, varid, 'long_name', 'gas fraction'))
      Call nf(nf90_enddef(ncid))
      Call nf(nf90_put_var(ncid, varid, gas_frac(1:n_gas_frac) ))
      Call nf(nf90_redef(ncid))
    END IF

    Call nf(nf90_def_var(ncid, 'p_calc', NF90_FLOAT, dimid2, varid))
    Call nf(nf90_put_att(ncid, varid, 'title', 'pressure' ))
    Call nf(nf90_put_att(ncid, varid, 'long_name', 'pressure' ))
    Call nf(nf90_put_att(ncid, varid, 'units', 'Pa' ))
    Call nf(nf90_enddef(ncid))
    Call nf(nf90_put_var(ncid, varid, p_calc(1:n_pt_pair) ))
    Call nf(nf90_redef(ncid))

    Call nf(nf90_def_var(ncid, 't_calc', NF90_FLOAT, dimid2, varid))
    Call nf(nf90_put_att(ncid, varid, 'title', 'temperature' ))
    Call nf(nf90_put_att(ncid, varid, 'long_name', 'temperature' ))
    Call nf(nf90_put_att(ncid, varid, 'units', 'K' ))
    Call nf(nf90_enddef(ncid))
    Call nf(nf90_put_var(ncid, varid, t_calc(1:n_pt_pair) ))
    Call nf(nf90_redef(ncid))

    Call nf(nf90_def_var(ncid, 'w_k', NF90_FLOAT, &
      (/dimid1,dimid4/), varid))
    Call nf(nf90_put_att(ncid, varid, 'title', 'weights' ))
    Call nf(nf90_put_att(ncid, varid, 'long_name', 'weights' ))
    Call nf(nf90_enddef(ncid))
    Call nf(nf90_put_var(ncid, varid, w_k(:,list_band(1:n_selected_band)) ))

    Call nf(nf90_redef(ncid))
    IF (l_self_broadening) THEN
      Call nf(nf90_def_var(ncid, 'kopt', NF90_FLOAT, &
        (/dimid1,dimid2,dimid3,dimid4/), varid))
      Call nf(nf90_enddef(ncid))
      Call nf(nf90_put_var(ncid, varid, kopt_all_sb ))
    ELSE
      Call nf(nf90_def_var(ncid, 'kopt', NF90_FLOAT, &
        (/dimid1,dimid2,dimid4/), varid))
      Call nf(nf90_enddef(ncid))
      Call nf(nf90_put_var(ncid, varid, kopt_all ))
    END IF
    Call nf(nf90_redef(ncid))
    Call nf(nf90_put_att(ncid, varid, 'title', 'k-term' ))
    Call nf(nf90_put_att(ncid, varid, 'long_name', 'k-term' ))
    IF (l_fit_cont_data) THEN
      Call nf(nf90_put_att(ncid, varid, 'units', 'm5 kg-2' ))
    ELSE
      Call nf(nf90_put_att(ncid, varid, 'units', 'm2 kg-1' ))
    END IF

    Call nf(nf90_close(ncid))

  END SUBROUTINE output_ck_cdf


  SUBROUTINE close_lbl_files

    use netcdf
    IMPLICIT NONE
    INTEGER :: len_base

!   Close files
    IF (l_lbl_exist) THEN
      CALL nf(nf90_close(ncidin_lbl))
    ELSE
      CALL nf(nf90_close(ncidout_lbl))
      ! Write empty file to indicate completion
      len_base=SCAN(file_lbl, ".", .TRUE.)-1
      IF (len_base < 0) len_base=LEN(TRIM(file_lbl))
      CALL nf(nf90_create(file_lbl(1:len_base)//'.done', &
              NF90_CLOBBER, ncidout_lbl))
      CALL nf(nf90_close(ncidout_lbl))
    ENDIF

  END SUBROUTINE close_lbl_files

  SUBROUTINE close_map_files
  
    use netcdf
    IMPLICIT NONE
  
!   Close files
    IF (l_load_map) THEN
      CALL nf(nf90_close(ncidin_map))
    ELSE IF (l_save_map) THEN
      CALL nf(nf90_close(ncidout_map))
    ENDIF
      
  END SUBROUTINE close_map_files

  Subroutine nf(status)
    USE netcdf
    Integer, Intent(IN):: status
    If (status /= NF90_NOERR) Then
       Write(*,*) 'netCDF-ERROR: ',nf90_strerror(status)
       Stop 'STOPPED!'
    End If
  End Subroutine nf


END SUBROUTINE corr_k_single
