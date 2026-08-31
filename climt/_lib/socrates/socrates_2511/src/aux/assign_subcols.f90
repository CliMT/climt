! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to calculate number of additional sub-columns to assign to 
!  each k-term when McICA is used
!
! Method:
!       Described in Hill et al, 2011 (DOI: 10.1002/qj.732)
!
!- ---------------------------------------------------------------------
PROGRAM assign_subcols

! Modules to set types of variables:
  USE def_spectrum
  USE def_mcica, ONLY: StrMcica, read_mcica_data
  USE realtype_rd
  USE rad_pcf

  IMPLICIT NONE

! Declaration of variables.

LOGICAL :: l_interactive
!       Switch for interactive use

! Spectral information:
  CHARACTER  (LEN=256) :: file_spectral_sw
!       Name of SW spectral file
  CHARACTER  (LEN=256) :: file_spectral_lw
!       Name of LW spectral file
  CHARACTER  (LEN=256) :: file_mcica_data = 'mcica_data'
!       Name of output MCICA data file
  CHARACTER  (LEN=256) :: file_mcica_xcw
!       Name of input MCICA data file for xcw array
  TYPE  (StrSpecData) :: Spectrum_sw
!       SW Spectral data
  TYPE  (StrSpecData) :: Spectrum_lw
!       LW Spectral data
  TYPE  (StrMcica) :: mcica_data
!       MCICA data
  INTEGER :: i_gas
!       Index of main gas
  REAL (RealK) :: ext_gas
!       extinction due to major gas
  REAL (RealK), ALLOCATABLE :: planck_flux(:)
!       mean Planck flux in each band

! Algorithm specfics
  INTEGER :: tot_subcol_gen
!       Number of sub-columns to generate
  INTEGER :: ioverlap = 3
!       Overlap method hardwired for legacy, no longer used from file
  INTEGER :: n_subcol_sw
!       No additional sub-columns in SW
  INTEGER :: n_subcol_lw
!       No additional sub-columns in LW

! Cloud optical properties
  REAL (RealK) :: k_ext_tot_cloud
!       total cloud extinction
  REAL  (RealK) :: absorb_cloud
!       cloud asborption
  REAL  (RealK), ALLOCATABLE :: cloud_parameter(:, :)
!       parameters in cloud optical properties parametrization

! Calculated 'importance'
  REAL (RealK), ALLOCATABLE :: iota_absorb_sw(:, :)
!       importance of k-term for absorption
  REAL  (RealK), ALLOCATABLE :: iota_extinct_sw(:, :)
!       importance of k-term for extinction
  REAL (RealK), ALLOCATABLE :: iota_absorb_lw(:, :)
!       importance of k-term for absorption
  REAL  (RealK), ALLOCATABLE :: iota_extinct_lw(:, :)
!       importance of k-term for extinction

  INTEGER, Parameter :: drop_param = 5
!       liquid water content doplet parametrization

  INTEGER, Parameter :: n_temperature = 200
!       no of temperatures used in calculation of mean Planckian flux

  REAL :: temperature
!       temperature for calculation of Planck function

  INTEGER :: i, j
!       Loop variables

  INTEGER :: ierr
!       Error flag

  INTEGER :: n_add_lw=6
!       Number of additional sub-columns in the LW

  INTEGER :: n_add_sw=10
!       Number of additional sub-columns in the SW

  INTEGER :: order_lw(99999), order_lw_single(99999)
!       Re-ordering of sub-columns required to match importance

  INTEGER :: total_k_lw
!       Total number of LW k-terms.

  INTEGER, ALLOCATABLE :: subcols_per_k_lw(:,:)
!       Number of sub-columns assigned to each LW k-term

  INTEGER, ALLOCATABLE :: subcols_per_k_sw(:,:)
!       Number of sub-columns assigned to each SW k-term

  LOGICAL :: l_exist
!       Flag for pre-existence of mcica_data file

! External functions:
  LOGICAL, EXTERNAL :: set_interactive
!       Function to set the flag for interactive operation

  INTEGER :: iu_mcd = 80

! Set the flag for interactive operation
  l_interactive=set_interactive()

! ------------------------------------------------------------------
! Spectral Data:
! ------------------------------------------------------------------
! Read in the SW spectral file.
  WRITE(*, "(a)") "Enter the name of the SW spectral file."
  READ(*, "(a)") file_spectral_sw
  CALL read_spectrum(file_spectral_sw, Spectrum_sw)
  WRITE(*, "(a)") &
    "Number of extra sub-columns in the SW for optimal sampling:"
  READ(*, *) n_add_sw

! Do SW importance calculation
  ALLOCATE(iota_absorb_sw(Spectrum_sw%Dim%nd_band                       &
                          , Spectrum_sw%Dim%nd_k_term))
  ALLOCATE(iota_extinct_sw(Spectrum_sw%Dim%nd_band                      &
                           , Spectrum_sw%Dim%nd_k_term))
  ALLOCATE(cloud_parameter(Spectrum_sw%Dim%nd_cloud_parameter           &
                           , Spectrum_sw%Basic%n_band))
  cloud_parameter(:, :)=Spectrum_sw%Drop%parm_list(:, :, drop_param)

  CALL k_term_importance(Spectrum_sw%Basic%n_band                       &
  , Spectrum_sw%Dim%nd_cloud_parameter                                  &
  , Spectrum_sw%Dim%nd_species                                          &
  , Spectrum_sw%Dim%nd_k_term                                           &
  , Spectrum_sw%Gas%index_absorb                                        &
  , cloud_parameter                                                     &
  , Spectrum_sw%Gas%i_band_k                                            &
  , Spectrum_sw%Solar%solar_flux_band                                   &
  , Spectrum_sw%Gas%w                                                   &
  , Spectrum_sw%Gas%k                                                   &
  , iota_extinct_sw                                                     &
  , iota_absorb_sw)

  IF (ALLOCATED(cloud_parameter)) DEALLOCATE(cloud_parameter)

  ALLOCATE(subcols_per_k_sw(Spectrum_sw%Dim%nd_band                     &
                            , Spectrum_sw%Dim%nd_k_term))
  n_add_sw = n_add_sw/2
  subcols_per_k_sw = 1

! Read LW spectral file
  WRITE(*, "(a)") "Enter the name of the LW spectral file."
  READ(*, "(a)") file_spectral_lw
  CALL read_spectrum(file_spectral_lw, Spectrum_lw)
  WRITE(*, "(a)") &
    "Number of extra sub-columns in the LW for optimal sampling:"
  READ(*, *) n_add_lw

! Calculate LW weights, assuming a Planck fn.
  ALLOCATE(planck_flux(Spectrum_lw%Dim%nd_band))
  planck_flux = 0.0

  DO i = 1, Spectrum_lw%Basic%n_band
    DO j = 1, n_temperature
      temperature = (j+129.0)/Spectrum_lw%Planck%t_ref_planck
      planck_flux(i) = planck_flux(i)                                   &
                    +Spectrum_lw%Planck%thermal_coeff(0,i)              &
                    +(Spectrum_lw%Planck%thermal_coeff(1,i)             &
                      *temperature)                                     &
                    +(Spectrum_lw%Planck%thermal_coeff(2,i)             &
                      *temperature**2)                                  &
                    +(Spectrum_lw%Planck%thermal_coeff(3,i)             &
                      *temperature**3)                                  &
                    +(Spectrum_lw%Planck%thermal_coeff(4,i)             &
                      *temperature**4)
    ENDDO
  END DO
  planck_flux(:) = planck_flux(:)/SUM(planck_flux)

! Calculate total number of LW k-terms for purposes of printing output
  total_k_lw=0
  DO i = 1, Spectrum_lw%Basic%n_band
    i_gas = Spectrum_lw%Gas%index_absorb(1, i)
    DO j=1, Spectrum_lw%Gas%i_band_k(i,i_gas)
      total_k_lw=total_k_lw+1
    ENDDO
  ENDDO

! Do importance calculation in the LW
  ALLOCATE(cloud_parameter(Spectrum_lw%Dim%nd_cloud_parameter           &
                           , Spectrum_lw%Basic%n_band))
  cloud_parameter(:, :) = Spectrum_lw%Drop%parm_list(:, :, drop_param)
  ALLOCATE(iota_absorb_lw(Spectrum_lw%Dim%nd_band                       &
                          , Spectrum_lw%Dim%nd_k_term))
  ALLOCATE(iota_extinct_lw(Spectrum_lw%Dim%nd_band                      &
                           , Spectrum_lw%Dim%nd_k_term))

  CALL k_term_importance(Spectrum_lw%Basic%n_band                       &
    , Spectrum_lw%Dim%nd_cloud_parameter                                &
    , Spectrum_lw%Dim%nd_species                                        &
    , Spectrum_lw%Dim%nd_k_term                                         &
    , Spectrum_lw%Gas%index_absorb                                      &
    , cloud_parameter                                                   &
    , Spectrum_lw%Gas%i_band_k                                          &
    , planck_flux                                                       &
    , Spectrum_lw%Gas%w                                                 &
    , Spectrum_lw%Gas%k                                                 &
    , iota_extinct_lw                                                   &
    , iota_absorb_lw)

  IF (ALLOCATED(cloud_parameter)) DEALLOCATE(cloud_parameter)

  ALLOCATE(subcols_per_k_lw(Spectrum_lw%Dim%nd_band                     &
                            , Spectrum_lw%Dim%nd_k_term))
  subcols_per_k_lw = 1

! Match SW and LW sub-columns by the absoroption importance of the 
! k-terms (before additional sub-columns are distributed to give the 
! single sampling reordering)
  CALL match_subcols(Spectrum_sw%Basic%n_band                           &
    , Spectrum_lw%Basic%n_band                                          &
    , Spectrum_sw%Dim%nd_species                                        &
    , Spectrum_lw%Dim%nd_species                                        &
    , Spectrum_sw%Dim%nd_k_term                                         &
    , Spectrum_lw%Dim%nd_k_term                                         &
    , Spectrum_sw%Gas%index_absorb                                      &
    , Spectrum_lw%Gas%index_absorb                                      &
    , Spectrum_sw%Gas%i_band_k                                          &
    , Spectrum_lw%Gas%i_band_k                                          &
    , iota_absorb_sw                                                    &
    , iota_absorb_lw                                                    &
    , subcols_per_k_sw                                                  &
    , subcols_per_k_lw                                                  &
    , order_lw_single)

! Assign subcolumns to k-terms depending on 'importance' for the SW 
! there are two calls to distribute_subcols, one for minimising one 
! where extinction importance is minimised, the other where absorption 
! importance is minimised.
  CALL distribute_subcols(Spectrum_sw%Basic%n_band                      &
    , Spectrum_sw%Dim%nd_k_term                                         &
    , n_add_sw                                                          &
    , iota_extinct_sw                                                   &
    , subcols_per_k_sw)

  CALL distribute_subcols(Spectrum_sw%Basic%n_band                      &
    , Spectrum_sw%Dim%nd_k_term                                         &
    , n_add_sw                                                          &
    , iota_absorb_sw                                                    &
    , subcols_per_k_sw)


! Assign subcolumns to k-terms depending on 'importance'
  CALL distribute_subcols(Spectrum_lw%Basic%n_band                      &
    , Spectrum_lw%Dim%nd_k_term                                         &
    , n_add_lw                                                          &
    , iota_absorb_lw                                                    &
    , subcols_per_k_lw)

! Match SW and LW sub-columns by the absoroption importance of the 
! k-terms
  CALL match_subcols(Spectrum_sw%Basic%n_band                           &
    , Spectrum_lw%Basic%n_band                                          &
    , Spectrum_sw%Dim%nd_species                                        &
    , Spectrum_lw%Dim%nd_species                                        &
    , Spectrum_sw%Dim%nd_k_term                                         &
    , Spectrum_lw%Dim%nd_k_term                                         &
    , Spectrum_sw%Gas%index_absorb                                      &
    , Spectrum_lw%Gas%index_absorb                                      &
    , Spectrum_sw%Gas%i_band_k                                          &
    , Spectrum_lw%Gas%i_band_k                                          &
    , iota_absorb_sw                                                    &
    , iota_absorb_lw                                                    &
    , subcols_per_k_sw                                                  &
    , subcols_per_k_lw                                                  &
    , order_lw)

! Read XCW array from an existing mcica_data file
  WRITE(*, "(a)") "Enter the name of a file containing the XCW array:"
  READ(*, "(a)") file_mcica_xcw
  CALL read_mcica_data(mcica_data, file_mcica_xcw)

! Read in the total number of generated sub-columns
  WRITE(*, "(a)") &
    "Total number of sub-columns to be generated"
  WRITE(*, "(a)") &
    "(or -1 to generate enough for single sampling, -2 for optimal sampling):"
  READ(*, *) tot_subcol_gen
  IF (tot_subcol_gen == -1) tot_subcol_gen = total_k_lw
  IF (tot_subcol_gen == -2) tot_subcol_gen = total_k_lw+n_add_lw

! Output data as a mcica_data file.
  INQUIRE(file=TRIM(file_mcica_data), exist=l_exist)
  IF (l_exist) THEN
     WRITE(*, '(3a)') &
          'Error: The file "',TRIM(file_mcica_data), &
          '" already exists: it will not be overwritten.'
  ELSE
    CALL get_free_unit(ierr, iu_mcd)
    OPEN(UNIT=iu_mcd, FILE=file_mcica_data, IOSTAT=ierr, STATUS='UNKNOWN')
    IF (ierr /= 0) THEN
     WRITE(*, '(3a)') &
          'Error: The file "',TRIM(file_mcica_data), &
          '" cannot be opened.'
      STOP
    ELSE
      WRITE(iu_mcd,'(a)') '+ ----'
      WRITE(iu_mcd,'(a)') 'Data file for the generation of sub-grid cloud.'
      WRITE(iu_mcd,'(4a)') 'Applicable for spectral files: ', &
        TRIM(file_spectral_sw(SCAN(file_spectral_sw,'/',.TRUE.)+1:)), ', ', &
        TRIM(file_spectral_lw(SCAN(file_spectral_lw,'/',.TRUE.)+1:))
      WRITE(iu_mcd,'(a/,a/)') '- ----','*DATA'
      WRITE(iu_mcd,'(a/,i5)') 'tot_subcol_gen', tot_subcol_gen
      WRITE(iu_mcd,'(a/,i5)') 'subcol_need_single', total_k_lw
      WRITE(iu_mcd,'(a/,i5)') 'subcol_need_optimal', total_k_lw+n_add_lw
      WRITE(iu_mcd,'(a/,i5)') 'ipph', mcica_data%ipph
      WRITE(iu_mcd,'(a/,i5,/)') 'ioverlap', ioverlap
      WRITE(iu_mcd,'(A24)') "lw_subcol_reorder_single"
      WRITE(iu_mcd,'(10i6)') order_lw_single(1:total_k_lw)
      WRITE(iu_mcd,'(/A11)') "sw_subcol_k"
      DO i=1, Spectrum_sw%Basic%n_band
        i_gas = Spectrum_sw%Gas%index_absorb(1, i)
        DO j=1, Spectrum_sw%Gas%i_band_k(i,i_gas)
          IF (subcols_per_k_sw(i,j) > 1) THEN 
            WRITE(iu_mcd,'(3i4)') i,j,subcols_per_k_sw(i,j)
          ENDIF
        ENDDO
      ENDDO
      WRITE(iu_mcd,'(A12/)') ' -99   0   0'
      WRITE(iu_mcd,'(A11)') "lw_subcol_k"
      DO i=1, Spectrum_lw%Basic%n_band
        i_gas = Spectrum_lw%Gas%index_absorb(1, i)
        DO j=1, Spectrum_lw%Gas%i_band_k(i,i_gas)
          IF (subcols_per_k_lw(i,j) > 1) THEN 
            WRITE(iu_mcd,'(3i4)') i,j,subcols_per_k_lw(i,j)
          ENDIF
        ENDDO
      ENDDO
      WRITE(iu_mcd,'(A12/)') ' -99   0   0'
      WRITE(iu_mcd,'(A25)') "lw_subcol_reorder_optimal"
      WRITE(iu_mcd,'(10i6)') order_lw(1:total_k_lw+n_add_lw)
      WRITE(iu_mcd,'(2(/a/,i6))') 'n1', mcica_data%n1, 'n2', mcica_data%n2
      WRITE(iu_mcd,'(a/,(5e15.8))') 'xcw', mcica_data%xcw
      WRITE(iu_mcd,'(a)') '*END'
      WRITE(*,'(3a)') "Written output to '",TRIM(file_mcica_data),"'" 
    ENDIF
    CLOSE(iu_mcd)
  END IF

  IF (ALLOCATED(subcols_per_k_lw)) DEALLOCATE(subcols_per_k_lw)
  IF (ALLOCATED(iota_extinct_lw)) DEALLOCATE(iota_extinct_lw)
  IF (ALLOCATED(iota_absorb_lw)) DEALLOCATE(iota_absorb_lw)
  IF (ALLOCATED(iota_extinct_sw)) DEALLOCATE(iota_extinct_sw)
  IF (ALLOCATED(iota_absorb_sw)) DEALLOCATE(iota_absorb_sw)

CONTAINS

  SUBROUTINE k_term_importance(n_band                                    &
    , n_cloud_parameter                                                  &
    , n_species                                                          &
    , n_k_term                                                           &
    , index_absorb                                                       &
    , cloud_parameter                                                    &
    , i_band_k                                                           &
    , weight_band                                                        &
    , w_esft                                                             &
    , k_esft                                                             &
    , iota_extinct                                                       &
    , iota_absorb)

    IMPLICIT NONE

!   Input variables
    INTEGER, Intent(IN) :: n_band
!         Number of bands

    INTEGER, Intent(IN) :: n_cloud_parameter
!         Size of cloud parameters array 
!         (depends on droplet optical properties parametrization)

    INTEGER, Intent(IN) :: n_species
!         Number of gases in each band

    INTEGER, Intent(IN) :: n_k_term
!         Size allocated for k-terms in each bands

    INTEGER, Intent(IN) :: index_absorb(n_species, n_band)
!         indeces of the gases in each band

    REAL (RealK), Intent(IN) :: cloud_parameter(n_cloud_parameter, n_band)
!         parameters for cloud droplet opltical properties parametrization

    INTEGER, Intent(IN) :: i_band_k(n_band, n_species)
!         Number of k-terms in each band

    REAL (RealK), Intent(IN) :: weight_band(n_band)
!         Weight of each band (from fraction of total solar flux in SW and
!         fraction of mean Planckian in LW)

    REAL (RealK), Intent(IN) :: w_esft(n_k_term, n_band, n_species)
!         Weight of each k-term as fraction of the band

    REAL (RealK), Intent(IN) :: k_esft(n_k_term, n_band, n_species)
!         Exponential esft terms

!   Output importances
    REAL (RealK), Intent(OUT) :: iota_extinct(n_band, n_k_term)
!         Importance calculated from extinction

    REAL (RealK), Intent(OUT) :: iota_absorb(n_band, n_k_term)
!         Importance calculated from absorption

!   Local variables
    REAL (RealK) :: weight_incr
!         Weight of each k-term as fraction of entire spetrum

    INTEGER :: i, j, i_gas
!         Loop variables

!   Local rather arbitrary Constants
    REAL (RealK), Parameter :: cloud_mass_thick = 0.203147793509528
!         Mass of optically thick cloud
    REAL (RealK), Parameter :: cloud_mass_thin = 0.00203147793509528
!         Mass of optically thin cloud
    REAL (RealK), Parameter :: h2o_mass = 25.0
!         Mass of H2O
    REAL (RealK), Parameter :: co2_mass = 5.0
!         Mass of CO2
    REAL (RealK), Parameter :: o3_mass = 0.008
!         Mass of O3
    REAL (RealK), Parameter :: radius_effect = 0.00001
!         effective cloud droplet radius

!   Initialize the calculated iota variables
    iota_extinct = 0.0
    iota_absorb = 0.0

    DO i = 1, n_band
!     Calculate cloud optical properties (constant in each band)
      i_gas = index_absorb(1, i)
      k_ext_tot_cloud = (cloud_parameter(1, i)+radius_effect            &
            *(cloud_parameter(2, i)+radius_effect                       &
            *cloud_parameter(3, i)))                                    &
            /(1.0e+00+radius_effect                                     &
            *(cloud_parameter(4, i)+radius_effect                       &
            *(cloud_parameter(5, i)+radius_effect                       &
            *cloud_parameter(6, i))))    
      absorb_cloud = k_ext_tot_cloud                                    &
            *(cloud_parameter(7, i)+radius_effect                       &
            *(cloud_parameter(8, i)+radius_effect                       &
            *cloud_parameter(9, i)))                                    &
            /(1.0e+00+radius_effect                                     &
            *(cloud_parameter(10, i)+radius_effect                      &
            *cloud_parameter(11, i)))

      DO j = 1, i_band_k(i, i_gas)
!       Calculate non-grey properties
        weight_incr = weight_band(i)*w_esft(j, i, i_gas)
        IF (i_gas == 1) ext_gas = exp(-(k_esft(j, i, i_gas)*h2o_mass))
        IF (i_gas == 2) ext_gas = exp(-(k_esft(j, i, i_gas)*co2_mass))
        IF (i_gas == 3) ext_gas = exp(-(k_esft(j, i, i_gas)*o3_mass))

!       Calculate iotas
        iota_absorb(i, j) = weight_incr*ext_gas                         &
                          *ABS(EXP(-absorb_cloud*cloud_mass_thick)      &
                               -EXP(-absorb_cloud*cloud_mass_thin))
        iota_extinct(i, j) = weight_incr*ext_gas                        &
                           *ABS(EXP(-k_ext_tot_cloud*cloud_mass_thick)  &
                                -EXP(-k_ext_tot_cloud*cloud_mass_thin))
      END DO
    END DO
  
  END SUBROUTINE


  SUBROUTINE distribute_subcols(n_band                                  &
    , n_k_term                                                          &
    , n_add                                                             &
    , iota                                                              &
    , subcols_per_k)

    IMPLICIT NONE

!   Input variables
    INTEGER, Intent(IN) :: n_band
!         Number of bands

    INTEGER, Intent(IN) :: n_k_term
!         Size allocated for k-terms in each bands

    INTEGER, Intent(IN) :: n_add
!         Number of additional sub-columns to assign

!   Modified variables
    REAL (RealK), Intent(INOUT) :: iota(n_band, n_k_term)
!         Importance of each k-term (importance/no sub-columns)

    INTEGER, Intent(INOUT) :: subcols_per_k(n_band, n_k_term)
!         Number of sub-columns assigned to each k-term


!   Local variables
    INTEGER :: ind(2)
!         indeces of most importantk-term

    REAL (RealK) :: iota_temp(n_band, n_k_term)
!         temporary variables used to hold importance


!   Assign additional sub-column to most important k-term, recalculating 
!   importance for each additional sub-column
    DO i = 1, n_add
      iota_temp = iota/subcols_per_k
      ind = MAXLOC(iota_temp)
      subcols_per_k(ind(1),ind(2)) = subcols_per_k(ind(1),ind(2))+1
    END DO
    iota = iota/subcols_per_k
    
  END SUBROUTINE


  SUBROUTINE match_subcols(n_band_sw                                    &
    , n_band_lw                                                         &
    , n_species_sw                                                      &
    , n_species_lw                                                      &
    , n_k_term_sw                                                       &
    , n_k_term_lw                                                       &
    , index_absorb_sw                                                   &
    , index_absorb_lw                                                   &
    , i_band_k_sw                                                       &
    , i_band_k_lw                                                       &
    , iota_sw                                                           &
    , iota_lw                                                           &
    , subcols_per_k_sw                                                  &
    , subcols_per_k_lw                                                  &
    , order_lw)

    IMPLICIT NONE

!   Input variables
    INTEGER, Intent(IN) :: n_band_sw
!         Number of SW bands

    INTEGER, Intent(IN) :: n_band_lw
!         Number of LW bands

    INTEGER, Intent(IN) :: n_species_sw
!         Number of gases in each SW band

    INTEGER, Intent(IN) :: n_species_lw
!         Number of gases in each LW band

    INTEGER, Intent(IN) :: n_k_term_sw
!         Size allocated for k-terms in each SW band

    INTEGER, Intent(IN) :: n_k_term_lw
!         Size allocated for k-terms in each LW band

    INTEGER, Intent(IN) :: index_absorb_sw(n_species_sw, n_band_sw)
!         indeces of the gases in each SW band

    INTEGER, Intent(IN) :: index_absorb_lw(n_species_lw, n_band_lw)
!         indeces of the gases in each LW band

    INTEGER, Intent(IN) :: i_band_k_sw(n_band_sw, n_species_sw)
!         Number of k-terms in each SW band

    INTEGER, Intent(IN) :: i_band_k_lw(n_band_lw, n_species_lw)
!         Number of k-terms in each LW band

    REAL (RealK), Intent(IN) :: iota_sw(n_band_sw, n_k_term_sw)
!         Importance of each SW k-term (importance/no sub-columns)

    REAL (RealK), Intent(IN) :: iota_lw(n_band_lw, n_k_term_lw)
!         Importance of each LW k-term (importance/no sub-columns)

    INTEGER, Intent(IN) :: subcols_per_k_sw(n_band_sw, n_k_term_sw)
!         Number of sub-columns assigned to each SW k-term

    INTEGER, Intent(IN) :: subcols_per_k_lw(n_band_lw, n_k_term_lw)
!         Number of sub-columns assigned to each LW k-term


!   Output variable
    INTEGER, Intent(OUT) :: order_lw(99999)
!       Re-ordering of sub-columns required to match importance

!   Local variables
    INTEGER :: i, j, k, i_gas
!       Loop variables

    INTEGER :: tot_subcol_sw
!       Total number of SW sub-columns

    INTEGER :: tot_subcol_lw
!       Total number of LW sub-columns

    INTEGER :: cumulative_subcol_sw(n_band_sw, n_k_term_sw)
!       Cumulative number of sub-columns per SW k-term  

    INTEGER :: cumulative_subcol_lw(n_band_lw, n_k_term_lw)
!       Cumulative number of sub-columns per LW k-term  

    REAL (RealK), ALLOCATABLE :: subcol_import_sw(:)
!       Importance of each sub-column in the SW
!       (equals the importance of k-term it's assigned to)

    REAL (RealK), ALLOCATABLE :: subcol_import_lw(:)
!       Importance of each sub-column in the LW
!       (equals the importance of k-term it's assigned to)

    REAL (RealK), ALLOCATABLE :: temp_order(:)
!       Array for temporarily storing the order of importance of the
!       sub-columns

    INTEGER, ALLOCATABLE :: order_sw(:)
!      Order of importance of the SW sub-columns.

!   Calculate total and cumulative number of sub-columns in the SW
    tot_subcol_sw = 0
    tot_subcol_lw = 0
    DO i = 1, n_band_sw
      i_gas = index_absorb_sw(1, i)
      DO j = 1, i_band_k_sw(i, i_gas)
        tot_subcol_sw = tot_subcol_sw+subcols_per_k_sw(i,j)
        cumulative_subcol_sw(i,j) = tot_subcol_sw 
      END DO
    END DO

!   Calculate importance of SW sub-columns
    ALLOCATE(subcol_import_sw(tot_subcol_sw))
    k = 1
    DO i = 1, n_band_sw
      i_gas = index_absorb_sw(1, i)
      DO j = 1, i_band_k_sw(i, i_gas)
        subcol_import_sw(k:cumulative_subcol_sw(i,j)) = iota_sw(i,j)
        k = cumulative_subcol_sw(i,j)+1
      END DO
    ENDDO

!   Calculate total and cumulative number of sub-columns in the LW
    DO i = 1, n_band_lw
      i_gas = index_absorb_lw(1, i)
      DO j = 1, i_band_k_lw(i, i_gas)
        tot_subcol_lw = tot_subcol_lw+subcols_per_k_lw(i,j)
        cumulative_subcol_lw(i,j) = tot_subcol_lw
      END DO
    END DO

!   Calculate importance of SW sub-columns
    ALLOCATE(subcol_import_lw(tot_subcol_lw))
    k = 1
    DO i = 1, n_band_lw
      i_gas = index_absorb_lw(1, i)
      DO j = 1, i_band_k_lw(i, i_gas)      
        subcol_import_lw(k:cumulative_subcol_lw(i,j)) = iota_lw(i,j)
        k = cumulative_subcol_lw(i,j)+1
      END DO
    ENDDO

!   Calculate the order of importance of the sub-columns in the SW
    ALLOCATE(order_sw(tot_subcol_sw))
    ALLOCATE(temp_order(tot_subcol_sw))
    temp_order = subcol_import_sw
    DO i = 1, tot_subcol_sw
      order_sw(i) = MAXLOC(temp_order,1)
      temp_order(order_sw(i)) = -999
    ENDDO
    IF (ALLOCATED(temp_order)) DEALLOCATE(temp_order)

!   Calculate the order of importance of the sub-columns in the LW
    ALLOCATE(temp_order(tot_subcol_lw))
    temp_order = subcol_import_lw
    DO i = 1, tot_subcol_lw
      order_lw(i) = MAXLOC(temp_order,1)
      temp_order(order_lw(i)) = -999
    ENDDO

!   Calculate the reordering of sub-columns necessary for matching rank
    temp_order = order_lw(1:tot_subcol_lw)
    k = tot_subcol_sw
    DO i = 1, tot_subcol_lw
      j = 1
      DO WHILE (temp_order(j) /= i)
        j = j+1
      ENDDO
      IF (j > tot_subcol_sw) THEN
        k = k+1
        order_lw(i) = k
      ELSE
        order_lw(i) = order_sw(j)
      ENDIF
    ENDDO

    IF (ALLOCATED(temp_order)) DEALLOCATE(temp_order)
    IF (ALLOCATED(order_sw)) DEALLOCATE(order_sw)
    
  END SUBROUTINE

END PROGRAM
