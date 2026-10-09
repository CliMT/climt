! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the pressure and temperature for c-k.
!
SUBROUTINE set_condition_ck_90 &
!
(l_interactive, &
 nd_pt, n_pt_pair, n_p, p, t, &
 nd_band, nd_species, &
 n_absorb, type_absorb, i_gas, i_index, &
 i_gas_1, i_index_1, i_gas_2, i_index_2, &
 n_band, n_band_absorb, index_absorb, &
 nd_continuum, n_band_continuum, index_continuum, &
 l_fit_line_data, l_fit_frn_continuum, l_fit_self_continuum, l_fit_cont_data, &
 umin_c, umax_c, n_path_c, n_pp, l_access_HITRAN, l_access_xsc, &
 l_access_cia, include_h2o_foreign_continuum, &
 l_use_h2o_frn_param, l_use_h2o_self_param, l_cont_line_abs_weight, &
 n_selected_band, list_band, &
 i_ck_fit, tol, max_path, max_path_wgt, n_k, nu_inc_0, line_cutoff, &
 l_ckd_cutoff, l_scale_pT, i_type_residual, i_scale_fnc, p_ref, t_ref, &
 l_load_map, l_load_wgt, l_save_map, file_map, &
 i_line_prof_corr, l_self_broadening, n_gas_frac, gas_frac, nd_gas_frac, &
 ierr &
)
!
! Description:
!   The pressures and temperatures and the amounts of gases
!   are set, either interactively or by reading from a file.
!
!
!
! Modules used
  USE realtype_rd
  USE def_std_io_icf
  USE ck_fit_pcf
  USE type_residual_pcf
  USE rad_pcf
  USE gas_list_pcf
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments.
!
! Sizes of arrays:
  INTEGER, Intent(IN) :: nd_pt
!   Size allocated for arrays of pressure and temperature at which
!   k-distributions are to be calculated
  INTEGER, Intent(IN) :: nd_band
!   Size allocated for spectral bands
  INTEGER, Intent(IN) :: nd_species
!   Size allocated for absorbing species
  INTEGER, Intent(IN) :: nd_continuum
!   Size allocated for continua
  INTEGER, Intent(IN) :: nd_gas_frac
!   Size allocated for gas fractions
!
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive operation
  INTEGER, Intent(IN) :: n_band
!   Number of bands
  INTEGER, Intent(IN) :: n_absorb
!   Number of gaseous absorbers
  INTEGER, Intent(IN) :: type_absorb(nd_species)
!   Identifiers of absorbers
  INTEGER, Intent(IN) :: n_band_absorb(nd_band)
!   Number of active absorbers in each band
  INTEGER, Intent(IN) :: index_absorb(nd_species, nd_band)
!   Lists of absorbers which are active in particular bands
!
  INTEGER, Intent(IN) :: n_band_continuum(nd_band)
!   Number of active continua in each band
  INTEGER, Intent(IN) :: index_continuum(nd_band, nd_continuum)
!   Lists of continua which are active in particular bands
  LOGICAL, Intent(IN) :: l_access_HITRAN
!   Flag for access to HITRAN database
  LOGICAL, Intent(IN) :: l_access_xsc
!   Flag for HITRAN cross-section file
  LOGICAL, Intent(IN) :: l_access_cia
!   Flag for HITRAN CIA file
!
  LOGICAL, Intent(OUT) :: l_fit_line_data
!   Flag requiring the fitting of line data
  LOGICAL, Intent(OUT) :: l_fit_frn_continuum
!   Flag requiring the fitting of foreign-broadened continuum data
  LOGICAL, Intent(OUT) :: l_fit_self_continuum
!   Flag requiring the fitting of self-broadened continuum data
  LOGICAL, Intent(OUT) :: l_fit_cont_data
!   Flag requiring generation of generalised continuum data
!
  INTEGER, Intent(OUT) :: i_gas
!   Identifiers of gas for which data are to be generated
  INTEGER, Intent(OUT) :: i_index
!   Identifiers of gas for which data are to be generated
  INTEGER, Intent(OUT) :: i_gas_1
!   Identifier of first continuum gas for which data are to be generated
  INTEGER, Intent(OUT) :: i_index_1
!   Identifier of first continuum gas for which data are to be generated
  INTEGER, Intent(OUT) :: i_gas_2
!   Identifier of second continuum gas for which data are to be generated
  INTEGER, Intent(OUT) :: i_index_2
!   Identifier of second continuum gas for which data are to be generated
  INTEGER, Intent(OUT) :: n_selected_band
!   Number of spectral bands in the range where the gas is active
  INTEGER, Intent(OUT) :: list_band(nd_band)
!   Number of spectral bands in the range where the gas is active
  LOGICAL, Intent(OUT) :: include_h2o_foreign_continuum
!   The foreign-broadened continuum is included with the H2O line
!   absorption
  LOGICAL, Intent(OUT) :: l_use_h2o_frn_param
!   Use foreign-broadened water vapor continuum parametrisation provided
!   in the code
  LOGICAL, Intent(OUT) :: l_use_h2o_self_param
!   Use self-broadened water vapor continuum parametrisation provided
!   in the code
  LOGICAL, Intent(OUT) :: l_cont_line_abs_weight
!   Use line absorption as weighting in continuum absorption transmissions
!
  REAL  (RealK), Intent(OUT) :: umin_c
!   Minimum pathlength for continuum absorption
  REAL  (RealK), Intent(OUT) :: umax_c
!   Maximum pathlength for continuum absorption
  INTEGER, Intent(OUT) :: n_path_c
!   Number of pathlengths for continuum absorption
  INTEGER, Intent(OUT) :: n_pp
!   Number of values of the partial pressure to be used
!
  INTEGER, Intent(OUT) :: i_ck_fit
!   Type of fitting to be carried out
  REAL  (RealK), Intent(OUT) :: tol
!   Tolerance for the fit
  REAL  (RealK), Intent(OUT) :: max_path
!   Maximum pathlength to be considered for the absorber (this
!   will be used to determine how weak absorption must be to be
!   considered as grey)
  REAL  (RealK), Intent(OUT) :: max_path_wgt
!   Maximum pathlength to be considered for the absorber used for weighting
!   in continuum transmissions
  INTEGER, Intent(OUT) :: n_k(nd_band)
!   Number of terms in the fit
!
  REAL  (RealK), Intent(OUT) :: nu_inc_0
!   Default spacing for integration in frequency
  REAL  (RealK), Intent(OUT) :: line_cutoff
!   Distance from line centre beyond which absorption is neglected
  LOGICAL, Intent(OUT) ::  l_ckd_cutoff
!   Line cutoff to be consistent with the CKD continuum
!
  INTEGER, Intent(OUT) :: n_pt_pair
!   Number of p,t pairs
  INTEGER, Intent(OUT) :: n_p
!   Number of unique pressures
  REAL  (RealK), Intent(OUT) :: p(nd_pt)
!   Pressures for generating data
  REAL  (RealK), Intent(OUT) :: t(nd_pt)
!   Temperatures for generating data
!
  LOGICAL, Intent(OUT) :: l_scale_pT
!   Flag indicating that a scaling with the pressure and temperature
!   is required
  INTEGER, Intent(OUT) :: i_type_residual
!   Type of residual used
  INTEGER, Intent(OUT) :: i_scale_fnc
!   Identifier of the scaling function used
  REAL  (RealK), Intent(OUT) :: p_ref(nd_band)
!   Reference pressures for scaling
  REAL  (RealK), Intent(OUT) :: t_ref(nd_band)
!   Reference temperatures for scaling
!
  INTEGER, Intent(OUT) :: i_line_prof_corr
!   Line profile correction type
!
  LOGICAL, Intent(OUT) :: l_self_broadening
!   Flags to include effects of self-broadening
  INTEGER, Intent(OUT) :: n_gas_frac
!   Number of gas fractions at which to tabulate ESFT terms
  REAL  (RealK), Intent(OUT) :: gas_frac(nd_gas_frac)
!   List of gas fractions at which to tabulate ESFT terms
!
  LOGICAL, Intent(OUT) :: l_load_map
!   Use pre-defined mapping of wavenumbers to g-space
  LOGICAL, Intent(OUT) :: l_load_wgt
!   Use pre-defined k-term weights
  LOGICAL, Intent(OUT) :: l_save_map
!   Save mapping of wavenumbers to g-space
  CHARACTER(LEN=*), Intent(OUT) :: file_map
!   Name of file with mapping
!
!
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
!
!
!
! Local variables.
!
  CHARACTER (LEN=1) :: char_yn
!   Interactive response
  CHARACTER (LEN=1) :: char_if
!   Interactive response
  INTEGER :: iu_file_in
!   Unit number for input from a file
  CHARACTER (LEN=80) :: line
!   Line of data
  INTEGER :: first_band
!   First band considered
  INTEGER :: last_band
!   Last band considered
  INTEGER :: ios
!   I/O error flag
  INTEGER :: i_swap
!   Swapping variable
  INTEGER :: ib
!   Loop variable
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  LOGICAL :: l_reverse
!   Flag for reversing order of bands
  LOGICAL :: l_finish
!   Flag for ending input of p and T
  LOGICAL :: l_file
!   Flag for data from a file
!
! Subroutines called:
  EXTERNAL &
    get_free_unit, open_file_in, set_extern_ckd_frn_data
!
  INTERFACE

    SUBROUTINE map_heap_func(a, map)
      USE realtype_rd
      REAL  (RealK), Intent(IN), Dimension(:) :: a
      INTEGER, Intent(OUT), Dimension(:) :: map
    END SUBROUTINE map_heap_func

    SUBROUTINE read_ref_pt_90(i_gas, &
      n_selected_band, list_band, &
      p_ref, t_ref, &
      ierr)
!
      USE realtype_rd
!
      INTEGER, Intent(IN) :: i_gas
      INTEGER, Intent(IN) :: n_selected_band
      INTEGER, Intent(IN) :: list_band(:)
      REAL  (RealK), Intent(OUT) :: p_ref(:)
      REAL  (RealK), Intent(OUT) :: t_ref(:)
      INTEGER, Intent(INOUT) :: ierr
!
    END SUBROUTINE read_ref_pt_90
!
  END INTERFACE
!
!
!
  CALL select_data_type_int
!
  IF (l_fit_line_data .OR. l_fit_cont_data) CALL select_gas_int
!
  l_use_h2o_frn_param = .FALSE.
  l_use_h2o_self_param = .FALSE.
  IF (l_fit_cont_data .AND. (.NOT. l_access_cia)) &
    CALL select_h2o_cont
!
  IF (l_fit_cont_data) CALL select_cont_weight
!
  CALL select_bands_int
!
  CALL select_generation_pT_int
!
  CALL select_scaling_int
!
  IF (l_access_HITRAN .OR. l_access_xsc .OR. l_fit_cont_data) &
    CALL select_line_details_int
!
  IF (l_fit_line_data) THEN
    CALL select_self_broadening
  ELSE
    l_self_broadening = .FALSE.
    n_gas_frac=1
    gas_frac=0.0_RealK
  END IF
!
  IF (l_fit_self_continuum .OR. l_fit_frn_continuum) &
    CALL select_cont_details_int
!
  IF (l_fit_line_data .OR. l_fit_cont_data) THEN
    CALL select_fit_type_int
  ELSE
    i_ck_fit = IP_ck_none
  END IF

  l_load_map = .FALSE.
  l_load_wgt = .FALSE.
  l_save_map = .FALSE.
  IF (i_ck_fit /= ip_ck_none) CALL select_mapping_wgt
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
  SUBROUTINE select_data_type_int
!
!
!
!   Select the mode of operation.
    WRITE(*, "(a)") "Are line absorption data to be generated? (Y/N)"
    DO
      READ(*, "(a)") char_yn
      IF ( (char_yn == 'Y') .OR. (char_yn == 'y') ) THEN
        l_fit_line_data=.TRUE.
        EXIT
      ELSE IF ( (char_yn == 'N') .OR. (char_yn == 'n') ) THEN
        l_fit_line_data=.FALSE.
        EXIT
      ELSE
        WRITE(iu_err, "(/a)") "Erroneous response"
        IF (l_interactive) THEN
          WRITE(*, "(a)") "Please re-enter."
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
    ENDDO
!
    IF (.NOT. l_fit_line_data) THEN
      WRITE(*, "(a)") &
        "Are generalised continuum data to be generated? (Y/N)"
      DO
        READ(*, "(a)") char_yn
        IF ( (char_yn == 'Y') .OR. (char_yn == 'y') ) THEN
          l_fit_cont_data=.TRUE.
          EXIT
        ELSE IF ( (char_yn == 'N') .OR. (char_yn == 'n') ) THEN
          l_fit_cont_data=.FALSE.
          EXIT
        ELSE
          WRITE(iu_err, "(/a)") "Erroneous response"
          IF (l_interactive) THEN
            WRITE(*, "(a)") "Please re-enter."
          ELSE
            ierr=i_err_fatal
            RETURN
          ENDIF
        ENDIF
      ENDDO
    ELSE
      l_fit_cont_data=.FALSE.
    END IF
!
    WRITE(*, "(a)") "Are foreign continuum data to be generated? (Y/N)"
    DO
      READ(*, "(a)") char_yn
      IF ( (char_yn == 'Y') .OR. (char_yn == 'y') ) THEN
        l_fit_frn_continuum=.TRUE.
        CALL set_extern_ckd_frn_data(ierr)
        EXIT
      ELSE IF ( (char_yn == 'N') .OR. (char_yn == 'n') ) THEN
        l_fit_frn_continuum=.FALSE.
        EXIT
      ELSE
        WRITE(iu_err, "(/a)") "Erroneous response"
        IF (l_interactive) THEN
          WRITE(*, "(a)") "Please re-enter."
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
    ENDDO
!
    WRITE(*, "(a)") &
      "Are self-broadened continuum data to be generated? (Y/N)"
    DO
      READ(*, "(a)") char_yn
      IF ( (char_yn == 'Y') .OR. (char_yn == 'y') ) THEN
        l_fit_self_continuum=.TRUE.
        CALL set_extern_ckd_self_data(ierr)
        EXIT
      ELSE IF ( (char_yn == 'N') .OR. (char_yn == 'n') ) THEN
        l_fit_self_continuum=.FALSE.
        EXIT
      ELSE
        WRITE(iu_err, "(/a)") "Erroneous response"
        IF (l_interactive) THEN
          WRITE(*, "(a)") "Please re-enter."
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
    ENDDO
!
    IF ( l_fit_frn_continuum .OR. l_fit_self_continuum ) THEN
      i_gas=IP_H2O
    ENDIF
!
!
!
  END SUBROUTINE select_data_type_int
!
!
!
  SUBROUTINE select_gas_int
!
!
!
!   Select the gas to be considered.
    WRITE(*, "(a)") "Enter the identifier(s) for the gas(es) to be considered."
    DO
      IF (l_fit_cont_data) THEN
        READ(*, *, IOSTAT=ios) i_gas_1, i_gas_2
        i_gas=0
      ELSE
        READ(*, *, IOSTAT=ios) i_gas
        i_gas_1=0
        i_gas_2=0
      END IF
      IF (ios /= 0) THEN
        WRITE(iu_err, "(/a)") "Erroneous response"
        IF (l_interactive) THEN
          WRITE(*, "(a)") "Please re-enter."
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ELSE
!       Check that the gas is in the spectral file.
        IF (l_fit_cont_data) THEN
          i_index_1=0
          i_index_2=0
          DO i=1,n_absorb
            IF (type_absorb(i) == i_gas_1) THEN
              i_index_1=i
              EXIT
            ENDIF
          ENDDO
          DO i=1,n_absorb
            IF (type_absorb(i) == i_gas_2) THEN
              i_index_2=i
              EXIT
            ENDIF
          ENDDO
        ELSE
          i_index=0
          DO i=1,n_absorb
            IF (type_absorb(i) == i_gas) THEN
              i_index=i
              EXIT
            ENDIF
          ENDDO
        END IF
        IF (l_fit_cont_data) THEN
          IF (i_index_1 == 0 .OR. i_index_2 == 0) THEN
            WRITE(iu_err, "(a)") &
              "One or both continuum gases is not in the spectral file."
            IF (l_interactive) THEN
              WRITE(*, "(a)") "Please re-enter."
            ELSE
              ierr=i_err_fatal
              RETURN
            ENDIF
          ELSE
            EXIT
          END IF
        ELSE IF (i_index == 0) THEN
          WRITE(iu_err, "(a)") "This gas is not in the spectral file."
          IF (l_interactive) THEN
            WRITE(*, "(a)") "Please re-enter."
          ELSE
            ierr=i_err_fatal
            RETURN
          ENDIF
        ELSE
          EXIT
        ENDIF
      ENDIF
    ENDDO
!
!
!
  END SUBROUTINE select_gas_int
!
!
!
  SUBROUTINE select_h2o_cont
!
!
!
!   Check if the selected continuum is one  of the two water vapour continua.
    IF ((i_gas_1 == ip_h2o .AND. i_gas_2 == ip_air) .OR. &
        (i_gas_1 == ip_air .AND. i_gas_2 == ip_h2o)) THEN
!     This is the foreign broadened water vapour continuum
!
      WRITE(*, '(a)') 'Do you wish to use the foreign-broadened water ' // &
        'vapour continuum included in the code? (Y/N)'
      DO
        READ(*, '(a)') char_yn
        IF ( (char_yn == 'Y') .OR. (char_yn == 'y') ) THEN
          l_use_h2o_frn_param=.TRUE.
          CALL set_extern_ckd_frn_data(ierr)
          EXIT
        ELSE IF ( (char_yn == 'N') .OR. (char_yn == 'n') ) THEN
          EXIT
        ELSE
          WRITE(iu_err, '(/a)') 'Erroneous response'
          IF (l_interactive) THEN
            WRITE(*, '(a)') 'Please re-enter.'
          ELSE
            ierr=i_err_fatal
            RETURN
          ENDIF
        ENDIF
      ENDDO
!
    ELSE IF (i_gas_1 == ip_h2o .AND. i_gas_2 == ip_h2o) THEN
!     This is the self-broadened water vapour continuum
!
      WRITE(*, '(a)') 'Do you wish to use the self-broadened water ' // &
        'vapour continuum included in the code? (Y/N)'
      DO
        READ(*, '(a)') char_yn
        IF ( (char_yn == 'Y') .OR. (char_yn == 'y') ) THEN
          l_use_h2o_self_param=.TRUE.
          CALL set_extern_ckd_self_data(ierr)
          EXIT
        ELSE IF ( (char_yn == 'N') .OR. (char_yn == 'n') ) THEN
          EXIT
        ELSE
          WRITE(iu_err, '(/a)') 'Erroneous response'
          IF (l_interactive) THEN
            WRITE(*, '(a)') 'Please re-enter.'
          ELSE
            ierr=i_err_fatal
            RETURN
          ENDIF
        ENDIF
      ENDDO
!
    END IF
!
!
!
  END SUBROUTINE select_h2o_cont
!
!
!
  SUBROUTINE select_cont_weight
!
!
!
!   Select weighting of continuum transmissions
    WRITE(*, '(a)') 'Do you wish to use line absorption transmissions as ' // &
        'weighting in continuum transmissions? (Y/N)'
    DO
      READ(*, '(a)') char_yn
      IF ( (char_yn == 'Y') .OR. (char_yn == 'y') ) THEN
        l_cont_line_abs_weight=.TRUE.
        
!       Select gas to use in weighting.
        WRITE(*, "(a)") "Enter the identifier for the weighting gas."
        DO
          READ(*, *, IOSTAT=ios) i_gas
          IF (ios /= 0) THEN
            WRITE(iu_err, "(/a)") "Erroneous response"
            IF (l_interactive) THEN
              WRITE(*, "(a)") "Please re-enter."
            ELSE
              ierr=i_err_fatal
              RETURN
            ENDIF
          ELSE
            WRITE(*, "(a)") "Enter the maximum pathlength for this gas."
            DO
              READ(*, *, IOSTAT=ios) max_path_wgt
              IF (ios /= 0) THEN
                WRITE(iu_err, "(/a)") "Erroneous response"
                IF (l_interactive) THEN
                  WRITE(*, "(a)") "Please re-enter."
                ELSE
                  ierr=i_err_fatal
                  RETURN
                ENDIF
              ELSE
                EXIT
              ENDIF
            ENDDO
            EXIT
          ENDIF
        ENDDO
        EXIT
      ELSE IF ( (char_yn == 'N') .OR. (char_yn == 'n') ) THEN
        l_cont_line_abs_weight=.FALSE.
        i_gas=0
        EXIT
      ELSE
        WRITE(iu_err, '(/a)') 'Erroneous response'
        IF (l_interactive) THEN
          WRITE(*, '(a)') 'Please re-enter.'
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
    ENDDO
!
!
!
  END SUBROUTINE select_cont_weight
!
!
!
  SUBROUTINE select_bands_int
!
!
!
!   Select the range of bands to be covered.
    WRITE(iu_stdout, '(/a)') 'Enter first and last bands to be considered.'
    DO 
      READ(iu_stdin, *, IOSTAT=ios) first_band, last_band
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ Erroneous response.'
        IF(l_interactive) THEN
          WRITE(iu_stdout, '(a)') 'please re-enter.'
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ELSE
        l_reverse=.false.
        IF (first_band > last_band) THEN
          i_swap=first_band
          first_band=last_band
          last_band=i_swap
          l_reverse=.true.
        ENDIF
!       Check limits.
        IF ( (first_band < 1) .OR. &
             (last_band > n_band) ) THEN
          WRITE(iu_err, '(a)') '+++ Response out of range.'
          IF (l_interactive) THEN
            WRITE(iu_stdout, '(a)') 'Please re-enter.'
          ELSE
            ierr=i_err_fatal
            RETURN
          ENDIF
        ELSE
          IF (l_fit_line_data .OR. l_fit_cont_data) THEN
!           Make the list of valid bands for this absorber.
            n_selected_band=0
            IF (l_reverse) THEN
              DO i=last_band, first_band, -1
                n_selected_band = n_selected_band+1
                list_band(n_selected_band) = i
              END DO
            ELSE
              DO i=first_band, last_band
                n_selected_band = n_selected_band+1
                list_band(n_selected_band) = i
              END DO
            END IF
          ELSE IF (l_fit_frn_continuum) THEN
!           Make the list of valid bands for this continuum.
            n_selected_band=0
            DO i=first_band, last_band
              n_selected_band = n_selected_band+1
              list_band(n_selected_band) = i
            END DO
          ELSE IF (l_fit_self_continuum) THEN
!           Make the list of valid bands for this continuum.
            n_selected_band=0
            DO i=first_band, last_band
              n_selected_band = n_selected_band+1
              list_band(n_selected_band) = i
            END DO
          END IF
          EXIT
        ENDIF
      ENDIF
    ENDDO
!
!
!
  END SUBROUTINE select_bands_int



  SUBROUTINE sorting_PT_points(P,T,map)

    IMPLICIT NONE

    ! Input P,T points
    REAL(RealK), DIMENSION(:), Intent(IN)   :: P,T
    ! Map indexing properly sorted P,T 
    INTEGER,DIMENSION(:), Intent(OUT):: map

    INTEGER,DIMENSION(SIZE(P)) :: Pmap,Pmap2
    INTEGER :: PuniqueNum ! Counts the number of unique pressure points
    INTEGER :: PTLen,i,j, PCount,idx,mapIdx
    REAL(RealK),DIMENSION(SIZE(P)) :: Pcopy, Tcopy,Punique
    INTEGER,DIMENSION(:), ALLOCATABLE :: PIdxSelection,Tmap
    REAL(RealK),DIMENSION(:), ALLOCATABLE :: TSelection
    LOGICAL :: PisPresent

    PTLen = SIZE(P)
    ! Sort the P points
    Pcopy=P
    Tcopy=T
    CALL map_heap_func(Pcopy,Pmap)
    Pcopy=Pcopy(Pmap)
    Tcopy=Tcopy(Pmap)

    ! Adds first element to Punique
    PuniqueNum = 1
    Punique(PuniqueNum)=P(1)
    ! Find unique P points
    DO idx =2,PTLen
      PisPresent=.FALSE.
      DO i = 1, PuniqueNum
        IF (P(idx) == Punique(i)) THEN
          PisPresent=.TRUE.
        END IF
      END DO
      IF (.NOT. PisPresent) THEN
        PuniqueNum =PuniqueNum+1
        Punique(PuniqueNum)=P(idx)
      END IF
    END DO
    !write(*,*) PuniqueNum, " unique Pressure Points"

    ! For each P, find every T, sort them, and output to map
    mapIdx=0
    DO idx=1,PuniqueNum
      PCount=0
      ! Find the Pressure indices corresponding to each P
      DO i = 1, PTLen
        IF (Punique(idx) == P(i)) THEN
          PCount=PCount+1
          Pmap2(Pcount)=i
        END IF
      END DO
      
      ! Find the Ts for each index
      ALLOCATE(TSelection(PCount))
      ALLOCATE(PIdxSelection(PCount))
      ALLOCATE(Tmap(PCount))
      ! Matching indices in the master lists and Ts
      PIdxSelection=Pmap2(1:PCount)
      TSelection=T(Pmap2(1:PCount))
      ! Sort T to find the proper order of indices
      CALL map_heap_func(TSelection,Tmap)
      PIdxSelection=PIdxSelection(Tmap)
      ! Write them to map for output
      DO j = 1, PCount
        mapIdx=mapIdx+1
        map(mapIdx) = PIdxSelection(j)
      END DO
      DEALLOCATE(TSelection)
      DEALLOCATE(PIdxSelection)
      DEALLOCATE(Tmap)
    END DO
  END SUBROUTINE sorting_PT_points
    
    
    
  SUBROUTINE select_generation_pT_int
!  
    INTEGER :: map_pt(nd_pt)
    
!  
!   Determine the pressures and temperatures at which data are to be
!   generated.
    WRITE(iu_stdout, '(/a)') &
      'Setting of pressures and temperatures:'
    WRITE(iu_stdout, '(/a)') &
      'Enter "f" to read from a file or "i" to set values interactively.'
    DO
      READ(iu_stdin, '(A)', iostat=ios) char_if
      IF  (ios /= 0) THEN
        WRITE(iu_err, '(/a)') '+++ Unrecognised reponse.'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(/a)') 'Please re-enter.'
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ELSE IF ( (char_if == "F").OR.(char_if == "f") ) THEN
        l_file=.true.
        CALL get_free_unit(ierr, iu_file_in)
        CALL open_file_in(ierr, iu_file_in, &
          'Enter name of file containing p/T-values.')
!       Read the lines of the file until the starting directive 
!       "*PTVAL" is encountered.
        DO
          READ(iu_file_in, '(A)', IOSTAT=ios) line
          IF  (ios /= 0) THEN
            WRITE(iu_err, '(/a)') '+++ Incorrect file: missing directive ?'
            ierr=i_err_fatal
            RETURN
          ELSE IF (line(1:6) == '*PTVAL') THEN
            EXIT
          ENDIF
        ENDDO
        EXIT
      ELSE IF ( (char_if.eq."I").or.(char_if.eq."i") ) THEN
        l_file=.false.
        EXIT
      ELSE
        WRITE(iu_err, '(/a)') '+++ Illegal input.'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(/a)') 'Please re-enter.'
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
    ENDDO
!  
!  
    n_pt_pair=0
    n_p=0
    l_finish=.false.
    DO
      CALL read_pt_line_90(ierr, l_interactive, l_file, &
        iu_file_in, &
        nd_pt, n_pt_pair, p, t, l_finish)
      IF (ierr /= i_normal) RETURN
      IF (l_finish) EXIT
      n_p = n_p + 1
    ENDDO
!  
    IF (n_pt_pair == 0) THEN
      WRITE(iu_err, '(a/)') &
        'At least one pair of values must be specified.'
      ierr=i_err_fatal
      RETURN
    ENDIF
!  
!   Close the file to allow the unit to be re-used.
    IF (l_file) CLOSE(iu_file_in)

!   Sort the P/T pairs into increasing P and T:
    Call sorting_PT_points( &
         p(1:n_pt_pair),t(1:n_pt_pair),map_pt(1:n_pt_pair))
    
    p(1:n_pt_pair)=p(map_pt(1:n_pt_pair))
    t(1:n_pt_pair)=t(map_pt(1:n_pt_pair))
!
!
!
  END SUBROUTINE select_generation_pT_int
!
!
!
  SUBROUTINE select_scaling_int
!  
!  
!  
    WRITE(*, "(a)") "Are scaling functions/lookup tables required? (Y/N)"
    DO
      READ(*, "(A)") char_yn
      IF ( (char_yn == 'Y') .OR. (char_yn == 'y') ) THEN
        l_scale_pT=.TRUE.
!  
!       Determine the scaling function.
        WRITE(*, "(/a)") &
          "Enter the type of scaling function."
        DO
          READ(*, *, IOSTAT=ios) i_scale_fnc
          IF ( (ios /= 0) .OR. &
               ( (i_scale_fnc /= IP_scale_power_law) .AND. &
                 (i_scale_fnc /= IP_scale_power_quad) .AND. &
                 (i_scale_fnc /= IP_scale_doppler_quad) .AND. &
                 (i_scale_fnc /= IP_scale_dbl_pow_law) .AND. &
                 (i_scale_fnc /= IP_scale_dbl_pow_quad) .AND. &
                 (i_scale_fnc /= IP_scale_dbl_dop_quad) .AND. &
                 (i_scale_fnc /= IP_scale_lookup) .AND. &
                 (i_scale_fnc /= IP_scale_t_lookup) ) ) THEN
            WRITE(iu_err, "(a)") "Erroneous response."
            IF (l_interactive) THEN
              WRITE(*, "(a)") "Please re-enter"
            ELSE
              ierr=i_err_fatal
              RETURN
            ENDIF
          ELSE
            EXIT
          ENDIF
        ENDDO

        IF (i_scale_fnc == IP_scale_lookup .OR. &
            i_scale_fnc == IP_scale_t_lookup) THEN
!         Set to arbitrary values (not used).
          p_ref=1.0_RealK
          t_ref=200.0_RealK
        ELSE
!         Determine reference conditions of scaling if required.
!         These may be set interactively or from a file.
          WRITE(*, "(/a)") &
            "Are the reference conditions to be set interactively "// &
            "or from a file? (I/F)"
          DO
            READ(*, "(a)") char_if
            IF ( (char_if == 'I') .OR. (char_if == 'i') ) THEN
              DO i=1, n_selected_band
                ib=list_band(i)
                WRITE(*, "(a, i5)") &
                  "Enter reference pressure and temperature in band ", ib
                DO
                  READ(*, *, IOSTAT=ios) p_ref(ib), t_ref(ib)
                  IF (ios /= 0) THEN
                    WRITE(iu_err, "(/a, i5)") &
                      "Erroneous specification of reference " // &
                      "conditions in band ", ib
                    IF (l_interactive) THEN
                      WRITE(*, "(a)") "Please re-enter."
                    ELSE
                      ierr=i_err_fatal
                      STOP
                    ENDIF
                  ELSE
                    EXIT
                  ENDIF
                ENDDO
              ENDDO
!         
            ELSE IF ( (char_if == 'F') .OR. (char_if == 'f') ) THEN
!         
              CALL read_ref_pt_90(i_gas, &
                n_selected_band, list_band, &
                p_ref, t_ref, ierr)
              IF (ierr /= i_normal) STOP
!         
            ELSE
!         
              WRITE(iu_err, "(a)") "Erroneous response"
              IF (l_interactive) THEN
                WRITE(*, "(a)") "Please re-enter."
              ELSE
                ierr=i_err_fatal
                STOP
              ENDIF
!         
            ENDIF
!         
            EXIT
!         
          ENDDO
        END IF
!  
!       Only the option of minimizing errors in the transmission is
!       currently permitted.
        i_type_residual=ip_scale_trans_residual
!  
        EXIT
!  
      ELSE IF ( (char_yn == 'N') .OR. (char_yn == 'n') ) THEN
!  
        l_scale_pT=.FALSE.
        i_scale_fnc=IP_scale_fnc_null
        EXIT
!  
      ELSE
!  
        WRITE(iu_err, "(a)") 'Invalid response.'
        IF (l_interactive) THEN
          WRITE(*, "(a)") "Please re-enter."
        ELSE
          ierr=i_err_fatal
          STOP
        ENDIF
      ENDIF
!  
    ENDDO
!
!
!
  END SUBROUTINE select_scaling_int
!
!
!
  SUBROUTINE select_line_details_int
!
    USE line_prof_corr_mod, ONLY: ip_lpc_unity
!
    l_ckd_cutoff = .FALSE.
    include_h2o_foreign_continuum=.FALSE.

    IF (l_access_HITRAN) THEN
!     Parameters for the integration over lines.
      IF (i_gas == ip_h2o) THEN
!
!       The CKD continuum model should be used only with a line cut-off 
!       of 25 cm-1 (2500 m-1): the line absorption within the cut-off 
!       is reduced by its value at the cut-off, and this absorption is
!       assigned to the continuum. This prescription keeps the continuum
!       contributions free of dicsontinuities.
        WRITE(iu_stdout, '(/a)') &
          'Are data to be adjusted for use with the CKD continuum?'
        READ(iu_stdin, "(A)") char_yn
        IF ( (char_yn == 'Y') .OR. (char_yn == 'y') ) THEN
          l_ckd_cutoff = .TRUE.
          WRITE(iu_stdout, '(/a)') &
            'Ensure that the cut-off is consistent with the CKD model.'
        ELSE
          l_ckd_cutoff = .FALSE.
        ENDIF
!
!       Although the foreign broadened continuum depends on (p-e), e/p
!       is generally small in regions of interest, so for fast application
!       in a GCM we can effectively treat the foreign continuum as depending
!       on the amount of water vapour with a pressure scaling, which means
!       that it may be counted along with line data.
        WRITE(iu_stdout, '(/a)') &
          'Do you wish to include the foreign continuum with the lines?'
        DO
          READ(iu_stdin, *) char_yn
          IF ( (char_yn == 'y') .OR. (char_yn == 'Y') ) THEN
            include_h2o_foreign_continuum=.TRUE.
            CALL set_extern_ckd_frn_data(ierr)
            IF (ierr /= i_normal) RETURN
            EXIT
          ELSE IF ( (char_yn == 'n') .OR. (char_yn == 'N') ) THEN
            include_h2o_foreign_continuum=.FALSE.
            EXIT
          ELSE
            WRITE(iu_err, '(a)') '+++ Erroneous response.'
            IF (l_interactive) THEN
              WRITE(iu_stdout, '(a)') 'Please re-enter.'
            ELSE
              ierr=i_err_fatal
              RETURN
            ENDIF
          ENDIF
        ENDDO
      ENDIF
!
      WRITE(iu_stdout, '(/a)') &
        'Enter the line-cutoff in m-1'
      DO
        READ(iu_stdin, *, IOSTAT=ios) line_cutoff
        IF (ios == 0) THEN
          EXIT
        ELSE
          CALL check_ios_int
        ENDIF
      ENDDO
!
      WRITE(iu_stdout, '(/a, i1, a)') &
        'Enter the type of line profile correction (enter ', &
        ip_lpc_unity, ' for no correction to the Voigt profile).'
      DO
        READ(iu_stdin, *, IOSTAT=ios) i_line_prof_corr
        IF (ios == 0) THEN
          EXIT
        ELSE
          CALL check_ios_int
        ENDIF
      ENDDO
    END IF
!
    WRITE(iu_stdout, '(/a)') &
      'Enter the frequency increment for integration in m-1'
    DO
      READ(iu_stdin, *, IOSTAT=ios) nu_inc_0
      IF (ios == 0) THEN
        EXIT
      ELSE
        CALL check_ios_int
      ENDIF
    ENDDO
!
!
!
  END SUBROUTINE select_line_details_int
!
!
!
  SUBROUTINE select_self_broadening
!
    LOGICAL :: l_data_region
!     Flag for input of data
!
    WRITE(iu_stdout, '(/a)') 'Is self-broadening required?'
    DO
      READ(iu_stdin, '(a)') char_yn
      IF ( (char_yn == 'Y') .OR. (char_yn.eq.'y') ) THEN
        l_self_broadening =  .TRUE. 
        WRITE(iu_stdout, '(a)') &
          'A list of gas fractions must also be supplied.'
        EXIT
      ELSE IF ( (char_yn == 'N') .OR. (char_yn.eq.'n') ) THEN
        l_self_broadening =  .FALSE. 
        n_gas_frac = 1
        EXIT
      ELSE
        WRITE(iu_err, '(a)') '+++ Erroneous response.'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(a)') 'Please re-enter.'
        ELSE
          STOP
        ENDIF
      ENDIF
    ENDDO
!
    WRITE(iu_stdout, '(/)')
!
    IF (l_self_broadening) THEN
!     Set a unit number for generic input
      CALL get_free_unit(ierr, iu_file_in)
      IF (ierr /= i_normal) RETURN

      CALL open_file_in(ierr, iu_file_in, &
        'Specify the file with gas fractions.')
      IF (ierr /= i_normal) RETURN
      l_data_region = .FALSE. 
      n_gas_frac = 0
      DO
        READ(iu_file_in, '(a)', IOSTAT=ios) line
        IF (ios /= 0) EXIT
        IF (l_data_region) THEN
          IF (line(1:4) /= '*END') THEN
            BACKSPACE(iu_file_in)
            n_gas_frac = n_gas_frac + 1
            IF (n_gas_frac > nd_gas_frac) THEN
              WRITE(iu_err, '(/a)')  &
                '*** Error: There are too many gas fractions: ' // &
                'increase npd_gas_frac and recompile.'
              ierr=i_err_fatal
              RETURN
            ENDIF
            READ(iu_file_in, *) gas_frac(n_gas_frac)
          ELSE
            l_data_region = .FALSE. 
          ENDIF
        ELSE
          IF (line(1:11) == '*BEGIN_DATA') THEN
            l_data_region= .TRUE. 
          ENDIF
        ENDIF
      ENDDO
!
      CLOSE(iu_file_in)
!
    ELSE
!     No self-broadening, set all gas fractions to zero
      gas_frac = 0.0_RealK
!
    END IF
!
!
!
  END SUBROUTINE select_self_broadening
!
!
!
  SUBROUTINE select_cont_details_int
!
!
!
    WRITE(iu_stdout, '(/a)') &
      'Enter the minimum pathlength for continuum absorption.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) umin_c
      IF (ios == 0) THEN
        EXIT
      ELSE
        CALL check_ios_int
      ENDIF
    ENDDO
!
    WRITE(iu_stdout, '(/a)') &
      'Enter the maximum pathlength for continuum absorption.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) umax_c
      IF (ios == 0) THEN
        EXIT
      ELSE
        CALL check_ios_int
      ENDIF
    ENDDO
!
    WRITE(iu_stdout, '(/a)') &
      'Enter the number of pathlengths for continuum absorption.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) n_path_c
      IF (ios == 0) THEN
        EXIT
      ELSE
        CALL check_ios_int
      ENDIF
    ENDDO
!
    WRITE(iu_stdout, '(/a)') &
      'Enter the number of partial pressures for continuum absorption.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) n_pp
      IF (ios == 0) THEN
        EXIT
      ELSE
        CALL check_ios_int
      ENDIF
    ENDDO
!
!
!
  END SUBROUTINE select_cont_details_int
!
!
!
  SUBROUTINE select_fit_type_int
!  
!  
!  
!   Selection of the style of fitting.
    WRITE(iu_stdout, '(/a, /a)') &
      'Enter the type of c-k fit required.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) i_ck_fit
      IF (ios == 0) THEN
        IF (i_ck_fit == IP_ck_tol .OR. i_ck_fit == IP_ck_bin) THEN
          WRITE(iu_stdout, '(/a, /a)') &
            'Enter the tolerance for the fit.'
          DO
            READ(iu_stdin, *, IOSTAT=ios) tol
            IF (ios == 0) THEN
              EXIT
            ELSE
              WRITE(iu_err, "(a)") 'Invalid response.'
              IF (l_interactive) THEN
                WRITE(*, "(a)") "Please re-enter."
              ELSE
                ierr=i_err_fatal
                RETURN
              ENDIF
            ENDIF
          ENDDO
        ELSE IF (i_ck_fit == IP_ck_fixed_n) THEN
!         Set a reasonable tolerance to be used in scaling function fits
          tol=1.0e-3_RealK
          WRITE(iu_stdout, '(/a, /a)') &
            'Enter the number of terms for the fit.'
          DO
            READ(iu_stdin, *, IOSTAT=ios) n_k(list_band(1))
            IF (ios == 0) THEN
              DO j=2, n_selected_band
                n_k(list_band(j)) = n_k(list_band(1))
              ENDDO
              EXIT
            ELSE
              WRITE(iu_err, "(a)") 'Invalid response.'
              IF (l_interactive) THEN
                WRITE(*, "(a)") "Please re-enter."
              ELSE
                ierr=i_err_fatal
                RETURN
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        EXIT
      ELSE
        WRITE(iu_err, "(a)") 'Invalid response.'
        IF (l_interactive) THEN
          WRITE(*, "(a)") "Please re-enter."
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
    ENDDO

    IF (i_ck_fit /= IP_ck_none) THEN
!     Setting of the maximum pathlength.
      WRITE(iu_stdout, '(/a, /a)') &
        'Enter the maximum pathlength for the absorber.'
      DO
        READ(iu_stdin, *, IOSTAT=ios) max_path
        IF (ios == 0) THEN
          EXIT
        ELSE
          CALL check_ios_int
        ENDIF
      ENDDO
    END IF
!
!
  END SUBROUTINE select_fit_type_int
!
!
!
  SUBROUTINE select_mapping_wgt
!
    WRITE(*, '(/a)') 'Do you wish to use a pre-defined mapping of ' // &
      'wavenumbers to g-space? (Y/N)'
    DO
      READ(*, '(a)') char_yn
      IF ( (char_yn == 'Y') .OR. (char_yn == 'y') ) THEN
        l_load_map=.TRUE.
        EXIT
      ELSE IF ( (char_yn == 'N') .OR. (char_yn == 'n') ) THEN
        l_load_map=.FALSE.
        EXIT
      ELSE
        WRITE(iu_err, "(/a)") "Erroneous response"
        IF (l_interactive) THEN
          WRITE(*, "(a)") "Please re-enter."
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
    ENDDO
!
    IF (.NOT.l_load_map) THEN
      WRITE(*, '(/a)') 'Do you wish to use pre-defined k-term weights? (Y/N)'
      DO
        READ(*, '(a)') char_yn
        IF ( (char_yn == 'Y') .OR. (char_yn == 'y') ) THEN
          l_load_wgt=.TRUE.
          EXIT
        ELSE IF ( (char_yn == 'N') .OR. (char_yn == 'n') ) THEN
          l_load_wgt=.FALSE.
          EXIT
        ELSE
          WRITE(iu_err, "(/a)") "Erroneous response"
          IF (l_interactive) THEN
            WRITE(*, "(a)") "Please re-enter."
          ELSE
            ierr=i_err_fatal
            RETURN
          ENDIF
        ENDIF
      ENDDO
    END IF
!
    IF (.NOT.l_load_map .AND. .NOT. l_load_wgt) THEN
      WRITE(*, '(/a)') 'Do you wish to save the mapping of ' // &
        'wavenumbers to g-space and k-term weights? (Y/N)'
      DO
        READ(*, '(a)') char_yn
        IF ( (char_yn == 'Y') .OR. (char_yn == 'y') ) THEN
          l_save_map=.TRUE.
          EXIT
        ELSE IF ( (char_yn == 'N') .OR. (char_yn == 'n') ) THEN
          l_save_map=.FALSE.
          EXIT
        ELSE
          WRITE(iu_err, "(/a)") "Erroneous response"
          IF (l_interactive) THEN
            WRITE(*, "(a)") "Please re-enter."
          ELSE
            ierr=i_err_fatal
            RETURN
          ENDIF
        ENDIF
      ENDDO
    END IF
!
    IF (l_load_map .OR. l_save_map .OR. l_load_wgt) THEN
      WRITE(iu_stdout, '(/a)') &
        'Give the name of the mapping file.'
      READ(iu_stdin, '(a)') file_map
    END IF
!
  END SUBROUTINE select_mapping_wgt
!
!
!
  SUBROUTINE check_ios_int
!
!
!
    WRITE(iu_err, "(a)") 'Invalid response.'
    IF (l_interactive) THEN
      WRITE(*, "(a)") "Please re-enter."
    ELSE
      ierr=i_err_fatal
      RETURN
    ENDIF
!
!
!
  END SUBROUTINE check_ios_int
!
!
!
END SUBROUTINE set_condition_ck_90
