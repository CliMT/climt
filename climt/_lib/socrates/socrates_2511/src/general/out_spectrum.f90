! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutines to write a spectral file.
!
SUBROUTINE out_spectrum(file_spectral, Spectrum, ierr)
!
! Method:
!   The output file is opened and a subroutine is called to
! write out each block of data.
!
!
! Modules used
  USE def_spectrum, ONLY: StrSpecData, StrSpecBasic, StrSpecSolar, &
    StrSpecRayleigh, StrSpecGas, StrSpecPlanck, StrSpecCont, StrSpecContGen, &
    StrSpecDrop, StrSpecAerosol, StrSpecIce, StrSpecVar
  USE rad_pcf, ONLY: i_normal, i_err_fatal, &
    ip_rayleigh_total, ip_rayleigh_custom, &
    ip_scale_lookup, ip_scale_t_lookup, &
    ip_scale_fnc_null, ip_scale_null, n_scale_variable, &
    ip_slingo_schrecker, ip_ackerman_stephens, ip_drop_pade_2, &
    ip_slingo_schr_phf, ip_drop_pade_2_phf, ip_ps_size_phf, &
    ip_clcmp_st_water, &
    ip_aerosol_param_dry, ip_aerosol_param_moist, name_aerosol_component, &
    ip_slingo_schrecker_ice, ip_ice_adt, ip_sun_shine_vn2_vis, &
    ip_sun_shine_vn2_ir, ip_ice_fu_solar, ip_ice_fu_ir, &
    ip_slingo_schr_ice_phf, ip_ice_fu_phf, ip_ice_adt_10, ip_ice_t_iwc, &
    ip_ice_iwc_only, ip_ice_baran, ip_ice_pade_2_phf, &
    ip_clcmp_st_ice
  USE def_std_io_icf, ONLY: iu_err
  USE missing_data_mod, ONLY: rmdi

  IMPLICIT NONE
!
!
!
! Dummy variables.
  CHARACTER (LEN=*), Intent(IN) :: file_spectral
!   Name of spectral file
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  TYPE (StrSpecData), Target :: Spectrum
!   Spectral data
!
! Local variables.
  INTEGER :: ios
!   I/O error flag
  CHARACTER (LEN=256) :: spectral_k
!   Name of extended spectral file
  CHARACTER (LEN=256) :: spectral_var
!   Name of spectral variability file
  INTEGER :: iu_spc, iu_spc1, iu_spc2
!   Unit numbers for output of the spectral data
  INTEGER :: i, ip, it
!   Loop variable
  LOGICAL :: l_exist_k
!   Existence flag for file



! Get a free unit number for the spectral file.
  CALL get_free_unit(ierr, iu_spc)
  IF (ierr /= i_normal) RETURN

! Open the file for writing.
  OPEN(UNIT=iu_spc, FILE=file_spectral, IOSTAT=ios, STATUS='UNKNOWN')
  IF (ios /= 0) THEN
    WRITE(iu_err, '(/a)') &
      '*** Error: File could not be opened.'
    ierr=i_err_fatal
    RETURN
  ENDIF

! Determine if an extended file is required and open it if so.
  l_exist_k = Spectrum%Planck%l_planck_tbl &
         .OR. Spectrum%Dim%nd_sub_band_gas > 1
  IF (ALLOCATED(Spectrum%Gas%i_scale_fnc)) THEN
    l_exist_k = l_exist_k .OR. &
     ANY(Spectrum%Gas%i_scale_fnc == ip_scale_lookup)
  END IF
  IF (ALLOCATED(Spectrum%ContGen%i_band_k_cont)) THEN
    l_exist_k = l_exist_k .OR. &
      ANY(Spectrum%ContGen%i_band_k_cont > 0)
  END IF
  IF (l_exist_k) THEN
    i=INDEX(file_spectral, ' ') - 1
    spectral_k = file_spectral(1:i) // '_k'
    CALL get_free_unit(ierr, iu_spc1)
    IF (ierr /= i_normal) RETURN
    OPEN(UNIT=iu_spc1, FILE=spectral_k, IOSTAT=ios, STATUS='UNKNOWN')
  END IF

! For each group of data which is present write out the 
! appropriate block, implicitly using the most recent type and
! version.
  IF (Spectrum%Basic%l_present(0)) CALL write_block_0_int
  IF (Spectrum%Basic%l_present(1)) &
    CALL write_block_1_int(Spectrum%Basic)
  IF (Spectrum%Basic%l_present(2)) &
    CALL write_block_2_int(Spectrum%Basic, Spectrum%Solar)
  IF (Spectrum%Basic%l_present(3)) &
    CALL write_block_3_int(Spectrum%Basic, Spectrum%Rayleigh)
  IF (Spectrum%Basic%l_present(4)) CALL write_block_4_int
  IF (Spectrum%Basic%l_present(5)) &
    CALL write_block_5_int(Spectrum%Basic, Spectrum%Gas)
  IF (Spectrum%Basic%l_present(6)) &
    CALL write_block_6_int(Spectrum%Basic, Spectrum%Planck)
! Block 7 is obsolete and omitted.
  IF (Spectrum%Basic%l_present(8)) &
    CALL write_block_8_int
  IF (Spectrum%Basic%l_present(9)) &
    CALL write_block_9_int(Spectrum%Basic, Spectrum%Cont)
  IF (Spectrum%Basic%l_present(10)) THEN
!   A separate  block is written for each type of drop.
    DO i=1, Spectrum%Dim%nd_drop_type
      IF (Spectrum%Drop%l_drop_type(i)) &
        CALL write_block_10_int(i, Spectrum%Basic, Spectrum%Drop)
    ENDDO
  ENDIF
  IF (Spectrum%Basic%l_present(11)) THEN
!   A different block is produced for each species.
    DO i=1, Spectrum%Aerosol%n_aerosol
      IF (Spectrum%Aerosol%l_aero_spec(i)) &
        CALL write_block_11_int(i, Spectrum%Basic, Spectrum%Aerosol)
    ENDDO
  ENDIF
  IF (Spectrum%Basic%l_present(12)) THEN
!   A separate  block is written for each type of crystal.
    DO i=1, Spectrum%Dim%nd_ice_type
      IF (Spectrum%Ice%l_ice_type(i)) &
        CALL write_block_12_int(i, Spectrum%Basic, Spectrum%Ice)
    ENDDO
  ENDIF
  IF (Spectrum%Basic%l_present(13)) &
    CALL write_block_13_int(Spectrum%Gas)
  IF (Spectrum%Basic%l_present(14)) &
    CALL write_block_14_int(Spectrum%Basic)
  IF (Spectrum%Basic%l_present(15)) THEN
!   A different block is produced for each species.
    DO i=1, Spectrum%Aerosol%n_aerosol
      IF (Spectrum%Aerosol%l_aero_spec(i)) &
        CALL write_block_15_int(i, Spectrum%Basic, Spectrum%Aerosol)
    ENDDO
  ENDIF
  IF (Spectrum%Basic%l_present(17)) &
    CALL write_block_17(Spectrum%Basic, Spectrum%Var)
  IF (Spectrum%Basic%l_present(18)) CALL write_block_18
  IF (Spectrum%Basic%l_present(19)) &
    CALL write_block_19(Spectrum%Basic, Spectrum%ContGen)
  IF (Spectrum%Basic%l_present(20)) CALL write_block_20

! Close the file.
  CLOSE(iu_spc)
  IF (l_exist_k) CLOSE(iu_spc1)


CONTAINS
!
!
!
  SUBROUTINE write_block_0_int

    USE gas_list_pcf

    IMPLICIT NONE


!   Local variables.
    INTEGER :: i, j
!     Loop variables
    INTEGER :: i_type, i_type_1, i_type_2
!     Identifier for gas
    INTEGER :: i_index_1, i_index_2
!     Indices of continuum gases
    INTEGER :: nd_k_term, nd_k_term_cont
!     Maximum number of k-terms in a band


    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =    0', ': SUBTYPE =    0', ': VERSION =    2'
    WRITE(iu_spc, '(a)') 'Summary of spectral data.'

    WRITE (iu_spc, '(a26, 1x, i5)') &
      'Number of spectral bands =', Spectrum%Basic%n_band
    WRITE(iu_spc, '(a35, 1x, i5)') &
      'Total number of gaseous absorbers =', Spectrum%Gas%n_absorb
    IF (ALLOCATED(Spectrum%Gas%i_band_k)) THEN
      nd_k_term = 0
      DO i=1, Spectrum%Basic%n_band
        DO j=1, Spectrum%Gas%n_band_absorb(i)
          nd_k_term = MAX(nd_k_term, &
            Spectrum%Gas%i_band_k(i, Spectrum%Gas%index_absorb(j, i)) )
        END DO
      END DO
      IF (nd_k_term > 0) WRITE(iu_spc, '(a37, 1x, i5)') &
        'Maximum number of k-terms in a band =', nd_k_term
    END IF

    IF (ALLOCATED(Spectrum%Gas%n_sub_band_gas)) THEN
      IF (ANY(Spectrum%Gas%n_sub_band_gas(1:Spectrum%Basic%n_band, &
                                          1:Spectrum%Gas%n_absorb) > 1)) THEN
        WRITE (iu_spc, '(a48, i6)') &
          'Maximum number of spectral sub-bands in a band =', &
          MAXVAL(Spectrum%Gas%n_sub_band_gas(1:Spectrum%Basic%n_band, &
                                             1:Spectrum%Gas%n_absorb))
      END IF
    END IF

    IF (Spectrum%ContGen%n_cont > 0) THEN
      WRITE(iu_spc, '(a, i5)') 'Total number of generalised continua = ', &
        Spectrum%ContGen%n_cont
      IF (ALLOCATED(Spectrum%ContGen%i_band_k_cont)) THEN
        nd_k_term_cont = 0
        DO i = 1, Spectrum%Basic%n_band
          DO j = 1, Spectrum%ContGen%n_band_cont(i)
            nd_k_term_cont = MAX(nd_k_term_cont, &
              Spectrum%ContGen%i_band_k_cont(i, &
                Spectrum%ContGen%index_cont(j, i)) )
          END DO
        END DO
        IF (nd_k_term_cont > 0) WRITE(iu_spc, '(a, i5)') &
            'Maximum number of continuum k-terms in a band = ', nd_k_term_cont
      END IF
    END IF

    IF (Spectrum%Aerosol%n_aerosol > 0) THEN
      WRITE(iu_spc, '(a26, 1x, i5)') &
        'Total number of aerosols =', Spectrum%Aerosol%n_aerosol
    END IF

    IF (Spectrum%Photol%n_pathway > 0) THEN
      WRITE(iu_spc, '(a37, 1x, i5)') &
        'Total number of photolysis pathways =', Spectrum%Photol%n_pathway
    END IF

    WRITE(iu_spc, '(a)') 'List of indexing numbers and absorbers.'
    WRITE(iu_spc, '(a)') &
      'Index       Absorber(identifier and name)'
    DO i=1, Spectrum%Gas%n_absorb
      i_type = Spectrum%Gas%type_absorb(i)
      WRITE(iu_spc, '(i5, 7x, i5, 7x, a)') &
        i, i_type, name_absorb(i_type)
    END DO

    IF (Spectrum%ContGen%n_cont > 0) THEN
      WRITE(iu_spc, '(a)') 'Listing of continuum indexing numbers and gases.'
      WRITE(iu_spc, '(a)') 'Index     ' // &
        'Gas 1(index and name)              ' // &
        'Gas 2(index and name)'
      DO j = 1, Spectrum%ContGen%n_cont
        i_index_1 = Spectrum%ContGen%index_cont_gas_1(j)
        i_index_2 = Spectrum%ContGen%index_cont_gas_2(j)
        i_type_1 = Spectrum%Gas%type_absorb(i_index_1)
        i_type_2 = Spectrum%Gas%type_absorb(i_index_2)
        WRITE(iu_spc, '(i5, 5x, i5, 5x, a, 5x, i5, 5x, a)') j, &
          i_index_1, name_absorb(i_type_1), &
          i_index_2, name_absorb(i_type_2)
      END DO
    END IF

    IF (Spectrum%Aerosol%n_aerosol > 0) THEN
      WRITE(iu_spc, '(a)') 'List of indexing numbers of aerosols.'
      WRITE(iu_spc, '(a)') &
        'Index       Aerosol(type number and name)'
      DO i=1, Spectrum%Aerosol%n_aerosol
        i_type=Spectrum%Aerosol%type_aerosol(i)
        WRITE(iu_spc, '(i5, 7x, i5, 7x, a)') &
          i, i_type, name_aerosol_component(i_type)
      END DO
    END IF

    IF (Spectrum%Photol%n_pathway > 0) THEN
      WRITE(iu_spc, '(a)') 'List of photolysis pathways.'
      WRITE(iu_spc, '(a)') &
       'Index  Absorber index  Reaction (Thermalisation indicator and Products)'
      DO i=1, Spectrum%Photol%n_pathway
        IF (Spectrum%Photol%pathway_products(i) > 0) THEN
          IF (Spectrum%Photol%l_thermalise(i)) THEN
            WRITE(iu_spc, '(2(i5, 7x), i5, l5, 2x, a)') i, &
              Spectrum%Photol%pathway_absorber(i), &
              Spectrum%Photol%pathway_products(i), &
              Spectrum%Photol%l_thermalise(i), &
              TRIM(photol_products(Spectrum%Photol%pathway_products(i), &
              Spectrum%Gas%type_absorb(Spectrum%Photol%pathway_absorber(i))))
          ELSE
            WRITE(iu_spc, '(3(i5, 7x), a)') i, &
              Spectrum%Photol%pathway_absorber(i), &
              Spectrum%Photol%pathway_products(i), &
              TRIM(photol_products(Spectrum%Photol%pathway_products(i), &
              Spectrum%Gas%type_absorb(Spectrum%Photol%pathway_absorber(i))))
          END IF
        ELSE
          IF (Spectrum%Photol%l_thermalise(i)) THEN
            WRITE(iu_spc, '(2(i5, 7x), i5, l5)') i, &
              Spectrum%Photol%pathway_absorber(i), &
              Spectrum%Photol%pathway_products(i), &
              Spectrum%Photol%l_thermalise(i)
          ELSE
            WRITE(iu_spc, '(3(i5, 7x))') i, &
              Spectrum%Photol%pathway_absorber(i), &
              Spectrum%Photol%pathway_products(i)
          END IF
        END IF
      END DO
    END IF

    WRITE(iu_spc, '(a4)') '*END'

  END SUBROUTINE write_block_0_int
!
!
!
  SUBROUTINE write_block_1_int(SpBasic) 
!
!
    TYPE (StrSpecBasic) :: SpBasic
!
!   Local variables.
    INTEGER :: i
!     Loop variable
!
!
    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =    1', ': SUBTYPE =    0', ': VERSION =    0'
    WRITE(iu_spc, '(a)') 'Specification of spectral intervals'
    WRITE(iu_spc, '(a)') &
      'Limits of spectral intervals (wavelengths in m.)'
    WRITE(iu_spc, '(a4, 8x, a11, 9x, a11)') &
      'Band', 'Lower limit', 'Upper limit'
!
    DO i=1, SpBasic%n_band
      WRITE(iu_spc, '(i5, 3x, 2(4x, 1pe16.9))' ) &
        i, SpBasic%wavelength_short(i), SpBasic%wavelength_long(i)
    ENDDO
!
    WRITE(iu_spc, '(a4)') '*END'
!
!
!
  END SUBROUTINE write_block_1_int
!
!
!
  SUBROUTINE write_block_2_int(SpBasic, SpSolar)

    TYPE  (StrSpecBasic) :: SpBasic
    TYPE  (StrSpecSolar) :: SpSolar

!   Local variables
    INTEGER :: i
!     Loop variable
    INTEGER :: version
!     Block version

    version = 0
    IF (ALLOCATED(SpSolar%weight_blue)) THEN
      IF (SpSolar%weight_blue(1) /= rmdi) version = 1
    END IF
    WRITE(iu_spc, '(a19, a16, a12, i4)') &
        '*BLOCK: TYPE =    2', ': SUBTYPE =    0', ': VERSION = ', version

    IF (version == 0) THEN
      WRITE(iu_spc, '(a  )') &
        'Normalized solar flux in each spectral interval.'
      WRITE(iu_spc, '(a4, 8x, a15)') 'Band', 'Normalized flux'
      
      DO i=1, SpBasic%n_band
        WRITE(iu_spc, '(i5, 7x, 1pe16.9)') i, SpSolar%solar_flux_band(i)
      ENDDO
    ELSE IF (version == 1) THEN
      WRITE(iu_spc, '(a  )') &
        'Normalized solar flux and blue fraction in each spectral interval.'
      WRITE(iu_spc, '(a4, 8x, a15, 5x, a11)') &
        'Band', 'Normalized flux', 'Weight blue'
      
      DO i=1, SpBasic%n_band
        WRITE(iu_spc, '(i5, 7x, 1pe16.9, 4x, 1pe16.9)') &
          i, SpSolar%solar_flux_band(i), SpSolar%weight_blue(i)
      ENDDO
    END IF

    WRITE(iu_spc, '(a4)') '*END'

  END SUBROUTINE write_block_2_int
!
!
!
  SUBROUTINE write_block_3_int(SpBasic, SpRayleigh)
!
!
!
    TYPE  (StrSpecBasic) :: SpBasic
    TYPE  (StrSpecRayleigh) :: SpRayleigh
!
!   Local variables.
    INTEGER :: i
!     Loop variable
!
!
    IF (SpRayleigh%i_rayleigh_scheme == ip_rayleigh_total) THEN
      WRITE(iu_spc, '(a19, a16, a16)') &
        '*BLOCK: TYPE =    3', ': SUBTYPE =    0', ': VERSION =    0'
      WRITE(iu_spc, '(a)') &
        'Rayleigh mass scattering coefficients at STP: unit m**2/kg'
      WRITE(iu_spc, '(a4, 8x, a20)') 'Band', 'Rayleigh coefficient'
      WRITE(iu_spc, '(17x, a5)') 'm2/kg'
!
      DO i=1, SpBasic%n_band
        WRITE(iu_spc, '(i5, 7x, 1pe16.9)') &
          i, SpRayleigh%rayleigh_coeff(i)
      ENDDO
!
    ELSE IF (SpRayleigh%i_rayleigh_scheme == ip_rayleigh_custom) THEN
      WRITE(iu_spc, '(a19, a16, a16)') &
        '*BLOCK: TYPE =    3', ': SUBTYPE =    1', ': VERSION =    0'
      WRITE(iu_spc, '(a)') &
        'Rayleigh mass scattering coefficients at STP: unit m**2/kg'
      WRITE(iu_spc, '(a37, 1x, i5)') &
        'Number of Rayleigh scattering gases =', SpRayleigh%n_gas_rayleigh
      WRITE(iu_spc, '(a)') &
        'Indexing numbers of gases'
      WRITE(iu_spc, '(15(2x, i3))') &
        SpRayleigh%index_rayleigh(1:SpRayleigh%n_gas_rayleigh)
!
      DO i=1, SpBasic%n_band
        WRITE(iu_spc, '(a7, i5, 5x, a45, /, (4(3x, 1pe16.9)))') &
          'Band = ', i, 'Rayleigh scattering ' // &
          'coefficient for each gas:', &
          SpRayleigh%rayleigh_coeff_gas(1:SpRayleigh%n_gas_rayleigh,i)
      ENDDO
!
    ELSE
      WRITE(iu_err, '(/a)') &
        '*** Error: Rayleigh scattering scheme not recognised.'
      ierr=i_err_fatal
      RETURN
!
    END IF
!
    WRITE(iu_spc, '(a4)') '*END'
!
!
!
  END SUBROUTINE write_block_3_int
!
!
!
  SUBROUTINE write_block_4_int
!
!
!
!   Local variables
    INTEGER :: i
!     Loop variable
    INTEGER :: j
!     Loop variable
!
!
!
    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =    4', ': SUBTYPE =    0', ': VERSION =    1'
    WRITE(iu_spc, '(a)') 'Gaseous absorbers in each interval'
    WRITE(iu_spc, '(a, /A)') &
      '(The number is the indexing number of the species as set out', &
      ' in the summary block 0.)'
    WRITE(iu_spc, '(a)') &
      'A zero indicates that there is no gaseous absorption ' // &
      'in the interval.'
    WRITE(iu_spc, '(a, a, a)') 'Band,', &
      ' number of absorbers and overlap method,', &
      ' followed by indexing numbers'

    DO i=1, Spectrum%Basic%n_band
      WRITE(iu_spc, '(i5, 2(7x, i5))') &
        i, Spectrum%Gas%n_band_absorb(i), Spectrum%Gas%i_overlap(i)
      IF (Spectrum%Gas%n_band_absorb(i) > 0) THEN
        WRITE(iu_spc, '(5x, 12(2x, i3))') &
          Spectrum%Gas%index_absorb(1:Spectrum%Gas%n_band_absorb(i), i)
      ENDIF
    ENDDO

    WRITE(iu_spc, '(a4)') '*END'
!
!
!
  END SUBROUTINE write_block_4_int
!
!
!
!
  SUBROUTINE write_block_5_int(SpBasic, SpGas)
!
!
!
    TYPE  (StrSpecBasic) :: SpBasic
    TYPE  (StrSpecGas) :: SpGas

!   Local variables.
    INTEGER :: i, j, k, igf, isb
!     Loop variables
    INTEGER :: i_index
!     Index of gas
    INTEGER :: i_index_sb
!     Index of gas in arrays with self-broadening


    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =    5', ': SUBTYPE =    0', ': version =    1'
    WRITE(iu_spc, '(a)') &
      'Exponential sum fiting coefficients: (exponents: m2/kg)'

    IF (ANY(SpGas%l_self_broadening)) THEN
      WRITE(iu_spc,'(a)') 'Self-broadened indexing numbers of all absorbers.'
      WRITE(iu_spc,'(a)') '(Zero means no self-broadening.)'
      WRITE(iu_spc,'(6(2x, i3))') SpGas%index_sb(1:SpGas%n_absorb)
    END IF

    WRITE(iu_spc, '(a4, 8x, a58, /, 12x, a48)') 'Band', &
      'Gas, Number of k-terms, Scaling type and scaling function,', &
      ' followed by reference pressure and temperature,'
    WRITE(iu_spc, '(12x, a47)') &
      ' k-terms, weights and scaling parameters.'

    DO i=1, SpBasic%n_band
      DO j=1, SpGas%n_band_absorb(i)
        i_index=SpGas%index_absorb(j, i)
        IF (SpGas%i_scale_fnc(i, i_index) == ip_scale_fnc_null) THEN
          SpGas%i_scale_k(i, i_index) = ip_scale_null
        END IF
        WRITE(iu_spc, '(i5, 4(7x, i5))') &
          i, i_index, SpGas%i_band_k(i, i_index), &
          SpGas%i_scale_k(i, i_index), SpGas%i_scale_fnc(i, i_index)
        IF (SpGas%i_scale_fnc(i, i_index) == ip_scale_lookup) THEN
          WRITE(iu_spc, '(2(7x, i5))') &
            SpGas%num_ref_p(i_index, i), SpGas%num_ref_t(i_index, i)
        ELSE IF (SpGas%i_scale_fnc(i, i_index) == ip_scale_t_lookup) THEN
          WRITE(iu_spc, '(a14, i5, a13)') 'Lookup table: ', &
            SpGas%num_ref_t(i_index, i), ' temperatures'
        ELSE
          WRITE(iu_spc, '(2(6x, 1pe16.9))') &
            SpGas%p_ref(i_index, i), SpGas%t_ref(i_index, i)
        END IF
        DO k=1, SpGas%i_band_k(i, i_index)
          IF (SpGas%i_scat(k, i, i_index)==0) THEN
            WRITE(iu_spc, '(2(3x, 1pe16.9), (t39, 2(3x, 1pe16.9)))') &
              SpGas%k(k, i, i_index), SpGas%w(k, i, i_index), &
              SpGas%scale(1:n_scale_variable(SpGas%i_scale_fnc(i, i_index)), &
                k, i, i_index)
          ELSE
            WRITE(iu_spc, '(2(3x, 1pe16.9),i3,(t42, 1pe16.9,3x,1pe16.9))') &
              SpGas%k(k, i, i_index), SpGas%w(k, i, i_index), &
              SpGas%i_scat(k, i, i_index), SpGas%scale(1:n_scale_variable &
                (SpGas%i_scale_fnc(i, i_index)), k, i, i_index)
          END IF
        END DO
      END DO
    END DO

    WRITE(iu_spc, '(a4)') '*END'


    IF (ANY(SpGas%i_scale_fnc == ip_scale_lookup)) THEN
      WRITE(iu_spc1,'(a6,a2,a12)') '*BLOCK', ': ', 'k-table     '
      WRITE(iu_spc1,'(/,2(a,i4),a)') 'Lookup table: ', &
        Spectrum%Dim%nd_pre, ' pressures, ', &
        Spectrum%Dim%nd_tmp, ' temperatures.'
      DO ip=1, Spectrum%Dim%nd_pre
        WRITE(iu_spc1,'(6(1PE13.6))') EXP(SpGas%p_lookup(ip)), &
          SpGas%t_lookup(1:Spectrum%Dim%nd_tmp,ip)
      END DO
      IF (ANY(SpGas%l_self_broadening)) THEN
        WRITE(iu_spc1,'(/,a,i4,a)') 'Lookup table: ', &
          SpGas%n_gas_frac, ' gas fractions.'
        WRITE(iu_spc1,'(6(1PE13.6))') &
          SpGas%gf_lookup(1:SpGas%n_gas_frac)
      END IF
      DO i=1, SpBasic%n_band
        DO j=1, SpGas%n_band_absorb(i)
          i_index=SpGas%index_absorb(j, i)
          IF (SpGas%i_scale_fnc(i, i_index) == ip_scale_lookup ) THEN
            WRITE(iu_spc1,'(/,3(a,i4))') 'Band: ',i,', gas: ',i_index, &
              ', k-terms: ',SpGas%i_band_k(i, i_index)
            IF (SpGas%l_self_broadening(i_index)) THEN
              i_index_sb = SpGas%index_sb(i_index)
              DO k=1, SpGas%i_band_k(i, i_index)
                DO igf=1, SpGas%n_gas_frac
                  DO ip=1, Spectrum%Dim%nd_pre
                    WRITE(iu_spc1,'(6(1PE13.6))') &
                      SpGas%k_lookup_sb(:,ip,igf,k,i_index_sb,i)
                  END DO
                END DO
              END DO
            ELSE
              DO k=1, SpGas%i_band_k(i, i_index)
                DO ip=1, Spectrum%Dim%nd_pre
                  WRITE(iu_spc1,'(6(1PE13.6))') &
                    SpGas%k_lookup(:,ip,k,i_index,i)
                END DO
              END DO
            END IF
          END IF
        END DO
      END DO
      WRITE(iu_spc1, '(a4)') '*END'
    END IF


    IF (ANY(SpGas%i_scale_fnc == ip_scale_t_lookup)) THEN
      WRITE(iu_spc1,'(//,a6,a2,a20)') '*BLOCK', ': ', 'temperature k-tables'
      WRITE(iu_spc1,'(/,a9,i4,a)') 'Maximum: ', &
        Spectrum%Dim%nd_t_lookup_gas, ' temperatures.'
      DO i_index=1, SpGas%n_absorb
        IF (ANY(SpGas%i_scale_fnc(:, i_index) == ip_scale_t_lookup)) THEN
          WRITE(iu_spc1,'(/,2(a,i4),a)') 'Lookup table for gas ', i_index, &
            ': ', SpGas%n_t_lookup_gas(i_index), ' temperatures.'
          WRITE(iu_spc1,'(6(1PE13.6))') &
            SpGas%t_lookup_gas(1:SpGas%n_t_lookup_gas(i_index), i_index)
        END IF
      END DO
      DO i=1, SpBasic%n_band
        DO j=1, SpGas%n_band_absorb(i)
          i_index=SpGas%index_absorb(j, i)
          IF (SpGas%i_scale_fnc(i, i_index) == ip_scale_t_lookup ) THEN
            WRITE(iu_spc1,'(/,3(a,i4))') 'Band: ',i,', gas: ',i_index, &
              ', k-terms: ',SpGas%i_band_k(i, i_index)
            DO k=1, SpGas%i_band_k(i, i_index)
              WRITE(iu_spc1,'(6(1PE13.6))') &
                SpGas%k_t_lookup_gas(1:SpGas%n_t_lookup_gas(i_index), &
                                     k, i_index, i)
            END DO
          END IF
        END DO
      END DO
      WRITE(iu_spc1, '(a4)') '*END'
    END IF


    IF (ANY(SpGas%n_sub_band_gas(1:SpBasic%n_band, 1:SpGas%n_absorb) > 1)) THEN
      WRITE(iu_spc1,'(//,a6,a2,a16)') '*BLOCK', ': ', 'sub-band mapping'
      DO i=1, SpBasic%n_band
        DO j=1, SpGas%n_band_absorb(i)
          i_index=SpGas%index_absorb(j, i)
          IF (SpGas%n_sub_band_gas(i, i_index) > 1) THEN
            WRITE(iu_spc1, '(/,2(a,i4),a,i6)') 'Band: ',i,', gas: ',i_index, &
              ', sub-bands: ', SpGas%n_sub_band_gas(i, i_index)
            WRITE(iu_spc1, '(a8,2x,a6,7x,a6,7x,a16,3x,a15)') &
              'Sub-band','k-term','weight','wavelength_short','wavelength_long'
            DO isb=1, SpGas%n_sub_band_gas(i, i_index)
              WRITE(iu_spc1, '(2i8, 3(2x,1PE16.9))') isb, &
                SpGas%sub_band_k(isb, i, i_index), &
                SpGas%sub_band_w(isb, i, i_index), &
                SpGas%wavelength_sub_band(:, isb, i, i_index)
            END DO
          END IF
        END DO
      END DO
      WRITE(iu_spc1, '(a4)') '*END'
    END IF
!
!
!
  END SUBROUTINE write_block_5_int
!
!
!
  SUBROUTINE write_block_6_int(SpBasic, SpPlanck)
!
!
!
    TYPE  (StrSpecBasic) :: SpBasic
    TYPE  (StrSpecPlanck) :: SpPlanck
!
!   Local variables.
    INTEGER :: i
!     Loop variable
!
!
!
    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =    6', ': SUBTYPE =    0', ': version =    1'
    WRITE(iu_spc, '(a)') 'Thermal source function.'
    IF (SpPlanck%l_planck_tbl) THEN
      WRITE(iu_spc, '(a14, a10)') 'Type of data: ', 'table     '
      
      WRITE(iu_spc, '(a23, 1x, i5, a1, 1x, a23, 1x, 1pe16.9)') &
        'Number of data points =', SpPlanck%n_deg_fit+1, ':', &
        'Reference temperature =', SpPlanck%t_ref_planck
!
!     Write Planck function table to extended spectral file
      WRITE(iu_spc1,'(a6, a2, a12)') '*BLOCK', ': ', 'Planck table'
      WRITE(iu_spc1,'(/,5(1PE16.8E3))') &
        SpPlanck%theta_planck_tbl(0:SpPlanck%n_deg_fit)
      DO i=1, SpBasic%n_band 
        WRITE(iu_spc1,'(/,3(a,i4))') 'Band: ',i
        WRITE(iu_spc1,'(5(1PE16.8E3))') &
          SpPlanck%thermal_coeff(0:SpPlanck%n_deg_fit,i)
      END DO
      WRITE(iu_spc1, '(a4)') '*END'
    ELSE
      WRITE(iu_spc, '(a14, a10)') 'Type of data: ', 'polynomial'
!
      WRITE(iu_spc, '(a22, 2x, i5, a1, 1x, a23, 1x, 1pe16.9)') &
        'Degree of polynomial =', SpPlanck%n_deg_fit, ':', &
        'Reference temperature =', SpPlanck%t_ref_planck
      WRITE(iu_spc, '(a)') 'Coefficients in each band.'
      WRITE(iu_spc, '(a4, 8x, 2x, a)') 'Band', 'w/m2'
      DO i=1, SpBasic%n_band
        WRITE(iu_spc, '(i5, 7x, (t13, 3(1pe16.9, 4x)))') &
          i, SpPlanck%thermal_coeff(0:SpPlanck%n_deg_fit, i)
      ENDDO
    END IF
!
    WRITE(iu_spc, '(a4)') '*END'
!
!
!
  END SUBROUTINE write_block_6_int
!
!
!
  SUBROUTINE write_block_8_int
!
!      
!
!   Local variables
    INTEGER :: i
!     Loop variable
    INTEGER :: j
!     Loop variable
!
!
!
    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =    8', ': SUBTYPE =    0', ': VERSION =    0'
    WRITE(iu_spc, '(a)') 'Continuum absorbers in each interval'
    WRITE(iu_spc, '(a)') &
      '(The number is the indexing number of each type as set out'
    WRITE(iu_spc, '(a)') &
      ' in the module rad_pcf.)'
    WRITE(iu_spc, '(a)') &
      'A zero indicates that there is no continuum absorption ' // &
      'in the interval.'
    WRITE(iu_spc, '(a4, 8x, a55)') &
      'Band', &
      'Number of active continua followed by indexing numbers'
!
    DO i=1, Spectrum%Basic%n_band
      WRITE(iu_spc, '(i5, 7x, i5)') i, Spectrum%Cont%n_band_continuum(i)
      IF (Spectrum%Cont%n_band_continuum(i) > 0) THEN
        WRITE(iu_spc, '(5x, 4(2x, i3))') &
          (Spectrum%Cont%index_continuum(i, j), &
          j=1, Spectrum%Cont%n_band_continuum(i))
      ENDIF
    ENDDO
!
    WRITE(iu_spc, '(a)') 'Indexing numbers of gases for continua:'
    WRITE(iu_spc, '(5x, a16, 1x, i5)') &
      'Index of water =', Spectrum%Cont%index_water
!
    WRITE(iu_spc, '(a4)') '*END'
!
!
!
  END SUBROUTINE write_block_8_int
!
!
!
  SUBROUTINE write_block_9_int(SpBasic, SpCont)
!
!
!
    TYPE  (StrSpecBasic) :: SpBasic
    TYPE  (StrSpecCont) :: SpCont
!
!   Local variables.
    INTEGER :: i
!     Loop variable
    INTEGER :: j
!     Loop variable
    INTEGER :: l
!     Loop variable
    INTEGER :: i_index
!     Index of gas
!
!
!
    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =    9', ': SUBTYPE =    0', ': VERSION =    0'
    WRITE(iu_spc, '(a)') &
      'Continuum absorption coefficients: (units: m5/(kg.mol))'
    WRITE(iu_spc, '(a4, 8x, a9, 3x, a51)') &
      'Band', 'Continuum', &
      'scaling function followed by reference pressure and'
    WRITE(iu_spc, '(12x, a48)') &
      ' Temperature, Extinction, and Scaling parameters.'
!
    DO i=1, SpBasic%n_band
      DO j=1, SpCont%n_band_continuum(i)
        i_index=SpCont%index_continuum(i, j)
        WRITE(iu_spc, '(i5, 2(7x, i5))') &
          i, i_index, SpCont%i_scale_fnc_cont(i, i_index)
        WRITE(iu_spc, '(2(6x, 1pe16.9))') &
          SpCont%p_ref_cont(i_index, i), SpCont%t_ref_cont(i_index, i)
        WRITE(iu_spc, '(6x, 1pe16.9, (t23, 2(6x, 1pe16.9)))') &
          SpCont%k_cont(i, i_index), &
          SpCont%scale_cont(1:n_scale_variable(SpCont%i_scale_fnc_cont &
                                              (i, i_index)), i, i_index)
      ENDDO
    ENDDO
!
    WRITE(iu_spc, '(a4)') '*END'
!
!
!
  END SUBROUTINE write_block_9_int
!
!
!
  SUBROUTINE write_block_10_int(i_d, SpBasic, SpDrop)
!
!
!
!   Dummy variables:
    INTEGER :: i_d
!     Type of droplet
    TYPE  (StrSpecBasic), Intent(IN) :: SpBasic
    TYPE  (StrSpecDrop), Intent(IN) :: SpDrop
!
!
!
!   Local variables.
    INTEGER :: n_parameter
!     Number of parameters
    INTEGER :: i
!     Loop variable
    INTEGER :: k
!     Loop variable

!   Functions called:
    INTEGER, EXTERNAL :: set_n_cloud_parameter


    IF ( (SpDrop%i_drop_parm(i_d) == IP_slingo_schrecker) .OR. &
         (SpDrop%i_drop_parm(i_d) == IP_ackerman_stephens) .OR. &
         (SpDrop%i_drop_parm(i_d) == IP_drop_pade_2).OR. &
         (SpDrop%i_drop_parm(i_d) == ip_slingo_schr_phf).OR. &
         (SpDrop%i_drop_parm(i_d) == IP_drop_pade_2_PHF).OR. &
         (SpDrop%i_drop_parm(i_d) == IP_ps_size_PHF) ) THEN
!
      WRITE(iu_spc, '(a19, a16, a16)') &
        '*BLOCK: TYPE =   10', ': SUBTYPE =    0', ': VERSION =    2'
!
!     Calculate the number of parameters for the scheme.
      n_parameter = set_n_cloud_parameter( &
        SpDrop%i_drop_parm(i_d), ip_clcmp_st_water, SpDrop%n_phf(i_d))
!
      WRITE(iu_spc, '(a42)') &
        'Parametrized scattering data for droplets.'
      WRITE(iu_spc, '(a27, i5)') &
        'Type number of droplets = ', i_d
      WRITE(iu_spc, '(a34, i5, a27, i5)') &
        'Index of parametrization scheme = ', SpDrop%i_drop_parm(i_d), &
        ':   Number of parameters = ', n_parameter
      WRITE(iu_spc, '(a42, i5)') &
        'Number of moments of the phase function = ', SpDrop%n_phf(i_d)
      WRITE(iu_spc, '(a39, 1pe12.5, a4, 1pe12.5)') &
        'Range of validity of parametrization = ', &
        SpDrop%parm_min_dim(i_d), ' -- ', SpDrop%parm_max_dim(i_d)
!
!     Write out the scattering parameters in each band.
      DO i=1, SpBasic%n_band
        WRITE(iu_spc, '(a7, i5, 5x, a19, /, (4(4x, 1pe12.5)))') &
          'Band = ', i, 'Fitting parameters:', &
          SpDrop%parm_list(1: n_parameter, i, i_d)
      ENDDO
!
    ENDIF
!
    WRITE(iu_spc, '(a4)') '*END'
!
!
!
  END SUBROUTINE write_block_10_int
!
!
!
  SUBROUTINE write_block_11_int(i_a, SpBasic, SpAerosol)
!
!
!
    INTEGER :: i_a
!     Species of aerosol
    TYPE  (StrSpecBasic), Intent(IN) :: SpBasic
    TYPE  (StrSpecAerosol), Intent(IN) :: SpAerosol
!
!   Local variables.
    INTEGER :: i
!     Loop variable
    INTEGER :: k
!     Loop variable
    INTEGER :: l
!     Loop variable
!
!
!
    IF (SpAerosol%i_aerosol_parm(i_a) == IP_aerosol_param_dry) THEN
!
      WRITE(iu_spc, '(a19, a16, a16)') &
        '*BLOCK: TYPE =   11', ': SUBTYPE =    0', ': VERSION =    2'
      WRITE(iu_spc, '(a39)') &
        'Scattering parameters for dry aerosols.'
      WRITE(iu_spc, '(a19, i5, 2x, a20)') &
        'Index of species = ' , i_a, &
        name_aerosol_component(SpAerosol%type_aerosol(i_a))
      WRITE(iu_spc, '(a36, i5)') &
        'Number of terms in phase function = ', &
        SpAerosol%n_aerosol_phf_term(i_a)
      WRITE(iu_spc, '(a4, 8x, 2(a10, 10x), a9)') &
        'Band', 'Absorption', 'Scattering', 'Phase fnc.'
      WRITE(iu_spc, '(12x, a9, 11x, a9)') '(m2.kg-1)', '(m2.kg-1)'
!
!     Write out the scattering parameters for each band
      DO i=1, SpBasic%n_band
        WRITE(iu_spc, fmt='(i5, 3(4x, 1pe16.9))') &
          i, SpAerosol%abs(1, i_a, i) , &
             SpAerosol%scat(1, i_a, i), &
             SpAerosol%phf_fnc(1, 1, i_a, i)
        IF (SpAerosol%n_aerosol_phf_term(i_a) > 1) THEN
          WRITE(iu_spc, fmt='((T50, 1PE16.9))') &
            (SpAerosol%phf_fnc(1, l, i_a, i), &
              l=2, SpAerosol%n_aerosol_phf_term(i_a))
        ENDIF
      ENDDO
!
    ELSE IF (SpAerosol%i_aerosol_parm(i_a) == IP_aerosol_param_moist) THEN
!
      WRITE(iu_spc, '(a19, a16, a16)') &
        '*BLOCK: TYPE =   11', ': SUBTYPE =    1', ': VERSION =    2'
      WRITE(iu_spc, '(a41)') &
        'Scattering parameters for moist aerosols.'
      WRITE(iu_spc, '(a19, i5, 2x, a20)') &
        'Index of species = ', i_a, &
        name_aerosol_component(SpAerosol%type_aerosol(i_a))
      WRITE(iu_spc, '(a28, i3)') &
        'Number of humidity values = ', SpAerosol%nhumidity(i_a)
      WRITE(iu_spc, '(a36, i5)') &
        'Number of terms in phase function = ', &
        SpAerosol%n_aerosol_phf_term(i_a)
!
!     Write out the scattering parameters for each band
      DO i=1, SpBasic%n_band
!
        WRITE(iu_spc, '(a7, i4)') 'Band = ', i
        WRITE(iu_spc, '(8x, a8, 12x, a10, 9x, a10, 10x, a10)') &
          'Humidity', 'Absorption', 'Scattering', 'Phase fnc.'
        WRITE(iu_spc, '(28x, a9, 10x, a9)') '(m2.kg-1)', '(m2.kg-1)'
!
        DO k=1, SpAerosol%nhumidity(i_a)
          WRITE(iu_spc, fmt='(3x, 1pe16.9, 3(4x, 1pe16.9))') &
            SpAerosol%humidities(k, i_a), &
            SpAerosol%abs(k, i_a, i), &
            SpAerosol%scat(k, i_a, i), &
            SpAerosol%phf_fnc(k, 1, i_a, i)
          IF (SpAerosol%n_aerosol_phf_term(i_a) > 1) THEN
            WRITE(iu_spc, fmt='((t63, 1pe16.9))') &
              (SpAerosol%phf_fnc(k, l, i_a, i), &
              l=2, SpAerosol%n_aerosol_phf_term(i_a))
          ENDIF
        ENDDO
!
      ENDDO
!
    ENDIF
!
    WRITE(iu_spc, '(a4)') '*END'
!
!
!
  END SUBROUTINE write_block_11_int
!
!
!
  SUBROUTINE write_block_12_int(i_ice, SpBasic, SpIce)
!
!
!
    INTEGER, Intent(IN) :: i_ice
!     Type of ice crystal
    TYPE  (StrSpecBasic), Intent(IN) :: SpBasic
    TYPE  (StrSpecIce), Intent(IN) :: SpIce
!
!
!
!   Local variables.
    INTEGER :: n_parameter
!     Number of parameters
    INTEGER :: i
!     Loop variable

!   Functions called:
    INTEGER, EXTERNAL :: set_n_cloud_parameter

    IF ( (SpIce%i_ice_parm(i_ice) == IP_slingo_schrecker_ice) .OR. &
         (SpIce%i_ice_parm(i_ice) == IP_ice_adt) .OR. &
         (SpIce%i_ice_parm(i_ice) == ip_sun_shine_vn2_vis) .OR. &
         (SpIce%i_ice_parm(i_ice) == ip_sun_shine_vn2_ir) .OR. &
         (SpIce%i_ice_parm(i_ice) == IP_ice_fu_solar) .OR. &
         (SpIce%i_ice_parm(i_ice) == IP_ice_fu_ir) .OR. &
         (SpIce%i_ice_parm(i_ice) == IP_slingo_schr_ice_phf) .OR. &
         (SpIce%i_ice_parm(i_ice) == IP_ice_fu_phf) .OR. &
         (SpIce%i_ice_parm(i_ice) == IP_ice_adt_10) .OR. &
         (SpIce%i_ice_parm(i_ice) == ip_ice_t_iwc) .OR. &
         (SpIce%i_ice_parm(i_ice) == ip_ice_iwc_only) .OR. &
         (SpIce%i_ice_parm(i_ice) == ip_ice_baran) .OR. &
         (SpIce%i_ice_parm(i_ice) == ip_ice_pade_2_phf) ) THEN
!
      WRITE(iu_spc, '(a19, a16, a16)') &
        '*BLOCK: TYPE =   12', ': SUBTYPE =    0', ': VERSION =    2'
!
!     Calculate the number of parameters for the scheme.
      n_parameter = set_n_cloud_parameter( &
        SpIce%i_ice_parm(i_ice), ip_clcmp_st_ice, SpIce%n_phf(i_ice))
!
      WRITE(iu_spc, '(a46)') &
        'Parametrized scattering data for ice crystals.'
      WRITE(iu_spc, '(a31, i5)') &
        'Type number of ice crystals = ', i_ice
      WRITE(iu_spc, '(a34, i5, a27, i5)') &
        'Index of parametrization scheme = ', &
        SpIce%i_ice_parm(i_ice), &
        ':   number of parameters = ', n_parameter
      WRITE(iu_spc, '(a42, i5)') &
        'Number of moments of the phase function = ', &
        SpIce%n_phf(i_ice)
      WRITE(iu_spc, '(a39, 1pe12.5, a4, 1pe12.5)') &
        'Range of validity of parametrization = ', &
        SpIce%parm_min_dim(i_ice), ' -- ', &
        SpIce%parm_max_dim(i_ice)
!
!     Write out the scattering parameters in each band.
      DO i=1, SpBasic%n_band
        WRITE(iu_spc, '(a7, i5, 5x, a19, /, (4(4x, 1pe12.5)))') &
          'Band = ', i, 'Fitting parameters:', &
          SpIce%parm_list(1:n_parameter, i, i_ice)
      ENDDO
!
    ENDIF
!
    WRITE(iu_spc, '(a4)') '*END'
!
!
!
  END SUBROUTINE write_block_12_int
!
!
!
  SUBROUTINE write_block_13_int(SpGas)

    TYPE  (StrSpecGas), Intent(IN) :: SpGas

    INTEGER :: i
!     Loop variable


    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =   13', ': SUBTYPE =    0', ': VERSION =    0'
    WRITE(iu_spc, '(a)') &
      'Doppler broadening coefficients for the gases'
    WRITE(iu_spc, '(a)') &
      'Species    Doppler flag   Doppler coefficient'
!   Write out the doppler terms for each species.
    DO i=1, SpGas%n_absorb
      WRITE(iu_spc, '(i5, 7x, l6, 9x, 1pe12.5)') &
        i, SpGas%l_doppler(i), SpGas%doppler_cor(i)
    ENDDO
    WRITE(iu_spc, '(a4)') '*END'

  END SUBROUTINE write_block_13_int



  SUBROUTINE write_block_14_int(SpBasic)
!
!
!
    TYPE  (StrSpecBasic), Intent(IN) :: SpBasic
!
!   Local variables
    INTEGER :: i
!     Loop variable
!
!
!
    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =   14', ': SUBTYPE =    0', ': VERSION =    0'
    WRITE(iu_spc, '(a)') 'Regions excluded from specific bands'
    WRITE(iu_spc, '(a)') &
      'A zero indicates that no bands have been excluded from the ' // &
      'particular band.'
    WRITE(iu_spc, '(a4, 8x, a55)') &
      'Band', &
      'Number of excluded bands followed by indices.'
!
    DO i=1, SpBasic%n_band
      WRITE(iu_spc, '(i5, 7x, i5)') i, SpBasic%n_band_exclude(i)
      IF (SpBasic%n_band_exclude(i) > 0) THEN
        WRITE(iu_spc, '(14x, 8(3x, i5))') &
          SpBasic%index_exclude(1:SpBasic%n_band_exclude(i), i)
      ENDIF
    ENDDO
!
    WRITE(iu_spc, '(a4)') '*END'
!
!
!
  END SUBROUTINE write_block_14_int



  SUBROUTINE write_block_15_int(i_a, SpBasic, SpAerosol)

    IMPLICIT NONE


    INTEGER :: i_a
!     Species of aerosol
    TYPE  (StrSpecBasic), Intent(IN) :: SpBasic
    TYPE  (StrSpecAerosol), Intent(IN) :: SpAerosol

!   Local variables.
    INTEGER :: i
!     Loop variable
    INTEGER :: k
!     Loop variable


    IF (SpAerosol%i_aerosol_parm(i_a) == IP_aerosol_param_dry) THEN
      WRITE(iu_spc, '(a19, a16, a16)') &
        '*BLOCK: TYPE =   15', ': SUBTYPE =    0', ': VERSION =    0'
      WRITE(iu_spc, '(a42)') &
        'Parameters for dry aerosol optical depths.'
      WRITE(iu_spc, '(a14,i5)') &
        'n_aod_wavel = ',SpAerosol%n_aod_wavel
      WRITE(iu_spc, '(a19, i5, a15, i5)') &
        'Index of species = ' , i_a, ' Type of AOD = ', &
        SpAerosol%i_aod_type(i_a)
      WRITE(iu_spc, '(a4, 8x, a10, 10x, a10)') &
        'Band', 'Absorption', 'Scattering'
!     Write out the scattering parameters for each band
      DO i=1, Spectrum%Aerosol%n_aod_wavel
        WRITE(iu_spc, fmt='(i5, 2(4x, 1pe12.5))') &
          i, SpAerosol%aod_abs(1, i_a, i) , &
             SpAerosol%aod_scat(1, i_a, i)
      END DO
    ELSE IF (SpAerosol%i_aerosol_parm(i_a) == IP_aerosol_param_moist) THEN
      WRITE(iu_spc, '(a19, a16, a16)') &
        '*BLOCK: TYPE =   15', ': SUBTYPE =    1', ': VERSION =    0'
      WRITE(iu_spc, '(a44)') &
        'Parameters for moist aerosol optical depths.'
      WRITE(iu_spc, '(a14,i5)') &
        'n_aod_wavel = ',SpAerosol%n_aod_wavel
      WRITE(iu_spc, '(a19, i5, a15, i5)') &
        'Index of species = ' , i_a, ' Type of AOD = ', &
        SpAerosol%i_aod_type(i_a)
!     Write out the scattering parameters for each band
      DO i = 1 , Spectrum%Aerosol%n_aod_wavel
        WRITE(iu_spc, '(a7, i5, /, a10, 10x, a10)') 'Band = ', i, &
          'Absorption', 'Scattering'
        DO k=1, SpAerosol%nhumidity(i_a)
          WRITE(iu_spc, fmt='(2(4x, 1pe12.5))') &
            SpAerosol%aod_abs(k, i_a, i) , &
            SpAerosol%aod_scat(k, i_a, i)
        END DO
      END DO
    ENDIF
    WRITE(iu_spc, '(a4)') '*END'


  END SUBROUTINE write_block_15_int



  SUBROUTINE write_block_17(SpBasic, SpVar)

    IMPLICIT NONE

    TYPE (StrSpecBasic), INTENT(IN) :: SpBasic
    TYPE (StrSpecVar),   INTENT(IN) :: SpVar

    ! Write block 17 in standard spectral file
    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =   17' , ': SUBTYPE =    0' , ': VERSION =    0'
    WRITE(iu_spc, '(a)') &
      'Specification of sub-bands for spectral variability data.'
    IF (SpVar%n_sub_band > SpBasic%n_band) WRITE(iu_spc, '(a)') &
      'Wavelength limits (m) and Rayleigh coefficients at STP (m2/kg).'
    WRITE(iu_spc, '(a, i0)') &
      'Number of spectral sub-bands = ', SpVar%n_sub_band
    IF (SpVar%n_sub_band > SpBasic%n_band) THEN
      WRITE(iu_spc, '(a)') &
        'Sub-band  Band  k-term    Lower limit         Upper limit' // &
        '       Rayleigh coeff'
      DO i=1, SpVar%n_sub_band
        WRITE(iu_spc, '(3i7, 2x, 1pe16.9, 2(4x, 1pe16.9))') &
          i, SpVar%index_sub_band(1:2,i), SpVar%wavelength_sub_band(1:2,i), &
          SpVar%rayleigh_coeff(i, 0)
      END DO
    END IF
    WRITE(iu_spc, '(a4)') '*END'

    IF (SpVar%n_times > 0) THEN
      ! Write spectral variability file
      i=INDEX(file_spectral, ' ') - 1
      spectral_var = file_spectral(1:i) // '_var'
      CALL get_free_unit(ierr, iu_spc2)
      IF (ierr /= i_normal) RETURN
      OPEN(UNIT=iu_spc2, FILE=spectral_var, IOSTAT=ios, STATUS='UNKNOWN')
      WRITE(iu_spc2, '(a, i0)') &
        'Number of times in look-up table = ', SpVar%n_times
      IF (SpVar%n_repeat_times > 0) WRITE(iu_spc2, '(a, i0)') &
        'Number of times for periodic repetition = ', SpVar%n_repeat_times
      IF (SpVar%n_rayleigh_coeff > 0) WRITE(iu_spc2, '(a, i0)') &
        'Number of Rayleigh coefficients given = ', SpVar%n_rayleigh_coeff
      WRITE(iu_spc2, '(a)') &
        'Year  Month  Day(of month)  Seconds(since midnight)  TSI(Wm-2 at 1 AU)'
      WRITE(iu_spc2, '(a)') &
        'Fraction of solar flux in each sub-band.'
      IF (SpVar%n_rayleigh_coeff > 0) WRITE(iu_spc2, '(a,i0,a)') &
        'Rayleigh coefficient in the first ', SpVar%n_rayleigh_coeff, &
        ' sub-bands.'
      WRITE(iu_spc2, '(a)') '*BEGIN: spectral variability data'
      DO i=1, SpVar%n_times
        WRITE(iu_spc2, '(4(i6),4x,1pe16.9)') &
          SpVar%time(1:4, i), SpVar%total_solar_flux(i)
        WRITE(iu_spc2, '(5(1pe16.9))') &
          SpVar%solar_flux_sub_band(1:SpVar%n_sub_band, i)
        IF (SpVar%n_rayleigh_coeff > 0) WRITE(iu_spc2, '(5(1pe16.9))') &
          SpVar%rayleigh_coeff(1:SpVar%n_rayleigh_coeff, i)
      END DO
      CLOSE(iu_spc2)
    END IF

  END SUBROUTINE write_block_17



  SUBROUTINE write_block_18

!   Local variables
    INTEGER :: i_band
!     Loop variable

    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =   18', ': SUBTYPE =    0', ': VERSION =    0'
    WRITE(iu_spc, '(a)') 'Generalised continua in each interval'
    WRITE(iu_spc, '(a, /A)') &
      '(The number is the indexing number of the continuum as set out', &
      ' in the summary block 0.)'
    WRITE(iu_spc, '(a)') &
      'A zero indicates that there in no generalised continuum absorption ' // &
      'in the interval.'
    WRITE(iu_spc, '(a, 8x, a, /, 12x, a)') 'Band', &
      'Number of active continua, flag for continuum being major absorber,', &
      'followed by indexing numbers'

    DO i_band = 1, Spectrum%Basic%n_band
        WRITE(iu_spc, '(i5, 7x, i5, 7x, l5)') i_band, &
          Spectrum%ContGen%n_band_cont(i_band), &
          Spectrum%ContGen%l_cont_major(i_band)
      IF (Spectrum%ContGen%n_band_cont(i_band) > 0) THEN
        WRITE(iu_spc, '(5x, 4(2x, i3))') &
          Spectrum%ContGen%index_cont(1:Spectrum%ContGen%n_band_cont(i_band), &
          i_band)
      ENDIF
    ENDDO

    WRITE(iu_spc, '(a4)') '*END'

  END SUBROUTINE write_block_18



  SUBROUTINE write_block_19(SpBasic, SpCont)

    TYPE  (StrSpecBasic),   INTENT(IN) :: SpBasic
    TYPE  (StrSpecContGen), INTENT(IN) :: SpCont

!   Local variables.
    INTEGER :: i, j, k
!     Loop variables
    INTEGER :: i_index
!     Index of continuum

    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =   19', ': SUBTYPE =    0', ': VERSION =    0'
    WRITE(iu_spc, '(a)') &
      'Exponential sum fiting continuum coefficients: (exponents: m5/kg2)'
    WRITE(iu_spc, '(a, i5)') &
      'Number of temperatures in look-up table = ', SpCont%n_t_lookup_cont
    WRITE(iu_spc, '(a4, 8x, a, /, 12x, a)') 'Band', &
      'Continuum, number of k-terms, overlap treatment, ', &
      'followed by, k-terms and weights.'

    WRITE(iu_spc1,'(/,a)') '*BLOCK: continuum k-table'
    WRITE(iu_spc1,'(/,2(a,i4),a)') 'Lookup table: ', &
      SpCont%n_t_lookup_cont, ' temperatures.'
    WRITE(iu_spc1,'(6(1PE13.6))') &
      SpCont%t_lookup_cont(1:SpCont%n_t_lookup_cont)

    DO i=1, SpBasic%n_band
      DO j=1, SpCont%n_band_cont(i)
        i_index = SpCont%index_cont(j, i)
        WRITE(iu_spc, '(i5, 3(7x, i5))') &
          i, i_index, SpCont%i_band_k_cont(i, i_index), &
          SpCont%i_cont_overlap_band(i, i_index)
        DO k=1, SpCont%i_band_k_cont(i, i_index)
          IF (SpCont%i_scat_cont(k, i, i_index)==0) THEN
            WRITE(iu_spc, '(2(3x, 1pe16.9))') &
              SpCont%k_cont(k, i, i_index), SpCont%w_cont(k, i, i_index)
          ELSE
            WRITE(iu_spc, '(2(3x, 1pe16.9),i3)') &
              SpCont%k_cont(k, i, i_index), SpCont%w_cont(k, i, i_index), &
              SpCont%i_scat_cont(k, i, i_index)
          END IF
        ENDDO
        WRITE(iu_spc1,'(/,3(a,i4))') 'Band: ',i,', continuum: ',i_index, &
          ', k-terms: ',SpCont%i_band_k_cont(i, i_index)
        DO k=1, SpCont%i_band_k_cont(i, i_index)
          WRITE(iu_spc1,'(6(1PE13.6))') &
            SpCont%k_lookup_cont(:,k,i_index,i)
        END DO
      ENDDO
    ENDDO

    WRITE(iu_spc, '(a4)') '*END'
    WRITE(iu_spc1, '(a4)') '*END'

  END SUBROUTINE write_block_19



  SUBROUTINE write_block_20

    IMPLICIT NONE

    INTEGER :: i, i_wl, n_t
    
    WRITE(iu_spc, '(a19, a16, a16)') &
      '*BLOCK: TYPE =   20', ': SUBTYPE =    0', ': VERSION =    0'
    WRITE(iu_spc, '(a)') &
      'Photolysis quantum yield lookup tables'
    WRITE(iu_spc, '(a,i5)') &
      'Max number of temperatures:', Spectrum%Dim%nd_t_lookup_photol
    WRITE(iu_spc, '(a,i6)') &
      'Max number of wavelengths:', Spectrum%Dim%nd_wl_lookup_photol
    DO i=1, Spectrum%Photol%n_pathway
      WRITE(iu_spc, '(/,a,i4,a,1pe16.9)') 'Pathway index:', &
        i, ', Threshold wavelength:', Spectrum%Photol%threshold_wavelength(i)
      n_t = Spectrum%Photol%n_t_lookup_photol(i)
      WRITE(iu_spc, '(a,i5,a,3(3x,1pe16.9))') 'Temperatures:', &
        n_t, ':', Spectrum%Photol%t_lookup_photol(1:MIN(n_t, 3), i)
      IF (n_t > 3) WRITE(iu_spc, '((19x, 3(3x,1pe16.9)))') &
        Spectrum%Photol%t_lookup_photol(4:n_t, i)
      WRITE(iu_spc, '(a,i6,a)') 'Wavelengths:', &
        Spectrum%Photol%n_wl_lookup_photol(i), ',     Quantum Yield:'
      DO i_wl=1, Spectrum%Photol%n_wl_lookup_photol(i)
        WRITE(iu_spc, '(4(3x,1pe16.9))') &
          Spectrum%Photol%wl_lookup_photol(i_wl, i), &
          Spectrum%Photol%quantum_yield(1:MIN(n_t, 3), i_wl, i)
        IF (n_t > 3) WRITE(iu_spc, '((19x, 3(3x,1pe16.9)))') &
          Spectrum%Photol%quantum_yield(4:n_t, i_wl, i)
      END DO
    END DO
    WRITE(iu_spc, '(a4)') '*END'

  END SUBROUTINE write_block_20

END SUBROUTINE out_spectrum
