! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 11.
!
! Method:
!   Scattering data are given in a parametrized form where the 
!   extinctions are proportion to the mass loading. For this,
!   a file of scattering properties is supplied. In the case of
!   dry aerosols, only one set of properties for each band is 
!   required. In the case of moist aerosols this file must
!   contain a number of blocks, each giving the single scattering
!   properties in each band at a particular value of the humidity.
!   A look-up table is constructed from these values.
!
!- ---------------------------------------------------------------------
SUBROUTINE make_block_11(Spectrum, ierr)

  USE realtype_rd
  USE def_spectrum
  USE rad_pcf
  USE scatter_pp_pcf
  USE file_type_pcf
  USE def_std_io_icf

  IMPLICIT NONE

  TYPE (StrSpecData), Intent(INOUT), TARGET :: Spectrum
!   Spectral file to be assigned
  INTEGER, Intent(INOUT) :: ierr
!   Error flag

! Local variables.
  CHARACTER (LEN=1) :: char_yn
!   Character response variable
  INTEGER :: i_scatter_type
!   Type of scatterer
  INTEGER :: iu_average
!   Unit number for reading averaged data
  INTEGER :: i_block
!   Block index
  INTEGER :: i_component
!   Aerosol component
  INTEGER :: i_species
!   Species of aerosol
  INTEGER :: idum
!   Dummy variable
  INTEGER :: i, k
!   Loop variables
  INTEGER :: i_input_type
!   Input type
  INTEGER :: n_phf_term
!   Temporary number of terms in the phase function
  INTEGER :: nhumidity
!   Temporary number of humidity values
  REAL (RealK) :: absorption(Spectrum%Dim%nd_band)
!   Absorption values of block
  REAL (RealK) :: scattering(Spectrum%Dim%nd_band)
!   Scattering values of block
  REAL (RealK), ALLOCATABLE :: phf_fnc(:, :, :, :)
  REAL (RealK), ALLOCATABLE :: phase_fnc(:, :)
!   Asymmetry values of block
  REAL (RealK) :: radius_eff, radius_eff_dry
!   Effective radius
  REAL (RealK) :: number_total
!   Total number density
  REAL (RealK) :: volume_fraction, volume_fraction_dry
!   Volume fraction
  REAL (RealK) :: humidity
!   Humidity

! Include data for aerosols
  INCLUDE 'aerosol_component.finc'


  IF (.NOT.Spectrum%Basic%l_present(11)) THEN
    IF (ALLOCATED(Spectrum%Aerosol%l_aero_spec)) &
       DEALLOCATE(Spectrum%Aerosol%l_aero_spec)
    ALLOCATE(Spectrum%Aerosol%l_aero_spec( &
                              Spectrum%Dim%nd_aerosol_species))
    Spectrum%Aerosol%l_aero_spec = .FALSE. 
    IF (ALLOCATED(Spectrum%Aerosol%i_aerosol_parm)) &
       DEALLOCATE(Spectrum%Aerosol%i_aerosol_parm)
    ALLOCATE(Spectrum%Aerosol%i_aerosol_parm( &
                              Spectrum%Dim%nd_aerosol_species))
!     Flag for humidity on-off
    IF (ALLOCATED(Spectrum%Aerosol%nhumidity)) &
       DEALLOCATE(Spectrum%Aerosol%nhumidity)
    ALLOCATE(Spectrum%Aerosol%nhumidity( &
                              Spectrum%Dim%nd_aerosol_species))
!     This is the number of humidity values for each 
!     Aerosol species. note that a value of zero indicates 
!     Dry aerosol (no humidity variation).
    IF (ALLOCATED(Spectrum%Aerosol%humidities)) &
       DEALLOCATE(Spectrum%Aerosol%humidities)
    ALLOCATE(Spectrum%Aerosol%humidities( &
                              Spectrum%Dim%nd_humidity, &
                              Spectrum%Dim%nd_aerosol_species))
!     These are the humidity values for each aerosol species.
!     They are assumed to occur in increasing order, and an 
!     error message should be printed if an attempt is made 
!     to break this rule.
    IF (ALLOCATED(Spectrum%Aerosol%n_aerosol_phf_term)) &
       DEALLOCATE(Spectrum%Aerosol%n_aerosol_phf_term)
    ALLOCATE(Spectrum%Aerosol%n_aerosol_phf_term( &
                              Spectrum%Dim%nd_aerosol_species))
    IF (ALLOCATED(Spectrum%Aerosol%abs)) &
       DEALLOCATE(Spectrum%Aerosol%abs)
    ALLOCATE(Spectrum%Aerosol%abs( &
                              Spectrum%Dim%nd_humidity, &
                              Spectrum%Dim%nd_aerosol_species, &
                              Spectrum%Dim%nd_band))
    IF (ALLOCATED(Spectrum%Aerosol%scat)) &
       DEALLOCATE(Spectrum%Aerosol%scat)
    ALLOCATE(Spectrum%Aerosol%scat( &
                              Spectrum%Dim%nd_humidity, &
                              Spectrum%Dim%nd_aerosol_species, &
                              Spectrum%Dim%nd_band))
    IF (ALLOCATED(Spectrum%Aerosol%phf_fnc)) &
       DEALLOCATE(Spectrum%Aerosol%phf_fnc)
    ALLOCATE(Spectrum%Aerosol%phf_fnc( &
                              Spectrum%Dim%nd_humidity, &
                              Spectrum%Dim%nd_phase_term, &
                              Spectrum%Dim%nd_aerosol_species, &
                              Spectrum%Dim%nd_band))
  END IF


! Obtain the input data.
  CALL get_free_unit(ierr, iu_average)

21 CALL open_file_in(ierr, iu_average, &
    'Enter the name of the file of averaged ' &
    //'scattering properties.')
  IF (ierr /= i_normal) THEN
    ierr=i_err_fatal
    RETURN
  ENDIF


! Read in and process each block.

! Read the input type. This is used to control the subsequent reading in.
2 READ(iu_average, '(13x,i5)', end=100) i_input_type
  READ(iu_average, '(///, 11x, i5, /, 23x, i5)', END=100) &
    i_block, i_scatter_type
  IF (i_scatter_type /= IP_type_aerosol) THEN
    WRITE(iu_err, '(a)') '+++ The scatterer is not an aerosol.'
    CLOSE(iu_average)
    goto 21
  ENDIF
  READ(iu_average, '(24x, i5)' ) i_component

! Find the species of aerosol which represents this component.
  i=0
20 i=i+1
  IF (Spectrum%Aerosol%type_aerosol(i) == i_component) THEN
    i_species=i
  ELSE IF (i < Spectrum%Aerosol%n_aerosol) THEN
    goto 20
  ELSE
    WRITE(iu_err, '(/a)') '+++ Illegal component'
    WRITE(iu_err, '(a)') 'This file contains data for an ' &
      //'aerosol which is not in the spectral file.'
    WRITE(iu_stdout, '(a)') 'Repeat this step.'
    CLOSE(iu_average)
    goto 21
  ENDIF

  IF ( (i_input_type  ==  IT_file_ave_mie_humid).OR. &
       (i_input_type  ==  IT_file_ave_phf_mie_humid) ) THEN
    READ(iu_average, '(/, 14x,f7.5)') humidity
    READ(iu_average, '(/, 22x, 1pe12.5, /, 22x, 1pe12.5, /, 39x, '// &
      '1pe12.5, /, 22x, 1pe12.5, /, 39x, 1pe12.5, /)' ) &
      number_total, radius_eff, radius_eff_dry, &
      volume_fraction, volume_fraction_dry
  ELSEIF ( (i_input_type == IT_file_ave_mie_dry).OR. &
           (i_input_type == IT_file_ave_phf_mie_dry) ) THEN
    humidity = 0.0_RealK
    READ(iu_average, '(/, 3(22x, 1pe12.5, /))' ) &
      number_total, radius_eff, volume_fraction
  ENDIF
  IF ( (i_input_type == IT_file_ave_phf_mie_humid).OR. &
       (i_input_type == IT_file_ave_phf_mie_dry) ) THEN
    READ(iu_average, '(43x, i3, //)') n_phf_term
    IF (n_phf_term > Spectrum%Dim%nd_phase_term) THEN
      IF (ALLOCATED(Spectrum%Aerosol%phf_fnc)) THEN
        ALLOCATE(phf_fnc(Spectrum%Dim%nd_humidity, &
                         Spectrum%Dim%nd_phase_term, &
                         Spectrum%Dim%nd_aerosol_species, &
                         Spectrum%Dim%nd_band))
        phf_fnc = Spectrum%Aerosol%phf_fnc
        DEALLOCATE(Spectrum%Aerosol%phf_fnc)
        ALLOCATE(Spectrum%Aerosol%phf_fnc( &
                         Spectrum%Dim%nd_humidity, &
                         n_phf_term, &
                         Spectrum%Dim%nd_aerosol_species, &
                         Spectrum%Dim%nd_band))
        Spectrum%Aerosol%phf_fnc(:,1:Spectrum%Dim%nd_phase_term,:,:) = phf_fnc
        DEALLOCATE(phf_fnc)
      END IF
      Spectrum%Dim%nd_phase_term = n_phf_term
    END IF
  ELSE
      n_phf_term=1
  ENDIF
  ALLOCATE(phase_fnc(n_phf_term, Spectrum%Dim%nd_band))
  READ(iu_average, '(/)')

  DO i=1, Spectrum%Basic%n_band
    READ(iu_average, '(i5,3(4x, 1pe16.9), /, (5x, 3(4x, 1pe16.9)))') &
      idum, absorption(i), scattering(i), &
      (phase_fnc(k, i), k=1, n_phf_term)
  ENDDO
  READ(iu_average, '(//)')

  IF (i_scatter_type /= IP_type_aerosol) THEN
    WRITE(iu_err, '(/a, i5, a)') '+++ Warning: This type ' &
      //'of scatterer is not aerosol: block ', &
      i_block, ' will be skipped.'
    goto 2
  ENDIF

! Check whether this component is already present.
! If it is present, the action to be taken depends on whether 
! the aerosol is dry or moist.
  IF (Spectrum%Aerosol%l_aero_spec(i_species)) THEN
    IF (Spectrum%Aerosol%nhumidity(i_species)  ==  0) THEN
!     We are trying to overwrite a dry aerosol. This should
!     never be done automatically, so stop the program or 
!     ask the user what to do.
      WRITE(iu_stdout, '(a, /a, i5, 1x, a1, a20, a1, 1x, a15)') &
        '+++ Warning: ', 'Data for (dry) aerosol component ', &
        i_component, '(', name_aerosol_component(i_component), ')', &
        'already exist.'
      WRITE(iu_stdout, '(a)') 'Do you wish to overwrite? (y/n)'
3     READ(iu_stdin, '(a)') char_yn
      IF ( (char_yn /= 'Y').AND.(char_yn /= 'y').and. &
           (char_yn /= 'N').AND.(char_yn /= 'n') ) THEN
        WRITE(iu_stdout, '(a)') '+++ Illegal response: ' &
          //'Please re-enter.'
        goto 3
      ENDIF
    ELSE
!     We are trying to overwrite a moist aerosol. In this
!     case, we might want to add data for a new, larger 
!     value of humidity. Check that the new value really is 
!     larger than the previous one.
      IF (humidity > Spectrum%Aerosol%humidities( &
        Spectrum%Aerosol%nhumidity(i_species), i_species) ) THEN
             char_yn = 'y'
      ELSE
        WRITE(iu_stdout, '(/a, i5, 1x, a1, a20, a1, /, a28, 1x, a48)') &
          '+++ Warning: Data for (moist) aerosol component ', &
          i_component, '(', name_aerosol_component(i_component), &
          ')', 'are not in ascending order of ', &
          'humidity, or moist aerosol is being overwritten.'
        WRITE(iu_stdout, '(a)') 'Do you wish to overwrite? (y/n)'
        DO
          READ(iu_stdin, '(a)') char_yn
          IF ( (char_yn == 'N').OR.(char_yn == 'n') ) THEN
            WRITE(iu_err, '(/a)') 'Program terminating.'
            ierr=i_err_fatal
            RETURN
          ELSE IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
            Spectrum%Aerosol%nhumidity(i_species) = 0
            EXIT
          ELSE
            WRITE(iu_stdout, '(a)') '+++ Illegal response: ' &
              //'Please re-enter.'
          END IF
        END DO
      ENDIF
    ENDIF
  ENDIF

! If this type of aerosol is new, initialise the number of
! humidities to be zero. this will be correct if it is a dry
! aerosol. If it is a moist aerosol, it will be incremented to
! 1 further down, before assigning to the main arrays.
  IF ( .NOT. Spectrum%Aerosol%l_aero_spec(i_species) ) THEN
    Spectrum%Aerosol%nhumidity(i_species) = 0
  ENDIF

  IF ( (.NOT.Spectrum%Aerosol%l_aero_spec(i_species)).OR. &
       (char_yn == 'Y').OR.(char_yn == 'y') ) THEN

!   Assign to the main arrays and scale out the volume fraction.
    i_component=Spectrum%Aerosol%type_aerosol(i_species)

    IF ( (i_input_type == IT_file_ave_mie_dry).OR. &
         (i_input_type == IT_file_ave_phf_mie_dry) ) THEN

!     Dry aerosol.
      Spectrum%Aerosol%i_aerosol_parm(i_species)=IP_aerosol_param_dry
      Spectrum%Aerosol%n_aerosol_phf_term(i_species)=n_phf_term
      DO i=1, Spectrum%Basic%n_band
        Spectrum%Aerosol%abs(1,i_species, i)=MAX(absorption(i) &
          /(volume_fraction*density_component(i_component)), 0.0_RealK)
        Spectrum%Aerosol%scat(1,i_species, i)=MAX(scattering(i) &
          /(volume_fraction*density_component(i_component)), 0.0_RealK)
        DO k=1, n_phf_term
          Spectrum%Aerosol%phf_fnc(1, k, i_species, i)=phase_fnc(k, i)
        ENDDO
      ENDDO

    ELSE IF ( (i_input_type == IT_file_ave_mie_humid).OR. &
              (i_input_type == IT_file_ave_phf_mie_humid) ) THEN

!     Moist aerosol.
      Spectrum%Aerosol%i_aerosol_parm(i_species)=IP_aerosol_param_moist
      nhumidity = Spectrum%Aerosol%nhumidity(i_species)+1
      Spectrum%Aerosol%nhumidity(i_species) = nhumidity
      Spectrum%Aerosol%n_aerosol_phf_term(i_species) = n_phf_term
      DO i=1, Spectrum%Basic%n_band
        Spectrum%Aerosol%abs(nhumidity, i_species, i)=MAX(absorption(i) &
          /(volume_fraction_dry*density_component(i_component)), 0.0_RealK)
        Spectrum%Aerosol%scat(nhumidity, i_species, i)=MAX(scattering(i) &
          /(volume_fraction_dry*density_component(i_component)), 0.0_RealK)
        DO k=1, n_phf_term
          Spectrum%Aerosol%phf_fnc(nhumidity, k, i_species, i) &
            = phase_fnc(k, i)
        ENDDO
      ENDDO
      Spectrum%Aerosol%humidities(nhumidity, i_species)=humidity
    ELSE
      WRITE(iu_err, '(/a, i5, 1x, a1, a20, a1, /, a26,1x,a39)') &
        '+++ Warning: Data for aerosol component ', i_component, &
        '(', name_aerosol_component(i_component), ')', &
        'have an illegal input type.', 'program terminating.'
      ierr=i_err_fatal
      RETURN
    ENDIF
!   This species now has data.
    Spectrum%Aerosol%l_aero_spec(i_species)=.TRUE.
  ENDIF
  DEALLOCATE(phase_fnc)
  goto 2


100 Spectrum%Basic%l_present(11)=.FALSE.
  DO i=1, Spectrum%Aerosol%n_aerosol
    Spectrum%Basic%l_present(11)=Spectrum%Basic%l_present(11).OR. &
      Spectrum%Aerosol%l_aero_spec(i)
  ENDDO

! Check for assignment of data.
  IF (.NOT.Spectrum%Basic%l_present(11)) THEN
    WRITE(iu_err, '(/a)') '+++ Warning from make_block_11'
    WRITE(iu_err, '(a)') 'No aerosol properties have been assigned.'
  ENDIF

  CLOSE(iu_average)

END SUBROUTINE make_block_11
