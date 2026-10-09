! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to fit parametrziations to cloud data.
!
SUBROUTINE cloud_fit_90 &
!
(l_interactive, fit_species, n_band, n_data, vol_frac, d, &
 particle_density, &
 absorption_ave, scattering_ave, n_phf_term, phf_fnc_ave, &
 ierr)
!
! Description:
!   The desired parametrization is selected. For each band, and for
! each optical property, the initial data are conditioned to give 
! a well-behaved fitting problem: this is solved using a conjugate
! gradient algorithm and the results are written to an output file.
! For clarity most work is carried out in internal subroutines.
!
!
!
! Modules used
  USE realtype_rd
  USE rad_pcf
  USE def_std_io_icf
  USE file_type_pcf
  USE cloud_fitting

  IMPLICIT NONE
!
!
!
! Dummy arguments.
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive operation
  CHARACTER  (LEN=5), Intent(IN) :: fit_species
!   Phase of condensate to be fitted
  INTEGER, Intent(IN) :: n_band
!   Number of spectral bands
  INTEGER, Intent(IN) :: n_data
!   Size of input data
!
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
!
  REAL  (RealK), Intent(IN), Dimension(:) :: vol_frac
!   Volume fractions for each block
  REAL  (RealK), Intent(IN), Dimension(:) :: d
!   Characteristic dimension of particles
  REAL  (RealK), Intent(IN), Dimension(:) :: particle_density
!   Internal density of particles
!
  REAL  (RealK), Intent(IN), Dimension(:, :) :: absorption_ave
!   Mean absorption for each block in each band
  REAL  (RealK), Intent(IN), Dimension(:, :) :: scattering_ave
!   Mean scattering for each block in each band
  INTEGER, Intent(IN) :: n_phf_term
!   Number of terms in the phase function
  REAL  (RealK), Intent(IN), Dimension(:, :, :) :: phf_fnc_ave
!   Mean phase function for each block in each band
!
!
!
!- End of header
!
! Local variables.
!
  INTEGER :: ios
!   I/O error flag
!
  INTEGER :: iu_cloud_fit
!   Unit number for input from a file
  INTEGER :: iu_monitor
!   Unit number for monitoring information
!
  INTEGER :: i_fit
!   Fitting scheme selected
  REAL  (RealK) :: density_particle(n_data)
!   Density of the scattering particle
!
  INTEGER :: n_property
!   Number of optical properties to be fitted using this scheme,
!   counint moments of the phase function as separate.
  CHARACTER  (LEN=10), Allocatable, Dimension(:) :: property
!   Name of the optical propety to be fitted: this varies from
!   scheme to scheme
!
  REAL  (RealK), Allocatable, Dimension(:) :: ds
!   Scaled sizes of particles
  REAL  (RealK) :: d_break
!   Break in dimension for ADT fit
  REAL  (RealK), Allocatable, Dimension(:) :: actual
!   The actual data to be fitted for each property
  REAL  (RealK), Allocatable, Dimension(:) :: scaled_fit
!   The fit to the scaled data
  REAL  (RealK), Allocatable, Dimension(:) :: parm_list
!   The list of fitting parameters
!
  REAL  (RealK) :: scaling_d
!   Scaling for the dimension of the particles
  REAL  (RealK) :: scaling
!   The scaling applied to the actual optical property to yield
!   the numerically conditioned equivalent
  INTEGER :: first_place
!   First position in the array for the current section of parameters
  INTEGER :: last_place
!   Last position in the array for the current section of parameters
  INTEGER :: i_term
!   Current term in the phase function
  INTEGER :: n_param
!   Total number of parameters in the scheme
!
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: i_b
!   Loop variable

! Functions called:
  INTEGER, EXTERNAL :: set_n_cloud_parameter


! Open the output file.
  CALL get_free_unit(ierr, iu_cloud_fit)
  CALL open_file_out(ierr, iu_cloud_fit, &
    "Enter the name of the file to contain the fitting parameters.")
  IF (ierr /= i_normal) RETURN
!
! Open the monitoring file.
  CALL get_free_unit(ierr, iu_monitor)
  CALL open_file_out(ierr, iu_monitor, &
    "Enter the name of the file to contain the monitoring information.")
  IF (ierr /= i_normal) RETURN
!
  WRITE(iu_stdout, '(/a)') &
    'Enter density of particle'
  DO 
    READ(iu_stdin, *, IOSTAT=ios) density_particle(n_data)
    IF (ios == 0) THEN
      EXIT
    ELSE
      WRITE(iu_stdout, '(/a)') &
        '*** Error: Invalid density '
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(/a)') &
          'Please re-enter '
      ELSE
       ierr=i_err_fatal
       RETURN
      ENDIF
    ENDIF
  ENDDO
  DO i = 1, n_data
    ! Use actual particle density where it is set
    IF (particle_density(i) > TINY(0.0_RealK)) THEN
      density_particle(i) = particle_density(i)
    ELSE
      density_particle(i) = density_particle(n_data)
    END IF
  END DO
!
! Choose the parametrization scheme
  CALL choose_scheme
!
! Construct a list of properties to be fitted.
  CALL set_properties_fit
!
  SELECT CASE(fit_species)
    CASE("Water")
      n_param = set_n_cloud_parameter(i_fit, ip_clcmp_st_water, n_phf_term)
    CASE("Ice  ")
      n_param = set_n_cloud_parameter(i_fit, ip_clcmp_st_ice, n_phf_term)
  END SELECT
  ALLOCATE(parm_list(n_param))
!
! Write headers to the output file.
  WRITE(iu_cloud_fit, '(a13, i5, /20x, a)') '*FILE_TYPE = ' &
    , IT_file_cloud_fit_phf, 'Parameters fitting single scattering'
  WRITE(iu_cloud_fit, '(/a38, 2(1pe12.5))') &
     'Min and max characteristic dimension: ', &
     MINVAL(d(1:n_data)), MAXVAL(d(1:n_data))
  WRITE(iu_cloud_fit, '(a33, i5)') &
     'Index of parametrization scheme =', i_fit
  WRITE(iu_cloud_fit, '(a31, i5, 1x, a14)') &
     'Number of fitting parameters = ', n_param, &
     'for each band.'
  WRITE(iu_cloud_fit, '(a31, i3)') &
     'Number of terms (Phase fnc.) = ', n_phf_term
!
  ALLOCATE(ds(n_data))
  ALLOCATE(actual(n_data))
!
!
  DO i_b=1, n_band
!
    DO i = 1, n_property
!
      WRITE(iu_monitor, '(a)') &
        '---------------------------------------------------------'
      WRITE(iu_monitor, '(/, a16, i3, a2, a10, a10, i3, //)') &
        'Fitted property ', i, ' (', property(i), ') in band ', i_b
!
!     Set the position in the full array of parameters.
      CALL set_position_param
!
!     Initialize for the fit.
      parm_list(first_place:last_place) = 0.0_RealK
!
!     For well-behaved fitting, the data must be conditioned so that O(1)
!     changes in the scaled size result in O(1) changes in the fitted
!     function.
      CALL condition_data
!
      CALL conjugate_gradient_cloud_90(ierr, iu_monitor, &
        fit_species, &
        n_data, ds, vol_frac, actual, &
        i_fit, property(i), last_place-first_place+1, &
        parm_list(first_place:last_place))
!
      ALLOCATE(scaled_fit(SIZE(ds)))
      WRITE(iu_monitor, '(/, a, /)') &
        'Scaled Actual data and Fit:'
      WRITE(iu_monitor, '(/, 3(3x, a12), /)') &
        'Scaled Size', 'Actual', 'Fitted'
      CALL cloud_fit_eval(fit_species, ds, i_fit, &
        property(i), parm_list(first_place:last_place), scaled_fit)
      DO j=1, SIZE(ds)
        WRITE(iu_monitor, '(3(3x, 1pe12.5))') &
          ds(j), actual(j), scaled_fit(j)
      ENDDO
      WRITE(iu_monitor, '(//)')
      DEALLOCATE(scaled_fit)
!
!     The parameters thus far yield the conditioned data and must now
!     be scaled to give the real fits.
      CALL rescale_parm(parm_list(first_place:last_place))
!
    ENDDO
!
!   Write the fit for this band to the output file.
    WRITE(iu_cloud_fit, '(a6, 1x, i5, 5x, a19)') &
      'Band =', i_b, 'Fitting Parameters:'
    WRITE(iu_cloud_fit, '(4(4x, 1pe12.5))') parm_list(1:n_param)
!
  ENDDO
!
  DEALLOCATE(actual)
  DEALLOCATE(ds)
  DEALLOCATE(parm_list)
!
  CLOSE(iu_cloud_fit)
  CLOSE(iu_monitor)
!
!
!
!
CONTAINS
!
!
!
  SUBROUTINE choose_scheme
!
!
!
    INTEGER :: ios
!     I/O error flag
!
!
!
    WRITE(iu_stdout, '(/a)') &
      'Enter the parametrization scheme to be used.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) i_fit
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ Erroneous response.'
        IF(l_interactive) THEN
          WRITE(iu_stdout, '(a)') 'Please re-enter.'
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ELSE IF ( ( (fit_species == "Water")             .AND. &
                  ( (i_fit /= IP_Slingo_Schrecker)     .AND. &
                    (i_fit /= IP_ackerman_stephens)    .AND. &
                    (i_fit /= IP_drop_pade_2)          .AND. &
                    (i_fit /= IP_Slingo_Schr_PHF)      .AND. &
                    (i_fit /= IP_drop_pade_2_PHF) ) )   .OR. &
                ( (fit_species == "Ice  ")             .AND. &
                  ( (i_fit /= IP_Slingo_Schrecker_ice) .AND. &
                    (i_fit /= IP_ice_adt)              .AND. &
                    (i_fit /= IP_ice_adt_10)           .AND. &
                    (i_fit /= IP_Slingo_Schr_ice_PHF)  .AND. &
                    (i_fit /= IP_ice_Fu_PHF)           .AND. &
                    (i_fit /= IP_ice_pade_2_phf)  ) ) ) THEN
        WRITE(iu_err, '(a)') '+++ Invalid parametrization scheme.'
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
!
!
  END SUBROUTINE choose_scheme
!
!
!
  SUBROUTINE set_properties_fit
!
!
!
    SELECT CASE (fit_species)
!
      CASE("Water")
        SELECT CASE (i_fit)
          CASE (IP_Slingo_Schrecker, IP_Ackerman_Stephens, &
            IP_Drop_Pade_2)
            n_property = 3
            ALLOCATE(property(n_property))
            property = (/ "Extinction", &
                          "Coalbedo  ", &
                          "Asymmetry " /)
          CASE (IP_Slingo_Schr_PHF, IP_drop_Pade_2_PHF)
            n_property = 2 + n_phf_term
            ALLOCATE(property(n_property))
            property(1:2) = (/ "Extinction", &
                               "Coalbedo  " /)
            property(3:2+n_phf_term) = &
              SPREAD("Moment   ", 1, n_phf_term)
        END SELECT
      CASE("Ice")
        SELECT CASE (i_fit)
          CASE (IP_Slingo_Schrecker_ice, IP_ice_adt, IP_ice_adt_10)
            n_property = 3
            ALLOCATE(property(n_property))
            property = (/ "Extinction", &
                          "Coalbedo  ", &
                          "Asymmetry " /)
          CASE (IP_Slingo_Schr_ice_PHF, IP_ice_Fu_PHF, IP_Ice_Pade_2_PHF)
            n_property = 2 + n_phf_term
            ALLOCATE(property(n_property))
            property(1:2) = (/ "Extinction", &
                               "Coalbedo  " /)
            property(3:2+n_phf_term) = &
              SPREAD("Moment   ", 1, n_phf_term)
      END SELECT
!
    END SELECT
!
!
!
  END SUBROUTINE set_properties_fit
!
!
!
  SUBROUTINE set_position_param
!
!
!
!   Set the position in the full array of parameters for this
!   property.
!
!
!
    SELECT CASE (fit_species)
      CASE("Water")
        SELECT CASE (i_fit)
          CASE (IP_Slingo_Schrecker)
            first_place = 2*i-1
            last_place  = 2*i
          CASE (IP_Ackerman_Stephens)
            first_place = 3*i-2
            last_place  = 3*i
          CASE (IP_Drop_Pade_2)
            SELECT CASE (i)
              CASE(1)
                first_place = 1
                last_place  = 6
              CASE(2)
                first_place = 7
                last_place  = 11
              CASE(3)
                first_place = 12
                last_place  = 16
            END SELECT
          CASE (IP_Slingo_Schr_PHF)
            first_place = 2*i-1
            last_place  = 2*i
            i_term      = i-2
          CASE (IP_Drop_Pade_2_PHF)
            SELECT CASE (i)
              CASE(1)
                first_place = 1
                last_place  = 6
              CASE(2)
                first_place = 7
                last_place  = 11
              CASE(3:)
                first_place = 5*i-3
                last_place  = 5*i+1
                i_term      = i-2
            END SELECT
        END SELECT
      CASE("Ice  ")
        SELECT CASE (i_fit)
          CASE (IP_Slingo_Schrecker_ice)
            first_place = 2*i-1
            last_place  = 2*i
          CASE (IP_ice_adt)
            first_place = 10*i-9
            last_place  = 10*i
          CASE (IP_ice_adt_10)
            first_place = 12*i-11
            last_place  = 12*i
          CASE (IP_Slingo_Schr_ice_PHF)
            first_place = 2*i-1
            last_place  = 2*i
            i_term       = i-2
          CASE (IP_ice_Fu_PHF)
            SELECT CASE (i)
              CASE(1)
                first_place = 1
                last_place  = 4
              CASE(2)
                first_place = 5
                last_place  = 9
              CASE(3:)
                first_place = 5*i-5
                last_place  = 5*i-1
                i_term      = i-2
            END SELECT
          CASE (IP_Ice_Pade_2_PHF)
            SELECT CASE (i)
              CASE(1)
                first_place = 1
                last_place  = 6
              CASE(2)
                first_place = 7
                last_place  = 11
              CASE(3:)
                first_place = 5*i-3
                last_place  = 5*i+1
                i_term      = i-2
            END SELECT
        END SELECT
    END SELECT
!
!
!
  END SUBROUTINE set_position_param
!
!
!
  SUBROUTINE condition_data
!
!
!
    SELECT CASE (fit_species)
!
      CASE("Water")
!
        SELECT CASE (i_fit)
!
          CASE (IP_Slingo_Schrecker)
!  
!           Scale the dimension:
            scaling_d  = REAL(n_data, RealK) / SUM(d(1:n_data))
            ds = scaling_d * d(1:n_data)
!           Repeat for the optical propeties
            SELECT CASE(property(i))
              CASE("Extinction")
                actual = (absorption_ave(i_b, 1:n_data) + &
                          scattering_ave(i_b, 1:n_data) ) / &
                         (vol_frac(1:n_data) * density_particle(1:n_data))
                scaling = REAL(n_data, RealK) / SUM(actual)
              CASE("Coalbedo  ")
                actual = absorption_ave(i_b, 1:n_data) / &
                  ( absorption_ave(i_b, 1:n_data) + &
                    scattering_ave(i_b, 1:n_data))
                scaling = MAXVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
              CASE("Asymmetry ")
                actual = phf_fnc_ave(i_b, 1, 1:n_data)
                scaling = 1.0_RealK
            END SELECT
            actual = scaling * actual
!
          CASE (IP_Ackerman_Stephens)
!  
!           Scale the dimension:
            scaling_d  = REAL(n_data, RealK) / SUM(d(1:n_data))
            ds = scaling_d * d(1:n_data)
!           Repeat for the optical propeties
            SELECT CASE(property(i))
              CASE("Extinction")
                actual = (absorption_ave(i_b, 1:n_data) + &
                          scattering_ave(i_b, 1:n_data) ) / &
                         (vol_frac(1:n_data) * density_particle(1:n_data))
                scaling = REAL(n_data, RealK) / SUM(actual)
              CASE("Coalbedo  ")
                actual = absorption_ave(i_b, :) / &
                  ( absorption_ave(i_b, 1:n_data) + &
                    scattering_ave(i_b, 1:n_data))
                scaling = MAXVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
              CASE("Asymmetry ")
                actual = phf_fnc_ave(i_b, 1, 1:n_data)
                scaling = 1.0_RealK
            END SELECT
            actual = scaling * actual
!
          CASE (IP_Drop_Pade_2)
!  
!           Scale the dimension:
            scaling_d  = REAL(n_data, RealK) / SUM(d(1:n_data))
            ds = scaling_d * d(1:n_data)
!           Repeat for the optical propeties
            SELECT CASE(property(i))
              CASE("Extinction")
                actual = (absorption_ave(i_b, 1:n_data) + &
                          scattering_ave(i_b, 1:n_data) ) / &
                         (vol_frac(1:n_data) * density_particle(1:n_data))
                scaling = REAL(n_data, RealK) / SUM(actual)
              CASE("Coalbedo  ")
                actual = absorption_ave(i_b, 1:n_data) / &
                  ( absorption_ave(i_b, 1:n_data) + &
                    scattering_ave(i_b, 1:n_data))
                scaling = MAXVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
              CASE("Asymmetry ")
                actual = phf_fnc_ave(i_b, 1, 1:n_data)
                scaling = 1.0_RealK
            END SELECT
            actual = scaling * actual
!
          CASE (IP_Slingo_Schr_PHF)
!  
!           Scale the dimension:
            scaling_d  = REAL(n_data, RealK) / SUM(d(1:n_data))
            ds = scaling_d * d(1:n_data)
!           Repeat for the optical propeties
            SELECT CASE(property(i))
              CASE("Extinction")
                actual = (absorption_ave(i_b, 1:n_data) + &
                          scattering_ave(i_b, 1:n_data) ) / &
                         (vol_frac(1:n_data) * density_particle(1:n_data))
                scaling = REAL(n_data, RealK) / SUM(actual)
              CASE("Coalbedo  ")
                actual = absorption_ave(i_b, 1:n_data) / &
                  ( absorption_ave(i_b, 1:n_data) + &
                    scattering_ave(i_b, 1:n_data))
                scaling = MAXVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
              CASE("Moment   ")
                actual = phf_fnc_ave(i_b, i_term, 1:n_data)
                scaling = 1.0_RealK
            END SELECT
            actual = scaling * actual
!
          CASE (IP_drop_Pade_2_PHF)
!  
!           Scale the dimension:
            scaling_d  = REAL(n_data, RealK) / SUM(d(1:n_data))
            ds = scaling_d * d(1:n_data)
!           Repeat for the optical propeties
            SELECT CASE(property(i))
              CASE("Extinction")
                actual = (absorption_ave(i_b, 1:n_data) + &
                          scattering_ave(i_b, 1:n_data) ) / &
                         (vol_frac(1:n_data) * density_particle(1:n_data))
                scaling = REAL(n_data, RealK) / SUM(actual)
              CASE("Coalbedo  ")
                actual = absorption_ave(i_b, 1:n_data) / &
                  ( absorption_ave(i_b, 1:n_data) + &
                    scattering_ave(i_b, 1:n_data))
                scaling = MAXVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
              CASE("Moment   ")
                actual = phf_fnc_ave(i_b, i_term, 1:n_data)
                scaling = 1.0_RealK
            END SELECT
            actual = scaling * actual
!
        END SELECT
!  
      CASE("Ice  ")
        SELECT CASE (i_fit)
          CASE (IP_Slingo_Schrecker_ice)
!  
!           Scale the dimension:
            scaling_d  = REAL(n_data, RealK) / SUM(d(1:n_data))
            ds = scaling_d * d(1:n_data)
!           Repeat for the optical propeties
            SELECT CASE(property(i))
              CASE("Extinction")
                actual = (absorption_ave(i_b, 1:n_data) + &
                          scattering_ave(i_b, 1:n_data) ) / &
                         (vol_frac(1:n_data) * density_particle(1:n_data))
                scaling = REAL(n_data, RealK) / SUM(actual)
              CASE("Coalbedo  ")
                actual = absorption_ave(i_b, 1:n_data) / &
                  ( absorption_ave(i_b, 1:n_data) + &
                    scattering_ave(i_b, 1:n_data))
                scaling = MAXVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
              CASE("Asymmetry ")
                actual = phf_fnc_ave(i_b, 1, 1:n_data)
                scaling = 1.0_RealK
            END SELECT
            actual = scaling * actual
!  
          CASE (IP_ice_adt)
!  
!           Scale the dimension:
            scaling_d = 1.0 / SQRT(MAXVAL(d(1:n_data)) * MINVAL(d(1:n_data)))
            ds = LOG(d(1:n_data) * scaling_d)
!           Repeat for the optical properties
            SELECT CASE(property(i))
              CASE("Extinction")
                actual = (absorption_ave(i_b, 1:n_data) + &
                          scattering_ave(i_b, 1:n_data) ) / &
                         (vol_frac(1:n_data) * density_particle(1:n_data))
                scaling = REAL(n_data, RealK) / SUM(actual)
              CASE("Coalbedo  ")
                actual = absorption_ave(i_b, 1:n_data) / &
                  ( absorption_ave(i_b, 1:n_data) + &
                    scattering_ave(i_b, 1:n_data))
                scaling = MAXVAL(actual) - MINVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
              CASE("Asymmetry ")
                actual = phf_fnc_ave(i_b, 1, 1:n_data)
                scaling = MAXVAL(actual) - MINVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
            END SELECT
            actual = scaling * actual
!  
          CASE (IP_ice_adt_10)
!  
!           Scale the dimension:
            scaling_d = 1.0 / MAXVAL(d(1:n_data))
            ds = d(1:n_data) * scaling_d
!           Repeat for the optical properties
            SELECT CASE(property(i))
              CASE("Extinction")
                actual = (absorption_ave(i_b, 1:n_data) + &
                          scattering_ave(i_b, 1:n_data) ) / &
                         (vol_frac(1:n_data) * density_particle(1:n_data))
                scaling = REAL(n_data, RealK) / SUM(actual)
              CASE("Coalbedo  ")
                actual = absorption_ave(i_b, 1:n_data) / &
                  ( absorption_ave(i_b, 1:n_data) + &
                    scattering_ave(i_b, 1:n_data))
                scaling = MAXVAL(actual) - MINVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
              CASE("Asymmetry ")
                actual = phf_fnc_ave(i_b, 1, 1:n_data)
                scaling = MAXVAL(actual) - MINVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
            END SELECT
            actual = scaling * actual
!  
          CASE (IP_ice_Fu_IR)
!  
!           Scale the dimension:
            scaling_d = 1.0 / MAXVAL(d(1:n_data))
            ds = d(1:n_data) *scaling_d
!           Repeat for the optical propeties
            SELECT CASE(property(i))
              CASE("Extinction")
                actual = (absorption_ave(i_b, 1:n_data) + &
                          scattering_ave(i_b, 1:n_data) ) / &
                         (vol_frac(1:n_data) * density_particle(1:n_data))
                scaling = REAL(n_data, RealK) / SUM(actual)
              CASE("Coalbedo  ")
                actual = absorption_ave(i_b, 1:n_data) / &
                  ( absorption_ave(i_b, 1:n_data) + &
                    scattering_ave(i_b, 1:n_data))
                scaling = MAXVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
              CASE("Asymmetry ")
                actual = phf_fnc_ave(i_b, 1, 1:n_data)
                scaling = MAXVAL(actual) - MINVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
            END SELECT
            actual = scaling * actual
!  
          CASE (IP_Slingo_Schr_ice_PHF)
!  
!           Scale the dimension:
            scaling_d = REAL(n_data, RealK)/ SUM(d(1:n_data))
            ds = scaling_d * d(1:n_data)
!           Repeat for the optical propeties
            SELECT CASE(property(i))
              CASE("Extinction")
                actual = (absorption_ave(i_b, 1:n_data) + &
                          scattering_ave(i_b, 1:n_data) ) / &
                         (vol_frac(1:n_data) * density_particle(1:n_data))
                scaling = REAL(n_data, RealK) / SUM(actual)
              CASE("Coalbedo  ")
                actual = absorption_ave(i_b, 1:n_data) / &
                  ( absorption_ave(i_b, 1:n_data) + &
                    scattering_ave(i_b, 1:n_data))
                scaling = MAXVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
              CASE("Moment   ")
                actual = phf_fnc_ave(i_b, i_term, 1:n_data)
                scaling = 1.0_RealK
            END SELECT
            actual = scaling * actual
!  
          CASE (IP_ice_Fu_PHF)
!  
!           Scale the dimension:
            scaling_d = 1.0 / MAXVAL(d(1:n_data))
            ds = d(1:n_data) *scaling_d
!           Repeat for the optical propeties
            SELECT CASE(property(i))
              CASE("Extinction")
                actual = (absorption_ave(i_b, 1:n_data) + &
                          scattering_ave(i_b, 1:n_data) ) / &
                         (vol_frac(1:n_data) * density_particle(1:n_data))
                scaling = REAL(n_data, RealK) / SUM(actual)
              CASE("Coalbedo  ")
                actual = absorption_ave(i_b, 1:n_data) / &
                  ( absorption_ave(i_b, 1:n_data) + &
                    scattering_ave(i_b, 1:n_data))
                scaling = MAXVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
              CASE("Moment   ")
                actual = phf_fnc_ave(i_b, i_term, 1:n_data)
                scaling = MAXVAL(actual) - MINVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
            END SELECT
            actual = scaling * actual
!
          CASE (IP_Ice_Pade_2_PHF)
!  
!           Scale the dimension:
            scaling_d  = REAL(n_data, RealK) / SUM(d(1:n_data))
            ds = scaling_d * d(1:n_data)
!           Repeat for the optical propeties
            SELECT CASE(property(i))
              CASE("Extinction")
                actual = (absorption_ave(i_b, 1:n_data) + &
                          scattering_ave(i_b, 1:n_data) ) / &
                         (vol_frac(1:n_data) * density_particle(1:n_data))
                scaling = REAL(n_data, RealK) / SUM(actual)
              CASE("Coalbedo  ")
                actual = absorption_ave(i_b, 1:n_data) / &
                  ( absorption_ave(i_b, 1:n_data) + &
                    scattering_ave(i_b, 1:n_data))
                scaling = MAXVAL(actual)
                IF (scaling > 0.0) THEN
                  scaling = 1.0 / scaling
                ELSE
                  scaling = 1.0
                ENDIF
              CASE("Moment   ")
                actual = phf_fnc_ave(i_b, i_term, 1:n_data)
                scaling = 1.0_RealK
            END SELECT
            actual = scaling * actual
!
      END SELECT
!
    END SELECT
!
!
!
  END SUBROUTINE condition_data
!
!
!
  SUBROUTINE rescale_parm(parm)

    USE PolynomialRoots, ONLY: CubicRoots, QuadraticRoots

    IMPLICIT NONE

    INTEGER, PARAMETER           :: dp=KIND(1.0D0)
    INTEGER                      :: i_root
    REAL (RealK), Dimension(:)   :: parm
    REAL (dp),    Dimension(0:3) :: poly
    COMPLEX (dp), Dimension(3)   :: roots
!     Variables used in the calculation of cubic and quadratic roots


    SELECT CASE(fit_species)
!
      CASE ("Water")
!
        SELECT CASE (i_fit)
!
          CASE (IP_Slingo_Schrecker)
!    
            SELECT CASE(property(i))
              CASE("Extinction")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) / (scaling * scaling_d)
              CASE("Coalbedo  ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
              CASE("Asymmetry ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
            END SELECT
!    
          CASE (IP_Ackerman_Stephens)
!    
            SELECT CASE(property(i))
              CASE("Extinction")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * EXP(parm(3) * LOG(scaling_d) ) / &
                  scaling
              CASE("Coalbedo  ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * EXP(parm(3) * LOG(scaling_d) ) / &
                  scaling
              CASE("Asymmetry ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * EXP(parm(3) * LOG(scaling_d) ) / &
                  scaling
            END SELECT
!    
          CASE (IP_drop_Pade_2)
!    
            SELECT CASE(property(i))
              CASE("Extinction")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
                parm(3) = parm(3) * scaling_d * scaling_d / scaling
                parm(4) = parm(4) * scaling_d
                parm(5) = parm(5) * scaling_d * scaling_d
                parm(6) = parm(6) * scaling_d * scaling_d * scaling_d
                ! Check these parameters will not produce NaNs
                poly(0)   = REAL(1.0, dp)
                poly(1:3) = REAL(parm(4:6), dp)
                CALL CubicRoots(poly, roots)
                DO i_root = 1, 3
                  IF (ABS(AIMAG(roots(i_root))) < TINY(1.0d0)) THEN
                    IF (DBLE(roots(i_root)) > MINVAL(d(1:n_data)) .AND. &
                        DBLE(roots(i_root)) < MAXVAL(d(1:n_data))) THEN
                      WRITE(iu_err, '(a,i0,a,e12.5)') &
                        'Warning: Extinction will be NaN for band ',i_b, &
                        ' with effective radius of',DBLE(roots(i_root))
                    END IF
                  END IF
                END DO
              CASE("Coalbedo  ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
                parm(3) = parm(3) * scaling_d * scaling_d / scaling
                parm(4) = parm(4) * scaling_d
                parm(5) = parm(5) * scaling_d * scaling_d
                poly(0)   = REAL(1.0, dp)
                poly(1:2) = REAL(parm(4:5), dp)
                CALL QuadraticRoots(poly, roots)
                DO i_root = 1, 2
                  IF (ABS(AIMAG(roots(i_root))) < TINY(1.0d0)) THEN
                    IF (DBLE(roots(i_root)) > MINVAL(d(1:n_data)) .AND. &
                        DBLE(roots(i_root)) < MAXVAL(d(1:n_data))) THEN
                      WRITE(iu_err, '(a,i0,a,e12.5)') &
                        'Warning: Coalbedo will be NaN for band ',i_b, &
                        ' with effective radius of',DBLE(roots(i_root))
                    END IF
                  END IF
                END DO
              CASE("Asymmetry ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
                parm(3) = parm(3) * scaling_d * scaling_d / scaling
                parm(4) = parm(4) * scaling_d
                parm(5) = parm(5) * scaling_d * scaling_d
                poly(0)   = REAL(1.0, dp)
                poly(1:2) = REAL(parm(4:5), dp)
                CALL QuadraticRoots(poly, roots)
                DO i_root = 1, 2
                  IF (ABS(AIMAG(roots(i_root))) < TINY(1.0d0)) THEN
                    IF (DBLE(roots(i_root)) > MINVAL(d(1:n_data)) .AND. &
                        DBLE(roots(i_root)) < MAXVAL(d(1:n_data))) THEN
                      WRITE(iu_err, '(a,i0,a,e12.5)') &
                        'Warning: Asymmetry will be NaN for band ',i_b, &
                        ' with effective radius of',DBLE(roots(i_root))
                    END IF
                  END IF
                END DO
            END SELECT
!    
          CASE (IP_Slingo_Schr_PHF)
!    
            SELECT CASE(property(i))
              CASE("Extinction")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) / (scaling * scaling_d)
              CASE("Coalbedo  ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
              CASE("Moment    ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
            END SELECT
!    
          CASE (IP_drop_Pade_2_PHF)
!    
            SELECT CASE(property(i))
              CASE("Extinction")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
                parm(3) = parm(3) * scaling_d * scaling_d / scaling
                parm(4) = parm(4) * scaling_d
                parm(5) = parm(5) * scaling_d * scaling_d
                parm(6) = parm(6) * scaling_d * scaling_d * scaling_d
              CASE("Coalbedo  ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
                parm(3) = parm(3) * scaling_d * scaling_d / scaling
                parm(4) = parm(4) * scaling_d
                parm(5) = parm(5) * scaling_d * scaling_d
              CASE("Moment    ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
                parm(3) = parm(3) * scaling_d * scaling_d / scaling
                parm(4) = parm(4) * scaling_d
                parm(5) = parm(5) * scaling_d * scaling_d
            END SELECT
!    
!    
        END SELECT
!
      CASE ("Ice  ")
!
        SELECT CASE (i_fit)
          CASE (IP_Slingo_Schrecker_ice)
!    
            SELECT CASE(property(i))
              CASE("Extinction")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) / (scaling * scaling_d)
              CASE("Coalbedo  ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
              CASE("Asymmetry ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
            END SELECT
!    
          CASE (IP_ice_adt)
!    
            SELECT CASE(property(i))
              CASE("Extinction")
                parm(1) = parm(1) - LOG(scaling)
                parm(10) = 1.0 / scaling_d
              CASE("Coalbedo  ")
                parm(1:9) = parm(1:9) / scaling
                parm(10) = 1.0 / scaling_d
              CASE("Asymmetry ")
                parm(1:9) = parm(1:9) / scaling
                parm(10) = 1.0 / scaling_d
            END SELECT
!    
          CASE (IP_ice_adt_10)
!    
            SELECT CASE(property(i))
              CASE("Extinction")
                parm(1) = parm(1) - LOG(scaling)
                parm(12) = 1.0 / scaling_d
              CASE("Coalbedo  ")
                parm(1:11) = parm(1:11) / scaling
                parm(12) = 1.0 / scaling_d
              CASE("Asymmetry ")
                parm(1:11) = parm(1:11) / scaling
                parm(12) = 1.0 / scaling_d
            END SELECT
!    
          CASE (IP_Slingo_Schr_ice_PHF)
!    
            SELECT CASE(property(i))
              CASE("Extinction")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) / (scaling * scaling_d)
              CASE("Coalbedo  ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
              CASE("Moment    ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
            END SELECT
!    
          CASE (IP_ice_Fu_PHF)
!    
            SELECT CASE(property(i))
              CASE("Extinction")
                parm(1:3) = parm(1:3) / scaling
                parm(4) = 1.0 / scaling_d
              CASE("Coalbedo  ")
                parm(1:4) = parm(1:4) / scaling
                parm(5) = 1.0 / scaling_d
              CASE("Moment   ")
                parm(1:4) = parm(1:4) / scaling
                parm(5) = 1.0 / scaling_d
            END SELECT
!    
          CASE (IP_Ice_Pade_2_PHF)
!    
            SELECT CASE(property(i))
              CASE("Extinction")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
                parm(3) = parm(3) * scaling_d * scaling_d / scaling
                parm(4) = parm(4) * scaling_d
                parm(5) = parm(5) * scaling_d * scaling_d
                parm(6) = parm(6) * scaling_d * scaling_d * scaling_d
                ! Check these parameters will not produce NaNs
                poly(0)   = REAL(1.0, dp)
                poly(1:3) = REAL(parm(4:6), dp)
                CALL CubicRoots(poly, roots)
                DO i_root = 1, 3
                  IF (ABS(AIMAG(roots(i_root))) < TINY(1.0d0)) THEN
                    IF (DBLE(roots(i_root)) > MINVAL(d(1:n_data)) .AND. &
                        DBLE(roots(i_root)) < MAXVAL(d(1:n_data))) THEN
                      WRITE(iu_err, '(a,i0,a,e12.5)') &
                        'Warning: Extinction will be NaN for band ',i_b, &
                        ' with effective radius of',DBLE(roots(i_root))
                    END IF
                  END IF
                END DO
              CASE("Coalbedo  ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
                parm(3) = parm(3) * scaling_d * scaling_d / scaling
                parm(4) = parm(4) * scaling_d
                parm(5) = parm(5) * scaling_d * scaling_d
                poly(0)   = REAL(1.0, dp)
                poly(1:2) = REAL(parm(4:5), dp)
                CALL QuadraticRoots(poly, roots)
                DO i_root = 1, 2
                  IF (ABS(AIMAG(roots(i_root))) < TINY(1.0d0)) THEN
                    IF (DBLE(roots(i_root)) > MINVAL(d(1:n_data)) .AND. &
                        DBLE(roots(i_root)) < MAXVAL(d(1:n_data))) THEN
                      WRITE(iu_err, '(a,i0,a,e12.5)') &
                        'Warning: Coalbedo will be NaN for band ',i_b, &
                        ' with effective radius of',DBLE(roots(i_root))
                    END IF
                  END IF
                END DO
              CASE("Asymmetry ","Moment    ")
                parm(1) = parm(1) / scaling
                parm(2) = parm(2) * scaling_d / scaling
                parm(3) = parm(3) * scaling_d * scaling_d / scaling
                parm(4) = parm(4) * scaling_d
                parm(5) = parm(5) * scaling_d * scaling_d
                poly(0)   = REAL(1.0, dp)
                poly(1:2) = REAL(parm(4:5), dp)
                CALL QuadraticRoots(poly, roots)
                DO i_root = 1, 2
                  IF (ABS(AIMAG(roots(i_root))) < TINY(1.0d0)) THEN
                    IF (DBLE(roots(i_root)) > MINVAL(d(1:n_data)) .AND. &
                        DBLE(roots(i_root)) < MAXVAL(d(1:n_data))) THEN
                      WRITE(iu_err, '(a,i0,a,e12.5)') &
                        'Warning: Asymmetry will be NaN for band ',i_b, &
                        ' with effective radius of',DBLE(roots(i_root))
                    END IF
                  END IF
                END DO
            END SELECT
!    
        END SELECT
    END SELECT
!
!
!
  END SUBROUTINE rescale_parm
!
!
!
END SUBROUTINE cloud_fit_90


