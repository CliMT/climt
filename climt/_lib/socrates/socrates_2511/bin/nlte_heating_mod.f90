! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculation of NLTE atmospheric heating rates

MODULE nlte_heating_mod

IMPLICIT NONE

! Description:
!   Calculation of the NLTE heating rates according to the Fomichev
!   matrix and recurrence parameterisations.
!
! Method:
!   Method follows these papers by Fomichev:
!
!   1. Fomichev, V. I., Blanchet J.-P., and Turner D. S. (1998), Matrix 
!   parameterization of the 15 $\mu$m CO$_2$ band cooling in the middle and 
!   upper atmosphere for variable CO$_2$ concentration, J. Geophys. Res., 
!   {\bf 103}, No. D10, 11505 - 11528.
! 
!   2. Fomichev, V. I., and Blanchet J.-P. (1995), Development of the new 
!   CCC/GCM radiation model for extension into the Middle Atmosphere, 
!   Atmosphere-Ocean, {\bf 33}, No. 3, 513-529.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NLTE_HEATING_MOD'

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE nlte_heating_lw(t, p, co2_mix_ratio, o3_mix_ratio, heat_rate, &
                           n_layer, n_profile)  

  USE realtype_rd, ONLY: RealK
  USE nlte_heating_data_mod, ONLY: n_level, sec_per_day, step_size, amm,       &
                                   matrix_min, o3_matrix_max, gas_const,       &
                                   molmass_o3, molmass_co2,                    &
                                   planck_exp_const_co2, planck_exp_const_o3,  &
                                   polyatom_max, vmr_o
  USE interpolate_p_mod, ONLY: ip_1_lin_lin, interpolate_p

  ! Imported variables
  INTEGER, INTENT(IN) :: n_layer, n_profile
  REAL (RealK), INTENT(IN) :: t(:, :), p(:, :)
  REAL (RealK), INTENT(IN) :: co2_mix_ratio(:, :)
  REAL (RealK), INTENT(IN) :: o3_mix_ratio(:, :)
  ! (Note: no mass mixing ratio for atomic oxygen available)

  ! Heating rate (passed in as LTE, returned with correction)
  REAL (RealK), INTENT(INOUT) :: heat_rate(:, :)

  ! Local variables

  ! Scale height, and pressure, temperature, gas mixing ratios and heating 
  ! rates corresponding to the scale height
  REAL (RealK), DIMENSION(n_level) :: x, p_x, t_x
  REAL (RealK), DIMENSION(n_level) :: gas_mix_ratio_x
  REAL (RealK), DIMENSION(n_level) :: hr_x = 0.0d0

  ! NLTE heating rate 
  REAL (RealK), DIMENSION(n_profile, n_layer) :: heat_rate_nlte

  ! Splining variables
  REAL (RealK), ALLOCATABLE, DIMENSION(:) :: xspline, yspline, yspline2
  REAL (RealK), DIMENSION(n_level) :: xspline_lev = 0.0d0, &
                                      yspline_lev = 0.0d0, yspline2_lev = 0.0d0

  ! Volume mixing ratios
  REAL (RealK), DIMENSION(n_level) :: vmr_co2 = 0.0d0, vmr_o3 = 0.0d0

  ! Base volume mixing ratio for CO2
  REAL (RealK) :: vmr_co2_base

  ! Planck exponent
  REAL (RealK), DIMENSION(n_level) :: planck_exp_co2 = 0.0d0, &
                                      planck_exp_o3 = 0.0d0

  ! Specific heat
  REAL (RealK) :: spec_heat

  ! Indexing variables and misc
  INTEGER :: i, j
  LOGICAL :: l_splined

  ALLOCATE(xspline(n_layer))
  ALLOCATE(yspline(n_layer))
  ALLOCATE(yspline2(n_layer))
  xspline(:) = 0.0d0
  yspline(:) = 0.0d0
  yspline2(:) = 0.0d0

  ! Create a scale height array 0 <= x <= 17.5 (top of atmosphere to bottom)
  ! in steps of delta-x = 0.25
  x = (/ (i*step_size, i = 0, n_level-1) /)

  ! Pressure corresponding to the scale height in Pascal
  DO i = 1, n_level
    p_x(i) = 100000./exp(x(i))
  END DO

  ! Reverse the pressure array
  p_x = p_x(n_level:1:-1)

  ! Loop over the profiles
  DO j = 1, n_profile

     ! Interpolate the temperature to the new scale height. 
     ! Note: needs some error handling here.
     DO i = 1, n_level
        l_splined = .false.
        CALL interpolate_p(n_layer, p(j, :), t(j, :), xspline, &
                           yspline, yspline2, p_x(i), t_x(i), IP_1_lin_lin, &
                           l_splined)
     END DO

     ! Reverse the temperature array
     t_x = t_x(n_level:1:-1)

     ! Interpolate the CO2 gas mixing ratios to the new scale height
     ! Note: needs some error handling here.
     DO i = 1, n_level
        l_splined = .false.
        CALL interpolate_p(n_layer, p(j, :), co2_mix_ratio(j, :), &
                           xspline, yspline, yspline2, p_x(i), &
                           gas_mix_ratio_x(i), IP_1_lin_lin, l_splined)
     END DO
     ! Reverse the gas mixing arrays
     gas_mix_ratio_x = gas_mix_ratio_x(n_level:1:-1)

     ! Convert mass mixing ratios (mmr) to volume mixing ratios (vmr)
     ! NOTE: Fomichev uses the fact that the mmr of CO2 is constant below 
     ! x=12.5. Here this uses the mmr directly from the UM.
     DO i = matrix_min, n_level
        vmr_co2(i) = gas_mix_ratio_x(i) * amm(i) / molmass_co2
     END DO

     ! A base CO2 volume mixing ratio is set up in ppm (needed for the 
     ! interpolation)
     vmr_co2_base = vmr_co2(matrix_min) * 1.0e6

     ! Interpolate the O3 gas mixing ratios to the new scale height
     DO i = 1, n_level
        l_splined = .false.
        CALL interpolate_p(n_layer, p(j, :), o3_mix_ratio(j, :), &
                           xspline, yspline, yspline2, p_x(i), &
                           gas_mix_ratio_x(i), IP_1_lin_lin, l_splined)
     END DO
     ! Reverse the gas mixing arrays
     gas_mix_ratio_x = gas_mix_ratio_x(n_level:1:-1)

     ! O3 needed for matrix parameterisation only, up to x=10.5
     DO i = matrix_min, o3_matrix_max
        vmr_o3(i) = gas_mix_ratio_x(i) * amm(i) / molmass_o3
     END DO
    
     ! Set up part of Planck equation
     ! Should I be splitting this to only pass in the values required?
     DO i = 1, n_level
        planck_exp_co2(i) = exp(-1.0 * planck_exp_const_co2 / t_x(i))
        planck_exp_o3(i) = exp(-1.0 * planck_exp_const_o3 / t_x(i))
     END DO

     ! Heating rates from the matrix parameterisation
     CALL hr_matrix_param(vmr_o3, planck_exp_co2, planck_exp_o3, hr_x, &
                          vmr_co2_base)

     ! Heating rates from the recurrence parameterisation
     CALL hr_recurr_param(x, t_x, vmr_co2, planck_exp_co2, hr_x)

     ! Convert heating rates into K/day
     ! For x <= 10.25 ignore monatomic gases 
     DO i = matrix_min, polyatom_max
        spec_heat = 3.5 * gas_const
        hr_x(i) = hr_x(i) * amm(i) * sec_per_day / spec_heat
     END DO
     ! x >= 10.5, consider the contribution of atomic Oxygen to the specific
     ! heat capacity
     DO i = polyatom_max+1, n_level
        spec_heat = gas_const * ((3.5 * (1.0 - vmr_o(i))) + (2.5 * vmr_o(i)))
        hr_x(i) = hr_x(i) * amm(i) * sec_per_day / spec_heat
     END DO

     ! Reverse the heating rate array
     hr_x = hr_x(n_level:1:-1)

     ! Interpolate heating rates to the original pressure levels
     DO i = 1, n_layer
        l_splined = .false.
        CALL interpolate_p(n_level, p_x, hr_x, xspline_lev, &
                           yspline_lev, yspline2_lev, p(j, i), &
                           heat_rate_nlte(j, i), IP_1_lin_lin, l_splined)
     END DO

     ! Blend the heating rates with the LTE rates
     DO i = 1, n_layer
        ! Boundary condition at 0.1 hPa = 10 Pa
        IF (p(j,i) < 10.) THEN 
           heat_rate(j, i) = heat_rate_nlte(j, i)
        END IF
     END DO

  END DO ! End of loop over n_profile

  DEALLOCATE(xspline)
  DEALLOCATE(yspline)
  DEALLOCATE(yspline2)

END SUBROUTINE nlte_heating_lw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE hr_matrix_param(vmr_o3, planck_exp_co2, planck_exp_o3, hr_x, &
                           vmr_co2_base)

  USE realtype_rd, ONLY: RealK
  USE nlte_heating_data_mod, ONLY: n_level, size_coeff, matrix_min, matrix_max,&
                                   matrix_coeff_co2_j, matrix_coeff_co2_conc,  &
                                   matrix_coeff_co2, o3_matrix_max,            &
                                   matrix_coeff_o3_j, delta_xj_co2_index,      &
                                   delta_xj_o3_index, matrix_coeff_o3
  USE interpolate_p_mod, ONLY: ip_1_lin_log, ip_3_lin_lin, interpolate_p
 
  ! Imported variables
  REAL (RealK), DIMENSION(n_level), INTENT(IN) :: vmr_o3
  REAL (RealK), DIMENSION(n_level), INTENT(IN) :: planck_exp_co2, planck_exp_o3
  REAL (RealK), DIMENSION(n_level), INTENT(INOUT) :: hr_x
  REAL (RealK), INTENT(IN) :: vmr_co2_base

  ! Variables for calculating the coefficients
  REAL (RealK), DIMENSION(4) :: coeff_forinterp
  REAL (RealK), DIMENSION(2) :: corr_coeff_co2

  ! Splining variables
  REAL (RealK), DIMENSION(size_coeff) :: xspline_coeff = 0.0d0, &
                                         yspline_coeff = 0.0d0, &
                                         yspline2_coeff = 0.0d0

  ! Indexing variables and misc
  INTEGER :: i, j, k, l
  LOGICAL :: l_splined
  INTEGER :: interp_mode, interp_mode_check

  ! Finding the heating rates per dimensionless scale height using the 
  ! matrix method for 2 < x < 12.5
  DO i = matrix_min, matrix_max

     ! Loop over the height levels for which the coeffs are calculated
     DO j = 1, size(matrix_coeff_co2_j)

        ! This is to take care of the coefficients of zero in the data
        ! tables which cause the interpolation routine to crash
        IF ((i <= 13) .AND. (j == 2)) THEN
           corr_coeff_co2(1) = 0.0d0
           corr_coeff_co2(2) = 0.0d0
        ELSE
       
           ! Loop over the coefficient types
           DO l = 1, 2

              ! Find the mode of interpolation (0 = linear negative, 1 = linear
              ! positive, 2 = spline)
              IF (matrix_coeff_co2(l, 1, j, i) < 0) THEN 
                 interp_mode = 0
              ELSE 
                 interp_mode = 1
              END IF
              DO k = 2, 4
                 IF (matrix_coeff_co2(l, k, j, i) < 0) THEN 
                    interp_mode_check = 0
                 ELSE 
                    interp_mode_check = 1
                 END IF
                 ! Use spline interpolation if the signs don't match
                 IF (interp_mode_check /= interp_mode) THEN 
                    interp_mode = 2
                    EXIT
                 END IF
              END DO

              ! Putting the coefficients in the correct form for interpolation
              DO k = 1, 4
                 coeff_forinterp(k) = matrix_coeff_co2(l, k, j, i) / &
                                      matrix_coeff_co2_conc(k)
              END DO

              ! If the signs are all the same, then use linear interpolation
              ! of the logs of the coefficients
              IF (interp_mode == 0 .OR. interp_mode == 1) THEN

                 IF (interp_mode == 0) THEN             
                    DO k = 1, 4
                       coeff_forinterp(k) = -1.0 * coeff_forinterp(k)
                    END DO
                 END IF

                 l_splined = .false.

                 CALL interpolate_p(size(matrix_coeff_co2_conc), &
                                    matrix_coeff_co2_conc, coeff_forinterp, &
                                    xspline_coeff, yspline_coeff, &
                                    yspline2_coeff, vmr_co2_base, &
                                    corr_coeff_co2(l), IP_1_lin_log, l_splined)
                 IF (interp_mode == 0) THEN
                    corr_coeff_co2(l) = -1.0 * corr_coeff_co2(l)
                 END IF

              ! If the signs are different, use a spline interpolation
              ELSE IF (interp_mode == 2) THEN

                 ! Is this actually a sensible interpolation routine to use?
                 l_splined = .false.
                 CALL interpolate_p(size(matrix_coeff_co2_conc), &
                                    matrix_coeff_co2_conc, coeff_forinterp, &
                                    xspline_coeff, yspline_coeff, &
                                    yspline2_coeff, vmr_co2_base, &
                                    corr_coeff_co2(l), IP_3_lin_lin, l_splined)
              END IF

              corr_coeff_co2(l) = corr_coeff_co2(l) * vmr_co2_base
           END DO ! End of loop over l
        END IF ! End of check i/j value check

        ! Calculate the CO2 heating rate contribution using the matrix approach
        ! If the height levels x_j are less than zero, these are taken as zero
        ! and do not contribute to the heating rate
        IF (i <= matrix_min + 4) THEN
          IF (j == 1) THEN
            hr_x(i) = planck_exp_co2(1)                                       &
                    * (corr_coeff_co2(1)                                      &
                       + corr_coeff_co2(2) * planck_exp_co2(i))
          END IF
          IF (j >= 3) THEN
            k = i + delta_xj_co2_index(j)
            hr_x(i) = hr_x(i) + planck_exp_co2(k)                             &
                              * (corr_coeff_co2(1)                            &
                                 + corr_coeff_co2(2) * planck_exp_co2(i))
          END IF
        ELSE IF (i <= matrix_min + 17) THEN
          IF (j == 1) THEN
            hr_x(i) = planck_exp_co2(1)                                       &
                    * (corr_coeff_co2(1)                                      &
                       + corr_coeff_co2(2) * planck_exp_co2(i))
          END IF
          IF (j >= 2) THEN
            k = i + delta_xj_co2_index(j)
            hr_x(i) = hr_x(i) + planck_exp_co2(k)                             &
                              * (corr_coeff_co2(1)                            &
                                 + corr_coeff_co2(2) * planck_exp_co2(i))
          END IF
        ELSE IF (i <= matrix_min + 42) THEN
          k = i + delta_xj_co2_index(j)
          hr_x(i) = hr_x(i) + planck_exp_co2(k)                               &
                            * (corr_coeff_co2(1)                              &
                               + corr_coeff_co2(2) * planck_exp_co2(i))
        END IF

     END DO ! End of loop over j

     ! Normalisation by 1.0e-4 because Fomichev calculates in units of 
     ! cm^2s^-3 , which is equivalent to 10e4 Jkg^-1s^-1
     hr_x(i) = hr_x(i) * 1.0e-4

  END DO ! End of loop over i

  ! Calculate the O3 heating rate contribution up to x=10.5
  DO i = matrix_min, o3_matrix_max
    DO j = 1, size(matrix_coeff_o3_j)
      IF (i <= matrix_min + 4) THEN
        IF (j == 1) THEN
          hr_x(i) = hr_x(i) &
            + (planck_exp_o3(1) * matrix_coeff_o3(1, i) * 1.0e-4 * vmr_o3(i))
        ELSE IF (j >= 3) THEN
          k = i + delta_xj_o3_index(j)
          hr_x(i) = hr_x(i) &
            + (planck_exp_o3(k) * matrix_coeff_o3(j, i) * 1.0e-4 * vmr_o3(i))
        END IF
      ELSE IF (i <= matrix_min + 17) THEN
        IF (j == 1) THEN
          hr_x(i) = hr_x(i) &
            + (planck_exp_o3(1) * matrix_coeff_o3(1, i) * 1.0e-4 * vmr_o3(i))
        ELSE IF (j >= 2) THEN
          k = i + delta_xj_o3_index(j)
          hr_x(i) = hr_x(i) &
            + (planck_exp_o3(k) * matrix_coeff_o3(j, i) * 1.0e-4 * vmr_o3(i))
        END IF
      ELSE IF (i <= matrix_min + 34) THEN
        k = i + delta_xj_o3_index(j)
        hr_x(i) = hr_x(i) &
          + (planck_exp_o3(k) * matrix_coeff_o3(j, i) * 1.0e-4 * vmr_o3(i))
      END IF
    END DO
  END DO

END SUBROUTINE hr_matrix_param


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE hr_recurr_param(x, t_x, vmr_co2, planck_exp_co2, hr_x)

  USE realtype_rd, ONLY: RealK
  USE nlte_heating_data_mod, ONLY: n_level, size_esc, size_alpha, boltz,       &
                                   cdr_o, co2_col_const, co2_col_top,          &
                                   einstein_a, amm, recurr_const, vmr_o,       &
                                   vmr_o2, matrix_esc_fcn, matrix_alpha_dep,   &
                                   vmr_n2
  USE interpolate_p_mod, ONLY: ip_1_lin_lin, interpolate_p

  ! Imported variables
  REAL (RealK), DIMENSION(n_level), INTENT(IN) :: x, t_x
  REAL (RealK), DIMENSION(n_level), INTENT(IN) :: vmr_co2
  REAL (RealK), DIMENSION(n_level), INTENT(IN) :: planck_exp_co2
  REAL (RealK), DIMENSION(n_level), INTENT(INOUT) :: hr_x

  ! Indexing variables
  INTEGER, PARAMETER :: recurr_arr_size = 17
  INTEGER, PARAMETER :: extend_arr_size = 21
  INTEGER :: recurr_skip = 50 ! Corresponds to x=12.25
  INTEGER :: top_index
  INTEGER :: alpha_arr_size = 6
  INTEGER :: high_skip

  ! CO2 column amount
  REAL (RealK), DIMENSION(recurr_arr_size) :: col_amount

  ! Variables for the recurrence coefficients
  REAL (RealK), DIMENSION(recurr_arr_size) :: recurr_coeff_dj
  REAL (RealK) :: col_dep_x_upper, col_dep_x_lower, delta_col_am
  REAL (RealK) :: alpha_dj
  REAL (RealK), DIMENSION(size_esc) :: xspline_esc = 0.0d0, &
                                       yspline_esc = 0.0d0, &
                                       yspline2_esc = 0.0d0
  REAL (RealK), DIMENSION(size_alpha) :: xspline_alpha = 0.0d0, &
                                         yspline_alpha = 0.0d0, &
                                         yspline2_alpha = 0.0d0

  ! Variables for calculating lambda
  REAL (RealK), DIMENSION(extend_arr_size) :: lambda
  REAL (RealK) :: number_density, t_cube_rt, cdr_n2, cdr_o2

  ! Variables for calculating gamma
  REAL (RealK) :: coeff_d_low, coeff_d_high, gamma_low, gamma_high
  
  ! Variables for calculating the heating rate
  REAL (RealK) :: up_flux

  ! Indexing variables and misc
  INTEGER :: i, j
  LOGICAL :: l_splined

  ! Calculate the CO2 column amount for 12.5 < x < 16.5. Note that I've
  ! still got a mystery factor in the co2_col_const (see 
  ! nlte_heating_data_mod.f90)
  col_amount(recurr_arr_size) = co2_col_top
  top_index = recurr_arr_size + recurr_skip
  col_dep_x_upper = vmr_co2(top_index) * exp(-x(top_index)) / amm(top_index) 
  DO i = recurr_arr_size-1, 1, -1
     j = i + recurr_skip
     col_dep_x_lower = vmr_co2(j) * exp(-x(j)) / amm(j) 
     delta_col_am = 0.5 * co2_col_const * (col_dep_x_upper + col_dep_x_lower) 
     col_amount(i) = col_amount(i+1) + delta_col_am
     col_dep_x_upper = col_dep_x_lower
  END DO

  ! Calculate the uncorrected escape function coefficients in the range
  ! 12.5 <= x <= 16.5
  DO i = 1, recurr_arr_size 
     l_splined = .false.
     CALL interpolate_p(size_esc, matrix_esc_fcn(1, :), &
                        matrix_esc_fcn(2, :), xspline_esc, yspline_esc, &
                        yspline2_esc, col_amount(i), recurr_coeff_dj(i), &
                        IP_1_lin_lin, l_splined)
  END DO

  ! Calculate the correction to the recurrence coefficients in range 
  ! 12.5 <= x <= 13.75 by interpolating to column density at x
  DO i = 1, alpha_arr_size 
     l_splined = .false.
     CALL interpolate_p(size_alpha, matrix_alpha_dep(1, :, i), &
                        matrix_alpha_dep(2, :, i), xspline_alpha, &
                        yspline_alpha, yspline2_alpha, col_amount(i), &
                        alpha_dj, IP_1_lin_lin, l_splined)
     ! These are added as they are in log form
     recurr_coeff_dj(i) = recurr_coeff_dj(i) + alpha_dj
  END DO
  ! Finally convert back from log form
  DO i = 1, recurr_arr_size 
     recurr_coeff_dj(i) = exp(recurr_coeff_dj(i))
  END DO
  
  ! Calculate the quantity lambda up to x=17.5
  DO i = 1, extend_arr_size 
     ! Calculate lambda
     number_density = 1.0e5 * exp(-x(i+recurr_skip)) / &
          (boltz * t_x(i+recurr_skip))
     t_cube_rt = t_x(i+recurr_skip)**(-1./3.)
     ! Collisional deactivation rate constants
     cdr_n2 = 1.0e-6 * ((5.5e-17 * sqrt(t_x(i+recurr_skip))) + &
          (6.7e-10 * exp(-83.8 * t_cube_rt)))
     cdr_o2 = 1.0e-6 * (1.0e-15 * exp(23.37 - (230.9 * t_cube_rt) + & 
          (564.0 * t_cube_rt * t_cube_rt)))
     lambda(i) = einstein_a / (einstein_a & 
          + number_density * ((vmr_n2(i+recurr_skip) * cdr_n2) &
          + (vmr_o2(i+recurr_skip) * cdr_o2) &
          + (vmr_o(i+recurr_skip) * cdr_o)))
  END DO
  
  ! Calculate the boundary condition at x = 12.5 (using the heating 
  ! rate found from the matrix method)
  gamma_low = (hr_x(1+recurr_skip) * amm(1+recurr_skip)) / & 
              (recurr_const * vmr_co2(1+recurr_skip) * (1.0 - lambda(1)))

  ! Calculate the heating rate using the recurrence approach for 
  ! 12.75 <= x <= 16.5
  DO i = 2, recurr_arr_size 
     coeff_d_low = 0.25 * ((3.0 * recurr_coeff_dj(i-1)) + recurr_coeff_dj(i))
     coeff_d_high = 0.25 * (recurr_coeff_dj(i-1) + (3.0 * recurr_coeff_dj(i)))

     gamma_high = (((1.0 - (lambda(i-1) * (1.0 - coeff_d_low))) * gamma_low) &
          + (coeff_d_low * planck_exp_co2(i-1+recurr_skip)) &
          - (coeff_d_high * planck_exp_co2(i+recurr_skip))) &
          / (1.0 - (lambda(i) * (1.0 - coeff_d_high)))

     hr_x(i+recurr_skip) = (recurr_const * vmr_co2(i+recurr_skip) & 
                            * (1.0 - lambda(i)) * gamma_high) & 
                            / amm(i+recurr_skip) 

     gamma_low = gamma_high
  END DO

  ! Calculate the heating rate for x >= 16.75
  high_skip = recurr_skip + recurr_arr_size
  ! Boundary condition at x = 16.5
  up_flux = gamma_high + planck_exp_co2(high_skip)
  DO i = 1, (extend_arr_size-recurr_arr_size)
     hr_x(i+high_skip) = (recurr_const * vmr_co2(i+high_skip) &
                          * (1.0 - lambda(i+recurr_arr_size)) &
                          * (up_flux - planck_exp_co2(i+high_skip))) &
                          / amm(i+high_skip) 
  END DO
     
END SUBROUTINE hr_recurr_param

END MODULE nlte_heating_mod
