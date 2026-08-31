! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculation of SW NLTE atmospheric heating rates

MODULE sw_nlte_heating_mod

USE realtype_rd, ONLY: RealK
IMPLICIT NONE

! Description:
!   Calculation of the SW NLTE heating rate corrections according to the
!   Fomichev parameterisation for CO2.
!
! Method:
!   Method follows these papers by Fomichev:
!
!   3. Ogibalov V. P., and V. I. Fomichev (2003), Parameterization of solar
!   heating by the near IR CO2 bands in the mesosphere, Adv. Space Res., 32,
!   No. 5, 759-764.
!
!   4. Fomichev V. I., V. P. Ogibalov, and S. R. Beagley (2004), Solar
!   heating by the near-IR CO2 bands in the mesosphere, Geophys. Res. Lett.,
!   31, L21102, doi:10.1029/2004GL020324.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SW_NLTE_HEATING_MOD'

! Ozone Hartley Band: 210nm - 320nm (requires NLTE correction)
REAL (RealK), PARAMETER :: hartley_wavelength_min = 2.1E-07_RealK
REAL (RealK), PARAMETER :: hartley_wavelength_max = 3.2E-07_RealK

! Extreme-UV: < 98.6nm (requires NLTE correction)
REAL (RealK), PARAMETER :: euv_wavelength_max = 9.86E-08_RealK

! LTE regions: < 1.1 microns, excluding above regions
REAL (RealK), PARAMETER :: lte_wavelength_max = 1.1E-06_RealK

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE nlte_heating_sw(trop_co2_mmr, coszin, p, heat_rate, &
                           heat_rate_lte_bands, heat_rate_hartley_band, &
                           heat_rate_euv_bands, &
                           n_layer, n_profile)

  USE nlte_heating_data_mod, ONLY: n_level, n_x_levs_min, n_co2_levs,          &
                                   n_x_levs_max, n_x_levs, param_co2_col_am,   &
                                   param_log_co2_col, matrix_min, molmass_co2, &
                                   o3e_c0_high, o3e_c0_low, o3e_c1_high,       &
                                   o3e_c1_low, o3e_c2_high, o3e_c2_low,        &
                                   polyatom_max, sec_per_day, step_size,       &
                                   param_heat_rate, o3e_c3_high, o3e_c3_low,   &
                                   amm, param_co2_vmr, vmr_o
  USE interpolate_p_mod, ONLY: ip_1_lin_lin, ip_1_log_lin, interpolate_p

  ! Imported variables
  INTEGER, INTENT(IN) :: n_layer, n_profile
  REAL (RealK), INTENT(IN) :: trop_co2_mmr(:)
  REAL (RealK), INTENT(IN) :: coszin(:)
  REAL (RealK), INTENT(IN) :: p(:, :)

  ! Heating rate (passed in as LTE, returned with correction)
  REAL (RealK), INTENT(INOUT) :: heat_rate(:, :)

  ! LTE heating rates for spectral regions that don't require correction
  REAL (RealK), INTENT(IN) :: heat_rate_lte_bands(:, :)

  ! LTE heating rate for ozone hartley band (needed for O3 correction)
  REAL (RealK), INTENT(IN) :: heat_rate_hartley_band(:, :)

  ! LTE heating rates for EUV bands (photoelectron heating efficiency)
  REAL (RealK), INTENT(IN) :: heat_rate_euv_bands(:, :)

  ! Local variables

  ! Scale height and pressure corresponding to the scale height
  REAL (RealK), DIMENSION(n_level) :: x, p_x

  ! NLTE heating rate
  REAL (RealK), DIMENSION(n_profile, n_layer) :: heat_rate_nlte
  REAL (RealK), DIMENSION(n_profile, n_layer) :: heat_rate_nlte_temp
  REAL (RealK), DIMENSION(n_x_levs_min:n_x_levs_max) :: heat_rate_co2

  ! Calculation variables
  REAL (RealK) :: mean_effect_path, co2_col_am, log_co2_col_am

  ! Splining variables
  REAL (RealK), DIMENSION(n_co2_levs) :: xspline = 0.0d0, &
                                         yspline = 0.0d0, yspline2 = 0.0d0
  REAL (RealK), DIMENSION(n_x_levs) :: xspline_lev = 0.0d0, &
                                       yspline_lev = 0.0d0, yspline2_lev = 0.0d0
  REAL (RealK), DIMENSION(n_x_levs) :: heat_rate_co2_interp, p_x_interp
  REAL (RealK), DIMENSION(n_level-9) :: xspline_o3 = 0.0d0, yspline_o3 = 0.0d0,&
                                        yspline2_o3 = 0.0d0

  ! Indexing variables and misc
  INTEGER :: i, j, idx
  LOGICAL :: l_splined

  REAL (RealK) :: trop_co2_vmr, sc_co2_col_am, sc_co2_vmr

  ! Specific heats
  REAL (RealK) :: spec_heat_conv

  ! Heating efficiency in the EUV bands is set to 0.05 following
  ! Roble et al 1987 DOI:10.1029/JA092iA08p08745
  REAL (RealK), PARAMETER :: euv_eff = 0.05_RealK

  ! O3 efficiency
  REAL (RealK) :: logp, z, topval
  REAL (RealK), ALLOCATABLE, DIMENSION(:) :: o3_eff
  REAL (RealK), DIMENSION(n_level) :: o3_eff_x
  REAL (RealK), DIMENSION(n_level-9) :: p_o3_interp, o3_eff_x_interp
  REAL (RealK), DIMENSION(n_level-9) :: x_interp
  REAL (RealK), ALLOCATABLE, DIMENSION(:) :: xin

  ALLOCATE(o3_eff(n_layer))
  ALLOCATE(xin(n_layer))

  ! Create a scale height array 0 <= x <= 16.5 (top of atmosphere to bottom)
  ! in steps of delta-x = 0.25
  x = (/ (i*step_size, i = 0, n_level-1) /)

  ! Pressure corresponding to the scale height in Pascal
  DO i = 1, n_level
    p_x(i) = 100000./exp(x(i))
  END DO

  ! Loop over the profiles
  DO j = 1, n_profile

     ! Calculate the path length
     mean_effect_path = 35. / sqrt((1224. * coszin(j) * coszin(j)) + 1.)

     ! Convert tropospheric CO2 mass mixing ratio to volume mixing ratio
     trop_co2_vmr = trop_co2_mmr(j) * amm(1) / molmass_co2

     ! Loop over the dimensionless scale height
     DO i = 1, n_x_levs
        idx = i + n_x_levs_min - 1

        ! Calculating the scaled CO2 column amount above the dimensionless
        ! scale height. The scaling using the tropospheric CO2 volume mixing
        ! ratios in the model and in the parameterisation data. The mean
        ! effect path is used to calculate the effective column amount.
        sc_co2_col_am = param_co2_col_am(idx) * trop_co2_vmr / param_co2_vmr(1)
        co2_col_am = sc_co2_col_am * mean_effect_path
        log_co2_col_am = log(co2_col_am)

        ! Interpolating the parameterised heating rate from the column
        ! amount
        IF (log_co2_col_am <= param_log_co2_col(1, idx)) THEN
           heat_rate_co2(idx) = param_heat_rate(1, idx)
        ELSE IF (log_co2_col_am >= param_log_co2_col(10, idx)) THEN
           heat_rate_co2(idx) = param_heat_rate(10, idx)
        ELSE
           l_splined = .false.
           CALL interpolate_p(n_co2_levs, param_log_co2_col(:, idx), &
                              param_heat_rate(:, idx), xspline, yspline, &
                              yspline2, log_co2_col_am, heat_rate_co2(idx), &
                              IP_1_lin_lin, l_splined)
        END IF

        ! Calculating the scaled CO2 volume mixing ratio and using this
        ! to find the heating rate in K/day
        sc_co2_vmr = param_co2_vmr(idx) * trop_co2_vmr / param_co2_vmr(1)
        heat_rate_co2(idx) = heat_rate_co2(idx) * sc_co2_vmr * sec_per_day

! Fomichev used a blending function, which I'm *fairly* sure isn't needed...
!        ! Calculating a blending value to apply this correction only over
!        ! the range up to x=15.5
!        IF (x(idx) < 7.0) THEN
!           blend = (x(idx) - 3.25) / 3.75
!        ELSE IF (x(idx) < 13.75) THEN
!           blend = 1.0
!        ELSE IF (x(idx) <= 15.5) THEN
!           blend = (15.5 - x(idx)) / 1.75
!        ELSE
!           blend = 0.0
!        END IF
!        heat_rate_co2(idx) = heat_rate_co2(idx) * blend

     END DO ! End of loop over x

     ! The UM heating rates are in Q(K/day)/Cp, the heating rate calculated
     ! here is in Q(K/day)/Cv
     ! For x <= 10.25 ignore monatomic gases
     DO i = matrix_min, polyatom_max
        spec_heat_conv = 5./7.
        heat_rate_co2(i) = heat_rate_co2(i) * spec_heat_conv
     END DO
     ! x >= 10.5, consider the contribution of atomic Oxygen to the specific
     ! heat capacity
     DO i = polyatom_max+1, n_x_levs_max
        spec_heat_conv = ((5./7.) * (1.0 - vmr_o(i))) + ((3./5.) * vmr_o(i))
        heat_rate_co2(i) = heat_rate_co2(i) * spec_heat_conv
     END DO

     ! Interpolate heating rates to the original pressure levels
     ! Reverse the pressure and heating rate arrays as required for the
     ! interpolation routine
     p_x_interp = p_x(n_x_levs_max:n_x_levs_min:-1)
     heat_rate_co2_interp = heat_rate_co2(n_x_levs_max:n_x_levs_min:-1)

     ! Interpolate to match the pressure levels of the model
     DO i = 1, n_layer
        l_splined = .false.
! There isn't a noticeable difference between IP_1_log_lin and IP_3_log_lin
        CALL interpolate_p(n_x_levs, p_x_interp, &
                           heat_rate_co2_interp, xspline_lev, yspline_lev, &
                           yspline2_lev, p(j, i), heat_rate_nlte(j, i), &
                           IP_1_log_lin, l_splined)
     END DO

     ! Calculating the O3 heating efficiency
     topval = 0.0
     DO i = 1, n_level
        logp = log10(p_x(i))
        IF (logp > 2) THEN
           o3_eff_x(i) = 1.0
        ELSE IF ((logp <= 2) .AND. (logp > 0)) THEN
           z = logp - 1.0
           o3_eff_x(i) = o3e_c0_high + (o3e_c1_high * z) &
             + (o3e_c2_high * z * z) + (o3e_c3_high * z * z * z)
        ELSE IF ((logp <= 0) .AND. (logp > -2)) THEN
           z = logp + 1.0
           o3_eff_x(i) = o3e_c0_low + (o3e_c1_low * z) &
             + (o3e_c2_low * z * z) + (o3e_c3_low * z * z * z)
           ! Finding the topmost value of the function to extend it - this
           ! is the last value to be set in the loop
           topval = o3_eff_x(i)
        ELSE
           o3_eff_x(i) = topval
        END IF
     END DO

     ! Interpolate over the discontinuity in the Hartley O3 efficiency
     ! function from 6 <= x <= 8
     p_o3_interp(1:24) = p_x(1:24)
     p_o3_interp(25:n_level-9) = p_x(34:n_level)
     o3_eff_x_interp(1:24) = o3_eff_x(1:24)
     o3_eff_x_interp(25:n_level-9) = o3_eff_x(34:n_level)

     x_interp(1:24) = x(1:24)
     x_interp(25:n_level-9) = x(34:n_level)
     DO i = 1, n_layer
        xin(i) = log(100000./p(j, i))
     END DO
     DO i = 1, n_layer
        l_splined = .false.
        CALL interpolate_p((n_level-9), x_interp, o3_eff_x_interp, &
                           xspline_o3, yspline_o3, yspline2_o3, xin(i), &
                           o3_eff(i), IP_1_lin_lin, l_splined)
     END DO

     ! The Hartley Band requires an NLTE correction.
     ! The Huggins/Chappuis O3 bands do not require an NLTE correction.
     ! The EUV bands require a correction for the heating efficiency of
     ! photoelectrons.
     DO i = 1, n_layer
        heat_rate_nlte_temp(j, i) = heat_rate_hartley_band(j, i) * o3_eff(i) &
          + heat_rate_euv_bands(j, i) * euv_eff &
          + heat_rate_lte_bands(j, i)
     END DO

     ! Adding the O3 to the CO2 heating
     DO i = 1, n_layer
        heat_rate_nlte(j, i) = heat_rate_nlte_temp(j, i) &
          + MAX(heat_rate_nlte(j, i), 0.0_RealK)
     END DO

     ! Blend the heating rates with the LTE rates (use the same level as LW)
     DO i = 1, n_layer
        IF (p(j,i) < 10.0) THEN
           heat_rate(j, i) = heat_rate_nlte(j, i)
        END IF
     END DO

  END DO ! End of loop over n_profile

  DEALLOCATE(xin)
  DEALLOCATE(o3_eff)

END SUBROUTINE nlte_heating_sw

END MODULE sw_nlte_heating_mod
