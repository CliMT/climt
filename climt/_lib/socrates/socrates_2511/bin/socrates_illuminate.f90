! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Calculate the parameters for external illumination of the atmosphere

module socrates_illuminate

use def_orbit, only: &
  ip_elements_user, &
  ip_elements_earth_fixed, &
  ip_elements_earth_secular_variation, &
  ip_spin_user, &
  ip_spin_earth_day, &
  ip_spin_fixed_sun

implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_ILLUMINATE'
contains

subroutine illuminate( &
    l_calendar_360, l_observer, l_stellar_position, l_stellar_angle, &
    l_slope, l_shading, &
    n_profile, n_horiz_layer, n_horiz_ang, i_elements, i_spin, &
    year, month, day_of_month, day_of_year, second_of_day, length_of_timestep, &
    epoch, eccentricity, eccentricity_inc, arg_periapsis, arg_periapsis_inc, &
    obliquity, obliquity_inc, semimajor_axis, semimajor_axis_inc, &
    mean_anomaly, mean_anomaly_inc, hour_angle, hour_angle_inc, &
    fixed_zenith_angle, fixed_azimuth_angle, observer_lon, observer_lat, &
    latitude, longitude, stellar_constant, &
    slope_aspect, slope_angle, horizon_aspect, horizon_angle, &
    sin_stellar_declination, stellar_eqn_of_time, stellar_constant_scaling, &
    sin_observer_declination, observer_eqn_of_time, observed_orbital_phase, &
    cos_zenith_angle, lit_fraction, stellar_irradiance, &
    observer_cos_zen_ang, observed_fraction, orographic_correction)

  use realtype_rd, only: RealK, RealExt
  use def_orbit, only: StrOrbit, set_orbit
  use solpos_mod, only: solpos
  use solang_mod, only: solang
  use solinc_mod, only: solinc

  implicit none

  logical, intent(in), optional :: l_calendar_360, l_observer
  logical, intent(in), optional :: l_stellar_position, l_stellar_angle
  logical, intent(in), optional :: l_slope, l_shading
  integer, intent(in), optional :: n_profile, n_horiz_layer, n_horiz_ang
  integer, intent(in), optional :: i_elements, i_spin
  integer, intent(in), optional :: year, month, day_of_month, day_of_year
  real(RealExt), intent(in), optional :: second_of_day, length_of_timestep
  real(RealExt), intent(in), optional :: &
    epoch, eccentricity, eccentricity_inc, arg_periapsis, arg_periapsis_inc, &
    obliquity, obliquity_inc, semimajor_axis, semimajor_axis_inc, &
    mean_anomaly, mean_anomaly_inc, hour_angle, hour_angle_inc, &
    fixed_zenith_angle, fixed_azimuth_angle, observer_lon, observer_lat
  real(RealExt), intent(in), optional :: latitude(:)
  real(RealExt), intent(in), optional :: longitude(:)
  real(RealExt), intent(in), optional :: stellar_constant
  real(RealExt), intent(in), optional :: slope_aspect(:), slope_angle(:)
  real(RealExt), intent(in), optional :: horizon_aspect(:), horizon_angle(:)

  real(RealExt), intent(inout), optional :: sin_stellar_declination
  real(RealExt), intent(inout), optional :: stellar_eqn_of_time
  real(RealExt), intent(inout), optional :: stellar_constant_scaling
  real(RealExt), intent(inout), optional :: sin_observer_declination
  real(RealExt), intent(inout), optional :: observer_eqn_of_time
  real(RealExt), intent(inout), optional :: observed_orbital_phase

  real(RealExt), intent(out), optional :: cos_zenith_angle(:)
  real(RealExt), intent(out), optional :: lit_fraction(:)
  real(RealExt), intent(out), optional :: stellar_irradiance(:)
  real(RealExt), intent(out), optional :: observer_cos_zen_ang(:)
  real(RealExt), intent(out), optional :: observed_fraction(:)
  real(RealExt), intent(out), optional :: orographic_correction(:)

  ! Local variables
  integer :: i, k, kk, l, ll
  type(StrOrbit) :: orbit
  real(RealK) :: eqt, sindec, scs, sindec_obs, eqt_obs, phase_obs
  real(RealK) :: second, dt
  real(RealK), allocatable :: cos_zen(:), lit_frac(:), lat(:), lon(:)
  real(RealK), allocatable :: sol_azimuth(:), cosz_beg(:), cosz_end(:)
  real(RealK), allocatable :: slope_asp(:), slope_ang(:), orog_corr(:)
  real(RealK), allocatable :: horiz_asp(:, :), horiz_ang(:, :, :)
  real(RealK), allocatable :: obs_cos_zen(:), obs_frac(:)


  call set_orbit(orbit, &
    i_elements          = i_elements, &
    i_spin              = i_spin, &
    epoch               = epoch, &
    eccentricity        = eccentricity, &
    eccentricity_inc    = eccentricity_inc, &
    arg_periapsis       = arg_periapsis, &
    arg_periapsis_inc   = arg_periapsis_inc, &
    obliquity           = obliquity, &
    obliquity_inc       = obliquity_inc, &
    semimajor_axis      = semimajor_axis, &
    semimajor_axis_inc  = semimajor_axis_inc, &
    mean_anomaly        = mean_anomaly, &
    mean_anomaly_inc    = mean_anomaly_inc, &
    hour_angle          = hour_angle, &
    hour_angle_inc      = hour_angle_inc, &
    fixed_zenith_angle  = fixed_zenith_angle, &
    fixed_azimuth_angle = fixed_azimuth_angle, &
    observer_lon        = observer_lon, &
    observer_lat        = observer_lat )


  if (present(second_of_day)) second = real(second_of_day, RealK)
  if (present(length_of_timestep)) dt = real(length_of_timestep, RealK)


  if (present(l_stellar_position)) then
  if (l_stellar_position) then
    ! Find stellar position from the planet surface
    call solpos(orbit, year, day_of_year, second, dt, &
      sindec, eqt, scs, l_observer, l_calendar_360, &
      sindec_obs, eqt_obs, phase_obs)

    if (present(sin_stellar_declination)) &
      sin_stellar_declination  = real(sindec, RealExt)
    if (present(stellar_eqn_of_time)) &
      stellar_eqn_of_time      = real(eqt, RealExt)
    if (present(stellar_constant_scaling)) &
      stellar_constant_scaling = real(scs, RealExt)
    if (present(sin_observer_declination)) &
      sin_observer_declination = real(sindec_obs, RealExt)
    if (present(observer_eqn_of_time)) &
      observer_eqn_of_time     = real(eqt_obs, RealExt)
    if (present(observed_orbital_phase)) &
      observed_orbital_phase   = real(phase_obs, RealExt)
  end if
  end if


  if (present(l_stellar_angle)) then
  if (l_stellar_angle) then
    ! Find stellar zenith angle and lit fraction of timestep
    if (present(sin_stellar_declination)) &
      sindec = real(sin_stellar_declination, RealK)
    if (present(stellar_eqn_of_time)) &
      eqt    = real(stellar_eqn_of_time, RealK)
    if (present(stellar_constant_scaling)) &
      scs    = real(stellar_constant_scaling, RealK)

    allocate(lat(n_profile), lon(n_profile))
    lat = real(latitude(1:n_profile), RealK)
    lon = real(longitude(1:n_profile), RealK)
    allocate(cos_zen(n_profile), lit_frac(n_profile))
    allocate(sol_azimuth(n_profile), cosz_beg(n_profile), cosz_end(n_profile))

    call solang(orbit, second, dt, sindec, eqt, lat, lon, n_profile, &
      lit_frac, cos_zen, sol_azimuth, cosz_beg, cosz_end)

    if (present(lit_fraction)) &
      lit_fraction(1:n_profile) = real(lit_frac, RealExt)
    if (present(cos_zenith_angle)) &
      cos_zenith_angle(1:n_profile) = real(cos_zen, RealExt)
    if (present(stellar_irradiance).and.present(stellar_constant)) &
      stellar_irradiance(1:n_profile) = stellar_constant &
        * real(scs, RealExt) * real(lit_frac, RealExt)

    deallocate(lit_frac, lon, lat)

    if (present(l_slope) .and. present(l_shading) .and. &
        present(orographic_correction)) then
    if (l_slope .or. l_shading) then
      allocate(orog_corr(n_profile))
      allocate(slope_ang(n_profile), slope_asp(n_profile))
      slope_ang = real(slope_angle(1:n_profile), RealK)
      slope_asp = real(slope_aspect(1:n_profile), RealK)
      if (l_shading) then
        allocate(horiz_ang(n_horiz_layer, n_horiz_ang, n_profile))
        allocate(horiz_asp(n_horiz_ang, n_profile))
        ll=0
        kk=0
        do l=1, n_profile
          do k=1, n_horiz_ang
            kk = kk + 1
            horiz_asp(k, l) = real(horizon_aspect(kk), RealK)
            do i=1, n_horiz_layer
              ll = ll + 1
              horiz_ang(i, k, l) = real(horizon_angle(ll), RealK)
            end do
          end do
        end do
        call solinc(n_profile, cos_zen, sol_azimuth, &
          slope_asp, slope_ang, orog_corr, l_shading, &
          n_horiz_ang, horiz_asp, horiz_ang, cosz_beg, cosz_end)
        deallocate(horiz_ang, horiz_asp)
      else
        call solinc(n_profile, cos_zen, sol_azimuth, &
          slope_asp, slope_ang, orog_corr)
      end if
      orographic_correction(1:n_profile) = real(orog_corr, RealExt)
      deallocate(slope_ang, slope_asp, orog_corr)
    else
      orographic_correction(1:n_profile) = 1.0_RealExt
    end if
    end if

    deallocate(cos_zen, sol_azimuth, cosz_beg, cosz_end)
  end if
  end if


  if (present(l_observer)) then
  if (l_observer) then
    ! Find observer zenith angle and observed fraction of timestep
    if (present(sin_observer_declination)) &
      sindec_obs = real(sin_observer_declination, RealK)
    if (present(observer_eqn_of_time)) &
      eqt_obs    = real(observer_eqn_of_time, RealK)

    allocate(lat(n_profile), lon(n_profile))
    lat = real(latitude(1:n_profile), RealK)
    lon = real(longitude(1:n_profile), RealK)
    allocate(obs_cos_zen(n_profile), obs_frac(n_profile))

    call solang(orbit, second, dt, &
      sindec_obs, eqt_obs, lat, lon, n_profile, &
      obs_frac, obs_cos_zen)

    if (present(observed_fraction)) &
      observed_fraction(1:n_profile) = real(obs_frac, RealExt)
    if (present(observer_cos_zen_ang)) &
      observer_cos_zen_ang(1:n_profile) = real(obs_cos_zen, RealExt)

    deallocate(obs_frac, obs_cos_zen, lon, lat)
  end if
  end if

end subroutine illuminate
end module socrates_illuminate
