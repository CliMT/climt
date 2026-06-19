! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the variables in the Socrates boundary conditions type
!
!------------------------------------------------------------------------------
module socrates_set_bound
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_BOUND'
contains

subroutine set_bound(bound, control, dimen, spectrum, &
  n_profile, profile_list, n_tile, &
  t_ground, flux_ground, flux_ground_1d, &
  cos_zenith_angle, solar_irrad, orog_corr, &
  l_grey_albedo, grey_albedo, albedo_diff, albedo_dir, &
  albedo_diff_1d, albedo_dir_1d, &
  frac_tile, t_tile, flux_tile, albedo_diff_tile, albedo_dir_tile, &
  frac_tile_1d, t_tile_1d, flux_tile_1d, albedo_diff_tile_1d, albedo_dir_tile_1d, &
  l_profile_last, l_debug, i_profile_debug)

use def_bound,    only: StrBound, allocate_bound
use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use realtype_rd,  only: RealK, RealExt
use rad_pcf,      only: &
  ip_solar, ip_infra_red, ip_surf_alb_diff, ip_surf_alb_dir

implicit none


! Boundary properties:
type(StrBound),    intent(out) :: bound

! Control options:
type(StrCtrl),     intent(in)  :: control

! Dimensions:
type(StrDim),      intent(in)  :: dimen

! Spectral data:
type(StrSpecData), intent(in)  :: spectrum

integer, intent(in) :: n_profile
!   Number of atmospheric profiles for radiation calculations
integer, intent(in), optional :: profile_list(:)
!   List of profiles to use from input fields
integer, intent(in), optional :: n_tile
!   Number of surface tiles for radiation calculations

real(RealExt), intent(in), optional :: t_ground(:)
!   Effective radiative temperature over whole grid-box
real(RealExt), intent(in), dimension(:,:), optional :: &
  flux_ground
!   Effective surface emission over whole grid-box
!   (n_profile,n_band)
real(RealExt), intent(in), dimension(:), optional :: &
  flux_ground_1d
!   Effective surface emission over whole grid-box
!   (n_profile*n_band)

real(RealExt), intent(in), optional :: cos_zenith_angle(:)
!   Cosine of solar zenith angle
real(RealExt), intent(in), optional :: solar_irrad(:)
!   Solar irradiance at top-of-atmosphere (mean over timestep)
real(RealExt), intent(in), optional :: orog_corr(:)
!   Orographic correction factor

logical, intent(in), optional :: l_grey_albedo
!   Set a single grey albedo for the surface
real(RealExt), intent(in), optional :: grey_albedo
!   Grey surface albedo

real(RealExt), intent(in), dimension(:, :), optional :: &
  albedo_diff, albedo_dir
!   Diffuse, direct albedo (n_profile, n_band)

real(RealExt), intent(in), dimension(:), optional :: &
  albedo_diff_1d, albedo_dir_1d
!   1d diffuse, direct albedo (n_band)

real(RealExt), intent(in), dimension(:, :), optional :: &
  frac_tile, t_tile
!   Tile fraction, temperature (n_profile, n_tile)

real(RealExt), intent(in), dimension(:, :, :), optional :: &
  flux_tile
!   Tile emission (n_profile, n_tile, n_band)

real(RealExt), intent(in), dimension(:, :, :), optional :: &
  albedo_diff_tile, albedo_dir_tile
!   Tile albedos (n_profile, n_tile, n_band)

real(RealExt), intent(in), dimension(:), optional :: &
  frac_tile_1d, t_tile_1d
!   1d tile fraction, temperature (n_tile)

real(RealExt), intent(in), dimension(:), optional :: &
  flux_tile_1d
!   Tile emission (n_profile*n_tile*n_band)

real(RealExt), intent(in), dimension(:), optional :: &
  albedo_diff_tile_1d, albedo_dir_tile_1d
!   1d tile albedos (n_tile*n_band)

logical, intent(in), optional :: l_profile_last
!   Loop over profiles is last in input fields

logical, intent(in), optional :: l_debug
integer, intent(in), optional :: i_profile_debug
!   Options for outputting debugging information

! Local variables.
integer :: l, ll, i_band, i_tile, list(n_profile)
logical :: l_grey_alb, l_last

! Allocate structure for the core radiation code interface
call allocate_bound(bound, dimen, spectrum)

if (present(profile_list)) then
  list = profile_list(1:n_profile)
else
  do l=1, n_profile
    list(l) = l
  end do
end if

if (present(l_profile_last)) then
  l_last = l_profile_last
else
  l_last = .false.
end if

if (present(l_grey_albedo)) then
  l_grey_alb = l_grey_albedo
else
  l_grey_alb = .false.
end if

! Surface albedo
if (l_grey_alb .and. present(grey_albedo)) then
  do i_band=1, spectrum%basic%n_band
    do l=1, n_profile
      bound%rho_alb(l, ip_surf_alb_diff, i_band) = real(grey_albedo, RealK)
    end do
    do l=1, n_profile
      bound%rho_alb(l, ip_surf_alb_dir,  i_band) = real(grey_albedo, RealK)
    end do
  end do
else
  if (present(albedo_diff)) then
    if (l_last) then
      do i_band=1, spectrum%basic%n_band
        do l=1, n_profile
          bound%rho_alb(l, ip_surf_alb_diff, i_band) &
            = real(albedo_diff(i_band, list(l)), RealK)
        end do
      end do
    else
      do i_band=1, spectrum%basic%n_band
        do l=1, n_profile
          bound%rho_alb(l, ip_surf_alb_diff, i_band) &
            = real(albedo_diff(list(l), i_band), RealK)
        end do
      end do
    end if
  else if (present(albedo_diff_1d)) then
    if (l_last) then
      do i_band=1, spectrum%basic%n_band
        do l=1, n_profile
          ll = spectrum%basic%n_band*(list(l)-1) + i_band
          bound%rho_alb(l, ip_surf_alb_diff, i_band) &
            = real(albedo_diff_1d(ll), RealK)
        end do
      end do
    else
      do i_band=1, spectrum%basic%n_band
        do l=1, n_profile
          ll = n_profile*(i_band-1) + list(l)
          bound%rho_alb(l, ip_surf_alb_diff, i_band) &
            = real(albedo_diff_1d(ll), RealK)
        end do
      end do
    end if
  else
    do i_band=1, spectrum%basic%n_band
      do l=1, n_profile
        bound%rho_alb(l, ip_surf_alb_diff, i_band) = 0.0_RealK
      end do
    end do
  end if
  if (present(albedo_dir)) then
    if (l_last) then
      do i_band=1, spectrum%basic%n_band
        do l=1, n_profile
          bound%rho_alb(l, ip_surf_alb_dir, i_band) &
            = real(albedo_dir(i_band, list(l)), RealK)
        end do
      end do
    else
      do i_band=1, spectrum%basic%n_band
        do l=1, n_profile
          bound%rho_alb(l, ip_surf_alb_dir, i_band) &
            = real(albedo_dir(list(l), i_band), RealK)
        end do
      end do
    end if
  else if (present(albedo_dir_1d)) then
    if (l_last) then
      do i_band=1, spectrum%basic%n_band
        do l=1, n_profile
          ll = spectrum%basic%n_band*(list(l)-1) + i_band
          bound%rho_alb(l, ip_surf_alb_dir, i_band) &
            = real(albedo_dir_1d(ll), RealK)
        end do
      end do
    else
      do i_band=1, spectrum%basic%n_band
        do l=1, n_profile
          ll = n_profile*(i_band-1) + list(l)
          bound%rho_alb(l, ip_surf_alb_dir, i_band) &
            = real(albedo_dir_1d(ll), RealK)
        end do
      end do
    end if
  else
    do i_band=1, spectrum%basic%n_band
      do l=1, n_profile
        bound%rho_alb(l, ip_surf_alb_dir, i_band) = 0.0_RealK
      end do
    end do
  end if
end if

! Surface temperature
if (present(t_ground)) then
  do l=1, n_profile
    bound%t_ground(l) = real(t_ground(list(l)), RealK)
  end do
else
  do l=1, n_profile
    bound%t_ground(l) = 0.0_RealK
  end do
end if

! Surface emission
if (control%l_flux_ground) then
  if (present(flux_ground)) then
    if (l_last) then
      do i_band=1, spectrum%basic%n_band
        do l=1, n_profile
          bound%flux_ground(l, i_band) &
            = real(flux_ground(i_band, list(l)), RealK)
        end do
      end do
    else
      do i_band=1, spectrum%basic%n_band
        do l=1, n_profile
          bound%flux_ground(l, i_band) &
            = real(flux_ground(list(l), i_band), RealK)
        end do
      end do
    end if
  else if (present(flux_ground_1d)) then
    if (l_last) then
      do i_band=1, spectrum%basic%n_band
        do l=1, n_profile
          ll = spectrum%basic%n_band*(list(l)-1) + i_band
          bound%flux_ground(l, i_band) &
            = real(flux_ground_1d(ll), RealK)
        end do
      end do
    else
      do i_band=1, spectrum%basic%n_band
        do l=1, n_profile
          ll = n_profile*(i_band-1) + list(l)
          bound%flux_ground(l, i_band) &
            = real(flux_ground_1d(ll), RealK)
        end do
      end do
    end if
  else
    do i_band=1, spectrum%basic%n_band
      do l=1, n_profile
        bound%flux_ground(l, i_band) = 0.0_RealK
      end do
    end do
  end if
end if

bound%n_point_tile=0
bound%n_tile=1
if (control%l_tile .and. present(n_tile) .and. &
    (present(frac_tile) .or. present(frac_tile_1d))) then

  ! Set up the surface tiling variables
  ! Treat all points as tiled when l_tile is true
  bound%n_tile=n_tile
  bound%n_point_tile = n_profile
  do l=1, n_profile
    bound%list_tile(l) = l
  end do

  if (present(frac_tile)) then
    if (l_last) then
      do i_tile=1, n_tile
        do l=1, n_profile
          bound%frac_tile(l, i_tile) = real(frac_tile(i_tile, list(l)), RealK)
        end do
      end do
    else
      do i_tile=1, n_tile
        do l=1, n_profile
          bound%frac_tile(l, i_tile) = real(frac_tile(list(l), i_tile), RealK)
        end do
      end do
    end if
  else
    if (l_last) then
      do i_tile=1, n_tile
        do l=1, n_profile
          ll = n_tile*(list(l)-1) + i_tile
          bound%frac_tile(l, i_tile) = real(frac_tile_1d(ll), RealK)
        end do
      end do
    else
      do i_tile=1, n_tile
        do l=1, n_profile
          ll = n_profile*(i_tile-1) + list(l)
          bound%frac_tile(l, i_tile) = real(frac_tile_1d(ll), RealK)
        end do
      end do
    end if
  end if

  ! Diffuse tile albedos
  if ((present(albedo_diff_tile) .or. present(albedo_diff_tile_1d)) &
      .and. .not. l_grey_alb) then
    if (present(albedo_diff_tile)) then
      if (l_last) then
        do i_band=1, spectrum%basic%n_band
          do i_tile=1, n_tile
            do l=1, n_profile
              bound%rho_alb_tile(l, ip_surf_alb_diff, i_tile, i_band) &
                = real(albedo_diff_tile(i_tile, i_band, list(l)), RealK)
            end do
          end do
        end do
      else
        do i_band=1, spectrum%basic%n_band
          do i_tile=1, n_tile
            do l=1, n_profile
              bound%rho_alb_tile(l, ip_surf_alb_diff, i_tile, i_band) &
                = real(albedo_diff_tile(list(l), i_tile, i_band), RealK)
            end do
          end do
        end do
      end if
    else
      if (l_last) then
        do i_band=1, spectrum%basic%n_band
          do i_tile=1, n_tile
            do l=1, n_profile
              ll = spectrum%basic%n_band*n_tile*(list(l)-1) &
                 + n_tile*(i_band-1) + i_tile
              bound%rho_alb_tile(l, ip_surf_alb_diff, i_tile, i_band) &
                = real(albedo_diff_tile_1d(ll), RealK)
            end do
          end do
        end do
      else
        do i_band=1, spectrum%basic%n_band
          do i_tile=1, n_tile
            do l=1, n_profile
              ll = n_profile*n_tile*(i_band-1) &
                 + n_profile*(i_tile-1) + list(l)
              bound%rho_alb_tile(l, ip_surf_alb_diff, i_tile, i_band) &
                = real(albedo_diff_tile_1d(ll), RealK)
            end do
          end do
        end do
      end if
    end if
    ! Ensure the total albedo is consistent with the tile albedos
    do i_band=1, spectrum%basic%n_band
      do l=1, n_profile
        bound%rho_alb(l, ip_surf_alb_diff, i_band) = 0.0_RealK
      end do
      do i_tile=1, n_tile
        do l=1, n_profile
          bound%rho_alb(l, ip_surf_alb_diff, i_band) &
            = bound%rho_alb(l, ip_surf_alb_diff, i_band) &
            + bound%rho_alb_tile(l, ip_surf_alb_diff, i_tile, i_band) &
            * bound%frac_tile(l, i_tile)
        end do
      end do
    end do
  else
    ! When not present just use the gridbox mean diffuse albedo 
    do i_band=1, spectrum%basic%n_band
      do i_tile=1, n_tile
        do l=1, n_profile
          bound%rho_alb_tile(l, ip_surf_alb_diff, i_tile, i_band) &
            = bound%rho_alb(l, ip_surf_alb_diff, i_band)
        end do
      end do
    end do
  end if

  ! Direct tile albedos
  if ((present(albedo_dir_tile) .or. present(albedo_dir_tile_1d)) &
      .and. .not. l_grey_alb) then
    if (present(albedo_dir_tile)) then
      if (l_last) then
        do i_band=1, spectrum%basic%n_band
          do i_tile=1, n_tile
            do l=1, n_profile
              bound%rho_alb_tile(l, ip_surf_alb_dir, i_tile, i_band) &
                = real(albedo_dir_tile(i_tile, i_band, list(l)), RealK)
            end do
          end do
        end do
      else
        do i_band=1, spectrum%basic%n_band
          do i_tile=1, n_tile
            do l=1, n_profile
              bound%rho_alb_tile(l, ip_surf_alb_dir, i_tile, i_band) &
                = real(albedo_dir_tile(list(l), i_tile, i_band), RealK)
            end do
          end do
        end do
      end if
    else
      if (l_last) then
        do i_band=1, spectrum%basic%n_band
          do i_tile=1, n_tile
            do l=1, n_profile
              ll = spectrum%basic%n_band*n_tile*(list(l)-1) &
                 + n_tile*(i_band-1) + i_tile
              bound%rho_alb_tile(l, ip_surf_alb_dir, i_tile, i_band) &
                = real(albedo_dir_tile_1d(ll), RealK)
            end do
          end do
        end do
      else
        do i_band=1, spectrum%basic%n_band
          do i_tile=1, n_tile
            do l=1, n_profile
              ll = n_profile*n_tile*(i_band-1) &
                 + n_profile*(i_tile-1) + list(l)
              bound%rho_alb_tile(l, ip_surf_alb_dir, i_tile, i_band) &
                = real(albedo_dir_tile_1d(ll), RealK)
            end do
          end do
        end do
      end if
    end if
    ! Ensure the total albedo is consistent with the tile albedos
    do i_band=1, spectrum%basic%n_band
      do l=1, n_profile
        bound%rho_alb(l, ip_surf_alb_dir, i_band) = 0.0_RealK
      end do
      do i_tile=1, n_tile
        do l=1, n_profile
          bound%rho_alb(l, ip_surf_alb_dir, i_band) &
            = bound%rho_alb(l, ip_surf_alb_dir, i_band) &
            + bound%rho_alb_tile(l, ip_surf_alb_dir, i_tile, i_band) &
            * bound%frac_tile(l, i_tile)
        end do
      end do
    end do
  else
    ! When not present just use the gridbox mean direct albedo 
    do i_band=1, spectrum%basic%n_band
      do i_tile=1, n_tile
        do l=1, n_profile
          bound%rho_alb_tile(l, ip_surf_alb_dir, i_tile, i_band) &
            = bound%rho_alb(l, ip_surf_alb_dir, i_band)
        end do
      end do
    end do
  end if

  if (present(t_tile)) then
    ! Set the tile temperatures (t_ground will not be used on these points)
    if (l_last) then
      do i_tile=1, n_tile
        do l=1, n_profile
          bound%t_tile(l, i_tile) = real(t_tile(i_tile, list(l)), RealK)
        end do
      end do
    else
      do i_tile=1, n_tile
        do l=1, n_profile
          bound%t_tile(l, i_tile) = real(t_tile(list(l), i_tile), RealK)
        end do
      end do
    end if
  else if (present(t_tile_1d)) then
    ! Set the tile temperatures from 1d input
    if (l_last) then
      do i_tile=1, n_tile
        do l=1, n_profile
          ll = n_tile*(list(l)-1) + i_tile
          bound%t_tile(l, i_tile) = real(t_tile_1d(ll), RealK)
        end do
      end do
    else
      do i_tile=1, n_tile
        do l=1, n_profile
          ll = n_profile*(i_tile-1) + list(l)
          bound%t_tile(l, i_tile) = real(t_tile_1d(ll), RealK)
        end do
      end do
    end if
  else
    ! When not present just use the gridbox mean surface temperature
    do i_tile=1, n_tile
      do l=1, n_profile
        bound%t_tile(l, i_tile) = bound%t_ground(l)
      end do
    end do
  end if

  if (present(flux_tile)) then
    ! Set the tile fluxes (flux_ground will not be used on these points)
    if (l_last) then
      do i_band=1, spectrum%basic%n_band
        do i_tile=1, n_tile
          if (control%l_flux_tile(i_tile)) then
            do l=1, n_profile
              bound%flux_tile(l, i_tile, i_band) &
                = real(flux_tile(i_tile, i_band, list(l)), RealK)
            end do
          end if
        end do
      end do
    else
      do i_band=1, spectrum%basic%n_band
        do i_tile=1, n_tile
          if (control%l_flux_tile(i_tile)) then
            do l=1, n_profile
              bound%flux_tile(l, i_tile, i_band) &
                = real(flux_tile(list(l), i_tile, i_band), RealK)
            end do
          end if
        end do
      end do
    end if
  else if (present(flux_tile_1d)) then
    if (l_last) then
      do i_band=1, spectrum%basic%n_band
        do i_tile=1, n_tile
          if (control%l_flux_tile(i_tile)) then
            do l=1, n_profile
              ll = spectrum%basic%n_band*n_tile*(list(l)-1) &
                 + n_tile*(i_band-1) + i_tile
              bound%flux_tile(l, i_tile, i_band) &
                = real(flux_tile_1d(ll), RealK)
          end do
          end if
        end do
      end do
    else
      do i_band=1, spectrum%basic%n_band
        do i_tile=1, n_tile
          if (control%l_flux_tile(i_tile)) then
            do l=1, n_profile
              ll = n_profile*n_tile*(i_band-1) &
                 + n_profile*(i_tile-1) + list(l)
              bound%flux_tile(l, i_tile, i_band) &
                = real(flux_tile_1d(ll), RealK)
            end do
          end if
        end do
      end do
    end if
  else
    ! When not present set to zero
    do i_band=1, spectrum%basic%n_band
      do i_tile=1, n_tile
        do l=1, n_profile
          bound%flux_tile(l, i_tile, i_band) = 0.0_RealK
        end do
      end do
    end do
  end if

end if


! Set the surface basis functions for a Lambertian surface.
bound%n_brdf_basis_fnc=1
! By defining F_{1,0,0,0} to be 4, rho_alb becomes equal to the diffuse albedo.
bound%f_brdf(1, 0, 0, 0)=4.0_RealK


! Orographic correction factor
if (present(orog_corr)) then
  do l=1, n_profile
    bound%orog_corr(l) = max(real(orog_corr(list(l)),RealK), epsilon(1.0_RealK))
  end do
end if


! Incident solar flux
if (present(cos_zenith_angle) .and. present(solar_irrad)) then
  do l=1, n_profile
    if (cos_zenith_angle(list(l)) > 0.0_RealExt) then
      bound%solar_irrad(l) = real(solar_irrad(list(l)), RealK)
      bound%zen_0(l)=1.0_RealK/real(cos_zenith_angle(list(l)), RealK)
    else
      bound%solar_irrad(l)=0.0_RealK
      bound%zen_0(l)=1.0_RealK
    end if
  end do
end if


if (present(l_debug)) then
  if (l_debug) then
    if (present(i_profile_debug)) then
      l = i_profile_debug
    else
      l = 1
    end if
    write(2000+l,'(A)') 'TEMP(K) SOLAR_IRRAD(WM-2) ZEN_0 OROG_CORR'
    write(2000+l,'(4(1pe16.8))') bound%t_ground(l), &
      bound%solar_irrad(l), bound%zen_0(l), bound%orog_corr(l)
    write(2000+l,'(A)') 'BAND DIFFUSE_ALBEDO DIRECT_ALBEDO'
    do i_band=1, spectrum%basic%n_band
      write(2000+l,'(i0, 2(1pe16.8))') i_band, &
        bound%rho_alb(l, ip_surf_alb_diff, i_band), &
        bound%rho_alb(l, ip_surf_alb_dir, i_band)
    end do
    if (control%l_tile .and. present(n_tile)) then
!      ll = findloc(bound%list_tile, l, 1)
      do ll=1, bound%n_point_tile
        if (bound%list_tile(ll) == l) then
          write(2000+l,'(A)') 'TILE FRAC TEMP(K)'
          do i_tile=1, n_tile
            write(2000+l,'(i0, 2(1pe16.8))') i_tile, &
              bound%frac_tile(ll, i_tile), bound%t_tile(ll, i_tile)
          end do
        end if
      end do
    end if
  end if
end if


end subroutine set_bound
end module socrates_set_bound
