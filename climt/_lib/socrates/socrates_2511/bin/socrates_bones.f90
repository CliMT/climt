! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Update radiative fluxes using simple correction schemes

module socrates_bones

implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_BONES'
contains

! The bare "bones" of radiative transfer calculations. In ancient China,
! the future weather was predicted by reading the cracks that appeared
! when oracle bones were heated. Not particularly accurate, but simple
! and cheap.

subroutine bones(n_profile, n_layer, n_layer_stride, n_tile, &
  l_cos_zen_correction, cos_zen_rts, lit_frac_rts, cos_zen_mts, lit_frac_mts, &
  l_trans_zen_correction, l_orog_corr_rts, orog_corr_rts, &
  l_grey_emis_correction, grey_albedo_tile, frac_tile, t_tile, &
  grey_albedo_tile_1d, frac_tile_1d, t_tile_1d, &
  l_profile_last, l_debug, i_profile_debug, &
  heating_rate_rts, flux_up_tile_rts, flux_up_blue_tile_rts, &
  flux_direct_toa_rts, flux_up_toa_rts, &
  flux_direct_surf_rts, flux_down_surf_rts, flux_up_surf_rts, &
  flux_direct_blue_surf_rts, flux_down_blue_surf_rts, &
  heating_rate_1d_rts, flux_up_tile_1d_rts, flux_up_blue_tile_1d_rts, &
  heating_rate_mts, flux_up_tile_mts, flux_up_blue_tile_mts, &
  flux_direct_toa_mts, flux_up_toa_mts, &
  flux_direct_surf_mts, flux_down_surf_mts, flux_up_surf_mts, &
  flux_direct_blue_surf_mts, flux_down_blue_surf_mts, &
  heating_rate_1d_mts, flux_up_tile_1d_mts, flux_up_blue_tile_1d_mts)

use realtype_rd, only: RealExt
use rad_ccf, only: stefan_boltzmann

implicit none

integer, intent(in) :: n_profile
!   Number of columns to operate on
integer, intent(in) :: n_layer
!   Number of layers for radiation
integer, intent(in), optional :: n_layer_stride
!   Number of layers in 1d arrays
integer, intent(in), optional :: n_tile
!   Number of surface tiles

logical, intent(in), optional :: l_cos_zen_correction
!   Apply simple solar zenith angle correction
real(RealExt), intent(in), optional :: cos_zen_rts(n_profile)
!   Mean cosine of solar zenith angle over lit fraction of radiation timestep
real(RealExt), intent(in), optional :: lit_frac_rts(n_profile)
!   Lit fraction of radiation timestep
real(RealExt), intent(in), optional :: cos_zen_mts(n_profile)
!   Mean cosine of solar zenith angle over lit fraction of model timestep
real(RealExt), intent(in), optional :: lit_frac_mts(n_profile)
!   Lit fraction of model timestep

logical, intent(in), optional :: l_trans_zen_correction
!   Apply transmission based solar zenith angle correction, DOI 10.1002/qj.385
logical, intent(in), optional :: l_orog_corr_rts
!   Orographic correction applied for the radiation timestep
real(RealExt), intent(in), optional :: orog_corr_rts(n_profile)
!   Orographic correction factor for the radiation timestep, DOI 10.1002/qj.956

logical, intent(in), optional :: l_grey_emis_correction
!   Apply surface temperature correction with grey emissivity per tile
real(RealExt), intent(in), optional :: grey_albedo_tile(:, :)
!   Grey albedo of tiles (n_profile, n_tile)
real(RealExt), intent(in), optional :: frac_tile(:, :)
!   Tile fractions (n_profile, n_tile)
real(RealExt), intent(in), optional :: t_tile(:, :)
!   Tile temperatures (n_profile, n_tile)
real(RealExt), intent(in), optional :: grey_albedo_tile_1d(:)
!   1d grey albedo of tiles (n_tile)
real(RealExt), intent(in), optional :: frac_tile_1d(:)
!   1d tile fractions (n_tile)
real(RealExt), intent(in), optional :: t_tile_1d(:)
!   1d tile temperatures (n_tile)

logical, intent(in), optional :: l_profile_last
!   Loop over profiles is last in 1d fields

logical, intent(in), optional :: l_debug
integer, intent(in), optional :: i_profile_debug
!   Options for outputting debugging information

! Input radiation timestep fields:
real(RealExt), intent(in), optional :: heating_rate_rts(n_profile, n_layer)
real(RealExt), intent(in), optional :: heating_rate_1d_rts(:)
!   Heating rate (Ks-1)
real(RealExt), intent(in), optional :: flux_up_tile_rts(:, :)
real(RealExt), intent(in), optional :: flux_up_tile_1d_rts(:)
!   Upwards flux on tiles (Wm-2) (n_profile, n_tile) and (n_tile)
real(RealExt), intent(in), optional :: flux_up_blue_tile_rts(:, :)
real(RealExt), intent(in), optional :: flux_up_blue_tile_1d_rts(:)
!   Upwards blue flux on tiles (Wm-2)
real(RealExt), intent(in), optional :: flux_direct_toa_rts(n_profile)
!   Direct flux at top-of-atmosphere
real(RealExt), intent(in), optional :: flux_up_toa_rts(n_profile)
!   Upward flux at top-of-atmosphere
real(RealExt), intent(in), optional :: flux_direct_surf_rts(n_profile)
!   Direct flux at the surface
real(RealExt), intent(in), optional :: flux_down_surf_rts(n_profile)
!   Total downward flux at the surface
real(RealExt), intent(in), optional :: flux_up_surf_rts(n_profile)
!   Upward flux at the surface
real(RealExt), intent(in), optional :: flux_direct_blue_surf_rts(n_profile)
!   Direct blue flux at the surface
real(RealExt), intent(in), optional :: flux_down_blue_surf_rts(n_profile)
!   Total downward blue flux at the surface

! Output model timestep fields:
real(RealExt), intent(out), optional :: heating_rate_mts(n_profile, n_layer)
real(RealExt), intent(out), optional :: heating_rate_1d_mts(:)
!   Heating rate (Ks-1)
real(RealExt), intent(out), optional :: flux_up_tile_mts(:, :)
real(RealExt), intent(out), optional :: flux_up_tile_1d_mts(:)
!   Upwards flux on tiles (Wm-2) (n_profile, n_tile) and (n_tile)
real(RealExt), intent(out), optional :: flux_up_blue_tile_mts(:, :)
real(RealExt), intent(out), optional :: flux_up_blue_tile_1d_mts(:)
!   Upwards blue flux on tiles (Wm-2)
real(RealExt), intent(out), optional :: flux_direct_toa_mts(n_profile)
!   Direct flux at top-of-atmosphere
real(RealExt), intent(out), optional :: flux_up_toa_mts(n_profile)
!   Upward flux at top-of-atmosphere
real(RealExt), intent(out), optional :: flux_direct_surf_mts(n_profile)
!   Direct flux at the surface
real(RealExt), intent(out), optional :: flux_down_surf_mts(n_profile)
!   Total downward flux at the surface
real(RealExt), intent(out), optional :: flux_up_surf_mts(n_profile)
!   Upward flux at the surface
real(RealExt), intent(out), optional :: flux_direct_blue_surf_mts(n_profile)
!   Direct blue flux at the surface
real(RealExt), intent(out), optional :: flux_down_blue_surf_mts(n_profile)
!   Total downward blue flux at the surface

! Local variables
integer :: i, l, ll, stride_layer
logical :: l_last
real(RealExt) :: cos_zen_scaling(n_profile), trans_zen_correction(n_profile)
real(RealExt) :: orog_corr(n_profile)
real(RealExt) :: scaling(n_profile)
real(RealExt) :: eps = epsilon(1.0_RealExt)


if (present(l_profile_last)) then
  l_last = l_profile_last
else
  l_last = .false.
end if

! Set the number of layers in the 1d arrays
if (present(n_layer_stride)) then
  stride_layer = n_layer_stride
else
  stride_layer = n_layer
end if


if (present(l_orog_corr_rts)) then
  if (l_orog_corr_rts) then
    ! If the orographic correction has been applied to the surface fluxes
    ! this needs to be accounted for in the transmission-based solar zenith
    ! angle correction.
    orog_corr = orog_corr_rts
  else
    orog_corr = 1.0_RealExt
  end if
else
  orog_corr = 1.0_RealExt
end if


trans_zen_correction = 1.0_RealExt
if (present(l_trans_zen_correction)) then
if (l_trans_zen_correction) then
  ! Transmission-based solar zenith angle correction for surface fluxes
  ! as described in Manners et al 2009, section 3.3 (DOI: 10.1002/qj.385).
  where (cos_zen_rts > eps .and. cos_zen_mts > eps .and. &
         flux_direct_surf_rts > eps .and. &
         flux_down_surf_rts > eps .and. &
         orog_corr*flux_direct_toa_rts > flux_direct_surf_rts .and. &
         orog_corr > sqrt(eps))
    trans_zen_correction = 1.0_RealExt &
      + (orog_corr - 0.5_RealExt) * (flux_direct_toa_rts &
      * (flux_direct_surf_rts/(orog_corr*flux_direct_toa_rts)) &
      **(cos_zen_rts/cos_zen_mts) &
      - flux_direct_surf_rts/orog_corr) / flux_down_surf_rts
  end where
end if
end if


if (present(l_cos_zen_correction)) then
if (l_cos_zen_correction) then
  ! A simple solar zenith angle correction that scales the fluxes and
  ! heating rates by the change in the cosine of the solar zenith angle.
  where (cos_zen_rts*lit_frac_rts > eps)
    cos_zen_scaling = cos_zen_mts*lit_frac_mts / (cos_zen_rts*lit_frac_rts)
  elsewhere
    cos_zen_scaling = 0.0_RealExt
  end where

  scaling = cos_zen_scaling
  call scale_field( heating_rate_rts, heating_rate_mts )
  call scale_field_1d_layer( heating_rate_1d_rts, heating_rate_1d_mts )
  call scale_field_surf( flux_direct_toa_rts, flux_direct_toa_mts )
  call scale_field_surf( flux_up_toa_rts, flux_up_toa_mts )

  ! Surface fields may also be adjusted for the transmission-based solar
  ! zenith angle correction. Note: we apply the same correction to the
  ! total and direct fluxes here to maintain the ratio of direct to diffuse
  ! flux over the radiation timestep. Using the separate correction to the
  ! direct flux as outlined in Manners et al 2009 can in some cases reduce
  ! the accuracy of this ratio.
  scaling = cos_zen_scaling * trans_zen_correction

  call scale_field( flux_up_tile_rts,      flux_up_tile_mts      )
  call scale_field( flux_up_blue_tile_rts, flux_up_blue_tile_mts )

  call scale_field_1d_tile( flux_up_tile_1d_rts,      flux_up_tile_1d_mts      )
  call scale_field_1d_tile( flux_up_blue_tile_1d_rts, flux_up_blue_tile_1d_mts )

  call scale_field_surf( flux_direct_surf_rts,      flux_direct_surf_mts      )
  call scale_field_surf( flux_down_surf_rts,        flux_down_surf_mts        )
  call scale_field_surf( flux_up_surf_rts,          flux_up_surf_mts          )
  call scale_field_surf( flux_direct_blue_surf_rts, flux_direct_blue_surf_mts )
  call scale_field_surf( flux_down_blue_surf_rts,   flux_down_blue_surf_mts   )

  ! Top-of-atmosphere outgoing flux is adjusted to reflect the change in
  ! net surface flux assuming that atmospheric absorption is unchanged.
  if (present(flux_up_toa_mts) .and. &
      present(flux_down_surf_rts) .and. present(flux_up_surf_rts)) then
    do i=1, n_profile
      flux_up_toa_mts(i) = max( 0.0_RealExt, flux_up_toa_mts(i) &
        - (trans_zen_correction(i) - 1.0_RealExt) * cos_zen_scaling(i) &
        * (flux_down_surf_rts(i) - flux_up_surf_rts(i)) )
    end do
  end if

  if (present(l_debug)) then
    if (l_debug) then
      write(9000,'(A)') 'heating_rate_1d_rts, heating_rate_1d_mts'
      do i=1, n_layer
        write(9000,'(2(1pe16.8))') &
          heating_rate_1d_rts(i), heating_rate_1d_mts(i)
      end do
    end if
  end if

end if
end if


if (present(l_grey_emis_correction) .and. &
    present(flux_down_surf_rts)) then
if (l_grey_emis_correction) then
  ! A surface temperature correction with grey emissivity per tile.
  ! Only the upward fluxes are corrected.
  scaling = 1.0_RealExt

  call scale_field( heating_rate_rts, heating_rate_mts )
  call scale_field_1d_layer( heating_rate_1d_rts, heating_rate_1d_mts )
  call scale_field_surf( flux_down_surf_rts, flux_down_surf_mts )

  if (present(flux_up_tile_mts) .and. present(n_tile) .and. &
      present(grey_albedo_tile) .and. present(t_tile)) then
    do i=1, n_tile
      flux_up_tile_mts(1:n_profile, i) &
        = flux_down_surf_rts(:) * grey_albedo_tile(1:n_profile, i) &
        + (1.0_RealExt - grey_albedo_tile(1:n_profile, i)) &
        * real(stefan_boltzmann, RealExt) * t_tile(1:n_profile, i)**4
    end do
    if (present(flux_up_surf_mts) .and. present(frac_tile)) then
      flux_up_surf_mts = sum(flux_up_tile_mts(1:n_profile, 1:n_tile) &
                       * frac_tile(1:n_profile, 1:n_tile), 2)
    end if
  end if
  if (present(flux_up_tile_1d_mts) .and. present(n_tile) .and. &
      present(grey_albedo_tile_1d) .and. present(t_tile_1d)) then
    if (l_last) then
      do l=1, n_profile
        do i=1, n_tile
          ll = n_tile*(l-1) + i
          flux_up_tile_1d_mts(ll) &
            = flux_down_surf_rts(l) * grey_albedo_tile_1d(ll) &
            + (1.0_RealExt - grey_albedo_tile_1d(ll)) &
            * real(stefan_boltzmann, RealExt) * t_tile_1d(ll)**4
        end do
      end do
    else
      do i=1, n_tile
        do l=1, n_profile
          ll = n_profile*(i-1) + l
          flux_up_tile_1d_mts(ll) &
            = flux_down_surf_rts(l) * grey_albedo_tile_1d(ll) &
            + (1.0_RealExt - grey_albedo_tile_1d(ll)) &
            * real(stefan_boltzmann, RealExt) * t_tile_1d(ll)**4
        end do
      end do
    end if
    if (present(flux_up_surf_mts) .and. present(frac_tile_1d)) then
      if (l_last) then
        do l=1, n_profile
          ll = n_tile*(l-1)
          flux_up_surf_mts(l) = sum(flux_up_tile_1d_mts(ll+1:ll+n_tile) &
                              * frac_tile_1d(ll+1:ll+n_tile))
        end do
      else
        ll = n_profile * n_tile
        do l=1, n_profile
          flux_up_surf_mts(l) = sum(flux_up_tile_1d_mts(l:ll:n_profile) &
                              * frac_tile_1d(l:ll:n_profile))
        end do
      end if
    end if
  end if

  ! Top-of-atmosphere outgoing flux is adjusted to reflect the change in
  ! upwards surface flux assuming that atmospheric absorption is unchanged.
  if (present(flux_up_toa_mts) .and. present(flux_up_surf_mts) .and. &
      present(flux_up_toa_rts) .and. present(flux_up_surf_rts)) then
    flux_up_toa_mts = flux_up_toa_rts + flux_up_surf_mts - flux_up_surf_rts
  end if

end if
end if

contains


subroutine scale_field(field_rts, field_mts)

  implicit none

  real(RealExt), intent(in), optional :: field_rts(:, :)
  real(RealExt), intent(out), optional :: field_mts(:, :)

  if (present(field_rts).and.present(field_mts)) then
    do i=1, size(field_rts, 2)
      field_mts(:, i) = field_rts(:, i) * scaling(:)
    end do
  end if

end subroutine scale_field


subroutine scale_field_1d_layer(field_rts, field_mts)

  implicit none

  real(RealExt), intent(in), optional :: field_rts(:)
  real(RealExt), intent(out), optional :: field_mts(:)

  if (present(field_rts).and.present(field_mts)) then
    if (l_last) then
      do l=1, n_profile
        do i=1, n_layer
          ll = stride_layer*(l-1) + i
          field_mts(ll) = field_rts(ll) * scaling(l)
        end do
      end do
    else
      do i=1, n_layer
        do l=1, n_profile
          ll = n_profile*(i-1) + l
          field_mts(ll) = field_rts(ll) * scaling(l)
        end do
      end do
    end if
  end if

end subroutine scale_field_1d_layer


subroutine scale_field_1d_tile(field_rts, field_mts)

  implicit none

  real(RealExt), intent(in), optional :: field_rts(:)
  real(RealExt), intent(out), optional :: field_mts(:)

  if (present(n_tile).and.present(field_rts).and.present(field_mts)) then
    if (l_last) then
      do l=1, n_profile
        do i=1, n_tile
          ll = n_tile*(l-1) + i
          field_mts(ll) = field_rts(ll) * scaling(l)
        end do
      end do
    else
      do i=1, n_tile
        do l=1, n_profile
          ll = n_profile*(i-1) + l
          field_mts(ll) = field_rts(ll) * scaling(l)
        end do
      end do
    end if
  end if

end subroutine scale_field_1d_tile


subroutine scale_field_surf(field_rts, field_mts)

  implicit none

  real(RealExt), intent(in), optional :: field_rts(:)
  real(RealExt), intent(out), optional :: field_mts(:)

  if (present(field_rts).and.present(field_mts)) then
    field_mts(:) = field_rts(:) * scaling(:)
  end if

end subroutine scale_field_surf


end subroutine bones
end module socrates_bones
