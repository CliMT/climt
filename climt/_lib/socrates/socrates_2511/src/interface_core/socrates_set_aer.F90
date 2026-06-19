! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the variables in the Socrates aerosol type
!
!------------------------------------------------------------------------------
module socrates_set_aer
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_AER'
contains

subroutine set_aer(aer, control, dimen, spectrum, &
  n_profile, n_layer, n_aer_mode, profile_list, n_layer_stride, n_aer_layer, &
  aer_mix_ratio, aer_absorption, aer_scattering, aer_asymmetry, &
  aer_mix_ratio_1d, aer_absorption_1d, aer_scattering_1d, aer_asymmetry_1d, &
  mean_rel_humidity, mean_rel_humidity_1d, &
  l_water_soluble, water_soluble, water_soluble_1d, &
  l_dust_like, dust_like, dust_like_1d, &
  l_oceanic, oceanic, oceanic_1d, &
  l_soot, soot, soot_1d, &
  l_ash, ash, ash_1d, &
  l_sulphuric, sulphuric, sulphuric_1d, &
  l_ammonium_sulphate, ammonium_sulphate, ammonium_sulphate_1d, &
  l_saharan_dust, saharan_dust, saharan_dust_1d, &
  l_accum_sulphate, accum_sulphate, accum_sulphate_1d, &
  l_aitken_sulphate, aitken_sulphate, aitken_sulphate_1d, &
  l_fresh_soot, fresh_soot, fresh_soot_1d, &
  l_aged_soot, aged_soot, aged_soot_1d, &
  l_sodium_chloride, sodium_chloride, sodium_chloride_1d, &
  l_seasalt_film, seasalt_film, seasalt_film_1d, &
  l_seasalt_jet, seasalt_jet, seasalt_jet_1d, &
  l_dust_div1, dust_div1, dust_div1_1d, &
  l_dust_div2, dust_div2, dust_div2_1d, &
  l_dust_div3, dust_div3, dust_div3_1d, &
  l_dust_div4, dust_div4, dust_div4_1d, &
  l_dust_div5, dust_div5, dust_div5_1d, &
  l_dust_div6, dust_div6, dust_div6_1d, &
  l_biomass_1, biomass_1, biomass_1_1d, &
  l_biomass_2, biomass_2, biomass_2_1d, &
  l_biogenic, biogenic, biogenic_1d, &
  l_ocff_fresh, ocff_fresh, ocff_fresh_1d, &
  l_ocff_aged, ocff_aged, ocff_aged_1d, &
  l_delta, delta, delta_1d, &
  l_murk, murk, murk_1d, &
  l_nitrate, nitrate, nitrate_1d, &
  l_twobindust_1, twobindust_1, twobindust_1_1d, &
  l_twobindust_2, twobindust_2, twobindust_2_1d, &
  l_invert, l_profile_last)

use def_aer,      only: StrAer, allocate_aer, allocate_aer_prsc
use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use realtype_rd,  only: RealK, RealExt
use rad_pcf,      only: ip_aersrc_classic_ron, ip_aersrc_classic_roff, &
  ip_water_soluble, ip_dust_like, ip_oceanic, ip_soot, ip_ash, ip_sulphuric, &
  ip_ammonium_sulphate, ip_saharan_dust, &
  ip_accum_sulphate, ip_aitken_sulphate, &
  ip_fresh_soot, ip_aged_soot, &
  ip_sodium_chloride, ip_seasalt_film, ip_seasalt_jet, &
  ip_dust_div1, ip_dust_div2, ip_dust_div3, &
  ip_dust_div4, ip_dust_div5, ip_dust_div6, &
  ip_biomass_1, ip_biomass_2, &
  ip_biogenic, &
  ip_ocff_fresh, ip_ocff_aged, &
  ip_delta, ip_murk, &
  ip_nitrate, &
  ip_twobindust_1, ip_twobindust_2

implicit none


! Aerosol properties:
type(StrAer),      intent(out) :: aer

! Control options:
type(StrCtrl),     intent(in)  :: control

! Dimensions:
type(StrDim),      intent(in)  :: dimen

! Spectral data:
type(StrSpecData), intent(in)  :: spectrum

integer, intent(in) :: n_profile
integer, intent(in) :: n_layer

integer, intent(in), optional :: n_aer_mode
!   Number of aerosol modes

integer, intent(in), optional :: profile_list(:)
!   List of profiles to use from input fields

integer, intent(in), optional :: n_layer_stride
!   Number of layers in input 1d arrays

integer, intent(in), optional :: n_aer_layer
!   Number of aerosol layers in 1d arrays (if different to other 1D fields)

real(RealExt), intent(in), optional :: aer_mix_ratio(:, :, :)
!   MODE aerosol mass-mixing ratio (n_profile, n_layer, n_mode)

real(RealExt), intent(in), optional :: aer_mix_ratio_1d(:)
!   1d MODE aerosol mass-mixing ratio (n_aer_layer*n_mode)

real(RealExt), intent(in), dimension(:, :, :, :), optional :: &
  aer_absorption, aer_scattering, aer_asymmetry
!   MODE aerosol optical properties (n_profile, n_layer, n_mode, n_band)

real(RealExt), intent(in), dimension(:), optional :: &
  aer_absorption_1d, aer_scattering_1d, aer_asymmetry_1d
!   1d MODE aerosol optical properties (n_aer_layer*n_mode*n_band)

real(RealExt), intent(in), optional :: mean_rel_humidity(:, :)
real(RealExt), intent(in), optional :: mean_rel_humidity_1d(:)
!   Mean relative humidity applicable for CLASSIC aerosols (clear-sky)

logical, intent(in), optional :: &
  l_water_soluble, l_dust_like, l_oceanic, l_soot, l_ash, l_sulphuric, &
  l_ammonium_sulphate, l_saharan_dust, &
  l_accum_sulphate, l_aitken_sulphate, &
  l_fresh_soot, l_aged_soot, &
  l_sodium_chloride, l_seasalt_film, l_seasalt_jet, &
  l_dust_div1, l_dust_div2, l_dust_div3, &
  l_dust_div4, l_dust_div5, l_dust_div6, &
  l_biomass_1, l_biomass_2, &
  l_biogenic, &
  l_ocff_fresh, l_ocff_aged, &
  l_delta, l_murk, &
  l_nitrate, &
  l_twobindust_1, l_twobindust_2
! Flags to include CLASSIC aerosols

real(RealExt), intent(in), dimension(:, :), optional :: &
  water_soluble, dust_like, oceanic, soot, ash, sulphuric, &
  ammonium_sulphate, saharan_dust, &
  accum_sulphate, aitken_sulphate, &
  fresh_soot, aged_soot, &
  sodium_chloride, seasalt_film, seasalt_jet, &
  dust_div1, dust_div2, dust_div3, &
  dust_div4, dust_div5, dust_div6, &
  biomass_1, biomass_2, &
  biogenic, &
  ocff_fresh, ocff_aged, &
  delta, murk, &
  nitrate, &
  twobindust_1, twobindust_2
! CLASSIC aerosol mass mixing ratios

real(RealExt), intent(in), dimension(:), optional :: &
  water_soluble_1d, dust_like_1d, oceanic_1d, soot_1d, ash_1d, sulphuric_1d, &
  ammonium_sulphate_1d, saharan_dust_1d, &
  accum_sulphate_1d, aitken_sulphate_1d, &
  fresh_soot_1d, aged_soot_1d, &
  sodium_chloride_1d, seasalt_film_1d, seasalt_jet_1d, &
  dust_div1_1d, dust_div2_1d, dust_div3_1d, &
  dust_div4_1d, dust_div5_1d, dust_div6_1d, &
  biomass_1_1d, biomass_2_1d, &
  biogenic_1d, &
  ocff_fresh_1d, ocff_aged_1d, &
  delta_1d, murk_1d, &
  nitrate_1d, &
  twobindust_1_1d, twobindust_2_1d
! 1d CLASSIC aerosol mass mixing ratios

logical, intent(in), optional :: l_invert
!   Flag to invert fields in the vertical
logical, intent(in), optional :: l_profile_last
!   Loop over profiles is last in input fields


! Local variables.
integer :: i, ii, j, l, ll, list(n_profile)
integer :: layer_offset, stride_layer
logical :: l_last


! Allocate structure for the core radiation code interface
call allocate_aer(aer, dimen, spectrum)
call allocate_aer_prsc(aer, dimen, spectrum)

if (present(profile_list)) then
  list = profile_list(1:n_profile)
else
  do l=1, n_profile
    list(l) = l
  end do
end if

layer_offset = 0
if (present(l_invert)) then
  if (l_invert) then
    ! The layer is indexed using an inverted loop counter
    layer_offset = n_layer + 1
  end if
end if

if (present(l_profile_last)) then
  l_last = l_profile_last
else
  l_last = .false.
end if

! Set the number of layers in the 1d aerosol arrays
if (present(n_aer_layer)) then
  stride_layer = n_aer_layer
else if (present(n_layer_stride)) then
  stride_layer = n_layer_stride
else
  stride_layer = n_layer
end if

! Set CLASSIC aerosol mixing ratios and relative humidity
if (control%l_aerosol) then
  call set_aer_field(aer%mean_rel_humidity, mean_rel_humidity, &
                                            mean_rel_humidity_1d)
  do j = 1, spectrum%aerosol%n_aerosol
    ! Here the values in the mixing ratio array will always be in
    ! the order set out in the spectral file:
    aer%mr_type_index(j) = j

    select case ( spectrum%aerosol%type_aerosol(j) )
    case(ip_water_soluble)
      call set_aer_mix_ratio(l_water_soluble, water_soluble, water_soluble_1d)
    case(ip_dust_like)
      call set_aer_mix_ratio(l_dust_like, dust_like, dust_like_1d)
    case(ip_oceanic)
      call set_aer_mix_ratio(l_oceanic, oceanic, oceanic_1d)
    case(ip_soot)
      call set_aer_mix_ratio(l_soot, soot, soot_1d)
    case(ip_ash)
      call set_aer_mix_ratio(l_ash, ash, ash_1d)
    case(ip_sulphuric)
      call set_aer_mix_ratio(l_sulphuric, sulphuric, sulphuric_1d)
    case(ip_ammonium_sulphate)
      call set_aer_mix_ratio(l_ammonium_sulphate, ammonium_sulphate, &
                                                  ammonium_sulphate_1d)
    case(ip_saharan_dust)
      call set_aer_mix_ratio(l_saharan_dust, saharan_dust, saharan_dust_1d)
    case(ip_accum_sulphate)
      call set_aer_mix_ratio(l_accum_sulphate, accum_sulphate, &
                                               accum_sulphate_1d)
    case(ip_aitken_sulphate)
      call set_aer_mix_ratio(l_aitken_sulphate, aitken_sulphate, &
                                                aitken_sulphate_1d)
    case(ip_fresh_soot)
      call set_aer_mix_ratio(l_fresh_soot, fresh_soot, fresh_soot_1d)
    case(ip_aged_soot)
      call set_aer_mix_ratio(l_aged_soot, aged_soot, aged_soot_1d)
    case(ip_sodium_chloride)
      call set_aer_mix_ratio(l_sodium_chloride, sodium_chloride, &
                                                sodium_chloride_1d)
    case(ip_seasalt_film)
      call set_aer_mix_ratio(l_seasalt_film, seasalt_film, seasalt_film_1d)
    case(ip_seasalt_jet)
      call set_aer_mix_ratio(l_seasalt_jet, seasalt_jet, seasalt_jet_1d)
    case(ip_dust_div1)
      call set_aer_mix_ratio(l_dust_div1, dust_div1, dust_div1_1d)
    case(ip_dust_div2)
      call set_aer_mix_ratio(l_dust_div2, dust_div2, dust_div2_1d)
    case(ip_dust_div3)
      call set_aer_mix_ratio(l_dust_div3, dust_div3, dust_div3_1d)
    case(ip_dust_div4)
      call set_aer_mix_ratio(l_dust_div4, dust_div4, dust_div4_1d)
    case(ip_dust_div5)
      call set_aer_mix_ratio(l_dust_div5, dust_div5, dust_div5_1d)
    case(ip_dust_div6)
      call set_aer_mix_ratio(l_dust_div6, dust_div6, dust_div6_1d)
    case(ip_biomass_1)
      call set_aer_mix_ratio(l_biomass_1, biomass_1, biomass_1_1d)
    case(ip_biomass_2)
      call set_aer_mix_ratio(l_biomass_2, biomass_2, biomass_2_1d)
    case(ip_biogenic)
      call set_aer_mix_ratio(l_biogenic, biogenic, biogenic_1d)
    case(ip_ocff_fresh)
      call set_aer_mix_ratio(l_ocff_fresh, ocff_fresh, ocff_fresh_1d)
    case(ip_ocff_aged)
      call set_aer_mix_ratio(l_ocff_aged, ocff_aged, ocff_aged_1d)
    case(ip_delta)
      call set_aer_mix_ratio(l_delta, delta, delta_1d)
    case(ip_murk)
      call set_aer_mix_ratio(l_murk, murk, murk_1d)
    case(ip_nitrate)
      call set_aer_mix_ratio(l_nitrate, nitrate, nitrate_1d)
    case(ip_twobindust_1)
      call set_aer_mix_ratio(l_twobindust_1, twobindust_1, twobindust_1_1d)
    case(ip_twobindust_2)
      call set_aer_mix_ratio(l_twobindust_2, twobindust_2, twobindust_2_1d)
    case default
      aer%mr_source(j) = ip_aersrc_classic_roff
    end select
  end do
end if

! Set MODE aerosol optical properties
if (present(n_aer_mode).and.control%l_aerosol_mode) then
  aer%n_mode = n_aer_mode
else
  aer%n_mode = 0
end if

if (aer%n_mode > 0) then
  call set_mode_field(aer%mode_mix_ratio, aer_mix_ratio, aer_mix_ratio_1d)
  call set_mode_opt_prop(aer%mode_absorption, aer_absorption, aer_absorption_1d)
  call set_mode_opt_prop(aer%mode_scattering, aer_scattering, aer_scattering_1d)
  call set_mode_opt_prop(aer%mode_asymmetry,  aer_asymmetry,  aer_asymmetry_1d)
end if


contains

  subroutine set_aer_field(out_field, full_field, oned_field)
    implicit none

    real(RealK), intent(inout) :: out_field(:, :)
!     Output field
    real(RealExt), intent(in), optional :: full_field(:, :)
!     Full field variable
    real(RealExt), intent(in), optional :: oned_field(:)
!     One-dimensional variable

    if (present(full_field)) then
      if (l_last) then
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            out_field(l, i) = real(full_field(ii, list(l)), RealK)
          end do
        end do
      else
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            out_field(l, i) = real(full_field(list(l), ii), RealK)
          end do
        end do
      end if
    else if (present(oned_field)) then
      if (l_last) then
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            ll = stride_layer*(list(l)-1) + ii
            out_field(l, i) = real(oned_field(ll), RealK)
          end do
        end do
      else
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            ll = n_profile*(ii-1) + list(l)
            out_field(l, i) = real(oned_field(ll), RealK)
          end do
        end do
      end if
    else
      out_field = 0.0_RealK
    end if

  end subroutine set_aer_field


  subroutine set_aer_mix_ratio(l_field, full_field, oned_field)
    implicit none

    logical, intent(in), optional :: l_field
!     Flag to use this aerosol
    real(RealExt), intent(in), optional :: full_field(:, :)
!     Full field variable
    real(RealExt), intent(in), optional :: oned_field(:)
!     One-dimensional variable

    logical :: l_use_field

    if (present(l_field)) then
      l_use_field = l_field
    else
      l_use_field = .false.
    end if

    if (l_use_field) then
      ! Indicate this aerosol is radiatively active
      aer%mr_source(j) = ip_aersrc_classic_ron
      if (present(full_field)) then
        if (l_last) then
          do i=1, n_layer
            ii = abs(layer_offset-i)
            do l=1, n_profile
              aer%mix_ratio(l, i, j) = real(full_field(ii, list(l)), RealK)
            end do
          end do
        else
          do i=1, n_layer
            ii = abs(layer_offset-i)
            do l=1, n_profile
              aer%mix_ratio(l, i, j) = real(full_field(list(l), ii), RealK)
            end do
          end do
        end if
      else if (present(oned_field)) then
        if (l_last) then
          do i=1, n_layer
            ii = abs(layer_offset-i)
            do l=1, n_profile
              ll = stride_layer*(list(l)-1) + ii
              aer%mix_ratio(l, i, j) = real(oned_field(ll), RealK)
            end do
          end do
        else
          do i=1, n_layer
            ii = abs(layer_offset-i)
            do l=1, n_profile
              ll = n_profile*(ii-1) + list(l)
              aer%mix_ratio(l, i, j) = real(oned_field(ll), RealK)
            end do
          end do
        end if
      else
        ! Turn off if no mixing ratio data provided
        aer%mr_source(j) = ip_aersrc_classic_roff
      end if
    else
      ! Indicate this aerosol is not radiatively active
      aer%mr_source(j) = ip_aersrc_classic_roff
    end if

  end subroutine set_aer_mix_ratio


  subroutine set_mode_field(out_field, full_field, oned_field)
    implicit none

    real(RealK), intent(inout) :: out_field(:, :, :)
!     Output field
    real(RealExt), intent(in), optional :: full_field(:, :, :)
!     Full field variable
    real(RealExt), intent(in), optional :: oned_field(:)
!     One-dimensional variable

    if (present(full_field)) then
      if (l_last) then
        do j=1, aer%n_mode
          do i=1, n_layer
            ii = abs(layer_offset-i)
            do l=1, n_profile
              out_field(l, i, j) = real(full_field(ii, j, list(l)), RealK)
            end do
          end do
        end do
      else
        do j=1, aer%n_mode
          do i=1, n_layer
            ii = abs(layer_offset-i)
            do l=1, n_profile
              out_field(l, i, j) = real(full_field(list(l), ii, j), RealK)
            end do
          end do
        end do
      end if
    else if (present(oned_field)) then
      if (l_last) then
        do j=1, aer%n_mode
          do i=1, n_layer
            ii = abs(layer_offset-i)
            do l=1, n_profile
              ll = aer%n_mode*stride_layer*(list(l)-1) &
                 + stride_layer*(j-1) + ii
              out_field(l, i, j) = real(oned_field(ll), RealK)
            end do
          end do
        end do
      else
        do j=1, aer%n_mode
          do i=1, n_layer
            ii = abs(layer_offset-i)
            do l=1, n_profile
              ll = n_profile*stride_layer*(j-1) &
                 + n_profile*(ii-1) + list(l)
              out_field(l, i, j) = real(oned_field(ll), RealK)
            end do
          end do
        end do
      end if
    else
      out_field = 0.0_RealK
    end if

  end subroutine set_mode_field


  subroutine set_mode_opt_prop(out_field, full_field, oned_field)
    implicit none

    real(RealK), intent(inout) :: out_field(:, :, :, :)
!     Output field
    real(RealExt), intent(in), optional :: full_field(:, :, :, :)
!     Full field variable
    real(RealExt), intent(in), optional :: oned_field(:)
!     One-dimensional variable

    integer :: band

    if (present(full_field)) then
      if (l_last) then
        do band=control%first_band, control%last_band
          do j=1, aer%n_mode
            do i=1, n_layer
              ii = abs(layer_offset-i)
              do l=1, n_profile
                out_field(l, i, j, band) &
                  = real(full_field(ii, j, band, list(l)), RealK)
              end do
            end do
          end do
        end do
      else
        do band=control%first_band, control%last_band
          do j=1, aer%n_mode
            do i=1, n_layer
              ii = abs(layer_offset-i)
              do l=1, n_profile
                out_field(l, i, j, band) &
                  = real(full_field(list(l), ii, j, band), RealK)
              end do
            end do
          end do
        end do
      end if
    else if (present(oned_field)) then
      if (l_last) then
        do band=control%first_band, control%last_band
          do j=1, aer%n_mode
            do i=1, n_layer
              ii = abs(layer_offset-i)
              do l=1, n_profile
                ll = spectrum%basic%n_band*aer%n_mode*stride_layer*(list(l)-1) &
                   + aer%n_mode*stride_layer*(band-1) &
                   + stride_layer*(j-1) + ii
                out_field(l, i, j, band) = real(oned_field(ll), RealK)
              end do
            end do
          end do
        end do
      else
        do band=control%first_band, control%last_band
          do j=1, aer%n_mode
            do i=1, n_layer
              ii = abs(layer_offset-i)
              do l=1, n_profile
                ll = aer%n_mode*stride_layer*n_profile*(band-1) &
                   + stride_layer*n_profile*(j-1) &
                   + n_profile*(ii-1) + list(l)
                out_field(l, i, j, band) = real(oned_field(ll), RealK)
              end do
            end do
          end do
        end do
      end if
    else
      out_field = 0.0_RealK
    end if

  end subroutine set_mode_opt_prop


end subroutine set_aer
end module socrates_set_aer
