! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the variables in the Socrates atm type
!
!------------------------------------------------------------------------------
module socrates_set_atm
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_ATM'
contains

subroutine set_atm(atm, control, dimen, spectrum, n_profile, n_layer, &
  profile_list, n_layer_stride, n_level_stride, &
  p_layer, t_layer, mass, density, p_level, t_level, r_layer, r_level, &
  p_layer_1d, t_layer_1d, mass_1d, density_1d, p_level_1d, t_level_1d, &
  r_layer_1d, r_level_1d, &
  h2o, co2, o3, n2o, co, ch4, o2, no, so2, no2, nh3, hno3, n2, &
  cfc11, cfc12, cfc113, hcfc22, hfc125, hfc134a, cfc114, tio, vo, h2, &
  he, ocs, na, k, feh, crh, li, rb, cs, ph3, c2h2, hcn, h2s, ar, o, n, no3, &
  n2o5, hono, ho2no2, h2o2, c2h6, ch3, h2co, ho2, &
  hdo, hcl, hf, cosso, tosso, yosos, &
  ch3cho, ch3ooh, ch3coch3, ch3cocho, chocho, c2h5cho, &
  hoch2cho, c2h5coch3, mvk, macr, pan, ch3ono2, &
  h2o_1d, co2_1d, o3_1d, n2o_1d, co_1d, ch4_1d, o2_1d, no_1d, so2_1d, no2_1d, &
  nh3_1d, hno3_1d, n2_1d, cfc11_1d, cfc12_1d, cfc113_1d, hcfc22_1d, hfc125_1d, &
  hfc134a_1d, cfc114_1d, tio_1d, vo_1d, h2_1d, he_1d, ocs_1d, na_1d, k_1d, &
  feh_1d, crh_1d, li_1d, rb_1d, cs_1d, ph3_1d, c2h2_1d, hcn_1d, h2s_1d, ar_1d, &
  o_1d, n_1d, no3_1d, n2o5_1d, hono_1d, ho2no2_1d, h2o2_1d, c2h6_1d, ch3_1d, &
  h2co_1d, ho2_1d, hdo_1d, hcl_1d, hf_1d, cosso_1d, tosso_1d, yosos_1d, &
  ch3cho_1d, ch3ooh_1d, ch3coch3_1d, ch3cocho_1d, chocho_1d, c2h5cho_1d, &
  hoch2cho_1d, c2h5coch3_1d, mvk_1d, macr_1d, pan_1d, ch3ono2_1d, &
  h2o_mix_ratio, co2_mix_ratio, o3_mix_ratio, n2o_mix_ratio, co_mix_ratio, &
  ch4_mix_ratio, o2_mix_ratio, no_mix_ratio, so2_mix_ratio, no2_mix_ratio, &
  nh3_mix_ratio, hno3_mix_ratio, n2_mix_ratio, cfc11_mix_ratio, &
  cfc12_mix_ratio, cfc113_mix_ratio, hcfc22_mix_ratio, hfc125_mix_ratio, &
  hfc134a_mix_ratio, cfc114_mix_ratio, tio_mix_ratio, vo_mix_ratio, &
  h2_mix_ratio, he_mix_ratio, ocs_mix_ratio, na_mix_ratio, k_mix_ratio, &
  feh_mix_ratio, crh_mix_ratio, li_mix_ratio, rb_mix_ratio, cs_mix_ratio, &
  ph3_mix_ratio, c2h2_mix_ratio, hcn_mix_ratio, h2s_mix_ratio, ar_mix_ratio, &
  o_mix_ratio, n_mix_ratio, no3_mix_ratio, n2o5_mix_ratio, hono_mix_ratio, &
  ho2no2_mix_ratio, h2o2_mix_ratio, c2h6_mix_ratio, ch3_mix_ratio, &
  h2co_mix_ratio, ho2_mix_ratio, hdo_mix_ratio, hcl_mix_ratio, hf_mix_ratio, &
  cosso_mix_ratio, tosso_mix_ratio, yosos_mix_ratio, ch3cho_mix_ratio, &
  ch3ooh_mix_ratio, ch3coch3_mix_ratio, ch3cocho_mix_ratio, chocho_mix_ratio, &
  c2h5cho_mix_ratio, hoch2cho_mix_ratio, c2h5coch3_mix_ratio, mvk_mix_ratio, &
  macr_mix_ratio, pan_mix_ratio, ch3ono2_mix_ratio, &
  l_h2o_well_mixed, l_co2_well_mixed, l_o3_well_mixed, l_n2o_well_mixed, &
  l_co_well_mixed, l_ch4_well_mixed, l_o2_well_mixed, l_no_well_mixed, &
  l_so2_well_mixed, l_no2_well_mixed, l_nh3_well_mixed, l_hno3_well_mixed, &
  l_n2_well_mixed, l_cfc11_well_mixed, l_cfc12_well_mixed, &
  l_cfc113_well_mixed, l_hcfc22_well_mixed, l_hfc125_well_mixed, &
  l_hfc134a_well_mixed, l_cfc114_well_mixed, l_tio_well_mixed, &
  l_vo_well_mixed, l_h2_well_mixed, l_he_well_mixed, l_ocs_well_mixed, &
  l_na_well_mixed, l_k_well_mixed, l_feh_well_mixed, l_crh_well_mixed, &
  l_li_well_mixed, l_rb_well_mixed, l_cs_well_mixed, l_ph3_well_mixed, &
  l_c2h2_well_mixed, l_hcn_well_mixed, l_h2s_well_mixed, l_ar_well_mixed, &
  l_o_well_mixed, l_n_well_mixed, l_no3_well_mixed, l_n2o5_well_mixed, &
  l_hono_well_mixed, l_ho2no2_well_mixed, l_h2o2_well_mixed, &
  l_c2h6_well_mixed, l_ch3_well_mixed, l_h2co_well_mixed, l_ho2_well_mixed, &
  l_hdo_well_mixed, l_hcl_well_mixed, l_hf_well_mixed, l_cosso_well_mixed, &
  l_tosso_well_mixed, l_yosos_well_mixed, l_ch3cho_well_mixed, &
  l_ch3ooh_well_mixed, l_ch3coch3_well_mixed, l_ch3cocho_well_mixed, &
  l_chocho_well_mixed, l_c2h5cho_well_mixed, l_hoch2cho_well_mixed, &
  l_c2h5coch3_well_mixed, l_mvk_well_mixed, l_macr_well_mixed, &
  l_pan_well_mixed, l_ch3ono2_well_mixed, &
  l_invert, l_profile_last, l_debug, i_profile_debug)

use def_atm,      only: StrAtm, allocate_atm
use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use realtype_rd,  only: RealK, RealExt
use gas_list_pcf, only: &
  ip_h2o, ip_co2, ip_o3, ip_n2o, ip_co, ip_ch4, ip_o2, ip_no, ip_so2, ip_no2, &
  ip_nh3, ip_hno3, ip_n2, ip_cfc11, ip_cfc12, ip_cfc113, ip_hcfc22, ip_hfc125, &
  ip_hfc134a, ip_cfc114, ip_tio, ip_vo, ip_h2, ip_he, ip_ocs, ip_na, ip_k, &
  ip_feh, ip_crh, ip_li, ip_rb, ip_cs, ip_ph3, ip_c2h2, ip_hcn, ip_h2s, ip_ar, &
  ip_o, ip_n, ip_no3, ip_n2o5, ip_hono, ip_ho2no2, ip_h2o2, ip_c2h6, ip_ch3, &
  ip_h2co, ip_ho2, ip_hdo, ip_hcl, ip_hf, ip_cosso, ip_tosso, ip_yosos, &
  ip_ch3cho, ip_ch3ooh, ip_ch3coch3, ip_ch3cocho, ip_chocho, ip_c2h5cho, &
  ip_hoch2cho, ip_c2h5coch3, ip_mvk, ip_macr, ip_pan, ip_ch3ono2

implicit none


! Atmospheric properties:
type(StrAtm),      intent(out) :: atm

! Control options:
type(StrCtrl),   intent(inout) :: control

! Dimensions:
type(StrDim),      intent(in)  :: dimen

! Spectral data:
type(StrSpecData), intent(in)  :: spectrum

integer, intent(in) :: n_profile
!   Number of atmospheric profiles for radiation calculations
integer, intent(in) :: n_layer
!   Number of atmospheric layers for radiation calculations
integer, intent(in), optional :: profile_list(:)
!   List of profiles to use from input fields
integer, intent(in), optional :: n_layer_stride
!   Number of layers in input 1d arrays
integer, intent(in), optional :: n_level_stride
!   Number of levels in input 1d arrays

real(RealExt), intent(in), optional :: p_layer(:, :), p_layer_1d(:)
!   Pressure at layer centres
real(RealExt), intent(in), optional :: t_layer(:, :), t_layer_1d(:)
!   Temperature at layer centres
real(RealExt), intent(in), optional :: mass(:, :), mass_1d(:)
!   Mass of layer (kg m-2)
real(RealExt), intent(in), optional :: density(:, :), density_1d(:)
!   Density of layer (kg m-3)
real(RealExt), intent(in), optional :: p_level(:, :), p_level_1d(:)
!   Pressure at layer boundaries
real(RealExt), intent(in), optional :: t_level(:, :), t_level_1d(:)
!   Temperature at layer boundaries
real(RealExt), intent(in), optional :: r_layer(:, :), r_layer_1d(:)
!   Radius (height from centre of planet) at layer centres
real(RealExt), intent(in), optional :: r_level(:, :), r_level_1d(:)
!   Radius (height from centre of planet) at layer boundaries

real(RealExt), intent(in), dimension(:, :), optional :: &
  h2o, co2, o3, n2o, co, ch4, o2, no, so2, no2, nh3, hno3, n2, &
  cfc11, cfc12, cfc113, hcfc22, hfc125, hfc134a, cfc114, tio, vo, h2, &
  he, ocs, na, k, feh, crh, li, rb, cs, ph3, c2h2, hcn, h2s, ar, o, n, no3, &
  n2o5, hono, ho2no2, h2o2, c2h6, ch3, h2co, ho2, &
  hdo, hcl, hf, cosso, tosso, yosos, &
  ch3cho, ch3ooh, ch3coch3, ch3cocho, chocho, c2h5cho, &
  hoch2cho, c2h5coch3, mvk, macr, pan, ch3ono2
!   Full field mass mixing ratios

real(RealExt), intent(in), dimension(:), optional :: &
  h2o_1d, co2_1d, o3_1d, n2o_1d, co_1d, ch4_1d, o2_1d, no_1d, so2_1d, no2_1d, &
  nh3_1d, hno3_1d, n2_1d, cfc11_1d, cfc12_1d, cfc113_1d, hcfc22_1d, hfc125_1d, &
  hfc134a_1d, cfc114_1d, tio_1d, vo_1d, h2_1d, he_1d, ocs_1d, na_1d, k_1d, &
  feh_1d, crh_1d, li_1d, rb_1d, cs_1d, ph3_1d, c2h2_1d, hcn_1d, h2s_1d, ar_1d, &
  o_1d, n_1d, no3_1d, n2o5_1d, hono_1d, ho2no2_1d, h2o2_1d, c2h6_1d, ch3_1d, &
  h2co_1d, ho2_1d, hdo_1d, hcl_1d, hf_1d, cosso_1d, tosso_1d, yosos_1d, &
  ch3cho_1d, ch3ooh_1d, ch3coch3_1d, ch3cocho_1d, chocho_1d, c2h5cho_1d, &
  hoch2cho_1d, c2h5coch3_1d, mvk_1d, macr_1d, pan_1d, ch3ono2_1d
!   1d mass mixing ratios

real(RealExt), intent(in), optional :: &
  h2o_mix_ratio, co2_mix_ratio, o3_mix_ratio, n2o_mix_ratio, co_mix_ratio, &
  ch4_mix_ratio, o2_mix_ratio, no_mix_ratio, so2_mix_ratio, no2_mix_ratio, &
  nh3_mix_ratio, hno3_mix_ratio, n2_mix_ratio, cfc11_mix_ratio, &
  cfc12_mix_ratio, cfc113_mix_ratio, hcfc22_mix_ratio, hfc125_mix_ratio, &
  hfc134a_mix_ratio, cfc114_mix_ratio, tio_mix_ratio, vo_mix_ratio, &
  h2_mix_ratio, he_mix_ratio, ocs_mix_ratio, na_mix_ratio, k_mix_ratio, &
  feh_mix_ratio, crh_mix_ratio, li_mix_ratio, rb_mix_ratio, cs_mix_ratio, &
  ph3_mix_ratio, c2h2_mix_ratio, hcn_mix_ratio, h2s_mix_ratio, ar_mix_ratio, &
  o_mix_ratio, n_mix_ratio, no3_mix_ratio, n2o5_mix_ratio, hono_mix_ratio, &
  ho2no2_mix_ratio, h2o2_mix_ratio, c2h6_mix_ratio, ch3_mix_ratio, &
  h2co_mix_ratio, ho2_mix_ratio, hdo_mix_ratio, hcl_mix_ratio, hf_mix_ratio, &
  cosso_mix_ratio, tosso_mix_ratio, yosos_mix_ratio, ch3cho_mix_ratio, &
  ch3ooh_mix_ratio, ch3coch3_mix_ratio, ch3cocho_mix_ratio, chocho_mix_ratio, &
  c2h5cho_mix_ratio, hoch2cho_mix_ratio, c2h5coch3_mix_ratio, mvk_mix_ratio, &
  macr_mix_ratio, pan_mix_ratio, ch3ono2_mix_ratio
!   Well mixed mass mixing ratios

logical, intent(in), optional :: &
  l_h2o_well_mixed, l_co2_well_mixed, l_o3_well_mixed, l_n2o_well_mixed, &
  l_co_well_mixed, l_ch4_well_mixed, l_o2_well_mixed, l_no_well_mixed, &
  l_so2_well_mixed, l_no2_well_mixed, l_nh3_well_mixed, l_hno3_well_mixed, &
  l_n2_well_mixed, l_cfc11_well_mixed, l_cfc12_well_mixed, &
  l_cfc113_well_mixed, l_hcfc22_well_mixed, l_hfc125_well_mixed, &
  l_hfc134a_well_mixed, l_cfc114_well_mixed, l_tio_well_mixed, &
  l_vo_well_mixed, l_h2_well_mixed, l_he_well_mixed, l_ocs_well_mixed, &
  l_na_well_mixed, l_k_well_mixed, l_feh_well_mixed, l_crh_well_mixed, &
  l_li_well_mixed, l_rb_well_mixed, l_cs_well_mixed, l_ph3_well_mixed, &
  l_c2h2_well_mixed, l_hcn_well_mixed, l_h2s_well_mixed, l_ar_well_mixed, &
  l_o_well_mixed, l_n_well_mixed, l_no3_well_mixed, l_n2o5_well_mixed, &
  l_hono_well_mixed, l_ho2no2_well_mixed, l_h2o2_well_mixed, &
  l_c2h6_well_mixed, l_ch3_well_mixed, l_h2co_well_mixed, l_ho2_well_mixed, &
  l_hdo_well_mixed, l_hcl_well_mixed, l_hf_well_mixed, l_cosso_well_mixed, &
  l_tosso_well_mixed, l_yosos_well_mixed, l_ch3cho_well_mixed, &
  l_ch3ooh_well_mixed, l_ch3coch3_well_mixed, l_ch3cocho_well_mixed, &
  l_chocho_well_mixed, l_c2h5cho_well_mixed, l_hoch2cho_well_mixed, &
  l_c2h5coch3_well_mixed, l_mvk_well_mixed, l_macr_well_mixed, &
  l_pan_well_mixed, l_ch3ono2_well_mixed
!   Flag to use the well mixed ratios

logical, intent(in), optional :: l_invert
!   Flag to invert fields in the vertical
logical, intent(in), optional :: l_profile_last
!   Loop over profiles is last in input fields

logical, intent(in), optional :: l_debug
integer, intent(in), optional :: i_profile_debug
!   Options for outputting debugging information

! Local variables.
integer :: i, ii, l, ll, i_gas, list(n_profile)
integer :: layer_offset, level_offset, stride_layer, stride_level
logical :: l_last

call allocate_atm(atm, dimen, spectrum)

! Set up atmosphere grid
atm%n_profile = n_profile
atm%n_layer   = n_layer

if (present(profile_list)) then
  list = profile_list(1:n_profile)
else
  do l=1, n_profile
    list(l) = l
  end do
end if

layer_offset = 0
level_offset = 0
if (present(l_invert)) then
  if (l_invert) then
    ! The layer is indexed using an inverted loop counter
    layer_offset = n_layer + 1
    level_offset = n_layer
  end if
end if

if (present(l_profile_last)) then
  l_last = l_profile_last
else
  l_last = .false.
end if

! Set the number of layers and levels in the 1d arrays
if (present(n_layer_stride)) then
  stride_layer = n_layer_stride
else
  stride_layer = n_layer
end if
if (present(n_level_stride)) then
  stride_level = n_level_stride
else
  stride_level = n_layer + 1
end if

! Set the pressures, temperatures, masses (per square metre) and densities
call set_layer_field(atm%p, p_layer, p_layer_1d)
call set_layer_field(atm%t, t_layer, t_layer_1d)
call set_layer_field(atm%mass, mass, mass_1d)
call set_layer_field(atm%density, density, density_1d)
call set_layer_field(atm%r_layer, r_layer, r_layer_1d)

call set_level_field(atm%p_level, p_level, p_level_1d)
call set_level_field(atm%t_level, t_level, t_level_1d)
call set_level_field(atm%r_level, r_level, r_level_1d)

! Set gas mass mixing ratios
do i_gas=1, spectrum%gas%n_absorb
  select case(spectrum%gas%type_absorb(i_gas))
  case(ip_h2o)
    call set_gas_mix_ratio(h2o, h2o_1d, h2o_mix_ratio, l_h2o_well_mixed)
  case(ip_co2)
    call set_gas_mix_ratio(co2, co2_1d, co2_mix_ratio, l_co2_well_mixed)
  case(ip_o3)
    call set_gas_mix_ratio(o3, o3_1d, o3_mix_ratio, l_o3_well_mixed)
  case(ip_n2o)
    call set_gas_mix_ratio(n2o, n2o_1d, n2o_mix_ratio, l_n2o_well_mixed)
  case(ip_co)
    call set_gas_mix_ratio(co, co_1d, co_mix_ratio, l_co_well_mixed)
  case(ip_ch4)
    call set_gas_mix_ratio(ch4, ch4_1d, ch4_mix_ratio, l_ch4_well_mixed)
  case(ip_o2)
    call set_gas_mix_ratio(o2, o2_1d, o2_mix_ratio, l_o2_well_mixed)
  case(ip_no)
    call set_gas_mix_ratio(no, no_1d, no_mix_ratio, l_no_well_mixed)
  case(ip_so2)
    call set_gas_mix_ratio(so2, so2_1d, so2_mix_ratio, l_so2_well_mixed)
  case(ip_no2)
    call set_gas_mix_ratio(no2, no2_1d, no2_mix_ratio, l_no2_well_mixed)
  case(ip_nh3)
    call set_gas_mix_ratio(nh3, nh3_1d, nh3_mix_ratio, l_nh3_well_mixed)
  case(ip_hno3)
    call set_gas_mix_ratio(hno3, hno3_1d, hno3_mix_ratio, l_hno3_well_mixed)
  case(ip_n2)
    call set_gas_mix_ratio(n2, n2_1d, n2_mix_ratio, l_n2_well_mixed)
  case(ip_cfc11)
    call set_gas_mix_ratio(cfc11, cfc11_1d, cfc11_mix_ratio, l_cfc11_well_mixed)
  case(ip_cfc12)
    call set_gas_mix_ratio(cfc12, cfc12_1d, cfc12_mix_ratio, l_cfc12_well_mixed)
  case(ip_cfc113)
    call set_gas_mix_ratio(cfc113, cfc113_1d, cfc113_mix_ratio, &
                           l_cfc113_well_mixed)
  case(ip_hcfc22)
    call set_gas_mix_ratio(hcfc22, hcfc22_1d, hcfc22_mix_ratio, &
                           l_hcfc22_well_mixed)
  case(ip_hfc125)
    call set_gas_mix_ratio(hfc125, hfc125_1d, hfc125_mix_ratio, &
                           l_hfc125_well_mixed)
  case(ip_hfc134a)
    call set_gas_mix_ratio(hfc134a, hfc134a_1d, hfc134a_mix_ratio, &
                           l_hfc134a_well_mixed)
  case(ip_cfc114)
    call set_gas_mix_ratio(cfc114, cfc114_1d, cfc114_mix_ratio, &
                           l_cfc114_well_mixed)
  case(ip_tio)
    call set_gas_mix_ratio(tio, tio_1d, tio_mix_ratio, l_tio_well_mixed)
  case(ip_vo)
    call set_gas_mix_ratio(vo, vo_1d, vo_mix_ratio, l_vo_well_mixed)
  case(ip_h2)
    call set_gas_mix_ratio(h2, h2_1d, h2_mix_ratio, l_h2_well_mixed)
  case(ip_he)
    call set_gas_mix_ratio(he, he_1d, he_mix_ratio, l_he_well_mixed)
  case(ip_ocs)
    call set_gas_mix_ratio(ocs, ocs_1d, ocs_mix_ratio, l_ocs_well_mixed)
  case(ip_na)
    call set_gas_mix_ratio(na, na_1d, na_mix_ratio, l_na_well_mixed)
  case(ip_k)
    call set_gas_mix_ratio(k, k_1d, k_mix_ratio, l_k_well_mixed)
  case(ip_feh)
    call set_gas_mix_ratio(feh, feh_1d, feh_mix_ratio, l_feh_well_mixed)
  case(ip_crh)
    call set_gas_mix_ratio(crh, crh_1d, crh_mix_ratio, l_crh_well_mixed)
  case(ip_li)
    call set_gas_mix_ratio(li, li_1d, li_mix_ratio, l_li_well_mixed)
  case(ip_rb)
    call set_gas_mix_ratio(rb, rb_1d, rb_mix_ratio, l_rb_well_mixed)
  case(ip_cs)
    call set_gas_mix_ratio(cs, cs_1d, cs_mix_ratio, l_cs_well_mixed)
  case(ip_ph3)
    call set_gas_mix_ratio(ph3, ph3_1d, ph3_mix_ratio, l_ph3_well_mixed)
  case(ip_c2h2)
    call set_gas_mix_ratio(c2h2, c2h2_1d, c2h2_mix_ratio, l_c2h2_well_mixed)
  case(ip_hcn)
    call set_gas_mix_ratio(hcn, hcn_1d, hcn_mix_ratio, l_hcn_well_mixed)
  case(ip_h2s)
    call set_gas_mix_ratio(h2s, h2s_1d, h2s_mix_ratio, l_h2s_well_mixed)
  case(ip_ar)
    call set_gas_mix_ratio(ar, ar_1d, ar_mix_ratio, l_ar_well_mixed)
  case(ip_o)
    call set_gas_mix_ratio(o, o_1d, o_mix_ratio, l_o_well_mixed)
  case(ip_n)
    call set_gas_mix_ratio(n, n_1d, n_mix_ratio, l_n_well_mixed)
  case(ip_no3)
    call set_gas_mix_ratio(no3, no3_1d, no3_mix_ratio, l_no3_well_mixed)
  case(ip_n2o5)
    call set_gas_mix_ratio(n2o5, n2o5_1d, n2o5_mix_ratio, l_n2o5_well_mixed)
  case(ip_hono)
    call set_gas_mix_ratio(hono, hono_1d, hono_mix_ratio, l_hono_well_mixed)
  case(ip_ho2no2)
    call set_gas_mix_ratio(ho2no2, ho2no2_1d, ho2no2_mix_ratio, &
                           l_ho2no2_well_mixed)
  case(ip_h2o2)
    call set_gas_mix_ratio(h2o2, h2o2_1d, h2o2_mix_ratio, l_h2o2_well_mixed)
  case(ip_c2h6)
    call set_gas_mix_ratio(c2h6, c2h6_1d, c2h6_mix_ratio, l_c2h6_well_mixed)
  case(ip_ch3)
    call set_gas_mix_ratio(ch3, ch3_1d, ch3_mix_ratio, l_ch3_well_mixed)
  case(ip_h2co)
    call set_gas_mix_ratio(h2co, h2co_1d, h2co_mix_ratio, l_h2co_well_mixed)
  case(ip_ho2)
    call set_gas_mix_ratio(ho2, ho2_1d, ho2_mix_ratio, l_ho2_well_mixed)
  case(ip_hdo)
    call set_gas_mix_ratio(hdo, hdo_1d, hdo_mix_ratio, l_hdo_well_mixed)
  case(ip_hcl)
    call set_gas_mix_ratio(hcl, hcl_1d, hcl_mix_ratio, l_hcl_well_mixed)
  case(ip_hf)
    call set_gas_mix_ratio(hf, hf_1d, hf_mix_ratio, l_hf_well_mixed)
  case(ip_cosso)
    call set_gas_mix_ratio(cosso, cosso_1d, cosso_mix_ratio, l_cosso_well_mixed)
  case(ip_tosso)
    call set_gas_mix_ratio(tosso, tosso_1d, tosso_mix_ratio, l_tosso_well_mixed)
  case(ip_yosos)
    call set_gas_mix_ratio(yosos, yosos_1d, yosos_mix_ratio, l_yosos_well_mixed)
  case(ip_ch3cho)
    call set_gas_mix_ratio(ch3cho, ch3cho_1d, ch3cho_mix_ratio, &
                           l_ch3cho_well_mixed)
  case(ip_ch3ooh)
    call set_gas_mix_ratio(ch3ooh, ch3ooh_1d, ch3ooh_mix_ratio, &
                           l_ch3ooh_well_mixed)
  case(ip_ch3coch3)
    call set_gas_mix_ratio(ch3coch3, ch3coch3_1d, ch3coch3_mix_ratio, &
                           l_ch3coch3_well_mixed)
  case(ip_ch3cocho)
    call set_gas_mix_ratio(ch3cocho, ch3cocho_1d, ch3cocho_mix_ratio, &
                           l_ch3cocho_well_mixed)
  case(ip_chocho)
    call set_gas_mix_ratio(chocho, chocho_1d, chocho_mix_ratio, &
                           l_chocho_well_mixed)
  case(ip_c2h5cho)
    call set_gas_mix_ratio(c2h5cho, c2h5cho_1d, c2h5cho_mix_ratio, &
                           l_c2h5cho_well_mixed)
  case(ip_hoch2cho)
    call set_gas_mix_ratio(hoch2cho, hoch2cho_1d, hoch2cho_mix_ratio, &
                           l_hoch2cho_well_mixed)
  case(ip_c2h5coch3)
    call set_gas_mix_ratio(c2h5coch3, c2h5coch3_1d, c2h5coch3_mix_ratio, &
                           l_c2h5coch3_well_mixed)
  case(ip_mvk)
    call set_gas_mix_ratio(mvk, mvk_1d, mvk_mix_ratio, l_mvk_well_mixed)
  case(ip_macr)
    call set_gas_mix_ratio(macr, macr_1d, macr_mix_ratio, l_macr_well_mixed)
  case(ip_pan)
    call set_gas_mix_ratio(pan, pan_1d, pan_mix_ratio, l_pan_well_mixed)
  case(ip_ch3ono2)
    call set_gas_mix_ratio(ch3ono2, ch3ono2_1d, ch3ono2_mix_ratio, &
                           l_ch3ono2_well_mixed)
  case default
    do i=1, n_layer
      do l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = 0.0_RealK
      end do
    end do
    call set_photol_only()
  end select
end do

if (present(l_debug)) then
  if (l_debug) then
    if (present(i_profile_debug)) then
      l = i_profile_debug
    else
      l = 1
    end if
    write(1000+l,'(A)') 'PRESS(PA) TEMP(K) MASS(KGM-2) DENSITY(KGM-3) GASES'
    do i=1, n_layer
      write(1000+l,'(7(1pe16.8))') &
        atm%p(l, i), atm%t(l, i), atm%mass(l, i), atm%density(l, i), &
        atm%gas_mix_ratio(l, i, 1:MIN(spectrum%gas%n_absorb,3))
    end do
    write(1000+l,'(A)') 'PLEV(PA) TEMP(K)'
    do i=0, n_layer
      write(1000+l,'(2(1pe16.8))') atm%p_level(l, i), atm%t_level(l, i)
    end do
  end if
end if

contains

  subroutine set_layer_field(out_field, full_field, oned_field)
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
    end if
  end subroutine set_layer_field


  subroutine set_level_field(out_field, full_field, oned_field)
    implicit none

    real(RealK), intent(inout) :: out_field(:, 0:)
!     Output field
    real(RealExt), intent(in), optional :: full_field(:, :)
!     Full field variable
    real(RealExt), intent(in), optional :: oned_field(:)
!     One-dimensional variable

    if (present(full_field)) then
      if (l_last) then
        do i=0, n_layer
          ii = abs(level_offset-i) + 1
          do l=1, n_profile
            out_field(l, i) = real(full_field(ii, list(l)), RealK)
          end do
        end do
      else
        do i=0, n_layer
          ii = abs(level_offset-i) + 1
          do l=1, n_profile
            out_field(l, i) = real(full_field(list(l), ii), RealK)
          end do
        end do
      end if
    else if (present(oned_field)) then
      if (l_last) then
        do i=0, n_layer
          ii = abs(level_offset-i) + 1
          do l=1, n_profile
            ll = stride_level*(list(l)-1) + ii
            out_field(l, i) = real(oned_field(ll), RealK)
          end do
        end do
      else
        do i=0, n_layer
          ii = abs(level_offset-i) + 1
          do l=1, n_profile
            ll = n_profile*(ii-1) + list(l)
            out_field(l, i) = real(oned_field(ll), RealK)
          end do
        end do
      end if
    end if
  end subroutine set_level_field


  subroutine set_gas_mix_ratio(full_field, oned_field, mix_ratio, l_well_mixed)
    implicit none

    real(RealExt), intent(in), optional :: full_field(:, :)
!     Full field mass mixing ratio
    real(RealExt), intent(in), optional :: oned_field(:)
!     One-dimensional mass mixing ratio
    real(RealExt), intent(in), optional :: mix_ratio
!     Well mixed mass mixing ratio
    logical, intent(in), optional :: l_well_mixed
!     Flag to use the well mixed ratio

    integer :: i_select
    integer, parameter :: ip_zero = 0
    integer, parameter :: ip_full_field = 1
    integer, parameter :: ip_oned_field = 2
    integer, parameter :: ip_well_mixed = 3

    if (present(full_field)) then
      i_select=ip_full_field
    else if (present(oned_field)) then
      i_select=ip_oned_field
    else if (present(mix_ratio)) then
      i_select=ip_well_mixed
    else
      i_select=ip_zero
    end if
    if ( present(mix_ratio) .and. present(l_well_mixed) ) then
      if (l_well_mixed) then
        i_select=ip_well_mixed
      end if
    end if

    select case(i_select)
    case(ip_full_field)
      if (l_last) then
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            atm%gas_mix_ratio(l, i, i_gas) &
              = max(real(full_field(ii, list(l)), RealK), 0.0_RealK)
          end do
        end do
      else
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            atm%gas_mix_ratio(l, i, i_gas) &
              = max(real(full_field(list(l), ii), RealK), 0.0_RealK)
          end do
        end do
      end if
    case(ip_oned_field)
      if (l_last) then
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            ll = stride_layer*(list(l)-1) + ii
            atm%gas_mix_ratio(l, i, i_gas) &
              = max(real(oned_field(ll), RealK), 0.0_RealK)
          end do
        end do
      else
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            ll = n_profile*(ii-1) + list(l)
            atm%gas_mix_ratio(l, i, i_gas) &
              = max(real(oned_field(ll), RealK), 0.0_RealK)
          end do
        end do
      end if
    case(ip_well_mixed)
      do i=1, n_layer
        do l=1, n_profile
          atm%gas_mix_ratio(l, i, i_gas) = real(mix_ratio, RealK)
        end do
      end do
      if (mix_ratio < epsilon(1.0_RealK)) call set_photol_only()
    case(ip_zero)
      do i=1, n_layer
        do l=1, n_profile
          atm%gas_mix_ratio(l, i, i_gas) = 0.0_RealK
        end do
      end do
      call set_photol_only()
    end select
  end subroutine set_gas_mix_ratio


  subroutine set_photol_only()
    implicit none

    do i=1, spectrum%photol%n_pathway
      if (spectrum%photol%pathway_absorber(i) == i_gas) then
        ! Mixing ratio is zero but photolysis rates are required
        control%l_photol_only(spectrum%gas%type_absorb(i_gas)) = .true.
      end if
    end do
  end subroutine set_photol_only

end subroutine set_atm
end module socrates_set_atm
