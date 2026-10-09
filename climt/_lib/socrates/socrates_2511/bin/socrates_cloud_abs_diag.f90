! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate cloud absorptivity averaged using the local temperature,
! or for a specified wavelength
!
!------------------------------------------------------------------------------
module socrates_cloud_abs_diag
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_CLOUD_ABS_DIAG'
contains

subroutine cloud_abs_diag(cloud_thermal_absorptivity, &
  cloud_absorptivity_wavelength, &
  control, dimen, spectrum, atm, cld, &
  list, layer_offset, l_last)

use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use def_atm,      only: StrAtm
use def_cld,      only: StrCld

use realtype_rd, only: RealK, RealExt
use rad_pcf, only: i_normal, i_err_fatal, ip_cloud_off, &
  ip_clcmp_st_water, ip_clcmp_cnv_water, ip_clcmp_st_ice, ip_clcmp_cnv_ice, &
  ip_drop_pade_2, ip_ice_pade_2_phf, ip_ice_fu_phf, ip_ice_baran
use ereport_mod, only: ereport
use errormessagelength_mod, only: errormessagelength
use set_cloud_geometry_mod, only: set_cloud_geometry
use planck_flux_band_mod, only: planck_flux_band
use opt_prop_pade_2_mod, only: opt_prop_pade_2
use opt_prop_fu_phf_mod, only: opt_prop_fu_phf
use opt_prop_baran_mod, only: opt_prop_baran

implicit none

! Cloud absorptivity averaged using the local Planckian in each band, or
! for the band containing a particular wavelength if a value > 0 is specified:
real(RealExt), intent(inout), pointer :: cloud_thermal_absorptivity(:, :)
real(RealExt), intent(in) :: cloud_absorptivity_wavelength

! Control options:
type(StrCtrl), intent(in) :: control

! Dimensions:
type(StrDim), intent(in) :: dimen

! Spectral data:
type (StrSpecData), intent(in) :: spectrum

! Atmospheric properties:
type(StrAtm), intent(in) :: atm

! Cloud properties:
type(StrCld), intent(in) :: cld

integer, intent(in) :: list(:)
!   List of profiles to fill in output fields
integer, intent(in) :: layer_offset
!   Offset to loop counter to allow indexing in inverted order
logical, intent(in) :: l_last
!   Flag to loop over profiles last in diagnostics

integer :: i_band, i_band_exclude, index_exclude, k, i, ii, l, ll
!   Loop variables
integer :: n_cloud_top
!   Topmost cloudy layer
integer :: n_cloud_profile(dimen%id_cloud_top:dimen%nd_layer)
!   Number of cloudy profiles
integer :: i_cloud_profile(dimen%nd_profile, dimen%id_cloud_top:dimen%nd_layer)
!   Profiles containing clouds
real(RealK) :: planck_sum(dimen%nd_profile, &
                          dimen%id_cloud_top:dimen%nd_layer)
!   Sum of Planckian over bands for local temperature
real(RealK) :: planck_fraction(dimen%nd_profile, &
                               dimen%id_cloud_top:dimen%nd_layer, &
                               control%first_band:control%last_band)
!   Fraction of Planckian in band for local temperature
real(RealK) :: k_ext_tot_cloud_comp(dimen%nd_profile, &
                                    dimen%id_cloud_top:dimen%nd_layer)
!   Total extinction for cloud component
real(RealK) :: k_ext_scat_cloud_comp(dimen%nd_profile, &
                                     dimen%id_cloud_top:dimen%nd_layer)
!   Scattering extinction for cloud component

real(RealK), parameter :: notzero = tiny(1.0_RealK)

logical :: l_found_band

integer :: ierr = i_normal
character (len=errormessagelength) :: cmessage
character (len=*), parameter :: RoutineName = 'CLOUD_ABS_DIAG'


! Only run if cloud is on and either the Planck flux is parametrised
! or a wavelength has been specified
if ( control%i_cloud_representation /= ip_cloud_off .and. &
     ( spectrum%basic%l_present(6) .or. &
       cloud_absorptivity_wavelength > 0.0_RealExt ) ) then

  call set_cloud_geometry(atm%n_profile, atm%n_layer, &
    .false., cld%w_cloud, &
    n_cloud_top, n_cloud_profile, i_cloud_profile, &
    dimen%nd_profile, dimen%nd_layer, dimen%id_cloud_top)

  if ( cloud_absorptivity_wavelength > 0.0_RealExt ) then
    l_found_band = .false.
  else
    do i=n_cloud_top, atm%n_layer
      do ll=1, n_cloud_profile(i)
        l=i_cloud_profile(ll, i)
        planck_fraction(l, i, control%first_band) &
          = planck_flux_band(spectrum, control%first_band, atm%t(l, i))
        planck_sum(l, i) = planck_fraction(l, i, control%first_band)
      end do
    end do
    do i_band=control%first_band+1, control%last_band
      do i=n_cloud_top, atm%n_layer
        do ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          planck_fraction(l, i, i_band) &
            = planck_flux_band(spectrum, i_band, atm%t(l, i))
          planck_sum(l, i) = planck_sum(l, i) &
            + planck_fraction(l, i, i_band)
        end do
      end do
    end do
  end if
  outer : do i_band=control%first_band, control%last_band
    if ( cloud_absorptivity_wavelength > 0.0_RealExt ) then
      if ( spectrum%basic%wavelength_short(i_band) > &
           cloud_absorptivity_wavelength ) cycle outer
      if ( spectrum%basic%wavelength_long(i_band) <= &
           cloud_absorptivity_wavelength ) cycle outer
      inner : do i_band_exclude = 1, spectrum%basic%n_band_exclude(i_band)
        index_exclude = spectrum%basic%index_exclude(i_band_exclude, i_band)
        if ( spectrum%basic%wavelength_short(index_exclude) <= &
             cloud_absorptivity_wavelength .and. &
             spectrum%basic%wavelength_long(index_exclude) > &
             cloud_absorptivity_wavelength ) cycle outer
      end do inner
      if (l_found_band) then
        cmessage = ' Error locating band'//&
                   ' for cloud_thermal_absorptivity diagnostic.'
        ierr=i_err_fatal
        call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
      end if
      l_found_band = .true.
      do i=n_cloud_top, atm%n_layer
        do ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          planck_fraction(l, i, i_band) = 1.0_RealK
        end do
      end do
    else
      do i=n_cloud_top, atm%n_layer
        do ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          planck_fraction(l, i, i_band) &
            = planck_fraction(l, i, i_band) / max(planck_sum(l, i), notzero)
        end do
      end do
    end if
    do k=1, cld%n_condensed
      select case (cld%type_condensed(k))
      case (ip_clcmp_st_water, ip_clcmp_cnv_water)
        select case (cld%i_condensed_param(k))
        case (ip_drop_pade_2)
          call opt_prop_pade_2(dimen%id_cloud_top, &
            n_cloud_top, atm%n_layer, n_cloud_profile, i_cloud_profile, &
            cld%condensed_param_list(:, k, i_band), &
            cld%condensed_mix_ratio(:, :, k), &
            cld%condensed_dim_char(:, :, k), &
            k_ext_tot_cloud_comp, k_ext_scat_cloud_comp)
        case default
          cmessage = ' Droplet parametrisation not implemented'//&
                     ' for cloud_thermal_absorptivity diagnostic.'
          ierr=i_err_fatal
          call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
        end select
      case (ip_clcmp_st_ice, ip_clcmp_cnv_ice)
        select case (cld%i_condensed_param(k))
        case (ip_ice_pade_2_phf)
          call opt_prop_pade_2(dimen%id_cloud_top, &
            n_cloud_top, atm%n_layer, n_cloud_profile, i_cloud_profile, &
            cld%condensed_param_list(:, k, i_band), &
            cld%condensed_mix_ratio(:, :, k), &
            cld%condensed_dim_char(:, :, k), &
            k_ext_tot_cloud_comp, k_ext_scat_cloud_comp)
        case (ip_ice_fu_phf)
          call opt_prop_fu_phf(dimen%id_cloud_top, &
            n_cloud_top, atm%n_layer, n_cloud_profile, i_cloud_profile, &
            cld%condensed_param_list(:, k, i_band), &
            cld%condensed_mix_ratio(:, :, k), &
            cld%condensed_dim_char(:, :, k), &
            k_ext_tot_cloud_comp, k_ext_scat_cloud_comp)
        case (ip_ice_baran)
          call opt_prop_baran(dimen%id_cloud_top, &
            n_cloud_top, atm%n_layer, n_cloud_profile, i_cloud_profile, &
            cld%condensed_param_list(:, k, i_band), &
            cld%condensed_mix_ratio(:, :, k), &
            cld%condensed_dim_char(:, :, k), atm%t, &
            k_ext_tot_cloud_comp, k_ext_scat_cloud_comp)
        case default
          cmessage = ' Ice-crystal parametrisation not implemented'//&
                     ' for cloud_thermal_absorptivity diagnostic.'
          ierr=i_err_fatal
          call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
        end select
      end select
      if (l_last) then
        do i=n_cloud_top, atm%n_layer
          ii = abs(layer_offset-i)
          do ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            cloud_thermal_absorptivity(ii, list(l)) &
              = cloud_thermal_absorptivity(ii, list(l)) + real( &
              ( k_ext_tot_cloud_comp(l, i) &
              - k_ext_scat_cloud_comp(l, i) ) &
              * cld%frac_cloud(l, i, cld%i_cloud_type(k)) &
              * planck_fraction(l, i, i_band), RealExt )
          end do
        end do
      else
        do i=n_cloud_top, atm%n_layer
          ii = abs(layer_offset-i)
          do ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            cloud_thermal_absorptivity(list(l), ii) &
              = cloud_thermal_absorptivity(list(l), ii) + real( &
              ( k_ext_tot_cloud_comp(l, i) &
              - k_ext_scat_cloud_comp(l, i) ) &
              * cld%frac_cloud(l, i, cld%i_cloud_type(k)) &
              * planck_fraction(l, i, i_band), RealExt )
          end do
        end do
      end if
    end do
  end do outer

  if ( cloud_absorptivity_wavelength > 0.0_RealExt .and. &
       .not. l_found_band ) then
    cmessage = ' Wavelength outside spectral file range'//&
                 ' for cloud_thermal_absorptivity diagnostic.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end if

end if

end subroutine cloud_abs_diag
end module socrates_cloud_abs_diag
