! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate an observed effective radius
!
! Purpose:
!   Calculates a 2-D cloud-top diagnostic from the input 3-D version
!   of the diagnostic.
!
! Method:
!   For each type of cloud containing water in any layer the input diagnostic
!   is weighted with the product of the area of the cloud and the probability
!   that light emitted from the cloud reaches the observing instrument.
!   (Modified from the UM routine: r2_cloud_level_diag)
!
!------------------------------------------------------------------------------
module socrates_cloud_level_diag
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_CLOUD_LEVEL_DIAG'
contains

subroutine cloud_level_diag(control, dimen, atm, cld, l_all_temps, &
  list, weighted_diag, sum_weight_diag)

use def_control, only: StrCtrl
use def_dimen,   only: StrDim
use def_atm,     only: StrAtm
use def_cld,     only: StrCld

use realtype_rd, only: RealK, RealExt
use ereport_mod, only: ereport
use errormessagelength_mod, only: errormessagelength
use rad_pcf, only: i_normal, i_err_fatal, &
  ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat, ip_cloud_csiw, &
  ip_cloud_combine_homogen, ip_cloud_combine_ice_water, &
  ip_cloud_split_homogen, ip_cloud_split_ice_water, &
  ip_clcmp_st_water, ip_clcmp_cnv_water, ip_cloud_mix_random, &
  ip_cloud_type_water, ip_cloud_type_strat, ip_cloud_type_conv, &
  ip_cloud_type_sw, ip_cloud_type_si, ip_cloud_type_cw, ip_cloud_type_ci

implicit none

! Control options:
type(StrCtrl), intent(in) :: control

! Dimensions:
type(StrDim), intent(in) :: dimen

! Atmospheric properties:
type(StrAtm), intent(in) :: atm

! Cloud properties:
type(StrCld), intent(in) :: cld

! Logical flags for diagnostics
logical, intent(in) :: l_all_temps
!   If TRUE, the routine has been called to obtain diagnostics for clouds
!   consisting of liquid water at any temperature (as done in MODIS retrievals).
!   If FALSE, only clouds with temperatures above freezing are to be diagnosed
!   (as done in AVHRR retrievals).

integer, intent(in) :: list(:)
!   List of profiles to fill in output fields

real(RealExt), intent(inout), pointer :: weighted_diag(:)
!   Weighted sum of cloud-top diagnostic and weighting function
real(RealExt), intent(inout), pointer :: sum_weight_diag(:)
!   Sum of weights for cloud-top diagnostic



integer :: i, l
!   Loop variables

real(RealK) :: &
  trans_overlying_space(atm%n_profile), &
!   Probability of a photon in clear air in the level above
!   the current one reaching space
  area_exposed(atm%n_profile), &
!   Total area of cloud in the current layer exposed to
!   clear air in the layer above
  area_exposed_st(atm%n_profile), &
!   Total area of stratiform cloud in the current layer
!   exposed to clear air in the layer above
  area_exposed_cnv(atm%n_profile), &
!   Total area of convective cloud in the current layer
!   exposed to clear air in the layer above
  area_clear_above(atm%n_profile), &
!   Area of the clear sky region in the layer above
  area_strat(atm%n_profile), &
!   Area of stratiform cloud in the current layer
  area_strat_above(atm%n_profile), &
!   Area of stratiform cloud in the layer above
  area_conv(atm%n_profile), &
!   Area of convective cloud in the current layer
  area_conv_above(atm%n_profile), &
!   Area of convective cloud in the layer above
  area_clear_clear(atm%n_profile), &
!   Area of boundary where clear sky overlies clear sky
  area_clear(atm%n_profile), &
!   Area of clear sky in the current layer down to a level
  area_uncorrelated(atm%n_profile)
!   Uncorrelated region on the interface

! Indicator function
real(RealK) :: chi_cnv(atm%n_profile)
!   Convective indicator function
real(RealK) :: chi_st(atm%n_profile)
!   Stratiform indicator function

integer                      :: ierr = i_normal
character (len=*), parameter :: RoutineName = 'CLOUD_LEVEL_DIAG'
character (len=errormessagelength) :: cmessage

real(RealK), parameter :: tm = 273.15_RealK
!   Minimum temperature for warm droplets consistent with AVHRR retrievals


! Initialization
do l=1, atm%n_profile
  weighted_diag(list(l))=0.0_RealExt
  sum_weight_diag(list(l))=0.0_RealExt
end do

if (.not.control%l_cloud) then
  return
end if

! Initialize the transmision above clouds
do l=1, atm%n_profile
  trans_overlying_space(l)=1.0_RealK
  area_clear_above(l)=1.0_RealK
end do
select case (control%i_cloud_representation)
case (ip_cloud_conv_strat, ip_cloud_split_homogen, &
      ip_cloud_csiw, ip_cloud_split_ice_water)
  do l=1, atm%n_profile
    area_strat_above(l)=0.0_RealK
    area_conv_above(l)=0.0_RealK
  end do
end select

! Step down through the atmosphere calculating contributions to the diagnostics
! and subsequently allowing for transmission through the current layer.
do i=dimen%id_cloud_top, atm%n_layer

  do l=1, atm%n_profile
    area_clear(l)=1.0_RealK-cld%w_cloud(l, i)
  end do

  ! Calculate the local area of cloud radiating into clear air.
  if (control%i_cloud == ip_cloud_mix_random) then
    ! Random overlap if explicitly set
    do l=1, atm%n_profile
      area_exposed(l)=cld%w_cloud(l, i)*area_clear_above(l)
    end do
  else
    ! Otherwise assume maximum overlap for this diagnostic
    do l=1, atm%n_profile
      area_exposed(l)=max(0.0_RealK, (cld%w_cloud(l, i) &
         +area_clear_above(l)-1.0_RealK))
    end do
  end if

  select case (control%i_cloud_representation)
  case (ip_cloud_homogen, ip_cloud_combine_homogen)
    ! Only exposed cloud containing water is included
    do l = 1, atm%n_profile
      if (cld%condensed_mix_ratio(l, i, ip_clcmp_st_water) > 0.0_RealK .and. &
         ((atm%t(l, i) > tm) .or. l_all_temps)) then
        weighted_diag(list(l)) = weighted_diag(list(l)) &
          + real(trans_overlying_space(l) * area_exposed(l) &
          * cld%condensed_dim_char(l, i, ip_clcmp_st_water), RealExt)
        sum_weight_diag(list(l)) = sum_weight_diag(list(l)) &
          + real(trans_overlying_space(l) * area_exposed(l), RealExt)
      end if
    end do

  case (ip_cloud_ice_water, ip_cloud_combine_ice_water)
    ! The exposed areas include only the parts of the clouds containing
    ! water droplets.
    do l = 1, atm%n_profile
      if ((atm%t(l, i) > tm) .or. l_all_temps) then
        weighted_diag(list(l)) = weighted_diag(list(l)) &
          + real(trans_overlying_space(l) * area_exposed(l) &
          * cld%frac_cloud(l, i, ip_cloud_type_water) &
          * cld%condensed_dim_char(l, i, ip_clcmp_st_water), RealExt)
        sum_weight_diag(list(l)) = sum_weight_diag(list(l)) &
          + real(trans_overlying_space(l) * area_exposed(l) &
          * cld%frac_cloud(l, i, ip_cloud_type_water), RealExt)
      end if
    end do

  case (ip_cloud_conv_strat, ip_cloud_split_homogen)
    ! Here, the different types of clouds overlap coherently so stratiform
    ! cloud will be exposed only if there is less stratiform cloud in the
    ! layer above and more clear air in the layer above:
    ! under these conditions the non-correlated areas overlap randomly.
    do l=1, atm%n_profile
      area_strat(l)=cld%w_cloud(l, i) &
         *cld%frac_cloud(l, i, ip_cloud_type_strat)
      area_conv(l)=cld%w_cloud(l, i) &
        *cld%frac_cloud(l, i, ip_cloud_type_conv)
      area_uncorrelated(l)=1.0_RealK &
        -min(area_clear(l), area_clear_above(l)) &
        -min(area_strat(l), area_strat_above(l)) &
        -min(area_conv(l), area_conv_above(l))
      ! First find the area of uncorrelated stratiform cloud.
      area_exposed_st(l)=max(0.0_RealK, area_strat(l)-area_strat_above(l))
      area_exposed_st(l)=max(0.0_RealK, area_exposed_st(l) &
         *(area_clear_above(l)-area_clear(l)))
      ! Now normalize within the uncorrelated region.
      ! If the uncorrelated area is 0 the exposed area must be 0,
      ! so no second branch of the IF-test is required.
      if (area_uncorrelated(l) > 0.0_RealK) then
        area_exposed_st(l)=area_exposed_st(l)/area_uncorrelated(l)
      end if
      area_exposed_cnv(l)=area_exposed(l)-area_exposed_st(l)
    end do

    ! The indicator functions for liquid water in convective or straiform
    ! clouds are set to 1 if there is any liquid water and to 0 otherwise.
    do l=1, atm%n_profile
      if (cld%condensed_mix_ratio(l, i, ip_clcmp_cnv_water) > 0.0_RealK) then
        chi_cnv(l)=1.0_RealK
      else
        chi_cnv(l)=0.0_RealK
      end if
      if (cld%condensed_mix_ratio(l, i, ip_clcmp_st_water) > 0.0_RealK) then
        chi_st(l)=1.0_RealK
      else
        chi_st(l)=0.0_RealK
      end if
    end do

    ! Include contributions from convective and stratiform water clouds.
    do l=1, atm%n_profile
      if ((atm%t(l, i) > tm) .or. l_all_temps) then
        weighted_diag(list(l))=weighted_diag(list(l)) &
           +real(trans_overlying_space(l) &
           *(area_exposed_cnv(l)*chi_cnv(l) &
           *cld%condensed_dim_char(l, i, ip_clcmp_cnv_water) &
           +area_exposed_st(l)*chi_st(l) &
           *cld%condensed_dim_char(l, i, ip_clcmp_st_water)), RealExt)
        sum_weight_diag(list(l))=sum_weight_diag(list(l)) &
           +real(trans_overlying_space(l) &
           *(area_exposed_cnv(l)*chi_cnv(l) &
           +area_exposed_st(l)*chi_st(l)), RealExt)
      end if
    end do

  case (ip_cloud_csiw, ip_cloud_split_ice_water)
    ! Here, the different types of clouds overlap coherently so stratiform
    ! cloud will be exposed only if there is less stratiform cloud in the
    ! layer above and more clear air in the layer above:
    ! under these conditions the non-correlated areas overlap randomly.
    ! The actual exposed areas of convective or stratiform cloud must then
    ! be weighted by factors representing the liquid portion of each cloud,
    ! since nothing is retrieved over ice. (The horizontal arrangement of
    ! ice and water within either type of cloud is random).
    do l=1, atm%n_profile

      area_strat(l)=cld%w_cloud(l, i) &
         *(cld%frac_cloud(l, i, ip_cloud_type_sw) &
         +cld%frac_cloud(l, i, ip_cloud_type_si))
      area_conv(l)=cld%w_cloud(l, i) &
         *(cld%frac_cloud(l, i, ip_cloud_type_cw) &
         +cld%frac_cloud(l, i, ip_cloud_type_ci))
      area_uncorrelated(l)=1.0_RealK &
         -min(area_clear(l), area_clear_above(l)) &
         -min(area_strat(l), area_strat_above(l)) &
         -min(area_conv(l), area_conv_above(l))
      area_exposed_st(l)=max(0.0_RealK, area_strat(l)-area_strat_above(l))
      if (area_uncorrelated(l) > 0.0_RealK) then
        area_exposed_st(l)=max(0.0_RealK, area_exposed_st(l) &
           *(area_clear_above(l)-area_clear(l))) &
           /area_uncorrelated(l)
      else
        area_exposed_st(l)=0.0_RealK
      end if
      area_exposed_cnv(l)=area_exposed(l)-area_exposed_st(l)

      if (cld%frac_cloud(l, i, ip_cloud_type_cw) > 0.0_RealK) then
        area_exposed_cnv(l)=area_exposed_cnv(l) &
           /(1.0_RealK &
           +cld%frac_cloud(l, i, ip_cloud_type_ci) &
           /cld%frac_cloud(l, i, ip_cloud_type_cw))
      else
        area_exposed_cnv(l)=0.0_RealK
      end if

      if (cld%frac_cloud(l, i, ip_cloud_type_sw) > 0.0_RealK) then
        area_exposed_st(l)=area_exposed_st(l) &
           /(1.0_RealK &
           +cld%frac_cloud(l, i, ip_cloud_type_si) &
           /cld%frac_cloud(l, i, ip_cloud_type_sw))
      else
        area_exposed_st(l)=0.0_RealK
      end if

    end do

    do l=1, atm%n_profile
      if ((atm%t(l, i) > tm) .or. l_all_temps) then
        weighted_diag(list(l))=weighted_diag(list(l)) &
           +real(trans_overlying_space(l) &
           *(area_exposed_cnv(l) &
           *cld%condensed_dim_char(l, i, ip_clcmp_cnv_water) &
           +area_exposed_st(l) &
           *cld%condensed_dim_char(l, i, ip_clcmp_st_water)), RealExt)
        sum_weight_diag(list(l))=sum_weight_diag(list(l)) &
           +real(trans_overlying_space(l) &
           *(area_exposed_cnv(l)+area_exposed_st(l)), RealExt)
      end if
    end do

  case default
    cmessage = 'The diagnostic of observed Re has not been' &
         //' implemented with this cloud representation.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end select


  ! Advance the stored quantities refferring to overlying layers.

  ! The transmission to space currently holds the probability that
  ! a photon travelling upwards in the clear air in the layer above
  ! will escape to space without encountering a cloud. To advance this
  ! to the current layer it must be multiplied by a factor representing
  ! the overlap assumption at the top of the present layer.
  if (control%i_cloud == ip_cloud_mix_random) then
    ! Random overlap if explicitly set
    do l=1, atm%n_profile
      trans_overlying_space(l)=trans_overlying_space(l)*area_clear_above(l)
    end do
  else
    ! Otherwise assume maximum overlap for this diagnostic
    do l=1, atm%n_profile
      area_clear_clear(l)=min(area_clear(l), area_clear_above(l))
      if (area_clear(l) > 0.0_RealK) then
        trans_overlying_space(l)=trans_overlying_space(l) &
           *area_clear_clear(l)/area_clear(l)
      else
        trans_overlying_space(l)=0.0_RealK
      end if
    end do
  end if

  ! Advance the areas of cloud.
  do l=1, atm%n_profile
    area_clear_above(l)=area_clear(l)
  end do
  select case (control%i_cloud_representation)
  case (ip_cloud_conv_strat, ip_cloud_split_homogen, &
        ip_cloud_csiw, ip_cloud_split_ice_water)
    do l=1, atm%n_profile
      area_strat_above(l)=area_strat(l)
      area_conv_above(l)=area_conv(l)
    end do
  end select

end do

end subroutine cloud_level_diag
end module socrates_cloud_level_diag
