! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module holds the interface to COSP2 which is called from LFRic
!
!------------------------------------------------------------------------------

module cosp_mod

use mod_cosp_config, only: n_backscatter_bins => sr_bins, &
                           n_isccp_tau_bins => numisccptaubins, &
                           n_isccp_pressure_bins => numisccppresbins
use cosp_input_mod, only: n_cloudsat_levels

implicit none
character(len=*), parameter, private :: ModuleName = 'cosp_mod'
contains

subroutine cosp( nlevels, &
        npoints, &
        ncolumns, &
        nclds, &
        ncldy, &
        p_full_levels, &
        p_half_levels, &
        hgt_full_levels, &
        hgt_half_levels, &
        d_mass, &
        t_n, &
        q_n, &
        w_cloud, &
        condensed_mix_ratio_water, &
        condensed_mix_ratio_ice, &
        cosp_crain, &
        cosp_csnow, &
        condensed_re_water, &
        condensed_re_ice, &
        cloud_extinction, &
        cloud_absorptivity, &
        frac_cloud_water, &
        frac_cloud_ice, &
        clw_sub_full, &
        t_surf, &
        p_surf, &
        hgt_surf, &
        cosp_sunlit, &
        x1r, &
        x1g, &
        x2r, &
        x2g, &
        x4g, &
        cosp_diag, &
        profile_list, &
        l_profile_last, &
        cosp_out_ext )

  use realtype_rd, only: RealExt
  use cosp_kinds, only: wp
  use cosp2_types_mod, only: &
      cosp2_config, cosp_inputs_host_model, &
      construct_cosp_inputs_host_model, &
      construct_cosp_optical_inputs, &
      construct_cosp_column_inputs, &
      construct_cosp_outputs, &
      destroy_cosp_optical_inputs, destroy_cosp_column_inputs, &
      destroy_cosp_outputs, destroy_cosp_inputs_host_model
  use cosp_input_mod, only: &
      cosp_radar_freq, cosp_k2, &
      cosp_use_gas_abs, cosp_do_ray, &
      cosp_isccp_topheight, &
      cosp_isccp_topheight_direction, &
      cosp_surface_radar, &
      cosp_nchannels, cosp_nlr, &
      cloudsat_micro_scheme, &
      cosp_emsfc_lw, &
      cosp_use_vgrid, cosp_csat_vgrid
  use cosp2_constants_mod, only: &
      cosp_qb_dist, cosp_hydroclass_init, &
      cosp_microphys_init, &
      i_lscliq, i_lscice, &
      i_cvrain, i_cvsnow
  use cosp_radiation_mod, only: &
      cosp_radiative_properties
  use mod_cosp, only: cosp_init, &
      cosp_column_inputs, cosp_simulator, &
      cosp_optical_inputs, cosp_outputs, &
      cosp_cleanup
  use mod_cosp_config, only: r_undef
  use mod_quickbeam_optics, only: &
      quickbeam_optics_init
  use cosp_def_diag, only: cospdiag
  use cosp2_diagnostics_mod, only: create_mask

  use ereport_mod, only: ereport
  use errormessagelength_mod, only: errormessagelength
  use rad_pcf, only: i_normal, i_err_fatal

  implicit none

  ! Input arguments
  integer, intent(in) :: nlevels
  integer, intent(in) :: npoints
  integer, intent(in) :: ncolumns
  integer, intent(in) :: nclds
  integer, intent(in), pointer :: ncldy(:)
  real(RealExt), intent(in), pointer :: p_full_levels(:,:)
  real(RealExt), intent(in), pointer :: p_half_levels(:,:)
  real(RealExt), intent(in), pointer :: hgt_full_levels(:,:)
  real(RealExt), intent(in), pointer :: hgt_half_levels(:,:)
  real(RealExt), intent(in), pointer :: d_mass(:,:)
  real(RealExt), intent(in), pointer :: T_n(:,:)
  real(RealExt), intent(in), pointer :: q_n(:,:)
  real(RealExt), intent(in), pointer :: w_cloud(:,:)
  real(RealExt), intent(in), pointer :: condensed_mix_ratio_water(:,:)
  real(RealExt), intent(in), pointer :: condensed_mix_ratio_ice(:,:)
  real(RealExt), intent(in), pointer :: cosp_crain(:,:)
  real(RealExt), intent(in), pointer :: cosp_csnow(:,:)
  real(RealExt), intent(in), pointer :: condensed_re_water(:,:)
  real(RealExt), intent(in), pointer :: condensed_re_ice(:,:)
  real(RealExt), intent(in), pointer :: cloud_extinction(:,:)
  real(RealExt), intent(in), pointer :: cloud_absorptivity(:,:)
  real(RealExt), intent(in), pointer :: frac_cloud_water(:,:)
  real(RealExt), intent(in), pointer :: frac_cloud_ice(:,:)
  real(RealExt), intent(in), pointer :: clw_sub_full(:,:,:)
  real(RealExt), intent(in), pointer :: t_surf(:), p_surf(:), hgt_surf(:)
  real(RealExt), intent(in), pointer :: cosp_sunlit(:)

  ! Microphysics
  real(RealExt), intent(in) :: x1r
  real(RealExt), intent(in) :: x1g
  real(RealExt), intent(in) :: x2r
  real(RealExt), intent(in) :: x2g
  real(RealExt), intent(in) :: x4g

  ! Output arguments
  type(cospdiag), intent(inout) :: cosp_diag

  ! Optional input arguments
  integer, intent(in), optional :: profile_list(:)
  logical, intent(in), optional :: l_profile_last

  ! Optional output arguments
  type(cosp_outputs), intent(out), target, optional :: cosp_out_ext

  ! Local variables
  integer :: list(npoints)
  integer :: i
  integer :: ii
  integer :: j
  integer :: l
  integer :: lll

  logical :: l_last

  real(wp) :: x1r_wp
  real(wp) :: x1g_wp
  real(wp) :: x2r_wp
  real(wp) :: x2g_wp
  real(wp) :: x4g_wp

  real(wp) :: cosp_temp
  real(wp) :: cosp_temp1
  real(wp) :: cosp_temp2

  ! COSP configuration options
  type(cosp2_config) :: cosp_cfg
  ! COSP inputs of host model
  type(cosp_inputs_host_model) :: cosp_hmodel
  ! COSP optical inputs
  type(cosp_optical_inputs) :: cosp_optical_in
  ! COSP column inputs
  type(cosp_column_inputs) :: cosp_column_in

  type(cosp_outputs), target :: cosp_out_int
  type(cosp_outputs), pointer :: cosp_out

  character(len=256) :: cosp_status(100)

  integer :: ierr = i_normal
  character (len=errormessagelength) :: cmessage
  character (len=*), parameter :: RoutineName = 'COSP'


  if (present(cosp_out_ext)) then
    cosp_out => cosp_out_ext
  else
    cosp_out => cosp_out_int
  end if

  ! Logicals for outputs
  ! ISCCP
  cosp_cfg%lalbisccp       = .false.
  cosp_cfg%ltauisccp       = .false.
  cosp_cfg%lpctisccp       = .false.
  cosp_cfg%lcltisccp       = .false.
  cosp_cfg%lmeantbisccp    = .false.
  cosp_cfg%lmeantbclrisccp = .false.
  cosp_cfg%lclisccp        =  &
    associated(cosp_diag%cosp_ctp_tau_histogram) ! 2337
  cosp_cfg%lboxptopisccp   = .false.
  cosp_cfg%lboxtauisccp    = .false.
  ! CALIPSO
  cosp_cfg%llidarbetamol532 = .false.
  cosp_cfg%latb532          =  &
    associated(cosp_diag%cosp_calipso_tot_backscatter) ! 2341
  cosp_cfg%latb532gbx       = .false.
  cosp_cfg%lcfadlidarsr532  =  &
    associated(cosp_diag%cosp_calipso_cfad_sr_40) ! 2370
  cosp_cfg%lclcalipso       = .false.
  cosp_cfg%lcllcalipso      =  &
    associated(cosp_diag%cosp_calipso_low_level_cl) .or. & ! 2344
    associated(cosp_diag%cosp_calipso_low_level_cl_mask)
  cosp_cfg%lclmcalipso      =  &
    associated(cosp_diag%cosp_calipso_mid_level_cl) .or. & ! 2345
    associated(cosp_diag%cosp_calipso_mid_level_cl_mask)
  cosp_cfg%lclhcalipso      =  &
    associated(cosp_diag%cosp_calipso_high_level_cl) .or. & ! 2346
    associated(cosp_diag%cosp_calipso_high_level_cl_mask)
  cosp_cfg%lcltcalipso      = .false.
  cosp_cfg%lparasolrefl     = .false.
  cosp_cfg%lclcalipsoliq    =  &
    associated(cosp_diag%cosp_calipso_cf_40_liq) .or. & ! 2473
    associated(cosp_diag%cosp_calipso_cf_40_mask)
  cosp_cfg%lclcalipsoice    =  &
    associated(cosp_diag%cosp_calipso_cf_40_ice) ! 2474
  cosp_cfg%lclcalipsoun     =  &
    associated(cosp_diag%cosp_calipso_cf_40_undet) ! 2475
  cosp_cfg%lcllcalipsoliq   = .false.
  cosp_cfg%lclmcalipsoliq   = .false.
  cosp_cfg%lclhcalipsoliq   = .false.
  cosp_cfg%lcltcalipsoliq   = .false.
  cosp_cfg%lcllcalipsoice   = .false.
  cosp_cfg%lclmcalipsoice   = .false.
  cosp_cfg%lclhcalipsoice   = .false.
  cosp_cfg%lcltcalipsoice   = .false.
  cosp_cfg%lcllcalipsoun    = .false.
  cosp_cfg%lclmcalipsoun    = .false.
  cosp_cfg%lclhcalipsoun    = .false.
  cosp_cfg%lcltcalipsoun    = .false.
  cosp_cfg%lclcalipsotmp    = .false.
  cosp_cfg%lclcalipsotmpliq = .false.
  cosp_cfg%lclcalipsotmpice = .false.
  cosp_cfg%lclcalipsotmpun  = .false.
  cosp_cfg%lclopaquecalipso = .false.
  cosp_cfg%lclthincalipso   = .false.
  cosp_cfg%lclzopaquecalipso = .false.
  cosp_cfg%lclcalipsoopaque = .false.
  cosp_cfg%lclcalipsothin   = .false.
  cosp_cfg%lclcalipsozopaque = .false.
  cosp_cfg%lclcalipsoopacity = .false.
  cosp_cfg%lclopaquetemp    = .false.
  cosp_cfg%lclthintemp      = .false.
  cosp_cfg%lclzopaquetemp   = .false.
  cosp_cfg%lclopaquemeanz   = .false.
  cosp_cfg%lclthinmeanz     = .false.
  cosp_cfg%lclthinemis      = .false.
  cosp_cfg%lclopaquemeanzse = .false.
  cosp_cfg%lclthinmeanzse   = .false.
  cosp_cfg%lclzopaquecalipsose = .false.
  ! Ground lidar
  cosp_cfg%llidarbetamol532gr = .false.
  cosp_cfg%lcfadlidarsr532gr = .false.
  cosp_cfg%latb532gr = .false.
  cosp_cfg%lclgrlidar532 = .false.
  cosp_cfg%lclhgrlidar532 = .false.
  cosp_cfg%lcllgrlidar532 = .false.
  cosp_cfg%lclmgrlidar532 = .false.
  cosp_cfg%lcltgrlidar532 = .false.
  cosp_cfg%llidarbetamol355 = .false.
  cosp_cfg%lcfadlidarsr355 = .false.
  ! ATLID
  cosp_cfg%latb355 = .false.
  cosp_cfg%lclatlid = .false.
  cosp_cfg%lclhatlid = .false.
  cosp_cfg%lcllatlid = .false.
  cosp_cfg%lclmatlid = .false.
  cosp_cfg%lcltatlid = .false.
  ! CloudSat
  cosp_cfg%lcfaddbze94  = .false.
  cosp_cfg%ldbze94      = .false.
  cosp_cfg%ldbze94gbx   = .false.
  cosp_cfg%lcloudsat_tcc = .false.
  cosp_cfg%lcloudsat_tcc2 = .false.
  ! CloudSat and CALIPSO
  cosp_cfg%lclcalipso2    = .false.
  cosp_cfg%lcltlidarradar = .false.
  cosp_cfg%lcllidarradar  = .false.
  ! RTTOV
  cosp_cfg%ltbrttov = .false.
  ! MISR
  cosp_cfg%lclmisr = .false.
  ! MODIS
  cosp_cfg%lclmodis      = .false.
  cosp_cfg%lcltmodis     = .false.
  cosp_cfg%lclwmodis     = .false.
  cosp_cfg%lclimodis     = .false.
  cosp_cfg%lclhmodis     = .false.
  cosp_cfg%lclmmodis     = .false.
  cosp_cfg%lcllmodis     = .false.
  cosp_cfg%ltautmodis    = .false.
  cosp_cfg%ltauwmodis    = .false.
  cosp_cfg%ltauimodis    = .false.
  cosp_cfg%ltautlogmodis = .false.
  cosp_cfg%ltauwlogmodis = .false.
  cosp_cfg%ltauilogmodis = .false.
  cosp_cfg%lreffclwmodis = .false.
  cosp_cfg%lreffclimodis = .false.
  cosp_cfg%lpctmodis     = .false.
  cosp_cfg%llwpmodis     = .false.
  cosp_cfg%liwpmodis     = .false.
  ! CloudSat and MODIS
  cosp_cfg%lptradarflag0 = .false.
  cosp_cfg%lptradarflag1 = .false.
  cosp_cfg%lptradarflag2 = .false.
  cosp_cfg%lptradarflag3 = .false.
  cosp_cfg%lptradarflag4 = .false.
  cosp_cfg%lptradarflag5 = .false.
  cosp_cfg%lptradarflag6 = .false.
  cosp_cfg%lptradarflag7 = .false.
  cosp_cfg%lptradarflag8 = .false.
  cosp_cfg%lptradarflag9 = .false.
  cosp_cfg%lradarpia     = .false.
  cosp_cfg%lwr_occfreq   = .false.
  cosp_cfg%lcfodd        = .false.
  ! Other
  cosp_cfg%lfracout = .false.

  ! Instrument flags
  cosp_cfg%lisccp =               &
         cosp_cfg%lalbisccp       &
    .or. cosp_cfg%ltauisccp       &
    .or. cosp_cfg%lpctisccp       &
    .or. cosp_cfg%lcltisccp       &
    .or. cosp_cfg%lmeantbisccp    &
    .or. cosp_cfg%lmeantbclrisccp &
    .or. cosp_cfg%lclisccp        &
    .or. cosp_cfg%lboxptopisccp   &
    .or. cosp_cfg%lboxtauisccp
  cosp_cfg%lcalipso  =                &
         cosp_cfg%llidarbetamol532    &
    .or. cosp_cfg%latb532             &
    .or. cosp_cfg%latb532gbx          &
    .or. cosp_cfg%lcfadlidarsr532     &
    .or. cosp_cfg%lclcalipso          &
    .or. cosp_cfg%lcllcalipso         &
    .or. cosp_cfg%lclmcalipso         &
    .or. cosp_cfg%lclhcalipso         &
    .or. cosp_cfg%lcltcalipso         &
    .or. cosp_cfg%lparasolrefl        &
    .or. cosp_cfg%lclcalipsoliq       &
    .or. cosp_cfg%lclcalipsoice       &
    .or. cosp_cfg%lclcalipsoun        &
    .or. cosp_cfg%lcllcalipsoliq      &
    .or. cosp_cfg%lclmcalipsoliq      &
    .or. cosp_cfg%lclhcalipsoliq      &
    .or. cosp_cfg%lcltcalipsoliq      &
    .or. cosp_cfg%lcllcalipsoice      &
    .or. cosp_cfg%lclmcalipsoice      &
    .or. cosp_cfg%lclhcalipsoice      &
    .or. cosp_cfg%lcltcalipsoice      &
    .or. cosp_cfg%lcllcalipsoun       &
    .or. cosp_cfg%lclmcalipsoun       &
    .or. cosp_cfg%lclhcalipsoun       &
    .or. cosp_cfg%lcltcalipsoun       &
    .or. cosp_cfg%lclcalipsotmp       &
    .or. cosp_cfg%lclcalipsotmpliq    &
    .or. cosp_cfg%lclcalipsotmpice    &
    .or. cosp_cfg%lclcalipsotmpun     &
    .or. cosp_cfg%lclopaquecalipso    &
    .or. cosp_cfg%lclthincalipso      &
    .or. cosp_cfg%lclzopaquecalipso   &
    .or. cosp_cfg%lclcalipsoopaque    &
    .or. cosp_cfg%lclcalipsothin      &
    .or. cosp_cfg%lclcalipsozopaque   &
    .or. cosp_cfg%lclcalipsoopacity   &
    .or. cosp_cfg%lclopaquetemp       &
    .or. cosp_cfg%lclthintemp         &
    .or. cosp_cfg%lclzopaquetemp      &
    .or. cosp_cfg%lclopaquemeanz      &
    .or. cosp_cfg%lclthinmeanz        &
    .or. cosp_cfg%lclthinemis         &
    .or. cosp_cfg%lclopaquemeanzse    &
    .or. cosp_cfg%lclthinmeanzse      &
    .or. cosp_cfg%lclzopaquecalipsose &
    .or. cosp_cfg%lclcalipso2         &
    .or. cosp_cfg%lcltlidarradar      &
    .or. cosp_cfg%lcllidarradar
  cosp_cfg%lcloudsat =           &
         cosp_cfg%lcfaddbze94    &
    .or. cosp_cfg%ldbze94        &
    .or. cosp_cfg%ldbze94gbx     &
    .or. cosp_cfg%lcloudsat_tcc  &
    .or. cosp_cfg%lcloudsat_tcc2 &
    .or. cosp_cfg%lclcalipso2    &
    .or. cosp_cfg%lcltlidarradar &
    .or. cosp_cfg%lcllidarradar  &
    .or. cosp_cfg%lptradarflag0  &
    .or. cosp_cfg%lptradarflag1  &
    .or. cosp_cfg%lptradarflag2  &
    .or. cosp_cfg%lptradarflag3  &
    .or. cosp_cfg%lptradarflag4  &
    .or. cosp_cfg%lptradarflag5  &
    .or. cosp_cfg%lptradarflag6  &
    .or. cosp_cfg%lptradarflag7  &
    .or. cosp_cfg%lptradarflag8  &
    .or. cosp_cfg%lptradarflag9  &
    .or. cosp_cfg%lradarpia      &
    .or. cosp_cfg%lwr_occfreq    &
    .or. cosp_cfg%lcfodd
  cosp_cfg%lmisr = cosp_cfg%lclmisr
  cosp_cfg%lmodis =             &
         cosp_cfg%lclmodis      &
    .or. cosp_cfg%lcltmodis     &
    .or. cosp_cfg%lclwmodis     &
    .or. cosp_cfg%lclimodis     &
    .or. cosp_cfg%lclhmodis     &
    .or. cosp_cfg%lclmmodis     &
    .or. cosp_cfg%lcllmodis     &
    .or. cosp_cfg%ltautmodis    &
    .or. cosp_cfg%ltauwmodis    &
    .or. cosp_cfg%ltauimodis    &
    .or. cosp_cfg%ltautlogmodis &
    .or. cosp_cfg%ltauwlogmodis &
    .or. cosp_cfg%ltauilogmodis &
    .or. cosp_cfg%lreffclwmodis &
    .or. cosp_cfg%lreffclimodis &
    .or. cosp_cfg%lpctmodis     &
    .or. cosp_cfg%llwpmodis     &
    .or. cosp_cfg%liwpmodis     &
    .or. cosp_cfg%lptradarflag0 &
    .or. cosp_cfg%lptradarflag1 &
    .or. cosp_cfg%lptradarflag2 &
    .or. cosp_cfg%lptradarflag3 &
    .or. cosp_cfg%lptradarflag4 &
    .or. cosp_cfg%lptradarflag5 &
    .or. cosp_cfg%lptradarflag6 &
    .or. cosp_cfg%lptradarflag7 &
    .or. cosp_cfg%lptradarflag8 &
    .or. cosp_cfg%lptradarflag9 &
    .or. cosp_cfg%lradarpia     &
    .or. cosp_cfg%lwr_occfreq   &
    .or. cosp_cfg%lcfodd
  cosp_cfg%lrttov = cosp_cfg%ltbrttov
  cosp_cfg%lgrlidar532 =             &
         cosp_cfg%llidarbetamol532gr &
    .or. cosp_cfg%lcfadlidarsr532gr  &
    .or. cosp_cfg%latb532gr          &
    .or. cosp_cfg%lclgrlidar532      &
    .or. cosp_cfg%lclhgrlidar532     &
    .or. cosp_cfg%lcllgrlidar532     &
    .or. cosp_cfg%lclmgrlidar532     &
    .or. cosp_cfg%lcltgrlidar532     &
    .or. cosp_cfg%llidarbetamol355   &
    .or. cosp_cfg%lcfadlidarsr355
  cosp_cfg%latlid =         &
         cosp_cfg%latb355   &
    .or. cosp_cfg%lclatlid  &
    .or. cosp_cfg%lclhatlid &
    .or. cosp_cfg%lcllatlid &
    .or. cosp_cfg%lclmatlid &
    .or. cosp_cfg%lcltatlid
  cosp_cfg%lparasol = cosp_cfg%lcalipso

  cosp_use_vgrid = cosp_cfg%lcloudsat .or. cosp_cfg%lcalipso


  ! Initialize PSDs
  x1r_wp = x1r
  x1g_wp = x1g
  x2r_wp = x2r
  x2g_wp = x2g
  x4g_wp = x4g
  call cosp_microphys_init(x1r_wp,x1g_wp,x2r_wp,x2g_wp,x4g_wp)

  ! Initialize Quickbeam
  call quickbeam_optics_init()

  ! Distributional parameters for hydrometeors in radar simulator
  call cosp_hydroclass_init(cosp_qb_dist)

  ! Initialize COSP simulator
  call cosp_init(cosp_cfg%lisccp, cosp_cfg%lmodis, cosp_cfg%lmisr, &
    cosp_cfg%lcloudsat, cosp_cfg%lcalipso, cosp_cfg%lgrlidar532, &
    cosp_cfg%latlid, cosp_cfg%lparasol, cosp_cfg%lrttov, &
    cosp_radar_freq, cosp_k2, cosp_use_gas_abs, cosp_do_ray, &
    cosp_isccp_topheight, cosp_isccp_topheight_direction, cosp_surface_radar, &
    cosp_optical_in%rcfg_cloudsat, cosp_use_vgrid, cosp_csat_vgrid, &
    cosp_nlr, nlevels, cloudsat_micro_scheme)

  ! Allocate memory for COSP types
  call construct_cosp_inputs_host_model(npoints, ncolumns, &
                        nlevels, cosp_hmodel)
  call construct_cosp_optical_inputs(cosp_cfg, npoints, ncolumns, &
                        nlevels, cosp_optical_in)
  call construct_cosp_column_inputs(npoints, ncolumns, &
                        nlevels, cosp_nchannels, cosp_column_in)
  call construct_cosp_outputs(cosp_cfg, npoints, ncolumns, &
                        nlevels, cosp_nlr, cosp_nchannels, cosp_out)

  ! Populate COSP types
  cosp_optical_in%emsfc_lw = cosp_emsfc_lw

  ! Profile list
  if (present(profile_list)) then
    list = profile_list(1:npoints)
  else
    do l=1, npoints
      list(l) = l
    end do
  end if

  ! Position of profile index in input arrays
  if (present(l_profile_last)) then
    l_last = l_profile_last
  else
    l_last = .false.
  end if

  if (l_last) then
    do i=1, nlevels
      ii = nlevels-i+1
        do l=1, npoints
          cosp_column_in%qv(l,i) = q_n(ii,list(l))
          cosp_column_in%pfull(l,i) = p_full_levels(ii,list(l))
          cosp_column_in%phalf(l,i+1) = p_half_levels(ii,list(l))
          cosp_column_in%hgt_matrix(l,i) = hgt_full_levels(ii,list(l))
          cosp_column_in%hgt_matrix_half(l,i+1) = hgt_half_levels(ii,list(l))
          cosp_column_in%at(l,i) = t_n(ii,list(l))
          cosp_hmodel%mr_gbx(l,i,i_cvrain) = cosp_crain(ii,list(l))
          cosp_hmodel%mr_gbx(l,i,i_cvsnow) = cosp_csnow(ii,list(l))
        end do ! l
    end do ! ii
  else
    do i=1, nlevels
      ii = nlevels-i+1
        do l=1, npoints
          cosp_column_in%qv(l,i) = q_n(list(l),ii)
          cosp_column_in%pfull(l,i) = p_full_levels(list(l),ii)
          cosp_column_in%phalf(l,i+1) = p_half_levels(list(l),ii)
          cosp_column_in%hgt_matrix(l,i) = hgt_full_levels(list(l),ii)
          cosp_column_in%hgt_matrix_half(l,i+1) = hgt_half_levels(list(l),ii)
          cosp_column_in%at(l,i) = t_n(list(l),ii)
          cosp_hmodel%mr_gbx(l,i,i_cvrain) = cosp_crain(list(l),ii)
          cosp_hmodel%mr_gbx(l,i,i_cvsnow) = cosp_csnow(list(l),ii)
        end do ! l
    end do ! ii
  end if ! l_last

  do l=1, npoints
    cosp_column_in%phalf(l,1) = 0.0_wp
    cosp_column_in%phalf(l,nlevels+1) = p_surf(list(l))
    cosp_column_in%surfelev(l) = hgt_surf(list(l))
    cosp_column_in%hgt_matrix_half(l,nlevels+1) = cosp_column_in%surfelev(l)
    cosp_column_in%hgt_matrix_half(l,1) = cosp_column_in%hgt_matrix(l,1) + &
      cosp_column_in%hgt_matrix(l,1) - cosp_column_in%hgt_matrix_half(l,2)
    cosp_column_in%skt(l) = t_surf(list(l))
    cosp_column_in%land(l) = 0.0_wp
  end do ! l

  ! Populate COSP types from SOCRATES radiation LW
  cosp_optical_in%ncolumns = ncolumns

  if (l_last) then
    do ii=1, nclds
      i = nlevels-ii+1
      do l=1, npoints
        cosp_temp1 = &
          condensed_mix_ratio_water(ii,list(l))*frac_cloud_water(ii,list(l))
        cosp_temp2 = &
          condensed_mix_ratio_ice(ii,list(l))*frac_cloud_ice(ii,list(l))
        do lll=1, ncldy(list(l))
          ! Cloud water SUBGRID mixing ratio (LIQUID)
          cosp_hmodel%mr_hydro(l,lll,i,i_lscliq) = &
            clw_sub_full(ii,lll,list(l))*cosp_temp1
          ! cloud water subgrid mixing ratio (ice)
          cosp_hmodel%mr_hydro(l,lll,i,i_lscice) = &
            clw_sub_full(ii,lll,list(l))*cosp_temp2
          ! cloud water effective radius (liquid)
          cosp_hmodel%reff_hydro(l,lll,i,i_lscliq) = &
            condensed_re_water(ii,list(l))
          ! cloud water effective dimension (ice)
          cosp_hmodel%reff_hydro(l,lll,i,i_lscice) = &
            condensed_re_ice(ii,list(l))
        end do ! lll
        ! total cloud gridbox fraction seen by radiation
        cosp_hmodel%tca_gbx(l,i) = w_cloud(ii,list(l))
        ! gridbox effective radii
        ! cloud water effective radius (liquid)
        cosp_hmodel%reff_gbx(l,i,i_lscliq) = &
          condensed_re_water(ii,list(l))
        ! cloud water effective dimension (ice)
        cosp_hmodel%reff_gbx(l,i,i_lscice) = &
          condensed_re_ice(ii,list(l))
      end do ! l
    end do ! ii
  else
    do ii=1, nclds
      i = nlevels-ii+1
      do l=1, npoints
        cosp_temp1 = &
          condensed_mix_ratio_water(list(l),ii)*frac_cloud_water(list(l),ii)
        cosp_temp2 = &
          condensed_mix_ratio_ice(list(l),ii)*frac_cloud_ice(list(l),ii)
        do lll=1, ncldy(list(l))
          ! Cloud water SUBGRID mixing ratio (LIQUID)
          cosp_hmodel%mr_hydro(l,lll,i,i_lscliq) = &
            clw_sub_full(list(l),ii,lll)*cosp_temp1
          ! cloud water subgrid mixing ratio (ice)
          cosp_hmodel%mr_hydro(l,lll,i,i_lscice) = &
            clw_sub_full(list(l),ii,lll)*cosp_temp2
          ! cloud water effective radius (liquid)
          cosp_hmodel%reff_hydro(l,lll,i,i_lscliq) = &
            condensed_re_water(list(l),ii)
          ! cloud water effective dimension (ice)
          cosp_hmodel%reff_hydro(l,lll,i,i_lscice) = &
            condensed_re_ice(list(l),ii)
        end do ! l
        ! total cloud gridbox fraction seen by radiation
        cosp_hmodel%tca_gbx(l,i) = w_cloud(list(l),ii)
        ! gridbox effective radii
        ! cloud water effective radius (liquid)
        cosp_hmodel%reff_gbx(l,i,i_lscliq) = &
          condensed_re_water(list(l),ii)
        ! cloud water effective dimension (ice)
        cosp_hmodel%reff_gbx(l,i,i_lscice) = &
          condensed_re_ice(list(l),ii)
      end do ! l
    end do ! ii
  end if ! l_last

  ! Populate COSP types from SOCRATES radiation SW and LW
  if (cosp_cfg%lmodis .or. cosp_cfg%lmisr .or. cosp_cfg%lisccp) then
    if (l_last) then
      do ii=1, nclds
        i = nlevels-ii+1
        do l= 1, npoints
          ! sub-grid cloud optical depth and emissivity for cosp assumes
          ! same generated sub-columns for liquid and ice.
          cosp_temp = cloud_extinction(ii,list(l)) * d_mass(ii,list(l))
          cosp_temp1 = cloud_absorptivity(ii,list(l)) * d_mass(ii,list(l))
          do lll=1, ncldy(list(l))
            cosp_optical_in%tau_067(l,lll,i) = &
              clw_sub_full(ii,lll,list(l))*cosp_temp
            cosp_optical_in%emiss_11(l,lll,i) = &
              1.0_wp-exp(-1.666_wp*cosp_temp1*clw_sub_full(ii,lll,list(l)))
          end do ! lll
        end do ! l
      end do ! ii
    else
      do ii=1, nclds
        i = nlevels-ii+1
        do l= 1, npoints
          ! sub-grid cloud optical depth and emissivity for cosp assumes
          ! same generated sub-columns for liquid and ice.
          cosp_temp = cloud_extinction(list(l),ii) * d_mass(list(l),ii)
          cosp_temp1 = cloud_absorptivity(list(l),ii) * d_mass(list(l),ii)
          do lll=1, ncldy(list(l))
            cosp_optical_in%tau_067(l,lll,i) = &
              clw_sub_full(list(l),ii,lll)*cosp_temp
            cosp_optical_in%emiss_11(l,lll,i) = &
              1.0_wp-exp(-1.666_wp*cosp_temp1*clw_sub_full(list(l),ii,lll))
          end do ! lll
        end do ! l
      end do ! ii
    end if ! l_last
  end if

  do l=1, npoints
    cosp_column_in%sunlit(l) = nint(cosp_sunlit(list(l)))
  end do

  ! Calculation of radiative properties for each instrument
  call cosp_radiative_properties(cosp_cfg,cosp_hmodel,cosp_column_in, &
    cosp_optical_in, cosp_qb_dist)

  ! Call to COSP
  cosp_status = cosp_simulator(cosp_optical_in, cosp_column_in, cosp_out, &
                               1, npoints, .false.)

  do i=1,size(cosp_status,1)
    if (cosp_status(i) /= '') then
      write(cmessage,'(i3,1x,a)') i, trim(cosp_status(i))
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    end if
  end do


! COSP: MASK FOR CALIPSO LOW-LEVEL CF (was 2321)
  if (associated(cosp_diag%cosp_calipso_low_level_cl_mask)) then
    do l=1, npoints
      call create_mask(r_undef, cosp_out%calipso_cldlayer(l,1), &
      cosp_diag%cosp_calipso_low_level_cl_mask(list(l)))
    end do
  end if

! COSP: CALIPSO LOW-LEVEL CLOUD (was 2344)
  if (associated(cosp_diag%cosp_calipso_low_level_cl)) then
    do l=1, npoints
      cosp_diag%cosp_calipso_low_level_cl(list(l)) &
        = 0.01_RealExt * cosp_out%calipso_cldlayer(l,1)
    end do
  end if

! COSP: MASK FOR CALIPSO MID-LEVEL CF (was 2322)
  if (associated(cosp_diag%cosp_calipso_mid_level_cl_mask)) then
    do l=1, npoints
      call create_mask(r_undef, cosp_out%calipso_cldlayer(l,2), &
      cosp_diag%cosp_calipso_mid_level_cl_mask(list(l)))
    end do
  end if

! COSP: CALIPSO MID-LEVEL CLOUD (was 2345)
  if (associated(cosp_diag%cosp_calipso_mid_level_cl)) then
    do l=1, npoints
      cosp_diag%cosp_calipso_mid_level_cl(list(l)) &
        = 0.01_RealExt * cosp_out%calipso_cldlayer(l,2)
    end do
  end if

! COSP: MASK FOR CALIPSO HIGH-LEVEL CF (was 2323)
  if (associated(cosp_diag%cosp_calipso_high_level_cl_mask)) then
    do l=1, npoints
      call create_mask(r_undef, cosp_out%calipso_cldlayer(l,3), &
      cosp_diag%cosp_calipso_high_level_cl_mask(list(l)))
    end do
  end if

! COSP: CALIPSO HIGH-LEVEL CLOUD (was 2346)
  if (associated(cosp_diag%cosp_calipso_high_level_cl)) then
    do l=1, npoints
      cosp_diag%cosp_calipso_high_level_cl(list(l)) &
        = 0.01_RealExt * cosp_out%calipso_cldlayer(l,3)
    end do
  end if

! COSP: ISCCP CTP-TAU HISTOGRAM (was 2330)
  if (associated(cosp_diag%cosp_cloud_weights)) then
    do l=1, npoints
      cosp_diag%cosp_cloud_weights(list(l)) &
        = real(cosp_column_in%sunlit(l), RealExt)
    end do
  end if

! COSP: ISCCP CTP-TAU HISTOGRAM (was 2337)
  if (associated(cosp_diag%cosp_ctp_tau_histogram)) then
    if (l_last) then
      do i=1, n_isccp_pressure_bins
        do j=1, n_isccp_tau_bins
          do l=1, npoints
            if (cosp_out%isccp_fq(l,j,i) == r_undef) then
               cosp_diag%cosp_ctp_tau_histogram(j,i,list(l)) = 0.0_RealExt
            else
               cosp_diag%cosp_ctp_tau_histogram(j,i,list(l)) &
               = 0.01_RealExt * cosp_out%isccp_fq(l,j,i)
            end if
          end do
        end do
      end do
    else
      do i=1, n_isccp_pressure_bins
        do j=1, n_isccp_tau_bins
          do l=1, npoints
            if (cosp_out%isccp_fq(l,j,i) == r_undef) then
               cosp_diag%cosp_ctp_tau_histogram(list(l),j,i) = 0.0_RealExt
            else
               cosp_diag%cosp_ctp_tau_histogram(list(l),j,i) &
               = 0.01_RealExt * cosp_out%isccp_fq(l,j,i)
            end if
          end do
        end do
      end do
    end if
  end if

! COSP: CALIPSO TOTAL BACKSCATTER (was 2341)
  if (associated(cosp_diag%cosp_calipso_tot_backscatter)) then
    if (l_last) then
      do i=1, nlevels
        ii=nlevels-i+1
        do lll=1, ncolumns
          do l=1, npoints
            cosp_diag%cosp_calipso_tot_backscatter(ii, lll, list(l)) &
              = cosp_out%calipso_beta_tot(l, lll, i)
          end do
        end do
      end do
      ! Fill uncalculated sub-columns with zeros
      do l=1, npoints
        do lll=ncolumns+1, ubound(cosp_diag%cosp_calipso_tot_backscatter, 2)
          do ii=1, nlevels
            cosp_diag%cosp_calipso_tot_backscatter(ii, lll, list(l)) &
              = 0.0_RealExt
          end do
        end do
      end do
    else
      do i=1, nlevels
        ii=nlevels-i+1
        do lll=1, ncolumns
          do l=1, npoints
            cosp_diag%cosp_calipso_tot_backscatter(list(l), ii, lll) &
              = cosp_out%calipso_beta_tot(l, lll, i)
          end do
        end do
      end do
      do lll=ncolumns+1, ubound(cosp_diag%cosp_calipso_tot_backscatter, 3)
        do ii=1, nlevels
          do l=1, npoints
            cosp_diag%cosp_calipso_tot_backscatter(list(l), ii, lll) &
              = 0.0_RealExt
          end do
        end do
      end do
    end if
  end if

! COSP CALIPSO CFAD SR 40 CSAT LEVELS (was 2370)
  if (associated(cosp_diag%cosp_calipso_cfad_sr_40)) then
    if (l_last) then
      do i=1, cosp_nlr
        do j=1, n_backscatter_bins
          do l=1, npoints
            cosp_diag%cosp_calipso_cfad_sr_40(j, i, list(l)) &
            = cosp_out%calipso_cfad_sr(l, j, i)
          end do
        end do
      end do
    else
      do i=1, cosp_nlr
        do j=1, n_backscatter_bins
          do l=1, npoints
            cosp_diag%cosp_calipso_cfad_sr_40(list(l), j, i) &
            = cosp_out%calipso_cfad_sr(l, j, i)
          end do
        end do
      end do
    end if
  end if

! COSP: MASK FOR CALIPSO CF 40 LVLS
  if (associated(cosp_diag%cosp_calipso_cf_40_mask)) then
    ! create_mask needs to be called for all cases to set r_undef to zero
    if (associated(cosp_diag%cosp_calipso_cf_40_ice)) then
      if (l_last) then
        do i=1, cosp_nlr
          do l=1, npoints
            call create_mask(r_undef, cosp_out%calipso_lidarcldphase(l, i, 1), &
            cosp_diag%cosp_calipso_cf_40_mask(i, list(l)))
          end do
        end do
      else
        do i=1, cosp_nlr
          do l=1, npoints
            call create_mask(r_undef, cosp_out%calipso_lidarcldphase(l, i, 1), &
            cosp_diag%cosp_calipso_cf_40_mask(list(l), i))
          end do
        end do
      end if
    end if
    if (associated(cosp_diag%cosp_calipso_cf_40_undet)) then
      if (l_last) then
        do i=1, cosp_nlr
          do l=1, npoints
            call create_mask(r_undef, cosp_out%calipso_lidarcldphase(l, i, 3), &
            cosp_diag%cosp_calipso_cf_40_mask(i, list(l)))
          end do
        end do
      else
        do i=1, cosp_nlr
          do l=1, npoints
            call create_mask(r_undef, cosp_out%calipso_lidarcldphase(l, i, 3), &
            cosp_diag%cosp_calipso_cf_40_mask(list(l), i))
          end do
        end do
      end if
    end if
    ! The liquid phase will always be calculated if the mask is requested
    if (l_last) then
      do i=1, cosp_nlr
        do l=1, npoints
          call create_mask(r_undef, cosp_out%calipso_lidarcldphase(l, i, 2), &
          cosp_diag%cosp_calipso_cf_40_mask(i, list(l)))
        end do
      end do
    else
      do i=1, cosp_nlr
        do l=1, npoints
          call create_mask(r_undef, cosp_out%calipso_lidarcldphase(l, i, 2), &
          cosp_diag%cosp_calipso_cf_40_mask(list(l), i))
        end do
      end do
    end if
  end if

! COSP: CALIPSO CF 40 LVLS (LIQ) (was 2473)
  if (associated(cosp_diag%cosp_calipso_cf_40_liq)) then
    if (l_last) then
      do i=1, cosp_nlr
        do l=1, npoints
          cosp_diag%cosp_calipso_cf_40_liq(i, list(l)) &
            = 0.01_RealExt * cosp_out%calipso_lidarcldphase(l, i, 2)
        end do
      end do
    else
      do i=1, cosp_nlr
        do l=1, npoints
          cosp_diag%cosp_calipso_cf_40_liq(list(l), i) &
            = 0.01_RealExt * cosp_out%calipso_lidarcldphase(l, i, 2)
        end do
      end do
    end if
  end if

! COSP: CALIPSO CF 40 LVLS (ICE) (was 2474)
  if (associated(cosp_diag%cosp_calipso_cf_40_ice)) then
    if (l_last) then
      do i=1, cosp_nlr
        do l=1, npoints
          cosp_diag%cosp_calipso_cf_40_ice(i, list(l)) &
            = 0.01_RealExt * cosp_out%calipso_lidarcldphase(l, i, 1)
        end do
      end do
    else
      do i=1, cosp_nlr
        do l=1, npoints
          cosp_diag%cosp_calipso_cf_40_ice(list(l), i) &
            = 0.01_RealExt * cosp_out%calipso_lidarcldphase(l, i, 1)
        end do
      end do
    end if
  end if

! COSP: CALIPSO CF 40 LVLS (UNDET) (was 2475)
  if (associated(cosp_diag%cosp_calipso_cf_40_undet)) then
    if (l_last) then
      do i=1, cosp_nlr
        do l=1, npoints
          cosp_diag%cosp_calipso_cf_40_undet(i, list(l)) &
            = 0.01_RealExt * cosp_out%calipso_lidarcldphase(l, i, 3)
        end do
      end do
    else
      do i=1, cosp_nlr
        do l=1, npoints
          cosp_diag%cosp_calipso_cf_40_undet(list(l), i) &
            = 0.01_RealExt * cosp_out%calipso_lidarcldphase(l, i, 3)
        end do
      end do
    end if
  end if

  if (.not.present(cosp_out_ext)) then
    call destroy_cosp_outputs(cosp_out)
    call cosp_cleanup
  end if
  call destroy_cosp_column_inputs(cosp_column_in)
  call destroy_cosp_optical_inputs(cosp_optical_in)
  call destroy_cosp_inputs_host_model(cosp_hmodel)

end subroutine cosp

end module cosp_mod
