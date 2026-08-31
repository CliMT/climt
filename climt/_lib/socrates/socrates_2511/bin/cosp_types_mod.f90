! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module holds derived types that are used for the COSP2 code
! configuration control structure as well as routines for the allocation
! and deallocation of input and output arrays.
!
!------------------------------------------------------------------------------

module cosp2_types_mod
use yomhook, only: lhook, dr_hook
use parkind1, only: jpim, jprb
use cosp_kinds, only: wp
use mod_cosp, only: cosp_optical_inputs, cosp_column_inputs, cosp_outputs
use mod_cosp_config, only: numisccptaubins, numisccppresbins,                  &
    nummisrtaubins, nummisrhgtbins, nummodistaubins, nummodispresbins,         &
    nummodistaubins, nummodisreffliqbins, nummodistaubins,                     &
    nummodisrefficebins, sr_bins, lidar_ncat, lidar_ntemp, lidar_ntype,        &
    parasol_nrefl, cloudsat_dbze_bins, wr_nregime, cfodd_ndbze,                &
    cfodd_nicod, cfodd_nclass, n_hydro
implicit none

character(len=*), private, parameter :: ModuleName='COSP_TYPES_MOD'

interface deallocate_and_nullify
  module procedure deallocate_and_nullify_1d, deallocate_and_nullify_2d,       &
                   deallocate_and_nullify_3d, deallocate_and_nullify_4d
end interface
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------- DERIVED TYPES ----------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Configuration choices
type :: cosp2_config
   ! Simulators
  logical :: lcloudsat,lcalipso,lisccp,lmodis,lmisr,lrttov,lgrlidar532,        &
          latlid,lparasol
  ! Output variables
  ! ISCCP
  logical :: lalbisccp, ltauisccp, lpctisccp, lcltisccp, lmeantbisccp,         &
             lmeantbclrisccp, lclisccp, lboxptopisccp, lboxtauisccp
  !   CALIPSO
  logical :: llidarbetamol532, latb532, latb532gbx, lcfadlidarsr532,           &
             lclcalipso, lcllcalipso, lclmcalipso, lclhcalipso, lcltcalipso,   &
             lparasolrefl, lclcalipsoliq, lclcalipsoice, lclcalipsoun,         &
             lcllcalipsoliq, lclmcalipsoliq, lclhcalipsoliq, lcltcalipsoliq,   &
             lcllcalipsoice, lclmcalipsoice, lclhcalipsoice, lcltcalipsoice,   &
             lcllcalipsoun, lclmcalipsoun, lclhcalipsoun, lcltcalipsoun,       &
             lclcalipsotmp, lclcalipsotmpliq, lclcalipsotmpice,                &
             lclcalipsotmpun, lclopaquecalipso, lclthincalipso,                &
             lclzopaquecalipso, lclcalipsoopaque, lclcalipsothin,              &
             lclcalipsozopaque, lclcalipsoopacity, lclopaquetemp,              &
             lclthintemp, lclzopaquetemp, lclopaquemeanz, lclthinmeanz,        &
             lclthinemis, lclopaquemeanzse, lclthinmeanzse,                    &
             lclzopaquecalipsose
  !   Ground lidar
  logical :: llidarbetamol532gr, lcfadlidarsr532gr, latb532gr,                 &
             lclgrlidar532, lclhgrlidar532, lcllgrlidar532, lclmgrlidar532,    &
             lcltgrlidar532, llidarbetamol355, lcfadlidarsr355
  !   ATLID
  logical :: latb355, lclatlid, lclhatlid, lcllatlid, lclmatlid, lcltatlid
  !   CloudSat
  logical :: lcfaddbze94, ldbze94, ldbze94gbx, lcloudsat_tcc, lcloudsat_tcc2
  !   CloudSat and CALIPSO
  logical :: lclcalipso2, lcltlidarradar, lcllidarradar
  !   RTTOV
  logical :: ltbrttov
  !   MISR
  logical :: lclmisr
  !   MODIS
  logical :: lclmodis, lcltmodis, lclwmodis, lclimodis, lclhmodis,             &
             lclmmodis, lcllmodis, ltautmodis, ltauwmodis, ltauimodis,         &
             ltautlogmodis, ltauwlogmodis, ltauilogmodis, lreffclwmodis,       &
             lreffclimodis, lpctmodis, llwpmodis, liwpmodis
  !   CloudSat and MODIS
  logical :: lptradarflag0, lptradarflag1, lptradarflag2, lptradarflag3,       &
             lptradarflag4, lptradarflag5, lptradarflag6, lptradarflag7,       &
             lptradarflag8, lptradarflag9, lradarpia, lwr_occfreq, lcfodd
  !   Other
  logical :: lfracout
end type cosp2_config

! Specific host-model inputs to optical calculations
type :: cosp_inputs_host_model
  integer :: npoints
  integer :: ncolumns
  integer :: nlevels
  ! points, columns, levels, hydrometeors
  real(wp), allocatable :: mr_hydro(:,:,:,:)
  real(wp), allocatable :: reff_hydro(:,:,:,:)
  ! points, levels, hydrometeors
  real(wp), allocatable :: mr_gbx(:,:,:)
  real(wp), allocatable :: reff_gbx(:,:,:)
  real(wp), allocatable :: numc_gbx(:,:,:)
  ! points, levels
  real(wp), allocatable :: tca_gbx(:,:)
  real(wp), allocatable :: cca_gbx(:,:)
  real(wp), allocatable :: dems_gbx(:,:)
  real(wp), allocatable :: demc_gbx(:,:)
  real(wp), allocatable :: dtaus_gbx(:,:)
  real(wp), allocatable :: dtauc_gbx(:,:)
end type cosp_inputs_host_model

contains

subroutine construct_cosp_inputs_host_model(np, nc, nl, y)
implicit none
! Inputs: number of horizontal gridpoints, subcolumns, and levels
integer, intent(in) :: np, nc, nl
! Outputs
type(cosp_inputs_host_model ), intent(out) :: y
y%npoints  = np
y%ncolumns = nc
y%nlevels  = nl
! 4D
allocate(y%mr_hydro(np, nc, nl, n_hydro))
allocate(y%reff_hydro(np, nc, nl, n_hydro))
! 3D
allocate(y%mr_gbx(np, nl, n_hydro))
allocate(y%reff_gbx(np, nl, n_hydro))
allocate(y%numc_gbx(np, nl, n_hydro))
! 2D
allocate(y%tca_gbx(np, nl))
allocate(y%cca_gbx(np, nl))
allocate(y%dems_gbx(np, nl))
allocate(y%demc_gbx(np, nl))
allocate(y%dtaus_gbx(np, nl))
allocate(y%dtauc_gbx(np, nl))
y%mr_hydro = 0.0_wp
y%reff_hydro = 0.0_wp
y%mr_gbx = 0.0_wp
y%reff_gbx = 0.0_wp
y%numc_gbx = 0.0_wp
y%tca_gbx = 0.0_wp
y%cca_gbx = 0.0_wp
y%dems_gbx = 0.0_wp
y%demc_gbx = 0.0_wp
y%dtaus_gbx = 0.0_wp
y%dtauc_gbx = 0.0_wp
end subroutine construct_cosp_inputs_host_model

subroutine destroy_cosp_inputs_host_model(y)
implicit none
type(cosp_inputs_host_model ), intent(in out) :: y
if (allocated(y%dtauc_gbx)) deallocate(y%demc_gbx)
if (allocated(y%dtaus_gbx)) deallocate(y%dems_gbx)
if (allocated(y%demc_gbx)) deallocate(y%demc_gbx)
if (allocated(y%dems_gbx)) deallocate(y%dems_gbx)
if (allocated(y%cca_gbx)) deallocate(y%cca_gbx)
if (allocated(y%tca_gbx)) deallocate(y%tca_gbx)
if (allocated(y%reff_gbx)) deallocate(y%reff_gbx)
if (allocated(y%numc_gbx)) deallocate(y%numc_gbx)
if (allocated(y%mr_gbx)) deallocate(y%mr_gbx)
if (allocated(y%reff_hydro)) deallocate(y%reff_hydro)
if (allocated(y%mr_hydro)) deallocate(y%mr_hydro)
end subroutine destroy_cosp_inputs_host_model

subroutine construct_cosp_optical_inputs(cosp_cfg, np, nc, nl, y)
implicit none
! Inputs: number of horizontal gridpoints, subcolumns, and levels
type(cosp2_config), intent(in) :: cosp_cfg
integer, intent(in) :: np, nc, nl
! Outputs
type(cosp_optical_inputs), intent(out) :: y

! Dimensions
y%npoints  = np
y%ncolumns = nc
y%nlevels  = nl
y%npart    = 4
y%nrefl    = parasol_nrefl
allocate(y%frac_out(np, nc, nl))
y%frac_out = 0.0_wp

if (cosp_cfg%lmodis .or. cosp_cfg%lmisr .or. cosp_cfg%lisccp) then
  allocate(y%tau_067(np, nc, nl), y%emiss_11(np, nc, nl))
  y%tau_067 = 0.0_wp
  y%emiss_11 = 0.0_wp
end if
if (cosp_cfg%lcalipso) then
  ! 3D
  allocate(y%betatot_calipso(np, nc, nl),                                      &
           y%betatot_ice_calipso(np, nc, nl),                                  &
           y%betatot_liq_calipso(np, nc, nl),                                  &
           y%tautot_calipso(np, nc, nl),                                       &
           y%tautot_ice_calipso(np, nc, nl),                                   &
           y%tautot_liq_calipso(np, nc, nl))
  ! 2D
  allocate(y%beta_mol_calipso(np, nl), y%tau_mol_calipso(np, nl),              &
            y%tautot_s_ice(np, nc), y%tautot_s_liq(np, nc))
  y%betatot_calipso = 0.0_wp
  y%betatot_ice_calipso = 0.0_wp
  y%betatot_liq_calipso = 0.0_wp
  y%tautot_calipso = 0.0_wp
  y%tautot_ice_calipso = 0.0_wp
  y%tautot_liq_calipso = 0.0_wp
  y%beta_mol_calipso = 0.0_wp
  y%tau_mol_calipso = 0.0_wp
  y%tautot_s_ice = 0.0_wp
  y%tautot_s_liq = 0.0_wp
end if

if (cosp_cfg%lgrlidar532) then
  ! 3D
  allocate(y%betatot_grlidar532(np, nc, nl),                                   &
            y%tautot_grlidar532(np, nc, nl))
  ! 2D
  allocate(y%beta_mol_grlidar532(np, nl), y%tau_mol_grlidar532(np, nl))
  y%betatot_grlidar532 = 0.0_wp
  y%tautot_grlidar532 = 0.0_wp
  y%beta_mol_grlidar532 = 0.0_wp
  y%tau_mol_grlidar532 = 0.0_wp
end if

if (cosp_cfg%latlid) then
  ! 3D
  allocate(y%betatot_atlid(np, nc, nl), y%tautot_atlid(np, nc, nl))
  ! 2D
  allocate(y%beta_mol_atlid(np, nl), y%tau_mol_atlid(np, nl))
  y%betatot_atlid = 0.0_wp
  y%tautot_atlid = 0.0_wp
  y%beta_mol_atlid = 0.0_wp
  y%tau_mol_atlid = 0.0_wp
end if

if (cosp_cfg%lcloudsat) then
  ! 3D
  allocate(y%z_vol_cloudsat(np, nc, nl), y%kr_vol_cloudsat(np, nc, nl),        &
           y%g_vol_cloudsat(np, nc, nl))
  ! 2D
  allocate(y%fracprecipice(np, nc))
  y%z_vol_cloudsat = 0.0_wp
  y%kr_vol_cloudsat = 0.0_wp
  y%g_vol_cloudsat = 0.0_wp
  y%fracprecipice = 0.0_wp
end if

if (cosp_cfg%lmodis) then
  ! 3D
  allocate(y%fracliq(np, nc, nl), y%asym(np, nc, nl), y%ss_alb(np, nc, nl))
  y%fracliq = 0.0_wp
  y%asym = 0.0_wp
  y%ss_alb = 0.0_wp
end if
end subroutine construct_cosp_optical_inputs

subroutine construct_cosp_column_inputs(np,nc,nl,nchan,y)
implicit none
! Inputs: number of horizontal gridpoints, ncolumns, nlevels, and channels
integer, intent(in) :: np, nc, nl, nchan
! Outputs
type(cosp_column_inputs), intent(out) :: y

! Dimensions
y%npoints  = np
y%ncolumns = nc
y%nlevels  = nl
! 1D
allocate(y%sunlit(np), y%skt(np), y%land(np), y%u_sfc(np), y%v_sfc(np),        &
         y%lat(np), y%lon(np), y%emis_sfc(nchan), y%seaice(np),                &
         y%surfelev(np))
! 2D
allocate(y%at(np,nl), y%pfull(np,nl), y%phalf(np,nl+1), y%qv(np,nl),           &
         y%o3(np,nl), y%hgt_matrix(np,nl), y%cloudice(np,nl),                  &
         y%cloudliq(np,nl), y%fl_snow(np,nl), y%fl_rain(np,nl),                &
         y%tca(np,nl), y%hgt_matrix_half(np,nl+1))
end subroutine construct_cosp_column_inputs

subroutine construct_cosp_outputs(cosp_cfg, np, nc, nl, nlvg, nchan, x)
implicit none
! Inputs
! Configuration logicals, number of horizontal gridpoints,
! subcolumns, levels, and levels in new grid
type(cosp2_config), intent(in) :: cosp_cfg
integer, intent(in) :: np, nc, nl, nlvg, nchan
! Outputs
type(cosp_outputs), intent(out) :: x

 ! ISCCP simulator outputs
if (cosp_cfg%lboxtauisccp)    allocate(x%isccp_boxtau(np, nc))
if (cosp_cfg%lboxptopisccp)   allocate(x%isccp_boxptop(np, nc))
if (cosp_cfg%lclisccp)        allocate(x%isccp_fq(np, numisccptaubins,         &
                                       numisccppresbins))
if (cosp_cfg%lcltisccp)       allocate(x%isccp_totalcldarea(np))
if (cosp_cfg%lpctisccp)       allocate(x%isccp_meanptop(np))
if (cosp_cfg%ltauisccp)       allocate(x%isccp_meantaucld(np))
if (cosp_cfg%lmeantbisccp)    allocate(x%isccp_meantb(np))
if (cosp_cfg%lmeantbclrisccp) allocate(x%isccp_meantbclr(np))
if (cosp_cfg%lalbisccp)       allocate(x%isccp_meanalbedocld(np))

! MISR simulator
if (cosp_cfg%lclmisr) then
  allocate(x%misr_fq(np, nummisrtaubins, nummisrhgtbins))
  ! These 3 fields are not output, but were part of the v1.4.0 cosp_misr, so
  ! they are still computed. should probably have a logical to control these
  ! outputs.
  allocate(x%misr_dist_model_layertops(np, nummisrhgtbins))
  allocate(x%misr_meanztop(np))
  allocate(x%misr_cldarea(np))
end if

! MODIS simulator
if (cosp_cfg%lcltmodis)     allocate(x%modis_cloud_fraction_total_mean(np))
if (cosp_cfg%lclwmodis)     allocate(x%modis_cloud_fraction_water_mean(np))
if (cosp_cfg%lclimodis)     allocate(x%modis_cloud_fraction_ice_mean(np))
if (cosp_cfg%lclhmodis)     allocate(x%modis_cloud_fraction_high_mean(np))
if (cosp_cfg%lclmmodis)     allocate(x%modis_cloud_fraction_mid_mean(np))
if (cosp_cfg%lcllmodis)     allocate(x%modis_cloud_fraction_low_mean(np))
if (cosp_cfg%ltautmodis) allocate(x%modis_optical_thickness_total_mean(np))
if (cosp_cfg%ltauwmodis) allocate(x%modis_optical_thickness_water_mean(np))
if (cosp_cfg%ltauimodis) allocate(x%modis_optical_thickness_ice_mean(np))
if (cosp_cfg%ltautlogmodis)                                                    &
  allocate(x%modis_optical_thickness_total_logmean(np))
if (cosp_cfg%ltauwlogmodis)                                                    &
  allocate(x%modis_optical_thickness_water_logmean(np))
if (cosp_cfg%ltauilogmodis)                                                    &
  allocate(x%modis_optical_thickness_ice_logmean(np))
if (cosp_cfg%lreffclwmodis)                                                    &
  allocate(x%modis_cloud_particle_size_water_mean(np))
if (cosp_cfg%lreffclimodis)                                                    &
  allocate(x%modis_cloud_particle_size_ice_mean(np))
if (cosp_cfg%lpctmodis) allocate(x%modis_cloud_top_pressure_total_mean(np))
if (cosp_cfg%llwpmodis) allocate(x%modis_liquid_water_path_mean(np))
if (cosp_cfg%liwpmodis) allocate(x%modis_ice_water_path_mean(np))
if (cosp_cfg%lclmodis) then
  allocate(x%modis_optical_thickness_vs_cloud_top_pressure(np,                 &
      nummodistaubins, nummodispresbins))
  allocate(x%modis_optical_thickness_vs_reffliq(np, nummodistaubins,           &
      nummodisreffliqbins))
  allocate(x%modis_optical_thickness_vs_reffice(np, nummodistaubins,           &
      nummodisrefficebins))
end if

! LIDAR simulator
if (cosp_cfg%llidarbetamol532) allocate(x%calipso_beta_mol(np, nl))
if (cosp_cfg%latb532)          allocate(x%calipso_beta_tot(np, nc, nl))
if (cosp_cfg%lcfadlidarsr532) then
  allocate(x%calipso_srbval(sr_bins+1))
  allocate(x%calipso_cfad_sr(np, sr_bins, nlvg))
  allocate(x%calipso_betaperp_tot(np, nc, nl))
end if
if (cosp_cfg%lclcalipso)       allocate(x%calipso_lidarcld(np, nlvg))
if (cosp_cfg%lclhcalipso .or. cosp_cfg%lclmcalipso .or.                        &
    cosp_cfg%lcllcalipso .or. cosp_cfg%lcltcalipso) then
  allocate(x%calipso_cldlayer(np,lidar_ncat))
end if
if (cosp_cfg%lclcalipsoice .or. cosp_cfg%lclcalipsoliq .or.                    &
    cosp_cfg%lclcalipsoun) then
  allocate(x%calipso_lidarcldphase(np, nlvg, 6))
end if
if (cosp_cfg%lclcalipsotmp .or. cosp_cfg%lclcalipsotmpliq .or.                 &
    cosp_cfg%lclcalipsoice .or. cosp_cfg%lclcalipsotmpun .or.                  &
    cosp_cfg%lclcalipsotmpice) then
  allocate(x%calipso_lidarcldtmp(np, lidar_ntemp, 5))
end if
if (cosp_cfg%lcllcalipsoice .or. cosp_cfg%lclmcalipsoice .or.                  &
    cosp_cfg%lclhcalipsoice .or. cosp_cfg%lcltcalipsoice .or.                  &
    cosp_cfg%lcllcalipsoliq .or. cosp_cfg%lclmcalipsoliq .or.                  &
    cosp_cfg%lclhcalipsoliq .or. cosp_cfg%lcltcalipsoliq .or.                  &
    cosp_cfg%lcllcalipsoun  .or. cosp_cfg%lclmcalipsoun  .or.                  &
    cosp_cfg%lclhcalipsoun  .or. cosp_cfg%lcltcalipsoun) then
  allocate(x%calipso_cldlayerphase(np, lidar_ncat, 6))
end if
if (cosp_cfg%lclopaquecalipso .or. cosp_cfg%lclthincalipso .or.                &
    cosp_cfg%lclzopaquecalipso) then
  allocate(x%calipso_cldtype(np, lidar_ntype))
end if
if (cosp_cfg%lclopaquetemp .or. cosp_cfg%lclthintemp .or.                      &
    cosp_cfg%lclzopaquetemp) then
  allocate(x%calipso_cldtypetemp(np, lidar_ntype))
end if
if (cosp_cfg%lclopaquemeanz .or. cosp_cfg%lclthinmeanz) then
  allocate(x%calipso_cldtypemeanz(np, 2))
end if
if (cosp_cfg%lclopaquemeanzse .or. cosp_cfg%lclthinmeanzse .or.                &
    cosp_cfg%lclzopaquecalipsose) then
  allocate(x%calipso_cldtypemeanzse(np, 3))
end if
if (cosp_cfg%lclthinemis) allocate(x%calipso_cldthinemis(np))
if (cosp_cfg%lclcalipsoopaque .or. cosp_cfg%lclcalipsothin .or.                &
    cosp_cfg%lclcalipsozopaque .or. cosp_cfg%lclcalipsoopacity) then
  allocate(x%calipso_lidarcldtype(np, nlvg, lidar_ntype+1))
end if

! GROUND LIDAR @ 532NM simulator
if (cosp_cfg%llidarbetamol532gr) allocate(x%grlidar532_beta_mol(np, nl))
if (cosp_cfg%latb532gr)          allocate(x%grlidar532_beta_tot(np, nc, nl))
if (cosp_cfg%lcfadlidarsr532gr) then
  allocate(x%grlidar532_srbval(sr_bins+1))
  allocate(x%grlidar532_cfad_sr(np, sr_bins, nlvg))
end if
if (cosp_cfg%lclgrlidar532)     allocate(x%grlidar532_lidarcld(np, nlvg))
if (cosp_cfg%lclhgrlidar532 .or. cosp_cfg%lclmgrlidar532 .or.                  &
    cosp_cfg%lcllgrlidar532 .or. cosp_cfg%lcltgrlidar532) then
  allocate(x%grlidar532_cldlayer(np, lidar_ncat))
end if

! ATLID simulator
if (cosp_cfg%llidarbetamol355) allocate(x%atlid_beta_mol(np, nl))
if (cosp_cfg%latb355)          allocate(x%atlid_beta_tot(np, nc, nl))
if (cosp_cfg%lcfadlidarsr355) then
  allocate(x%atlid_srbval(sr_bins+1))
  allocate(x%atlid_cfad_sr(np, sr_bins, nlvg))
end if
if (cosp_cfg%lclatlid)     allocate(x%atlid_lidarcld(np, nlvg))
if (cosp_cfg%lclhatlid .or. cosp_cfg%lclmatlid .or.                            &
    cosp_cfg%lcllatlid .or. cosp_cfg%lcltatlid) then
  allocate(x%atlid_cldlayer(np, lidar_ncat))
end if

! PARASOL
if (cosp_cfg%lparasolrefl) then
  allocate(x%parasolpix_refl(np, nc, parasol_nrefl))
  allocate(x%parasolgrid_refl(np, parasol_nrefl))
end if

! Cloudsat simulator
if (cosp_cfg%ldbze94)        allocate(x%cloudsat_ze_tot(np, nc, nl))
if (cosp_cfg%lcfaddbze94)    allocate(x%cloudsat_cfad_ze(np,                   &
                                      cloudsat_dbze_bins, nlvg))
if (cosp_cfg%lptradarflag0 .or. cosp_cfg%lptradarflag1 .or.                    &
    cosp_cfg%lptradarflag2 .or. cosp_cfg%lptradarflag3 .or.                    &
    cosp_cfg%lptradarflag4 .or. cosp_cfg%lptradarflag5 .or.                    &
    cosp_cfg%lptradarflag6 .or. cosp_cfg%lptradarflag7 .or.                    &
    cosp_cfg%lptradarflag8 .or. cosp_cfg%lptradarflag9) then
  allocate(x%cloudsat_precip_cover(np, cloudsat_dbze_bins))
end if
if (cosp_cfg%lradarpia) allocate(x%cloudsat_pia(np))

! Combined CALIPSO/CLOUDSAT fields
if (cosp_cfg%lclcalipso2)    allocate(x%lidar_only_freq_cloud(np, nlvg))
if (cosp_cfg%lcltlidarradar) allocate(x%radar_lidar_tcc(np))
if (cosp_cfg%lcloudsat_tcc)  allocate(x%cloudsat_tcc(np))
if (cosp_cfg%lcloudsat_tcc2) allocate(x%cloudsat_tcc2(np))

! RTTOV
if (cosp_cfg%ltbrttov) allocate(x%rttov_tbs(np, nchan))

! Joint MODIS/CloudSat Statistics
if (cosp_cfg%lwr_occfreq)  allocate(x%wr_occfreq_ntotal(np, wr_nregime))
if (cosp_cfg%lcfodd)       allocate(x%cfodd_ntotal(np, cfodd_ndbze,            &
                                    cfodd_nicod, cfodd_nclass))
end subroutine construct_cosp_outputs

subroutine destroy_cosp_optical_inputs(y)
implicit none
type(cosp_optical_inputs), intent(in out) :: y

if (allocated(y%tau_067))             deallocate(y%tau_067)
if (allocated(y%emiss_11))            deallocate(y%emiss_11)
if (allocated(y%frac_out))            deallocate(y%frac_out)
if (allocated(y%beta_mol_calipso))    deallocate(y%beta_mol_calipso)
if (allocated(y%tau_mol_calipso))     deallocate(y%tau_mol_calipso)
if (allocated(y%betatot_calipso))     deallocate(y%betatot_calipso)
if (allocated(y%betatot_ice_calipso)) deallocate(y%betatot_ice_calipso)
if (allocated(y%betatot_liq_calipso)) deallocate(y%betatot_liq_calipso)
if (allocated(y%tautot_calipso))      deallocate(y%tautot_calipso)
if (allocated(y%tautot_ice_calipso))  deallocate(y%tautot_ice_calipso)
if (allocated(y%tautot_liq_calipso))  deallocate(y%tautot_liq_calipso)
if (allocated(y%tautot_s_liq))        deallocate(y%tautot_s_liq)
if (allocated(y%tautot_s_ice))        deallocate(y%tautot_s_ice)
if (allocated(y%z_vol_cloudsat))      deallocate(y%z_vol_cloudsat)
if (allocated(y%kr_vol_cloudsat))     deallocate(y%kr_vol_cloudsat)
if (allocated(y%g_vol_cloudsat))      deallocate(y%g_vol_cloudsat)
if (allocated(y%asym))                deallocate(y%asym)
if (allocated(y%ss_alb))              deallocate(y%ss_alb)
if (allocated(y%fracliq))             deallocate(y%fracliq)
if (allocated(y%beta_mol_grlidar532)) deallocate(y%beta_mol_grlidar532)
if (allocated(y%betatot_grlidar532))  deallocate(y%betatot_grlidar532)
if (allocated(y%tau_mol_grlidar532))  deallocate(y%tau_mol_grlidar532)
if (allocated(y%tautot_grlidar532))   deallocate(y%tautot_grlidar532)
if (allocated(y%beta_mol_atlid))      deallocate(y%beta_mol_atlid)
if (allocated(y%betatot_atlid))       deallocate(y%betatot_atlid)
if (allocated(y%tau_mol_atlid))       deallocate(y%tau_mol_atlid)
if (allocated(y%tautot_atlid))        deallocate(y%tautot_atlid)
if (allocated(y%fracprecipice))       deallocate(y%fracprecipice)
if (allocated(y%rcfg_cloudsat%n_scale_flag)) &
   deallocate(y%rcfg_cloudsat%n_scale_flag)
if (allocated(y%rcfg_cloudsat%z_scale_flag)) &
   deallocate(y%rcfg_cloudsat%z_scale_flag)
if (allocated(y%rcfg_cloudsat%z_scale_added_flag)) &
   deallocate(y%rcfg_cloudsat%z_scale_added_flag)
if (allocated(y%rcfg_cloudsat%ze_scaled)) &
   deallocate(y%rcfg_cloudsat%ze_scaled)
if (allocated(y%rcfg_cloudsat%zr_scaled)) &
   deallocate(y%rcfg_cloudsat%zr_scaled)
if (allocated(y%rcfg_cloudsat%kr_scaled)) &
   deallocate(y%rcfg_cloudsat%kr_scaled)
if (allocated(y%rcfg_cloudsat%fc)) &
   deallocate(y%rcfg_cloudsat%fc)
if (allocated(y%rcfg_cloudsat%rho_eff)) &
   deallocate(y%rcfg_cloudsat%rho_eff)
if (allocated(y%rcfg_cloudsat%base_list)) &
   deallocate(y%rcfg_cloudsat%base_list)
if (allocated(y%rcfg_cloudsat%step_list)) &
   deallocate(y%rcfg_cloudsat%step_list)
end subroutine destroy_cosp_optical_inputs

subroutine destroy_cosp_column_inputs(y)
implicit none
type(cosp_column_inputs), intent(in out) :: y

if (allocated(y%sunlit))          deallocate(y%sunlit)
if (allocated(y%skt))             deallocate(y%skt)
if (allocated(y%land))            deallocate(y%land)
if (allocated(y%at))              deallocate(y%at)
if (allocated(y%pfull))           deallocate(y%pfull)
if (allocated(y%phalf))           deallocate(y%phalf)
if (allocated(y%qv))              deallocate(y%qv)
if (allocated(y%o3))              deallocate(y%o3)
if (allocated(y%hgt_matrix))      deallocate(y%hgt_matrix)
if (allocated(y%u_sfc))           deallocate(y%u_sfc)
if (allocated(y%v_sfc))           deallocate(y%v_sfc)
if (allocated(y%lat))             deallocate(y%lat)
if (allocated(y%lon))             deallocate(y%lon)
if (allocated(y%emis_sfc))        deallocate(y%emis_sfc)
if (allocated(y%cloudice))        deallocate(y%cloudice)
if (allocated(y%cloudliq))        deallocate(y%cloudliq)
if (allocated(y%seaice))          deallocate(y%seaice)
if (allocated(y%fl_rain))         deallocate(y%fl_rain)
if (allocated(y%fl_snow))         deallocate(y%fl_snow)
if (allocated(y%tca))             deallocate(y%tca)
if (allocated(y%hgt_matrix_half)) deallocate(y%hgt_matrix_half)
if (allocated(y%surfelev))        deallocate(y%surfelev)
end subroutine destroy_cosp_column_inputs

subroutine destroy_cosp_outputs(y)
implicit none
type(cosp_outputs), intent(in out) :: y

! deallocate and nullify
call deallocate_and_nullify(y%calipso_beta_mol)
call deallocate_and_nullify(y%calipso_temp_tot)
call deallocate_and_nullify(y%calipso_betaperp_tot)
call deallocate_and_nullify(y%calipso_beta_tot)
call deallocate_and_nullify(y%calipso_tau_tot)
call deallocate_and_nullify(y%calipso_lidarcldphase)
call deallocate_and_nullify(y%calipso_lidarcldtype)
call deallocate_and_nullify(y%calipso_cldlayerphase)
call deallocate_and_nullify(y%calipso_lidarcldtmp)
call deallocate_and_nullify(y%calipso_cldlayer)
call deallocate_and_nullify(y%calipso_cldtype)
call deallocate_and_nullify(y%calipso_cldtypetemp)
call deallocate_and_nullify(y%calipso_cldtypemeanz)
call deallocate_and_nullify(y%calipso_cldtypemeanzse)
call deallocate_and_nullify(y%calipso_cldthinemis)
call deallocate_and_nullify(y%calipso_lidarcld)
call deallocate_and_nullify(y%calipso_srbval)
call deallocate_and_nullify(y%calipso_cfad_sr)
call deallocate_and_nullify(y%grlidar532_beta_mol)
call deallocate_and_nullify(y%grlidar532_beta_tot)
call deallocate_and_nullify(y%grlidar532_cldlayer)
call deallocate_and_nullify(y%grlidar532_lidarcld)
call deallocate_and_nullify(y%grlidar532_cfad_sr)
call deallocate_and_nullify(y%grlidar532_srbval)
call deallocate_and_nullify(y%atlid_beta_mol)
call deallocate_and_nullify(y%atlid_beta_tot)
call deallocate_and_nullify(y%atlid_cldlayer)
call deallocate_and_nullify(y%atlid_lidarcld)
call deallocate_and_nullify(y%atlid_cfad_sr)
call deallocate_and_nullify(y%atlid_srbval)
call deallocate_and_nullify(y%parasolpix_refl)
call deallocate_and_nullify(y%parasolgrid_refl)
call deallocate_and_nullify(y%cloudsat_ze_tot)
call deallocate_and_nullify(y%cloudsat_cfad_ze)
call deallocate_and_nullify(y%cloudsat_precip_cover)
call deallocate_and_nullify(y%cloudsat_pia)
call deallocate_and_nullify(y%cloudsat_tcc)
call deallocate_and_nullify(y%cloudsat_tcc2)
call deallocate_and_nullify(y%radar_lidar_tcc)
call deallocate_and_nullify(y%cloudsat_tcc)
call deallocate_and_nullify(y%cloudsat_tcc2)
call deallocate_and_nullify(y%lidar_only_freq_cloud)
call deallocate_and_nullify(y%isccp_totalcldarea)
call deallocate_and_nullify(y%isccp_meantb)
call deallocate_and_nullify(y%isccp_meantbclr)
call deallocate_and_nullify(y%isccp_meanptop)
call deallocate_and_nullify(y%isccp_meantaucld)
call deallocate_and_nullify(y%isccp_meanalbedocld)
call deallocate_and_nullify(y%isccp_boxtau)
call deallocate_and_nullify(y%isccp_boxptop)
call deallocate_and_nullify(y%isccp_fq)
call deallocate_and_nullify(y%misr_fq)
call deallocate_and_nullify(y%misr_dist_model_layertops)
call deallocate_and_nullify(y%misr_meanztop)
call deallocate_and_nullify(y%misr_cldarea)
call deallocate_and_nullify(y%rttov_tbs)
call deallocate_and_nullify(y%modis_cloud_fraction_total_mean)
call deallocate_and_nullify(y%modis_cloud_fraction_ice_mean)
call deallocate_and_nullify(y%modis_cloud_fraction_water_mean)
call deallocate_and_nullify(y%modis_cloud_fraction_high_mean)
call deallocate_and_nullify(y%modis_cloud_fraction_mid_mean)
call deallocate_and_nullify(y%modis_cloud_fraction_low_mean)
call deallocate_and_nullify(y%modis_optical_thickness_total_mean)
call deallocate_and_nullify(y%modis_optical_thickness_water_mean)
call deallocate_and_nullify(y%modis_optical_thickness_ice_mean)
call deallocate_and_nullify(y%modis_optical_thickness_total_logmean)
call deallocate_and_nullify(y%modis_optical_thickness_water_logmean)
call deallocate_and_nullify(y%modis_optical_thickness_ice_logmean)
call deallocate_and_nullify(y%modis_cloud_particle_size_water_mean)
call deallocate_and_nullify(y%modis_cloud_particle_size_ice_mean)
call deallocate_and_nullify(y%modis_cloud_top_pressure_total_mean)
call deallocate_and_nullify(y%modis_liquid_water_path_mean)
call deallocate_and_nullify(y%modis_ice_water_path_mean)
call deallocate_and_nullify( &
  y%modis_optical_thickness_vs_cloud_top_pressure)
call deallocate_and_nullify(y%modis_optical_thickness_vs_reffliq)
call deallocate_and_nullify(y%modis_optical_thickness_vs_reffice)
call deallocate_and_nullify(y%cfodd_ntotal)
call deallocate_and_nullify(y%wr_occfreq_ntotal)
end subroutine destroy_cosp_outputs

subroutine deallocate_and_nullify_1d(x)
implicit none
real(wp), pointer :: x(:)
if (associated(x)) then
  deallocate(x)
  nullify(x)
end if
end subroutine deallocate_and_nullify_1d

subroutine deallocate_and_nullify_2d(x)
implicit none
real(wp), pointer :: x(:,:)
if (associated(x)) then
  deallocate(x)
  nullify(x)
end if
end subroutine deallocate_and_nullify_2d

subroutine deallocate_and_nullify_3d(x)
implicit none
real(wp), pointer :: x(:,:,:)
if (associated(x)) then
  deallocate(x)
  nullify(x)
end if
end subroutine deallocate_and_nullify_3d

subroutine deallocate_and_nullify_4d(x)
implicit none
real(wp), pointer :: x(:,:,:,:)
if (associated(x)) then
  deallocate(x)
  nullify(x)
end if
end subroutine deallocate_and_nullify_4d
end module cosp2_types_mod
