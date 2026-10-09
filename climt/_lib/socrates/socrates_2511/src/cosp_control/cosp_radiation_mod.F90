! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module calculates optical properties that are required by some of
! the COSP2 simulators
!
!------------------------------------------------------------------------------

module cosp_radiation_mod
use mod_quickbeam_optics, only: size_distribution, quickbeam_optics, gases
use cosp_optics, only: lidar_optics
use mod_cosp_stats, only: cosp_change_vertical_grid
use cosp_kinds, only: wp
use cosp_input_mod, only: cosp_lidar_ice_type
use cosp2_types_mod, only: cosp2_config, cosp_inputs_host_model
use mod_cosp_config, only: r_undef, n_hydro, nlvgrid, vgrid_zl, vgrid_zu, &
  cloudsat_preclvl, use_vgrid
use mod_cosp, only: cosp_optical_inputs, cosp_column_inputs
use cosp2_constants_mod, only: i_lscliq, i_lscice, i_lsrain, i_lsiagg, &
  i_cvcliq, i_cvcice, i_cvrain, i_cvsnow, i_lsgrpl

use ereport_mod

implicit none

character(len=*), parameter, private :: ModuleName='COSP_RADIATION_MOD'

integer, parameter :: calipso_npart = 4
integer, parameter :: calipso_frequency = 532

contains
subroutine cosp_radiative_properties(cosp_cfg,cosp_hmodel,cosp_column_in, &
  cosp_optical_in,quickbeam_size_dist)

use errormessagelength_mod, only: errormessagelength

implicit none
!----Input arguments
! configuration options
type(cosp2_config),intent(in) :: cosp_cfg
! cosp host model inputs
type(cosp_inputs_host_model), intent(in) :: cosp_hmodel
! cosp column inputs
type(cosp_column_inputs), intent(in) :: cosp_column_in
! cosp optical inputs
type(cosp_optical_inputs), intent(in out) :: cosp_optical_in
! quickbeam size distribution
type(size_distribution), intent(in out) :: quickbeam_size_dist
!----Local variables
integer :: i, j
integer :: npoints, ncolumns, nlevels
integer, allocatable :: cloudsat_preclvl_index(:)
real(wp), allocatable :: g_vol(:,:)
real(wp), allocatable :: fracprecipice(:,:,:)
real(wp), allocatable :: fracprecipice_statgrid(:,:,:)
real(wp), allocatable :: Np(:,:,:)
character(len=25), parameter :: RoutineName='COSP_RADIATIVE_PROPERTIES'

! Use local variables to make the code more readable
npoints  = cosp_column_in%npoints
ncolumns = cosp_column_in%ncolumns
nlevels  = cosp_column_in%nlevels

! Optical properties for CALIPSO simulator
if (cosp_cfg%lcalipso) then
  call lidar_optics(npoints, ncolumns, nlevels, calipso_npart, &
    cosp_lidar_ice_type, calipso_frequency, .false., &
    cosp_hmodel%mr_hydro(:,:,:,i_lscliq), &
    cosp_hmodel%mr_hydro(:,:,:,i_lscice), &
    cosp_hmodel%mr_hydro(:,:,:,i_cvcliq), &
    cosp_hmodel%mr_hydro(:,:,:,i_cvcice), &
    cosp_hmodel%reff_gbx(:,:,i_lscliq), &
    cosp_hmodel%reff_gbx(:,:,i_lscice), &
    cosp_hmodel%reff_gbx(:,:,i_cvcliq), &
    cosp_hmodel%reff_gbx(:,:,i_cvcice), &
    cosp_column_in%pfull, cosp_column_in%phalf, cosp_column_in%at, &
    cosp_optical_in%beta_mol_calipso, cosp_optical_in%betatot_calipso, &
    cosp_optical_in%tau_mol_calipso, cosp_optical_in%tautot_calipso, &
    cosp_optical_in%tautot_s_liq, cosp_optical_in%tautot_s_ice, &
    cosp_optical_in%betatot_ice_calipso, &
    cosp_optical_in%betatot_liq_calipso, &
    cosp_optical_in%tautot_ice_calipso, &
    cosp_optical_in%tautot_liq_calipso)
end if

! Radiative properties for CloudSat simulator
if (cosp_cfg%lcloudsat) then
  ! Compute gaseous absorption (assume identical for each subcolun)
  allocate(g_vol(npoints, nlevels))
  g_vol(:,:)=0.0_wp
  do j = 1,nlevels
    do i = 1,npoints
      if (cosp_optical_in%rcfg_cloudsat%use_gas_abs==1 .or. &
        (cosp_optical_in%rcfg_cloudsat%use_gas_abs==2 .and. j==1)) then
        g_vol(i,j) = gases(cosp_column_in%pfull(i,j), &
        cosp_column_in%at(i,j), cosp_column_in%qv(i,j), &
        cosp_optical_in%rcfg_cloudsat%freq)
      end if
      cosp_optical_in%g_vol_cloudsat(i,:,j)=g_vol(i,j)
    end do
  end do
  deallocate(g_vol)

  ! Loop over all subcolumns
  allocate(fracprecipice(npoints, ncolumns, nlevels))
  ! np not used in the um. if needed, then it will have to be dimensioned
  ! as (npoints, ncolumns, nlevels, n_hydro)
  allocate(np(npoints, nlevels, n_hydro))
  fracprecipice(:,:,:) = 0.0_wp
  np(:,:,:) = 0.0_wp
  do i = 1,ncolumns
    call quickbeam_optics(quickbeam_size_dist,                                 &
      cosp_optical_in%rcfg_cloudsat, npoints, nlevels, r_undef,                &
      cosp_hmodel%mr_hydro(:,i,:,:)*1000.0_wp,                                 &
      cosp_hmodel%reff_hydro(:,i,:,:)*1.0e6_wp,                                &
      np(:,:,:), cosp_column_in%pfull, cosp_column_in%at,                      &
      cosp_column_in%qv, cosp_optical_in%z_vol_cloudsat(:,i,:),                &
      cosp_optical_in%kr_vol_cloudsat(:,i,:))

    ! At each model level, what fraction of the precipitation is frozen?
    where (cosp_hmodel%mr_hydro(:,i,:,i_lsrain) > 0.0_wp .or.                  &
          cosp_hmodel%mr_hydro(:,i,:,i_lsiagg) > 0.0_wp .or.                   &
          cosp_hmodel%mr_hydro(:,i,:,i_cvrain) > 0.0_wp .or.                   &
          cosp_hmodel%mr_hydro(:,i,:,i_cvsnow) > 0.0_wp .or.                   &
          cosp_hmodel%mr_hydro(:,i,:,i_lsgrpl) > 0)
      fracprecipice(:,i,:) = (cosp_hmodel%mr_hydro(:,i,:,i_lsiagg) +           &
        cosp_hmodel%mr_hydro(:,i,:,i_cvsnow) +                                 &
        cosp_hmodel%mr_hydro(:,i,:,i_lsgrpl)) /                                &
        (cosp_hmodel%mr_hydro(:,i,:,i_lsiagg) +                                &
        cosp_hmodel%mr_hydro(:,i,:,i_cvsnow) +                                 &
        cosp_hmodel%mr_hydro(:,i,:,i_lsgrpl) +                                 &
        cosp_hmodel%mr_hydro(:,i,:,i_lsrain) +                                 &
        cosp_hmodel%mr_hydro(:,i,:,i_cvrain))
    end where
  end do
  deallocate(np)

  ! Regrid frozen fraction to Cloudsat/Calipso statistical grid
  if (use_vgrid) then
    allocate(fracprecipice_statgrid(npoints, ncolumns, nlvgrid))
    fracprecipice_statgrid(:,:,:) = 0.0_wp
    call cosp_change_vertical_grid(npoints, ncolumns, nlevels,                 &
      cosp_column_in%hgt_matrix(:,nlevels:1:-1),                               &
      cosp_column_in%hgt_matrix_half(:,nlevels:1:-1),                          &
      fracprecipice(:,:,nlevels:1:-1), nlvgrid, vgrid_zl(nlvgrid:1:-1),        &
      vgrid_zu(nlvgrid:1:-1), fracprecipice_statgrid(:,:,nlvgrid:1:-1))

    ! Find proper layer above de surface elevation to compute precip flags
    ! in cloudsat/calipso statistical grid
    allocate(cloudsat_preclvl_index(npoints))
    cloudsat_preclvl_index(:) = 0.0_wp
    ! Computing altitude index for precip flags calculation
    ! (one layer above surfelev layer)
    cloudsat_preclvl_index(:) = cloudsat_preclvl -                             &
      floor( cosp_column_in%surfelev(:)/(vgrid_zl(1)-vgrid_zl(2)) )
    ! For near-surface diagnostics, only need the
    ! frozen fraction at one layer
    do i=1,npoints
      cosp_optical_in%fracprecipice(i,:) =                                     &
      fracprecipice_statgrid(i,:,cloudsat_preclvl_index(i))
    end do

    ! Free memory of local arrays
    deallocate(cloudsat_preclvl_index)
    deallocate(fracprecipice_statgrid)
  end if
  deallocate(fracprecipice)
end if
end subroutine cosp_radiative_properties
end module cosp_radiation_mod
