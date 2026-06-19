! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This main program is used to drive the COSP2 code calling the same
! interface that is called within LFRic for offline tests.
!
!------------------------------------------------------------------------------

program cosp_driver

use realtype_rd, only: RealExt
use cosp_kinds, only: wp
use mod_cosp_io, only: &
    nc_read_input_file, &
    write_cosp2_output
use mod_cosp, only: cosp_outputs
use cosp_mod, only: &
    cosp, &
    n_backscatter_bins, &
    n_isccp_tau_bins, &
    n_isccp_pressure_bins, &
    n_cloudsat_levels
use cosp_input_mod, only: n_hydro
use cosp_def_diag, only: cospdiag

implicit none

  ! Test data
  character(len=600) :: filein
  character(len=256) :: foutput

! variables as used in github test code
  integer :: &
       nlon,nlat,geomode
  real(wp) :: &
       emsfc_lw
  real(wp),dimension(:),allocatable,target:: &
       lon,       & ! longitude (deg)
       lat,       & ! latitude (deg)
       skt,       & ! skin temperature (k)
       surfelev,  & ! surface elevation (m)
       landmask,  & ! land/sea mask (0/1)
       u_wind,    & ! u-component of wind (m/s)
       v_wind,    & ! v-component of wind (m/s)
       sunlit       ! sunlit flag
  real(wp),dimension(:,:),allocatable,target :: &
       p,         & ! model pressure levels (pa)
       ph,        & ! moddel pressure @ half levels (pa)
       zlev,      & ! model level height (m)
       zlev_half, & ! model level height @ half-levels (m)
       t,         & ! temperature (k)
       sh,        & ! specific humidity (kg/kg)
       rh,        & ! relative humidity (1)
       tca,       & ! total cloud fraction (1)
       cca,       & ! convective cloud fraction (1) 
       mr_lsliq,  & ! mass mixing ratio for stratiform cloud liquid (kg/kg)
       mr_lsice,  & ! mass mixing ratio for stratiform cloud ice (kg/kg)
       mr_ccliq,  & ! mass mixing ratio for convective cloud liquid (kg/kg)
       mr_ccice,  & ! mass mixing ratio for convective cloud ice (kg/kg)
       mr_ozone,  & ! mass mixing ratio for ozone (kg/kg)
       fl_lsrain, & ! precipitation flux (rain) for stratiform cloud (kg/m^2/s)
       fl_lssnow, & ! precipitation flux (snow) for stratiform cloud (kg/m^2/s)
       fl_lsgrpl, & ! precipitation flux (groupel) for stratiform cloud (kg/m^2/s)
       fl_ccrain, & ! precipitation flux (rain) for convective cloud (kg/m^2/s)
       fl_ccsnow, & ! precipitation flux (snow) for convective cloud (kg/m^2/s)
       dtau_s,    & ! 0.67micron optical depth (stratiform cloud) (1)
       dtau_c,    & ! 0.67micron optical depth (convective cloud) (1)
       dem_s,     & ! 11micron emissivity (stratiform cloud) 
       dem_c        ! 11microm emissivity (convective cloud)
  real(wp),dimension(:,:,:),allocatable,target :: &
       reff         ! subcolumn effective radius

integer :: nlevels
integer :: n_profile_full
integer :: n_profile_list
integer :: ncolumns
integer :: nclds
integer :: lun

integer, pointer :: ncldy(:)
integer, allocatable :: profile_list(:)

real(RealExt), pointer :: p_full_levels(:,:)
real(RealExt), pointer :: p_half_levels(:,:)
real(RealExt), pointer :: hgt_full_levels(:,:)
real(RealExt), pointer :: hgt_half_levels(:,:)
real(RealExt), pointer :: d_mass(:,:)
real(RealExt), pointer :: T_n(:,:)
real(RealExt), pointer :: q_n(:,:)
real(RealExt), pointer :: w_cloud(:,:)
real(RealExt), pointer :: condensed_mix_ratio_water(:,:)
real(RealExt), pointer :: condensed_mix_ratio_ice(:,:)
real(RealExt), pointer :: cosp_crain(:,:)
real(RealExt), pointer :: cosp_csnow(:,:)
real(RealExt), pointer :: condensed_re_water(:,:)
real(RealExt), pointer :: condensed_re_ice(:,:)
real(RealExt), pointer :: cloud_extinction(:,:)
real(RealExt), pointer :: cloud_absorptivity(:,:)
real(RealExt), pointer :: frac_cloud_water(:,:)
real(RealExt), pointer :: frac_cloud_ice(:,:)
real(RealExt), pointer :: clw_sub_full(:,:,:)
real(RealExt), pointer :: t_surf(:), p_surf(:), hgt_surf(:)
real(RealExt), pointer :: cosp_sunlit(:)

real(RealExt), pointer :: d_p_full_levels(:,:)
real(RealExt), pointer :: d_p_half_levels(:,:)
real(RealExt), pointer :: d_hgt_full_levels(:,:)
real(RealExt), pointer :: d_hgt_half_levels(:,:)
real(RealExt), pointer :: d_d_mass(:,:)
real(RealExt), pointer :: d_T_n(:,:)
real(RealExt), pointer :: d_q_n(:,:)
real(RealExt), pointer :: d_w_cloud(:,:)
real(RealExt), pointer :: d_condensed_mix_ratio_water(:,:)
real(RealExt), pointer :: d_condensed_mix_ratio_ice(:,:)
real(RealExt), pointer :: d_cosp_crain(:,:)
real(RealExt), pointer :: d_cosp_csnow(:,:)
real(RealExt), pointer :: d_condensed_re_water(:,:)
real(RealExt), pointer :: d_condensed_re_ice(:,:)
real(RealExt), pointer :: d_cloud_extinction(:,:)
real(RealExt), pointer :: d_cloud_absorptivity(:,:)
real(RealExt), pointer :: d_frac_cloud_water(:,:)
real(RealExt), pointer :: d_frac_cloud_ice(:,:)
real(RealExt), pointer :: d_clw_sub_full(:,:,:)

! Microphysics constants
real(RealExt) :: x1r
real(RealExt) :: x1g
real(RealExt) :: x2r
real(RealExt) :: x2g
real(RealExt) :: x4g

type(cosp_outputs) :: cosp_out_ext

type(cospdiag) :: cosp_diag

integer :: i, j, k, l, rep
integer :: ierr

real(RealExt), pointer :: cosp_calipso_low_level_cl_mask(:) ! 2321
real(RealExt), pointer :: cosp_calipso_mid_level_cl_mask(:) ! 2322
real(RealExt), pointer :: cosp_calipso_high_level_cl_mask(:) ! 2323

real(RealExt), pointer :: cosp_calipso_low_level_cl(:) ! 2344
real(RealExt), pointer :: cosp_calipso_mid_level_cl(:) ! 2345
real(RealExt), pointer :: cosp_calipso_high_level_cl(:) ! 2346

real(RealExt), pointer :: cosp_cloud_weights(:) ! 2330
real(RealExt), pointer :: cosp_ctp_tau_histogram(:) ! 2337

real(RealExt), pointer :: cosp_calipso_tot_backscatter(:) ! 2341

real(RealExt), pointer :: cosp_calipso_cf_40_mask(:) ! 2325

real(RealExt), pointer :: cosp_calipso_cfad_sr_40(:) ! 2370

real(RealExt), pointer :: cosp_calipso_cf_40_liq(:) ! 2473
real(RealExt), pointer :: cosp_calipso_cf_40_ice(:) ! 2474
real(RealExt), pointer :: cosp_calipso_cf_40_undet(:) ! 2475

real(wp), allocatable :: cosp_lon(:)
real(wp), allocatable :: cosp_lat(:)

logical :: l_profile_last


l_profile_last = .true.


! Reading test data
call get_command_argument(1, filein)
filein=trim(filein)

n_profile_full = 153
nlevels = 38

allocate(lon(n_profile_full),lat(n_profile_full),p(n_profile_full,nlevels),ph(n_profile_full,nlevels),             &
         zlev(n_profile_full,nlevels),zlev_half(n_profile_full,nlevels),t(n_profile_full,nlevels),          &
         sh(n_profile_full,nlevels),rh(n_profile_full,nlevels),tca(n_profile_full,nlevels),                 &
         cca(n_profile_full,nlevels),mr_lsliq(n_profile_full,nlevels),mr_lsice(n_profile_full,nlevels),     &
         mr_ccliq(n_profile_full,nlevels),mr_ccice(n_profile_full,nlevels),                          &
         fl_lsrain(n_profile_full,nlevels),fl_lssnow(n_profile_full,nlevels),                        &
         fl_lsgrpl(n_profile_full,nlevels),fl_ccrain(n_profile_full,nlevels),                        &
         fl_ccsnow(n_profile_full,nlevels),reff(n_profile_full,nlevels,n_hydro),                     &
         dtau_s(n_profile_full,nlevels),dtau_c(n_profile_full,nlevels),dem_s(n_profile_full,nlevels),       &
         dem_c(n_profile_full,nlevels),skt(n_profile_full),landmask(n_profile_full),                        &
         mr_ozone(n_profile_full,nlevels),u_wind(n_profile_full),v_wind(n_profile_full),sunlit(n_profile_full),    &
         surfelev(n_profile_full))

call nc_read_input_file(filein,n_profile_full,nlevels,n_hydro,lon,lat,p,ph,zlev,zlev_half,    &
                        t,sh,rh,tca,cca,mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain, &
                        fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,reff,dtau_s,dtau_c,    &
                        dem_s,dem_c,skt,landmask,mr_ozone,u_wind,v_wind,sunlit,        &
                        emsfc_lw,geomode,nlon,nlat,surfelev)


allocate(ncldy(n_profile_full))

n_profile_list = 153
allocate(profile_list(n_profile_list))

do l=1, n_profile_list
  profile_list(l) = l
end do ! n_profile_list


ncolumns = 1

if (l_profile_last) then
  allocate(p_full_levels(nlevels,n_profile_full))
  allocate(p_half_levels(nlevels,n_profile_full))
  allocate(hgt_full_levels(nlevels,n_profile_full))
  allocate(hgt_half_levels(nlevels,n_profile_full))
  allocate(d_mass(nlevels,n_profile_full))
  allocate(t_n(nlevels,n_profile_full))
  allocate(q_n(nlevels,n_profile_full))
  allocate(w_cloud(nlevels,n_profile_full))
  allocate(condensed_mix_ratio_water(nlevels,n_profile_full))
  allocate(condensed_mix_ratio_ice(nlevels,n_profile_full))
  allocate(cosp_crain(nlevels,n_profile_full))
  allocate(cosp_csnow(nlevels,n_profile_full))
  allocate(condensed_re_water(nlevels,n_profile_full))
  allocate(condensed_re_ice(nlevels,n_profile_full))
  allocate(cloud_extinction(nlevels,n_profile_full))
  allocate(cloud_absorptivity(nlevels,n_profile_full))
  allocate(frac_cloud_water(nlevels,n_profile_full))
  allocate(frac_cloud_ice(nlevels,n_profile_full))
  allocate(clw_sub_full(nlevels,ncolumns,n_profile_full))
else
  allocate(p_full_levels(n_profile_full,nlevels))
  allocate(p_half_levels(n_profile_full,nlevels))
  allocate(hgt_full_levels(n_profile_full,nlevels))
  allocate(hgt_half_levels(n_profile_full,nlevels))
  allocate(d_mass(n_profile_full,nlevels))
  allocate(t_n(n_profile_full,nlevels))
  allocate(q_n(n_profile_full,nlevels))
  allocate(w_cloud(n_profile_full,nlevels))
  allocate(condensed_mix_ratio_water(n_profile_full,nlevels))
  allocate(condensed_mix_ratio_ice(n_profile_full,nlevels))
  allocate(cosp_crain(n_profile_full,nlevels))
  allocate(cosp_csnow(n_profile_full,nlevels))
  allocate(condensed_re_water(n_profile_full,nlevels))
  allocate(condensed_re_ice(n_profile_full,nlevels))
  allocate(cloud_extinction(n_profile_full,nlevels))
  allocate(cloud_absorptivity(n_profile_full,nlevels))
  allocate(frac_cloud_water(n_profile_full,nlevels))
  allocate(frac_cloud_ice(n_profile_full,nlevels))
  allocate(clw_sub_full(n_profile_full,nlevels,ncolumns))
end if ! l_profile_last

nclds = nlevels
ncldy = ncolumns

p_full_levels = 0.0
p_half_levels = 0.0
hgt_full_levels = 0.0
hgt_half_levels = 0.0
t_n = 0.0
q_n = 0.0
w_cloud = 0.0
condensed_mix_ratio_water = 0.0
condensed_mix_ratio_ice = 0.0
cosp_crain = 0.0
cosp_csnow = 0.0
condensed_re_water = 0.0
condensed_re_ice = 0.0
cloud_extinction = 0.0
cloud_absorptivity = 0.0

if (l_profile_last) then
  do i=1, nlevels
    do j=1, n_profile_full
      p_full_levels(i,j) = p(j,i)
      p_half_levels(i,j) = ph(j,i)
      hgt_full_levels(i,j) = zlev(j,i)
      hgt_half_levels(i,j) = zlev_half(j,i)
      t_n(i,j) = t(j,i)
      q_n(i,j) = sh(j,i)
      w_cloud(i,j) = tca(j,i)
      condensed_mix_ratio_water(i,j) = mr_lsliq(j,i) + mr_ccliq(j,i)
      condensed_mix_ratio_ice(i,j) = mr_lsice(j,i) + mr_ccice(j,i)
      cosp_crain(i,j) = fl_lsrain(j,i) + fl_ccrain(j,i)
      cosp_csnow(i,j) = fl_lssnow(j,i) + fl_ccsnow(j,i)
      condensed_re_water(i,j) = reff(j,i,1)
      condensed_re_ice(i,j) = reff(j,i,1)
      cloud_extinction(i,j) = dtau_s(j,i) + dtau_c(j,i)
      cloud_absorptivity(i,j) = dem_s(j,i) + dem_c(j,i)
    end do ! j
  end do ! i
else
  do i=1, nlevels
    do j=1, n_profile_full
      p_full_levels(j,i) = p(j,i)
      p_half_levels(j,i) = ph(j,i)
      hgt_full_levels(j,i) = zlev(j,i)
      hgt_half_levels(j,i) = zlev_half(j,i)
      t_n(j,i) = t(j,i)
      q_n(j,i) = sh(j,i)
      w_cloud(j,i) = tca(j,i)
      condensed_mix_ratio_water(j,i) = mr_lsliq(j,i) + mr_ccliq(j,i)
      condensed_mix_ratio_ice(j,i) = mr_lsice(j,i) + mr_ccice(j,i)
      cosp_crain(j,i) = fl_lsrain(j,i) + fl_ccrain(j,i)
      cosp_csnow(j,i) = fl_lssnow(j,i) + fl_ccsnow(j,i)
      condensed_re_water(j,i) = reff(j,i,1)
      condensed_re_ice(j,i) = reff(j,i,1)
      cloud_extinction(j,i) = dtau_s(j,i) + dtau_c(j,i)
      cloud_absorptivity(j,i) = dem_s(j,i) + dem_c(j,i)
    end do ! j
  end do ! i
endif ! l_profile_last

deallocate(p)
deallocate(zlev_half)
deallocate(t)
deallocate(sh)
deallocate(rh)
deallocate(tca)
deallocate(cca)
deallocate(mr_lsliq)
deallocate(mr_lsice)
deallocate(mr_ccliq)
deallocate(mr_ccice)
deallocate(mr_ozone)
deallocate(fl_lsrain)
deallocate(fl_lssnow)
deallocate(fl_lsgrpl)
deallocate(fl_ccrain)
deallocate(fl_ccsnow)
deallocate(reff)
deallocate(dtau_s)
deallocate(dtau_c)
deallocate(dem_s)
deallocate(dem_c)

allocate(t_surf(n_profile_full),p_surf(n_profile_full),hgt_surf(n_profile_full))
allocate(cosp_sunlit(n_profile_full))

allocate(cosp_lon(n_profile_full))
allocate(cosp_lat(n_profile_full))

do j=1, n_profile_full
  t_surf(j) = skt(j)
  p_surf(j) = ph(j,1)
  hgt_surf(j) = 0.0
  cosp_sunlit(j) = sunlit(j)
  cosp_lon(j) = lon(j)
  cosp_lat(j) = lat(j)
end do

deallocate(ph)

d_mass=1.0
frac_cloud_water=1.0
frac_cloud_ice=1.0
clw_sub_full=0.1

! For rain taken from [namelist:run_precip]
x1r = 2.2000e-1
x2r = 2.2000
! For graupel taken from mphys_psd_mod.f90
x1g = 7.9e9
x2g = -2.58
x4g = 0.0

d_p_full_levels => p_full_levels
d_p_half_levels => p_half_levels
d_hgt_full_levels => hgt_full_levels
d_hgt_half_levels => hgt_half_levels
d_d_mass => d_mass
d_t_n => t_n
d_q_n => q_n
d_w_cloud => w_cloud
d_condensed_mix_ratio_water => condensed_mix_ratio_water
d_condensed_mix_ratio_ice => condensed_mix_ratio_ice
d_cosp_crain => cosp_crain
d_cosp_csnow => cosp_csnow
d_condensed_re_water => condensed_re_water
d_condensed_re_ice => condensed_re_ice
d_cloud_extinction => cloud_extinction
d_cloud_absorptivity => cloud_absorptivity
d_frac_cloud_water => frac_cloud_water
d_frac_cloud_ice => frac_cloud_ice
d_clw_sub_full => clw_sub_full


allocate(cosp_calipso_low_level_cl_mask(n_profile_full)) ! 2321
allocate(cosp_calipso_mid_level_cl_mask(n_profile_full)) ! 2322
allocate(cosp_calipso_high_level_cl_mask(n_profile_full)) ! 2323

allocate(cosp_calipso_low_level_cl(n_profile_full)) ! 2344
allocate(cosp_calipso_mid_level_cl(n_profile_full)) ! 2345
allocate(cosp_calipso_high_level_cl(n_profile_full)) ! 2346

allocate(cosp_cloud_weights(n_profile_full)) ! 2330
allocate(cosp_ctp_tau_histogram(n_profile_full*7*7)) ! 2337

allocate(cosp_calipso_tot_backscatter(n_profile_full*ncolumns*nlevels)) ! 2341

allocate(cosp_calipso_cfad_sr_40(n_profile_full*n_backscatter_bins*n_cloudsat_levels)) ! 2370

allocate(cosp_calipso_cf_40_mask(n_profile_full*n_cloudsat_levels)) ! 2325
allocate(cosp_calipso_cf_40_liq(n_profile_full*n_cloudsat_levels)) ! 2473
allocate(cosp_calipso_cf_40_ice(n_profile_full*n_cloudsat_levels)) ! 2474
allocate(cosp_calipso_cf_40_undet(n_profile_full*n_cloudsat_levels)) ! 2475

! 2321
cosp_diag%cosp_calipso_low_level_cl_mask(1:n_profile_full) &
=> cosp_calipso_low_level_cl_mask(1:n_profile_full)

! 2322
cosp_diag%cosp_calipso_mid_level_cl_mask(1:n_profile_full) &
=> cosp_calipso_mid_level_cl_mask(1:n_profile_full)

! 2323
cosp_diag%cosp_calipso_high_level_cl_mask(1:n_profile_full) &
=> cosp_calipso_high_level_cl_mask(1:n_profile_full)

! 2325
if (l_profile_last) then
cosp_diag%cosp_calipso_cf_40_mask(1:n_cloudsat_levels,1:n_profile_full) &
=> cosp_calipso_cf_40_mask(1:n_profile_full*n_cloudsat_levels)
else
cosp_diag%cosp_calipso_cf_40_mask(1:n_profile_full,1:n_cloudsat_levels) &
=> cosp_calipso_cf_40_mask(1:n_profile_full*n_cloudsat_levels)
end if ! l_profile_last

! 2330
cosp_diag%cosp_cloud_weights(1:n_profile_full) &
=> cosp_cloud_weights(1:n_profile_full)

! 2337
if (l_profile_last) then
cosp_diag%cosp_ctp_tau_histogram(1:n_isccp_tau_bins,1:n_isccp_pressure_bins,1:n_profile_full) &
=> cosp_ctp_tau_histogram(1:n_profile_full*n_isccp_tau_bins*n_isccp_pressure_bins)
else
cosp_diag%cosp_ctp_tau_histogram(1:n_profile_full,1:n_isccp_tau_bins,1:n_isccp_pressure_bins) &
=> cosp_ctp_tau_histogram(1:n_profile_full*n_isccp_tau_bins*n_isccp_pressure_bins)
end if ! l_profile_last

! 2341
if (l_profile_last) then
cosp_diag%cosp_calipso_tot_backscatter(1:nlevels,1:ncolumns,1:n_profile_full) &
=> cosp_calipso_tot_backscatter(1:n_profile_full*ncolumns*nlevels)
else
cosp_diag%cosp_calipso_tot_backscatter(1:n_profile_full,1:nlevels,1:ncolumns) &
=> cosp_calipso_tot_backscatter(1:n_profile_full*ncolumns*nlevels)
end if ! l_profile_last

! 2344
cosp_diag%cosp_calipso_low_level_cl(1:n_profile_full) &
=> cosp_calipso_low_level_cl(1:n_profile_full)

! 2345
cosp_diag%cosp_calipso_mid_level_cl(1:n_profile_full) &
=> cosp_calipso_mid_level_cl(1:n_profile_full)

! 2346
cosp_diag%cosp_calipso_high_level_cl(1:n_profile_full) &
=> cosp_calipso_high_level_cl(1:n_profile_full)

! 2370
if (l_profile_last) then
cosp_diag%cosp_calipso_cfad_sr_40(1:n_backscatter_bins,1:n_cloudsat_levels,1:n_profile_full) &
=> cosp_calipso_cfad_sr_40(1:n_profile_full*n_backscatter_bins*n_cloudsat_levels)
else
cosp_diag%cosp_calipso_cfad_sr_40(1:n_profile_full,1:n_backscatter_bins,1:n_cloudsat_levels) &
=> cosp_calipso_cfad_sr_40(1:n_profile_full*n_backscatter_bins*n_cloudsat_levels)
end if ! l_profile_last

! 2473
if (l_profile_last) then
cosp_diag%cosp_calipso_cf_40_liq(1:n_cloudsat_levels,1:n_profile_full) &
=> cosp_calipso_cf_40_liq(1:n_profile_full*n_cloudsat_levels)
else
cosp_diag%cosp_calipso_cf_40_liq(1:n_profile_full,1:n_cloudsat_levels) &
=> cosp_calipso_cf_40_liq(1:n_profile_full*n_cloudsat_levels)
end if ! l_profile_last

! 2474
if (l_profile_last) then
cosp_diag%cosp_calipso_cf_40_ice(1:n_cloudsat_levels,1:n_profile_full) &
=> cosp_calipso_cf_40_ice(1:n_profile_full*n_cloudsat_levels)
else
cosp_diag%cosp_calipso_cf_40_ice(1:n_profile_full,1:n_cloudsat_levels) &
=> cosp_calipso_cf_40_ice(1:n_profile_full*n_cloudsat_levels)
end if ! l_profile_last

! 2475
if (l_profile_last) then
cosp_diag%cosp_calipso_cf_40_undet(1:n_cloudsat_levels,1:n_profile_full) &
=> cosp_calipso_cf_40_undet(1:n_profile_full*n_cloudsat_levels)
else
cosp_diag%cosp_calipso_cf_40_undet(1:n_profile_full,1:n_cloudsat_levels) &
=> cosp_calipso_cf_40_undet(1:n_profile_full*n_cloudsat_levels)
end if ! l_profile_last

call cosp( nlevels, &
           n_profile_list, &
           ncolumns, &
           nclds, &
           ncldy, &
           d_p_full_levels, &
           d_p_half_levels, &
           d_hgt_full_levels, &
           d_hgt_half_levels, &
           d_d_mass, &
           d_t_n, &
           d_q_n, &
           d_w_cloud, &
           d_condensed_mix_ratio_water, &
           d_condensed_mix_ratio_ice, &
           d_cosp_crain, &
           d_cosp_csnow, &
           d_condensed_re_water, &
           d_condensed_re_ice, &
           d_cloud_extinction, &
           d_cloud_absorptivity, &
           d_frac_cloud_water, &
           d_frac_cloud_ice, &
           d_clw_sub_full, &
           t_surf, &
           p_surf, &
           hgt_surf, &
           cosp_sunlit, &
           x1r, &
           x1g, &
           x2r, &
           x2g, &
           x4g, &
           cosp_diag = cosp_diag, &
           profile_list = profile_list, &
           l_profile_last = l_profile_last, &
           cosp_out_ext = cosp_out_ext )

deallocate(ncldy)
deallocate(p_full_levels)
deallocate(p_half_levels)
deallocate(hgt_full_levels)
deallocate(hgt_half_levels)
deallocate(d_mass)
deallocate(t_n)
deallocate(q_n)
deallocate(w_cloud)
deallocate(condensed_mix_ratio_water)
deallocate(condensed_mix_ratio_ice)
deallocate(cosp_crain)
deallocate(cosp_csnow)
deallocate(condensed_re_water)
deallocate(condensed_re_ice)
deallocate(cloud_extinction)
deallocate(cloud_absorptivity)
deallocate(frac_cloud_water)
deallocate(frac_cloud_ice)
deallocate(clw_sub_full)
deallocate(t_surf, p_surf, hgt_surf)
deallocate(cosp_sunlit)

ierr=0
call get_free_unit(ierr,lun)
if (ierr /= 0) stop
open(lun,file='output.txt')

do l=1,n_profile_full

  if (l==10*int(l/10)+1) then

    write(lun,*) 'profile: ',l

    ! 2330
    write(lun,*) 'cloud_weights (2330)'
    write(lun,*) cosp_diag%cosp_cloud_weights(l)

    ! 2337
    write(lun,*) 'ctp_tau_histogram (2337)'
    if (l_profile_last) then
      do i=1, n_isccp_pressure_bins
        write(lun,'(10f5.1)') cosp_diag%cosp_ctp_tau_histogram(1:n_isccp_tau_bins, i, l)
      end do
    else
      do i=1, n_isccp_pressure_bins
        write(lun,'(10f5.1)') cosp_diag%cosp_ctp_tau_histogram(l, 1:n_isccp_tau_bins, i)
      end do
    end if ! l_profile_last

    ! 2341
    write(lun,*) 'calipso_tot_backscatter (2341)'
    if (l_profile_last) then
      do i=1, nlevels
        write(lun,'(i4,1000e16.6)') i, cosp_diag%cosp_calipso_tot_backscatter(i, 1:ncolumns, l)
      end do
    else
      do i=1, nlevels
        write(lun,'(i4,1000e16.6)') i, cosp_diag%cosp_calipso_tot_backscatter(l, i, 1:ncolumns)
      end do
    end if ! l_profile_last

    ! 2321-2323 and 2344-2346
    write(lun,*) 'low/mid/high-level cloud mask/value (2321-2323, 2344-2346)'
    write(lun,'(10f5.1)') cosp_diag%cosp_calipso_low_level_cl_mask(l),cosp_diag%cosp_calipso_low_level_cl(l), &
                          cosp_diag%cosp_calipso_mid_level_cl_mask(l),cosp_diag%cosp_calipso_mid_level_cl(l), &
                          cosp_diag%cosp_calipso_high_level_cl_mask(l),cosp_diag%cosp_calipso_high_level_cl(l)

    ! 2325, 2473-2475 and 2370
    if (l_profile_last) then
      write(lun,*) 'on 40 levels: '
      write(lun,*) 'cf mask/liq/ice/undet and cfad_sr_backscatter_bins (2325, 2473-2475 and 2370)'
      do i=1,n_cloudsat_levels
        write(lun,'(i4,100f5.1)') i,cosp_diag%cosp_calipso_cf_40_mask(i,l),cosp_diag%cosp_calipso_cf_40_liq(i,l), &
                                  cosp_diag%cosp_calipso_cf_40_ice(i,l),cosp_diag%cosp_calipso_cf_40_undet(i,l), &
                                  cosp_diag%cosp_calipso_cfad_sr_40(1:n_backscatter_bins, i, l)
      end do
    else
      do i=1,n_cloudsat_levels
        write(lun,'(i4,100f5.1)') i,cosp_diag%cosp_calipso_cf_40_mask(l,i),cosp_diag%cosp_calipso_cf_40_liq(l,i), &
                                  cosp_diag%cosp_calipso_cf_40_ice(l,i),cosp_diag%cosp_calipso_cf_40_undet(l,i), &
                                  cosp_diag%cosp_calipso_cfad_sr_40(l, 1:n_backscatter_bins, i)
    end do
    end if ! l_profile_last

  end if ! profile selection

end do

close(lun)


foutput='cosp_out.nc'
call write_cosp2_output(n_profile_full, ncolumns, nlevels, &
                        zlev(1,nlevels:1:-1), &
                        cosp_lon, cosp_lat, cosp_out_ext, foutput)

deallocate(zlev)

deallocate(cosp_lon)
deallocate(cosp_lat)

deallocate(cosp_calipso_low_level_cl_mask) ! 2321
deallocate(cosp_calipso_mid_level_cl_mask) ! 2322
deallocate(cosp_calipso_high_level_cl_mask) ! 2323

deallocate(cosp_calipso_low_level_cl) ! 2344
deallocate(cosp_calipso_mid_level_cl) ! 2345
deallocate(cosp_calipso_high_level_cl) ! 2346

deallocate(cosp_calipso_cf_40_mask) ! 2325

deallocate(cosp_calipso_cfad_sr_40) ! 2370

deallocate(cosp_calipso_cf_40_liq) ! 2473
deallocate(cosp_calipso_cf_40_ice) ! 2474
deallocate(cosp_calipso_cf_40_undet) ! 2475

end program cosp_driver
