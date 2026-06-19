! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module holds parameter settings related to COSP2
!
!------------------------------------------------------------------------------

module cosp_input_mod

use cosp_kinds, only: wp
use errormessagelength_mod, only: errormessagelength

implicit none

!- Control variables set in run_cosp namelist
integer :: i_cosp_version    = 0
logical :: cosp_cloudsat_sim = .false.
logical :: cosp_lidar_sim    = .false.
logical :: cosp_isccp_sim    = .false.
logical :: cosp_misr_sim     = .false.
logical :: cosp_modis_sim    = .false.
logical :: cosp_rttov_sim    = .false.
logical :: cosp_use_vgrid    = .false.

!- Defaults for vertical grid if cosp_use_vgrid is set to true
integer, parameter :: n_cloudsat_levels = 40
integer :: cosp_nlr        = n_cloudsat_levels
logical :: cosp_csat_vgrid = .true.

!- Other control variables
integer,parameter :: cosp_ncolumns_max = 64
integer,parameter :: n_hydro = 9

!- Inputs related to radar simulations
real(wp)    :: cosp_radar_freq     = 94.0
integer :: cosp_surface_radar  = 0
integer :: cosp_use_mie_tables = 0
integer :: cosp_use_gas_abs    = 1
integer :: cosp_do_ray         = 0
integer :: cosp_melt_lay       = 0
real(wp) :: cosp_k2            = -1
logical :: cosp_use_reff       = .true.
logical :: cosp_use_precipitation_fluxes = .false.
character(len=64) :: cloudsat_micro_scheme = 'um'
integer :: cosp_qb_dtype(n_hydro),cosp_qb_phase(n_hydro)
real(wp) :: cosp_qb_dmin(n_hydro),cosp_qb_dmax(n_hydro), cosp_qb_apm(n_hydro), &
        cosp_qb_bpm(n_hydro),cosp_qb_rho(n_hydro), cosp_qb_p1(n_hydro), &
        cosp_qb_p2(n_hydro),cosp_qb_p3(n_hydro)
!- Inputs related to lidar simulations
integer :: cosp_nprmts_max_hydro = 12
integer :: cosp_naero            = 1
integer :: cosp_nprmts_max_aero  = 1
integer :: cosp_lidar_ice_type   = 0
!- Inputs related to isccp simulator
integer :: cosp_overlap                   = 3
integer :: cosp_isccp_topheight           = 1
integer :: cosp_isccp_topheight_direction = 2
real(wp)    :: cosp_emsfc_lw                  = 0.99
!- Inputs related to rttov
integer :: cosp_platform   = 1
integer :: cosp_satellite  = 15
integer :: cosp_instrument = 0
integer,parameter :: cosp_nchannels = 8
integer :: cosp_channels(cosp_nchannels) = [1,3,5,6,8,10,11,13]
real(wp) :: cosp_surfem(cosp_nchannels) = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
real(wp) :: cosp_zenang = 50.0
real(wp) :: cosp_co2    = 5.241e-04
real(wp) :: cosp_ch4    = 9.139e-07
real(wp) :: cosp_n2o    = 4.665e-07
real(wp) :: cosp_co     = 2.098e-07

! HCLASS table for quickbeam.
! CP is not used in the version of quickbeam included in COSP
!                  lsl     lsi   lsr     lsa   cvl     cvi   cvr     cvs     lsg
data cosp_qb_dtype/ 1,      1,    1,     -1,    1,      1,    1,      1,     -1/
data cosp_qb_phase/ 0,      1,    0,      1,    0,      1,    0,      1,      1/
data cosp_qb_dmin/ -1,     -1,   -1,     -1,   -1,     -1,   -1,     -1,     -1/
data cosp_qb_dmax/ -1,     -1,   -1,     -1,   -1,     -1,   -1,     -1,     -1/
data cosp_qb_apm/  -1, 0.0257,   -1, 0.0185,   -1, 0.0185,   -1, 0.0444,  261.8/
data cosp_qb_bpm/  -1,    2.0,   -1,    1.9,   -1,    1.9,   -1,    2.1,      3/
data cosp_qb_rho/1000,     -1, 1000,     -1, 1000,     -1, 1000,     -1,     -1/
data cosp_qb_p1/   -1,     -1,   -1,     -1,   -1,     -1,   -1,     -1,     -1/
data cosp_qb_p2/   10,     40, 1000,     40,   10,     40, 1000,    120,   1000/
data cosp_qb_p3/    3,      1,    1,      1,    3,      1,    1,      1,    3.5/

end module cosp_input_mod
