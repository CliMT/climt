! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module holds properties related to the different hydrometeor types
! used in COSP2
!
!------------------------------------------------------------------------------

module cosp2_constants_mod
use mod_quickbeam_optics, only: size_distribution
use cosp_kinds, only: wp
use cosp_input_mod, only:cosp_qb_dtype, cosp_qb_phase, cosp_qb_dmin, &
  cosp_qb_dmax, cosp_qb_apm, cosp_qb_bpm, cosp_qb_rho, cosp_qb_p1, &
  cosp_qb_p2, cosp_qb_p3
use mod_cosp_config, only: n_hydro
implicit none

character(len=16) :: cosp_version='COSP v2'
character(len=*), parameter, private :: ModuleName='COSP2_CONSTANTS_MOD'

! Indices to address arrays of LS and CONV hydrometeors
integer,parameter :: i_lscliq = 1
integer,parameter :: i_lscice = 2
integer,parameter :: i_lsrain = 3
integer,parameter :: i_lsiagg = 4
integer,parameter :: i_cvcliq = 5
integer,parameter :: i_cvcice = 6
integer,parameter :: i_cvrain = 7
integer,parameter :: i_cvsnow = 8
integer,parameter :: i_lsgrpl = 9

! Stratiform and convective clouds in frac_out
integer, parameter :: i_lsc = 1 ! Large-scale clouds
integer, parameter :: i_cvc = 2 ! Convective clouds
! Microphysical settings
real(wp) :: n_ax(n_hydro)
real(wp) :: n_bx(n_hydro)
real(wp) :: alpha_x(n_hydro)
real(wp) :: c_x(n_hydro)
real(wp) :: d_x(n_hydro)
real(wp) :: g_x(n_hydro)
real(wp) :: a_x(n_hydro)
real(wp) :: b_x(n_hydro)
real(wp) :: gamma_1(n_hydro)
real(wp) :: gamma_2(n_hydro)
real(wp) :: gamma_3(n_hydro)
real(wp) :: gamma_4(n_hydro)
real(wp) :: x_1(n_hydro)
real(wp) :: x_2(n_hydro)
real(wp) :: x_3(n_hydro)
real(wp) :: x_4(n_hydro)

real,parameter :: icarus_missing_value = -1.0e30

! Quickbeam size distribution
type(size_distribution) :: cosp_qb_dist

contains

subroutine cosp_microphys_init(x1r,x1g,x2r,x2g,x4g)
implicit none
real(wp), intent(in) :: x1r
real(wp), intent(in) :: x1g
real(wp), intent(in) :: x2r
real(wp), intent(in) :: x2g
real(wp), intent(in) :: x4g
! Settings for the precipitation flux to mixing ratio conversion
!          LSL      LSI      LSR           LSA      CVL      CVI      CVR          CVS          LSG
n_ax = [   -1.0_wp, -1.0_wp,     0.22_wp,  -1.0_wp, -1.0_wp, -1.0_wp,     0.22_wp,    4.0e6_wp,      -1.0_wp]
n_bx = [   -1.0_wp, -1.0_wp,     2.20_wp,  -1.0_wp, -1.0_wp, -1.0_wp,     2.20_wp,      0.0_wp,      -4.0_wp]
alpha_x = [-1.0_wp, -1.0_wp,      0.0_wp,  -1.0_wp, -1.0_wp, -1.0_wp,      0.0_wp,      0.0_wp,       2.5_wp]
c_x = [    -1.0_wp, -1.0_wp,    386.8_wp,  -1.0_wp, -1.0_wp, -1.0_wp,    386.8_wp,     14.3_wp,     253.0_wp]
d_x = [    -1.0_wp, -1.0_wp,     0.67_wp,  -1.0_wp, -1.0_wp, -1.0_wp,     0.67_wp,    0.416_wp,     0.734_wp]
g_x = [    -1.0_wp, -1.0_wp,      0.4_wp,  -1.0_wp, -1.0_wp, -1.0_wp,      0.4_wp,      0.4_wp,       0.4_wp]
a_x = [    -1.0_wp, -1.0_wp,    523.6_wp,  -1.0_wp, -1.0_wp, -1.0_wp,    523.6_wp,   0.0444_wp,     261.8_wp]
b_x = [    -1.0_wp, -1.0_wp,      3.0_wp,  -1.0_wp, -1.0_wp, -1.0_wp,      3.0_wp,      2.1_wp,       3.0_wp]
gamma_1 = [-1.0_wp, -1.0_wp, 14.78119_wp,  -1.0_wp, -1.0_wp, -1.0_wp, 14.78119_wp, 3.382827_wp, 1120.6197_wp]
gamma_2 = [-1.0_wp, -1.0_wp,      6.0_wp,  -1.0_wp, -1.0_wp, -1.0_wp,      6.0_wp, 2.197659_wp,  287.8853_wp]
gamma_3 = [-1.0_wp, -1.0_wp,      2.0_wp,  -1.0_wp, -1.0_wp, -1.0_wp,      2.0_wp,      2.0_wp,  52.34278_wp]
gamma_4 = [-1.0_wp, -1.0_wp,      6.0_wp,  -1.0_wp, -1.0_wp, -1.0_wp,      6.0_wp,      6.0_wp,  287.8853_wp]
x_1 = [    -1.0_wp, -1.0_wp,      x1r,     -1.0_wp, -1.0_wp, -1.0_wp,      x1r,    3.382827_wp,       x1g]
x_2 = [    -1.0_wp, -1.0_wp,      x2r,     -1.0_wp, -1.0_wp, -1.0_wp,      x2r,    2.197659_wp,       x2g]
x_3 = [    -1.0_wp, -1.0_wp,      0.0_wp,  -1.0_wp, -1.0_wp, -1.0_wp,      0.0_wp,      2.0_wp,       0.0_wp]
x_4 = [    -1.0_wp, -1.0_wp,      0.0_wp,  -1.0_wp, -1.0_wp, -1.0_wp,      0.0_wp,      6.0_wp,       x4g]
end subroutine cosp_microphys_init

subroutine cosp_hydroclass_init(quickbeam_size_dist)
implicit none
!-----Input arguments
! Quickbeam size distribution
type(size_distribution), intent(out) :: quickbeam_size_dist
!-----Local variables
character(len=*), parameter :: RoutineName='COSP_HYDROCLASS_INIT'

! Initialise all elements
quickbeam_size_dist%p1 = -1.0
quickbeam_size_dist%p2 = -1.0
quickbeam_size_dist%p3 = -1.0
quickbeam_size_dist%dmin = -1.0
quickbeam_size_dist%dmax = -1.0
quickbeam_size_dist%apm = -1.0
quickbeam_size_dist%bpm = -1.0
quickbeam_size_dist%rho = -1.0
quickbeam_size_dist%dtype = -1
quickbeam_size_dist%phase = -1
! Fill in with hydrometeor parameters
quickbeam_size_dist%p1(1:n_hydro) = cosp_qb_p1
quickbeam_size_dist%p2(1:n_hydro) = cosp_qb_p2
quickbeam_size_dist%p3(1:n_hydro) = cosp_qb_p3
quickbeam_size_dist%dmin(1:n_hydro) = cosp_qb_dmin
quickbeam_size_dist%dmax(1:n_hydro) = cosp_qb_dmax
quickbeam_size_dist%apm(1:n_hydro) = cosp_qb_apm
quickbeam_size_dist%bpm(1:n_hydro) = cosp_qb_bpm
quickbeam_size_dist%rho(1:n_hydro) = cosp_qb_rho
quickbeam_size_dist%dtype(1:n_hydro) = cosp_qb_dtype
quickbeam_size_dist%phase(1:n_hydro) = cosp_qb_phase

return
end subroutine cosp_hydroclass_init
end module cosp2_constants_mod
