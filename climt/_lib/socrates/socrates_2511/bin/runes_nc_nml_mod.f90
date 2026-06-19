! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! To run the offline driver (runes_nc) it requires the relevant parts
! of the LFRic configuration namelist
!
! Description: The relevant namelists are defined in this module and
!              given default values which may be overwritten 
!              with actual values
!
! ---------------------------------------------------------------------

MODULE runes_nc_nml_mod

USE realtype_rd, ONLY : RealExt

IMPLICIT NONE

! Namelist elements potentially required for offline driver

!AEROSOL
LOGICAL :: easyaerosol_lw
LOGICAL :: easyaerosol_sw
LOGICAL :: l_radaer
LOGICAL :: sulphuric_strat_climatology

REAL(RealExt) :: sulphuric_strat_column

! RADIATION
CHARACTER(LEN=200) :: cloud_entrapment
CHARACTER(LEN=200) :: cloud_inhomogeneity
CHARACTER(LEN=200) :: cloud_overlap
CHARACTER(LEN=200) :: cloud_representation
CHARACTER(LEN=200) :: droplet_effective_radius
CHARACTER(LEN=200) :: mcica_data_file
CHARACTER(LEN=200) :: scatter_method_lw
CHARACTER(LEN=200) :: spectral_file_lw
CHARACTER(LEN=200) :: spectral_file_sw
CHARACTER(LEN=200) :: topography

INTEGER :: i_cloud_ice_type_lw
INTEGER :: i_cloud_ice_type_sw
INTEGER :: i_cloud_liq_type_lw
INTEGER :: i_cloud_liq_type_sw

LOGICAL :: l_rayleigh_sw
LOGICAL :: l_trans_zen_correction

REAL(RealExt) :: cloud_vertical_decorr
REAL(RealExt) :: constant_droplet_effective_radius
REAL(RealExt) :: liu_aparam
REAL(RealExt) :: liu_bparam

! RADIATIVE GASES
CHARACTER(LEN=200) :: cfc113_rad_opt
CHARACTER(LEN=200) :: cfc11_rad_opt
CHARACTER(LEN=200) :: cfc12_rad_opt
CHARACTER(LEN=200) :: ch4_rad_opt
CHARACTER(LEN=200) :: co2_rad_opt
CHARACTER(LEN=200) :: co_rad_opt
CHARACTER(LEN=200) :: cs_rad_opt
CHARACTER(LEN=200) :: h2_rad_opt
CHARACTER(LEN=200) :: h2o_rad_opt
CHARACTER(LEN=200) :: hcfc22_rad_opt
CHARACTER(LEN=200) :: hcn_rad_opt
CHARACTER(LEN=200) :: he_rad_opt
CHARACTER(LEN=200) :: hfc134a_rad_opt
CHARACTER(LEN=200) :: k_rad_opt
CHARACTER(LEN=200) :: li_rad_opt
CHARACTER(LEN=200) :: n2_rad_opt
CHARACTER(LEN=200) :: n2o_rad_opt
CHARACTER(LEN=200) :: na_rad_opt
CHARACTER(LEN=200) :: nh3_rad_opt
CHARACTER(LEN=200) :: o2_rad_opt
CHARACTER(LEN=200) :: o3_rad_opt
CHARACTER(LEN=200) :: rb_rad_opt
CHARACTER(LEN=200) :: so2_rad_opt
CHARACTER(LEN=200) :: tio_rad_opt
CHARACTER(LEN=200) :: vo_rad_opt

REAL(RealExt) :: cfc113_mix_ratio
REAL(RealExt) :: cfc11_mix_ratio
REAL(RealExt) :: cfc12_mix_ratio
REAL(RealExt) :: ch4_mix_ratio
REAL(RealExt) :: co2_mix_ratio
REAL(RealExt) :: co_mix_ratio
REAL(RealExt) :: cs_mix_ratio
REAL(RealExt) :: h2_mix_ratio
REAL(RealExt) :: h2o_mix_ratio
REAL(RealExt) :: hcfc22_mix_ratio
REAL(RealExt) :: hcn_mix_ratio
REAL(RealExt) :: he_mix_ratio
REAL(RealExt) :: hfc134a_mix_ratio
REAL(RealExt) :: k_mix_ratio
REAL(RealExt) :: li_mix_ratio
REAL(RealExt) :: n2_mix_ratio
REAL(RealExt) :: n2o_mix_ratio
REAL(RealExt) :: na_mix_ratio
REAL(RealExt) :: nh3_mix_ratio
REAL(RealExt) :: o2_mix_ratio
REAL(RealExt) :: o3_mix_ratio
REAL(RealExt) :: rb_mix_ratio
REAL(RealExt) :: so2_mix_ratio
REAL(RealExt) :: tio_mix_ratio
REAL(RealExt) :: vo_mix_ratio

NAMELIST /aerosol/ &
easyaerosol_lw, &
easyaerosol_sw, &
l_radaer, &
sulphuric_strat_climatology, &
sulphuric_strat_column

NAMELIST /radiation/ &
cloud_entrapment, &
cloud_inhomogeneity, &
cloud_overlap, &
cloud_representation, &
droplet_effective_radius, &
mcica_data_file, &
scatter_method_lw, &
spectral_file_lw, &
spectral_file_sw, &
topography, &
i_cloud_ice_type_lw, &
i_cloud_ice_type_sw, &
i_cloud_liq_type_lw, &
i_cloud_liq_type_sw, &
l_rayleigh_sw, &
l_trans_zen_correction, &
cloud_vertical_decorr, &
constant_droplet_effective_radius, &
liu_aparam, &
liu_bparam

NAMELIST /radiative_gases/ &
cfc113_rad_opt, &
cfc11_rad_opt, &
cfc12_rad_opt, &
ch4_rad_opt, &
co2_rad_opt, &
co_rad_opt, &
cs_rad_opt, &
h2_rad_opt, &
h2o_rad_opt, &
hcfc22_rad_opt, &
hcn_rad_opt, &
he_rad_opt, &
hfc134a_rad_opt, &
k_rad_opt, &
li_rad_opt, &
n2_rad_opt, &
n2o_rad_opt, &
na_rad_opt, &
nh3_rad_opt, &
o2_rad_opt, &
o3_rad_opt, &
rb_rad_opt, &
so2_rad_opt, &
tio_rad_opt, &
vo_rad_opt, &
cfc113_mix_ratio, &
cfc11_mix_ratio, &
cfc12_mix_ratio, &
ch4_mix_ratio, &
co2_mix_ratio, &
co_mix_ratio, &
cs_mix_ratio, &
h2_mix_ratio, &
h2o_mix_ratio, &
hcfc22_mix_ratio, &
hcn_mix_ratio, &
he_mix_ratio, &
hfc134a_mix_ratio, &
k_mix_ratio, &
li_mix_ratio, &
n2_mix_ratio, &
n2o_mix_ratio, &
na_mix_ratio, &
nh3_mix_ratio, &
o2_mix_ratio, &
o3_mix_ratio, &
rb_mix_ratio, &
so2_mix_ratio, &
tio_mix_ratio, &
vo_mix_ratio


CONTAINS

! Set default values for the namelist entries
SUBROUTINE set_nml_default

IMPLICIT NONE

!Aerosol
easyaerosol_lw=.false.
easyaerosol_sw=.false.
l_radaer=.false.
sulphuric_strat_climatology=.false.
sulphuric_strat_column=0.0_RealExt

!Radiation

!Radiative gases
cfc113_rad_opt='off'
cfc11_rad_opt='off'
cfc12_rad_opt='off'
ch4_rad_opt='off'
co2_rad_opt='off'
co_rad_opt='off'
cs_rad_opt='off'
h2_rad_opt='off'
h2o_rad_opt='off'
hcfc22_rad_opt='off'
hcn_rad_opt='off'
he_rad_opt='off'
hfc134a_rad_opt='off'
k_rad_opt='off'
li_rad_opt='off'
n2_rad_opt='off'
n2o_rad_opt='off'
na_rad_opt='off'
nh3_rad_opt='off'
o2_rad_opt='off'
o3_rad_opt='off'
rb_rad_opt='off'
so2_rad_opt='off'
tio_rad_opt='off'
vo_rad_opt='off'
cfc113_mix_ratio=0.0_RealExt
cfc11_mix_ratio=0.0_RealExt
cfc12_mix_ratio=0.0_RealExt
ch4_mix_ratio=0.0_RealExt
co2_mix_ratio=0.0_RealExt
co_mix_ratio=0.0_RealExt
cs_mix_ratio=0.0_RealExt
h2_mix_ratio=0.0_RealExt
h2o_mix_ratio=0.0_RealExt
hcfc22_mix_ratio=0.0_RealExt
hcn_mix_ratio=0.0_RealExt
he_mix_ratio=0.0_RealExt
hfc134a_mix_ratio=0.0_RealExt
k_mix_ratio=0.0_RealExt
li_mix_ratio=0.0_RealExt
n2_mix_ratio=0.0_RealExt
n2o_mix_ratio=0.0_RealExt
na_mix_ratio=0.0_RealExt
nh3_mix_ratio=0.0_RealExt
o2_mix_ratio=0.0_RealExt
o3_mix_ratio=0.0_RealExt
rb_mix_ratio=0.0_RealExt
so2_mix_ratio=0.0_RealExt
tio_mix_ratio=0.0_RealExt
vo_mix_ratio=0.0_RealExt

END SUBROUTINE set_nml_default

END MODULE runes_nc_nml_mod
