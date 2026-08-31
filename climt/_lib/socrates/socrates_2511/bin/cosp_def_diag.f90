! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module defines a type (CospDiag) which holds a pointer for each
! COSP2 diagnostic
!
!------------------------------------------------------------------------------

module cosp_def_diag

use realtype_rd, only: RealExt

implicit none

type :: cospdiag

real(RealExt), pointer :: cosp_calipso_low_level_cl_mask(:) => null()
! COSP: MASK FOR CALIPSO LOW-LEVEL CF (was 2321)

real(RealExt), pointer :: cosp_calipso_mid_level_cl_mask(:) => null()
! COSP: MASK FOR CALIPSO MID-LEVEL CF (was 2322)

real(RealExt), pointer :: cosp_calipso_high_level_cl_mask(:) => null()
! COSP: MASK FOR CALIPSO HIGH-LEVEL CF (was 2323)

real(RealExt), pointer :: cosp_calipso_cf_40_mask(:,:) => null()
! COSP: MASK FOR CALIPSO CF 40 LVLS

real(RealExt), pointer :: cosp_cloud_weights(:) => null()
! COSP: ISCCP/MISR/MODIS CLOUD WEIGHTS (was 2330)

real(RealExt), pointer :: cosp_ctp_tau_histogram(:,:,:) => null()
! COSP: ISCCP CTP-TAU HISTOGRAM (was 2337)

real(RealExt), pointer :: cosp_calipso_tot_backscatter(:,:,:) => null()
! COSP: CALIPSO TOTAL BACKSCATTER (was 2341)

real(RealExt), pointer :: cosp_calipso_low_level_cl(:) => null()
! COSP: CALIPSO LOW-LEVEL CLOUD (was 2344)

real(RealExt), pointer :: cosp_calipso_mid_level_cl(:) => null()
! COSP: CALIPSO MID-LEVEL CLOUD (was 2345)

real(RealExt), pointer :: cosp_calipso_high_level_cl(:) => null()
! COSP: CALIPSO HIGH-LEVEL CLOUD (was 2346)

real(RealExt), pointer :: cosp_calipso_cfad_sr_40(:,:,:) => null()
! COSP: CALIPSO CFAD SR 40 CSAT LEVELS (was 2370)

real(RealExt), pointer :: cosp_calipso_cf_40_liq(:,:) => null()
! COSP: CALIPSO CF 40 LVLS (LIQ) (was 2473)

real(RealExt), pointer :: cosp_calipso_cf_40_ice(:,:) => null()
! COSP: CALIPSO CF 40 LVLS (ICE) (was 2474)

real(RealExt), pointer :: cosp_calipso_cf_40_undet(:,:) => null()
! COSP: CALIPSO CF 40 LVLS (UNDET) (was 2475)

end type cospdiag

end module cosp_def_diag
