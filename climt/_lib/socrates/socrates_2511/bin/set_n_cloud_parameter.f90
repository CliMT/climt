! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Function to set number of cloudy parameters.
!
! Method:
!       Straightforward
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------
FUNCTION set_n_cloud_parameter(i_scheme, i_component, n_phase_term)


  USE rad_pcf, ONLY: ip_ackerman_stephens, ip_clcmp_cnv_ice,            &
                     ip_clcmp_cnv_water, ip_clcmp_st_ice,               &
                     ip_clcmp_st_water, ip_drop_pade_2,                 &
                     ip_drop_pade_2_phf, ip_ice_adt, ip_ice_adt_10,     &
                     ip_ice_baran, ip_ice_fu_ir,ip_ice_fu_phf,          &
                     ip_ice_fu_solar, ip_ice_iwc_only, ip_ice_t_iwc,    &
                     ip_ps_size_phf, ip_slingo_schr_ice_phf,            &
                     ip_slingo_schr_phf, ip_slingo_schrecker,           &
                     ip_slingo_schrecker_ice, ip_sun_shine_vn2_ir,      &
                     ip_sun_shine_vn2_vis, ip_ice_pade_2_phf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      i_scheme                                                          &
!       Parametrization scheme
    , i_component                                                       &
!       Component in cloud
    , n_phase_term
!       Number of terms in the phase function

  INTEGER ::                                                            &
      set_n_cloud_parameter
!       Returned number of coefficients in parametrization

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_N_CLOUD_PARAMETER'


  IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

  IF ( (i_component == ip_clcmp_st_water).OR.                           &
       (i_component == ip_clcmp_cnv_water) ) THEN

    IF (i_scheme == ip_slingo_schrecker) THEN
      set_n_cloud_parameter=6
    ELSE IF (i_scheme == ip_ackerman_stephens) THEN
      set_n_cloud_parameter=9
    ELSE IF (i_scheme == ip_drop_pade_2) THEN
      set_n_cloud_parameter=16
    ELSE IF (i_scheme == ip_slingo_schr_phf) THEN
      set_n_cloud_parameter=4+2*n_phase_term
    ELSE IF (i_scheme == IP_drop_pade_2_PHF) THEN
      set_n_cloud_parameter=11+5*n_phase_term
    ELSE IF (i_scheme == IP_ps_size_PHF) THEN
      set_n_cloud_parameter=11+5*n_phase_term
    END IF

  ELSE IF ( (i_component == ip_clcmp_st_ice).OR.                        &
            (i_component == ip_clcmp_cnv_ice) ) THEN

    IF (i_scheme == ip_slingo_schrecker_ice) THEN
      set_n_cloud_parameter=6
    ELSE IF (i_scheme == ip_ice_adt) THEN
      set_n_cloud_parameter=30
    ELSE IF (i_scheme == ip_ice_adt_10) THEN
      set_n_cloud_parameter=36
    ELSE IF (i_scheme == ip_sun_shine_vn2_vis) THEN
      set_n_cloud_parameter=6
    ELSE IF (i_scheme == ip_sun_shine_vn2_ir) THEN
      set_n_cloud_parameter=0
    ELSE IF (i_scheme == ip_ice_fu_solar) THEN
      set_n_cloud_parameter=14
    ELSE IF (i_scheme == ip_ice_fu_ir) THEN
      set_n_cloud_parameter=10
    ELSE IF (i_scheme == ip_slingo_schr_ice_phf) THEN
      set_n_cloud_parameter=4+2*n_phase_term
    ELSE IF (i_scheme == ip_ice_fu_phf) THEN
      set_n_cloud_parameter=9+5*n_phase_term
    ELSE IF (i_scheme == ip_ice_t_iwc) THEN
      set_n_cloud_parameter=9
    ELSE IF (i_scheme == ip_ice_iwc_only) THEN
      set_n_cloud_parameter=4
    ELSE IF (i_scheme == ip_ice_baran) THEN
      set_n_cloud_parameter=5
    ELSE IF (i_scheme == ip_ice_pade_2_phf) THEN
      set_n_cloud_parameter=11+5*n_phase_term
    END IF

  END IF


  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
END FUNCTION set_n_cloud_parameter
