! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT****************************************
!
! Subroutine to correct optical scattering properties for cloud
! inhomogeneities with ip_cairns.
!
! Method:
!   Straight forward. Applying Eq. (11) from Cairns et al., Journal of
!   Atmospheric Sciences, Volume 57, 2000.
!
!   Note that this routine modifies the optical propeties of a
!   single condensed component of the cloud.
!
!- -----------------------------------------------------------------------------
MODULE opt_prop_inhom_corr_cairns_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'OPT_PROP_INHOM_CORR_CAIRNS_MOD'
CONTAINS
SUBROUTINE opt_prop_inhom_corr_cairns(                                         &
      n_layer, n_cloud_top, n_cloud_profile, i_cloud_profile                   &
    , l_rescale, n_order_forward                                               &
    , rel_var_dens, k_ext_tot_cloud, k_ext_scat_cloud                          &
    , phase_fnc_cloud, forward_scatter_cloud                                   &
    , nd_profile, nd_layer, id_ct                                              &
    )


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                       &
      nd_profile                                                               &
!       Size allocated for profiles
    , nd_layer                                                                 &
!       Size allocated for layers
    , id_ct
!       Topmost declared cloudy layer

! Dummy variables.
  INTEGER, INTENT(IN) ::                                                       &
      n_layer                                                                  &
!       Number of layers
    , n_order_forward                                                          &
!       Order used in forming the forward scattering parameter
    , n_cloud_top                                                              &
!       Topmost cloudy layer
    , n_cloud_profile(id_ct: nd_layer)                                         &
!       Number of cloudy profiles
    , i_cloud_profile(nd_profile, id_ct: nd_layer)
!       Profiles containing clouds
  LOGICAL, INTENT(IN) ::                                                       &
      l_rescale
!       Delta-rescaling required

  REAL (RealK), INTENT(IN) ::                                                  &
      rel_var_dens(nd_profile, id_ct: nd_layer)
!       Relative variance of cloud density
  REAL (RealK), INTENT(INOUT) ::                                               &
      k_ext_scat_cloud(nd_profile, id_ct: nd_layer)                            &
!       Scattering extinction
    , k_ext_tot_cloud(nd_profile, id_ct: nd_layer)                             &
!       Total extinction
    , phase_fnc_cloud(nd_profile, id_ct: nd_layer)                             &
!       First order of cloudy phase function
    , forward_scatter_cloud(nd_profile, id_ct: nd_layer)
!       Cloudy forward scattering

! Local variables.
  INTEGER                                                                      &
      l                                                                        &
!       Loop variable
    , ll                                                                       &
!       Loop variable
    , i
!       Loop variable
  REAL (RealK) ::                                                              &
      asymmetry_process                                                        &
!       Asymmetry of current process
    , omega_process                                                            &
!       Single scattering albedo of current process
    , k_ext_tot_corr                                                           &
!       Correction to total extinction
    , omega_process_corr                                                       &
!       Correction to single scattering albedo
    , asymmetry_process_corr
!       Correction to asymmetry

  REAL (RealK), PARAMETER :: tiny_k=TINY(k_ext_scat_cloud)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'OPT_PROP_INHOM_CORR_CAIRNS'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Apply inhomogeneity corrections to optical scattering parameters

  DO i=n_cloud_top, n_layer
!CDIR NODEP
    DO ll=1, n_cloud_profile(i)
      l=i_cloud_profile(ll, i)

      k_ext_tot_corr=1.0_RealK/(1.0_RealK + rel_var_dens(l, i))
      IF (k_ext_scat_cloud(l, i) > tiny_k) THEN
        asymmetry_process=phase_fnc_cloud(l, i)/k_ext_scat_cloud(l, i)
        omega_process=k_ext_scat_cloud(l, i)/k_ext_tot_cloud(l, i)

        omega_process_corr=1.0_RealK/(1.0_RealK + rel_var_dens(l, i)           &
          *(1.0_RealK - omega_process))
        asymmetry_process_corr                                                 &
          =(1.0_RealK + rel_var_dens(l, i)*(1.0_RealK - omega_process))        &
          /(1.0_RealK + rel_var_dens(l, i)                                     &
          *(1.0_RealK - omega_process*asymmetry_process))

        k_ext_scat_cloud(l, i)=k_ext_scat_cloud(l, i)                          &
          *k_ext_tot_corr*omega_process_corr
        phase_fnc_cloud(l, i)=phase_fnc_cloud(l, i)                            &
          *k_ext_tot_corr*omega_process_corr*asymmetry_process_corr
      
        IF (l_rescale) THEN
          forward_scatter_cloud(l, i)=forward_scatter_cloud(l, i)              &
            *k_ext_tot_corr*omega_process_corr                                 &
            *asymmetry_process_corr**n_order_forward
        END IF
      END IF
      k_ext_tot_cloud(l, i)=k_ext_tot_cloud(l, i)*k_ext_tot_corr
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE opt_prop_inhom_corr_cairns
END MODULE opt_prop_inhom_corr_cairns_mod
