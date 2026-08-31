! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to solve for mixed fluxes scattering without a matrix.
!
! Method:
!   Gaussian elimination in an upward direction is employed to
!   determine effective albedos for lower levels of the atmosphere.
!   This allows a downward pass of back-substitution to be carried
!   out to determine the upward and downward fluxes.
!
!   This version has been modified by Robin Hogan to allow
!   shadowing, as documented in Shonk & Hogan, 2007, J. Climate.
!
!- ---------------------------------------------------------------------
MODULE solver_mix_direct_hogan_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLVER_MIX_DIRECT_HOGAN_MOD'
CONTAINS
SUBROUTINE solver_mix_direct_hogan(n_profile, n_layer, n_cloud_top      &
     , t, r, s_down, s_up                                               &
     , t_cloud, r_cloud, s_down_cloud, s_up_cloud                       &
     , v11, v21, v12, v22                                               &
     , u11, u12, u21, u22                                               &
     , flux_inc_down                                                    &
     , source_ground_free, source_ground_cloud, albedo_surface_diff     &
     , flux_total                                                       &
     , nd_profile, nd_layer, id_ct                                      &
     )


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , id_ct
!       Topmost declared cloudy layer


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_cloud_top
!       Topmost cloudy layer
  REAL (RealK), INTENT(IN) ::                                           &
      t(nd_profile, nd_layer)                                           &
!       Clear-sky transmission
    , r(nd_profile, nd_layer)                                           &
!       Clear-sky reflection
    , s_down(nd_profile, nd_layer)                                      &
!       Clear-sky downward source function
    , s_up(nd_profile, nd_layer)                                        &
!       Clear-sky upward source function
    , t_cloud(nd_profile, nd_layer)                                     &
!       Cloudy transmission
    , r_cloud(nd_profile, nd_layer)                                     &
!       Cloudy reflection
    , s_down_cloud(nd_profile, nd_layer)                                &
!       Downward cloudy source function
    , s_up_cloud(nd_profile, nd_layer)
!       Upward cloudy source function
  REAL (RealK), INTENT(IN) ::                                           &
      v11(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , v21(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , v12(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , v22(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , u11(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , u12(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , u21(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , u22(nd_profile, id_ct-1: nd_layer)
!       Energy transfer coefficient
  REAL (RealK), INTENT(IN) ::                                           &
      flux_inc_down(nd_profile)                                         &
!       Incident flux
    , source_ground_free(nd_profile)                                    &
!       Source from ground (clear sky)
    , source_ground_cloud(nd_profile)                                   &
!       Source from ground (cloudy region)
    , albedo_surface_diff(nd_profile)
!       Diffuse albedo
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_total(nd_profile, 2*nd_layer+2)
!       Total flux

! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable

! Effective coupling albedos and source functions:
  REAL (RealK) ::                                                       &
      alpha11(nd_profile)                                               &
    , alpha22(nd_profile)                                               &
    , g1(nd_profile)                                                    &
    , g2(nd_profile)
! Terms for downward propagation:
  REAL (RealK) ::                                                       &
      gamma11(nd_profile, nd_layer)                                     &
    , gamma22(nd_profile, nd_layer)                                     &
    , beta11_inv(nd_profile, nd_layer)                                  &
    , beta22_inv(nd_profile, nd_layer)                                  &
    , h1(nd_profile, nd_layer)                                          &
    , h2(nd_profile, nd_layer)

! Auxilairy numerical variables required only in the current layer:
  REAL (RealK) ::                                                       &
      theta11                                                           &
    , theta22                                                           &
    , lambda1                                                           &
    , lambda2                                                           &
    , lambda                                                            &
    , dum1                                                              &
    , dum2

! Temporary fluxes
  REAL (RealK) ::                                                       &
      flux_down_1(nd_profile)                                           &
!       Downward fluxes outside clouds just below I'th level
    , flux_down_2(nd_profile)                                           &
!       Downward fluxes inside clouds just below I'th level
    , flux_up_1(nd_profile)                                             &
!       Upward fluxes outside clouds just above I'th level
    , flux_up_2(nd_profile)
!       Upward fluxes inside clouds just above I'th level

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SOLVER_MIX_DIRECT_HOGAN'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Upward elimination through the cloudy layers.
  DO i=n_layer+1, 1, -1

    IF ( i < n_layer+1 .AND. i >= n_cloud_top ) THEN

      DO l=1, n_profile

        theta11 = alpha11(l)*v11(l, i)+alpha22(l)*v21(l, i)
        theta22 = alpha11(l)*v12(l, i)+alpha22(l)*v22(l, i)

        beta11_inv(l, i)=1.0_RealK/(1.0_RealK-theta11*r(l, i))
        gamma11(l, i)=theta11*t(l, i)
        h1(l, i)=g1(l)+theta11*s_down(l, i)

        beta22_inv(l, i)=1.0_RealK/(1.0_RealK-theta22*r_cloud(l, i))
        gamma22(l, i)=theta22*t_cloud(l, i)
        h2(l, i)=g2(l)+theta22*s_down_cloud(l, i)

        lambda1 = s_up(l, i)+h1(l, i)*t(l, i)*beta11_inv(l, i)
        lambda2 = s_up_cloud(l, i)+h2(l, i)*t_cloud(l, i)               &
             * beta22_inv(l, i)

        alpha11(l)=r(l, i)                                              &
             + theta11*t(l, i)*t(l, i)*beta11_inv(l, i)
        g1(l)=u11(l, i-1)*lambda1 + u12(l, i-1)*lambda2

        alpha22(l)=r_cloud(l, i)                                        &
             + theta22*t_cloud(l, i)*t_cloud(l, i)*beta22_inv(l, i)
        g2(l)=u21(l, i-1)*lambda1 + u22(l, i-1)*lambda2

      END DO

    ELSE IF ( i < n_cloud_top-1 ) THEN

      DO l=1, n_profile

        beta11_inv(l, i)=1.0e+00_RealK                                  &
          /(1.0e+00_RealK-alpha11(l)*r(l, i))
        gamma11(l, i)=alpha11(l)*t(l, i)
        h1(l, i)=g1(l)+alpha11(l)*s_down(l, i)

        lambda=t(l, i)*beta11_inv(l, i)
        alpha11(l)=r(l, i)+lambda*gamma11(l, i)
        g1(l)=s_up(l, i)+lambda*h1(l, i)

      END DO

    ELSE IF ( i == n_cloud_top-1 ) THEN

!     The layer above the cloud: only one set of alphas is now needed.
!     This will not be presented if there is cloud in the top layer.
      DO l=1, n_profile

        IF (n_cloud_top < n_layer) THEN
!         If there is no cloud in the column the V's will not be
!         assigned so an if test is required.
          theta11=alpha11(l)*v11(l, i)+alpha22(l)*v21(l, i)
        ELSE
          theta11=alpha11(l)
        END IF

        beta11_inv(l, i)=1.0e+00_RealK/(1.0e+00_RealK-theta11*r(l, i))
        gamma11(l, i)=theta11*t(l, i)
        h1(l, i)=g1(l)+theta11*s_down(l, i)

        lambda=t(l, i)*beta11_inv(l, i)
        alpha11(l)=r(l, i)+lambda*gamma11(l, i)
        g1(l)=s_up(l, i)+lambda*h1(l, i)

      END DO

    ELSE IF ( I == n_layer+1 ) THEN

!     Initialize at the bottom of the column for upward elimination.
      DO l=1, n_profile
        alpha11(l)=albedo_surface_diff(l)
        alpha22(l)=albedo_surface_diff(l)
        g1(l)=source_ground_free(l)
        g2(l)=source_ground_cloud(l)
      END DO

    END IF

  END DO


! Initialize for downward back-substitution.
  DO l=1, n_profile
    flux_total(l, 2)=flux_inc_down(l)
  END DO
  IF (n_cloud_top > 1) THEN
    DO l=1, n_profile
      flux_total(l, 1)=alpha11(l)*flux_total(l, 2)+g1(l)
    END DO
  ELSE
    DO l=1, n_profile
      flux_total(l, 1)=g1(l)+flux_inc_down(l)                           &
        *(v11(l, 0)*alpha11(l)+v21(l, 0)*alpha22(l))
    END DO
  END IF

! Sweep downward through the clear-sky region, finding the downward
! flux at the top of the layer and the upward flux at the bottom.
  DO i=1, n_cloud_top-1
    DO l=1, n_profile
      flux_total(l, 2*i+1)=(gamma11(l, i)*flux_total(l, 2*i)            &
        +h1(l, i))*beta11_inv(l, i)
      flux_total(l, 2*i+2)=t(l, i)*flux_total(l, 2*i)                   &
        +r(l, i)*flux_total(l, 2*i+1)+s_down(l, i)
    END DO
  END DO

! Pass into the top cloudy layer. Use FLUX_DOWN_[1,2] to hold,
! provisionally, the downward fluxes just below the top of the
! layer, then calculate the upward fluxes at the bottom and
! finally the downward fluxes at the bottom of the layer.

! The main loop of back-substitution. The provisional use of the
! downward fluxes is as above.
  DO i=n_cloud_top, n_layer

    IF ( i > n_cloud_top ) THEN
      DO l=1, n_profile
        dum1 = flux_down_1(l)
        dum2 = flux_down_2(l)
        flux_down_1(l)=v11(l, i-1)*dum1+v12(l, i-1)*dum2
        flux_down_2(l)=v21(l, i-1)*dum1+v22(l, i-1)*dum2
        flux_up_1(l)=(gamma11(l, i)*flux_down_1(l)                      &
           +h1(l, i))*beta11_inv(l, i)
        flux_up_2(l)=(gamma22(l, i)*flux_down_2(l)                      &
           +h2(l, i))*beta22_inv(l, i)
        flux_down_1(l)=t(l, i)*flux_down_1(l)                           &
           +r(l, i)*flux_up_1(l)+s_down(l, i)
        flux_down_2(l)=t_cloud(l, i)*flux_down_2(l)                     &
           +r_cloud(l, i)*flux_up_2(l)+s_down_cloud(l, i)
      END DO
    ELSE
      DO l=1, n_profile
        flux_down_1(l)=v11(l, i-1)*flux_total(l, 2*i)
        flux_down_2(l)=v21(l, i-1)*flux_total(l, 2*i)
        flux_up_1(l)=(gamma11(l, i)*flux_down_1(l)                      &
           +h1(l, i))*beta11_inv(l, i)
        flux_up_2(l)=(gamma22(l, i)*flux_down_2(l)                      &
             +h2(l, i))*beta22_inv(l, i)
        flux_down_1(l)=t(l, i)*flux_down_1(l)                           &
           +r(l, i)*flux_up_1(l)+s_down(l, i)
        flux_down_2(l)=t_cloud(l, i)*flux_down_2(l)                     &
           +r_cloud(l, i)*flux_up_2(l)+s_down_cloud(l, i)
      END DO
    END IF


!   Calculate the overall flux.
    DO l=1, n_profile
      flux_total(l, 2*i+1)=flux_up_1(l)+flux_up_2(l)
      flux_total(l, 2*i+2)=flux_down_1(l)+flux_down_2(l)
    END DO

  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE solver_mix_direct_hogan
END MODULE solver_mix_direct_hogan_mod
