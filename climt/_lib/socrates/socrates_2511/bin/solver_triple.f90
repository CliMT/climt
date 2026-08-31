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
!- ---------------------------------------------------------------------
MODULE solver_triple_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLVER_TRIPLE_MOD'
CONTAINS
SUBROUTINE solver_triple(n_profile, n_layer, n_cloud_top                &
     , t, r, s_down, s_up                                               &
     , t_strat, r_strat, s_down_strat, s_up_strat                       &
     , t_conv, r_conv, s_down_conv, s_up_conv                           &
     , v11, v12, v13, v21, v22, v23, v31, v32, v33                      &
     , u11, u12, u13, u21, u22, u23, u31, u32, u33                      &
     , flux_inc_down                                                    &
     , source_ground_free, source_ground_strat                          &
     , source_ground_conv, albedo_surface_diff                          &
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
    , t_strat(nd_profile, nd_layer)                                     &
!       Stratfiform transmission
    , r_strat(nd_profile, nd_layer)                                     &
!       Stratfiform reflection
    , s_down_strat(nd_profile, nd_layer)                                &
!       Downward stratfiform source function
    , s_up_strat(nd_profile, nd_layer)                                  &
!       Upward stratfiform source function
    , t_conv(nd_profile, nd_layer)                                      &
!       Convective transmission
    , r_conv(nd_profile, nd_layer)                                      &
!       Convective reflection
    , s_down_conv(nd_profile, nd_layer)                                 &
!       Downward convective source function
    , s_up_conv(nd_profile, nd_layer)
!       Upward convective source function
  REAL (RealK), INTENT(IN) ::                                           &
      v11(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v12(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v13(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v21(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v22(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v23(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v31(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v32(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v33(nd_profile, id_ct-1: nd_layer)
!       Energy transfer coefficient for downward radiation
  REAL (RealK) ::                                                       &
      u11(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u12(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u13(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u21(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u22(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u23(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u31(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u32(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u33(nd_profile, id_ct-1: nd_layer)
!       Energy transfer coefficient for upward radiation
  REAL (RealK), INTENT(IN) ::                                           &
      flux_inc_down(nd_profile)                                         &
!       Incident flux
    , source_ground_free(nd_profile)                                    &
!       Source from ground (clear sky)
    , source_ground_strat(nd_profile)                                   &
!       Source from ground (cloudy region)
    , source_ground_conv(nd_profile)                                    &
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
    , alpha12(nd_profile)                                               &
    , alpha13(nd_profile)                                               &
    , alpha21(nd_profile)                                               &
    , alpha22(nd_profile)                                               &
    , alpha23(nd_profile)                                               &
    , alpha31(nd_profile)                                               &
    , alpha32(nd_profile)                                               &
    , alpha33(nd_profile)                                               &
    , g1(nd_profile)                                                    &
    , g2(nd_profile)                                                    &
    , g3(nd_profile)
! Terms for downward propagation:
  REAL (RealK) ::                                                       &
      gamma11(nd_profile, nd_layer)                                     &
    , gamma12(nd_profile, nd_layer)                                     &
    , gamma13(nd_profile, nd_layer)                                     &
    , gamma21(nd_profile, nd_layer)                                     &
    , gamma22(nd_profile, nd_layer)                                     &
    , gamma23(nd_profile, nd_layer)                                     &
    , gamma31(nd_profile, nd_layer)                                     &
    , gamma32(nd_profile, nd_layer)                                     &
    , gamma33(nd_profile, nd_layer)                                     &
    , beta11_inv(nd_profile, nd_layer)                                  &
    , beta21(nd_profile, nd_layer)                                      &
    , beta22_inv(nd_profile, nd_layer)                                  &
    , beta31(nd_profile, nd_layer)                                      &
    , beta32(nd_profile, nd_layer)                                      &
    , beta33_inv(nd_profile, nd_layer)                                  &
    , h1(nd_profile, nd_layer)                                          &
    , h2(nd_profile, nd_layer)                                          &
    , h3(nd_profile, nd_layer)

! Auxilairy numerical variables required only in the current layer:
  REAL (RealK) ::                                                       &
      theta11                                                           &
    , theta12                                                           &
    , theta13                                                           &
    , theta21                                                           &
    , theta22                                                           &
    , theta23                                                           &
    , theta31                                                           &
    , theta32                                                           &
    , theta33                                                           &
    , lambda3                                                           &
    , lambda2                                                           &
    , lambda1                                                           &
    , lambda                                                            &
    , dum1, dum2, dum3
!           dummy variables

! Temporary fluxes
  REAL (RealK) ::                                                       &
      flux_down_1(nd_profile)                                           &
!       Downward fluxes outside clouds just below i'th level
    , flux_down_2(nd_profile)                                           &
!       Downward fluxes inside clouds just below i'th level
    , flux_down_3(nd_profile)                                           &
!       Downward fluxes inside clouds just below i'th level
    , flux_up_1(nd_profile)                                             &
!       Upward fluxes outside clouds just above i'th level
    , flux_up_2(nd_profile)                                             &
!       Upward fluxes inside clouds just above i'th level
    , flux_up_3(nd_profile)
!       Upward fluxes inside clouds just above i'th level

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SOLVER_TRIPLE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! This routine is specific to cases of three regions and it is
! assumed that 1 represents clear skies, 2 represents startiform
! clouds and 3 represents convective cloud.

! Upward elimination through the cloudy layers.
  DO i=n_layer+1, 1, -1
    IF ( i < n_layer+1 .AND. i >= n_cloud_top ) THEN
      DO l=1, n_profile

        theta11=alpha11(l)*v11(l, i)+alpha12(l)*v21(l, i)               &
          +alpha13(l)*v31(l, i)
        theta12=alpha11(l)*v12(l, i)+alpha12(l)*v22(l, i)               &
          +alpha13(l)*v32(l, i)
        theta13=alpha11(l)*v13(l, i)+alpha12(l)*v23(l, i)               &
          +alpha13(l)*v33(l, i)
        theta21=alpha21(l)*v11(l, i)+alpha22(l)*v21(l, i)               &
          +alpha23(l)*v31(l, i)
        theta22=alpha21(l)*v12(l, i)+alpha22(l)*v22(l, i)               &
          +alpha23(l)*v32(l, i)
        theta23=alpha21(l)*v13(l, i)+alpha22(l)*v23(l, i)               &
          +alpha23(l)*v33(l, i)
        theta31=alpha31(l)*v11(l, i)+alpha32(l)*v21(l, i)               &
          +alpha33(l)*v31(l, i)
        theta32=alpha31(l)*v12(l, i)+alpha32(l)*v22(l, i)               &
          +alpha33(l)*v32(l, i)
        theta33=alpha31(l)*v13(l, i)+alpha32(l)*v23(l, i)               &
          +alpha33(l)*v33(l, i)
        beta31(l, i)=-theta31*r(l, i)
        beta32(l, i)=-theta32*r_strat(l, i)
        beta33_inv(l, i)=1.0e+00_RealK                                  &
          /(1.0e+00_RealK-theta33*r_conv(l, i))
        gamma31(l, i)=theta31*t(l, i)
        gamma32(l, i)=theta32*t_strat(l, i)
        gamma33(l, i)=theta33*t_conv(l, i)
        h3(l, i)=g3(l)+theta31*s_down(l, i)                             &
          +theta32*s_down_strat(l, i)                                   &
          +theta33*s_down_conv(l, i)

        lambda3=theta23*r_conv(l, i)*beta33_inv(l, i)
        beta22_inv(l, i)=1.0e+00_RealK                                  &
          /(1.0e+00_RealK-theta22*r_strat(l, i)+lambda3*beta32(l, i))
        beta21(l, i)=-theta21*r(l, i)+lambda3*beta31(l, i)
        gamma21(l, i)=theta21*t(l, i)+lambda3*gamma31(l, i)
        gamma22(l, i)=theta22*t_strat(l, i)+lambda3*gamma32(l, i)
        gamma23(l, i)=theta23*t_conv(l, i)+lambda3*gamma33(l, i)
        h2(l, i)=g2(l)+theta21*s_down(l, i)                             &
          +theta22*s_down_strat(l, i)+theta23*s_down_conv(l, i)         &
          +lambda3*h3(l, i)

        lambda3=theta13*r_conv(l, i)*beta33_inv(l, i)
        lambda2=(theta12*r_strat(l, i)-lambda3*beta32(l, i))            &
          *beta22_inv(l, i)
        beta11_inv(l, i)=1.0e+00_RealK                                  &
          /(1.0e+00_RealK-theta11*r(l, i)+lambda3*beta31(l, i)          &
          +lambda2*beta21(l, i))
        gamma11(l, i)=theta11*t(l, i)+lambda3*gamma31(l, i)             &
          +lambda2*gamma21(l, i)
        gamma12(l, i)=theta12*t_strat(l, i)+lambda3*gamma32(l, i)       &
          +lambda2*gamma22(l, i)
        gamma13(l, i)=theta13*t_conv(l, i)+lambda3*gamma33(l, i)        &
          +lambda2*gamma23(l, i)
        h1(l, i)=g1(l)+theta11*s_down(l, i)                             &
          +theta12*s_down_strat(l, i)+theta13*s_down_conv(l, i)         &
          +lambda3*h3(l, i)+lambda2*h2(l, i)

        lambda3=u33(l, i-1)*t_conv(l, i)*beta33_inv(l, i)
        lambda2=(u32(l, i-1)*t_strat(l, i)+lambda3*beta32(l, i))        &
          *beta22_inv(l, i)
        lambda1=(u31(l, i-1)*t(l, i)+lambda3*beta31(l, i)               &
          +lambda2*beta21(l, i))*beta11_inv(l, i)
        alpha31(l)=u31(l, i-1)*r(l, i)+lambda3*gamma31(l, i)            &
          +lambda2*gamma21(l, i)+lambda1*gamma11(l, i)
        alpha32(l)=u32(l, i-1)*r_strat(l, i)                            &
          +lambda3*gamma32(l, i)+lambda2*gamma22(l, i)                  &
          +lambda1*gamma12(l, i)
        alpha33(l)=u33(l, i-1)*r_conv(l, i)                             &
          +lambda3*gamma33(l, i)+lambda2*gamma23(l, i)                  &
          +lambda1*gamma13(l, i)
        g3(l)=u31(l, i-1)*s_up(l, i)+u32(l, i-1)*s_up_strat(l, i)       &
          +u33(l, i-1)*s_up_conv(l, i)                                  &
          +lambda3*h3(l, i)+lambda2*h2(l, i)+lambda1*h1(l, i)

        lambda3=u23(l, i-1)*t_conv(l, i)*beta33_inv(l, i)
        lambda2=(u22(l, i-1)*t_strat(l, i)+lambda3*beta32(l, i))        &
          *beta22_inv(l, i)
        lambda1=(u21(l, i-1)*t(l, i)+lambda3*beta31(l, i)               &
          +lambda2*beta21(l, i))*beta11_inv(l, i)
        alpha21(l)=u21(l, i-1)*r(l, i)+lambda3*gamma31(l, i)            &
          +lambda2*gamma21(l, i)+lambda1*gamma11(l, i)
        alpha22(l)=u22(l, i-1)*r_strat(l, i)                            &
          +lambda3*gamma32(l, i)+lambda2*gamma22(l, i)                  &
          +lambda1*gamma12(l, i)
        alpha23(l)=u23(l, i-1)*r_conv(l, i)                             &
          +lambda3*gamma33(l, i)+lambda2*gamma23(l, i)                  &
          +lambda1*gamma13(l, i)
        g2(l)=u21(l, i-1)*s_up(l, i)+u22(l, i-1)*s_up_strat(l, i)       &
          +u23(l, i-1)*s_up_conv(l, i)                                  &
          +lambda3*h3(l, i)+lambda2*h2(l, i)+lambda1*h1(l, i)

        lambda3=u13(l, i-1)*t_conv(l, i)*beta33_inv(l, i)
        lambda2=(u12(l, i-1)*t_strat(l, i)+lambda3*beta32(l, i))        &
          *beta22_inv(l, i)
        lambda1=(u11(l, i-1)*t(l, i)+lambda3*beta31(l, i)               &
          +lambda2*beta21(l, i))*beta11_inv(l, i)
        alpha11(l)=u11(l, i-1)*r(l, i)+lambda3*gamma31(l, i)            &
          +lambda2*gamma21(l, i)+lambda1*gamma11(l, i)
        alpha12(l)=u12(l, i-1)*r_strat(l, i)                            &
          +lambda3*gamma32(l, i)+lambda2*gamma22(l, i)                  &
          +lambda1*gamma12(l, i)
        alpha13(l)=u13(l, i-1)*r_conv(l, i)                             &
          +lambda3*gamma33(l, i)+lambda2*gamma23(l, i)                  &
          +lambda1*gamma13(l, i)
        g1(l)=u11(l, i-1)*s_up(l, i)+u12(l, i-1)*s_up_strat(l, i)       &
          +u13(l, i-1)*s_up_conv(l, i)                                  &
          +lambda3*h3(l, i)+lambda2*h2(l, i)+lambda1*h1(l, i)

      END DO

    ELSE IF ( i < n_cloud_top-1 ) THEN

!     Continue through the cloud free region: if there is no such
!     region the DO-loop will not be executed.
      DO l=1, n_profile

        beta11_inv(l, i)=1.0e+00_RealK                                  &
          /(1.0e+00_RealK-alpha11(l)*r(l, i))
        gamma11(l, i)=alpha11(l)*t(l, i)
        h1(l, i)=g1(l)+alpha11(l)*s_down(l, i)

        lambda1=t(l, i)*beta11_inv(l, i)
        alpha11(l)=r(l, i)+lambda1*gamma11(l, i)
        g1(l)=s_up(l, i)+lambda1*h1(l, i)

      END DO

    ELSE IF ( i == n_cloud_top-1 ) THEN

!     The layer above the cloud, if the cloud does not reach to the top
!     of the column: only one set of alphas is now needed.
      DO l=1, n_profile

        IF (n_cloud_top <  n_layer) THEN
!         If there is no cloud in the column the V's will not be
!         assigned so an IF-test is required.
          theta11=alpha11(l)*v11(l, i)+alpha12(l)*v21(l, i)             &
            +alpha13(l)*v31(l, i)
        ELSE
          theta11=alpha11(l)
        END IF

        beta11_inv(l,i)=1.0e+00_RealK/(1.0e+00_RealK-theta11*r(l,i))
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
        alpha12(l)=0.0e+00_RealK
        alpha13(l)=0.0e+00_RealK
        alpha21(l)=0.0e+00_RealK
        alpha22(l)=albedo_surface_diff(l)
        alpha23(l)=0.0e+00_RealK
        alpha31(l)=0.0e+00_RealK
        alpha32(l)=0.0e+00_RealK
        alpha33(l)=albedo_surface_diff(l)
        g1(l)=source_ground_free(l)
        g2(l)=source_ground_strat(l)
        g3(l)=source_ground_conv(l)
      END DO

    END IF
  END DO


! Initialize for downward back-substitution.
  DO l=1, n_profile
    flux_total(l, 2)=flux_inc_down(l)
  END DO
  IF (n_cloud_top >  1) THEN
    DO l=1, n_profile
      flux_total(l, 1)=alpha11(l)*flux_total(l, 2)+g1(l)
    END DO
  ELSE
    DO l=1, n_profile
      flux_total(l, 1)=g1(l)+flux_inc_down(l)                           &
        *(v11(l, 0)*alpha11(l)+v21(l, 0)*alpha12(l)                     &
        + v31(l, 0)*alpha13(l))
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

! Pass into the top cloudy layer. use flux_down_[1,2,3] to hold,
! provisionally, the downward fluxes just below the top of the
! layer, then calculate the upward fluxes at the bottom and
! finally the downward fluxes at the bottom of the layer.

! The main loop of back-substitution. the provisional use of the
! downward fluxes is as above.
  DO i=n_cloud_top, n_layer

    IF ( i > n_cloud_top ) THEN
      DO l=1, n_profile
        dum1 = flux_down_1(l)
        dum2 = flux_down_2(l)
        dum3 = flux_down_3(l)
        flux_down_1(l)=v11(l, i-1)*dum1 + v12(l, i-1)*dum2              &
                      +v13(l, i-1)*dum3                   
        flux_down_2(l)=v21(l, i-1)*dum1 + v22(l, i-1)*dum2              &
                      +v23(l, i-1)*dum3                   
        flux_down_3(l)=v31(l, i-1)*dum1 + v32(l, i-1)*dum2              &
                      +v33(l, i-1)*dum3
        flux_up_1(l)=(gamma11(l, i)*flux_down_1(l)                      &
          +gamma12(l, i)*flux_down_2(l)                                 &
          +gamma13(l, i)*flux_down_3(l) + h1(l, i))*beta11_inv(l, i)
        flux_up_2(l)=(gamma21(l, i)*flux_down_1(l)                      &
          +gamma22(l, i)*flux_down_2(l)                                 &
          +gamma23(l, i)*flux_down_3(l)+h2(l, i)                        &
          -beta21(l, i)*flux_up_1(l))*beta22_inv(l, i)
        flux_up_3(l)=(gamma31(l, i)*flux_down_1(l)                      &
          +gamma32(l, i)*flux_down_2(l)                                 &
          +gamma33(l, i)*flux_down_3(l)+h3(l, i)                        &
          -beta31(l, i)*flux_up_1(l)                                    &
          -beta32(l, i)*flux_up_2(l))                                   &
          *beta33_inv(l, i)
        flux_down_1(l)=t(l, i)*flux_down_1(l)                           &
          +r(l, i)*flux_up_1(l)+s_down(l, i)
        flux_down_2(l)=t_strat(l, i)*flux_down_2(l)                     &
          +r_strat(l, i)*flux_up_2(l)+s_down_strat(l, i)
        flux_down_3(l)=t_conv(l, i)*flux_down_3(l)                      &
          +r_conv(l, i)*flux_up_3(l)+s_down_conv(l, i)
      END DO

    ELSE

!     If there are no clouds n_cloud_top may be out-of-bounds for
!     these arrays so an if test is required.
      DO l=1, n_profile
        flux_down_1(l)=v11(l, i-1)*flux_total(l, 2*i)
        flux_down_2(l)=v21(l, i-1)*flux_total(l, 2*i)
        flux_down_3(l)=v31(l, i-1)*flux_total(l, 2*i)
        flux_up_1(l)=(gamma11(l, i)*flux_down_1(l)                      &
          +gamma12(l, i)*flux_down_2(l)                                 &
          +gamma13(l, i)*flux_down_3(l)                                 &
          +h1(l, i))*beta11_inv(l, i)
        flux_up_2(l)=(gamma21(l, i)*flux_down_1(l)                      &
          +gamma22(l, i)*flux_down_2(l)                                 &
          +gamma23(l, i)*flux_down_3(l)+h2(l, i)                        &
          -beta21(l, i)*flux_up_1(l))*beta22_inv(l, i)
        flux_up_3(l)=(gamma31(l, i)*flux_down_1(l)                      &
          +gamma32(l, i)*flux_down_2(l)                                 &
          +gamma33(l, i)*flux_down_3(l)+h3(l, i)                        &
          -beta31(l, i)*flux_up_1(l)-beta32(l, i)*flux_up_2(l))         &
          *beta33_inv(l, i)
        flux_down_1(l)=t(l, i)*flux_down_1(l)                           &
          +r(l, i)*flux_up_1(l)+s_down(l, i)
        flux_down_2(l)=t_strat(l, i)*flux_down_2(l)                     &
          +r_strat(l, i)*flux_up_2(l)+s_down_strat(l, i)
        flux_down_3(l)=t_conv(l, i)*flux_down_3(l)                      &
          +r_conv(l, i)*flux_up_3(l)+s_down_conv(l, i)
      END DO
    END IF

!   Calculate the overall flux.
    DO l=1, n_profile
      flux_total(l, 2*i+1)=flux_up_1(l)+flux_up_2(l)                    &
        +flux_up_3(l)
      flux_total(l, 2*i+2)=flux_down_1(l)+flux_down_2(l)                &
        +flux_down_3(l)
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE solver_triple
END MODULE solver_triple_mod
