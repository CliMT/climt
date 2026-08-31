! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate two-stream coefficients in cloudy regions.
!
! Method:
!   The coefficients for each region are determined and averaged.
!
!- ---------------------------------------------------------------------
MODULE two_coeff_region_fast_lw_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TWO_COEFF_REGION_FAST_LW_MOD'
CONTAINS
SUBROUTINE two_coeff_region_fast_lw(ierr                                &
    , n_profile, n_layer, n_cloud_top                                   &
    , l_ir_source_quad, n_source_coeff                                  &
    , n_cloud_type, frac_cloud                                          &
    , n_region, i_region_cloud, frac_region                             &
    , tau_clr, tau                                                      &
    , isolir                                                            &
    , trans, reflect, source_coeff                                      &
    , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_source_coeff        &
    , nd_cloud_type, nd_region                                          &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_fatal, ip_infra_red, ip_region_clear
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength
  USE two_coeff_fast_lw_mod, ONLY: two_coeff_fast_lw

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer                                                          &
!       Maximum number of layers
    , nd_layer_clr                                                      &
!       Maximum number of completely clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_source_coeff                                                   &
!       Size allocated for source coefficients
    , nd_cloud_type                                                     &
!       Maximum number of types of cloud
    , nd_region
!       Maximum number of cloudy regions

! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_cloud_top                                                       &
!       Topmost cloudy layer
    , isolir                                                            &
!       Spectral region
    , n_cloud_type                                                      &
!       Number of types of clouds
    , n_source_coeff
!       Number of source coefficients

  INTEGER, INTENT(IN) ::                                                &
      n_region                                                          &
!       Number of cloudy regions
    , i_region_cloud(nd_cloud_type)
!       Regions in which types of clouds fall

  LOGICAL, INTENT(IN) ::                                                &
      l_ir_source_quad
!       Use a quadratic source in the infra-red

! optical properties of layer:
  REAL (RealK), INTENT(IN) ::                                           &
      frac_cloud(nd_profile, nd_layer, nd_cloud_type)                   &
!       Fractions of different types of clouds
    , frac_region(nd_profile, nd_layer, nd_region)                      &
!       Fractions of total cloud occupied by each region
    , tau_clr(nd_profile, nd_layer_clr)                                 &
!       Clear-sky optical depth
    , tau(nd_profile, id_ct: nd_layer, 0: nd_cloud_type)
!       Optical depth

! Coefficients in the two-stream equations:
  REAL (RealK), INTENT(OUT) ::                                          &
      trans(nd_profile, nd_layer, nd_region)                            &
!       Diffuse transmission coefficient
    , reflect(nd_profile, nd_layer, nd_region)                          &
!       Diffuse reflection coefficient
    , source_coeff(nd_profile, nd_layer                                 &
      , nd_source_coeff, nd_region)
!       Source coefficients in two-stream equations

! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , j                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , i_region
!       Loop variable over regions

! Coefficients in the two-stream equations:
  REAL (RealK) ::                                                       &
      trans_temp(nd_profile, nd_layer)                                  &
!       Temporary diffuse transmission coefficient
    , source_coeff_temp(nd_profile, nd_layer, nd_source_coeff)
!       Temporary source coefficients in two-stream equations

! Variables for gathering:
  INTEGER                                                               &
      n_list                                                            &
!       Number of points in list
    , l_list(nd_profile)                                                &
!       List of collected points
    , ll
  REAL (RealK) ::                                                       &
      tau_gathered(nd_profile, nd_layer)                                &
!       Gathered optical depth
    , tmp_inv(nd_profile)
!       Temporary work array

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'TWO_COEFF_REGION_FAST_LW'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! This routine should not be used outside the IR.
  IF (isolir /= ip_infra_red) THEN
    cmessage = '*** Erroneous use of non-scattering code.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)
  END IF

! Determine the optical properties of the clear-sky regions of
! the layers.

  CALL two_coeff_fast_lw(n_profile, 1, n_cloud_top-1                    &
    , l_ir_source_quad, tau_clr                                         &
    , trans(1, 1, ip_region_clear)                                      &
    , source_coeff(1, 1, 1, ip_region_clear)                            &
    , nd_profile, nd_layer, 1, nd_layer_clr, nd_source_coeff            &
    )
  CALL two_coeff_fast_lw(n_profile, n_cloud_top, n_layer                &
    , l_ir_source_quad, tau_clr                                         &
    , trans(1, 1, ip_region_clear)                                      &
    , source_coeff(1, 1, 1, ip_region_clear)                            &
    , nd_profile, nd_layer, id_ct, nd_layer, nd_source_coeff            &
    )
  DO i=1, n_layer
    DO l=1, n_profile
      reflect(l, i, ip_region_clear)=0.0e+00_RealK
    END DO
  END DO


! Now deal with clouds.

! Initialize the full arrays for cloudy regions.

  DO i_region=1, n_region
    IF (i_region /= ip_region_clear) THEN
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          trans(l, i, i_region)=0.0e+00_RealK
          reflect(l, i, i_region)=0.0e+00_RealK
        END DO
      END DO
      DO j=1, n_source_coeff
        DO i=n_cloud_top, n_layer
          DO l=1, n_profile
            source_coeff(l, i, j, i_region)=0.0e+00_RealK
          END DO
        END DO
      END DO

    END IF

  END DO



! Consider each type of cloud in turn, checking which region it
! contrubutes to and form weighted sums of cloud properties.

  DO k=1, n_cloud_type


!   Set the region in which clouds of this type are included.
    i_region=i_region_cloud(k)

    DO i=n_cloud_top, n_layer

!     Form a list of points where cloud of this type exists
!     on this row for gathering.
      n_list=0
      DO l=1, n_profile
        IF (frac_cloud(l, i, k) >  0.0e+00_RealK) THEN
          n_list=n_list+1
          l_list(n_list)=l
        END IF
      END DO


      IF (n_list >  0) THEN

!       Gather the optical properties. Though we consider only
!       one layer at a time the lower routines will operate on
!       arrays with vertical structure, so the gathered arrays
!       are two-dimensional.

        DO l=1, n_list
          tau_gathered(l, i)=tau(l_list(l), i, k)
        END DO

        CALL two_coeff_fast_lw(n_list, i, i                             &
          , l_ir_source_quad, tau_gathered                              &
          , trans_temp                                                  &
          , source_coeff_temp                                           &
          , nd_profile, nd_layer, id_ct, nd_layer, nd_source_coeff      &
          )

        DO ll=1, n_list
          l=l_list(ll)
          trans(l, i, i_region)=trans(l, i, i_region)                   &
            +frac_cloud(l, i, k)*trans_temp(ll, i)
        END DO
        DO j=1, n_source_coeff
          DO ll=1, n_list
            l=l_list(ll)
            source_coeff(l, i, j, i_region)                             &
              =source_coeff(l, i, j, i_region)                          &
              +frac_cloud(l, i, k)                                      &
              *source_coeff_temp(ll, i, j)
          END DO
        END DO

      END IF

    END DO
  END DO


! Finally, scale the weighted sums by the cloud fractions.
  DO i_region=1, n_region
    IF (i_region /= ip_region_clear) THEN
      DO i=n_cloud_top, n_layer

!       Gather points within this region.
        n_list=0
        DO l=1,n_profile
          IF (frac_region(l, i, i_region) >  0.0e+00_RealK) THEN
            n_list=n_list+1
            l_list(n_list)=l
          END IF
        END DO
        DO ll=1, n_list
          l=l_list(ll)
          tmp_inv(ll)=1.0e+00_RealK/frac_region(l, i, i_region)
          trans(l, i, i_region)=trans(l, i, i_region)                   &
            *tmp_inv(ll)
        END DO
        DO j=1, n_source_coeff
          DO ll=1, n_list
            l=l_list(ll)
            source_coeff(l, i, j, i_region)                             &
              =source_coeff(l, i, j, i_region)                          &
              *tmp_inv(ll)
          END DO
        END DO
      END DO
    END IF
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE two_coeff_region_fast_lw
END MODULE two_coeff_region_fast_lw_mod
