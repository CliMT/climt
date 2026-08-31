! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to find energy transfer coefficients for coupled overlap.
!
! Method:
!   Energy transfer coefficients for upward and downward radiation
!   at the edges of the layers are calculated assuming maximal
!   overlap of regions of the same nature and random overlap of
!   regions of a different nature.
!
!   Storage and Indexing: Now that solvers for the net flux are no
!   longer supported, the overlap coefficients can be stored more
!   easily. The coefficient referring to downward transfer from the
!   kth to the jth region is stored with a third index of
!   K+N_REGION*(J-1): the coeffieint for upward transfer is stored
!   with an index of N_REGION*(N_REGION+J-1)+K, so that in both
!   cases the originating region changes most frequently with the
!   index.
!
!- ---------------------------------------------------------------------
MODULE overlap_coupled_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'OVERLAP_COUPLED_MOD'
CONTAINS
SUBROUTINE overlap_coupled(n_profile, n_layer, n_cloud_top              &
    , w_cloud, w_free, n_region, type_region, frac_region, p            &
    , i_cloud                                                           &
    , cloud_overlap                                                     &
    , nd_profile, nd_layer, nd_overlap_coeff, nd_region                 &
    , id_ct, dp_corr_strat, dp_corr_conv, tot_cloud_cover               &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_cloud_mix_max, ip_cloud_mix_random,             &
                     ip_cloud_part_corr, ip_cloud_part_corr_cnv,        &
                     ip_cloud_triple, ip_region_clear, ip_region_conv,  &
                     ip_region_strat
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer                                                          &
!       Maximum number of layers
    , nd_overlap_coeff                                                  &
!       Maximum number of overlap coefficients
    , nd_region                                                         &
!       Maximum number of regions
    , id_ct
!       Topmost declared cloudy layer


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_cloud_top                                                       &
!       Topmost cloudy layer
    , n_region                                                          &
!       Number of cloudy regions
    , type_region(nd_region)                                            &
!       Array holding the type of each region
    , i_cloud
!       Cloud scheme selected
  REAL (RealK), INTENT(IN) ::                                           &
      w_cloud(nd_profile, id_ct: nd_layer)                              &
!       Cloud amounts
    , frac_region(nd_profile, id_ct: nd_layer, nd_region)               &
!       Fractions of total cloud amount occupied by
!       different regions
    , p(nd_profile, nd_layer)                                           &
!       Pressures at the middles of layers
    , dp_corr_strat                                                     &
!       Decorrelation pressure scale for large scale cloud
    , dp_corr_conv
!       Decorrelation pressure scale for convective cloud

  REAL (RealK), INTENT(OUT) ::                                          &
      w_free(nd_profile, id_ct: nd_layer)                               &
!       Cloud-free amounts
    , cloud_overlap(nd_profile, id_ct-1: nd_layer                       &
        , nd_overlap_coeff)
!       Coefficients for transfer of energy at interface

  REAL (RealK), INTENT(OUT) ::                                          &
       tot_cloud_cover(nd_profile)
!         Total cloud cover

! Local arguments.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , j                                                                 &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , k
!       Loop variable


! Fixed local values:
  REAL (RealK) ::                                                       &
      dp_corr                                                           &
!       Pressure scale over which correlation between cloudy
!       layers is lost
    , corr_factor(nd_profile, nd_region)
!       Correlation factors for each region across the boundary
!       between layers: this represents the fraction of the
!       potentially maximally overlapped region that is actually
!       maximally overlapped.
  REAL (RealK) ::                                                       &
      area_lower(nd_profile, nd_region)                                 &
!       Areas of regions in lower layer
    , area_upper(nd_profile, nd_region)                                 &
!       Areas of regions in lower layer
    , area_overlap(nd_profile, nd_region, nd_region)                    &
!       Areas of overlap between the different regions:
!       the first index refers to the upper layer
    , area_random_upper(nd_profile, nd_region)                          &
!       Areas of each region in the upper layer
!       to be overlapped randomly
    , area_random_lower(nd_profile, nd_region)                          &
!       Areas of each region in the lower layer
!       to be overlapped randomly
    , area_random_tot(nd_profile)                                       &
!       Total randomly overlapped area
    , tol_cloud
!       Tolerance used to detect cloud amounts of 0

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='OVERLAP_COUPLED'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  tol_cloud=1.0e+02_RealK*EPSILON(tol_cloud)

! Set the free fractions in each layer.
  DO i=n_cloud_top, n_layer
    DO l=1, n_profile
      w_free(l, i)=1.0e+00_RealK-w_cloud(l, i)
    END DO
  END DO

! Use the total cloud cover temporarily to hold the clear-sky
! fraction and convert back to cloud cover later.
  DO l=1, n_profile
    IF (w_free(l, n_cloud_top) > tol_cloud) THEN
      tot_cloud_cover(l)=w_free(l, n_cloud_top)
    ELSE
      tot_cloud_cover(l) = 0.0e+00_RealK
    END IF
  END DO


! We consider each boundary in turn, comparing the fractions
! of each region in the layers above and below the boundary.

! Initialize for the layer above the clouds: here the clear
! region will cover the grid-box.
  DO k=1, n_region
    IF (type_region(k) == ip_region_clear) THEN
      DO l=1, n_profile
        area_upper(l, k)=1.0e+00_RealK
      END DO
    ELSE
      DO l=1, n_profile
        area_upper(l, k)=0.0e+00_RealK
      END DO
    END IF
  END DO

  DO i=n_cloud_top-1, n_layer

!   Set the correlations between like regions at each interface.

    IF ( (i_cloud == ip_cloud_triple).OR.                               &
         (i_cloud == ip_cloud_mix_max) ) THEN

      DO k=1, n_region
        DO l=1, n_profile
          corr_factor(l, k)=1.0e+00_RealK
        END DO
      END DO

    ELSE IF (i_cloud == ip_cloud_mix_random) THEN

      DO k=1, n_region
        DO l=1, n_profile
          corr_factor(l, k)=0.0e+00_RealK
        END DO
      END DO

    ELSE IF ( (i_cloud == ip_cloud_part_corr).OR.                       &
              (i_cloud == ip_cloud_part_corr_cnv) ) THEN

      DO k=1, n_region

!       Experimental version: set the pressure scales over
!       which decorrelation occurs.
        IF (type_region(k) == ip_region_clear) THEN
          dp_corr=1.0e+00_RealK
        ELSE IF (type_region(k) == ip_region_strat) THEN
          dp_corr=dp_corr_strat
        ELSE IF (type_region(k) == ip_region_conv) THEN
          dp_corr=dp_corr_conv
        END IF

        IF ( (i <  n_layer).AND.(i >  1) ) THEN
          DO l=1, n_profile
            corr_factor(l, k)=EXP((p(l, i)-p(l, i+1))/dp_corr)
          END DO
        ELSE
!         At the surface and the top of the atmosphere
!         the correlation factor is irrelevant.
          DO l=1, n_profile
            corr_factor(l, k)=1.0e+00_RealK
          END DO
        END IF

      END DO

    END IF

!   Set areas of the regions in the lower layer.
    DO k=1, n_region
      IF (i <  n_layer) THEN
        IF (type_region(k) == ip_region_clear) THEN
          DO l=1, n_profile
            area_lower(l, k)=w_free(l, i+1)
          END DO
        ELSE
          DO l=1, n_profile
            area_lower(l, k)=w_cloud(l, i+1)                            &
              *frac_region(l, i+1, k)
          END DO
        END IF
      ELSE
!       At the very bottom of the column we imagine a notional
!       clear layer below the ground surface.
        IF (type_region(k) == ip_region_clear) THEN
          DO l=1, n_profile
            area_lower(l, k)=1.0e+00_RealK
          END DO
        ELSE
          DO l=1, n_profile
            area_lower(l, k)=0.0e+00_RealK
          END DO
        END IF
      END IF

!     Begin by setting the maximally overlapped parts of the
!     atmospheric column. The area of common overlap betwen
!     like regions may be incremented by randomly overlapped
!     fractions later.

      DO l=1, n_profile
        area_overlap(l, k, k)=corr_factor(l, k)                         &
          *MIN(area_lower(l, k), area_upper(l, k))
      END DO

    END DO

!   Calculate total column cloud cover for use as a diagnostic.
!   (Done at this stage to make use of the maximally overlapped
!   fractions before they are incremented.)
!   We calculate this quantity by imagining a totally transparent
!   atmosphere containing totally opaque clouds and finding the
!   transmission.

    IF ( (i_cloud == ip_cloud_mix_max)    .OR.                          &
         (i_cloud == ip_cloud_mix_random) .OR.                          &
         (i_cloud == ip_cloud_part_corr) ) THEN

      IF ( (i >= n_cloud_top).AND.(i < n_layer) ) THEN

        DO l=1, n_profile
          IF (w_free(l, i+1) > tol_cloud) THEN
            tot_cloud_cover(l)=tot_cloud_cover(l)*w_free(l, i+1) /      &
              (1.0e+00_RealK - area_overlap(l, 2, 2))
          ELSE
            tot_cloud_cover(l) = 0.0e+00_RealK
          END IF
        END DO

      END IF

    ELSE IF ( (i_cloud == ip_cloud_triple)         .OR.                 &
              (i_cloud == ip_cloud_part_corr_cnv) ) THEN

      IF ( (i >= n_cloud_top).AND.(i < n_layer) ) THEN

        DO l=1, n_profile
          IF (w_free(l, i+1) > tol_cloud) THEN
            tot_cloud_cover(l)=tot_cloud_cover(l)*w_free(l, i+1) /      &
              (1.0e+00_RealK - area_overlap(l, 2, 2)                    &
                             - area_overlap(l, 3, 3))
          ELSE
            tot_cloud_cover(l) = 0.0e+00_RealK
          END IF
        END DO

      END IF

    END IF

!   Find the remaining areas of overlap on the assumption that
!   the overlap is random. We initialize the areas of overlap to
!   0 and reset later when such an area is present.
    DO k=1, n_region
      DO j=1, k-1
        DO l=1, n_profile
          area_overlap(l, k, j)=0.0e+00_RealK
          area_overlap(l, j, k)=0.0e+00_RealK
        END DO
      END DO
    END DO

    DO l=1, n_profile
      area_random_tot(l)=1.0e+00_RealK-area_overlap(l, 1, 1)
    END DO
    DO k=2, n_region
      DO l=1, n_profile
        area_random_tot(l)=area_random_tot(l)-area_overlap(l, k, k)
      END DO
    END DO
    DO k=1, n_region
      DO l=1, n_profile
        area_random_upper(l, k)                                         &
          =area_upper(l, k)-area_overlap(l, k, k)
        area_random_lower(l, k)                                         &
          =area_lower(l, k)-area_overlap(l, k, k)
      END DO
    END DO
!   To calculate the contributions of random overlap to the
!   areas of overlap we take the randomly overlapped portion
!   of the kth region in the upper layer. The probability that
!   this is overalpped with the randomly overlapped portion of
!   the jth region in the lower layer will be equal to
!   the randomly overlapped area of the lower jth region divided
!   by the total randomly overalpped area. The ratio might become
!   ill-conditioned for small amounts of cloud, the but this
!   should not be an issue as the randomly overalpped area would
!   then be small.
    DO k=1, n_region
      DO j=1, n_region
        DO l=1, n_profile
          IF (area_random_tot(l) >  tol_cloud) THEN
            area_overlap(l, k, j)=area_overlap(l, k, j)                 &
              +area_random_upper(l, k)                                  &
              *area_random_lower(l, j)/area_random_tot(l)
          END IF
        END DO
      END DO
    END DO

!   Now proceed to find the energy transfer coefficients
!   between the various regions.

!   Coefficients for the downward transfer of energy:

!   To avoid division by 0 we initialize to default values
!   and reset.
    DO k=1, n_region
      DO l=1, n_profile
        cloud_overlap(l, i, n_region*(k-1)+k)=1.0e+00_RealK
      END DO
      DO j=1, k-1
        DO l=1, n_profile
          cloud_overlap(l, i, n_region*(j-1)+k)=0.0e+00_RealK
          cloud_overlap(l, i, n_region*(k-1)+j)=0.0e+00_RealK
        END DO
      END DO
    END DO

    DO k=1, n_region
      DO l=1, n_profile
        IF (area_upper(l, k) >  tol_cloud) THEN
          DO j=1, n_region
            cloud_overlap(l, i, n_region*(j-1)+k)                       &
              =area_overlap(l, k, j)/area_upper(l, k)
          END DO
        END IF
      END DO
    END DO


!   Coefficients for upward flow of energy:

!   To avoid division by 0 we initialize to default values
!   and reset.
    DO k=1, n_region
      DO l=1, n_profile
        cloud_overlap(l, i, n_region*(n_region+k-1)+k)                  &
          =1.0e+00_RealK
      END DO
      DO j=1, k-1
        DO l=1, n_profile
          cloud_overlap(l, i, n_region*(n_region+j-1)+k)                &
            =0.0e+00_RealK
          cloud_overlap(l, i, n_region*(n_region+k-1)+j)                &
            =0.0e+00_RealK
        END DO
      END DO
    END DO

    DO k=1, n_region
      DO l=1, n_profile
        IF (area_lower(l, k) >  tol_cloud) THEN
          DO j=1, n_region
            cloud_overlap(l, i, n_region*(n_region+j-1)+k)              &
              =area_overlap(l, j, k)/area_lower(l, k)
          END DO
        END IF
      END DO
    END DO


!   Reassign the fractions in the upper layer to step down
!   through the atmosphere.
    IF (i <  n_layer) THEN
      DO k=1, n_region
        DO l=1, n_profile
          area_upper(l, k)=area_lower(l, k)
        END DO
      END DO
    END IF

  END DO

  DO l=1, n_profile
     tot_cloud_cover(l)=1.0e+00_RealK-tot_cloud_cover(l)
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE overlap_coupled
END MODULE overlap_coupled_mod
