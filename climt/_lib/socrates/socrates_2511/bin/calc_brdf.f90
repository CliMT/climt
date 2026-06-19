! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate the arrays of BRDF terms.
!
! Purpose:
!   This routine is called to calculate a set of arrays related to
!   the BRDF for later efficiency.
!
! Method:
!   As this routine is called only once speed is not too critical
!   so direct calculation is used.
!
! Symmetries of the BRDF and storage:
!   Since the BRDF is defined only for downwelling incident
!   radiances and upwelling reflected radiances it cannot be
!   uniquely defined as a double expension in spherical harmonics.
!   To make a unique expansion we stipulate that only harmonics of
!   even parity will be used: if odd harmonics were chosen we would
!   get into difficulties with the Gibb's phenomenon, as all odd
!   harmonics vanish on the horizontal.
!   F(j, l, l', m) will therefore be 0 unless l+m and l'+m are
!   both even, so the indices of storage for l and l' are set to
!   l/2 and l'/2. There are more efficient schemes of storage
!   that depend on the fact that F vanishes if l<m or l'<m, but
!   such a scheme has not been selected because of its extra
!   complexity.
!
!- ---------------------------------------------------------------------
MODULE calc_brdf_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CALC_BRDF_MOD'
CONTAINS
SUBROUTINE calc_brdf(isolir, ms_min, ms_max                             &
    , ia_sph_mm                                                         &
    , uplm_sol, uplm_zero                                               &
    , n_brdf_basis_fnc, ls_brdf_trunc, f_brdf                           &
    , n_profile, n_direction, direction                                 &
    , brdf_sol, brdf_hemi                                               &
    , nd_profile, nd_radiance_profile, nd_direction                     &
    , nd_max_order, nd_sph_coeff                                        &
    , nd_brdf_basis_fnc, nd_brdf_trunc                                  &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_solar, ip_infra_red
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE eval_uplm_mod, ONLY: eval_uplm

  IMPLICIT NONE


! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_radiance_profile                                               &
!       Size allocated for profiles where radiances are calculated
    , nd_direction                                                      &
!       Size allocated for viewing directions
    , nd_max_order                                                      &
!       Size allcoated for polar orders
    , nd_sph_coeff                                                      &
!       Size allocated for spherical coefficients
    , nd_brdf_basis_fnc                                                 &
!       Size allocated for BRDF basis functions
    , nd_brdf_trunc
!       Size allocated for truncation of BRDFs


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile
!       Number of atmospheric profiles
  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Spectral region
  INTEGER, INTENT(IN) ::                                                &
      ms_min                                                            &
!       Lowest azimuthal order calculated
    , ms_max                                                            &
!       Highest azimuthal order calculated
    , ia_sph_mm(0: nd_max_order)
!       Address of spherical coefficient for (m, m) for each m
  REAL (RealK), INTENT(IN) ::                                           &
      uplm_zero(nd_sph_coeff)
!       Array of Upsilon_l^m and derivatives at polar angles of pi/2
  INTEGER, INTENT(IN) ::                                                &
      n_brdf_basis_fnc                                                  &
!       Number of basis functions for BRDFs
    , ls_brdf_trunc
!       Order of truncation applied to BRDFs
  REAL (RealK), INTENT(IN) ::                                           &
      f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                      &
        , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)
!       Array of moments of BRDF basis functions
  INTEGER, INTENT(IN) ::                                                &
      n_direction
!       Number of viewing directions
  REAL (RealK), INTENT(IN) ::                                           &
      direction(nd_radiance_profile, nd_direction, 2)                   &
!       Cosines of polar viewing angles and actual azimuthal
!       viewing angles
    , uplm_sol(nd_profile, nd_sph_coeff)
!       Upsilon terms for solar radiation

  REAL (RealK), INTENT(OUT) ::                                          &
      brdf_sol(nd_profile, nd_brdf_basis_fnc, nd_direction)             &
!       The BRDF evaluated for scattering from the solar
!       beam into the viewing direction
    , brdf_hemi(nd_profile, nd_brdf_basis_fnc, nd_direction)
!       The BRDF evaluated for scattering from isotropic
!       radiation into the viewing direction


! Local variables
  INTEGER                                                               &
      ls                                                                &
!       Polar order of harmonic
    , lsr                                                               &
!       Reduced polar order of harmonic
    , ls_p                                                              &
!       Polar order of harmonic
    , lsr_p                                                             &
!       Reduced polar order of harmonic
    , ms                                                                &
!       Azimuthal order of spherical harmonic
    , j                                                                 &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , id
!       Loop variable

  REAL (RealK) ::                                                       &
      up_lm(nd_profile, nd_brdf_trunc+1, nd_direction)                  &
!       Polar parts of spherical harmonics in the viewing
!       directions
    , ss1(nd_profile, 0: nd_brdf_trunc)                                 &
!       Products of the BRDF and the solar harmonics
    , azim_factor(nd_profile, nd_direction)                             &
!       Azimuthal factors
    , kappa(nd_brdf_trunc+1)                                            &
!       Hemispherical quadratic integrals of spherical harmonics
!       (reduced storage does not seem worth the effort here)
    , fk(nd_brdf_trunc+1)
!       Sum of products of the BRDF and KAPPA over l'

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_BRDF'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (isolir == ip_solar) THEN

!   Initialize the BRDF for solar radiation
    DO id=1, n_direction
      DO j=1, n_brdf_basis_fnc
        DO l=1, n_profile
          brdf_sol(l, j, id)=0.0e+00_RealK
        END DO
      END DO
    END DO

!   Loop over azimuthal orders.
    DO ms=ms_min, MIN(ms_max, ls_brdf_trunc)

!     Caclulate the azimuthal factors.
      IF (ms == 0) THEN
        DO id=1, n_direction
          DO l=1, n_profile
            azim_factor(l, id)=1.0e+00_RealK
          END DO
        END DO
      ELSE
        DO id=1, n_direction
          DO l=1, n_profile
            azim_factor(l, id)                                          &
            =2.0e+00_RealK*COS(REAL(ms,RealK)*direction(l,id,2))
          END DO
        END DO
      END IF

!     Calculate spherical harmonics in the viewing directions
!     at this azimuthal order.
      DO id=1, n_direction
        CALL eval_uplm(ms, ls_brdf_trunc                                &
          , n_profile, direction(1, id, 1), up_lm(1, 1, id)             &
          , nd_profile)
      END DO

!     Now loop over basis functions.
      DO j=1, n_brdf_basis_fnc

!       The array SS1 pulls in the solar dependence, which is
!       independent of the viewing direction. At this stage both
!       MS and J are fixed.
        DO ls=ms, ls_brdf_trunc, 2
          DO l=1, n_profile
            ss1(l, ls)=f_brdf(j, ls/2, ms/2, ms)                        &
             *uplm_sol(l, ia_sph_mm(ms))
          END DO

          DO ls_p=ms+2, ls_brdf_trunc, 2
            DO l=1, n_profile
              ss1(l, ls)=f_brdf(j, ls/2, ls_p/2, ms)                    &
               *uplm_sol(l, ia_sph_mm(ms)+ls_p-ms)
            END DO
          END DO
        END DO

!       Now consider each direction, incrementing the solar
!       BRDF.
        DO id=1, n_direction
          DO ls=ms, ls_brdf_trunc, 2
            DO l=1, n_profile
              brdf_sol(l, j, id)=brdf_sol(l, j, id)                     &
                +ss1(l, ls)*up_lm(l, ls+1-ms, id)                       &
                *azim_factor(l, id)
            END DO
          END DO
        END DO

      END DO

    END DO

  ELSE IF (isolir == ip_infra_red) THEN

!   Initialize.
    DO id=1, n_direction
      DO j=1, n_brdf_basis_fnc
        DO l=1, n_profile
          brdf_hemi(l, j, id)=0.0e+00_RealK
        END DO
      END DO
    END DO

!   Only azimuthally symmetric terms contribute.
    DO ms=ms_min, 0

      DO lsr_p=1, ls_brdf_trunc-ms+1, 2
        kappa(lsr_p)=2.0e+00_RealK*pi                                   &
          *uplm_zero(ia_sph_mm(0)+lsr_p-1)*uplm_zero(2)                 &
          /REAL((lsr_p-2)*(lsr+1+2*ms), RealK)
      END DO

!     Calculate spherical harmonics in the viewing directions
!     at this azimuthal order.
      DO id=1, n_direction
        CALL eval_uplm(ms, ls_brdf_trunc                                &
          , n_profile, direction(1, id, 1), up_lm(1, 1, id)             &
          , nd_profile)
      END DO

!     Now loop over basis functions.
      DO j=1, n_brdf_basis_fnc

        DO lsr=1, ls_brdf_trunc-ms+1, 2
          fk(lsr)=0.0e+00_RealK
          DO lsr_p=1, ls_brdf_trunc-ms+1, 2
            fk(lsr)=fk(lsr)+kappa(lsr_p)                                &
              *f_brdf(j, (lsr-1+ms)/2, (lsr_p-1+ms)/2, ms)
          END DO
          DO id=1, n_direction
            DO l=1, n_profile
              brdf_hemi(l, j, id)=brdf_hemi(l, j, id)                   &
                +fk(lsr)*up_lm(l, lsr, id)
            END DO
          END DO
        END DO

      END DO

    END DO

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_brdf
END MODULE calc_brdf_mod
