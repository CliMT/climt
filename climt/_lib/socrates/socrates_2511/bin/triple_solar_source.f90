! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set the solar solar terms in a triple column.
!
! Method:
!   The direct beam is calculated by propagating down through
!   the column. These direct fluxes are used to define the
!   source terms in each layer.
!
!- ---------------------------------------------------------------------
MODULE triple_solar_source_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TRIPLE_SOLAR_SOURCE_MOD'
CONTAINS
SUBROUTINE triple_solar_source(control, bound                           &
     , n_profile, n_layer, n_cloud_top                                  &
     , n_region, flux_inc_direct                                        &
     , l_scale_solar, adjust_solar_ke                                   &
     , trans_0_dir, trans_0, source_coeff                               &
     , v11, v12, v13, v21, v22, v23, v31, v32, v33                      &
     , flux_direct                                                      &
     , flux_direct_ground                                               &
     , s_up, s_down                                                     &
     , nd_profile, nd_layer, id_ct, nd_source_coeff, nd_region          &
     )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_bound,   ONLY: StrBound
  USE rad_pcf, ONLY: ip_direct_csr_scaling, ip_direct_noscaling,        &
                     ip_region_clear, ip_region_conv, ip_region_strat,  &
                     ip_scf_solar_down, ip_scf_solar_up
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Boundary conditions:
  TYPE(StrBound),     INTENT(IN)    :: bound

! Sizes of dummy arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_source_coeff                                                   &
!       Size allocated for source coefficients
    , nd_region
!       Maximum number of cloudy regions


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_cloud_top                                                       &
!       Top cloudy layer
    , n_region
!       Number of cloudy regions

! Special arrays for equivalent extinction:
  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar
!       Scaling applied to solar flux
  REAL (RealK), INTENT(IN) ::                                           &
      adjust_solar_ke(nd_profile, nd_layer)
!       Adjustment to solar fluxes with equivalent extinction

  REAL (RealK), INTENT(IN) ::                                           &
      flux_inc_direct(nd_profile)
!       Incident direct solar flux

! Optical properties:
  REAL (RealK), INTENT(IN) ::                                           &
      trans_0(nd_profile, nd_layer, nd_region)                          &
!       Direct transmission
    , trans_0_dir(nd_profile, nd_layer, nd_region)                      &
!       Direct transmission using direct tau
    , source_coeff(nd_profile, nd_layer                                 &
      , nd_source_coeff, nd_region)
!       Source coefficients

! Energy transfer coefficients:
  REAL (RealK), INTENT(IN) ::                                           &
      v11(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , v12(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , v13(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , v21(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , v22(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , v23(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , v31(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , v32(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient
    , v33(nd_profile, id_ct-1: nd_layer)
!       Energy transfer coefficient

! Calculated direct flux and source terms:
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_profile, 0: nd_layer)                              &
!       Overall direct flux
    , flux_direct_ground(nd_profile, nd_region)                         &
!       Direct fluxes at ground beneath each region
    , s_up(nd_profile, nd_layer, nd_region)                             &
!       Upward source functions
    , s_down(nd_profile, nd_layer, nd_region)
!       Downward source functions


! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , k
!       Loop variable

  REAL (RealK) ::                                                       &
      solar_top(nd_profile, nd_region)                                  &
!       Solar fluxes at top of layer
    , solar_base(nd_profile, nd_region)                                 &
!       Solar fluxes at base of layer
    , flux_direct_dir(nd_profile, 0: nd_layer)
!       Overall unscaled direct flux

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIPLE_SOLAR_SOURCE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The clear and cloudy direct fluxes are calculated separately
! and added together to form the total direct flux.

! Set incident fluxes.
  DO l=1, n_profile
    flux_direct(l, 0)=flux_inc_direct(l)
  END DO
  IF (control%i_direct_tau == ip_direct_noscaling .OR.                 &
      control%i_direct_tau == ip_direct_csr_scaling) THEN
    DO l=1, n_profile
      flux_direct_dir(l, 0)=flux_inc_direct(l)
    END DO
  END IF

! With equivalent extinction the direct solar flux must be corrected.
  IF (l_scale_solar) THEN

    DO i=1, n_cloud_top-1
      DO l=1, n_profile
        flux_direct(l, i)                                               &
          =flux_direct(l, i-1)*trans_0(l, i, ip_region_clear)           &
          *adjust_solar_ke(l, i)
        s_up(l, i, ip_region_clear)                                     &
          =source_coeff(l, i, ip_scf_solar_up, ip_region_clear)         &
          *flux_direct(l, i-1)
        s_down(l, i, ip_region_clear)                                   &
          =(source_coeff(l, i, ip_scf_solar_down, ip_region_clear)      &
          -trans_0(l, i, ip_region_clear))*flux_direct(l, i-1)          &
          +flux_direct(l, i)
      END DO
    END DO

    IF (control%i_direct_tau == ip_direct_noscaling .OR.                &
        control%i_direct_tau == ip_direct_csr_scaling) THEN
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          flux_direct_dir(l, i)                                         &
            =flux_direct_dir(l, i-1)                                    &
            *trans_0_dir(l, i, ip_region_clear)                         &
            *adjust_solar_ke(l, i)
        END DO
      END DO
    END IF

  ELSE

    DO i=1, n_cloud_top-1
      DO l=1, n_profile
        flux_direct(l, i)                                               &
          =flux_direct(l, i-1)*trans_0(l, i, ip_region_clear)
        s_up(l, i, ip_region_clear)                                     &
          =source_coeff(l, i, ip_scf_solar_up, ip_region_clear)         &
          *flux_direct(l, i-1)
        s_down(l, i, ip_region_clear)                                   &
          =source_coeff(l, i, ip_scf_solar_down, ip_region_clear)       &
          *flux_direct(l, i-1)
      END DO
    END DO

    IF (control%i_direct_tau == ip_direct_noscaling .OR.                &
        control%i_direct_tau == ip_direct_csr_scaling) THEN
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          flux_direct_dir(l, i)                                         &
            =flux_direct_dir(l, i-1)                                    &
            *trans_0_dir(l, i, ip_region_clear)
        END DO
      END DO
    END IF

  END IF



! Clear and cloudy region.
! Initialize partial fluxes:
  DO l=1, n_profile
    solar_base(l, ip_region_clear)=flux_direct(l, n_cloud_top-1)
    solar_base(l, ip_region_strat)=0.0e+00_RealK
    solar_base(l, ip_region_conv)=0.0e+00_RealK
  END DO


  DO i=n_cloud_top, n_layer

!   Transfer fluxes across the interface.

    DO l=1, n_profile
      solar_top(l, ip_region_clear)                                     &
        =v11(l, i-1)*solar_base(l, ip_region_clear)                     &
        +v12(l, i-1)*solar_base(l, ip_region_strat)                     &
        +v13(l, i-1)*solar_base(l, ip_region_conv)
      solar_top(l, ip_region_strat)                                     &
        =v21(l, i-1)*solar_base(l, ip_region_clear)                     &
        +v22(l, i-1)*solar_base(l, ip_region_strat)                     &
        +v23(l, i-1)*solar_base(l, ip_region_conv)
      solar_top(l, ip_region_conv)                                      &
        =v31(l, i-1)*solar_base(l, ip_region_clear)                     &
        +v32(l, i-1)*solar_base(l, ip_region_strat)                     &
        +v33(l, i-1)*solar_base(l, ip_region_conv)
    END DO


!   Propagate the fluxes through the layer:
    IF (l_scale_solar) THEN

      DO k=1, n_region
        DO l=1, n_profile
          solar_base(l, k)                                              &
            =solar_top(l, k)                                            &
            *trans_0(l, i, k)*adjust_solar_ke(l, i)
          s_up(l, i, k)                                                 &
            =source_coeff(l, i, ip_scf_solar_up, k)                     &
            *solar_top(l, k)
          s_down(l, i, k)                                               &
            =(source_coeff(l, i, ip_scf_solar_down, k)                  &
            -trans_0(l, i, k))*solar_top(l, k)                          &
            +solar_base(l, k)
        END DO
      END DO

    ELSE

      DO k=1, n_region
        DO l=1, n_profile
          solar_base(l, k)=solar_top(l, k)                              &
            *trans_0(l, i, k)
          s_up(l, i, k)                                                 &
            =source_coeff(l, i, ip_scf_solar_up, k)                     &
            *solar_top(l, k)
          s_down(l, i, k)                                               &
            =source_coeff(l, i, ip_scf_solar_down, k)                   &
            *solar_top(l, k)
        END DO
      END DO

    END IF


!   Calculate the total direct flux.
    DO l=1, n_profile
      flux_direct(l, i)=solar_base(l, ip_region_clear)                  &
        +solar_base(l, ip_region_strat)                                 &
        +solar_base(l, ip_region_conv)
    END DO

  END DO


! Repeat for direct flux using trans without scaling
  IF (control%i_direct_tau == ip_direct_noscaling .OR.                  &
      control%i_direct_tau == ip_direct_csr_scaling) THEN

!   Clear and cloudy region.
!   Initialize partial fluxes:
    DO l=1, n_profile
      solar_base(l, ip_region_clear)=flux_direct_dir(l, n_cloud_top-1)
      solar_base(l, ip_region_strat)=0.0e+00_RealK
      solar_base(l, ip_region_conv)=0.0e+00_RealK
    END DO
    
    
    DO i=n_cloud_top, n_layer
    
!     Transfer fluxes across the interface.
    
      DO l=1, n_profile
        solar_top(l, ip_region_clear)                                   &
          =v11(l, i-1)*solar_base(l, ip_region_clear)                   &
          +v12(l, i-1)*solar_base(l, ip_region_strat)                   &
          +v13(l, i-1)*solar_base(l, ip_region_conv)
        solar_top(l, ip_region_strat)                                   &
          =v21(l, i-1)*solar_base(l, ip_region_clear)                   &
          +v22(l, i-1)*solar_base(l, ip_region_strat)                   &
          +v23(l, i-1)*solar_base(l, ip_region_conv)
        solar_top(l, ip_region_conv)                                    &
          =v31(l, i-1)*solar_base(l, ip_region_clear)                   &
          +v32(l, i-1)*solar_base(l, ip_region_strat)                   &
          +v33(l, i-1)*solar_base(l, ip_region_conv)
      END DO
    
    
!     Propagate the fluxes through the layer:
      IF (l_scale_solar) THEN
    
        DO k=1, n_region
          DO l=1, n_profile
            solar_base(l, k)                                            &
              =solar_top(l, k)                                          &
              *trans_0_dir(l, i, k)*adjust_solar_ke(l, i)
          END DO
        END DO
    
      ELSE
    
        DO k=1, n_region
          DO l=1, n_profile
            solar_base(l, k)=solar_top(l, k)                            &
              *trans_0_dir(l, i, k)
          END DO
        END DO
    
      END IF
    
    
!     Calculate the total direct flux.
      DO l=1, n_profile
        flux_direct_dir(l, i)=solar_base(l, ip_region_clear)            &
          +solar_base(l, ip_region_strat)                               &
          +solar_base(l, ip_region_conv)
      END DO
    
    END DO


!   From this point, use the unscaled direct flux as the direct component.
    DO i= 0, n_layer
      DO l=1, n_profile
        flux_direct(l, i)=flux_direct_dir(l, i)
      END DO
    END DO

  END IF


! Pass the last value at the base of the cloud out.
  DO k=1, n_region
    DO l=1, n_profile
      flux_direct_ground(l, k)=solar_base(l, k)
    END DO
  END DO


! Correct the direct flux at the ground for sloping terrain
  IF (control%l_orog) THEN
     flux_direct(1:n_profile, n_layer) =                                &
        flux_direct(1:n_profile, n_layer) *                             &
        bound%orog_corr(1:n_profile)

     DO k=1, n_region
        flux_direct_ground(1:n_profile, k) =                            &
           flux_direct_ground(1:n_profile, k) *                         &
           bound%orog_corr(1:n_profile)

        s_down(1:n_profile, n_layer, k) =                               &
           s_down(1:n_profile, n_layer, k) +                            &
           solar_base(1:n_profile, k) *                                 &
           (bound%orog_corr(1:n_profile) - 1.0_RealK)

     END DO
  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE triple_solar_source
END MODULE triple_solar_source_mod
