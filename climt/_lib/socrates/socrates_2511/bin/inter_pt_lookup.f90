! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate pressure, temperature and gas fraction interpolation
! factor for k-terms.
!
!-----------------------------------------------------------------------

MODULE inter_pt_lookup_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INTER_PT_LOOKUP_MOD'
CONTAINS
SUBROUTINE inter_pt_lookup(nd_profile, nd_layer, nd_pre, nd_tmp         &
     , nd_gas_frac, nd_species, nd_species_sb                           &
     , n_profile, n_layer, n_gas_frac, n_absorb                         &
     , type_absorb, l_self_broadening, index_sb                         &
     , l_water, l_mixing_ratio, i_pointer_water                         &
     , p, t, gas_mix_ratio                                              &
     , p_lookup, t_lookup, gf_lookup                                    &
     , fac00, fac01, fac10, fac11, jp, jt, jtt, fgf, jgf, jgfp1)

  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: mol_weight_air
  USE gas_list_pcf, ONLY: molar_weight, ip_h2o
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


  INTEGER, INTENT(IN) ::                                                &
       nd_profile                                                       &
!       Max number of profile
     , nd_layer                                                         &
!       Max number of layer
     , nd_pre                                                           &
!       Number of lookup pressures
     , nd_tmp                                                           &
!       Number of lookup temperatures
     , nd_gas_frac                                                      &
!       Maximum number of gas fractions
     , nd_species                                                       &
!       Maximum number of species
     , nd_species_sb                                                    &
!       Maximum number of species with self-broadening
     , n_layer                                                          &
!       Number of layers
     , n_profile                                                        &
!       Number of profiles
     , n_gas_frac                                                       &
!       Number of gas fractions
     , n_absorb                                                         &
!       Number of absorbers
     , type_absorb(nd_species)                                          &
!       Actual types of absorbers
     , index_sb(nd_species)                                             &
!       Gas index numbers in arrays with self-broadening
     , i_pointer_water
!       Pointer to water vapour

  LOGICAL, INTENT(IN) :: &
       l_self_broadening(nd_species)                                    &
!       Flag for self-broadening of species
     , l_water                                                          &
!       Flag for water present in spectral file
     , l_mixing_ratio
!       True if mixing ratios are with respect to dry mass

  REAL (RealK), INTENT(IN) ::                                           &
       p(nd_profile, nd_layer)                                          &
!          Actual pressure
     , t(nd_profile, nd_layer)                                          &
!          Layer temperature
     , gas_mix_ratio(nd_profile, nd_layer, nd_species)                  &
!          Mass mixing ratios of gases
     , p_lookup(nd_pre)                                                 &
!          Lookup table pressures
     , t_lookup(nd_tmp, nd_pre)                                         &
!          Lookup table temperatures
     , gf_lookup(nd_gas_frac)
!          Lookup table gas fractions

  REAL (RealK), INTENT(OUT) ::                                          &
       fac00(nd_profile, nd_layer)                                      &
     , fac01(nd_profile, nd_layer)                                      &
     , fac10(nd_profile, nd_layer)                                      &
     , fac11(nd_profile, nd_layer)                                      &
!          Multiplication factors for P & T interpolation
     , fgf(nd_profile, nd_layer, nd_species_sb)
!          Multiplication factors for gas fraction interpolation

  INTEGER, INTENT(OUT) ::                                               &
       jp(nd_profile, nd_layer)                                         &
!       Index of reference pressure level such that the actual
!       pressure is between JP and JP+1
     , jt(nd_profile, nd_layer)                                         &
!       Index of reference temperature at level I such that the actual
!       temperature is between JT and JT+1
     , jtt(nd_profile, nd_layer)                                        &
!       Index of reference temperature at level I+1 such that the actual
!       temperature is between JTT and JTT+1
     , jgf(nd_profile, nd_layer, nd_species_sb)                         &
!       Index of reference gas fraction such that the actual gas
!       fraction is between jgf and jgf+1
     , jgfp1(nd_profile, nd_layer, nd_species_sb)
!       Equal to jgf + 1, except if n_gas_frac == 1, in which case it
!       is equal to jgf to avoid out of bounds errors
!       in interpolation.

! Local variables
  REAL (RealK) ::                                                       &
       fp, ft, ftt, compfp                                              &
!          Fraction factor for P & T interpolation
     , p_layer, t_layer                                                 &
!          Inter-medium variable
     , p_lookup_min, p_lookup_max                                       &
!          Minimum and maximum values for look-up pressure
     , gf_lookup_min, gf_lookup_max                                     &
!          Minimum and maximum values for look-up gas fraction
     , gf_layer(nd_profile, nd_layer)                                   &
!          Gas fraction of gas in layer
     , gf_layer_loc                                                     &
!          Temporary storage of layer gas fraction
     , ratio_gas                                                        &
!          Ratio of molar weight of air to molar weight of gas
     , ratio_water
!          Ratio of molar weight of water to molar weight of gas

  INTEGER ::                                                            &
       i, l, j                                                          &
!       Vertical, horizontal and gas loop index
     , j_sb
!       Gas loop index in arrays with self-broadened

  REAL (RealK), PARAMETER :: eps = EPSILON(1.0_RealK)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='INTER_PT_LOOKUP'

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  p_lookup_min = p_lookup(1)      + MAX(ABS(p_lookup(1)     )*eps, eps)
  p_lookup_max = p_lookup(nd_pre) - MAX(ABS(p_lookup(nd_pre))*eps, eps)

  DO i=1, n_layer
    DO l=1, n_profile
!     Find the reference pressure on the lower side of the layer
!     pressure and store in JP. Store in FP the fraction
!     of the difference (in ln(pressure)) between JP and JP+1
!     that the layer pressure lies for gas absorption
!     coefficient interpolation.
      p_layer = MIN( MAX( LOG(p(l,i)), p_lookup_min ), p_lookup_max )

      jp(l,i) = MINLOC( p_layer - p_lookup, 1, &
                        p_layer - p_lookup >= 0.0_RealK)

      fp = ( p_layer             - p_lookup(jp(l,i)) ) &
         / ( p_lookup(jp(l,i)+1) - p_lookup(jp(l,i)) )

!     Find the reference temperature on the lower side of the
!     layer temperature for each of the layers JP and JP+1. Store
!     these indices in JT and JTT, resp.
!     Store in FT (resp. FTT) the fraction of the way between JT
!     (JTT) and the next highest reference temperature that the
!     layer temperature falls.
      t_layer=MIN( MAX( t(l,i), t_lookup(1,jp(l,i))*(1.0_RealK+eps) ), &
                   t_lookup(nd_tmp,jp(l,i))*(1.0_RealK-eps) )

      jt(l,i) = MINLOC( t_layer - t_lookup(:,jp(l,i)), 1, &
                        t_layer - t_lookup(:,jp(l,i)) >= 0.0_RealK)

      ft = ( t_layer                     - t_lookup(jt(l,i),jp(l,i)) ) &
         / ( t_lookup(jt(l,i)+1,jp(l,i)) - t_lookup(jt(l,i),jp(l,i)) )

      t_layer=MIN( MAX( t(l,i), t_lookup(1,jp(l,i)+1)*(1.0_RealK+eps) ), &
                   t_lookup(nd_tmp,jp(l,i)+1)*(1.0_RealK-eps) )

      jtt(l,i) = MINLOC( t_layer - t_lookup(:,jp(l,i)+1), 1, &
                         t_layer - t_lookup(:,jp(l,i)+1) >= 0.0_RealK)

      ftt=(t_layer                       -t_lookup(jtt(l,i),jp(l,i)+1)) &
        / (t_lookup(jtt(l,i)+1,jp(l,i)+1)-t_lookup(jtt(l,i),jp(l,i)+1))

!     Multiply the pressure fraction with the appropriate temperature
!     fraction for use in the interpolations
      compfp=1.0_RealK-fp
      fac00(l,i)=compfp*(1.0_RealK-ft)
      fac10(l,i)=fp*(1.0_RealK-ftt)
      fac01(l,i)=compfp*ft
      fac11(l,i)=fp*ftt
    END DO
  END DO

! Find the gas fraction interpolation parameters for all gases with
! self-broadening.
  IF (ANY(l_self_broadening(1:n_absorb))) THEN
    gf_lookup_min = gf_lookup(1) + MAX(gf_lookup(1)*eps, eps)
    gf_lookup_max = gf_lookup(n_gas_frac) - MAX(gf_lookup(n_gas_frac)*eps, eps)
    DO j=1, n_absorb
      IF (l_self_broadening(j)) THEN
        j_sb=index_sb(j)

    !   Find gas fraction in the layer.
        ratio_gas = mol_weight_air/(molar_weight(type_absorb(j))*1.0E-03_RealK)
        IF (l_water) THEN
          ratio_water = mol_weight_air/(molar_weight(ip_h2o)*1.0E-03_RealK)
          IF (l_mixing_ratio) THEN
            DO i=1, n_layer
              DO l=1, n_profile
                gf_layer(l,i) = gas_mix_ratio(l,i,j)*ratio_gas &
                  /(1.0_RealK + gas_mix_ratio(l,i,i_pointer_water) &
                  * ratio_water)
              END DO
            END DO
          ELSE
            DO i=1, n_layer
              DO l=1, n_profile
                gf_layer(l,i) = gas_mix_ratio(l,i,j)*ratio_gas &
                  /(1.0_RealK + gas_mix_ratio(l,i,i_pointer_water) &
                  * (ratio_water - 1.0_RealK))
              END DO
            END DO
          END IF
        ELSE
          DO i=1, n_layer
            DO l=1, n_profile
              gf_layer(l,i) = gas_mix_ratio(l,i,j)*ratio_gas
            END DO
          END DO
        END IF

    !   Find the interpolation parameters.
        IF (n_gas_frac == 1) THEN
!         Use only gas fraction available.
          DO i=1, n_layer
            DO l=1, n_profile
              jgf(l,i,j_sb) = 1
              jgfp1(l,i,j_sb) = 1
              fgf(l,i,j_sb) = 1.0_RealK
            END DO
          END DO
        ELSE
!         Linear interpolation.
          DO i=1, n_layer
            DO l=1, n_profile
              gf_layer_loc = MIN( MAX( gf_layer(l,i), gf_lookup_min), &
                                  gf_lookup_max)

              jgf(l,i,j_sb) = &
                MINLOC( gf_layer_loc - gf_lookup(1:n_gas_frac), 1, &
                        gf_layer_loc - gf_lookup(1:n_gas_frac) >= 0.0_RealK)
              jgfp1(l,i,j_sb) = jgf(l,i,j_sb) + 1
              fgf(l,i,j_sb) = 1.0_RealK - &
                ( gf_layer_loc - gf_lookup(jgf(l,i,j_sb)) ) &
                / ( gf_lookup(jgf(l,i,j_sb)+1) - gf_lookup(jgf(l,i,j_sb)) )
            END DO
          END DO
        END IF
      END IF
    END DO
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE inter_pt_lookup
END MODULE inter_pt_lookup_mod
