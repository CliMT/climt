! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Interpolation routine for ESFT onto model grid
!
! Method:
!   Interpolation of gas and continuum absorption coefficients.
!
!-----------------------------------------------------------------------
MODULE inter_k_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INTER_K_MOD'
CONTAINS
SUBROUTINE inter_k(n_profile, n_layer, n_band_absorb                    &
  , mix_gas_band, n_mix_gas, index_mix_gas                              &
  , i_band_esft, i_band                                                 &
  , n_continuum, k_continuum, l_continuum, index_continuum              &
  , k_h2oc, fac00, fac01, fac10, fac11                                  &
  , fac00c, fac01c, fac10c, fac11c, facc00, facc01                      &
  , jp, jph2oc, jt, jtt, jto2c, jtswo3                                  &
  , k_esft, k_mix_gas                                                   &
  , f_mix, gas_mix_ratio, gas_mix_amt                                   &
  , k_esft_layer, k_mix_gas_layer, k_contm_layer                        &
! dimensions
  , nd_profile, nd_layer                                                &
  , nd_band, nd_species, nd_continuum                                   &
  , nd_esft_term, nd_mix, nd_tmp, nd_pre                                &
  , nd_band_mix_gas                                                     &
  )

  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_cont_h2o, ip_cont_o2n2_2
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

! Dummy array sizes
  INTEGER, INTENT(IN) ::                                                &
       nd_profile                                                       &
!         Maximum number of profiles
     , nd_layer                                                         &
!         Maximum number of layers
     , nd_band                                                          &
!         Number of bands
     , nd_species                                                       &
!         Number of species
     , nd_continuum                                                     &
!         Number of continua
     , nd_mix                                                           &
!         Number of binary interpolation parameter for mixture species
     , nd_esft_term                                                     &
!         maximum number of esft terms
     , nd_tmp                                                           &
!         number of reference temperatures for k-distribution
     , nd_pre                                                           &
!         number of reference pressure for k-distribution
     , nd_band_mix_gas
!         Number of bands where mixture gases exist
  REAL (RealK) ::                                                       &
       k_esft(nd_tmp, nd_pre, nd_esft_term, nd_species)                 &
!         Absorption coefficients at reference conditions
     , k_h2oc(nd_pre,nd_tmp,nd_esft_term)                               &
!         H2O continuum Absorption coefficients at reference conditions
     , k_mix_gas(nd_pre, nd_tmp, nd_mix, nd_esft_term, nd_band_mix_gas) &
!         Absorption coefficients of minor species at reference conditions
     , f_mix                                                            &
!         Parameter for determining mixed gas amount
     , k_esft_layer(nd_profile, nd_layer, nd_esft_term, nd_species)     &
!         Interpolated absorption coefficients
     , k_mix_gas_layer(nd_profile, nd_esft_term, nd_layer)              &
!         Interpolated absorption coefficients of mixed species
     , k_contm_layer(nd_profile, nd_esft_term, nd_layer, nd_continuum)  &
!         Interpolated continuum absorption coefficients
     , k_continuum(nd_esft_term, 5, nd_band, nd_continuum)              &
!         Reference continuum absorption coefficients
     , gas_mix_ratio(nd_profile, nd_layer, nd_species)                  &
!        Gaseous mass mixing ratios
     , gas_mix_amt(nd_profile, nd_layer)                                &
!        Gas mixing ratio for AER overlap scheme
     , fac00(nd_profile, nd_layer)                                      &
     , fac01(nd_profile, nd_layer)                                      &
     , fac10(nd_profile, nd_layer)                                      &
     , fac11(nd_profile, nd_layer)                                      &
!         Multiplication factors for P & T interpolation
     , fac00c(nd_profile, nd_layer)                                     &
     , fac01c(nd_profile, nd_layer)                                     &
     , fac10c(nd_profile, nd_layer)                                     &
     , fac11c(nd_profile, nd_layer)                                     &
!         Multiplication factors for P & T interpolation
     , facc00(nd_profile, nd_layer)                                     &
     , facc01(nd_profile, nd_layer)
!         Multiplication factors for continuum T interpolation

  INTEGER ::                                                            &
       n_layer                                                          &
!       Number of layers
     , n_profile                                                        &
!       Number of profiles
     , n_band_absorb(nd_band)                                           &
!       Number of gases in each band
     , mix_gas_band(nd_band)                                            &
!       Sequence band number (not real band number) of mixed species
     , index_mix_gas(2,nd_band_mix_gas)                                 &
!       Index of mixed species
     , n_mix_gas                                                        &
!       Indexing band where mixed species occurs
     , n_continuum                                                      &
!       Number of continuum
     , index_continuum(nd_band, nd_continuum)                           &
!       Indices of continua
     , i_band                                                           &
!        band number
     , i_band_esft                                                      &
!        Number of ESFT term in each band
     , jp(nd_profile, nd_layer), jp1                                    &
!        Index of reference pressure level such that the actual
!        pressure is between JP and JP1
     , jph2oc(nd_profile, nd_layer)                                     &
!        same as jp but for water vapour pressure
     , jt(nd_profile, nd_layer)                                         &
!        Index of reference temperature at level i such that the actual
!        temperature is between JT and JT+1
     , jtt(nd_profile, nd_layer)                                        &
!        Index of reference temperature i+1 such that the actual
!        temperature is between JT and JT+1
     , jto2c(nd_profile, nd_layer)                                      &
!        Index of continuum  reference temperature at level I
!        such that the actual temperature is between JTO2C and JTO2C+1
     , jtswo3(nd_profile, nd_layer)
!        Index of sw o3 reference temp

  LOGICAL ::                                                            &
       l_continuum
!        water continuum absorption if set true

! Local variables
  REAL (RealK), PARAMETER :: oneminus = .999999_RealK
!        Uplimit of binary parameter. If interpolated binary
!        parameter >= 1.0 then ONEMINUS is used in stead

  REAL (RealK) ::                                                       &
       fac000(nd_profile, nd_layer)                                     &
     , fac100(nd_profile, nd_layer)                                     &
     , fac010(nd_profile, nd_layer)                                     &
     , fac110(nd_profile, nd_layer)                                     &
     , fac001(nd_profile, nd_layer)                                     &
     , fac101(nd_profile, nd_layer)                                     &
     , fac011(nd_profile, nd_layer)                                     &
     , fac111(nd_profile, nd_layer)                                     &
!        Fractional factor for interpolation
     , fu, compfu                                                       &
!        fractional factor for binary parameter interpolation
     , specparm, specmult

  INTEGER ::                                                            &
       i, ig, k, l                                                   &
!        Loop index
     , ju(nd_profile, nd_layer)                                         &
!        Index of amount at level JP and JP1 such that
     , jtt1, jto2c1, jt1, ju1, ic                                       &
     , ib
!        Band index for mixed gases

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='INTER_K'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


  IF (l_continuum .AND. ( n_continuum  /=  0 ) ) THEN

    DO ic = 1, n_continuum
      IF (index_continuum(i_band, ic)  ==  ip_cont_h2o ) THEN
        DO k=1, i_band_esft
          DO i=1, n_layer
            DO l=1, n_profile
              jp1=jph2oc(l,i)+1
              jt1=jtswo3(l,i)+1
              k_contm_layer(l,k,i,ic) = MAX(0.0_RealK,              &
                 fac00c(l,i)*k_h2oc(jph2oc(l,i),jtswo3(l,i),k)      &
                +fac10c(l,i)*k_h2oc(jp1,jtswo3(l,i),k)              &
                +fac01c(l,i)*k_h2oc(jph2oc(l,i),jt1,k)              &
                +fac11c(l,i)*k_h2oc(jp1,jt1,k)                      &
                )
            END DO
          END DO
        END DO
      ELSE IF(index_continuum(i_band, ic)  ==  ip_cont_o2n2_2) THEN

!       These coefficients have temperature dependency.
!       o2n2_2 for 1.27 nu
!vdir noloopchg
        DO k=1, i_band_esft
          DO i=1, n_layer
            DO l=1, n_profile
              jto2c1=jto2c(l,i)+1
              k_contm_layer(l,k,i,ic) = MAX(0.0_RealK,              &
                 facc00(l,i)*k_continuum(k,jto2c(l,i),i_band,ic)    &
                +facc01(l,i)*k_continuum(k,jto2c1,    i_band,ic)    &
                )
            END DO
          END DO
        END DO
      ELSE

!       O2 continuum: no temperature dependency
        DO k=1, i_band_esft
          DO i=1, n_layer
            DO l=1, n_profile
              k_contm_layer(l,k,i,ic)=k_continuum(k,1,i_band,ic)
            END DO
          END DO
        END DO
      END IF
    END DO
  END IF


  DO ig=1, n_band_absorb(i_band)
    DO k=1, i_band_esft
      DO i=1, n_layer
        DO l=1, n_profile
          jp1=jp(l,i)+1
          jt1=jt(l,i)+1
          jtt1=jtt(l,i)+1
          k_esft_layer(l,i,k,ig) = MAX(0.0_RealK,                   &
             fac00(l,i)*k_esft(jt(l,i),jp(l,i),k,ig)                &
            +fac10(l,i)*k_esft(jtt(l,i),jp1,k,ig)                   &
            +fac01(l,i)*k_esft(jt1,jp(l,i),k,ig)                    &
            +fac11(l,i)*k_esft(jtt1,jp1,k,ig)                       &
            )
        END DO
      END DO
    END DO
  END DO


  IF ( n_mix_gas  /=  0) THEN

!   Mixed species in this band

    ib=mix_gas_band(i_band)
    DO i=1, n_layer
      DO l=1, n_profile
        gas_mix_amt(l,i)                                            &
          =gas_mix_ratio(l, i, index_mix_gas(1,ib))                 &
                                                  ! key species
          +f_mix*gas_mix_ratio(l, i, index_mix_gas(2,ib))
                                                  ! minor species


        specparm=                                                   &
          MIN(gas_mix_ratio(l, i, index_mix_gas(1,ib))              &
            /gas_mix_amt(l,i), oneminus)

        specmult=8.0_RealK*specparm
        ju(l,i)=1 + INT(specmult)
        fu=1.0_RealK-MOD(specmult, 1.0_RealK)
        compfu=1. -fu
        fac000(l,i)=fac00(l,i)*fu
        fac100(l,i)=fac10(l,i)*fu
        fac010(l,i)=fac01(l,i)*fu
        fac110(l,i)=fac11(l,i)*fu
        fac001(l,i)=fac00(l,i)*compfu
        fac101(l,i)=fac10(l,i)*compfu
        fac011(l,i)=fac01(l,i)*compfu
        fac111(l,i)=fac11(l,i)*compfu
      END DO
    END DO


!cdir noloopchg
    DO k=1, i_band_esft
      DO i=1, n_layer
        DO l=1, n_profile
          jp1=jp(l,i)+1
          jt1=jt(l,i)+1
          jtt1=jtt(l,i)+1
          ju1=ju(l,i)+1
          k_mix_gas_layer(l,k,i)                                    &
            =MAX(0.0_RealK,                                         &
            fac000(l,i)*k_mix_gas(jp(l,i),jt(l,i),ju(l,i),k,ib)     &
           +fac100(l,i)*k_mix_gas(jp1,jtt(l,i),ju(l,i),k,ib)        &
           +fac010(l,i)*k_mix_gas(jp(l,i),jt1,ju(l,i),k,ib)         &
           +fac110(l,i)*k_mix_gas(jp1,jtt1,ju(l,i),k,ib)            &
           +fac001(l,i)*k_mix_gas(jp(l,i),jt(l,i),ju1,k,ib)         &
           +fac101(l,i)*k_mix_gas(jp1,jtt(l,i),ju1,k,ib)            &
           +fac011(l,i)*k_mix_gas(jp(l,i),jt1,ju1,k,ib)             &
           +fac111(l,i)*k_mix_gas(jp1,jtt1,ju1,k,ib)                &
            )
        END DO
      END DO
    END DO

  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE inter_k
END MODULE inter_k_mod
