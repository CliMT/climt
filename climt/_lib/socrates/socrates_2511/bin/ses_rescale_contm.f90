! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Perform path-length scaling for the different continuum absorbers.
!
!- ---------------------------------------------------------------------
MODULE ses_rescale_contm_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SES_RESCALE_CONTM_MOD'
CONTAINS
SUBROUTINE ses_rescale_contm(nd_profile, nd_layer                       &
     , i_continuum, n_profile, n_layer                                  &
     , p, t, gas_mix_ratio, amount_continuum                            &
     )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_cont_h2o, ip_cont_o4
  USE rad_ccf, ONLY: mol_weight_air, n_avogadro, k_boltzmann
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


  INTEGER, INTENT(IN) ::                                                &
       nd_profile                                                       &
!         Max number of profile
     , nd_layer                                                         &
!         Max number of layer
     , n_profile                                                        &
!         Number of profiles
     , n_layer                                                          &
!         Number of layers
     , i_continuum
!         gas id index of continuum
  REAL (RealK), INTENT(IN) ::                                           &
       gas_mix_ratio(nd_profile, nd_layer)                              &
!         Gaseous mass mixing ratios
     , t(nd_profile, nd_layer)                                          &
!         Temperature
     , p(nd_profile, nd_layer)
!         Pressure
  REAL (RealK), INTENT(OUT) ::                                          &
       amount_continuum(nd_profile, nd_layer)
!         Amount of continuum

! Local variables.
  INTEGER ::                                                            &
       l                                                                &
!         Loop variable
     , i
!         Loop variable
   REAL (RealK) ::                                                      &
       o2mix, n2mix                                                     &
!         O2 & N2 mixing ratio
     , o2n2, o2o2                                                       &
!         produce of mixing ratio and avog number
     , k_m_air
!         K * air molecular weight

   PARAMETER(o2mix=0.2315e-23, n2mix=0.7552e-23                         &
     , o2n2=1.e-10*n_avogadro*o2mix*n2mix                               &
     , o2o2=1.e-10*n_avogadro*o2mix*o2mix)

!  Note: E-23 in o2mix n2mix is for unit conversion
!        E-10 for o2n2 o2o2 also for unit conversion

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SES_RESCALE_CONTM'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  k_m_air=k_boltzmann*mol_weight_air
  IF (i_continuum  ==  ip_cont_h2o ) THEN
    DO i=1, n_layer
      DO l=1, n_profile
        amount_continuum(l, i) = gas_mix_ratio(l, i)
      END DO
    END DO
  ELSE IF (i_continuum  ==  ip_cont_o4  ) THEN

! (i_continuum=2 --> O4)
! tau = k * p * o2mix * o2mix * N / (K*T*m_air) * d_mass

    DO i=1, n_layer
      DO l=1, n_profile
        amount_continuum(l, i)=p(l,i)*o2o2/(t(l,i)*k_m_air)
      END DO
    END DO
  ELSE

! o2n2 continuum at 1.27 & 1.06 microns
! tau = k * p * o2mix * n2mix * N / (K*T*m_air) * d_mass
! convert k to m5 molec-2 requires a factor 1.e-10

    DO i=1, n_layer
      DO l=1, n_profile
        amount_continuum(l, i)=p(l,i) * o2n2/(t(l,i)*k_m_air)
      END DO
    END DO
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE ses_rescale_contm
END MODULE ses_rescale_contm_mod
