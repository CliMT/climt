! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set moist aerosol properties independent of bands.
!
! Method:
!   The clear-sky mean relative humidities are calculated and
!   pointers to the lookup tables are set.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------
SUBROUTINE set_moist_aerosol_properties(ierr                            &
    , n_profile, n_layer                                                &
    , n_aerosol, i_aerosol_parametrization, nhumidity                   &
    , water_mix_ratio, t, p, w_cloud, delta_humidity                    &
    , mean_rel_humidity, i_humidity_pointer                             &
    , nd_profile, nd_layer, id_ct, nd_aerosol_species                   &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer                                                          &
!       Maximum number of layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_aerosol_species
!       Maximum number of aerosols

! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag

  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_aerosol                                                         &
!       Number of aerosol species
    , i_aerosol_parametrization(nd_aerosol_species)                     &
!       Parametrizations of aerosol species
    , nhumidity(nd_aerosol_species)
!       Number of humidity values
  INTEGER, INTENT(OUT) ::                                               &
      i_humidity_pointer(nd_profile, nd_layer)
!       Pointers to look-up tables
  REAL (RealK), INTENT(IN) ::                                           &
      water_mix_ratio(nd_profile, nd_layer)                             &
!       Mixing ratio of water vapour
    , t(nd_profile, nd_layer)                                           &
!       Temperatures
    , p(nd_profile, nd_layer)                                           &
!       Pressures
    , w_cloud(nd_profile, id_ct: nd_layer)
!       Cloud fraction
  REAL (RealK), INTENT(OUT) ::                                          &
      mean_rel_humidity(nd_profile, nd_layer)                           &
!       Mean humidities of layers
    , delta_humidity
!       Increment in humidity


! local variables.
  INTEGER ::                                                            &
      i                                                                 &
!       Loop variable
    , j                                                                 &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , nhumidity_common
!       Common number of humidities for moist aerosols
  REAL (RealK) ::                                                       &
      mix_ratio_sat(nd_profile, nd_layer)                               &
!       Saturated humidity mixing ratio
    , mix_ratio_clr(nd_profile, nd_layer)
!       Humidity mixing ratio of the clear-sky part of the grid-box

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=256)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'SET_MOIST_AEROSOL_PROPERTIES'


  IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Set up array of pointers to `include' the effects of humidity.
! Calculate the saturated mixing ratio.
  DO i=1, n_layer
! DEPENDS ON: qsat_wat
    CALL qsat_wat(mix_ratio_sat(1, i), t(1, i), p(1, i)                 &
      , n_profile)
  END DO
! Calculate the clear-sky mixing ratio
  mix_ratio_clr(1:n_profile, 1:id_ct-1)                                 &
     = water_mix_ratio(1:n_profile, 1:id_ct-1)
  DO i=id_ct, n_layer
    DO l = 1, n_profile
      mix_ratio_clr(l, i) =                                             &
        (water_mix_ratio(l, i) - (w_cloud(l, i) *                       &
         mix_ratio_sat(l, i))) /                                        &
        MAX( (1.00e+00_RealK - w_cloud(l, i)), 0.001_RealK)
      ! MAKE SURE THAT 0 <= MIX_RATIO_CLR <= WATER_MIX_RATIO
      mix_ratio_clr(l, i) = MIN(mix_ratio_clr(l, i),                    &
                                water_mix_ratio(l, i))
      mix_ratio_clr(l, i) = MAX(mix_ratio_clr(l, i),                    &
                                0.00e+00_RealK)
    END DO
  END DO

! Determine the number of humidities to be used for moist
! aerosols. This must be the same for all moist aerosols
! in the current version of the code.
  nhumidity_common=0
  DO j=1, n_aerosol
    IF (i_aerosol_parametrization(j) == ip_aerosol_param_moist)         &
         THEN
      IF (nhumidity(j) >  0) THEN
!       Set the actual common value.
        IF (nhumidity_common == 0) THEN
          nhumidity_common=nhumidity(j)
        ELSE IF (nhumidity(j) /= nhumidity_common) THEN
!         There is an inconsistency.
          cmessage = '***Error: The look-up tables for moist aerosols ' &
            //'are of different sizes. This is not permitted.'
          ierr=i_err_fatal
          CALL ereport(RoutineName, ierr, cmessage)
        END IF
      END IF
    END IF
  END DO
! The look-up table is assumed to be uniform in humidity.
  delta_humidity=1.0e+00_RealK                                          &
    /(REAL(nhumidity_common, RealK)-1.0e+00_RealK)

! Compute the clear-sky mean relative humidity
  DO i = 1, n_layer
    DO l = 1, n_profile
      mean_rel_humidity(l,i)                                            &
         = mix_ratio_clr(l, i)*(1.0e+00_RealK-mix_ratio_sat(l, i))      &
         /((1.0e+00_RealK-mix_ratio_clr(l, i))*mix_ratio_sat(l, i))
!     Check that the clear-sky mean relative humidity
!     lies between 0.0 (incl.) and 1.0 (excl.)
      mean_rel_humidity(l, i)=MIN(mean_rel_humidity(l, i)               &
        , 0.99999_RealK)
      mean_rel_humidity(l, i)=MAX(mean_rel_humidity(l, i)               &
        , 0.00000_RealK)
      i_humidity_pointer(l, i)=1                                        &
         +INT(mean_rel_humidity(l, i)*(nhumidity_common-1))
    END DO
  END DO


  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_moist_aerosol_properties
