! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to scale amounts of absorbers.
!
! Method:
!   A scaling factor for the gas mixing ratio is determined
!   by the type of scaling selected.
!
!- ---------------------------------------------------------------------
MODULE scale_absorb_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SCALE_ABSORB_MOD'
CONTAINS
SUBROUTINE scale_absorb(ierr, n_profile, n_layer                        &
    , p, t                                                              &
    , gas_frac_rescaled                                                 &
    , i_fnc, p_reference, t_reference, scale_parameter                  &
    , iex, i_band                                                       &
    , l_doppler, doppler_correction                                     &
    , nd_profile, nd_layer                                              &
    , nd_scale_variable                                                 &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_fatal, ip_scale_power_law,                   &
                     ip_scale_power_quad, ip_scale_doppler_quad,        &
                     ip_scale_dbl_pow_law, ip_scale_fnc_null,           &
                     ip_scale_dbl_pow_quad, ip_scale_wenyi
  USE vectlib_mod, ONLY : rtor_v
  USE scale_wenyi, ONLY: plg, ttb, tto, gk250b, gk4, gk6
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_layer                                                          &
!       Size allocated for layers
    , nd_scale_variable
!       Size allocated for of scaling variables

! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , i_fnc                                                             &
!       Type of scaling function
    , iex                                                               &
!       Index of ESFT term
    , i_band
!       Band being considered
  LOGICAL, INTENT(IN) ::                                                &
      l_doppler
!       Flag for Doppler term
  REAL (RealK), INTENT(IN) ::                                           &
      p(nd_profile, nd_layer)                                           &
!       Pressure
    , t(nd_profile, nd_layer)                                           &
!       Temperature
    , p_reference                                                       &
!       Reference pressure
    , t_reference                                                       &
!       Reference temperature
    , scale_parameter(nd_scale_variable)                                &
!       Scaling paramters
    , doppler_correction
!       Doppler-broadening correction
  REAL (RealK), INTENT(OUT) ::                                          &
      gas_frac_rescaled(nd_profile, nd_layer)
!       Rescaling factor for mass fraction of gas

! Local variables.
  INTEGER ::                                                            &
      l, i, jp, jp1, jt, jt1
!       Loop variables
  REAL (RealK) ::                                                       &
      pressure_offset
!       Offset to pressure

  REAL (RealK) :: cgp, gkpb, gkpc
  REAL (RealK) :: pwk_in(n_profile,n_layer) ! Workspace
  REAL (RealK) :: pwk(n_profile,n_layer)    ! Workspace
  REAL (RealK) :: twk_in(n_profile,n_layer) ! Workspace
  REAL (RealK) :: twk(n_profile,n_layer)    ! Workspace
  REAL (RealK) :: sp1(n_profile,n_layer)    ! Workspace
  REAL (RealK) :: sp2(n_profile,n_layer)    ! Workspace
  INTEGER :: n_input      ! No. of inputs for rtor_v function
  REAL (RealK) :: tmp, t_inv, p_ref_off_inv

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'SCALE_ABSORB'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the offset to the pressure for the Doppler correction.
  IF (l_doppler) THEN
    pressure_offset=doppler_correction
  ELSE
    pressure_offset=0.0e+00_RealK
  END IF

  IF ((i_fnc == ip_scale_power_law)  .OR.                               &
      (i_fnc == ip_scale_power_quad) .OR.                               &
      (i_fnc == ip_scale_doppler_quad)) THEN
    DO i=1, n_layer
      DO l=1, n_profile
        sp1(l,i)=scale_parameter(1)
        sp2(l,i)=scale_parameter(2)
      END DO
    END DO
    n_input=(n_layer)*n_profile
  END IF

! The array gas_frac_rescaled is used initially to hold only the
! scaling functions, and only later is it multiplied by the
! mixing ratios
  IF (i_fnc == ip_scale_power_law) THEN

    t_inv = 1.0_RealK/t_reference
    p_ref_off_inv = 1.0_RealK/(p_reference+pressure_offset)
    DO i=1, n_layer
      DO l=1, n_profile
        pwk_in(l,i)=(p(l,i)+pressure_offset)*p_ref_off_inv
        twk_in(l,i)=t(l,i)*t_inv
      END DO
    END DO
    CALL rtor_v(n_input,pwk_in,sp1,pwk)
    CALL rtor_v(n_input,twk_in,sp2,twk)
    DO i=1, n_layer
      DO l=1, n_profile
        gas_frac_rescaled(l, i)=pwk(l,i)*twk(l,i)
      END DO
    END DO

  ELSE IF (i_fnc == ip_scale_dbl_pow_law) THEN

    DO i=1, n_layer
      DO l=1, n_profile
        IF (p(l, i) > scale_parameter(5)) THEN
          gas_frac_rescaled(l, i)=                                      &
            EXP( scale_parameter(3)*LOG( p(l, i)/scale_parameter(5) )   &
               + scale_parameter(4)*LOG( t(l, i)/scale_parameter(6) ) )
        ELSE
          gas_frac_rescaled(l, i)=                                      &
            EXP( scale_parameter(1)*LOG( p(l, i)/scale_parameter(5) )   &
               + scale_parameter(2)*LOG( t(l, i)/scale_parameter(6) ) )
        END IF
      END DO
    END DO

  ELSE IF (i_fnc == ip_scale_fnc_null) THEN

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN

  ELSE IF (i_fnc == ip_scale_power_quad) THEN

    p_ref_off_inv = 1.0_RealK/(p_reference+pressure_offset)
    DO i=  1, n_layer
      DO l=1, n_profile
        pwk_in(l,i)=(p(l,i)+pressure_offset)*p_ref_off_inv
      END DO
    END DO
    CALL rtor_v(n_input,pwk_in,sp1,pwk)
    t_inv = 1.0_RealK/t_reference
    DO i=1, n_layer
      DO l=1, n_profile
        tmp = t(l,i)*t_inv - 1.0_RealK
        gas_frac_rescaled(l, i)=pwk(l,i) &
          *(1.0e+00_RealK+tmp*scale_parameter(2) &
          +scale_parameter(3)*tmp*tmp)
      END DO
    END DO

  ELSE IF (i_fnc == ip_scale_dbl_pow_quad) THEN

    DO i=1, n_layer
      DO l=1, n_profile
        IF (p(l, i) > scale_parameter(7)) THEN
          gas_frac_rescaled(l, i)=                                      &
            EXP( scale_parameter(4)*LOG( p(l, i)/scale_parameter(7) ) ) &
            *( 1.0e+00_RealK + scale_parameter(5)                       &
            *( t(l, i)/scale_parameter(8)-1.0e+00_RealK )               &
            +scale_parameter(6)                                         &
            *( t(l, i)/scale_parameter(8)-1.0e+00_RealK )**2 )
        ELSE
          gas_frac_rescaled(l, i)=                                      &
            EXP( scale_parameter(1)*LOG( p(l, i)/scale_parameter(7) ) ) &
            *( 1.0e+00_RealK + scale_parameter(2)                       &
            *( t(l, i)/scale_parameter(8)-1.0e+00_RealK )               &
            +scale_parameter(3)                                         &
            *( t(l, i)/scale_parameter(8)-1.0e+00_RealK )**2 )
        END IF
      END DO
    END DO

  ELSE IF (i_fnc == ip_scale_doppler_quad) THEN

!   There is no Doppler term here since it is implicitly included
!   in the scaling.
    DO i=  1, n_layer
      DO l=1, n_profile
        pwk_in(l,i)=(p(l,i)+scale_parameter(2)) &
                   /(p_reference+scale_parameter(2))
      END DO
    END DO
    CALL rtor_v(n_input,pwk_in,sp1,pwk)
    t_inv = 1.0_RealK/t_reference
    DO i=1, n_layer
      DO l=1, n_profile
        tmp = t(l,i)*t_inv - 1.0_RealK
        gas_frac_rescaled(l, i)=pwk(l,i) &
          *(1.0e+00_RealK+tmp*scale_parameter(3) &
          +scale_parameter(4)*tmp*tmp)
      END DO
    END DO

  ELSE IF (i_fnc == ip_scale_wenyi) THEN

    IF (i_band  ==  4) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          cgp  = MAX(-5.5_RealK, LOG(p(l, i)/100.0_RealK))
          jp   = INT((5.5+cgp)*2.)+1
          jp   = MIN(MAX(jp,1),26)
          jp1  = MIN(jp+1,26)
          jt   = INT((t(l,i)-240.0)/20.+7.)
          jt   = MIN(MAX(jt,1),10)
          jt1  = MIN(jt+1,10)
          gkpb = gk4(jt,jp,iex)+(gk4(jt,jp1,iex)-                       &
             gk4(jt,jp,iex))*(cgp-plg(jp))*2.
          gkpc = gk4(jt1,jp,iex)+(gk4(jt1,jp1,iex)-                     &
             gk4(jt1,jp,iex))*(cgp-plg(jp))*2.
          gas_frac_rescaled(l, i) = (gkpb+(gkpc-gkpb)*                  &
                     (t(l,i)-ttb(jt))/20.0)/gk250b(iex)
        END DO
      END DO
    ELSE IF(i_band == 6)THEN
      DO i=1, n_layer
        DO l=1, n_profile
          cgp  = MAX(-5.5_RealK, LOG(p(l, i)/100.0_RealK))
          jp   = INT((5.5+cgp)*2.)+1
          jp   = MIN(MAX(jp,1),26)
          jp1  = MIN(jp+1,26)
          jt   = INT((t(l,i)-240.0)/20.+5.)
          jt   = MIN(MAX(jt, 1),9)
          jt1  = MIN(jt+1,9)
          gkpb = gk6(jt,jp,iex)+(gk6(jt,jp1,iex)-                       &
                  gk6(jt,jp,iex))*(cgp-plg(jp))*2.
          gkpc = gk6(jt1,jp,iex)+(gk6(jt1,jp1,iex)-                     &
                  gk6(jt1,jp,iex))*(cgp-plg(jp))*2.
          gas_frac_rescaled(l, i) = (gkpb+(gkpc-gkpb)*                  &
                               (t(l,i)-tto(jt))/20.0)
        END DO
      END DO
    END IF

  ELSE
    cmessage = '*** Error: an illegal type of scaling has been given.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE scale_absorb
END MODULE scale_absorb_mod
