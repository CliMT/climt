! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate the singly scattered solar radiance.
!
! Purpose:
!   This subroutine is used to increment the radiances in the
!   required directions on the viewing levels with the singly
!   scattered solar radiance.
!
! Method:
!   Each direction is considered in turn. For each layer of the
!   atmosphere an angular factor involving the phase function
!   and for each viewing level a geometric factor involving the
!   optical depth between the layer in question and the viewing
!   level is calculated. The product of these with the solar
!   beam gives the contribution of that layer to the increment
!   to the radiance.
!
!- ---------------------------------------------------------------------
MODULE single_scat_sol_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SINGLE_SCAT_SOL_MOD'
CONTAINS
SUBROUTINE single_scat_sol(n_profile, n_layer                           &
    , n_direction, direction                                            &
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
    , i_direct, mu_0                                                    &
    , tau, omega, phase_fnc_solar                                       &
    , radiance                                                          &
    , nd_profile, nd_radiance_profile                                   &
    , nd_layer, nd_direction, nd_viewing_level                          &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Dummy arguments:
! Sizes of arrays:
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_radiance_profile                                               &
!       Size allocated for atmospheric profiles where radiances
!       are calculated
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_viewing_level                                                  &
!       Size allocated for viewing levels
    , nd_direction
!       Size allocated for viewing directions

! Atmospheric structure:
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of atmospheric profiles
    , n_layer
!       Number of atmospheric layers

! Viewing geometry:
  INTEGER, INTENT(IN) ::                                                &
      n_direction                                                       &
!       Number of directions
    , n_viewing_level                                                   &
!       Number of levels where the radiance is calculated
    , i_rad_layer(nd_viewing_level)
!       Indices of layers containing viewing levels
  REAL (RealK), INTENT(IN) ::                                           &
      direction(nd_radiance_profile, nd_direction, 2)                   &
!       Cosines of polar viewing angles
    , frac_rad_layer(nd_viewing_level)
!       Fraction optical depth into its layer of the viewing level

! Solar Radiances
  REAL (RealK), INTENT(IN) ::                                           &
      i_direct(nd_profile, 0: nd_layer)                                 &
!       Direct solar radiances
    , mu_0(nd_profile)
!       Cosines of solar zenith angles

! Optical properties of the atmosphere
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile, nd_layer)                                         &
!       Optical depths
    , omega(nd_profile, nd_layer)                                       &
!       Albedos of single scattering
    , phase_fnc_solar(nd_radiance_profile, nd_layer, nd_direction)
!       Phase function


  REAL (RealK), INTENT(INOUT) ::                                        &
      radiance(nd_radiance_profile, nd_viewing_level, nd_direction)
!       Radiances


! Local variables
  INTEGER                                                               &
      l                                                                 &
!       Loop variable (points)
    , ll                                                                &
!       Loop variable
    , i                                                                 &
!       Loop variable
    , ii                                                                &
!       Loop variable
    , ir                                                                &
!       Loop variable (radiative levels)
    , id
!       Loop variable (directions)
  INTEGER                                                               &
      n_list_up                                                         &
!       Numbers of points where the current viewing direction
!       points up
    , list_up(nd_profile)                                               &
!       List up points with where the current viewing direction
!       points up
    , n_list_down                                                       &
!       Numbers of points where the current viewing direction
!       points down
    , list_down(nd_profile)
!       List up points with where the current viewing direction
!       points up
  REAL (RealK) ::                                                       &
      geom(nd_profile)                                                  &
!       Geometrical factor
    , m_slant_depth_near(nd_profile)                                    &
!       Minus slantwise optical distance between the radiance
!       level and the nearer boundary of the current layer
    , m_slant_depth_far(nd_profile)                                     &
!       Minus slantwise optical distance between the radiance
!       level and the farther boundary of the current layer
    , tau_i(nd_profile)                                                 &
!       Optical depth of the relevant part of the current layer
    , trans_d(nd_profile)                                               &
!       Direct transmission from the layer containing the viewing
!       level to the viewung level
    , d_mu
!       Difference in cosines of directions

! Variables related to the treatment of ill-conditioning
  REAL (RealK) ::                                                       &
      eps_r                                                             &
!       The smallest real number such that 1.0-EPS_R is not 1
!       to the computer's precision
    , sq_eps_r                                                          &
!       The square root of the above
    , eta                                                               &
!       The conditioning weight
    , eta_nm
!       The conditioning multiplier applied in the numerator

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SINGLE_SCAT_SOL'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the tolerances used in avoiding ill-conditioning, testing
! on any variable.
  eps_r=EPSILON(mu_0(1))
  sq_eps_r=SQRT(eps_r)

! Consider each direction in turn. Collect points where the
! viewing direction is upward and points where it is downward,
! then calculate the angular and geometric factors.
  DO id=1, n_direction

!   Collect points where the viewing direction is upward:
!   (Horizontal directions are not allowed).
    n_list_up=0
    DO l=1, n_profile
      IF (direction(l, id, 1) >  0.0e+00_RealK) THEN
        n_list_up=n_list_up+1
        list_up(n_list_up)=l
      END IF
    END DO

!   Collect points where the viewing direction is downward:
    n_list_down=0
    DO l=1, n_profile
      IF (direction(l, id, 1) <  0.0e+00_RealK) THEN
        n_list_down=n_list_down+1
        list_down(n_list_down)=l
      END IF
    END DO

!   Go through each atmospheric layer calculating the radiance
!   at each observing level.
    DO i=1, n_layer

!     Calculate the geometric factors:
      DO ir=1, n_viewing_level

!       Upward Radiances:
!       Contributions arise only from layers below the viewing
!       level.
        IF (i >= i_rad_layer(ir)) THEN

!         Calculate minus the slantwise optical depths to the
!         boundaries of the layer. If the level where the radiance
!         is required lies in the current layer we perform the
!         calculation for a temporary layer reaching from the
!         viewing level to the bottom of the current layer.
          IF (i >  i_rad_layer(ir)) THEN
!           Full layers are required.
            DO ll=1, n_list_up
              l=list_up(ll)
              m_slant_depth_near(l)                                     &
                =(1.0e+00_RealK-frac_rad_layer(ir))                     &
                *tau(l, i_rad_layer(ir))
            END DO
            DO ii=i_rad_layer(ir)+1, i-1
              DO ll=1, n_list_up
                l=list_up(ll)
                m_slant_depth_near(l)                                   &
                 =m_slant_depth_near(l)+tau(l, ii)
              END DO
            END DO
            DO ll=1, n_list_up
              l=list_up(ll)
              m_slant_depth_near(l)                                     &
                =-m_slant_depth_near(l)/direction(l, id, 1)
              m_slant_depth_far(l)=m_slant_depth_near(l)                &
                -tau(l, i)/direction(l, id, 1)
!             Collect the local optical depth to allow the use of
!             generic code later.
              tau_i(l)=tau(l, i)
              trans_d(l)=1.0e+00_RealK
            END DO
          ELSE IF (i == i_rad_layer(ir)) THEN
!           The viewing level lies in the current layer.
            DO ll=1, n_list_up
              l=list_up(ll)
              m_slant_depth_near(l)=0.0e+00_RealK
              m_slant_depth_far(l)                                      &
                =-(1.0e+00_RealK-frac_rad_layer(ir))*tau(l, i)          &
                /direction(l, id, 1)
              tau_i(l)=(1.0e+00_RealK-frac_rad_layer(ir))*tau(l, i)
              trans_d(l)                                                &
                =EXP(-frac_rad_layer(ir)*tau(l, i)/mu_0(l))
            END DO
          END IF


!         Set the geometrical term and increment the radiance.
          DO ll=1, n_list_up
            l=list_up(ll)
            geom(l)=(mu_0(l)/(mu_0(l)+direction(l, id, 1)))             &
              *(EXP(m_slant_depth_near(l))                              &
              -EXP(m_slant_depth_far(l)-tau_i(l)/mu_0(l)))
            radiance(l, ir, id)=radiance(l, ir, id)                     &
              +i_direct(l, i-1)*trans_d(l)*geom(l)*(omega(l, i)         &
              /(4.0e+00_RealK*pi))*phase_fnc_solar(l, i, id)
          END DO

        END IF


!       Downward Radiances:
!       Contributions arise only from layers above the viewing
!       level.
        IF (i <= i_rad_layer(ir)) THEN

!         Calculate the slantwise optical depths to the
!         boundaries of the layer. If the observing level lies
!         within the current layer we perform the calculation for
!         a layer reaching from the top of the current layer to
!         the observing level.
          IF (i <  i_rad_layer(ir)) THEN
            DO ll=1, n_list_down
              l=list_down(ll)
              m_slant_depth_near(l)                                     &
              =frac_rad_layer(ir)*tau(l, i_rad_layer(ir))
            END DO
            DO ii=i_rad_layer(ir)-1, i+1, -1
              DO ll=1, n_list_down
                l=list_down(ll)
                m_slant_depth_near(l)                                   &
                  =m_slant_depth_near(l)+tau(l, ii)
              END DO
            END DO
            DO ll=1, n_list_down
              l=list_down(ll)
              m_slant_depth_near(l)                                     &
                =m_slant_depth_near(l)/direction(l, id, 1)
              m_slant_depth_far(l)=m_slant_depth_near(l)                &
                +tau(l, i)/direction(l, id, 1)
              tau_i(l)=tau(l, i)
            END DO
          ELSE
!           The viewing level lies in the current layer.
            DO ll=1, n_list_down
              l=list_down(ll)
              tau_i(l)=frac_rad_layer(ir)*tau(l, i)
              m_slant_depth_near(l)=0.0e+00_RealK
              m_slant_depth_far(l)=tau_i(l)/direction(l, id, 1)
            END DO
          END IF


!         Set the geometrical terms for the solar integral.
          DO ll=1, n_list_down
            l=list_down(ll)
!           This may exhibit ill-conditioning, so it is perturbed
!           using L'Hopital's rule.
            d_mu=mu_0(l)+direction(l, id, 1)
            eta=eps_r/(d_mu+SIGN(sq_eps_r, d_mu))
            eta_nm=(1.0e+00_RealK-eta*tau_i(l)                          &
              /(mu_0(l)*direction(l, id, 1)))
            geom(l)=(mu_0(l)/(d_mu+eta))                                &
              *(EXP(m_slant_depth_near(l)-tau_i(l)/mu_0(l))             &
              *eta_nm                                                   &
              -EXP(m_slant_depth_far(l)))
            radiance(l, ir, id)=radiance(l, ir, id)                     &
              +i_direct(l, i-1)*geom(l)                                 &
              *(omega(l, i)/(4.0e+00_RealK*pi))                         &
              *phase_fnc_solar(l, i, id)
          END DO

        END IF

      END DO
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE single_scat_sol
END MODULE single_scat_sol_mod
