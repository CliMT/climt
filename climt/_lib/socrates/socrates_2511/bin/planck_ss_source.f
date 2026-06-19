! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calcaluate Planckian fluxes or radiances.
!
! Purpose:
!   To calculate the Planckian flux or radiance integrated across
!   a given spectral region.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE planck_ss_source(ierr, n_profile, n_layer
     &  , i_angular_integration
     &  , t_level, l_ir_source_quad, t, t_ground
     &  , wavelength_short, wavelength_long
     &  , n_viewing_level, i_rad_layer, frac_rad_layer
     &  , diff_planck, diff_planck_2, d_planck_flux_surface
     &  , planck_flux, planck_radiance
     &  , nd_profile, nd_layer
     &  , nd_radiance_profile, nd_viewing_level
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE rad_pcf
      USE rad_ccf, ONLY: pi, h_planck, c_light, k_boltzmann
!
!
      IMPLICIT NONE
!
!
!
!
!
!     Declaration of variables.
!
!     Sizes of arrays:
      INTEGER, Intent(IN) ::
     &    nd_profile
!           Size allocated for atmospheric profiles
     &  , nd_layer
!           Size allocated for atmospheric layers
     &  , nd_radiance_profile
!           Number of profiles in arrays of radiances
     &  , nd_viewing_level
!           Number of layers in arrays of radiances
!
!     Control of errors:
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
!
!     Algorithmic control:
      INTEGER, Intent(IN) ::
     &    i_angular_integration
!           Type of angular integration
      LOGICAL, Intent(IN) ::
     &    l_ir_source_quad
!           Flag to use a quadratic source function in the IR
!
!
!
!     General atmospheric properties:
      INTEGER, Intent(IN) ::
     &    n_profile
!           Number of profiles
     &  , n_layer
!           Number of layers
      REAL  (RealK), Intent(IN) ::
     &    t(nd_profile, nd_layer)
!           Temperatures at the mid-points of layers
     &  , t_level(nd_profile, 0: nd_layer)
!           Temperatures at the edges of layers
     &  , t_ground(nd_profile)
!           Temperature of the surface
!
!     Infra-red properties:
      REAL  (RealK), Intent(IN) ::
     &    wavelength_short
!           Shorter wavelength of spectral region
     &  , wavelength_long
!           Longer wavelength of spectral region
!
      INTEGER, Intent(IN) ::
     &    n_viewing_level
!           Number of levels on which radiances are calculated
     &  , i_rad_layer(nd_viewing_level)
!           Layers in which radiances are calculated
      REAL  (RealK), Intent(IN) ::
     &    frac_rad_layer(nd_viewing_level)
!           Fractions below the tops of the layers
!
!     Planckian variables:
      REAL  (RealK), Intent(OUT) ::
     &    diff_planck(nd_profile, nd_layer)
!           Difference in pi*Planckian function
     &  , diff_planck_2(nd_profile, nd_layer)
!           2x2nd differences of Planckian
     &  , d_planck_flux_surface(nd_profile)
!           Differential Planckian flux from the surface
     &  , planck_flux(nd_profile, 0: nd_layer)
!           Total Planckian flux
     &  , planck_radiance(nd_radiance_profile, nd_viewing_level)
!           Total Planckian radiance
!
!
!
!     Local variables:
!
!     Controlling variables
      INTEGER
     &    i
!           Loop variable
     &  , l
!           Loop variable
!
!     Local folded constants:
      REAL  (RealK) ::
     &    cf1
!           Folded constant
     &  , cf2
!           Folded constant
      parameter(
     &    cf1=h_planck*c_light/k_boltzmann
     &  , cf2=2.0_RealK*c_light*k_boltzmann/cf1**3
     &  )
      REAL  (RealK) ::
     &    arg_l
!           Argument to Planckian function at the longer wavelength
     &  , arg_s
!           Argument to Planckian function at the longer wavelength
     &  , tol
!           Tolerance used in calculating the Planckian
     &  , t_temp
!           Temporary temperature
     &  , planck_radiance_top
!           Planckian radiance at the top of a layer
     &  , planck_radiance_bottom
!           Planckian radiance at the bottom of a layer
!
      REAL  (RealK) ::
     &    planck_cumul
!           Function to calculate the cumulative 
!           radiant Planck function
!
!     Functions called:
      EXTERNAL
     &    planck_cumul
!
!
!
!     Crude tolerance. The question of narrow bands has not yet been
!     addressed.
      tol=1.0e-06_RealK
!
!
!     The integral of the Planckian across the region is calculated
!     as the difference between the cumulative Planckians from infinite
!     wavelengths up to this point. This cumulative function is
!     expressed in a form independent of physical constants.
      DO i=0, n_layer
        DO l=1, n_profile
          arg_l=cf1/(t_level(l, i)*wavelength_long)
          arg_s=cf1/(t_level(l, i)*wavelength_short)
          planck_flux(l, i)=pi*cf2*t_level(l, i)**4
     &      *(planck_cumul(arg_s, tol)-planck_cumul(arg_l, tol))
        ENDDO
      ENDDO
!
!     Form the differences.
      DO i=1, n_layer
        DO l=1, n_profile
          diff_planck(l, i)=planck_flux(l, i)-planck_flux(l, i-1)
        ENDDO
      ENDDO
!
!     Second differences are required for quadratic variations in the
!     source function.
      IF (l_ir_source_quad) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            arg_l=cf1/(t(l, i)*wavelength_long)
            arg_s=cf1/(t(l, i)*wavelength_short)
            diff_planck_2(l, i)=pi*cf2*t_level(l, i)**4
     &        *(planck_cumul(arg_s, tol)-planck_cumul(arg_l, tol))
            diff_planck_2(l, i)=2.0_RealK*(planck_flux(l, i)
     &        +planck_flux(l, i-1)-2.0_RealK*diff_planck_2(l, i))
          ENDDO
        ENDDO
      ENDIF
!
!     Surface terms:
      DO l=1, n_profile
        arg_l=cf1/(t_ground(l)*wavelength_long)
        arg_s=cf1/(t_ground(l)*wavelength_short)
        d_planck_flux_surface(l)=pi*cf2*t_ground(l)**4
     &    *(planck_cumul(arg_s, tol)-planck_cumul(arg_l, tol))
     &    -planck_flux(l, n_layer)
      ENDDO
!
!     Calculate radiances at the levels.
      IF (i_angular_integration == IP_spherical_harmonic) THEN
        DO i=1, n_viewing_level
          IF (l_ir_source_quad) THEN
            WRITE(iu_err, '(/a)')
     &        'a quadratic source function is not yet available '
     &        , 'in calculations of radiances.'
            ierr=i_err_fatal
            RETURN
          ELSE
            DO l=1, n_profile
!             We must interpolate the Planckian radiance consistently
!             with the rest of the program: interpolating the
!             temperature will not work.
              t_temp=t_level(l, i_rad_layer(i)-1)
              arg_l=cf1
     &          /(t_temp*wavelength_long)
              arg_s=cf1
     &          /(t_temp*wavelength_short)
              planck_radiance_top=cf2*t_temp**4
     &          *(planck_cumul(arg_s, tol)-planck_cumul(arg_l, tol))
              t_temp=t_level(l, i_rad_layer(i))
              arg_l=cf1
     &          /(t_temp*wavelength_long)
              arg_s=cf1
     &          /(t_temp*wavelength_short)
              planck_radiance_bottom=cf2*t_temp**4
     &          *(planck_cumul(arg_s, tol)-planck_cumul(arg_l, tol))
              planck_radiance(l, i)=planck_radiance_top
     &          +(planck_radiance_bottom-planck_radiance_top)
     &          *frac_rad_layer(i)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!
!
!
      RETURN
      END
