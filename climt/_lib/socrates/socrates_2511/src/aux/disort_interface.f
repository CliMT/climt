! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to interface to DISORT.
!
! Purpose:
!   To transfer between the sets of variables normally used and
!   those required by DISORT.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE disort_interface(ierr
     &  , i_angular_integration, l_rescale, n_order_forward
     &  , i_2stream, i_solver, n_order_gauss
     &  , i_truncation, ls_global_trunc, ls_brdf_trunc
     &  , l_henyey_greenstein_pf
     &  , isolir, zen_0, azim_0, solar_irrad
     &  , t_level, l_ir_source_quad, t, t_ground
     &  , wavelength_short, wavelength_long
     &  , l_iso_inc, iso_inc
     &  , n_profile, n_layer
     &  , tau, omega, phase_fnc, n_order_phase, forward_scatter
     &  , n_direction, direction
     &  , n_viewing_level, viewing_level
     &  , n_brdf_basis_fnc, f_brdf, rho_alb
     &  , flux_direct, flux_down, flux_up, radiance
     &  , nd_profile, nd_layer
     &  , nd_phase_term, nd_sph_coeff
     &  , nd_flux_profile, nd_radiance_profile
     &  , nd_viewing_level, nd_direction
     &  , nd_brdf_basis_fnc, nd_brdf_trunc
     &  , nd_source_coeff
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE rad_pcf
      USE rad_ccf, ONLY: pi
      USE set_rad_layer_mod, ONLY: set_rad_layer
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
!     Treatment of errors:
      INTEGER, Intent(IN) ::
     &    ierr
!           Error flag
!
!     Sizes of arrays:
      INTEGER, Intent(IN) ::
     &    nd_profile
!           Size allocated for atmospheric profiles
     &  , nd_layer
!           Size allocated for atmospheric layers
     &  , nd_phase_term
!           Size allocated for terms in the phase function
     &  , nd_sph_coeff
!           Size of array for spherical coefficients
     &  , nd_flux_profile
!           Number of profiles in arrays of output fluxes
     &  , nd_direction
!           Size allowed for directional radiances
     &  , nd_brdf_basis_fnc
!           Size allowed for BRDF basis functions
     &  , nd_brdf_trunc
!           Size allowed for orders of BRDFs
     &  , nd_radiance_profile
!           Number of profiles in arrays of radiances
     &  , nd_viewing_level
!           Number of layers in arrays of radiances
     &  , nd_source_coeff
!           Size allocated for two-stream source coefficients
!
!     Angular integration
      INTEGER, Intent(IN) ::
     &    i_angular_integration
!           Type of angular integration
     &  , i_2stream
!           Type of two-stream scheme
     &  , i_solver
!           Type of solver
     &  , n_order_gauss
!           Order of Gaussian integration (IR only)
     &  , i_truncation
!           Type of spherical truncation
     &  , ls_global_trunc
!           Overall order of spherical truncation
     &  , ls_brdf_trunc
!           Order of truncation applied to BRDFs
     &  , n_order_forward
!           Order of term used to `define' the 
!           forward scattering fraction
      LOGICAL, Intent(IN) ::
     &    l_rescale
!           Flag for rescaling of the phase function
     &  , l_ir_source_quad
!           Flag to use a quadratic source function in the IR
     &  , l_henyey_greenstein_pf
!           Flag for Henyey-Greenstein phase functions
!
!
!     Optical Properties:
      INTEGER, Intent(IN) ::
     &    n_order_phase
!           Order of terms in the phase function
!     Optical properties will be altered by rescaling and resetting for
!     Henyey-Greenstein phase functions.
      REAL  (RealK), Intent(INOUT) ::
     &    tau(nd_profile, nd_layer)
!           Optical depths
     &  , omega(nd_profile, nd_layer)
!           Albedos of single scattering
     &  , phase_fnc(nd_profile, nd_layer, nd_phase_term)
!           Moments of the phase function
     &  , forward_scatter(nd_profile, nd_layer)
!           Forward-scattering fractions
!
!     Specification of radiances required:
      INTEGER, Intent(IN) ::
     &    n_direction
!           Number of directions at which to calculate radiances
     &  , n_viewing_level
!           Number of levels where radiances are required
      REAL  (RealK), Intent(IN) ::
     &    direction(nd_radiance_profile, nd_direction, 2)
!           Directions in which to calculate radiances
     &  , viewing_level(nd_layer+1)
!           List of required levels: in this driver radiances are
!           calculated at all possible levels in the atmosphere
!
!     General atmospheric properties
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
!     Spectral region
      INTEGER, Intent(IN) ::
     &    isolir
!           Switch for spectral region
!
!     Solar Fields
      REAL  (RealK), Intent(IN) ::
     &    zen_0(nd_profile)
!           Secants or cosines of solar zenith angles
     &  , azim_0(nd_profile)
!           Solar azimuthal angles
     &  , solar_irrad(nd_profile)
!           Solar irradiance at the top of the atmosphere
!
!     Infra-red properties:
      REAL  (RealK), Intent(IN) ::
     &    wavelength_short
!           Shorter wavelength of region (IR only)
     &  , wavelength_long
!           Longer wavelength of region (IR only)
!
!     Isotropic incident radiation:
      LOGICAL, Intent(IN) ::
     &    l_iso_inc
!           Flag for incident isotropic radiation
      REAL  (RealK), Intent(IN) ::
     &    iso_inc(nd_profile)
!           Incident isotropic radiance
!
!     Surface properties:
      INTEGER, Intent(IN) ::
     &    n_brdf_basis_fnc
!           Number of basis functions for BRDFs
      REAL  (RealK), Intent(IN) ::
     &    rho_alb(nd_profile, nd_brdf_basis_fnc)
!           Weights of the basis functions
     &  , f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2
     &      , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)
!           Array of BRDF basis terms
!
!
!     Fluxes and radiances calculated:
      REAL  (RealK), Intent(OUT) ::
     &    flux_direct(nd_flux_profile, 0: nd_layer)
!           Direct flux
     &  , flux_down(nd_flux_profile, 0: nd_layer)
!           Totol downward flux
     &  , flux_up(nd_flux_profile, 0: nd_layer)
!           Totol downward flux
     &  , radiance(nd_radiance_profile, nd_viewing_level, nd_direction)
!           Calculated radiances
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
     &  , ls
!           Loop variable
     &  , j
!           Loop variable
     &  , k
!           Loop variable
     &  , k1
!           Loop variable
      LOGICAL
     &    l_included
!           Flag to test whether the angle is included
!
      INTEGER
     &    i_rad_layer(nd_viewing_level)
!           Layers in which radiances are calculated
      REAL  (RealK) ::
     &    frac_rad_layer(nd_viewing_level)
!           Fractions below the tops of the layers
!
!     Variables required by DISORT:
      INCLUDE 'dec_disort.finc'
!
!
!
!     Subroutines called:
      EXTERNAL
     &    disort
!
!
!
!       Determine which layers will be required to give radiances.
        CALL set_rad_layer(ierr
     &    , n_layer, n_viewing_level, viewing_level
     &    , i_rad_layer, frac_rad_layer
     &    , nd_viewing_level
     &    )
!
      INCLUDE 'call_disort.finc'
!
!
!
      RETURN
      END
