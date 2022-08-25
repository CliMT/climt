=======
History
=======

v.0.16.11
---------

* New component BucketHydrology that implements Manabe first generation land model
* BucketHydrology calculates the sensible and latent heat flux within the component
* Conservation test for the component also added
* Moving CI to Github Actions

v.0.16.8
--------

* Fix timeout for all MAC builds

v.0.16.6
--------

* Prevent MAC OS builds from timing out

v.0.16.5
--------

* Fix formatting errors which prevent pypi deployment

v.0.16.4
--------

* Fix MCICA for the shortwave component of RRTMG
* Revise random number generation for MCICA
* Improvement of the user interface to control MCICA

v.0.16.3
--------

* update numpy requirement to avoid binary incompatibility error
* Fix error in documentation

v.0.16.2
--------

* Fix wheel build on Mac

v.0.16.1
--------

* Fixed issue with Mac build
* Few changes in the dry convection component. Significantly improves the performance.
* Changed logo!
* Fixed failing docs build

v0.16.0
-------

* Added some documentation for using RRTMG with McICA
* CI Testing for Mac and py37 added.
* Refactored initialisation code
* Enable the McICA version of RRTMG Longwave for consistency
  with the Shortwave component.
* Fix bugs in IceSheet
* Add tests to verify conservation of quantities
* Fix bugs in initialisation
* Fix energy conservation in surface flux scheme
* Enable the McICA version of RRTMG Shortwave,
  so that partial cloud fractions can be used.
* Add GMD example scripts to repository.
* Fix docs to reflect API changes after refactor.
* Fix wrong initialisation to use sigma values instead of pressure values
  of optical depth for GrayLongwaveRadiation

Breaking Changes
----------------

* The flux outputs of GrayLongwaveRadiation have been renamed to eliminate
  `on_interface_levels` to keep consistency with other components.
* All arrays are now 3/2d by default based on their expected dimensions.
* horizontal dimensions are now `lon`, `lat`, but inputs
  used by components remain the same (`latitude`, `longitude`).



v.0.14.8
--------

Many of the changes in this version come from changes in Sympl 0.4.0. We recommend
reading those changes in the Sympl documentation.

* Updated component APIs to work with Sympl 0.4.0
* Many components which previously required horizontal dimensions now use
  wildcard matches for column dimensions.
* Switched many print statements to logging calls.
* Fixed bugs in some components

Breaking Changes
----------------

* get_constant and set_constant have been removed, use the ones in Sympl.
* Emanuel convection scheme can no longer be set to perform dry adiabatic
  adjustment to the boundary layer. This has been implemented in a separate
  component.
* ClimtPrognostic, ClimtImplicitPrognostic, ClimtDiagnostic, ClimtImplicit have
  been removed. Use the base types in Sympl.
* State initialization has been entirely re-worked. get_default_state now takes in
  an optional grid state instead of options to do with the state grid. A function
  get_grid is provided which can create a grid state, or one can be created manually.
  A grid state is a state containing air pressure and sigma on mid and interface
  levels, as well as surface pressure.
* Replaced references to "thermal_capacity" with references to "heat_capacity" in
  component quantity names.

v.0.14.7
--------

* Fix issue with pip v10 and pandas 0.22 conflicts

v.0.14.3
--------

* Fix release issue because of pip API change

v.0.14.1
--------
* Fix appveyor fail due to pip changes

v.0.14.0
--------

* Fixed broken version numbers

v.0.12.0
--------

* new release to fix version numbers and create zenodo ID

v.0.9.4
-------

* Added attributes to inputs/outputs/ etc., to work with ScalingWrapper
  Added tests as well.
* Added tests for constants functions
* Fixed requirements to ensure this version of climt installs
  the correct versions of sympl and numpy.

v.0.9.3
-------

* Released because of a labelling issue. See 0.9.2 for details.

v.0.9.2
--------
* Updated documentation
* Cleaned up examples
* Added (*)_properties as a property to all components
* The gas constant for dry air in the Emanuel scheme is now renamed _Rdair
* RRTMG LW and SW are now OpenMP parallel
* Added Instellation component to calculate zenith angle
* Added tests to increase coverage
* New constants handling functionality added
* Travis builds now use stages
* Appveyor CI up and running
* Pre-installation of cython and numpy no longer necessary for source builds
* Added snow-ice component
* Ozone profiles do not need to be specified externally
* Now also tested on Python 3.6

Breaking Changes
----------------

* API for constants setting changed to `set_constant_from_dict` and `add_constants_from_dict`
* `GfsDynamicalCore` renamed to `GFSDynamicalCore` for consistency
* `get_prognostic_version` method of `ClimtImplicit` renamed to `prognostic_version`, and
  no longer accepts timestep as an argument. The current timestep should be set in
  `ClimtImplicit.current_time_step` during each iteration.
* `RRTMGShortwave` now uses sympl's solar constant by default instead of from fortran.

v.0.9.1
-------
* Held-Suarez and moist GCM with grey radiation work!
* Added DCMIP initial conditions, test 4 tried out.
* Dynamical core integrated now.
* BIG change in the build system. Tests pass on Mac as well
* Arrays can now have arbitrary dtype (to use qualitative, string, quantities)
* Added Emanuel Convection, surface energy balance model and ice sheet energy balance
* 2D coordinates are now supported for horizontal coordinates
* Replaced create_output_arrays() with a more general
  get_state_dict_for() and get_numpy_arrays_from_state()
  combination.
* State arrays now have coordinates
* Updated documentation
* RTD finally working, phew!
* Added RRTMG Longwave, Simple Physics
* Added helper functions to reduce boilerplate code in components

Breaking Changes
----------------

Latest
-------

* method to obtain piecewise constant prognostic has been renamed to
  :code:`piecewise_constant_version`
* Ozone profile has been modified
* Heating rate for RRTMG top-of-atmosphere is no longer manually set to zero
* Components no longer accept constants during initialisation. All constant handling
  is done internally.

v.0.9
------
* SlabSurface no longer uses depth_slab_surface as input
* changed order of outputs of GfsDynamicalCore and SimplePhysics to conform
  to TimeStepper order of diagnostics, new_state
* get_default_state now accepts mid_levels and interface_levels instead of z
  to specify vertical coordinates.
* mass_to_volume_mixing_ratio now uses numpy arrays instead of DataArrays.
