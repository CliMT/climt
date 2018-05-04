=======
History
=======

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
