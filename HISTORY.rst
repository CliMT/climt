=======
History
=======

Latest
------
* Added snow-ice component

v.0.9
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
* SlabSurface no longer uses depth_slab_surface as input
* changed order of outputs of GfsDynamicalCore and SimplePhysics to conform
  to TimeStepper order of diagnostics, new_state
* get_default_state now accepts mid_levels and interface_levels instead of z
  to specify vertical coordinates.
* mass_to_volume_mixing_ratio now uses numpy arrays instead of DataArrays.
