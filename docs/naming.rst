=================
Naming Quantities
=================

* Vertical coordinates should specify whether they are on full or half levels
  using the form "<variable>_on_<full/half>_levels",
  e.g. "pressure_on_full_levels", "height_on_half_levels", etc.
  This must be in the name rather than as an attribute so that the same
  coordinate may be present on both half and full levels (with different names).
* Each dimension name (e.g. "latitude" or "height") should have a coordinate
  variable corresponding to that dimension. This means if the model is using
  a vertical coordinate of "air_pressure_on_full_levels", the vertical dimension
  on that DataArray should be named "air_pressure_on_full_levels" and not
  "height" or "vertical" or "z".
