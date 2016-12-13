=================
Naming Quantities
=================

* Any quantity that is on vertical interface levels should be
  named using the form "<variable>_on_interface_levels", while quantities
  on mid levels should not have anything appended to the variable name.
  Quantities that do not specify "_on_interface_levels" are assumed to be on
  mid levels.
  This must be in the name rather than as an attribute so that the same
  coordinate may be present on both half and full levels (with different names).
