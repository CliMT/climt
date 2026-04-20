Generating correlated-k tables
==============================

climt ships pre-built low-resolution correlated-k tables for a handful of
planetary scenarios. They are produced offline from
`G. Chaverot's Zenodo dataset <https://zenodo.org/records/16795590>`_
(CC-BY 4.0) using `exo_k` and the converter in
``scripts/generate_picket_fence_tables.py``.

The Chaverot tables are ``Ktable5d`` objects whose x-axis is H2O VMR. Pass
``--vmr`` to select the target mixing ratio (log-linear interpolation). Typical
values: Earth ~0.01, Mars ~1e-5, Venus ~3e-5.

Running the converter::

    python scripts/generate_picket_fence_tables.py \\
        --input /path/to/Earth_376ppmCO2_R500_10-30k.ktable.SI.h5 \\
        --output climt/_data/picket_fence/correlated_k/earth_low_res_lw.nc \\
        --kind lw \\
        --bands 10,500,1250,2500,3250 \\
        --ngpt 2 \\
        --vmr 0.01

See also the ``linepyline`` workflow for per-gas tables, which supports
gas-forcing experiments (e.g., CO₂ doubling).
