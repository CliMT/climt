Generating correlated-k tables
==============================

climt ships pre-built low-resolution correlated-k tables for a handful of
planetary scenarios. They are produced offline from
`G. Chaverot's Zenodo dataset <https://zenodo.org/records/16795590>`_
(CC-BY 4.0) using `exo_k` and the converter in
``scripts/generate_cork_tables.py``.

The Chaverot tables are ``Ktable5d`` objects whose x-axis is H2O VMR. Pass
``--vmr`` to select the target mixing ratio (log-linear interpolation). Typical
values: Earth ~0.01, Mars ~1e-5, Venus ~3e-5.

Running the converter::

    python scripts/generate_cork_tables.py \\
        --input /path/to/Earth_376ppmCO2_R500_10-30k.ktable.SI.h5 \\
        --output climt/_data/cork/correlated_k/earth_low_res_lw.nc \\
        --kind lw \\
        --bands 10,500,1250,2500,3250 \\
        --ngpt 2 \\
        --vmr 0.01

See also the ``linepyline`` workflow for per-gas tables, which supports
gas-forcing experiments (e.g., CO₂ doubling).

linepyline-based tables
-----------------------

For scenarios outside Chaverot's Zenodo coverage (Titan, TRAPPIST-1e,
ad-hoc compositions), climt ships a second offline pipeline using
`linepyline <https://github.com/exoclime/linepyline>`_ and HITRAN 2024 line
data. The driver is ``scripts/generate_cork_tables_linepyline.py``.

Building a Titan table (run from the ``linepyline`` conda env)::

    conda run -n linepyline python scripts/generate_cork_tables_linepyline.py \\
        --scenario titan --kind lw \\
        --output climt/_data/cork/correlated_k/titan_lw.nc

Building a custom scenario: copy one of the entries in the ``SCENARIOS``
dict at the top of the driver, edit the absorber dict / VMR grid / band
edges, and re-run. The output netCDF drops straight into climt's
``correlated_k`` data directory and is picked up by name.

Titan in particular requires HITRAN CIA flat files in
``climt/_data/cork/cia/``; see the README there for the download
list.
