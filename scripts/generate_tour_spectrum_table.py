#!/usr/bin/env python
"""Generate the Modelling Tour's high-resolution longwave spectrum table.

Pages 1-3 of ``docs/modelling-tour/`` draw an OLR spectrum and read a radiating
level off it. The shipped 14-band ``earth_low_res_lw`` table is the right tool
for anything quantitative or swept -- it is what every test pins -- but 14 bars
do not look like a spectrum, and a 180 cm^-1 window band cannot show the window
*narrowing*. This script builds a ~100-band sibling for the picture.

It is a thin wrapper over ``scripts/generate_cork_tables_linepyline.py``'s
``earth_hifi`` scenario, changing exactly two things:

* **band structure** -- each of the 14 shipped bands is split into ``--nsub``
  equal-wavenumber pieces, so every existing edge survives. The hi-res grid is a
  strict *refinement* of the reference grid, which is what lets
  ``tests/test_spectrum_table.py`` validate band-by-band by plain summation,
  with no interpolation error anywhere in the comparison.
* **CO2 axis** -- thinned from 10 nodes to ``--co2-nodes`` over the same
  10-10000 ppm range. The table is a browser download, and the CO2 axis is a
  whole multiplicative factor on its size; pages 1-3 never sweep CO2 (page 5
  does, on the shipped table), so they only need enough nodes to sit at 280 ppm
  without interpolation error. The gate in ``tests/test_spectrum_table.py``
  checks the 280 -> 560 ppm forcing survives the thinning anyway.

Everything else -- the (T, p, X_H2O) grids, the pseudovoigt line shape, the
decoupled MT_CKD continuum -- is inherited from ``earth_hifi``, so the two
tables differ only in resolution and not in physics.

Usage (needs the ``linepyline`` conda env, which owns the HITRAN line data):

    conda run -n linepyline python scripts/generate_tour_spectrum_table.py \\
        --output /tmp/earth_spectrum_lw.nc

The defaults reproduce the shipped table exactly: 8 g-points over 4 sub-bands
each, i.e. 56 bands. ``--ngpt 4 --nsub 7`` (98 bands) was tried first and is
equally valid on every test except the per-band OLR one, where 4 g-points cannot
resolve a band's spread of absorption strengths; see Task 13's log in
``docs/superpowers/plans/2026-08-12-modelling-tour-radiation.md``.

then convert to the numpy-native format the browser reads (``climt`` env, since
the converter goes through climt's netCDF reader):

    conda run -n climt python scripts/convert_ck_table_to_npz.py \\
        /tmp/earth_spectrum_lw.nc
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys

import numpy as np

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)
from generate_cork_tables_linepyline import SCENARIOS, build_table  # noqa: E402

# The scenario the shipped `earth_low_res_lw` table was built from -- see its
# `source` attribute, `linepyline:earth_hifi`. Not `earth`, which is a 6-band
# structure with CO2 baked in at 376 ppm and no runtime axis.
SCENARIO = "earth_hifi"


def refine_band_edges(edges, nsub):
    """Split every band into ``nsub`` equal-wavenumber pieces.

    Args:
        edges: (nband+1,) reference band edges in cm^-1, ascending.
        nsub: sub-divisions per reference band.

    Returns:
        (nband*nsub + 1,) edges whose set contains every input edge exactly.
    """
    edges = np.asarray(edges, dtype=float)
    pieces = [np.linspace(lo, hi, nsub + 1)[:-1]
              for lo, hi in zip(edges[:-1], edges[1:])]
    return np.concatenate(pieces + [edges[-1:]])


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output", required=True,
                    help="output .nc path")
    ap.add_argument("--ngpt", type=int, default=8,
                    help="g-points per band (default 8, as shipped; 4 is too "
                         "few to resolve a band's spread of absorption)")
    ap.add_argument("--nsub", type=int, default=4,
                    help="sub-divisions per shipped band (default 4 -> 56 "
                         "bands, as shipped)")
    ap.add_argument("--co2-nodes", type=int, default=5,
                    help="log-spaced CO2 nodes over 10-10000 ppm (default 5; "
                         "the shipped table uses 10)")
    ap.add_argument("--dnu", type=float, default=0.1,
                    help="line-by-line wavenumber resolution, cm^-1")
    args = ap.parse_args()

    cfg = SCENARIOS[SCENARIO]
    band_edges = refine_band_edges(cfg["lw_band_edges"], args.nsub)
    nband = len(band_edges) - 1

    # Same 10-10000 ppm span as the shipped table, fewer nodes.
    co2_grid = np.logspace(-5, -2, args.co2_nodes)

    print(f"[tour-spectrum] {nband} bands "
          f"({len(cfg['lw_band_edges']) - 1} x {args.nsub}), "
          f"ngpt={args.ngpt}, {args.co2_nodes} CO2 nodes "
          f"({co2_grid[0] * 1e6:.0f}-{co2_grid[-1] * 1e6:.0f} ppm)")
    print(f"[tour-spectrum] narrowest band = "
          f"{np.diff(band_edges).min():.1f} cm^-1 "
          f"({np.diff(band_edges).min() / args.dnu:.0f} line-by-line points)")

    # build_table reads the CO2 axis off the scenario registry, so thinning it
    # means editing the registry entry for the duration of this call.
    original_co2 = cfg["co2_vmr_grid"]
    cfg["co2_vmr_grid"] = co2_grid
    try:
        build_table(
            SCENARIO, "lw", args.output,
            ngpt=args.ngpt,
            dnu=args.dnu,
            line_shape="pseudovoigt",
            binning=False,
            band_edges=band_edges,
            include_mtckd_continuum=True,
            # The shipped table stores the H2O continuum band-grey and separate
            # from the line k-distribution; matching that is what makes the two
            # tables comparable at all.
            decouple_continuum=True,
        )
    finally:
        cfg["co2_vmr_grid"] = original_co2

    size = os.path.getsize(args.output)
    print(f"[tour-spectrum] {args.output}  {size / 1e6:.2f} MB")
    print(f"[tour-spectrum] sha256 {_sha256(args.output)}")


if __name__ == "__main__":
    main()
