"""Generate a single-band cork LW table that mimics climt's *default* gray
longwave optical depth.

``get_default_state`` initializes ``longwave_optical_depth_on_interface_levels``
to the linear profile ``tau(p) = tau0 * (1 - p/ps)`` with ``tau0 = 1`` (see the
``longwave_optical_depth`` entry in ``climt/_core/initialization.py``). Because
CO2 is well mixed, a *constant* mass-absorption coefficient makes cork's
per-layer optical depth proportional to the layer mass (proportional to dp),
whose surface->TOA cumulative integral is proportional to ``(1 - p/ps)`` -
exactly the shape of the default gray tau. So a single constant ``k`` scaled to
the right total reproduces the default gray optical depth.

Matching convention (see tests/test_grey_limit.py): cork applies a diffusivity
factor D=1.66 inside its two-stream kernel, whereas ``GrayLongwaveRadiation``
expects tau already scaled (D=1). So the target is

    D * sum_layers(k * m_layer) = tau0        (tau0 = 1 by default)

We calibrate ``k`` in two passes: build the table with k=1 to read the column
mass path ``sum_layers(m_layer)`` that cork sees in the default state, then set
``k = tau0 / (D * sum_m)``.

    python scripts/generate_gray_default_table.py \\
        --output climt/_data/cork/correlated_k/single_band_gray_lw.nc

This reproduces the *scipy-free* .nc; run scripts/convert_ck_table_to_npz.py on
it to ship the .npz the pure wheel loads.
"""
import argparse
import os

import numpy as np

DIFFUSIVITY_FACTOR = 1.66  # matches climt/_components/cork/lw/kernels.py


def write_single_band_table(out_path, k_value, band_wn=(10.0, 3250.0)):
    """Write a single-band, single-g-point LW k-table with constant opacity."""
    from scipy.io import netcdf_file

    T_grid = np.linspace(100.0, 400.0, 8).astype("f4")
    log_p_grid = np.log(np.logspace(-1, 6, 8)).astype("f4")
    nband, ngpt, nT, nP = 1, 1, T_grid.size, log_p_grid.size

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with netcdf_file(out_path, "w") as nc:
        nc.createDimension("band", nband)
        nc.createDimension("gpoint", ngpt)
        nc.createDimension("temperature", nT)
        nc.createDimension("pressure", nP)
        nc.createDimension("gas", 1)
        nc.createDimension("bounds", 2)

        k = nc.createVariable(
            "k_coefficients", "f4",
            ("gas", "band", "gpoint", "temperature", "pressure"))
        k[:] = np.full((1, nband, ngpt, nT, nP), k_value, dtype="f4")

        w = nc.createVariable("gpoint_weights", "f4", ("band", "gpoint"))
        w[:] = np.ones((nband, ngpt), dtype="f4")

        tg = nc.createVariable("temperature_grid", "f4", ("temperature",))
        tg[:] = T_grid
        pg = nc.createVariable("pressure_grid_log", "f4", ("pressure",))
        pg[:] = log_p_grid

        bl = nc.createVariable("band_wavenumber_limits", "f4", ("band", "bounds"))
        bl[:] = np.asarray([band_wn], dtype="f4")

        pf = nc.createVariable("planck_fraction", "f4",
                               ("band", "gpoint", "temperature"))
        pf[:] = np.ones((nband, ngpt, nT), dtype="f4")

        nc.gas_names = "co2"
        nc.overlap_method = "additive"
        nc.resolution = "gray_default"
        nc.source = "synthetic_gray_default_fixture"


def _column_mass_path(nc_path):
    """Return sum_layers(m_layer) cork sees in the default state (k=1 table)."""
    import sympl
    from climt import get_default_state, get_grid
    from climt._components.cork import CorkLongwaveRadiation

    sympl.set_backend(sympl.DataArrayBackend())
    cork = CorkLongwaveRadiation(
        optics="correlated_k",
        table=nc_path,
    )
    state = get_default_state([cork], grid_state=get_grid(nx=1, ny=1, nz=28))
    _, diag = cork(state)
    tau_layer = diag["longwave_optical_depth_per_band"].values[:, 0, 0, 0]
    return float(tau_layer.sum())  # k=1 => tau_layer == m_layer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--total-optical-depth", type=float, default=1.0,
                    help="Target D-scaled column optical depth (default gray tau0=1).")
    ap.add_argument("--diffusivity", type=float, default=DIFFUSIVITY_FACTOR,
                    help="Diffusivity factor D the component will run with. "
                         "The calibration target is D * sum(k * m) == "
                         "--total-optical-depth, so this MUST match the D "
                         "passed to CorkLongwaveRadiation. The EC2213 notes "
                         "use D=2; cork's default is 1.66.")
    args = ap.parse_args()

    # Pass 1: write with k=1 and read the column mass path cork integrates.
    write_single_band_table(args.output, k_value=1.0)
    sum_m = _column_mass_path(args.output)

    # Pass 2: calibrate k so D * sum(k * m) == target total.
    k_value = args.total_optical_depth / (args.diffusivity * sum_m)
    write_single_band_table(args.output, k_value=k_value)
    print(f"wrote {args.output}")
    print(f"  column mass path sum_m = {sum_m:.6f}")
    print(f"  calibrated k = {k_value:.6g}  "
          f"(raw total tau = {k_value*sum_m:.6f}, "
          f"D-scaled at D={args.diffusivity} = "
          f"{args.diffusivity*k_value*sum_m:.6f})")


if __name__ == "__main__":
    main()
