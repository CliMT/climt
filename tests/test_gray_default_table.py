"""``single_band_gray_lw`` reproduces climt's *default* gray longwave radiation.

``get_default_state`` seeds ``longwave_optical_depth_on_interface_levels`` with
the linear profile ``tau0 * (1 - p/ps)`` (tau0 = 1). The ``single_band_gray_lw``
cork table is a constant-k, single-band fixture calibrated (see
``scripts/generate_gray_default_table.py``) so that cork's diffusivity-scaled
column optical depth equals that default gray tau. Running the two components on
the default state must then give the same heating rates and OLR.
"""
import numpy as np

from climt import GrayLongwaveRadiation, get_default_state, get_grid
from climt._components.cork import CorkLongwaveRadiation

DIFFUSIVITY_FACTOR = 1.66  # matches climt/_components/cork/lw/kernels.py


def test_gray_default_table_matches_default_gray_radiation():
    gray = GrayLongwaveRadiation()
    cork = CorkLongwaveRadiation(optics="correlated_k", table="single_band_gray_lw")
    state = get_default_state([gray, cork], grid_state=get_grid(nx=1, ny=1, nz=28))

    tend_gray, diag_gray = gray(state)   # consumes default-state tau (tau0=1)
    tend_cork, diag_cork = cork(state)

    # cork's diffusivity-scaled cumulative optical depth reproduces default tau.
    tau_layer = diag_cork["longwave_optical_depth_per_band"].values[:, 0, 0, 0]
    tau_cork_int = np.concatenate([[0.0], np.cumsum(DIFFUSIVITY_FACTOR * tau_layer)])
    tau_default = state["longwave_optical_depth_on_interface_levels"].values[:, 0, 0]
    np.testing.assert_allclose(tau_cork_int, tau_default, atol=2e-3)

    # Same longwave heating rates and OLR.
    np.testing.assert_allclose(
        tend_cork["air_temperature"].values,
        tend_gray["air_temperature"].values,
        rtol=1e-2, atol=1e-7,
    )
    olr_gray = diag_gray["upwelling_longwave_flux_in_air"].values.flat[-1]
    olr_cork = diag_cork["upwelling_longwave_flux_in_air"].values.flat[-1]
    np.testing.assert_allclose(olr_cork, olr_gray, rtol=1e-3)


def test_generator_honours_diffusivity_and_absolute_paths(tmp_path):
    """Calibration must use the same D the component will run with.

    The EC2213 notes' tau_inf is the diffusivity-SCALED column optical depth,
    so a table built for D=2 must satisfy 2 * sum(k * m) == tau_inf.
    """
    import subprocess
    import sys

    import numpy as np
    import sympl

    import climt
    from climt import CorkLongwaveRadiation, get_default_state, get_grid

    out = tmp_path / "tour_gray_lw.nc"
    tau_inf, D = 4.0, 2.0
    subprocess.run(
        [sys.executable, "scripts/generate_gray_default_table.py",
         "--output", str(out), "--total-optical-depth", str(tau_inf),
         "--diffusivity", str(D)],
        check=True,
    )

    sympl.set_backend(climt.UnytBackend())
    lw = CorkLongwaveRadiation(optics="correlated_k", table=str(out),
                               diffusivity_factor=D)
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=28))
    _, diag = lw(state)
    tau_col = float(diag["longwave_optical_depth_per_band"].values[:, 0, 0, 0].sum())
    np.testing.assert_allclose(D * tau_col, tau_inf, rtol=1e-3)
