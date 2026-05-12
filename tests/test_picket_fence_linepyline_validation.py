"""LBL-vs-picket-fence consistency on one column per scenario."""
import numpy as np
import pytest

linepyline = pytest.importorskip("linepyline")


@pytest.mark.slow
def test_olr_matches_linepyline_for_hab2():
    """TRAPPIST-1e Hab2 OLR from picket-fence within 20% of linepyline LBL."""
    from climt import load_atmospheric_properties, reset_atmospheric_properties
    from climt._components.picket_fence import PicketFenceLongwave

    # Pick one (T, p) profile
    p = np.logspace(2, 5, 30)               # 100 Pa -> 1e5 Pa
    T = np.linspace(160, 250, 30)
    Ts, ps = 260.0, 1e5

    # Reference: linepyline LBL
    rtm = linepyline.rtm(background_gas=None, surface_gravity=3.721)
    ds = rtm.radiative_transfer(
        10.0, 3250.0, 0.5, p, ps, T, Ts,
        absorbers={"CO2": 1.0}, line_shape="lorentz", binning=True,
    )
    olr_lbl = float(ds.olr.integrate("nu"))

    # Picket-fence
    load_atmospheric_properties("trappist1e")
    try:
        lw = PicketFenceLongwave(optics="correlated_k",
                                 table="trappist1e_hab2_lw")
        # ... build a state with these (T, p), call lw, extract OLR ...
        from climt import get_default_state
        state = get_default_state([lw])
        # Overwrite the default profile with ours
        state["air_temperature"].values[:, 0] = T
        state["air_pressure"].values[:, 0] = p
        state["surface_temperature"].values[:] = Ts
        tend, diag = lw(state)
        olr_pf = float(diag["upwelling_longwave_flux_in_air"].values[-1, 0])
    finally:
        reset_atmospheric_properties()

    rel_err = abs(olr_pf - olr_lbl) / olr_lbl
    assert rel_err < 0.20, f"OLR mismatch {rel_err:.1%}: pf={olr_pf:.1f}, lbl={olr_lbl:.1f}"
