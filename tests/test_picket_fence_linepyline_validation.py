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

    # Load planet constants before building either the LBL or PF state
    load_atmospheric_properties("trappist1e")
    try:
        from sympl import get_constant
        g = get_constant("gravitational_acceleration", "m/s^2")

        # Reference: linepyline LBL
        rtm = linepyline.rtm(background_gas=None, surface_gravity=g)
        ds = rtm.radiative_transfer(
            10.0, 3250.0, 0.5, p, ps, T, Ts,
            absorbers={"CO2": 1.0}, line_shape="lorentz", binning=True,
        )
        olr_lbl = float(ds.olr.integrate("nu"))

        # Picket-fence
        lw = PicketFenceLongwave(optics="correlated_k",
                                 table="trappist1e_hab2_lw")
        from climt import get_default_state, get_grid
        state = get_default_state([lw], grid_state=get_grid(nz=30))
        # climt convention is index 0 = surface (pressure decreases with
        # index); linepyline's p/T arrays here are bottom-up the opposite way,
        # so flip them before assigning.
        p_climt = p[::-1]
        T_climt = T[::-1]
        state["air_temperature"].values[:, 0, 0] = T_climt
        state["air_pressure"].values[:, 0, 0] = p_climt
        state["surface_temperature"].values[0, 0] = Ts
        # Build interface pressures consistent with the layer profile.
        p_int = np.empty(len(p_climt) + 1)
        p_int[1:-1] = np.sqrt(p_climt[:-1] * p_climt[1:])
        p_int[0] = ps
        p_int[-1] = p_climt[-1] ** 2 / p_int[-2]
        state["air_pressure_on_interface_levels"].values[:, 0, 0] = p_int
        tend, diag = lw(state)
        olr_pf = float(diag["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    finally:
        reset_atmospheric_properties()

    rel_err = abs(olr_pf - olr_lbl) / olr_lbl
    assert rel_err < 0.20, f"OLR mismatch {rel_err:.1%}: pf={olr_pf:.1f}, lbl={olr_lbl:.1f}"
