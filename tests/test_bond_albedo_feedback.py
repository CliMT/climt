import numpy as np
import sympl

from climt._components.picket_fence import PicketFenceShortwave


def _make_state_with_high_albedo(sw):
    from climt import get_default_state, get_grid
    sympl.set_backend(sympl.DataArrayBackend())
    grid = get_grid(nx=1, ny=1, nz=20)
    state = get_default_state([sw], grid_state=grid)
    state["zenith_angle"].values[:] = np.pi / 4
    state["surface_albedo_for_direct_shortwave"].values[:] = 0.8
    # Non-zero T_irr so A_B shifts T_eff and therefore optics
    state["irradiation_temperature"].values[:] = 1500.0
    return state


def test_bond_albedo_feedback_changes_olr():
    """Enabling bond_albedo_feedback shifts the TOA upward/downward fluxes."""
    sw_no_feedback = PicketFenceShortwave(
        optics="parmentier", bond_albedo_feedback=False
    )
    sw_feedback = PicketFenceShortwave(
        optics="parmentier", bond_albedo_feedback=True
    )

    s1 = _make_state_with_high_albedo(sw_no_feedback)
    s2 = _make_state_with_high_albedo(sw_feedback)

    _, d1 = sw_no_feedback(s1)
    _, d2 = sw_feedback(s2)

    up_no = d1["upwelling_shortwave_flux_in_air"].values[-1, :]
    up_yes = d2["upwelling_shortwave_flux_in_air"].values[-1, :]
    # Feedback must shift the TOA flux — equal to four decimals would be a bug.
    assert not np.allclose(up_no, up_yes, rtol=1e-4)
