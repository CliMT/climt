"""G-point resolution convergence for the Earth LW option-B correlated-k tables.

Verifies that as g-points per band increase (2 → 4 → 8), the broadband OLR
converges monotonically toward the high-resolution limit. All three tables use
identical band edges (10/500/800/1800/3250 cm⁻¹); only the quadrature order
changes, so any non-monotone behaviour would indicate a table-generation bug.
"""
import numpy as np
import pytest

from climt import get_default_state
from climt._components.picket_fence import PicketFenceLongwave


@pytest.fixture(scope="module")
def olr_by_gpt():
    """OLR (W/m²) for the 2-, 4-, and 8-gpt Earth LW option-B tables."""
    results = {}
    for label, table in [
        ("2gpt", "earth_low_res_lw"),
        ("4gpt", "earth_low_res_lw_4gpt"),
        ("8gpt", "earth_low_res_lw_8gpt"),
    ]:
        lw = PicketFenceLongwave(optics="correlated_k", table=table)
        state = get_default_state([lw])
        _, diag = lw(state)
        results[label] = float(
            diag["upwelling_longwave_flux_in_air"].values.flat[-1]
        )
    return results


def test_olr_converges_monotonically(olr_by_gpt):
    """OLR difference must shrink as g-points double: |4gpt−2gpt| > |8gpt−4gpt|."""
    diff_2_4 = abs(olr_by_gpt["4gpt"] - olr_by_gpt["2gpt"])
    diff_4_8 = abs(olr_by_gpt["8gpt"] - olr_by_gpt["4gpt"])
    assert diff_4_8 < diff_2_4, (
        f"OLR not converging: |4-2gpt|={diff_2_4:.3f} W/m², "
        f"|8-4gpt|={diff_4_8:.3f} W/m²"
    )


def test_olr_values_reported(olr_by_gpt, capsys):
    """Print OLR at each resolution so the numbers appear in the test log."""
    for label, olr in olr_by_gpt.items():
        print(f"  OLR ({label}): {olr:.2f} W/m²")
    diff_2_4 = abs(olr_by_gpt["4gpt"] - olr_by_gpt["2gpt"])
    diff_4_8 = abs(olr_by_gpt["8gpt"] - olr_by_gpt["4gpt"])
    print(f"  |4-2gpt|: {diff_2_4:.3f} W/m²  |8-4gpt|: {diff_4_8:.3f} W/m²")
