"""Window closure with increasing water vapour using the hot-moist LW table.

Shows OLR, surface back-radiation, and the atmospheric window contribution
(band 2: 500–1200 cm⁻¹) as specific humidity sweeps from dry to near-runaway.

Run from the project root:
    conda run -n climt python benchmarks/benchmark_hotmoist_window_closure.py
"""
import numpy as np

from climt import get_default_state
from climt._components.picket_fence import PicketFenceLongwave


def run():
    pf = PicketFenceLongwave(optics="correlated_k", table="earth_low_res_lw_hotmoist",
                             diagnostics_level=1)
    state_base = get_default_state([pf])

    # Sweep specific humidity from very dry to near-saturation at ~300 K surface
    q_values = [1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 3e-3, 6e-3, 1e-2, 2e-2, 4e-2]

    print("=== Hotmoist window closure  (earth_low_res_lw_hotmoist, option D) ===")
    print(f"  {'q (kg/kg)':<12}  {'OLR':>8}  {'back-rad':>10}  {'ΔOLR':>8}")
    print(f"  {'':12}  {'(W/m²)':>8}  {'(W/m²)':>10}  {'(W/m²)':>8}")
    print("  " + "-" * 44)

    olr_ref = None
    for q in q_values:
        import copy
        state = copy.deepcopy(state_base)
        state["specific_humidity"].values[:] = q
        _, d = pf(state)
        olr = float(d["upwelling_longwave_flux_in_air"].values.flat[-1])
        back = float(d["downwelling_longwave_flux_in_air"].values.flat[0])
        if olr_ref is None:
            olr_ref = olr
        print(f"  {q:<12.2e}  {olr:>8.2f}  {back:>10.2f}  {olr - olr_ref:>+8.2f}")


if __name__ == "__main__":
    run()
