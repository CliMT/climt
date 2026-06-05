"""Window closure with increasing water vapour using the hot-moist LW table.

Shows broadband OLR / back-radiation and a per-band breakdown as specific
humidity sweeps from dry to near-runaway.  Band structure (option D):
  Band 0:  10– 500 cm⁻¹  H₂O rotation (far IR, always opaque)
  Band 1: 500–1200 cm⁻¹  CO₂ 15µm + atmospheric window  ← closes with H₂O
  Band 2: 1200–2000 cm⁻¹  H₂O vibration-rotation
  Band 3: 2000–3250 cm⁻¹  H₂O near-IR (controls runaway limit)

Run from the project root:
    conda run -n climt python benchmarks/benchmark_hotmoist_window_closure.py
"""
import copy

import numpy as np

from climt import get_default_state
from climt._components.picket_fence import PicketFenceLongwave

BAND_LABELS = [
    "10–500",
    "500–1200",
    "1200–2000",
    "2000–3250",
]

Q_VALUES = [1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 3e-3, 6e-3, 1e-2, 2e-2, 4e-2]


def run():
    pf = PicketFenceLongwave(optics="correlated_k", table="earth_low_res_lw_hotmoist")
    state_base = get_default_state([pf])

    rows = []
    for q in Q_VALUES:
        state = copy.deepcopy(state_base)
        state["specific_humidity"].values[:] = q
        _, d = pf(state)

        olr  = float(d["upwelling_longwave_flux_in_air"].values.flat[-1])
        back = float(d["downwelling_longwave_flux_in_air"].values.flat[0])

        # per-band: shape (levels, lat, lon, band) — take TOA for OLR, sfc for back-rad
        up_band = d["upwelling_longwave_flux_in_air_per_band"].values
        dn_band = d["downwelling_longwave_flux_in_air_per_band"].values
        olr_band  = up_band[-1, 0, 0, :]   # TOA, squeezed
        back_band = dn_band[0,  0, 0, :]   # surface, squeezed

        rows.append(dict(q=q, olr=olr, back=back,
                         olr_band=olr_band, back_band=back_band))

    olr_ref  = rows[0]["olr"]
    back_ref = rows[0]["back"]

    # ── broadband summary ─────────────────────────────────────────────────────
    print("=== Broadband  (option D hotmoist, Earth std atm) ===")
    hdr = f"  {'q (kg/kg)':<12}  {'OLR':>8}  {'ΔOLR':>7}  {'back-rad':>9}  {'Δback':>7}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        print(f"  {r['q']:<12.2e}"
              f"  {r['olr']:>8.2f}"
              f"  {r['olr']-olr_ref:>+7.2f}"
              f"  {r['back']:>9.2f}"
              f"  {r['back']-back_ref:>+7.2f}")

    # ── per-band OLR ─────────────────────────────────────────────────────────
    print()
    print("=== Per-band OLR at TOA  (W/m²) ===")
    band_hdr = f"  {'q (kg/kg)':<12}" + "".join(f"  {lb:>12}" for lb in BAND_LABELS)
    print(band_hdr)
    print("  " + "-" * (len(band_hdr) - 2))
    olr_band_ref = rows[0]["olr_band"]
    for r in rows:
        deltas = "".join(
            f"  {v:>6.1f}({v-ref:>+5.1f})"
            for v, ref in zip(r["olr_band"], olr_band_ref)
        )
        print(f"  {r['q']:<12.2e}{deltas}")

    # ── per-band back-radiation ───────────────────────────────────────────────
    print()
    print("=== Per-band back-radiation at surface  (W/m²) ===")
    print(band_hdr)
    print("  " + "-" * (len(band_hdr) - 2))
    back_band_ref = rows[0]["back_band"]
    for r in rows:
        deltas = "".join(
            f"  {v:>6.1f}({v-ref:>+5.1f})"
            for v, ref in zip(r["back_band"], back_band_ref)
        )
        print(f"  {r['q']:<12.2e}{deltas}")


if __name__ == "__main__":
    run()
