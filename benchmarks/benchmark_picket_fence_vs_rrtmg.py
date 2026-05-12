"""Accuracy metrics for PicketFence correlated-k vs RRTMG on a standard Earth column.

Prints OLR, back-radiation, column RMSE, and max absolute/relative errors
for both LW and SW.  Run from the project root:

    conda run -n climt python benchmarks/benchmark_picket_fence_vs_rrtmg.py
"""
import numpy as np


def _rrtmg_available():
    try:
        from climt import RRTMGLongwave, RRTMGShortwave  # noqa: F401
        return True
    except (ImportError, OSError):
        return False


def lw_metrics():
    from climt import RRTMGLongwave, get_default_state
    from climt._components.picket_fence import PicketFenceLongwave

    rrtmg = RRTMGLongwave()
    state_rr = get_default_state([rrtmg])
    _, d_rr = rrtmg(state_rr)
    up_rr = d_rr["upwelling_longwave_flux_in_air"].values.squeeze()
    dn_rr = d_rr["downwelling_longwave_flux_in_air"].values.squeeze()

    pf = PicketFenceLongwave(optics="correlated_k", table="earth_low_res_lw")
    state_pf = get_default_state([pf])
    _, d_pf = pf(state_pf)
    up_pf = d_pf["upwelling_longwave_flux_in_air"].values.squeeze()
    dn_pf = d_pf["downwelling_longwave_flux_in_air"].values.squeeze()

    print("=== LW  (option B 4-band 2gpt  vs  RRTMG-LW 16-band) ===")
    print(f"  OLR (TOA upwelling)")
    print(f"    RRTMG : {up_rr[-1]:.2f} W/m²")
    print(f"    PF    : {up_pf[-1]:.2f} W/m²")
    print(f"    diff  : {up_pf[-1] - up_rr[-1]:+.2f} W/m²")
    print(f"  Back-radiation (surface downwelling)")
    print(f"    RRTMG : {dn_rr[0]:.2f} W/m²")
    print(f"    PF    : {dn_pf[0]:.2f} W/m²")
    print(f"    diff  : {dn_pf[0] - dn_rr[0]:+.2f} W/m²")
    print(f"  Column upwelling   RMSE {np.sqrt(np.mean((up_pf-up_rr)**2)):.2f} W/m²"
          f"  max {np.max(np.abs(up_pf-up_rr)):.2f} W/m²"
          f"  max-rel {np.max(np.abs(up_pf-up_rr)/np.maximum(np.abs(up_rr),1e-3)):.1%}")
    print(f"  Column downwelling RMSE {np.sqrt(np.mean((dn_pf-dn_rr)**2)):.2f} W/m²"
          f"  max {np.max(np.abs(dn_pf-dn_rr)):.2f} W/m²"
          f"  max-rel {np.max(np.abs(dn_pf-dn_rr)/np.maximum(np.abs(dn_rr),1e-3)):.1%}")


def sw_metrics():
    from climt import RRTMGShortwave, get_default_state
    from climt._components.picket_fence import PicketFenceShortwave

    rrtmg = RRTMGShortwave()
    state_rr = get_default_state([rrtmg])
    state_rr["zenith_angle"].values[:] = np.pi / 4
    _, d_rr = rrtmg(state_rr)
    dn_rr = d_rr["downwelling_shortwave_flux_in_air"].values.squeeze()
    up_rr = d_rr["upwelling_shortwave_flux_in_air"].values.squeeze()

    pf = PicketFenceShortwave(optics="correlated_k", table="earth_low_res_sw")
    state_pf = get_default_state([pf])
    state_pf["zenith_angle"].values[:] = np.pi / 4
    _, d_pf = pf(state_pf)
    dn_pf = d_pf["downwelling_shortwave_flux_in_air"].values.squeeze()
    up_pf = d_pf["upwelling_shortwave_flux_in_air"].values.squeeze()

    print()
    print("=== SW  (3-band 2gpt  vs  RRTMG-SW 14-band,  zenith=45°) ===")
    print(f"  Surface downwelling")
    print(f"    RRTMG : {dn_rr[0]:.2f} W/m²")
    print(f"    PF    : {dn_pf[0]:.2f} W/m²")
    print(f"    diff  : {dn_pf[0] - dn_rr[0]:+.2f} W/m²")
    print(f"  TOA upwelling (reflected)")
    print(f"    RRTMG : {up_rr[-1]:.2f} W/m²")
    print(f"    PF    : {up_pf[-1]:.2f} W/m²")
    print(f"    diff  : {up_pf[-1] - up_rr[-1]:+.2f} W/m²")
    print(f"  Column downwelling RMSE {np.sqrt(np.mean((dn_pf-dn_rr)**2)):.2f} W/m²"
          f"  max {np.max(np.abs(dn_pf-dn_rr)):.2f} W/m²"
          f"  max-rel {np.max(np.abs(dn_pf-dn_rr)/np.maximum(np.abs(dn_rr),1e-3)):.1%}")
    print(f"  Column upwelling   RMSE {np.sqrt(np.mean((up_pf-up_rr)**2)):.2f} W/m²"
          f"  max {np.max(np.abs(up_pf-up_rr)):.2f} W/m²"
          f"  max-rel {np.max(np.abs(up_pf-up_rr)/np.maximum(np.abs(up_rr),1e-3)):.1%}")


def lw_convergence():
    from climt import get_default_state
    from climt._components.picket_fence import PicketFenceLongwave

    print()
    print("=== LW g-point convergence  (option B, Earth standard atm) ===")
    olrs = {}
    back = {}
    for label, table in [
        ("2gpt", "earth_low_res_lw"),
        ("4gpt", "earth_low_res_lw_4gpt"),
        ("8gpt", "earth_low_res_lw_8gpt"),
    ]:
        lw = PicketFenceLongwave(optics="correlated_k", table=table)
        state = get_default_state([lw])
        _, d = lw(state)
        olrs[label] = float(d["upwelling_longwave_flux_in_air"].values.flat[-1])
        back[label] = float(d["downwelling_longwave_flux_in_air"].values.flat[0])
        print(f"  OLR ({label}): {olrs[label]:.2f} W/m²   back-rad: {back[label]:.2f} W/m²")

    for name, vals in [("OLR", olrs), ("back-rad", back)]:
        d24 = abs(vals["4gpt"] - vals["2gpt"])
        d48 = abs(vals["8gpt"] - vals["4gpt"])
        print(f"  {name}: |4-2gpt| {d24:.3f} W/m²   |8-4gpt| {d48:.3f} W/m²"
              f"   ratio: {d24/d48:.1f}×")


if __name__ == "__main__":
    if not _rrtmg_available():
        print("RRTMG not available — skipping.")
    else:
        lw_metrics()
        sw_metrics()
        lw_convergence()
