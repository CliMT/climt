"""Accuracy metrics for Cork correlated-k vs RRTMG on a standard Earth column.

Prints OLR, back-radiation, column RMSE, and max absolute/relative errors
for both LW and SW.  Run from the project root:

    conda run -n climt python benchmarks/benchmark_cork_vs_rrtmg.py
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
    from climt._components.cork import CorkLongwaveRadiation

    rrtmg = RRTMGLongwave()
    state_rr = get_default_state([rrtmg])
    _, d_rr = rrtmg(state_rr)
    up_rr = d_rr["upwelling_longwave_flux_in_air"].values.squeeze()
    dn_rr = d_rr["downwelling_longwave_flux_in_air"].values.squeeze()

    cork = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state_cork = get_default_state([cork])
    _, d_cork = cork(state_cork)
    up_cork = d_cork["upwelling_longwave_flux_in_air"].values.squeeze()
    dn_cork = d_cork["downwelling_longwave_flux_in_air"].values.squeeze()

    print("=== LW  (option B 4-band 2gpt  vs  RRTMG-LW 16-band) ===")
    print(f"  OLR (TOA upwelling)")
    print(f"    RRTMG : {up_rr[-1]:.2f} W/m²")
    print(f"    CORK    : {up_cork[-1]:.2f} W/m²")
    print(f"    diff  : {up_cork[-1] - up_rr[-1]:+.2f} W/m²")
    print(f"  Back-radiation (surface downwelling)")
    print(f"    RRTMG : {dn_rr[0]:.2f} W/m²")
    print(f"    CORK    : {dn_cork[0]:.2f} W/m²")
    print(f"    diff  : {dn_cork[0] - dn_rr[0]:+.2f} W/m²")
    print(f"  Column upwelling   RMSE {np.sqrt(np.mean((up_cork-up_rr)**2)):.2f} W/m²"
          f"  max {np.max(np.abs(up_cork-up_rr)):.2f} W/m²"
          f"  max-rel {np.max(np.abs(up_cork-up_rr)/np.maximum(np.abs(up_rr),1e-3)):.1%}")
    print(f"  Column downwelling RMSE {np.sqrt(np.mean((dn_cork-dn_rr)**2)):.2f} W/m²"
          f"  max {np.max(np.abs(dn_cork-dn_rr)):.2f} W/m²"
          f"  max-rel {np.max(np.abs(dn_cork-dn_rr)/np.maximum(np.abs(dn_rr),1e-3)):.1%}")


def sw_metrics():
    from climt import RRTMGShortwave, get_default_state
    from climt._components.cork import CorkShortwaveRadiation

    rrtmg = RRTMGShortwave()
    state_rr = get_default_state([rrtmg])
    state_rr["zenith_angle"].values[:] = np.pi / 4
    _, d_rr = rrtmg(state_rr)
    dn_rr = d_rr["downwelling_shortwave_flux_in_air"].values.squeeze()
    up_rr = d_rr["upwelling_shortwave_flux_in_air"].values.squeeze()

    cork = CorkShortwaveRadiation(optics="correlated_k", table="earth_low_res_sw")
    state_cork = get_default_state([cork])
    state_cork["zenith_angle"].values[:] = np.pi / 4
    _, d_cork = cork(state_cork)
    dn_cork = d_cork["downwelling_shortwave_flux_in_air"].values.squeeze()
    up_cork = d_cork["upwelling_shortwave_flux_in_air"].values.squeeze()

    print()
    print("=== SW  (3-band 2gpt  vs  RRTMG-SW 14-band,  zenith=45°) ===")
    print(f"  Surface downwelling")
    print(f"    RRTMG : {dn_rr[0]:.2f} W/m²")
    print(f"    CORK    : {dn_cork[0]:.2f} W/m²")
    print(f"    diff  : {dn_cork[0] - dn_rr[0]:+.2f} W/m²")
    print(f"  TOA upwelling (reflected)")
    print(f"    RRTMG : {up_rr[-1]:.2f} W/m²")
    print(f"    CORK    : {up_cork[-1]:.2f} W/m²")
    print(f"    diff  : {up_cork[-1] - up_rr[-1]:+.2f} W/m²")
    print(f"  Column downwelling RMSE {np.sqrt(np.mean((dn_cork-dn_rr)**2)):.2f} W/m²"
          f"  max {np.max(np.abs(dn_cork-dn_rr)):.2f} W/m²"
          f"  max-rel {np.max(np.abs(dn_cork-dn_rr)/np.maximum(np.abs(dn_rr),1e-3)):.1%}")
    print(f"  Column upwelling   RMSE {np.sqrt(np.mean((up_cork-up_rr)**2)):.2f} W/m²"
          f"  max {np.max(np.abs(up_cork-up_rr)):.2f} W/m²"
          f"  max-rel {np.max(np.abs(up_cork-up_rr)/np.maximum(np.abs(up_rr),1e-3)):.1%}")


if __name__ == "__main__":
    if not _rrtmg_available():
        print("RRTMG not available — skipping.")
    else:
        lw_metrics()
        sw_metrics()
