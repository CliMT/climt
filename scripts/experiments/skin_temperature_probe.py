"""Does cork's non-grey LW admit an isothermal skin-temperature solution?

Physics under test: for an optically thin layer the energy balance is

    sum_b eps_b * F_up_b  =  2 * sum_b eps_b * pi*B_b(T),    eps_b = k_b * m

The layer mass m cancels, so every thin layer that sees the same upwelling
spectrum must equilibrate to the SAME temperature -- gray or non-grey. Non-grey
changes what that temperature is (weighting toward the strong bands), not that
it is height-independent.

Two probes, neither of which needs a converged time integration:

  1. --budget: per-band optical depth and layer energy budget of the top
     layers in the default state, so we can see which layers are actually thin.

  2. --sweep: hold the top region isothermal at T_top, sweep T_top, and record
     each layer's LW heating rate. If the model admits a skin temperature, all
     thin layers cross zero heating at the SAME T_top. If their crossings fan
     out with height, cork has no isothermal solution and the RCE run will keep
     cooling the lid forever.

    conda run -n climt python scripts/experiments/skin_temperature_probe.py --budget
    conda run -n climt python scripts/experiments/skin_temperature_probe.py --sweep
"""
import argparse

import numpy as np
import sympl

from climt import (CorkLongwaveRadiation, get_default_state, get_grid)

SIGMA = 5.670374419e-8
NZ = 18
TABLE = "earth_low_res_lw"


def build(diagnostics_level=1, co2=None):
    lw = CorkLongwaveRadiation(optics="correlated_k", table=TABLE,
                               diagnostics_level=diagnostics_level)
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=NZ))
    if co2 is not None:
        # The table's co2_vmr axis bottoms out at 1e-5 mol/mol; staying on-grid
        # keeps this an opacity change rather than an extrapolation artifact.
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = co2
    return lw, state


def budget(co2=None):
    lw, state = build(co2=co2)
    _, diag = lw(state)
    p = state["air_pressure"].values[:, 0, 0] / 100.0
    T = state["air_temperature"].values[:, 0, 0]
    q = state["specific_humidity"].values[:, 0, 0]
    tau = diag["longwave_optical_depth_per_band"].values[:, 0, 0, :]     # (nlev, nband)
    up = diag["upwelling_longwave_flux_in_air_per_band"].values[:, 0, 0, :]
    dn = diag["downwelling_longwave_flux_in_air_per_band"].values[:, 0, 0, :]

    print(f"{'lev':>3} {'p(hPa)':>9} {'T(K)':>7} {'q(kg/kg)':>10} "
          f"{'tau_tot':>9} {'max_tau_b':>10} {'thin?':>6}")
    for k in range(NZ - 1, -1, -1):
        tt = tau[k].sum()
        print(f"{k:3d} {p[k]:9.2f} {T[k]:7.2f} {q[k]:10.3e} "
              f"{tt:9.4f} {tau[k].max():10.4f} "
              f"{'yes' if tau[k].max() < 0.1 else 'no':>6}")

    print("\nper-band optical depth, top 4 layers (band: 1..14)")
    for k in range(NZ - 1, NZ - 5, -1):
        print(f"  lev {k:2d}: " + " ".join(f"{v:7.4f}" for v in tau[k]))

    print("\nspectral SHAPE of the layer opacity (tau_b / sum_b tau_b), top 4")
    for k in range(NZ - 1, NZ - 5, -1):
        s = tau[k] / tau[k].sum()
        print(f"  lev {k:2d}: " + " ".join(f"{v:7.4f}" for v in s))

    print("\nupwelling flux per band entering each layer from below, top 4")
    for k in range(NZ - 1, NZ - 5, -1):
        print(f"  lev {k:2d}: " + " ".join(f"{v:7.3f}" for v in up[k]))
    print("\nband-integrated: up(bottom of layer), dn(top of layer)")
    for k in range(NZ - 1, NZ - 5, -1):
        print(f"  lev {k:2d}: up={up[k].sum():8.3f}  dn={dn[k + 1].sum():8.3f}")


def sweep(k_lo=10, n=61, t_min=100.0, t_max=280.0, co2=None):
    """Hold levels k_lo.. top isothermal at T_top; find each layer's zero-heating T.

    Note the confound this design carries when the block is optically THICK:
    the flux arriving at the top of the block is then emitted by the block
    itself, so it tracks B(T_top) and no crossing can exist. The test is only
    meaningful once the block is thin -- which is exactly the regime the
    isothermal-skin argument is about, reachable here via --co2.
    """
    lw, state = build(co2=co2)
    T0 = state["air_temperature"].values[:, 0, 0].copy()
    p = state["air_pressure"].values[:, 0, 0] / 100.0

    temps = np.linspace(t_min, t_max, n)
    H = np.zeros((n, NZ))
    for i, t in enumerate(temps):
        state["air_temperature"].values[k_lo:, 0, 0] = t
        _, diag = lw(state)
        H[i] = diag["air_temperature_tendency_from_longwave"].values[:, 0, 0]
    state["air_temperature"].values[:, 0, 0] = T0

    print(f"isothermal block: levels {k_lo}..{NZ - 1}   "
          f"(sweeping T_top over {t_min:.0f}-{t_max:.0f} K)\n")
    print(f"{'lev':>3} {'p(hPa)':>9} {'T_zero-heating(K)':>19}")
    roots = {}
    for k in range(NZ - 1, k_lo - 1, -1):
        h = H[:, k]
        s = np.where(np.diff(np.sign(h)) != 0)[0]
        if len(s) == 0:
            print(f"{k:3d} {p[k]:9.2f} {'no crossing':>19}   "
                  f"(H spans {h.min():+.3g} .. {h.max():+.3g} K/day)")
            continue
        j = s[0]
        root = np.interp(0.0, [h[j], h[j + 1]], [temps[j], temps[j + 1]]) \
            if h[j] < h[j + 1] else \
            np.interp(0.0, [h[j + 1], h[j]], [temps[j + 1], temps[j]])
        roots[k] = root
        print(f"{k:3d} {p[k]:9.2f} {root:19.2f}")

    if len(roots) > 1:
        v = np.array(list(roots.values()))
        print(f"\nspread across the block: {v.min():.2f} - {v.max():.2f} K "
              f"(range {v.max() - v.min():.2f} K)")
        print("A skin temperature exists only if this spread is small.")
    print(f"\nanalytic gray reference: (OLR/2sigma)^(1/4) with OLR=240 "
          f"-> {(240.0 / (2 * SIGMA)) ** 0.25:.2f} K")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget", action="store_true")
    ap.add_argument("--sweep", action="store_true")
    ap.add_argument("--k-lo", type=int, default=10)
    ap.add_argument("--co2", type=float, default=None,
                    help="override CO2 mole fraction (default state: 3.3e-4; "
                         "table axis bottoms out at 1e-5)")
    args = ap.parse_args()
    sympl.set_backend(sympl.DataArrayBackend())
    if args.budget:
        budget(co2=args.co2)
    if args.sweep:
        sweep(k_lo=args.k_lo, co2=args.co2)
