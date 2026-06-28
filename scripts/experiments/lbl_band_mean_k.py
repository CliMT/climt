"""Run linepyline LBL on a column profile (isothermal or dry adiabat).

Saves kappa(p, nu) and T(p) so a downstream script (in the climt env) can
integrate Planck-weighted band-mean opacities over arbitrary band edges
and compare to CORK correlated-k tables.

Usage (linepyline env):
    conda run -n linepyline python examples/lbl_band_mean_k.py --profile iso
    conda run -n linepyline python examples/lbl_band_mean_k.py --profile adiabat
"""
import argparse

import numpy as np
import xarray as xr

import linepyline as lpl

CO2_VMR = 376e-6
NU_MIN = 10.0
NU_MAX = 3250.0
DNU = 0.1

# Dry adiabat parameters
T_SURF_ADIABAT = 288.0
P_SURF_ADIABAT = 101325.0
KAPPA = 2.0 / 7.0  # R/cp for diatomic ideal gas
T_FLOOR = 80.0  # clip to keep within CORK table T_grid (min 50 K)

# Isothermal
T_ISO = 250.0


def build_profile(profile):
    p = np.logspace(np.log10(10.0), np.log10(101325.0), 60)  # Pa
    if profile == "iso":
        T = np.full_like(p, T_ISO)
        Ts = T_ISO
    elif profile == "adiabat":
        T = np.maximum(T_FLOOR, T_SURF_ADIABAT * (p / P_SURF_ADIABAT) ** KAPPA)
        Ts = T_SURF_ADIABAT
    else:
        raise ValueError(profile)
    return p, T, Ts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", choices=["iso", "adiabat"], required=True)
    args = ap.parse_args()

    p, T, Ts = build_profile(args.profile)
    ps = p[-1]

    p_da = xr.DataArray(p, dims=("p",), coords={"p": p})
    T_da = xr.DataArray(T, dims=("p",), coords={"p": p})

    rtm = lpl.rtm(background_gas="air", use_numba=True)
    ds = rtm.radiative_transfer(
        NU_MIN, NU_MAX, DNU, p_da, ps, T_da, Ts,
        absorbers={"CO2": CO2_VMR}, line_shape="pseudovoigt",
    )

    # Linepyline's output p coord is in hPa AND may double the level count
    # (staggering NaN-filled half-levels). Reduce to one finite row per input
    # level by dropping NaN rows, then verify and reassign our input p (Pa).
    kappa_raw = ds["kappa"].values
    if ds["kappa"].dims[0] != "p":
        kappa_raw = kappa_raw.T
    good = np.isfinite(kappa_raw).all(axis=1)
    kappa_good = kappa_raw[good]
    if kappa_good.shape[0] != len(p):
        raise SystemExit(
            f"finite kappa rows ({kappa_good.shape[0]}) != input p levels ({len(p)}); "
            "linepyline output structure differs from expectation"
        )
    nu_vals = ds["nu"].values
    kappa_da = xr.DataArray(
        kappa_good, dims=("p", "nu"),
        coords={"p": p, "nu": nu_vals},
    )
    T_da_out = xr.DataArray(T, dims=("p",), coords={"p": p})
    keep = xr.Dataset(
        {
            "kappa": kappa_da,
            "T": T_da_out,
        },
        attrs={
            "profile": args.profile,
            "Ts_K": Ts,
            "CO2_vmr": CO2_VMR,
            "background_gas": "air",
            "line_shape": "pseudovoigt",
            "dnu_cm-1": DNU,
        },
    )
    out = f"examples/lbl_k_{args.profile}.nc"
    keep.to_netcdf(out)
    print(f"wrote {out}  kappa shape={keep['kappa'].shape}")
    print(f"  nu: {NU_MIN}..{NU_MAX} cm-1 dnu={DNU}")
    print(f"  p:  {p[0]:.2f}..{p[-1]:.1f} Pa, {len(p)} levels")
    print(f"  T range: {T.min():.1f}..{T.max():.1f} K")


if __name__ == "__main__":
    main()
