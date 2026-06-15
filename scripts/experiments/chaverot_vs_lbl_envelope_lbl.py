"""LBL (linepyline) extraction across the same (T, p) envelope.

Companion to chaverot_vs_lbl_envelope_native.py. Computes line-by-line CO2 kappa
(pseudovoigt, binning=False, dnu=0.1, MT_CKD-free dry CO2) at the same (T, p)
nodes, writing band-mean k over the three diagnostic 15-um bands and the fine
kappa(nu) spectrum over 480-820 cm^-1.

Run from the linepyline env:
    /Users/joymonteiro/miniconda3/envs/linepyline/bin/python \
        scripts/experiments/chaverot_vs_lbl_envelope_lbl.py
"""
from __future__ import annotations

import os

import numpy as np
import xarray as xr

import linepyline as lpl

OUT = os.path.join(os.path.dirname(__file__), "..", "..",
                   "debug_data", "chaverot_vs_lbl_lbl.npz")

CO2_VMR = 376e-6
T_NODES = [110.0, 170.0, 230.0, 290.0]
P_NODES = [1.0e2, 1.0e3, 1.0e4, 1.0e5]   # ascending, used as the column profile
NU_MIN, NU_MAX, DNU = 300.0, 1000.0, 0.1  # wide enough for 25 cm^-1 line cutoff
SPEC_LO, SPEC_HI = 480.0, 820.0
BANDS = {"500-630": (500.0, 630.0),
         "630-700": (630.0, 700.0),
         "700-800": (700.0, 800.0)}


def main():
    rtm = lpl.rtm(background_gas="air", use_numba=True)
    p_prof = np.array(P_NODES)
    nT, nP = len(T_NODES), len(P_NODES)

    band_k = {b: np.zeros((nT, nP)) for b in BANDS}
    nu_ref = None
    spec = None

    for a, T in enumerate(T_NODES):
        p_da = xr.DataArray(p_prof, dims=("p",), coords={"p": p_prof})
        T_da = xr.DataArray(np.full_like(p_prof, T), dims=("p",),
                            coords={"p": p_prof})
        ds = rtm.radiative_transfer(
            NU_MIN, NU_MAX, DNU, p_da, float(p_prof[-1]), T_da, float(T),
            absorbers={"CO2": CO2_VMR}, line_shape="pseudovoigt",
        )
        kappa = ds["kappa"].values
        if ds["kappa"].dims[0] != "p":
            kappa = kappa.T                       # -> (p, nu)
        nu = ds["nu"].values
        if nu_ref is None:
            nu_ref = nu
            spec_sel = (nu >= SPEC_LO) & (nu < SPEC_HI)
            spec = np.zeros((nT, nP, int(spec_sel.sum())))
        for b in range(nP):
            row = np.maximum(kappa[b], 0.0)
            spec[a, b, :] = row[spec_sel]
            for name, (lo, hi) in BANDS.items():
                sel = (nu >= lo) & (nu < hi)
                band_k[name][a, b] = float(np.mean(row[sel]))

    np.savez(
        OUT,
        T_nodes=np.array(T_NODES), p_nodes=np.array(P_NODES),
        wns_spec=nu_ref[(nu_ref >= SPEC_LO) & (nu_ref < SPEC_HI)], spec=spec,
        **{f"band_{b.replace('-', '_')}": band_k[b] for b in BANDS},
    )
    print(f"wrote {OUT}")
    for b in BANDS:
        print(f"  LBL band-mean k {b} (m^2/kg), rows=T cols=p:")
        print(band_k[b])


if __name__ == "__main__":
    main()
