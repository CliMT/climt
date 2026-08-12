"""Native Chaverot R=500 source extraction across a (T, p) envelope.

Companion to chaverot_vs_lbl_envelope_lbl.py + plot_chaverot_vs_lbl_envelope.py.
Reads the native R=500 Earth CO2 k-table via exo_k at its own (T, p) grid nodes
(no interpolation), and writes (a) band-mean k over the three diagnostic 15-um
bands and (b) the per-bin g-integrated opacity spectrum over 480-820 cm^-1, both
converted m^2/molecule -> m^2/kg-air, for a grid of (T, p) nodes.

Run from the radiation env (exo_k + numpy):
    /Users/joymonteiro/miniconda3/envs/radiation/bin/python \
        scripts/experiments/chaverot_vs_lbl_envelope_native.py
"""
from __future__ import annotations

import os

import numpy as np
import exo_k

SRC = ("/Users/joymonteiro/github/chaverot/data/correlated-k_tables/"
       "Earth_376ppmCO2_R500_10-30k.ktable.SI.h5")
OUT = os.path.join(os.path.dirname(__file__), "..", "..",
                   "debug_data", "chaverot_vs_lbl_native.npz")

N_A = 6.02214076e23
M_AIR = 0.028964
CONV = N_A / M_AIR              # m^2/molecule -> m^2/kg-air

T_NODES = [110.0, 170.0, 230.0, 290.0]
P_NODES = [1.0e2, 1.0e3, 1.0e4, 1.0e5]
X_DRY = 1e-6
SPEC_LO, SPEC_HI = 480.0, 820.0
BANDS = {"500-630": (500.0, 630.0),
         "630-700": (630.0, 700.0),
         "700-800": (700.0, 800.0)}


def main():
    kt = exo_k.Ktable5d(filename=SRC)
    wns = np.asarray(kt.wns)
    bw = np.diff(np.asarray(kt.wnedges))
    w = np.asarray(kt.weights)
    Tg = np.asarray(kt.tgrid)
    pg = np.asarray(kt.pgrid)
    xg = np.asarray(kt.xgrid)
    ix = int(np.argmin(np.abs(xg - X_DRY)))

    iTs = [int(np.argmin(np.abs(Tg - T))) for T in T_NODES]
    iPs = [int(np.argmin(np.abs(pg - P))) for P in P_NODES]
    assert all(abs(Tg[i] - T) < 1 for i, T in zip(iTs, T_NODES))
    assert all(abs(pg[i] - P) < 1 for i, P in zip(iPs, P_NODES))

    nT, nP = len(T_NODES), len(P_NODES)
    band_k = {b: np.zeros((nT, nP)) for b in BANDS}

    spec_sel = (wns >= SPEC_LO) & (wns < SPEC_HI)
    wns_spec = wns[spec_sel]
    spec = np.zeros((nT, nP, spec_sel.sum()))

    for a, iT in enumerate(iTs):
        for b, iP in enumerate(iPs):
            # per-bin g-integrated k (m^2/molecule), all bins
            perbin = (kt.kdata[iP, iT, ix, :, :] * w[None, :]).sum(axis=1)
            spec[a, b, :] = perbin[spec_sel] * CONV
            for name, (lo, hi) in BANDS.items():
                sel = (wns >= lo) & (wns < hi)
                band_k[name][a, b] = np.average(perbin[sel],
                                                weights=bw[sel]) * CONV

    np.savez(
        OUT,
        T_nodes=np.array(T_NODES), p_nodes=np.array(P_NODES),
        wns_spec=wns_spec, spec=spec,
        **{f"band_{b.replace('-', '_')}": band_k[b] for b in BANDS},
    )
    print(f"wrote {OUT}")
    print(f"T_nodes={T_NODES} p_nodes={P_NODES} (Pa), x={xg[ix]:.0e}")
    for b in BANDS:
        print(f"  native band-mean k {b} (m^2/kg), rows=T cols=p:")
        print(band_k[b])


if __name__ == "__main__":
    main()
