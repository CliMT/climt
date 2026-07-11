"""Generate the line-by-line .npz inputs for the linepyline-based figures.

``linepyline`` (Rodrigo Caballero, GPL-3) wraps a HITRAN line-by-line opacity
calculator. It is an optional, environment-specific dependency, so the heavy LBL
computation happens here, once, and is frozen into small ``.npz`` files committed
under ``_artifacts/``. ``_figures.py`` then plots from those arrays in the
``climt`` env (which has matplotlib), so figures 02 and 04 are reproducible in CI
without linepyline.

Run this in an environment that has ``linepyline`` installed::

    python docs/radiative-transfer/_lbl_data.py \
        --out-dir docs/radiative-transfer/_artifacts

Produces:
  02_lbl_H2O_1000_1200.npz   nu, kappa            (H2O, T=296 K, p=1 bar)
  04_correlation_across_T.npz  k_lo, k_hi, T_lo, T_hi  (sorted H2O kappa)
"""
import argparse
import os

import numpy as np


def _kappa(r, nu_min, nu_max, dnu, p, T):
    """H2O mass absorption coefficient (m^2/kg) on the wavenumber grid."""
    k = r.get_kappa_hitran("H2O", nu_min, nu_max, dnu, p=p, T=T)
    return np.asarray(k)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True, help="directory for the .npz files")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    import linepyline as lpl

    r = lpl.rtm()

    # Figure 2.1: full-resolution H2O spectrum, 1000-1200 cm^-1, T=296 K, p=1 bar.
    nu = np.asarray(r.get_nu_grid(1000.0, 1200.0, 0.01))
    kappa = _kappa(r, 1000.0, 1200.0, 0.01, p=1.0e5, T=296.0)
    out02 = os.path.join(args.out_dir, "02_lbl_H2O_1000_1200.npz")
    np.savez_compressed(out02, nu=nu, kappa=kappa, T=296.0, p=1.0e5)
    print(f"wrote {out02}  ({nu.size} points)")

    # Figure 4.1: same band at two temperatures, sorted, to show the
    # correlated-k rank-order assumption. Coarser grid keeps the file small.
    k_lo = np.sort(_kappa(r, 1000.0, 1200.0, 0.05, p=1.0e5, T=250.0))
    k_hi = np.sort(_kappa(r, 1000.0, 1200.0, 0.05, p=1.0e5, T=350.0))
    out04 = os.path.join(args.out_dir, "04_correlation_across_T.npz")
    np.savez_compressed(out04, k_lo=k_lo, k_hi=k_hi, T_lo=250.0, T_hi=350.0)
    print(f"wrote {out04}  ({k_lo.size} points)")


if __name__ == "__main__":
    main()
