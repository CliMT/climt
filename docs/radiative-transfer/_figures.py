"""Walkthrough figure functions for the radiative-transfer chapters.

Ported from the phase-4 plan's ``_generate_figures.py`` and adapted to the
CORK module layout (``climt._components.cork.*`` — the scheme formerly named
"picket fence"). Driven by ``sources.yml`` through ``scripts/build_experiments.py``
(``make experiments``): each artifact invokes this module with
``--figure NAME --out PATH`` and one PNG is written per call.

Two figures (``02_lbl_H2O_1000_1200`` and ``04_correlation_across_T``) are
line-by-line and depend on the optional ``linepyline`` package. To keep them
reproducible here and in CI without that dependency, the heavy LBL arrays are
frozen into ``.npz`` files (committed beside the PNGs) by ``_lbl_data.py``, which
is run once in a ``linepyline`` env. These two functions just plot from the
frozen ``.npz``; regenerate the inputs with::

    python docs/radiative-transfer/_lbl_data.py \
        --out-dir docs/radiative-transfer/_artifacts
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def _save(fig, out):
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


def fig_01_mean_of_exp(out):
    """Chapter 1: mean-of-exponential vs exponential-of-mean."""
    sigma = np.concatenate([0.01 * np.ones(90), 5.0 * np.ones(10)])
    L = np.linspace(0, 1, 200)
    T_true = np.mean(np.exp(-sigma[:, None] * L[None, :]), axis=0)
    T_naive = np.exp(-np.mean(sigma) * L)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(L, T_true, "k-", label=r"$\langle e^{-\sigma L}\rangle$ (true)")
    ax.plot(L, T_naive, "r--", label=r"$e^{-\langle\sigma\rangle L}$ (naive)")
    ax.set_xlabel("path length L")
    ax.set_ylabel("transmission")
    ax.set_title("Mean of an exponential $\\neq$ exponential of the mean")
    ax.legend()
    _save(fig, out)


def _load_npz(out):
    """Load the frozen LBL .npz that sits beside the target PNG."""
    npz = os.path.splitext(out)[0] + ".npz"
    if not os.path.exists(npz):
        raise FileNotFoundError(
            f"{npz} not found. Generate it once in a linepyline env:\n"
            f"  python docs/radiative-transfer/_lbl_data.py "
            f"--out-dir {os.path.dirname(out)}")
    return np.load(npz)


def fig_02_lbl_spectrum(out):
    """Chapter 2: H2O LBL absorption near 1000-1200 cm^-1 (from frozen .npz)."""
    d = _load_npz(out)
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.semilogy(d["nu"], d["kappa"], lw=0.4)
    ax.set_xlabel("wavenumber (cm$^{-1}$)")
    ax.set_ylabel(r"$\kappa$ (m$^{2}$/kg)")
    ax.set_title("H$_2$O LBL absorption, T=296 K, p=1 bar")
    _save(fig, out)


def fig_03_k_distribution_construction(out):
    """Chapter 3: build a CDF of σ and drop 2-point Gauss-Legendre markers."""
    rng = np.random.default_rng(0)
    # Toy band: mix of weak and strong lines
    sigma = np.concatenate([rng.uniform(1e-24, 1e-22, 800),
                            rng.uniform(1e-20, 1e-18, 200)])
    sigma_sorted = np.sort(sigma)
    g = np.linspace(0, 1, sigma.size)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.semilogy(g, sigma_sorted, "k-")
    for gi, wi in zip([0.211, 0.789], [0.5, 0.5]):
        ki = np.interp(gi, g, sigma_sorted)
        ax.plot(gi, ki, "ro")
        ax.annotate(f"g={gi:.2f}, w={wi:.2f}", (gi, ki),
                    textcoords="offset points", xytext=(6, 4))
    ax.set_xlabel("cumulative fraction g")
    ax.set_ylabel(r"k(g), cm$^{2}$/molecule")
    ax.set_title("k-distribution and 2-point Gauss quadrature")
    _save(fig, out)


def fig_04_correlation_across_T(out):
    """Chapter 4: sorted κ at two different T share rank order (from frozen .npz)."""
    d = _load_npz(out)
    s1, s2 = d["k_lo"], d["k_hi"]
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.loglog(s1, s2, ".", ms=1)
    ax.plot([s1.min(), s1.max()], [s1.min(), s1.max()], "k--", alpha=0.4)
    ax.set_xlabel(r"sorted $\kappa$(T=250 K), m$^2$/kg")
    ax.set_ylabel(r"sorted $\kappa$(T=350 K), m$^2$/kg")
    ax.set_title("Correlated-k assumption: rank order is approximately preserved")
    _save(fig, out)


def fig_06_picket_fence_opacity(out):
    """Chapter 6: Parmentier LW opacity structure for a hot-Jupiter column."""
    from climt._components.cork.optics.parmentier import (
        compute_rosseland_mean_opacity,
        load_freedman2014_coefficients,
        load_parmentier_coefficients,
        lookup_ratio_coefficients,
    )
    coeffs = load_freedman2014_coefficients()
    parm = load_parmentier_coefficients("solar_composition")
    T_eff = 1450.0
    lookup_ratio_coefficients(parm, T_eff)
    p = np.geomspace(1.0, 1e7, 50)
    T = 1000.0 * np.ones_like(p)
    kR = np.array([compute_rosseland_mean_opacity(t, pi, coeffs)
                   for t, pi in zip(T, p)])
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.loglog(kR, p, "k-", label=r"$\kappa_R$")
    ax.invert_yaxis()
    ax.set_xlabel("opacity (m$^2$/kg)")
    ax.set_ylabel("pressure (Pa)")
    ax.set_title("Rosseland mean opacity (Freedman 2014)")
    ax.legend()
    _save(fig, out)


def fig_07_two_stream_phases(out):
    """Chapter 7: diffuse R and T as functions of optical depth for ssa=1, g=0.5."""
    from climt._components.cork.sw.kernels import _sw_dif_and_source
    taus = np.geomspace(1e-3, 1e2, 40)
    Rdif, Tdif = [], []
    for t in taus:
        rdif, tdif, _, _, _ = _sw_dif_and_source(t, 1.0, 0.5, 0.7)
        Rdif.append(rdif)
        Tdif.append(tdif)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.semilogx(taus, Rdif, label="$R_{dif}$")
    ax.semilogx(taus, Tdif, label="$T_{dif}$")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel("amplitude")
    ax.set_title("Two-stream diffuse R, T (ssa=1, g=0.5, $\\mu_0$=0.7)")
    ax.legend()
    _save(fig, out)


FIGURES = {
    "01_mean_of_exp": fig_01_mean_of_exp,
    "02_lbl_H2O_1000_1200": fig_02_lbl_spectrum,
    "03_k_distribution_construction": fig_03_k_distribution_construction,
    "04_correlation_across_T": fig_04_correlation_across_T,
    "06_picket_fence_opacity": fig_06_picket_fence_opacity,
    "07_two_stream_phases": fig_07_two_stream_phases,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figure", required=True, choices=sorted(FIGURES),
                    help="which walkthrough figure to generate")
    ap.add_argument("--out", required=True, help="output PNG path")
    args = ap.parse_args()
    FIGURES[args.figure](args.out)


if __name__ == "__main__":
    main()
