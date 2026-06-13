"""Convert a Chaverot (Exo_k-format) correlated-k table into climt's picket-fence
netCDF format and re-bin to picket-fence resolution.

Usage:
    python scripts/generate_picket_fence_tables.py \\
        --input /path/to/chaverot_table.h5 \\
        --output climt/_data/picket_fence/correlated_k/earth_low_res_lw.nc \\
        --kind lw \\
        --bands 10,500,1250,2500,3250 \\
        --ngpt 2 \\
        --mixture-molar-mass 0.028964

WHY WE KEEP THE H2O VMR AXIS
----------------------------
Chaverot tables are ``Ktable5d`` objects: (P, T, X_H2O, ν, g). The X axis is
H2O volume mixing ratio — it exists specifically so H2O can be a *runtime*
variable. An earlier version of this script collapsed that axis at
table-build time by log-interpolating to a single representative VMR (e.g.
0.01 for Earth). That was wrong for any realistic simulation: RRTMG (and any
real model) drives H2O from a dynamic ``specific_humidity`` field, but a
collapsed table is frozen at whatever fraction the generator picked. Headline
fluxes (OLR, surface SW) can *accidentally* match RRTMG across a range of
humidity profiles if the bands chosen are CO2-dominated, but the table
stops responding to moisture entirely — it cannot distinguish a dry desert
from a moist tropical column.

This version preserves the H2O axis in the output netCDF. The component
reads ``specific_humidity`` from state, converts to mole fraction, and does
trilinear interpolation in (T, log P, log X_H2O) at runtime. The "gas" in
``gas_names = ("h2o",)`` is therefore honest: H2O really is the only
runtime-variable constituent. Everything else (CO2, N2, continua) stays
baked into the k-coefficients at the composition Chaverot built them for.

The netCDF attribute ``background_is_premixed = "true"`` tells climt's
picket-fence component that this is a bulk-air table: per-layer optical
depth is computed as ``tau = k(T, p, X_H2O) * column_mass_of_AIR`` (not
``k * column_mass_of_H2O``). The k is already m^2/kg-of-mixture, so it
scales with air mass, not H2O mass.

OTHER NOTES
-----------
    * Requires ``exo_k`` (GPL-3), not a climt runtime dep. Run from a dev env.
    * Chaverot tables (Zenodo 16795590, CC-BY 4.0) cover the planets / mixtures
      listed in ``regenerate_shipped_picket_fence_tables.py``.
    * ``bands`` argument is N+1 band edges in wavenumber (cm^-1). For LW on
      Earth, 10, 500, 1250, 2500, 3250 cm^-1 gives rotational / window /
      CO2 / H2O-vib splits.
    * Chaverot tables store k in m^2/molecule; climt expects m^2/kg-of-mixture.
      ``--mixture-molar-mass`` (kg/mol) is used for the N_A/M conversion.
      Typical values: Earth (dry-air) 0.028964, Mars (CO2-dominant) 0.04355,
      Venus (CO2-dominant) 0.04344.
    * ``planck_fraction`` is written as the true per-band Planck integral,
      B_band(T)/(sigma T^4), evaluated at each temperature grid point. It is
      replicated across g-points (correlated-k assumes uniform distribution
      of Planck emission within a band).
"""
import argparse
import os

import numpy as np

# Avogadro constant (1/mol)
_N_AVOGADRO = 6.02214076e23
# Planck's constant (J s), speed of light (m/s), Boltzmann constant (J/K),
# Stefan-Boltzmann (W/m^2/K^4).
_H_PLANCK = 6.62607015e-34
_C_LIGHT = 2.99792458e8
_K_BOLTZ = 1.380649e-23
_SIGMA_SB = 5.670374419e-8


def _planck_B_wavenumber(wn_cm, T):
    """Planck function B_nu evaluated per wavenumber (cm^-1), units W/m^2/sr/(cm^-1).

    Uses B_wn = 2 h c^2 nu^3 / (exp(h c nu / kT) - 1), where nu is wavenumber
    in cm^-1 converted to m^-1 via factor 100.
    """
    nu_m = 100.0 * np.asarray(wn_cm, dtype=np.float64)  # m^-1
    x = _H_PLANCK * _C_LIGHT * nu_m / (_K_BOLTZ * T)
    with np.errstate(over="ignore"):
        denom = np.expm1(np.clip(x, 0.0, 700.0))
    denom = np.where(denom == 0, np.inf, denom)
    B = 2.0 * _H_PLANCK * _C_LIGHT ** 2 * nu_m ** 3 / denom
    # Convert per m^-1 to per cm^-1 (differential) by multiplying by 100:
    return B * 100.0


def _band_planck_fraction(wn_lo_cm, wn_hi_cm, T, n_quad=64):
    """Fraction of total Planck emission pi*integral(B_wn) dwn in [wn_lo, wn_hi]
    relative to sigma*T^4. Midpoint rule on a linear grid (sufficient for the
    smooth Planck spectrum over typical ~500 cm^-1 bands).
    """
    grid = np.linspace(wn_lo_cm, wn_hi_cm, n_quad + 1)
    mid = 0.5 * (grid[1:] + grid[:-1])
    dwn = np.diff(grid)
    B = _planck_B_wavenumber(mid, T)
    band_flux = np.pi * np.sum(B * dwn)
    return band_flux / (_SIGMA_SB * T ** 4)


def _solar_band_irradiance(band_edges_wn, S0=1361.0, T_sun=5778.0, n_quad=512):
    """Per-band solar irradiance (W/m²) from a T_sun blackbody scaled to S0.

    Band fractions are computed from the full blackbody integral (10–100000
    cm⁻¹) and then rescaled so they sum exactly to S0.  This ensures the
    total solar constant is conserved regardless of which spectral range the
    bands cover, while preserving the correct relative distribution across
    bands (near-IR vs visible vs UV).

    Values are replicated across g-points by the caller; the component weights
    them by g-point weights that sum to 1, so the per-g-point value equals
    the full band solar flux.
    """
    wn_ref = np.linspace(10.0, 100000.0, 100 * n_quad)
    B_ref = np.pi * _planck_B_wavenumber(wn_ref, T_sun)
    total = np.trapz(B_ref, wn_ref)

    fractions = []
    for lo, hi in zip(band_edges_wn[:-1], band_edges_wn[1:]):
        mask = (wn_ref >= lo) & (wn_ref <= hi)
        fractions.append(np.trapz(B_ref[mask], wn_ref[mask]) / total)

    fractions = np.array(fractions)
    # Renormalize so bands sum to S0 — conserves solar constant.
    return S0 * fractions / fractions.sum()


# Rayleigh cross-section reference: Bodhaine et al. (1999) at 0.55 µm.
_SIGMA_RAYLEIGH_REF = 4.95e-31   # m²/molecule at ν_ref
_WN_RAYLEIGH_REF = 18182.0        # cm⁻¹ (0.55 µm)


def _rayleigh_per_band(band_edges_wn, solar_per_band, molar_mass=0.028964,
                       n_quad=512):
    """Solar-spectrum-weighted Rayleigh coefficient per band (m²/kg-of-air).

    σ_R(ν) = σ_ref × (ν/ν_ref)⁴  (Bodhaine 1999).
    Band average is weighted by the solar irradiance spectrum so that the
    scattering optical depth matches the irradiance-weighted band opacity.
    Returned coefficient r satisfies τ_R = r × (dp/g).
    """
    wn_fine = np.linspace(10.0, 100000.0, 100 * n_quad)
    sigma_fine = _SIGMA_RAYLEIGH_REF * (wn_fine / _WN_RAYLEIGH_REF) ** 4

    # solar-spectrum weight for weighting within each band
    B_fine = np.pi * _planck_B_wavenumber(wn_fine, 5778.0)

    out = []
    for lo, hi in zip(band_edges_wn[:-1], band_edges_wn[1:]):
        mask = (wn_fine >= lo) & (wn_fine <= hi)
        wts = B_fine[mask]
        sigma_avg = np.trapz(sigma_fine[mask] * wts, wn_fine[mask]) / np.trapz(wts, wn_fine[mask])
        # Convert m²/molecule → m²/kg-of-air
        out.append(sigma_avg * _N_AVOGADRO / molar_mass)
    return np.array(out)


def _two_stretch_gl_nodes(ngpt, g_split):
    """Two-stretch Gauss-Legendre nodes on [0, 1].

    Places ngpt//2 GL nodes on [0, g_split] and ngpt - ngpt//2 on [g_split, 1].
    The upper stretch densely samples the strong-line tail (g→1) that plain
    Gauss-Legendre under-resolves for molecular bands with steep line cores.

    Weights sum to 1.0 over the combined node set.
    """
    n_lo = ngpt // 2
    n_hi = ngpt - n_lo
    xi_lo, wi_lo = np.polynomial.legendre.leggauss(n_lo)
    xi_hi, wi_hi = np.polynomial.legendre.leggauss(n_hi)
    # Map [-1, 1] -> [0, g_split] and [g_split, 1]
    g_lo = 0.5 * (xi_lo + 1.0) * g_split
    g_hi = g_split + 0.5 * (xi_hi + 1.0) * (1.0 - g_split)
    w_lo = 0.5 * wi_lo * g_split
    w_hi = 0.5 * wi_hi * (1.0 - g_split)
    g_nodes = np.concatenate([g_lo, g_hi])
    weights = np.concatenate([w_lo, w_hi])
    return g_nodes, weights


def _rebin_to_bands(ktable, band_edges_wn, ngpt,
                    quadrature="gauss-legendre", g_split=0.95):
    """Re-bin an Exo_k Ktable/Ktable5d to new bands with ``ngpt`` g-points per band.

    For a Ktable5d (Chaverot), the H2O VMR (X) axis is preserved end-to-end
    so the component can interpolate in log(X) at runtime. See the module
    docstring for why we no longer collapse it.
    """
    k_binned = ktable.copy()
    k_binned.bin_down(np.asarray(band_edges_wn))

    nT = k_binned.Nt
    nP = k_binned.Np
    nband = k_binned.Nw
    original_g = k_binned.ggrid

    is_5d = hasattr(k_binned, 'Nx')
    if is_5d:
        nX = k_binned.Nx
        kdata = k_binned.kdata  # (Np, Nt, Nx, Nw, Ng)
        x_grid = np.asarray(k_binned.xgrid)
    else:
        nX = 1
        kdata = k_binned.kdata[:, :, np.newaxis, :, :]
        x_grid = np.asarray([1.0])  # dummy, not written for non-H2O tables

    if quadrature == "gauss-legendre":
        xi, wi = np.polynomial.legendre.leggauss(ngpt)
        xi = 0.5 * (xi + 1.0)
        wi = 0.5 * wi
    elif quadrature == "two-stretch":
        if ngpt < 2:
            raise ValueError("two-stretch quadrature requires ngpt >= 2")
        xi, wi = _two_stretch_gl_nodes(ngpt, g_split)
    else:
        raise ValueError(f"unknown quadrature: {quadrature!r}")

    # Output layout: (nband, ngpt, T, P, X). Transposed from Chaverot's
    # (P, T, X, …) layout to match climt's column-oriented k-lookup.
    k_out = np.zeros((nband, ngpt, nT, nP, nX))
    for b in range(nband):
        for t in range(nT):
            for p in range(nP):
                for ix in range(nX):
                    kline = kdata[p, t, ix, b, :]
                    k_out[b, :, t, p, ix] = np.interp(xi, original_g, kline)

    return {
        "k": k_out,  # (nband, ngpt, T, P, X)
        "g_weights": np.broadcast_to(wi, (nband, ngpt)).copy(),
        "T_grid": np.asarray(k_binned.tgrid),
        "log_p_grid": np.log(np.asarray(k_binned.pgrid)),
        "h2o_vmr_grid": x_grid,
        "has_h2o_axis": is_5d,
        "wn_edges": np.asarray(band_edges_wn),
    }


def _write_climt_netcdf(out_path, rebinned, kind, gas_names=("h2o",),
                        rayleigh=None, solar_per_gpt=None):
    from scipy.io import netcdf_file

    nband, ngpt, nT, nP, nX = rebinned["k"].shape
    T_grid = rebinned["T_grid"]
    edges = rebinned["wn_edges"]
    has_h2o_axis = rebinned["has_h2o_axis"]

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with netcdf_file(out_path, "w") as nc:
        nc.createDimension("band", nband)
        nc.createDimension("gpoint", ngpt)
        nc.createDimension("temperature", nT)
        nc.createDimension("pressure", nP)
        nc.createDimension("gas", len(gas_names))
        nc.createDimension("bounds", 2)
        if has_h2o_axis:
            nc.createDimension("h2o_vmr", nX)

        k_dims = ("gas", "band", "gpoint", "temperature", "pressure")
        if has_h2o_axis:
            k_dims = k_dims + ("h2o_vmr",)
            k_arr = rebinned["k"][np.newaxis]  # (1, nband, ngpt, T, P, X)
        else:
            k_arr = rebinned["k"][np.newaxis, :, :, :, :, 0]  # drop trivial X
        k = nc.createVariable("k_coefficients", "f4", k_dims)
        k[:] = k_arr

        w = nc.createVariable("gpoint_weights", "f4", ("band", "gpoint"))
        w[:] = rebinned["g_weights"]

        tg = nc.createVariable("temperature_grid", "f4", ("temperature",))
        tg[:] = T_grid

        pg = nc.createVariable("pressure_grid_log", "f4", ("pressure",))
        pg[:] = rebinned["log_p_grid"]

        if has_h2o_axis:
            xg = nc.createVariable("h2o_vmr_grid", "f4", ("h2o_vmr",))
            xg[:] = rebinned["h2o_vmr_grid"]

        limits = np.column_stack([edges[:-1], edges[1:]])
        bl = nc.createVariable("band_wavenumber_limits", "f4", ("band", "bounds"))
        bl[:] = limits

        if kind == "lw":
            pf = nc.createVariable("planck_fraction", "f4",
                                   ("band", "gpoint", "temperature"))
            pf_arr = np.zeros((nband, ngpt, nT), dtype=np.float32)
            for ib in range(nband):
                for it, T in enumerate(T_grid):
                    f_band = _band_planck_fraction(
                        float(edges[ib]), float(edges[ib + 1]), float(T)
                    )
                    pf_arr[ib, :, it] = f_band
            pf[:] = pf_arr
        else:
            sf = nc.createVariable("solar_source_per_gpoint", "f4",
                                   ("band", "gpoint"))
            # solar_per_gpt shape: (nband,) — same value replicated across
            # g-points. The component weights by g-point weights that sum to 1,
            # so the full band solar flux equals solar_per_gpt[b].
            if solar_per_gpt is not None:
                solar_arr = np.asarray(solar_per_gpt)
            else:
                raise ValueError(
                    "solar_per_gpt must be supplied for SW tables; "
                    "compute via _solar_band_irradiance()"
                )
            sf[:] = solar_arr[:, np.newaxis]  # broadcast to (band, gpoint)
            r = nc.createVariable("rayleigh_coefficient", "f4", ("band",))
            r[:] = rayleigh if rayleigh is not None else 0.0

        nc.gas_names = ",".join(gas_names)
        nc.overlap_method = "additive"
        nc.resolution = "low"
        nc.source = "chaverot_zenodo_16795590_rebinned"
        # Flag the component: the non-H2O background (CO2, N2, continua) is
        # pre-mixed into k, so optical depth scales with column mass of AIR,
        # not of H2O. H2O VMR selects the X-axis slice via live interpolation.
        if has_h2o_axis:
            nc.background_is_premixed = "true"


def main():
    import exo_k

    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Chaverot Exo_k table (.h5 or .nc)")
    ap.add_argument("--output", required=True, help="climt .nc path")
    ap.add_argument("--kind", choices=("lw", "sw"), required=True)
    ap.add_argument("--bands", required=True,
                    help="comma-separated band edges in wavenumber (cm^-1)")
    ap.add_argument("--ngpt", type=int, default=2)
    ap.add_argument("--quadrature", default="gauss-legendre",
                    choices=("gauss-legendre", "two-stretch"),
                    help="g-node placement. 'gauss-legendre' (default) is the "
                         "shipped behavior. 'two-stretch' splits [0,1] at "
                         "--g-split and applies GL to each subinterval, "
                         "densely sampling the strong-line tail g->1.")
    ap.add_argument("--g-split", type=float, default=0.95,
                    help="Split point for two-stretch quadrature (default 0.95). "
                         "Ignored for gauss-legendre.")
    ap.add_argument("--mixture-molar-mass", type=float, default=0.028964,
                    help="Molar mass of the gas mixture in kg/mol (used to "
                         "convert Chaverot k from m^2/molecule to m^2/kg). "
                         "Earth (dry air) 0.028964, Mars 0.04355, Venus 0.04344.")
    args = ap.parse_args()

    edges = [float(s) for s in args.bands.split(",")]
    try:
        ktable = exo_k.Ktable(filename=args.input)
    except ValueError:
        ktable = exo_k.Ktable5d(filename=args.input)

    rebinned = _rebin_to_bands(ktable, edges, args.ngpt,
                               quadrature=args.quadrature, g_split=args.g_split)

    # Convert k from m^2/molecule to m^2/kg-of-mixture using Avogadro / molar mass.
    conversion = _N_AVOGADRO / args.mixture_molar_mass
    rebinned["k"] = rebinned["k"] * conversion

    solar_per_gpt = rayleigh = None
    if args.kind == "sw":
        solar_per_gpt = _solar_band_irradiance(edges)
        rayleigh = _rayleigh_per_band(edges, solar_per_gpt,
                                      molar_mass=args.mixture_molar_mass)
        print(f"  solar per band (W/m²): {np.round(solar_per_gpt, 2)}")
        print(f"  rayleigh coeff (m²/kg): {rayleigh}")

    _write_climt_netcdf(args.output, rebinned, args.kind,
                        solar_per_gpt=solar_per_gpt, rayleigh=rayleigh)
    print(f"wrote {args.output}  shape={rebinned['k'].shape} "
          f"has_h2o_axis={rebinned['has_h2o_axis']} "
          f"(k x{conversion:.3e} m^2/kg per m^2/molec)")


if __name__ == "__main__":
    main()
