# scripts/experiments/tropopause.py
"""Robust cross-climate tropopause locator for the picket-fence RCE columns.

Replaces the constant-theta+10K marker (formerly in rce_moist_pf_vs_rrtmg.py)
that mislocated the tropopause for runaway-warm columns -- the 583 hPa artifact.
Pure numpy, so it is Pyodide-safe and importable from the experiment notebook.

Primary criterion: theta-curvature. The tropopause is the kink where potential
temperature transitions from the convective, near-constant-theta troposphere to
the rapidly increasing stratosphere -- the interior level of maximum |d^2 theta /
d(ln p)^2|. Secondary cross-check: the cold-point (temperature minimum), which
rides high (low pressure) in the moist tropics.
"""
import numpy as np

_RD_OVER_CP = 287.0 / 1004.0   # dry-air R/cp, matches the RCE scripts
_P0_HPA = 1000.0


def potential_temperature(temperature, pressure):
    """theta = T (p0/p)^(Rd/cp); pressure in hPa, same shape as temperature."""
    temperature = np.asarray(temperature, dtype=float)
    pressure = np.asarray(pressure, dtype=float)
    return temperature * (_P0_HPA / pressure) ** _RD_OVER_CP


def find_tropopause(temperature, pressure):
    """Locate the tropopause from one column's T(p) profile.

    Parameters
    ----------
    temperature : (nz,) array of K, ordered surface -> top (index 0 = surface).
    pressure : (nz,) array of hPa, same ordering, strictly decreasing.

    Returns
    -------
    dict:
        p_curvature : tropopause pressure (hPa), the headline value (interior
            peak of |d^2 theta / d(ln p)^2|, i.e. max absolute curvature).
        p_coldpoint : pressure (hPa) of the temperature minimum (cross-check).
        k_curvature, k_coldpoint : the corresponding integer level indices.

    Degenerate columns (fewer than 4 levels, or no usable interior curvature)
    return p_curvature = pressure[-1] (model top) as a documented sentinel.
    """
    temperature = np.asarray(temperature, dtype=float)
    pressure = np.asarray(pressure, dtype=float)
    nz = pressure.size

    k_cold = int(np.argmin(temperature))

    if nz < 4:
        return {
            "p_curvature": float(pressure[-1]),
            "p_coldpoint": float(pressure[k_cold]),
            "k_curvature": nz - 1,
            "k_coldpoint": k_cold,
        }

    theta = potential_temperature(temperature, pressure)
    lnp = np.log(pressure)

    # Second derivative on the (non-uniform) ln-p grid; np.gradient twice handles
    # unequal spacing. The stratospheric kink is the interior extremum of largest
    # magnitude; skip the lowest nz//5 levels to avoid boundary-layer kinks, and
    # exclude the two end levels where np.gradient is one-sided. The index-space
    # floor assumes the tropospheric base is near 1000 hPa (true for these RCE
    # columns); shallow or low-pressure-bottom columns may need it adjusted.
    dtheta = np.gradient(theta, lnp)
    d2theta = np.gradient(dtheta, lnp)
    interior = np.arange(max(1, nz // 5), nz - 1)
    k_curv = int(interior[np.argmax(np.abs(d2theta[interior]))])

    return {
        "p_curvature": float(pressure[k_curv]),
        "p_coldpoint": float(pressure[k_cold]),
        "k_curvature": k_curv,
        "k_coldpoint": k_cold,
    }
