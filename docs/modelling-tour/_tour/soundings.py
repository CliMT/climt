"""Prescribed vertical profiles for the Modelling Tour pages.

Nothing here integrates in time. Every tour page states a profile outright and
asks the radiation code what it makes of it, which keeps each cell to a single
component call and keeps the causal structure visible.

Kept importable and natively testable: the pages import these functions rather
than defining profiles inline, and ``tests/test_modelling_tour.py`` guards them.
"""
import numpy as np

RD = 287.0        # gas constant for dry air, J/(kg K)
G = 9.81          # gravitational acceleration, m/s^2
EPSILON = 0.622   # ratio of molar masses, water vapour to dry air


def saturation_vapour_pressure(T):
    """Saturation vapour pressure over liquid water, in Pa (Bolton 1980)."""
    Tc = np.asarray(T) - 273.15
    return 611.2 * np.exp(17.67 * Tc / (Tc + 243.5))


def lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8, gamma=6.5e-3,
                        T_strat=200.0, q_floor=1e-7):
    """A troposphere at a constant lapse rate under an isothermal stratosphere.

    Hydrostatic balance with a constant lapse rate gives T(p) = T_surf *
    (p/ps)**(RD*gamma/G); the profile is then clipped at ``T_strat``. Humidity
    is set at fixed relative humidity ``rh`` and floored so the stratosphere
    stays inside the k-table's H2O axis.

    Args:
        p: (nz,) mid-level pressure, Pa, surface first.
        ps: surface pressure, Pa.
        T_surf: surface temperature, K.
        rh: relative humidity, dimensionless (0-1).
        gamma: lapse rate, K/m.
        T_strat: isothermal stratosphere temperature, K.
        q_floor: minimum specific humidity, kg/kg.

    Returns:
        (T, q) each (nz,) — air temperature in K, specific humidity in kg/kg.
    """
    p = np.asarray(p, dtype=float)
    T = T_surf * (p / ps) ** (RD * gamma / G)
    T = np.maximum(T, T_strat)
    e = rh * saturation_vapour_pressure(T)
    q = EPSILON * e / np.maximum(p - e, 1.0)
    return T, np.maximum(q, q_floor)


def analytic_gray_equilibrium(p, ps, tau_inf=4.0, T_e=255.0):
    """The course notes' closed-form gray radiative-equilibrium profile.

    Chapter 8 of *Principles of Planetary Climate* derives

        T(tau)  = T_e * [(1 + tau_inf - tau) / 2]**0.25
        T_ground = T_e * (1 + tau_inf / 2)**0.25

    with ``tau`` measured UP from the surface, so ``tau_star = tau_inf - tau``
    is measured down from space. **``tau_inf`` is the diffusivity-SCALED column
    optical depth** — a component evaluating this profile must be built with the
    same ``diffusivity_factor`` the table was calibrated against, or the heating
    rate will not vanish.

    Optical depth is taken linear in pressure (a well-mixed absorber):
    ``tau_star = tau_inf * p / ps``.

    Args:
        p: (nz,) mid-level pressure, Pa, surface first.
        ps: surface pressure, Pa.
        tau_inf: diffusivity-scaled column optical depth.
        T_e: emission temperature, K.

    Returns:
        (T, T_surf, tau_star) — (nz,) air temperature in K, scalar surface
        temperature in K, and (nz,) optical depth measured down from space.
    """
    p = np.asarray(p, dtype=float)
    tau_star = tau_inf * p / ps
    T = T_e * ((1.0 + tau_star) / 2.0) ** 0.25
    T_surf = T_e * (1.0 + tau_inf / 2.0) ** 0.25
    return T, T_surf, tau_star


def apply_sounding(state, T, q=None, T_surf=None):
    """Write a prescribed profile into a climt state, in place.

    Args:
        state: a climt state dict from ``get_default_state``.
        T: (nz,) air temperature, K.
        q: (nz,) specific humidity in kg/kg, or None to leave it alone.
        T_surf: scalar surface temperature in K, or None to leave it alone.
    """
    state["air_temperature"].values[:, 0, 0] = T
    if T_surf is not None:
        state["surface_temperature"].values[:] = T_surf
    if q is not None and "specific_humidity" in state:
        state["specific_humidity"].values[:, 0, 0] = q
