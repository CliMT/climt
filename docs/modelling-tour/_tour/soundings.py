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


def saturation_specific_humidity(T, p):
    """Saturation specific humidity, kg/kg, at temperature ``T`` (K) and
    pressure ``p`` (Pa).

    Written the way ``climt.GridScaleCondensation`` writes it -- Bolton's
    vapour pressure, ``epsilon = Rd / Rv`` from sympl's constants, and
    ``q = epsilon e / (p - (1 - epsilon) e)`` -- so that "relative humidity
    1" here is exactly the threshold at which that component starts to
    condense. (``lapse_rate_sounding`` above uses the rounder 0.622 and the
    mixing-ratio form; the two differ by well under 1 % in the troposphere,
    which is harmless for a prescribed profile but not for a test asserting
    that the condensation holds relative humidity at 1.)
    """
    from sympl import get_constant

    epsilon = (float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
               / float(get_constant("gas_constant_of_vapor_phase",
                                    "J/kg/degK")))
    e_sat = saturation_vapour_pressure(T)
    return epsilon * e_sat / np.maximum(
        np.asarray(p, dtype=float) - (1.0 - epsilon) * e_sat, 1.0)


def relative_humidity(state):
    """A column's relative humidity profile, (nz,) as a fraction, bottom
    first -- ``specific_humidity`` over :func:`saturation_specific_humidity`.
    """
    def column(name):
        values = np.asarray(state[name].values, dtype=float)
        return values.reshape(values.shape[0], -1)[:, 0]

    return (column("specific_humidity")
            / saturation_specific_humidity(column("air_temperature"),
                                           column("air_pressure")))


def _moist_constants():
    from sympl import get_constant

    def constant(name, units):
        return float(get_constant(name, units))

    Rd = constant("gas_constant_of_dry_air", "J/kg/degK")
    return dict(
        Rd=Rd,
        epsilon=Rd / constant("gas_constant_of_vapor_phase", "J/kg/degK"),
        Cp=constant("heat_capacity_of_dry_air_at_constant_pressure",
                    "J/kg/degK"),
        g=constant("gravitational_acceleration", "m/s^2"),
        Lv=constant("latent_heat_of_condensation", "J/kg"),
    )


def moist_adiabatic_lapse_rate(T, p):
    """How fast a saturated parcel cools as it rises, K/km, at ``T`` (K) and
    ``p`` (Pa) -- the pseudo-adiabatic lapse rate.

        Gamma_m = g (1 + Lv r / (Rd T)) / (cp + Lv^2 r epsilon / (Rd T^2))

    with ``r`` the saturation mixing ratio. The numerator's extra term is the
    parcel's lower density; the denominator's is the latent heat released as
    it condenses, and that is the one that matters: warm air holds a lot of
    vapour and releases a lot of heat per kelvin of cooling, so its rate is
    far below the dry g/cp. Cold air holds little, and its rate tends back to
    g/cp. Same constants, and same q_sat, as ``GridScaleCondensation``.
    """
    return _moist_lapse(np.asarray(T, dtype=float), p, _moist_constants())


def _moist_lapse(T, p, c):
    """``moist_adiabatic_lapse_rate`` with the constants already looked up
    (sympl's ``get_constant`` is too slow to call inside an integration)."""
    e_sat = saturation_vapour_pressure(T)
    q = c["epsilon"] * e_sat / np.maximum(
        np.asarray(p, dtype=float) - (1.0 - c["epsilon"]) * e_sat, 1.0)
    r = q / (1.0 - q)
    return 1e3 * c["g"] * (1.0 + c["Lv"] * r / (c["Rd"] * T)) / (
        c["Cp"] + c["Lv"] ** 2 * r * c["epsilon"] / (c["Rd"] * T ** 2))


def moist_adiabat(T_start, p_start, p, substeps=40):
    """The temperature of a saturated parcel lifted from ``(T_start,
    p_start)``, at each pressure in ``p`` (Pa) -- a moist adiabat.

    Integrates ``dT/dln p = Gamma_m R_d T / g`` upward, ``substeps`` steps
    per level. Levels below ``p_start`` come back as NaN: the parcel starts
    where it starts.
    """
    c = _moist_constants()
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    T, p_now = float(T_start), float(p_start)
    for k in np.argsort(-p):                       # bottom up
        if p[k] > p_start:
            continue
        step = np.log(p[k] / p_now) / substeps     # negative: going up
        for _ in range(substeps):
            rate = _moist_lapse(T, p_now, c) / 1e3             # K/m
            T += rate * c["Rd"] * T / c["g"] * step
            p_now *= np.exp(step)
        out[k] = T
    return out


def lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8, gamma=6.5e-3,
                        T_strat=200.0, q_floor=1e-7, gamma_strat=0.0):
    """A troposphere at a constant lapse rate under a settable stratosphere.

    Hydrostatic balance with a constant lapse rate gives T(p) = T_surf *
    (p/ps)**(RD*gamma/G); the profile is then clipped at ``T_strat``. Humidity
    is set at fixed relative humidity ``rh`` and floored so the stratosphere
    stays inside the k-table's H2O axis.

    ``gamma_strat`` gives the stratosphere a temperature gradient instead of
    leaving it isothermal: positive values WARM with height, as Earth's real
    ozone-heated stratosphere does. Page 5 needs this, because the 15 um core's
    forcing is set by the temperature structure its emission level sits in, and
    an isothermal cap is the one structure in which that forcing has nowhere to
    come from. Height above the tropopause is taken as ``H * log(p_trop / p)``
    with the isothermal scale height ``H = RD * T_strat / G``, which is exact
    for the isothermal case and a good approximation for the small gradients
    this is used with.

    Args:
        p: (nz,) mid-level pressure, Pa, surface first.
        ps: surface pressure, Pa.
        T_surf: surface temperature, K.
        rh: relative humidity, dimensionless (0-1).
        gamma: lapse rate, K/m.
        T_strat: stratosphere temperature at the tropopause, K.
        q_floor: minimum specific humidity, kg/kg.
        gamma_strat: stratospheric temperature gradient, K/m, positive upward.
            0.0 (the default) leaves the stratosphere isothermal.

    Returns:
        (T, q) each (nz,) — air temperature in K, specific humidity in kg/kg.
    """
    p = np.asarray(p, dtype=float)
    T = T_surf * (p / ps) ** (RD * gamma / G)
    T = np.maximum(T, T_strat)

    if gamma_strat:
        # The tropopause is where the tropospheric profile first reaches T_strat.
        p_trop = ps * (T_strat / T_surf) ** (G / (RD * gamma))
        scale_height = RD * T_strat / G
        above = p < p_trop
        z = scale_height * np.log(p_trop / p[above])
        T[above] = T_strat + gamma_strat * z

    e = rh * saturation_vapour_pressure(T)
    # Vapour pressure cannot be a large fraction of the total pressure. Without
    # this the fixed-RH formula returns nonsense wherever the air is warm and
    # the pressure low -- a stratosphere warmed by gamma_strat reaches ~273 K at
    # 1 hPa, where saturation alone is 600 Pa against a total of 100, and q came
    # out at 300 kg/kg. Clipping at half the total pressure bounds q by
    # EPSILON = 0.622 and is a no-op everywhere else: the isothermal 200 K cap
    # has e ~ 0.2 Pa, and the 288 K surface ~1400 Pa against 1000 hPa.
    e = np.minimum(e, 0.5 * p)
    q = EPSILON * e / (p - e)
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
