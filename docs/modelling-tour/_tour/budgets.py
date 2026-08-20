"""Read a column's energy and water budgets out of its state.

Pages 11 and 12 do not decide a column has equilibrated by squinting at a
curve; they read these numbers. The same functions are the acceptance criterion
in ``scripts/generate_tour_equilibria.py`` and the residual test that keeps the
shipped equilibrium states from going stale, so "converged" means one thing
across the pages, the generator and CI.

**Sign convention: positive is into the box.** ``toa_imbalance`` positive means
the planet is gaining energy and will warm; ``surface_imbalance`` positive
means the surface is gaining energy.

**Level indexing.** Index 0 is the bottom (the surface interface, or the
lowest model layer); index -1 is the top. climt states are ordered this way
throughout; a page that gets it backwards reports the OLR as the downwelling
longwave at the ground.

**The shortwave.** No page in this tranche runs a shortwave component -- the
absorbed shortwave is written into ``downwelling_shortwave_flux_in_air`` at the
surface interface and the upwelling shortwave is zero. So the atmosphere
absorbs no sunlight, and the planet's absorbed shortwave is exactly the
surface's. That is why :func:`toa_imbalance` reads it out of the state instead
of taking it as an argument: a page cannot then quote a budget computed against
a ``SOLAR`` it did not actually use.
"""
import numpy as np
from sympl import get_constant

DAY_SECONDS = 86400.0


def _column(state, name):
    """A quantity's single column as a plain 1-D float array."""
    values = np.asarray(state[name].values, dtype=float)
    return values.reshape(values.shape[0], -1)[:, 0]


def _scalar(state, name, default=0.0):
    """A surface scalar as a plain float, or ``default`` if absent."""
    if name not in state:
        return float(default)
    return float(np.asarray(state[name].values, dtype=float).ravel()[0])


def pressure_thickness(state):
    """Layer mass in pressure units, (nz,) in Pa, bottom first."""
    p_int = _column(state, "air_pressure_on_interface_levels")
    return p_int[:-1] - p_int[1:]


def column_enthalpy(state):
    """Column-integrated moist enthalpy, J m^-2.

    ``Cp_dry * T + Lv * q``, mass-weighted by ``dp / g`` -- the same integral
    ``tests/test_conservation.py`` uses, so page 10's claim that dry convective
    adjustment conserves it is guarded by a test that already exists.
    """
    Cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure",
                       "J/kg/degK")
    Lv = get_constant("latent_heat_of_condensation", "J/kg")
    g = get_constant("gravitational_acceleration", "m/s^2")

    T = _column(state, "air_temperature")
    q = (_column(state, "specific_humidity")
         if "specific_humidity" in state else np.zeros_like(T))
    return float(np.sum((float(Cpd) * T + float(Lv) * q)
                        * pressure_thickness(state) / float(g)))


def absorbed_shortwave(state):
    """Net downward shortwave at the surface, W m^-2 -- and at the TOA.

    Equal at both, because no page here puts a shortwave absorber in the air.
    """
    return float(_column(state, "downwelling_shortwave_flux_in_air")[0]
                 - _column(state, "upwelling_shortwave_flux_in_air")[0])


def olr(state):
    """Outgoing longwave radiation: upwelling longwave at the top, W m^-2."""
    return float(_column(state, "upwelling_longwave_flux_in_air")[-1])


def toa_imbalance(state):
    """Net downward energy flux at the top of the atmosphere, W m^-2.

    Positive means the planet is gaining energy. At equilibrium this is zero,
    and how close to zero is the convergence test.
    """
    return float(absorbed_shortwave(state) - olr(state))


def surface_imbalance(state):
    """Net downward energy flux into the surface, W m^-2.

    Shortwave in, longwave down in, longwave up out, and the two turbulent
    fluxes out. The turbulent terms are zero on pages that run no boundary
    layer, which is what makes those pages' surface--air discontinuity the
    thing page 8 then goes and erases.
    """
    return float(
        absorbed_shortwave(state)
        + _column(state, "downwelling_longwave_flux_in_air")[0]
        - _column(state, "upwelling_longwave_flux_in_air")[0]
        - _scalar(state, "surface_upward_sensible_heat_flux")
        - _scalar(state, "surface_upward_latent_heat_flux")
    )


def evaporation_rate(state):
    """Evaporation implied by the surface latent heat flux, mm day^-1.

    ``LHF / Lv`` is a mass flux in kg m^-2 s^-1, which over water is mm s^-1.
    Page 12 closes its moisture budget by comparing this with
    :func:`precipitation_rate`. A page with no boundary layer has no latent
    heat flux in its state at all; that reads as zero rather than raising, so
    the budget table is the same code on every page.
    """
    Lv = get_constant("latent_heat_of_condensation", "J/kg")
    return (_scalar(state, "surface_upward_latent_heat_flux")
            / float(Lv) * DAY_SECONDS)


def precipitation_rate(state, timestep):
    """Total precipitation, mm day^-1: convective plus grid-scale.

    ``EmanuelConvectionPython`` reports ``convective_precipitation_rate``
    already in mm day^-1. ``GridScaleCondensation`` reports
    ``precipitation_amount`` in kg m^-2 accumulated over the step it was
    called with, so it needs the timestep to become a rate. Either may be
    absent; a page with neither gets 0.0.

    Args:
        state: a climt state after at least one step.
        timestep: the ``UnytTimeDelta`` the step was taken with.
    """
    convective = _scalar(state, "convective_precipitation_rate")
    grid_scale = _scalar(state, "precipitation_amount")
    seconds = float(timestep.total_seconds())
    return float(convective + grid_scale * DAY_SECONDS / seconds)


def summary(state):
    """The five numbers pages 11 and 12 print beside every experiment."""
    return dict(
        surface_temperature=_scalar(state, "surface_temperature"),
        olr=olr(state),
        absorbed_shortwave=absorbed_shortwave(state),
        toa_imbalance=toa_imbalance(state),
        surface_imbalance=surface_imbalance(state),
    )
