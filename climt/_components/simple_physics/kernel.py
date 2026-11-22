"""numba port of the Reed and Jablonowski (2012) simple physics.

A line-by-line port of ``simple_physics_func`` in
``fortran/_lib/simple_physics/simple_physics_custom.f90``, with the level
reversal and layer thickness of the Cython wrapper ``_simple_physics.pyx`` of
the climt-fortran package.
The Fortran module's constants become fields of :class:`SimplePhysicsParams`
(one value per column) and its switches become :class:`SimplePhysicsFlags`.
"""

from typing import NamedTuple

import numpy as np

from ..._core.backend import HAS_NUMBA, _renamed, prange
from ..._core.parameters import PARALLEL_MIN_COLUMNS

if HAS_NUMBA:
    from numba import njit


class SimplePhysicsParams(NamedTuple):
    gravitational_acceleration: float              # m s-2
    heat_capacity_of_dry_air_at_constant_pressure: float   # J kg-1 K-1
    gas_constant_of_dry_air: float                 # J kg-1 K-1
    gas_constant_of_vapor_phase: float             # J kg-1 K-1
    latent_heat_of_condensation: float             # J kg-1
    planetary_radius: float                        # m
    planetary_rotation_rate: float                 # s-1
    density_of_liquid_water: float                 # kg m-3
    top_of_boundary_layer: float                   # Pa
    boundary_layer_influence_height: float         # Pa
    drag_coefficient_heat_fluxes: float
    base_momentum_drag_coefficient: float
    wind_dependent_momentum_drag_coefficient: float
    maximum_momentum_drag_coefficient: float


class SimplePhysicsFlags(NamedTuple):
    simulate_cyclone: bool
    large_scale_condensation: bool
    boundary_layer: bool
    surface_fluxes: bool
    use_external_surface_temperature: bool
    use_external_surface_specific_humidity: bool


# Literal constants of the Fortran. Literals written there without the _r8
# kind suffix are default (single precision) reals, so they carry their
# float32 value; that includes pi = 4.*atan(1.).
_f32 = np.float32
_PI = float(_f32(4.0) * np.arctan(_f32(1.0)))
_SST_TC = 302.15       # SST for the tropical cyclone test (K)
# Latitude arrives in degrees; the conversion is not in the Fortran, so it
# uses the exact float64 pi (as numpy.deg2rad does).
_DEG_TO_RAD = np.pi / 180.0
_T0 = 273.16           # control temperature for qsat (K)
_E0 = 610.78           # saturation vapour pressure at T0 (Pa)
_V20 = 20.0            # wind speed above which Cd = Cm (m s-1)
_P0 = 100000.0         # reference pressure for potential temperature (Pa)
_T00 = 288.0           # moist baroclinic wave: mean surface T (K)
_U0 = 35.0             # moist baroclinic wave: zonal wind constant (m s-1)
_LATW = 2.0 * _PI / 9.0            # moist baroclinic wave: halfwidth
_ETA0 = 0.252                      # moist baroclinic wave: jet centre
_ETAV = (1.0 - _ETA0) * 0.5 * _PI  # moist baroclinic wave: auxiliary
_Q0 = float(_f32(0.021))           # moist baroclinic wave: maximum q
# Surface saturation vapour pressure (Buck 1981) over water, Tsurf > 271 K,
_W_A, _W_B, _W_C = float(_f32(1.0007)), float(_f32(3.46e-8)), float(_f32(611.21))
_W_D, _W_E = float(_f32(17.966)), float(_f32(247.15))
# ... and over ice.
_I_A, _I_B, _I_C = float(_f32(1.0003)), float(_f32(4.18e-8)), float(_f32(611.15))
_I_D, _I_E = float(_f32(22.452)), 272.5
_EPS_FACTOR = float(_f32(0.378))

# Rows of the per-column work array. The column profiles are stored top
# level first (index 0 at the model top), the Fortran's order.
(_T, _Q, _U, _V, _PMID, _PINT, _PDEL, _RPDEL, _KM, _KE,
 _CA, _CC, _CE, _CAM, _CCM, _CEM, _CFU, _CFV, _CFT, _CFQ) = range(20)
_NWORK = 20


def _simple_physics_columns(u, v, t, q, p, p_int, ps, ts, qs, lat, dt,
                            params, flags):
    nlev, ncol = t.shape
    kb = nlev - 1   # lowest full level, Fortran pver
    ks = nlev       # surface interface, Fortran pver+1
    t_out = np.empty((nlev, ncol))
    u_out = np.empty((nlev, ncol))
    v_out = np.empty((nlev, ncol))
    q_out = np.empty((nlev, ncol))
    precip = np.zeros(ncol)
    sensible_heat_flux = np.zeros(ncol)
    latent_heat_flux = np.zeros(ncol)

    for i in prange(ncol):
        gravit = params.gravitational_acceleration[i]
        cpair = params.heat_capacity_of_dry_air_at_constant_pressure[i]
        rair = params.gas_constant_of_dry_air[i]
        rh2o = params.gas_constant_of_vapor_phase[i]
        latvap = params.latent_heat_of_condensation[i]
        a = params.planetary_radius[i]
        omega = params.planetary_rotation_rate[i]
        rhow = params.density_of_liquid_water[i]
        pbltop = params.top_of_boundary_layer[i]
        pblconst = params.boundary_layer_influence_height[i]
        c_flux = params.drag_coefficient_heat_fluxes[i]
        cd0 = params.base_momentum_drag_coefficient[i]
        cd1 = params.wind_dependent_momentum_drag_coefficient[i]
        cm = params.maximum_momentum_drag_coefficient[i]
        epsilo = rair / rh2o
        zvir = (rh2o / rair) - 1.0

        # Zeroed: Km and Ke are read by the boundary layer even when the
        # surface fluxes, which set them, are off.
        w = np.zeros((_NWORK, nlev + 1))
        tc = w[_T]
        qc = w[_Q]
        uc = w[_U]
        vc = w[_V]
        pmid = w[_PMID]
        pint = w[_PINT]
        pdel = w[_PDEL]
        rpdel = w[_RPDEL]
        km = w[_KM]
        ke = w[_KE]
        ca = w[_CA]
        cc = w[_CC]
        ce = w[_CE]
        cam = w[_CAM]
        ccm = w[_CCM]
        cem = w[_CEM]
        cfu = w[_CFU]
        cfv = w[_CFV]
        cft = w[_CFT]
        cfq = w[_CFQ]

        # Reverse the level axis (climt has index 0 at the surface) and take
        # the layer thickness and its inverse, as _simple_physics.pyx does.
        for k in range(nlev):
            r = kb - k
            tc[k] = t[r, i]
            qc[k] = q[r, i]
            uc[k] = u[r, i]
            vc[k] = v[r, i]
            pmid[k] = p[r, i]
        for k in range(nlev + 1):
            pint[k] = p_int[ks - k, i]
        for k in range(nlev):
            pdel[k] = pint[k + 1] - pint[k]
            rpdel[k] = 1.0 / pdel[k]

        # Height of the lowest full level.
        dlnpint = np.log(ps[i]) - np.log(pint[kb])
        za = rair / gravit * tc[kb] * (1.0 + zvir * qc[kb]) * 0.5 * dlnpint

        # Surface temperature: external; the constant SST of the tropical
        # cyclone test (5-1); or the latitude-dependent SST of the moist
        # baroclinic wave (test 4-3), whose formula takes radians.
        if flags.use_external_surface_temperature:
            tsurf = ts[i]
        elif flags.simulate_cyclone:
            tsurf = _SST_TC
        else:
            latr = lat[i] * _DEG_TO_RAD
            sinlat = np.sin(latr)
            coslat = np.cos(latr)
            tsurf = (
                (_T00 + _PI * _U0 / rair * 1.5 * np.sin(_ETAV) * (np.cos(_ETAV)) ** 0.5 *
                 ((-2.0 * sinlat ** 6 * (coslat ** 2 + 1.0 / 3.0) + 10.0 / 63.0) *
                  _U0 * (np.cos(_ETAV)) ** 1.5 +
                  (8.0 / 5.0 * coslat ** 3 * (sinlat ** 2 + 2.0 / 3.0) - _PI / 4.0) *
                  a * omega * 0.5)) /
                (1.0 + zvir * _Q0 * np.exp(-(latr / _LATW) ** 4)))

        # Large-scale condensation.
        precl = 0.0
        if flags.large_scale_condensation:
            for k in range(nlev):
                qsat = epsilo * _E0 / pmid[k] * np.exp(
                    -latvap / rh2o * ((1.0 / tc[k]) - 1.0 / _T0))
                if qc[k] > qsat:
                    tmp = 1.0 / dt * (qc[k] - qsat) / (
                        1.0 + (latvap / cpair) *
                        (epsilo * latvap * qsat / (rair * tc[k] ** 2)))
                    dtdt = latvap / cpair * tmp
                    dqdt = -tmp
                    precl = precl + tmp * pdel[k] / (gravit * rhow)
                    tc[k] = tc[k] + dtdt * dt
                    qc[k] = qc[k] + dqdt * dt

        # Surface fluxes of momentum, sensible heat and moisture.
        if flags.surface_fluxes:
            wind = np.sqrt(uc[kb] ** 2 + vc[kb] ** 2)
            ke[ks] = c_flux * wind * za
            if wind < _V20:
                cd = cd0 + cd1 * wind
                km[ks] = cd * wind * za
            else:
                cd = cm
                km[ks] = cm * wind * za
            for k in range(nlev):
                if pint[k] >= pbltop:
                    km[k] = km[ks]
                    ke[k] = ke[ks]
                else:
                    decay = np.exp(-(pbltop - pint[k]) ** 2 / (pblconst) ** 2)
                    km[k] = km[ks] * decay
                    ke[k] = ke[ks] * decay

            uc[kb] = uc[kb] / (1.0 + cd * wind * dt / za)
            vc[kb] = vc[kb] / (1.0 + cd * wind * dt / za)

            rho = pmid[kb] / (rair * tc[kb])
            temp_flux = c_flux * wind * (tsurf - tc[kb])
            sensible_heat_flux[i] = rho * cpair * temp_flux
            dtdt = temp_flux * (rho * gravit) / (pint[ks] - pint[kb])
            tc[kb] = tc[kb] + dtdt * dt

            if flags.use_external_surface_specific_humidity:
                qsats = qs[i]
            else:
                if tsurf > 271.0:
                    esats = (_W_A + _W_B * ps[i]) * _W_C * np.exp(
                        _W_D * (tsurf - 273.0) / (_W_E + (tsurf - 273.0)))
                else:
                    esats = (_I_A + _I_B * ps[i]) * _I_C * np.exp(
                        _I_D * (tsurf - 273.0) / (_I_E + (tsurf - 273.0)))
                qsats = epsilo * esats / (ps[i] - _EPS_FACTOR * esats)
            rho = pmid[kb] / (rair * tc[kb])
            temp_flux = c_flux * wind * (qsats - qc[kb])
            latent_heat_flux[i] = latvap * rho * temp_flux
            dqdt = temp_flux * (rho * gravit) / (pint[ks] - pint[kb])
            qc[kb] = qc[kb] + dqdt * dt

        # Boundary layer: implicit vertical diffusion of u, v, potential
        # temperature and q (tridiagonal solve, top to bottom).
        if flags.boundary_layer:
            for k in range(nlev - 1):
                rho = (pint[k + 1] / (rair * (tc[k + 1] + tc[k]) / 2.0))
                dpmid = pmid[k + 1] - pmid[k]
                cam[k] = rpdel[k] * dt * gravit * gravit * km[k + 1] * rho * rho / dpmid
                ccm[k + 1] = rpdel[k + 1] * dt * gravit * gravit * km[k + 1] * rho * rho / dpmid
                ca[k] = rpdel[k] * dt * gravit * gravit * ke[k + 1] * rho * rho / dpmid
                cc[k + 1] = rpdel[k + 1] * dt * gravit * gravit * ke[k + 1] * rho * rho / dpmid
            cam[kb] = 0.0
            ccm[0] = 0.0
            cem[ks] = 0.0
            ca[kb] = 0.0
            cc[0] = 0.0
            ce[ks] = 0.0
            cfu[ks] = 0.0
            cfv[ks] = 0.0
            cft[ks] = 0.0
            cfq[ks] = 0.0
            kappa = rair / cpair
            for k in range(kb, -1, -1):
                denom_m = 1.0 + cam[k] + ccm[k] - cam[k] * cem[k + 1]
                denom_e = 1.0 + ca[k] + cc[k] - ca[k] * ce[k + 1]
                ce[k] = cc[k] / denom_e
                cem[k] = ccm[k] / denom_m
                cfu[k] = (uc[k] + cam[k] * cfu[k + 1]) / denom_m
                cfv[k] = (vc[k] + cam[k] * cfv[k + 1]) / denom_m
                cft[k] = ((_P0 / pmid[k]) ** kappa * tc[k] + ca[k] * cft[k + 1]) / denom_e
                cfq[k] = (qc[k] + ca[k] * cfq[k + 1]) / denom_e

            uc[0] = cfu[0]
            vc[0] = cfv[0]
            tc[0] = cft[0] * (pmid[0] / _P0) ** kappa
            qc[0] = cfq[0]
            for k in range(1, nlev):
                uc[k] = cem[k] * uc[k - 1] + cfu[k]
                vc[k] = cem[k] * vc[k - 1] + cfv[k]
                tc[k] = (ce[k] * tc[k - 1] * (_P0 / pmid[k - 1]) ** kappa + cft[k]) * \
                    (pmid[k] / _P0) ** kappa
                qc[k] = ce[k] * qc[k - 1] + cfq[k]

        # Back to climt's level order.
        precip[i] = precl
        for k in range(nlev):
            r = kb - k
            t_out[r, i] = tc[k]
            q_out[r, i] = qc[k]
            u_out[r, i] = uc[k]
            v_out[r, i] = vc[k]

    return (t_out, u_out, v_out, q_out, precip, sensible_heat_flux,
            latent_heat_flux)


if HAS_NUMBA:
    _simple_physics_serial = njit(cache=True)(
        _renamed(_simple_physics_columns, "_simple_physics_serial"))
    _simple_physics_parallel = njit(parallel=True, cache=True)(
        _renamed(_simple_physics_columns, "_simple_physics_parallel"))
else:
    _simple_physics_serial = _simple_physics_columns
    _simple_physics_parallel = _simple_physics_columns


def simple_physics_kernel(u, v, t, q, p, p_int, ps, ts, qs, lat, dt, params,
                          flags):
    """Advance one timestep of the simple physics.

    Args:
        u, v, t, q, p: Arrays of shape (nlev, ncol): eastward and northward
            wind (m s-1), air temperature (K), specific humidity (kg kg-1)
            and air pressure (Pa).
        p_int: Air pressure on interface levels (Pa), shape (nlev+1, ncol).
        ps, ts, qs, lat: Arrays of shape (ncol,): surface air pressure (Pa),
            surface temperature (K), surface specific humidity (kg kg-1) and
            latitude (degrees_north).
        dt (float): Timestep (s).
        params (SimplePhysicsParams): Fields are arrays of shape (ncol,);
            see :func:`climt._core.parameters.per_column`.
        flags (SimplePhysicsFlags): Switches.

    Level index 0 is nearest the surface (climt's order).

    Returns:
        tuple: ``(t_out, u_out, v_out, q_out, precip, sensible_heat_flux,
        latent_heat_flux)``, new arrays; the inputs are not modified.
        ``latent_heat_flux`` is not clipped at zero (the component clips it).
    """
    if t.shape[1] >= PARALLEL_MIN_COLUMNS:
        return _simple_physics_parallel(u, v, t, q, p, p_int, ps, ts, qs, lat,
                                        float(dt), params, flags)
    return _simple_physics_serial(u, v, t, q, p, p_int, ps, ts, qs, lat,
                                  float(dt), params, flags)
