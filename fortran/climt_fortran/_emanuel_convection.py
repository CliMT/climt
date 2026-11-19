"""Emanuel convection (``convect43c.f90``) through cffi.

Same functions, arguments and in-place behaviour as the 1.0 Cython module:
``init_emanuel_convection`` sets the Fortran module's constants once,
``convect`` runs the scheme column by column and fills the caller's output
arrays.

``convect`` checks its arguments in Python, then loops over columns in a
cached numba function that calls the Fortran routine directly. Converting
the 24 arguments to cffi pointers in Python instead made a one-column call
about 1.5 times slower than the Cython version.
"""

import operator

import numba
import numpy as np

from . import _emanuel_convection_cffi as _compiled
from ._arrays import as_double_scalar_pointer, as_int_pointer
from ._numba_cffi import external_function

ffi, lib = _compiled.ffi, _compiled.lib
_emanuel_convection = external_function(_compiled, "emanuel_convection")
_FLOAT64 = np.dtype(np.float64)
_INT32 = np.dtype(np.int32)
_NO_TRACERS = np.zeros((0, 0, 0))

_initialised = False


def init_emanuel_convection(
        py_least_conv_level,
        py_thresh_water_level,
        py_crit_temp, py_entrain_coeff,
        py_downdraft_frac_area,
        py_precip_frac_outside_cloud,
        py_rain_speed, py_snow_speed,
        py_rain_evap_coeff, py_snow_evap_coeff,
        py_mom_tran_coeff, py_max_neg_temp_pert,
        py_beta, py_alpha, py_damp_amp,
        py_Cpd, py_Cpv, py_Cl,
        py_gas_const_vapour, py_gas_const_air,
        py_lat_heat, py_grav,
        py_density_water,
        py_reference_timescale):
    """Set the scheme's parameters and physical constants in the Fortran module."""
    global _initialised
    _initialised = False
    constants = (
        py_thresh_water_level, py_crit_temp, py_entrain_coeff, py_downdraft_frac_area,
        py_precip_frac_outside_cloud, py_rain_speed, py_snow_speed, py_rain_evap_coeff,
        py_snow_evap_coeff, py_mom_tran_coeff, py_max_neg_temp_pert, py_beta, py_alpha,
        py_damp_amp, py_Cpd, py_Cpv, py_Cl, py_gas_const_vapour, py_gas_const_air,
        py_lat_heat, py_grav, py_density_water, py_reference_timescale,
    )
    # The dry convective PBL adjustment modifies inputs, so it has been moved
    # to a separate component. Here we disable it (pbl = 0).
    lib.init_emanuel_convection_fortran(
        as_int_pointer(ffi, 0),
        as_int_pointer(ffi, py_least_conv_level),
        *[as_double_scalar_pointer(ffi, value) for value in constants])
    _initialised = True


def _check(name, array, dtype, shape=None, output=False):
    """Raise unless ``array`` is a numpy array of ``dtype`` (and of ``shape``,
    and writeable if ``output``). Any memory layout is accepted."""
    if not isinstance(array, np.ndarray) or array.dtype != dtype:
        found = array.dtype if isinstance(array, np.ndarray) else type(array).__name__
        raise TypeError("%s: expected a %s array, got %s" % (name, dtype.name, found))
    if shape is not None and array.shape != shape:
        raise ValueError("%s: expected shape %s, got %s" % (name, shape, array.shape))
    if output and not array.flags.writeable:
        raise ValueError("%s: array is read-only" % name)


def convect(
        nlevs,
        ncols,
        max_conv_lev,
        num_tracers,
        dt,
        temp,
        q,
        qs,
        u,
        v,
        pmid,
        pint,
        conv_state,
        precip,
        downdraft_vel_scale,
        downdraft_temp_scale,
        downdraft_q_scale,
        cloud_base_mass_flux,
        cape,
        dtemp,
        dq,
        du,
        dv,
        dtracers=None,
        tracers=None):
    """Run the scheme on each of ``ncols`` columns of ``nlevs`` levels.

    Inputs ``temp``, ``q``, ``qs``, ``u``, ``v``, ``pmid`` (``(ncols, nlevs)``)
    and ``pint`` (``(ncols, nlevs + 1)``) are float64 and are not modified.
    ``conv_state`` (int32), ``precip``, ``downdraft_*_scale``, ``cape``
    (``(ncols,)``) and ``dtemp``, ``dq``, ``du``, ``dv`` (``(ncols, nlevs)``)
    are filled in place; ``cloud_base_mass_flux`` is read and updated in place.
    Arrays may have any memory layout. ``tracers`` (``(ncols, num_tracers,
    nlevs)`` or broadcastable to it; zeros if None) is read, and ``dtracers``
    is left unchanged, as in 1.0.
    """
    if not _initialised:
        raise ValueError('Emanuel scheme not initialised.')
    nlevs, ncols = operator.index(nlevs), operator.index(ncols)
    max_conv_lev, num_tracers = operator.index(max_conv_lev), operator.index(num_tracers)
    dt = float(dt)
    levels = (ncols, nlevs)
    columns = (ncols,)
    for name, array, shape in (
            ("temp", temp, levels), ("q", q, levels), ("qs", qs, levels), ("u", u, levels),
            ("v", v, levels), ("pmid", pmid, levels), ("pint", pint, (ncols, nlevs + 1))):
        _check(name, array, _FLOAT64, shape)
    _check("conv_state", conv_state, _INT32, columns, output=True)
    for name, array, shape in (
            ("precip", precip, columns), ("downdraft_vel_scale", downdraft_vel_scale, columns),
            ("downdraft_temp_scale", downdraft_temp_scale, columns),
            ("downdraft_q_scale", downdraft_q_scale, columns),
            ("cloud_base_mass_flux", cloud_base_mass_flux, columns), ("cape", cape, columns),
            ("dtemp", dtemp, levels), ("dq", dq, levels), ("du", du, levels), ("dv", dv, levels)):
        _check(name, array, _FLOAT64, shape, output=True)
    if num_tracers == 0:
        tracers = _NO_TRACERS
    elif tracers is None:
        tracers = np.zeros((ncols, num_tracers, nlevs))
    else:
        _check("tracers", tracers, _FLOAT64)
        tracers = np.ascontiguousarray(np.broadcast_to(tracers, (ncols, num_tracers, nlevs)))
    _convect_columns(
        nlevs, ncols, max_conv_lev, num_tracers, dt, temp, q, qs, u, v, pmid, pint,
        conv_state, precip, downdraft_vel_scale, downdraft_temp_scale, downdraft_q_scale,
        cloud_base_mass_flux, cape, dtemp, dq, du, dv, tracers)


@numba.njit(cache=True)
def _convect_columns(
        nlevs, ncols, max_conv_lev, num_tracers, dt, temp, q, qs, u, v, pmid, pint,
        conv_state, precip, downdraft_vel_scale, downdraft_temp_scale, downdraft_q_scale,
        cloud_base_mass_flux, cape, dtemp, dq, du, dv, tracers):
    """The 1.0 Cython column loop: Fortran takes every argument by pointer and
    may write to all of them, so each column goes through private contiguous
    buffers, and results are copied to the caller's arrays after each call."""
    integers = np.empty(3, np.int32)  # nlevs, max_conv_lev, num_tracers
    integers[0] = nlevs
    integers[1] = max_conv_lev
    integers[2] = num_tracers
    timestep = np.empty(1)
    timestep[0] = dt
    profiles = np.empty((6, nlevs))  # temp, q, qs, u, v, pmid
    interfaces = np.empty(nlevs + 1)  # pint
    tendencies = np.zeros((4, nlevs))  # dtemp, dq, du, dv
    state = np.empty(1, np.int32)  # conv_state
    # precip, downdraft_vel_scale, downdraft_temp_scale, downdraft_q_scale,
    # cloud_base_mass_flux, cape
    diagnostics = np.empty(6)
    # Fortran sees a column's tracers as TRA(nlevs, num_tracers): C order
    # (num_tracers, nlevs).
    tracer = np.zeros((num_tracers, nlevs))
    dtracer = np.zeros((num_tracers, nlevs))

    for col in range(ncols):
        profiles[0] = temp[col]
        profiles[1] = q[col]
        profiles[2] = qs[col]
        profiles[3] = u[col]
        profiles[4] = v[col]
        profiles[5] = pmid[col]
        interfaces[:] = pint[col]
        if num_tracers > 0:
            tracer[:, :] = tracers[col]
        # The 1.0 module passed pointers into these arrays: start from the
        # caller's values (cloud_base_mass_flux is read; the others are kept
        # wherever Fortran does not write them).
        state[0] = conv_state[col]
        diagnostics[0] = precip[col]
        diagnostics[1] = downdraft_vel_scale[col]
        diagnostics[2] = downdraft_temp_scale[col]
        diagnostics[3] = downdraft_q_scale[col]
        diagnostics[4] = cloud_base_mass_flux[col]
        diagnostics[5] = cape[col]

        _emanuel_convection(
            ffi.from_buffer(profiles[0]), ffi.from_buffer(profiles[1]), ffi.from_buffer(profiles[2]),
            ffi.from_buffer(profiles[3]), ffi.from_buffer(profiles[4]), ffi.from_buffer(profiles[5]),
            ffi.from_buffer(interfaces),
            ffi.from_buffer(integers[0:1]), ffi.from_buffer(integers[1:2]), ffi.from_buffer(integers[2:3]),
            ffi.from_buffer(timestep), ffi.from_buffer(state),
            ffi.from_buffer(tendencies[0]), ffi.from_buffer(tendencies[1]),
            ffi.from_buffer(tendencies[2]), ffi.from_buffer(tendencies[3]),
            ffi.from_buffer(diagnostics[0:1]), ffi.from_buffer(diagnostics[1:2]),
            ffi.from_buffer(diagnostics[2:3]), ffi.from_buffer(diagnostics[3:4]),
            ffi.from_buffer(diagnostics[4:5]), ffi.from_buffer(diagnostics[5:6]),
            ffi.from_buffer(tracer), ffi.from_buffer(dtracer))

        conv_state[col] = state[0]
        precip[col] = diagnostics[0]
        downdraft_vel_scale[col] = diagnostics[1]
        downdraft_temp_scale[col] = diagnostics[2]
        downdraft_q_scale[col] = diagnostics[3]
        cloud_base_mass_flux[col] = diagnostics[4]
        cape[col] = diagnostics[5]
        dtemp[col] = tendencies[0]
        dq[col] = tendencies[1]
        du[col] = tendencies[2]
        dv[col] = tendencies[3]
