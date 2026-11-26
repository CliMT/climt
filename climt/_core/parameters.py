"""Per-column parameters for numba column kernels.

A kernel takes its parameters as a NamedTuple whose fields are arrays of
shape ``(ncol,)``, so that every column can have its own value.
:func:`per_column` builds that from a NamedTuple whose fields may be
scalars.
"""

import copy
import typing

import numpy as np

PARALLEL_MIN_COLUMNS = 4096
"""Kernels use their parallel variant at or above this many columns.

Set from benchmarks/scm_build.py (benchmarks/results/scm_build.json): the
smallest measured count (of 16, 256, 4096) at which NUMBA_NUM_THREADS threads
beat one thread for the GMD column with both Emanuel schemes. The gain is
about 1-2%, because the serial RRTMG Fortran dominates the step.
"""


def _field_dtype(cls, name):
    try:
        annotations = typing.get_type_hints(cls)
    except Exception:
        annotations = {}
    return np.int64 if annotations.get(name) is int else np.float64


def per_column(params, ncol):
    """Return ``params`` with every field as an array of shape ``(ncol,)``.

    Args:
        params (NamedTuple): Fields are scalars or arrays.
        ncol (int): Number of columns.

    Returns:
        NamedTuple: Same type as ``params``. Fields annotated ``int`` in the
        NamedTuple class become int64 arrays, all others float64. A scalar
        or shape-(1,) field is broadcast (a read-only view); a shape-(ncol,)
        field is passed through.

    Raises:
        ValueError: If a field has any other shape. The message names the
        field.
    """
    cls = type(params)
    fields = []
    for name, value in zip(params._fields, params):
        array = np.asarray(value, dtype=_field_dtype(cls, name))
        if array.shape == (ncol,):
            fields.append(array)
        elif array.shape in ((), (1,)):
            fields.append(np.broadcast_to(array.reshape(()), (ncol,)))
        else:
            raise ValueError(
                "parameter {!r} has shape {}; expected a scalar, shape (1,) "
                "or shape ({},)".format(name, array.shape, ncol))
    return cls(*fields)


def with_parameters(component, values, ncol):
    """A shallow copy of ``component`` with per-column parameter values.

    Args:
        component: A component with a NamedTuple ``_params`` and a
            ``per_column_parameters`` class attribute.
        values (dict): Parameter name to value (scalar or shape (ncol,)).
        ncol (int): Number of columns.

    Returns:
        A shallow copy whose ``_params`` is ``per_column(...)`` of the
        original's, with the fields named in ``values`` replaced. The
        original component is not modified.

    Raises:
        ValueError: If a name is not in
        ``type(component).per_column_parameters``. The message names the
        component class, the parameter and the declared names.
    """
    cls = type(component)
    declared = tuple(getattr(cls, "per_column_parameters", ()))
    for name in values:
        if name not in declared:
            raise ValueError(
                "{} does not declare {!r} as a per-column parameter; "
                "declared: {}".format(cls.__name__, name, ", ".join(declared)))
    new = copy.copy(component)
    new._params = per_column(component._params._replace(**values), ncol)
    if hasattr(new, "_column_params"):
        new._column_params = None
    return new
