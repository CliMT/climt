"""Per-column parameters for numba column kernels.

A kernel takes its parameters as a NamedTuple whose fields are arrays of
shape ``(ncol,)``, so that every column can have its own value.
:func:`per_column` builds that from a NamedTuple whose fields may be
scalars.
"""

import numpy as np

PARALLEL_MIN_COLUMNS = 32
"""Kernels use their parallel variant at or above this many columns."""


def _field_dtype(cls, name):
    annotations = getattr(cls, "__annotations__", {})
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
