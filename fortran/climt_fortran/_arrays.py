"""Check numpy arrays and turn them, and Python scalars, into the pointers the
Fortran routines take.

Shared by the cffi wrappers. ``check_array`` validates an argument's dtype,
shape and writability. ``as_double_pointer`` checks before handing out
a pointer, so a wrong dtype or layout raises instead of letting Fortran read
the wrong memory. It never copies: callers make arrays C-contiguous first, as
the ``.pyx`` files did.
"""

import numpy as np


def as_double_pointer(ffi, array):
    """cffi double* to a C-contiguous float64 array. Raises TypeError naming the
    dtype if not float64, ValueError if not C-contiguous (callers make arrays
    contiguous first, as the .pyx did)."""
    if not isinstance(array, np.ndarray) or array.dtype != np.float64:
        found = array.dtype if isinstance(array, np.ndarray) else type(array).__name__
        raise TypeError("expected a float64 array, got %s" % found)
    if not array.flags.c_contiguous:
        raise ValueError("expected a C-contiguous array")
    # require_writable: Fortran may write through any pointer it is given.
    return ffi.from_buffer("double[]", array, require_writable=True)


def as_int_pointer(ffi, value):
    """A new int32_t* holding value."""
    return ffi.new("int32_t *", value)


def as_double_scalar_pointer(ffi, value):
    """A new double* holding value."""
    return ffi.new("double *", value)


def check_array(name, array, dtype, shape=None, output=False):
    """Raise unless ``array`` is a numpy array of ``dtype`` (and of ``shape``,
    and writeable if ``output``). Any memory layout is accepted."""
    if not isinstance(array, np.ndarray) or array.dtype != dtype:
        found = array.dtype if isinstance(array, np.ndarray) else type(array).__name__
        raise TypeError("%s: expected a %s array, got %s" % (name, dtype.name, found))
    if shape is not None and array.shape != shape:
        raise ValueError("%s: expected shape %s, got %s" % (name, shape, array.shape))
    if output and not array.flags.writeable:
        raise ValueError("%s: array is read-only" % name)
