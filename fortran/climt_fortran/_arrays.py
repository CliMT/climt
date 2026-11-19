"""Turn numpy arrays and Python scalars into the pointers the Fortran routines take.

Shared by the cffi wrappers. ``as_double_pointer`` checks before handing out
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
