"""Check numpy arrays and turn them, and Python scalars, into the pointers the
Fortran routines take.

Shared by the cffi wrappers. ``check_array`` validates an argument's dtype,
shape and writability. ``as_double_pointer`` checks before handing out
a pointer, so a wrong dtype or layout raises instead of letting Fortran read
the wrong memory. It never copies: callers make arrays C-contiguous first, as
the ``.pyx`` files did.

``input_pointer`` and ``output_pointer`` do both for one float64 argument of
a known shape, in any memory layout: an array that is not C-contiguous is
passed to Fortran as a contiguous copy, and an output copy is written back
by ``copy_back_outputs`` after the call. They are the per-argument step of a
plain-Python wrapper and run tens of times per call, so the common case is
checked inline, and ``check_array`` (with its messages) runs only when that
check fails.
"""

import numpy as np

_FLOAT64 = np.dtype(np.float64)


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


def input_pointer(ffi, name, array, shape):
    """double* to input ``array``, checked to be float64 of ``shape``: the
    array itself if C-contiguous, else a C-contiguous copy."""
    if type(array) is not np.ndarray or array.dtype != _FLOAT64 or array.shape != shape:
        check_array(name, array, _FLOAT64, shape)  # raises, or accepts an ndarray subclass
    if not array.flags.c_contiguous:
        array = np.ascontiguousarray(array)
    # Writable (the third argument): Fortran may write through any pointer.
    return ffi.from_buffer("double[]", array, True)


def output_pointer(ffi, name, array, shape, copy_back):
    """double* to output ``array``, checked to be writeable float64 of
    ``shape``: the array itself if C-contiguous, else a C-contiguous copy,
    recorded in the list ``copy_back`` for ``copy_back_outputs``."""
    if type(array) is not np.ndarray or array.dtype != _FLOAT64 or array.shape != shape:
        check_array(name, array, _FLOAT64, shape, output=True)
    if not array.flags.c_contiguous:
        check_array(name, array, _FLOAT64, shape, output=True)  # writeable, before copying
        buffer = np.ascontiguousarray(array)  # keeps the caller's values
        copy_back.append((array, buffer))
        array = buffer
    return ffi.from_buffer("double[]", array, True)


def copy_back_outputs(copy_back):
    """Copy each output computed in a contiguous copy back to the caller's array."""
    for array, buffer in copy_back:
        array[...] = buffer
