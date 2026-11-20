"""Check numpy arrays and turn Python scalars into the pointers the Fortran
routines take.

Shared by the cffi wrappers. ``check_array`` validates an argument's dtype,
shape and writability, so a wrong dtype or shape raises, naming the argument,
instead of letting Fortran read the wrong memory. ``contiguous_arguments``
checks a call's arrays and makes C-contiguous writeable copies of those that
are not, for the numba kernels that pass them to Fortran;
``copy_back_outputs`` writes the output copies back after the call.
"""

import numpy as np

_FLOAT64 = np.dtype(np.float64)


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


def contiguous_arguments(arrays, specs, copy_back):
    """``arrays`` checked against ``specs`` (a ``(name, shape, output)`` each,
    see ``check_array``), as C-contiguous, aligned, writeable arrays: each
    array that is not one is replaced by such a copy, and output copies are
    recorded in the list ``copy_back`` for ``copy_back_outputs``."""
    result = []
    for array, (name, shape, output) in zip(arrays, specs, strict=True):
        check_array(name, array, _FLOAT64, shape, output)
        flags = array.flags
        if not (flags.c_contiguous and flags.aligned and flags.writeable):
            copy = np.array(array, order="C")  # keeps the caller's values
            if output:
                copy_back.append((array, copy))
            array = copy
        result.append(array)
    return result


def copy_back_outputs(copy_back):
    """Copy each output computed in a contiguous copy back to the caller's array."""
    for array, buffer in copy_back:
        array[...] = buffer
