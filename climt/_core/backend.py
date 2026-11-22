# -*- coding: utf-8 -*-
import types

import numpy as np

try:
    from numba import njit, prange

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    def njit(x, **kwargs):
        return x

    prange = range


def jit_compile(func=None, backend=np, **kwargs):
    if func is None or not callable(func):
        actual_backend = func if func is not None else backend

        def decorator(f):
            return jit_compile(f, backend=actual_backend, **kwargs)

        return decorator
    if HAS_NUMBA and backend is np:
        return njit(func, **kwargs)
    return func


def _renamed(func, name):
    """A copy of ``func`` under another name.

    numba's on-disk cache files entries under the function's name and first
    line and keys them by signature, not by compile options, so a serial and
    a parallel dispatcher of one function would load each other's code.
    """
    clone = types.FunctionType(func.__code__, func.__globals__, name,
                               func.__defaults__, func.__closure__)
    clone.__qualname__ = name
    return clone
