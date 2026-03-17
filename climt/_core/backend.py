# -*- coding: utf-8 -*-
import numpy as np

try:
    from numba import njit, prange

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    njit = lambda x, **kwargs: x
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
