# -*- coding: utf-8 -*-
from datetime import timedelta

import numpy as np

try:
    import jax
    import jax.numpy as jnp

    HAS_JAX = True
except ImportError:
    HAS_JAX = False
    jnp = None

try:
    from numba import extending, njit, prange

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    njit = lambda x, **kwargs: x
    prange = range
    extending = None

from sympl import StateBackend


def get_array_namespace(*arrays):
    """
    Returns the appropriate array namespace (numpy or jax.numpy)
    based on the input arrays.
    """
    if HAS_JAX:
        for arr in arrays:
            # Check for JaxStateContainer explicitly
            if hasattr(arr, "__jax_array__"):
                return jnp
            # Only return jnp if it's a real JAX array or tracer
            # memoryviews should return np
            if isinstance(arr, (jax.Array, jnp.ndarray)):
                return jnp

            # Handle DataArray if its underlying data is JAX
            if hasattr(arr, "values"):
                data = arr.values
                if isinstance(data, (jax.Array, jnp.ndarray)):
                    return jnp

    return np


def jit_compile(func=None, backend=np, **kwargs):
    if func is None or not callable(func):
        actual_backend = func if func is not None else backend

        def decorator(f):
            return jit_compile(f, backend=actual_backend, **kwargs)

        return decorator
    if HAS_JAX and backend is jnp:
        return jax.jit(func, **kwargs)
    if HAS_NUMBA and backend is np:
        return njit(func, **kwargs)
    return func


def set_item(arr, idx, val):
    if HAS_JAX and (
        isinstance(arr, (jax.Array, jnp.ndarray)) or hasattr(arr, "__jax_array__")
    ):
        data = getattr(arr, "data", arr)
        return data.at[idx].set(val)
    arr[idx] = val
    return arr


if HAS_NUMBA:

    @extending.overload(set_item)
    def overload_set_item(arr, idx, val):
        def impl(arr, idx, val):
            arr[idx] = val
            return arr

        return impl


def vectorize_component(single_column_func, backend=np, **kwargs):
    use_jit = kwargs.get("jit", True)
    if HAS_JAX and backend is jnp and use_jit:
        in_axes = kwargs.get("in_axes", 0)
        out_axes = kwargs.get("out_axes", 0)
        return jax.vmap(single_column_func, in_axes=in_axes, out_axes=out_axes)

    def manual_vectorization(*args, **kwargs):
        first_array = args[0]
        data = getattr(first_array, "values", first_array)
        data = getattr(data, "data", data)
        ncol = data.shape[-1]
        results = []
        for i in range(ncol):
            col_args = []
            for arg in args:
                if hasattr(arg, "shape") and arg.ndim > 0:
                    col_args.append(arg[..., i])
                else:
                    col_args.append(arg)
            results.append(single_column_func(*col_args, **kwargs))
        if not results:
            return None
        xp = get_array_namespace(*args)
        if isinstance(results[0], tuple):
            return tuple(xp.stack(res, axis=-1) for res in zip(*results))
        else:
            return xp.stack(results, axis=-1)

    return manual_vectorization
