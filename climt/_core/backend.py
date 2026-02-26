# -*- coding: utf-8 -*-
import numpy as np

try:
    import jax
    import jax.numpy as jnp
    HAS_JAX = True
except ImportError:
    HAS_JAX = False
    jnp = None

try:
    from numba import njit, prange, extending
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    njit = lambda x, **kwargs: x
    prange = range
    extending = None

def get_array_namespace(*arrays):
    """
    Returns the appropriate array namespace (numpy or jax.numpy)
    based on the input arrays.
    """
    if HAS_JAX:
        for arr in arrays:
            data = getattr(arr, "values", arr)
            data = getattr(data, "data", data)
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
    """
    A backend-agnostic way to set items in an array.
    Uses functional updates for JAX, in-place for NumPy.
    """
    if HAS_JAX and isinstance(arr, (jax.Array, jnp.ndarray)):
        return arr.at[idx].set(val)
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
    if HAS_JAX and backend is jnp:
        in_axes = kwargs.get("in_axes", 0)
        out_axes = kwargs.get("out_axes", 0)
        return jax.vmap(single_column_func, in_axes=in_axes, out_axes=out_axes)
        
    def manual_vectorization(*args, **kwargs):
        first_array = args[0]
        ncol = first_array.shape[-1]
        results = []
        for i in range(ncol):
            col_args = [
                arg[..., i] if isinstance(arg, np.ndarray) and arg.ndim > 0 else arg 
                for arg in args
            ]
            results.append(single_column_func(*col_args, **kwargs))
        if not results: return None
        if isinstance(results[0], tuple):
            return tuple(np.stack(res, axis=-1) for res in zip(*results))
        else:
            return np.stack(results, axis=-1)
    return manual_vectorization
