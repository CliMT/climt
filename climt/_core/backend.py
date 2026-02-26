# -*- coding: utf-8 -*-
import numpy as np
from datetime import timedelta

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
    if HAS_JAX and (isinstance(arr, (jax.Array, jnp.ndarray)) or hasattr(arr, "__jax_array__")):
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
        in_axes = kwargs.get("in_axes", 0); out_axes = kwargs.get("out_axes", 0)
        return jax.vmap(single_column_func, in_axes=in_axes, out_axes=out_axes)
    def manual_vectorization(*args, **kwargs):
        first_array = args[0]; data = getattr(first_array, "values", first_array); data = getattr(data, "data", data); ncol = data.shape[-1]
        results = []
        for i in range(ncol):
            col_args = []
            for arg in args:
                if hasattr(arg, "shape") and arg.ndim > 0: col_args.append(arg[..., i])
                else: col_args.append(arg)
            results.append(single_column_func(*col_args, **kwargs))
        if not results: return None
        xp = get_array_namespace(*args)
        if isinstance(results[0], tuple): return tuple(xp.stack(res, axis=-1) for res in zip(*results))
        else: return xp.stack(results, axis=-1)
    return manual_vectorization

class JaxTimeDelta(timedelta):
    def __new__(cls, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], timedelta):
            td = args[0]
            return super().__new__(cls, days=td.days, seconds=td.seconds, microseconds=td.microseconds)
        return super().__new__(cls, *args, **kwargs)
    def total_seconds(self):
        if not HAS_JAX: return super().total_seconds()
        return jnp.array(super().total_seconds())

class JaxStateContainer:
    def __init__(self, data, dims, units="dimensionless"):
        self.data = data; self.dims = tuple(dims); self.units = units
    def __repr__(self): return f"JaxStateContainer(data={self.data}, dims={self.dims}, units={self.units})"
    def __array__(self, dtype=None): return np.asarray(self.data, dtype=dtype)
    def __jax_array__(self): return self.data
    __array_priority__ = 100.0
    @property
    def shape(self): return self.data.shape
    @property
    def values(self): return self.data
    @property
    def attrs(self): return {"units": self.units}
    def to_units(self, units): return self
    def _align_and_op(self, other, op_func):
        if hasattr(other, "data") and hasattr(other, "dims"):
            if self.dims != other.dims:
                other_indices = {dim: i for i, dim in enumerate(other.dims)}; perm = [other_indices[dim] for dim in self.dims]; aligned_other_data = other.data.transpose(perm)
                return JaxStateContainer(op_func(self.data, aligned_other_data), self.dims, self.units)
            return JaxStateContainer(op_func(self.data, other.data), self.dims, self.units)
        return JaxStateContainer(op_func(self.data, other), self.dims, self.units)
    def __add__(self, other): return self._align_and_op(other, lambda x, y: x + y)
    def __radd__(self, other): return self.__add__(other)
    def __sub__(self, other): return self._align_and_op(other, lambda x, y: x - y)
    def __rsub__(self, other): return JaxStateContainer(other - self.data, self.dims, self.units)
    def __mul__(self, other):
        if isinstance(other, timedelta): other = JaxTimeDelta(other).total_seconds()
        return self._align_and_op(other, lambda x, y: x * y)
    def __rmul__(self, other): return self.__mul__(other)
    def __truediv__(self, other):
        if isinstance(other, timedelta): other = JaxTimeDelta(other).total_seconds()
        return self._align_and_op(other, lambda x, y: x / y)
    def __rtruediv__(self, other): return JaxStateContainer(other / self.data, self.dims, self.units)

if HAS_JAX:
    from jax.tree_util import register_pytree_node
    register_pytree_node(JaxStateContainer, lambda container: ((container.data,), (container.dims, container.units)), lambda aux, children: JaxStateContainer(children[0], aux[0], aux[1]))

class JaxBackend(StateBackend):
    def __init__(self):
        if not HAS_JAX: raise ImportError("JAX is required for JaxBackend")
    def get_array(self, state_value, name, target_units, target_dims, dim_lengths):
        data = getattr(state_value, "data", state_value); current_dims = getattr(state_value, "dims", target_dims)
        if current_dims == target_dims: return data
        src_indices = {dim: i for i, dim in enumerate(current_dims)}; perm = [src_indices[dim] for dim in target_dims if dim in current_dims]; aligned_data = data.transpose(perm)
        reshape_shape = []; current_axis_idx = 0
        for dim in target_dims:
            if dim in current_dims: reshape_shape.append(aligned_data.shape[current_axis_idx]); current_axis_idx += 1
            else: reshape_shape.append(1)
        return jnp.broadcast_to(aligned_data.reshape(reshape_shape), tuple([dim_lengths[dim] for dim in target_dims]))
    def create_quantity(self, data, name, units, dims, reference_state=None): return JaxStateContainer(jnp.array(data), dims, units)
    def get_dims(self, state_value): return getattr(state_value, "dims", None)
    def get_shape(self, state_value): return getattr(state_value, "shape", jnp.shape(state_value))
