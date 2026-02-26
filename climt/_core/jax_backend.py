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

from sympl import StateBackend

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
        self.data = jnp.array(data) if HAS_JAX else np.array(data)
        self.dims = tuple(dims)
        self.units = units

    def __repr__(self):
        return f"JaxStateContainer(data={self.data}, dims={self.dims}, units={self.units})"

    def __array__(self, dtype=None):
        return np.asarray(self.data, dtype=dtype)

    def __jax_array__(self):
        return self.data

    __array_priority__ = 100.0

    @property
    def shape(self): return self.data.shape
    @property
    def values(self): return self.data
    @property
    def attrs(self): return {"units": self.units}

    def _align_and_op(self, other, op_func):
        if isinstance(other, JaxStateContainer):
            if self.dims != other.dims:
                other_indices = {dim: i for i, dim in enumerate(other.dims)}
                perm = [other_indices[dim] for dim in self.dims]
                aligned_other_data = other.data.transpose(perm)
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
    register_pytree_node(
        JaxStateContainer,
        lambda container: ((container.data,), (container.dims, container.units)),
        lambda aux, children: JaxStateContainer(children[0], aux[0], aux[1])
    )

class JaxBackend(StateBackend):
    def __init__(self):
        if not HAS_JAX: raise ImportError("JAX is required for JaxBackend")

    def get_array(self, state_value, name, target_units, target_dims, dim_lengths):
        if not isinstance(state_value, JaxStateContainer):
            if isinstance(state_value, (jax.Array, jnp.ndarray)): return state_value
            raise TypeError(f"Expected JaxStateContainer for '{name}', got {type(state_value)}")
        data = state_value.data
        if state_value.dims == target_dims: return data
        src_indices = {dim: i for i, dim in enumerate(state_value.dims)}
        perm = [src_indices[dim] for dim in target_dims if dim in state_value.dims]
        aligned_data = data.transpose(perm)
        reshape_shape = []
        current_axis_idx = 0
        for dim in target_dims:
            if dim in state_value.dims:
                reshape_shape.append(aligned_data.shape[current_axis_idx]); current_axis_idx += 1
            else:
                reshape_shape.append(1)
        return jnp.broadcast_to(aligned_data.reshape(reshape_shape), tuple([dim_lengths[dim] for dim in target_dims]))

    def create_quantity(self, data, name, units, dims, reference_state=None):
        return JaxStateContainer(data, dims, units)

    def get_dims(self, state_value):
        if hasattr(state_value, "dims"): return state_value.dims
        return None

    def get_shape(self, state_value):
        if hasattr(state_value, "shape"): return state_value.shape
        return getattr(state_value, "data", state_value).shape
