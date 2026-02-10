import numpy as np
import jax
import jax.numpy as jnp
from sympl import StateBackend
from datetime import timedelta
import unyt

# Force CPU usage to avoid GPU/TPU discovery overhead and ensure zero-copy with NumPy
jax.config.update("jax_platform_name", "cpu")
jax.config.update("jax_enable_x64", True)


class JaxTimeDelta(timedelta):
    """
    A timedelta subclass that returns a JAX scalar from total_seconds().
    """

    def __new__(cls, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], timedelta):
            td = args[0]
            return super().__new__(
                cls, days=td.days, seconds=td.seconds, microseconds=td.microseconds
            )
        return super().__new__(cls, *args, **kwargs)

    def total_seconds(self):
        return jnp.array(super().total_seconds(), dtype=jnp.float64)


class JaxStateContainer:
    """
    A lightweight wrapper around a jax.numpy.ndarray to store dimension names and units.

    This object is intended to be immutable to match JAX's functional paradigm.
    """

    def __init__(self, data, dims, units=None):
        if not isinstance(data, (jnp.ndarray, np.ndarray, jax.Array)):
            data = jnp.asarray(data)
        self.data = data
        self.dims = tuple(dims)
        self.units = units

    def __repr__(self):
        return (
            f"JaxStateContainer(data={self.data}, dims={self.dims}, units={self.units})"
        )

    def __array__(self, dtype=None):
        # On CPU, this should be zero-copy for JAX arrays
        return np.asarray(self.data, dtype=dtype)

    __array_priority__ = 100.0

    @property
    def shape(self):
        return self.data.shape

    @property
    def values(self):
        # Return numpy view. Note: JAX arrays are immutable, so this view
        # will not allow in-place modification unless we convert to numpy.
        # But for 'climt', we want to encourage functional updates.
        return np.asarray(self.data)

    @property
    def attrs(self):
        return {"units": str(self.units)} if self.units else {}

    def to_units(self, target_units):
        if self.units is None or self.units == target_units:
            return self

        # Simple unit conversion logic using unyt for the factor
        # This keeps the compute JAX-friendly by just using a scalar multiplication
        try:
            # We use unyt to get the conversion factor
            # Sanitization similar to UnytBackend
            u_src = unyt.Unit(self.units.replace("^", "**").replace(" ", "*"))
            u_dst = unyt.Unit(target_units.replace("^", "**").replace(" ", "*"))
            factor = (1.0 * u_src).to(u_dst).value

            return JaxStateContainer(self.data * factor, self.dims, target_units)
        except Exception as e:
            raise ValueError(
                f"Unit conversion failed from {self.units} to {target_units}: {e}"
            )

    def _align_and_op(self, other, op):
        if isinstance(other, JaxStateContainer):
            if self.dims != other.dims:
                if set(self.dims) != set(other.dims):
                    raise ValueError(
                        f"Incompatible dimensions: {self.dims} vs {other.dims}"
                    )

                other_indices = {dim: i for i, dim in enumerate(other.dims)}
                perm = [other_indices[dim] for dim in self.dims]
                aligned_other_data = other.data.transpose(perm)
                return JaxStateContainer(
                    op(self.data, aligned_other_data), self.dims, self.units
                )

            return JaxStateContainer(op(self.data, other.data), self.dims, self.units)

        # Handle scalar/numpy/jax arrays
        if isinstance(other, (timedelta, JaxTimeDelta)):
            other = other.total_seconds()

        return JaxStateContainer(op(self.data, other), self.dims, self.units)

    def __add__(self, other):
        return self._align_and_op(other, lambda x, y: x + y)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self._align_and_op(other, lambda x, y: x - y)

    def __rsub__(self, other):
        return JaxStateContainer(other - self.data, self.dims, self.units)

    def __mul__(self, other):
        return self._align_and_op(other, lambda x, y: x * y)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return self._align_and_op(other, lambda x, y: x / y)

    def __rtruediv__(self, other):
        return JaxStateContainer(other / self.data, self.dims, self.units)


class JaxBackend(StateBackend):
    """
    A Sympl StateBackend that uses JAX for computation.
    """

    def __init__(self):
        self._unit_cache = {}

    def get_array(self, state_value, name, target_units, target_dims, dim_lengths):
        if not isinstance(state_value, JaxStateContainer):
            raise TypeError(
                f"Expected JaxStateContainer for '{name}', got {type(state_value)}"
            )

        # 1. Unit Conversion
        converted_container = state_value.to_units(target_units)
        data = converted_container.data
        current_dims = converted_container.dims

        # 2. Dimension alignment (Transpose, Reshape, Broadcast)
        # Check for extra dimensions
        extra_dims = set(current_dims) - set(target_dims)
        if extra_dims:
            raise ValueError(f"Extra dimensions in '{name}': {extra_dims}")

        # Construct permutation for transpose
        src_indices = {dim: i for i, dim in enumerate(current_dims)}
        perm = [src_indices[dim] for dim in target_dims if dim in current_dims]
        aligned_data = data.transpose(perm)

        # Construct shape for reshape (insert singletons)
        reshape_shape = []
        current_axis_idx = 0
        for dim in target_dims:
            if dim in current_dims:
                reshape_shape.append(aligned_data.shape[current_axis_idx])
                current_axis_idx += 1
            else:
                reshape_shape.append(1)

        reshaped_data = aligned_data.reshape(reshape_shape)

        # Broadcast to target shape
        final_shape = tuple([dim_lengths[dim] for dim in target_dims])
        broadcasted_data = jnp.broadcast_to(reshaped_data, final_shape)

        # Return a WRITABLE numpy array.
        # JAX arrays are read-only when viewed as numpy.
        # For components that need to modify the array in-place (like many climt wrappers),
        # we must provide a copy or a writable buffer.
        return np.array(broadcasted_data).copy()

    def create_quantity(self, data, name, units, dims, reference_state=None):
        return JaxStateContainer(data, dims, units)

    def get_dims(self, state_value):
        return state_value.dims

    def get_shape(self, state_value):
        return state_value.shape
