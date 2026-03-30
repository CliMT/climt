from copy import deepcopy

import numpy as np
import unyt
from sympl import StateBackend

# Define aliases for units used in climt that are not in unyt by default

try:
    unyt.define_unit("degrees_east", 1.0 * unyt.degree)

except RuntimeError as e:
    if "already defined" not in str(e):
        print(f"Error defining degrees_east: {e}")

        # Optionally raise, but for now print

    pass


try:
    unyt.define_unit("degrees_north", 1.0 * unyt.degree)

except RuntimeError as e:
    if "already defined" not in str(e):
        print(f"Error defining degrees_north: {e}")

    pass


try:
    unyt.define_unit("degK", 1.0 * unyt.Kelvin)

except RuntimeError as e:
    if "already defined" not in str(e):
        print(f"Error defining degK: {e}")

    pass

try:
    unyt.define_unit("degrees_N", 1.0 * unyt.degree)
except RuntimeError:
    pass

try:
    unyt.define_unit("degrees_E", 1.0 * unyt.degree)
except RuntimeError:
    pass

try:
    unyt.define_unit("degrees_W", 1.0 * unyt.degree)
except RuntimeError:
    pass

try:
    unyt.define_unit("degrees_S", 1.0 * unyt.degree)
except RuntimeError:
    pass

try:
    unyt.define_unit("radians", 1.0 * unyt.radian)
except RuntimeError:
    pass


from datetime import timedelta


class UnytTimeDelta(timedelta):
    """
    A timedelta subclass that returns a unyt_quantity with units of seconds
    from total_seconds().
    """

    def __new__(cls, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], timedelta):
            td = args[0]
            return super().__new__(
                cls, days=td.days, seconds=td.seconds, microseconds=td.microseconds
            )
        return super().__new__(cls, *args, **kwargs)

    def total_seconds(self):
        return super().total_seconds() * unyt.s


class UnytStateContainer:
    """
    A lightweight wrapper around a unyt_array to store dimension names.

    Attributes:
        data (unyt.unyt_array): The underlying data with units.
        dims (tuple of str): The names of the dimensions.
    """

    def __init__(self, data, dims):
        if not isinstance(data, (unyt.unyt_array, np.ndarray)):
            # Ideally strict, but helpful to be flexible if easy.
            # For now, strict to match design.
            raise TypeError(
                f"Data must be a unyt.unyt_array or numpy.ndarray, got {type(data)}"
            )
        self.data = data
        self.dims = tuple(dims)

    def __repr__(self):
        return f"UnytStateContainer(data={self.data}, dims={self.dims})"

    def __array__(self, dtype=None):
        return np.asarray(self.data, dtype=dtype)

    __array_priority__ = 100.0

    def copy(self, deep=False):
        if deep:
            return UnytStateContainer(deepcopy(self.data), self.dims)
        else:
            return UnytStateContainer(self.data.copy(), self.dims)

    @property
    def shape(self):
        return self.data.shape

    @property
    def values(self):
        return np.asarray(self.data)

    @property
    def attrs(self):
        if hasattr(self.data, "units"):
            return {"units": str(self.data.units)}
        return {}

    def to_units(self, units):
        if not isinstance(self.data, unyt.unyt_array):
            return self

        # Sanitization logic similar to backend
        sanitized_units = units.replace("^", "**").replace(" ", "*")
        return UnytStateContainer(self.data.to(sanitized_units), self.dims)

    def __add__(self, other):
        if isinstance(other, UnytStateContainer):
            # Align other to self
            if self.dims != other.dims:
                # Basic alignment: check if sets match
                if set(self.dims) != set(other.dims):
                    raise ValueError(
                        f"Cannot add UnytStateContainers with incompatible dimensions: {self.dims} vs {other.dims}"
                    )

                # Transpose other to match self
                # map dim name to axis index in other
                other_indices = {dim: i for i, dim in enumerate(other.dims)}
                perm = [other_indices[dim] for dim in self.dims]
                aligned_other_data = other.data.transpose(perm)
                return UnytStateContainer(self.data + aligned_other_data, self.dims)

            return UnytStateContainer(self.data + other.data, self.dims)

        if not isinstance(
            other, (unyt.unyt_array, np.ndarray, float, int, np.floating, np.integer)
        ):
            return NotImplemented

        other_data = self._to_numpy_if_needed(other)
        return UnytStateContainer(self.data + other_data, self.dims)

    def __radd__(self, other):
        return self.__add__(other)

    def _to_numpy_if_needed(self, obj):
        if isinstance(obj, memoryview):
            return np.asarray(obj)
        return obj

    def __mul__(self, other):
        from datetime import timedelta

        if isinstance(other, timedelta):
            # If it's UnytTimeDelta, total_seconds() already returns unyt.s
            # If it's standard timedelta, it returns float.
            other = other.total_seconds()
            if not isinstance(other, unyt.unyt_array):
                other = other * unyt.s

        if isinstance(other, UnytStateContainer):
            # For multiplication, alignment is also needed if we want element-wise
            if self.dims != other.dims and set(self.dims) == set(other.dims):
                other_indices = {dim: i for i, dim in enumerate(other.dims)}
                perm = [other_indices[dim] for dim in self.dims]
                aligned_other_data = other.data.transpose(perm)
                return UnytStateContainer(self.data * aligned_other_data, self.dims)

            return UnytStateContainer(self.data * other.data, self.dims)

        if not isinstance(
            other, (unyt.unyt_array, np.ndarray, float, int, np.floating, np.integer)
        ):
            return NotImplemented

        other_data = self._to_numpy_if_needed(other)
        return UnytStateContainer(self.data * other_data, self.dims)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __sub__(self, other):
        if isinstance(other, UnytStateContainer):
            if self.dims != other.dims:
                if set(self.dims) != set(other.dims):
                    raise ValueError(
                        f"Cannot subtract UnytStateContainers with incompatible dimensions: {self.dims} vs {other.dims}"
                    )

                other_indices = {dim: i for i, dim in enumerate(other.dims)}
                perm = [other_indices[dim] for dim in self.dims]
                aligned_other_data = other.data.transpose(perm)
                return UnytStateContainer(self.data - aligned_other_data, self.dims)

            return UnytStateContainer(self.data - other.data, self.dims)

        if not isinstance(
            other, (unyt.unyt_array, np.ndarray, float, int, np.floating, np.integer)
        ):
            return NotImplemented

        other_data = self._to_numpy_if_needed(other)
        return UnytStateContainer(self.data - other_data, self.dims)

    def __rsub__(self, other):
        if not isinstance(
            other, (unyt.unyt_array, np.ndarray, float, int, np.floating, np.integer)
        ):
            return NotImplemented
        other_data = self._to_numpy_if_needed(other)
        return UnytStateContainer(other_data - self.data, self.dims)

    def __truediv__(self, other):
        from datetime import timedelta

        if isinstance(other, timedelta):
            other = other.total_seconds()
            if not isinstance(other, unyt.unyt_array):
                other = other * unyt.s

        if isinstance(other, UnytStateContainer):
            return UnytStateContainer(self.data / other.data, self.dims)

        if not isinstance(
            other, (unyt.unyt_array, np.ndarray, float, int, np.floating, np.integer)
        ):
            return NotImplemented

        other_data = self._to_numpy_if_needed(other)
        return UnytStateContainer(self.data / other_data, self.dims)

    def __rtruediv__(self, other):
        if not isinstance(
            other, (unyt.unyt_array, np.ndarray, float, int, np.floating, np.integer)
        ):
            return NotImplemented
        other_data = self._to_numpy_if_needed(other)
        return UnytStateContainer(other_data / self.data, self.dims)


class UnytBackend(StateBackend):
    """
    A Sympl StateBackend that uses unyt for unit handling.
    """

    def __init__(self):
        self._unit_cache = {}

    def _get_cached_unit(self, unit_str):
        if unit_str in self._unit_cache:
            return self._unit_cache[unit_str]

        if unit_str != unit_str.strip():
            raise ValueError(
                f"Unit string '{unit_str}' contains leading or trailing whitespace. This is not allowed."
            )

        # Sanitization logic
        sanitized_units = unit_str.replace("^", "**").replace(" ", "*")

        try:
            unit_obj = unyt.Unit(sanitized_units)
            self._unit_cache[unit_str] = unit_obj
            return unit_obj
        except Exception as e:
            # Fallback or re-raise with context if needed, but for now just let it fail or return string if unyt.Unit fails (though unlikely if valid)
            # If unyt.Unit fails, we can't really do much.
            raise ValueError(
                f"Failed to parse unit '{unit_str}' (sanitized: '{sanitized_units}'): {e}"
            )

    def get_array(self, state_value, name, target_units, target_dims, dim_lengths):
        """
        Extract a raw numpy array from the state container, converting units and aligning dimensions.
        """
        if not isinstance(state_value, UnytStateContainer):
            raise TypeError(
                f"Expected UnytStateContainer for '{name}', got {type(state_value)}"
            )

        # Handle non-numeric data (stored as plain numpy array)
        if not isinstance(state_value.data, unyt.unyt_array):
            # Assume alignment is needed if dims differ, but unit conversion is skipped
            converted_data = state_value.data
        else:
            # 1. Unit Conversion
            try:
                target_unit_obj = self._get_cached_unit(target_units)
                # .to() returns a new unyt_array in the target units
                converted_data = state_value.data.to(target_unit_obj)
            except Exception as e:
                raise ValueError(
                    f"Unit conversion failed for '{name}' from {state_value.data.units} to {target_units}: {e}"
                )

        current_dims = state_value.dims

        # 2. Check for extra dimensions in source that are not in target
        # This implies we are asking for a slice or projection, which this backend doesn't support implicitly yet
        # or it's an error.
        # Typically get_numpy_arrays expects to return the full state including all target dims.
        extra_dims = set(current_dims) - set(target_dims)
        if extra_dims:
            raise ValueError(
                f"State value for '{name}' has dimensions {current_dims} which are not in target dimensions {target_dims}. Extra: {extra_dims}"
            )

        # 3. Transpose to align existing dimensions with target order
        # We find the permutation of current_dims that matches their relative order in target_dims

        # indices of current dims in the source array
        src_indices = {dim: i for i, dim in enumerate(current_dims)}

        # Construct the permutation list:
        # iterate target_dims, if it exists in current_dims, add its source index
        perm = [src_indices[dim] for dim in target_dims if dim in current_dims]

        # Transpose data
        aligned_data = converted_data.transpose(perm)

        # 4. Reshape to insert singleton dimensions for missing dims
        # We need to construct the shape tuple for reshape
        # If dim is in current_dims, use its size (which is now in correct relative order in aligned_data)
        # If dim is NOT in current_dims, use 1

        reshape_shape = []
        # We keep track of which axis of aligned_data we are consuming
        current_axis_idx = 0

        for dim in target_dims:
            if dim in current_dims:
                reshape_shape.append(aligned_data.shape[current_axis_idx])
                current_axis_idx += 1
            else:
                reshape_shape.append(1)

        reshaped_data = aligned_data.reshape(reshape_shape)

        # 5. Broadcast to full target shape
        # Construct full shape from dim_lengths
        try:
            final_shape = tuple([dim_lengths[dim] for dim in target_dims])
        except KeyError as e:
            raise ValueError(
                f"Missing dimension length for '{e}' in dim_lengths. Available: {dim_lengths.keys()}"
            )

        broadcasted_data = np.broadcast_to(reshaped_data, final_shape)

        # Return raw numpy array
        return np.array(broadcasted_data)

    def create_quantity(self, data, name, units, dims, reference_state=None):
        """
        Wrap raw data into a UnytStateContainer.
        """
        # Ensure data is a numpy array (handle lists, etc.)
        data_arr = np.asarray(data)

        # Check if numeric
        if not np.issubdtype(data_arr.dtype, np.number):
            return UnytStateContainer(data_arr, dims)

        # Create unyt array
        try:
            unit_obj = self._get_cached_unit(units)
            unyt_arr = unyt.unyt_array(data_arr, unit_obj)
        except Exception as e:
            print(f"Error creating unyt_array for '{name}' with units '{units}': {e}")
            raise e

        return UnytStateContainer(unyt_arr, dims)

    def get_dims(self, state_value):
        if not isinstance(state_value, UnytStateContainer):
            raise TypeError(f"Expected UnytStateContainer, got {type(state_value)}")
        return state_value.dims

    def get_shape(self, state_value):
        if not isinstance(state_value, UnytStateContainer):
            raise TypeError(f"Expected UnytStateContainer, got {type(state_value)}")
        return state_value.data.shape

    def get_container_type(self):
        return UnytStateContainer
