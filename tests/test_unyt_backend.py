import numpy as np
import pytest
import unyt

from climt._core.unyt_backend import UnytBackend, UnytStateContainer


def test_unyt_backend_create_quantity():
    backend = UnytBackend()
    data = np.array([1, 2, 3])
    name = "temperature"
    units = "degree_Celsius"
    dims = ("z",)

    container = backend.create_quantity(data, name, units, dims)

    assert isinstance(container, UnytStateContainer)
    assert container.dims == dims
    assert container.data.units == unyt.degree_Celsius
    np.testing.assert_array_equal(container.data.v, data)


def test_unyt_backend_get_array_simple():
    backend = UnytBackend()
    data = np.array([273.15, 283.15])
    unyt_arr = unyt.unyt_array(data, "K")
    container = UnytStateContainer(unyt_arr, ("z",))

    # Test valid extraction
    result = backend.get_array(container, "temp", "degC", ("z",), {})
    expected = np.array([0.0, 10.0])

    np.testing.assert_allclose(result, expected)


def test_unyt_backend_get_array_transpose():
    backend = UnytBackend()
    # shape (2, 3) -> (y, x)
    data = np.array([[1, 2, 3], [4, 5, 6]])
    unyt_arr = unyt.unyt_array(data, "m")
    container = UnytStateContainer(unyt_arr, ("y", "x"))

    # Request (x, y) -> shape (3, 2)
    result = backend.get_array(container, "height", "m", ("x", "y"), {})

    expected = data.T
    np.testing.assert_array_equal(result, expected)


def test_unyt_backend_get_array_unit_conversion_error():
    backend = UnytBackend()
    unyt_arr = unyt.unyt_array([1], "m")
    container = UnytStateContainer(unyt_arr, ("x",))

    with pytest.raises(ValueError, match="Unit conversion failed"):
        backend.get_array(container, "height", "kg", ("x",), {})


def test_unyt_backend_get_array_missing_dim_error():
    backend = UnytBackend()
    unyt_arr = unyt.unyt_array([1], "m")
    container = UnytStateContainer(unyt_arr, ("x",))

    with pytest.raises(ValueError, match="Missing dimension"):
        backend.get_array(container, "height", "m", ("x", "y"), {})


if __name__ == "__main__":
    # Manually run if executed as script
    try:
        test_unyt_backend_create_quantity()
        test_unyt_backend_get_array_simple()
        test_unyt_backend_get_array_transpose()
        test_unyt_backend_get_array_unit_conversion_error()
        test_unyt_backend_get_array_missing_dim_error()
        print("All manual tests passed!")
    except Exception as e:
        print(f"Test failed: {e}")
        raise
