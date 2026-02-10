import pytest
import sympl
import unyt

import climt
from climt import (
    EmanuelConvection,
    UnytBackend,
    UnytStateContainer,
    get_default_state,
    get_grid,
)


def test_unyt_backend_initialization():
    # Set the backend
    sympl.set_backend(UnytBackend())

    # Create a component list
    components = [EmanuelConvection()]

    # Get default state
    state = get_default_state(components, grid_state=get_grid(nx=4, ny=4, nz=10))

    # Verify objects in state
    for name, value in state.items():
        if name == "time":
            continue
        assert isinstance(value, UnytStateContainer), (
            f"'{name}' is not UnytStateContainer"
        )
        assert isinstance(value.data, unyt.unyt_array), (
            f"'{name}'.data is not unyt_array"
        )

        # Verify specific units to ensure mapping worked
        if name == "longitude":
            # unyt might normalize units, but it should be compatible with degrees
            assert value.data.units.same_dimensions_as(unyt.degree)
        if name == "air_temperature":
            assert value.data.units.same_dimensions_as(unyt.Kelvin)


def test_unyt_backend_unit_parsing():
    backend = UnytBackend()
    # Test space handling
    c = backend.create_quantity([1], "test", "kg m^-2 s^-1", ("x",))
    assert c.data.units == unyt.kg / (unyt.m**2 * unyt.s)

    # Test caret handling
    c = backend.create_quantity([1], "test", "m^2", ("x",))
    assert c.data.units == unyt.m**2


if __name__ == "__main__":
    test_unyt_backend_initialization()
    test_unyt_backend_unit_parsing()
