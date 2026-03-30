import pytest
import sympl
from sympl._core.backend import DataArrayBackend


@pytest.fixture(autouse=True)
def reset_sympl_backend():
    """Reset sympl backend to DataArrayBackend and reset constants after each test.

    Prevents backend state leaks (e.g. UnytBackend set in one test
    breaking sympl.set_constant in a later test that runs with a
    backend whose get_container_type() returns None).
    Also resets constants to ensure isolation.
    """
    yield
    sympl.set_backend(DataArrayBackend())
    sympl.reset_constants()
