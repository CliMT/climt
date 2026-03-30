import pytest
import sympl
from sympl._core.backend import DataArrayBackend


@pytest.fixture(autouse=True)
def reset_sympl_backend():
    """Save and restore sympl backend and constants after each test.

    Prevents backend state leaks (e.g. UnytBackend set in one test
    breaking sympl.set_constant in a later test that runs with a
    backend whose get_container_type() returns None).
    """
    saved_backend = sympl.get_backend()
    yield
    sympl.set_backend(saved_backend)
    sympl.reset_constants()
    # Re-set climt's custom constants that reset_constants() clears
    sympl.set_constant("top_of_model_pressure", 20.0, "Pa")
