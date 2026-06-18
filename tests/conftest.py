import pytest
import sympl
from sympl._core.backend import DataArrayBackend


@pytest.fixture(autouse=True)
def reset_sympl_backend():
    """Save and restore sympl backend, constants, and climt globals after each test.

    Prevents leaks of:
    - sympl backend state (e.g. UnytBackend set in one test
      breaking sympl.set_constant in a later test that runs with a
      backend whose get_container_type() returns None).
    - climt's atmospheric-profile snapshot stack and configurable
      radiation band counts (set by Cork* schemes at init,
      which otherwise cause shape mismatches in later RRTMG tests).
    """
    saved_backend = sympl.get_backend()
    yield
    sympl.set_backend(saved_backend)
    sympl.reset_constants()
    # Re-set climt's custom constants that reset_constants() clears
    sympl.set_constant("top_of_model_pressure", 20.0, "Pa")
    
    try:
        from climt._core.condensibles import _H2O_DEFAULTS
        for _name, (_val, _units) in _H2O_DEFAULTS.items():
            sympl.set_constant(_name, _val, _units)
    except ImportError:
        pass

    # Reset climt-side global state that survives sympl.reset_constants().
    from climt._core import initialization as _init
    _init._num_longwave_bands = None
    _init._num_shortwave_bands = None
    from climt._core import atmospheric_properties as _ap
    _ap._snapshot_stack.clear()
    _ap._active_profile["name"] = None
    _ap._active_profile["path"] = None
