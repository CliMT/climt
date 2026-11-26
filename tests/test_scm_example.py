"""examples/scm_ensemble.py runs and gives finite output."""
import importlib.util
import os

import numpy as np

_PATH = os.path.join(os.path.dirname(__file__), os.pardir, "examples", "scm_ensemble.py")


def _load():
    spec = importlib.util.spec_from_file_location("scm_ensemble", _PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_example_runs_small():
    scm_ensemble = _load()
    result = scm_ensemble.main(n=4, n_steps=20)
    assert result.arrays["air_temperature"].shape[0] >= 1
    assert np.all(np.isfinite(result.arrays["air_temperature"]))
