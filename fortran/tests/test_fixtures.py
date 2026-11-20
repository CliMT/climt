"""The recorded Cython fixtures hold every call the parity and speed tests use."""

import json
import os

import numpy as np
import pytest

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
CALLS = ("emanuel", "rrtmg_lw", "rrtmg_lw_mcica", "rrtmg_sw", "rrtmg_sw_mcica")


@pytest.mark.parametrize("call", CALLS)
def test_reference_has_inputs_and_outputs(call):
    with np.load(os.path.join(FIXTURES, "cython_reference.npz")) as reference:
        keys = reference.files
    assert any(key.startswith(call + "/in/") for key in keys)
    assert any(key.startswith(call + "/out/") for key in keys)


@pytest.mark.parametrize("call", CALLS)
def test_timings_have_every_call(call):
    with open(os.path.join(FIXTURES, "cython_timings.json")) as timings_file:
        assert call in json.load(timings_file)
