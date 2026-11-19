"""The Fortran reference generator refuses to run on the numba ports."""

import importlib.util
from pathlib import Path

import pytest

SCRIPT = (Path(__file__).resolve().parent.parent / "scripts"
          / "generate_fortran_reference.py")


def _load_script():
    spec = importlib.util.spec_from_file_location(
        "generate_fortran_reference", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_refuses_when_the_components_run_numba(tmp_path, monkeypatch):
    script = _load_script()
    # If the guard ever fails, only tmp_path is written, never the fixtures.
    out_dir = tmp_path / "fortran_reference"
    monkeypatch.setattr(script, "OUT_DIR", str(out_dir))
    with pytest.raises(SystemExit) as excinfo:
        script.main()
    message = str(excinfo.value.code)
    assert "SimplePhysics" in message
    assert "DcmipInitialConditions" in message
    assert "b7c3302e" in message
    assert not out_dir.exists()
