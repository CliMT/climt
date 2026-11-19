import importlib
import importlib.metadata
import os

import pytest

import climt_fortran

EXPECTED = ("_dcmip", "_emanuel_convection", "_rrtmg_lw", "_rrtmg_sw", "_simple_physics")
REQUIRED = os.environ.get("CLIMT_REQUIRE_FORTRAN") == "1"
ENTRY_POINTS = {
    "_dcmip": ("get_baroclinic_wave_ics", "get_tropical_cyclone_ics"),
    "_emanuel_convection": ("init_emanuel_convection", "convect"),
    "_rrtmg_lw": ("set_constants", "initialise_rrtm_radiation", "rrtm_calculate_longwave_fluxes"),
    "_rrtmg_sw": ("set_constants", "initialise_rrtm_radiation", "rrtm_calculate_shortwave_fluxes"),
    "_simple_physics": ("init_simple_physics",),
}


def test_metadata():
    assert climt_fortran.INTERFACE_VERSION == 1
    assert climt_fortran.MODULES == EXPECTED
    assert climt_fortran.__version__ == importlib.metadata.version("climt-fortran")


def test_available_reports_every_module_without_importing():
    result = climt_fortran.available()
    assert tuple(sorted(result)) == EXPECTED
    assert all(isinstance(v, bool) for v in result.values())
    if REQUIRED:
        assert all(result.values()), result


@pytest.mark.parametrize("name", EXPECTED)
def test_extension_exposes_the_functions_climt_calls(name):
    if not climt_fortran.available()[name]:
        if REQUIRED:
            pytest.fail(name + " is not built")
        pytest.skip("stub install")
    module = importlib.import_module("climt_fortran." + name)
    for function in ENTRY_POINTS[name]:
        assert callable(getattr(module, function))


def test_stub_wheel_has_no_binaries(tmp_path):
    import glob, subprocess, sys, zipfile
    from pathlib import Path
    project = Path(__file__).resolve().parent.parent
    env = dict(os.environ, CLIMT_FORTRAN_STUB="1")
    subprocess.run([sys.executable, "-m", "pip", "wheel", str(project), "--no-deps", "-w", str(tmp_path)],
                   env=env, check=True)
    wheels = glob.glob(str(tmp_path / "climt_fortran-*-py3-none-any.whl"))
    assert len(wheels) == 1
    names = zipfile.ZipFile(wheels[0]).namelist()
    assert "climt_fortran/__init__.py" in names
    assert not [n for n in names if n.endswith((".so", ".pyd", ".pyx", ".c", ".a", ".f90", ".F90"))]
