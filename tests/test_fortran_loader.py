import importlib.machinery
import sys

import pytest

from climt._core import fortran

ABI_MESSAGE = "numpy.core.multiarray failed to import"


@pytest.fixture(autouse=True)
def isolated(monkeypatch):
    for name in [n for n in sys.modules if n == "climt_fortran" or n.startswith("climt_fortran.")]:
        monkeypatch.delitem(sys.modules, name)
    fortran._reset_cache()
    yield
    for name in [n for n in sys.modules if n == "climt_fortran" or n.startswith("climt_fortran.")]:
        sys.modules.pop(name, None)
    fortran._reset_cache()


def install_fake(monkeypatch, tmp_path, interface=1, present=fortran.MODULES, failing=(), wrappers_only=()):
    """A fake climt_fortran package on sys.path. Each module in ``present`` or
    ``failing`` gets a wrapper ``<name>.py`` and a stand-in ``<name>_cffi.py``
    for its compiled module; ``wrappers_only`` modules get the wrapper only."""
    package = tmp_path / "climt_fortran"
    package.mkdir()
    lines = ['__version__ = "9.9.9"']
    if interface is not None:
        lines.append("INTERFACE_VERSION = %d" % interface)
    (package / "__init__.py").write_text("\n".join(lines) + "\n")
    for name in present + wrappers_only:
        (package / (name + ".py")).write_text("LOADED = True\n")
    for name in failing:
        (package / (name + ".py")).write_text("raise ImportError(%r)\n" % ABI_MESSAGE)
    for name in present + failing:
        (package / (name + "_cffi.py")).write_text("")
    monkeypatch.syspath_prepend(str(tmp_path))


def test_missing(monkeypatch):
    monkeypatch.setitem(sys.modules, "climt_fortran", None)
    status = fortran.fortran_status()
    assert status["state"] == "missing" and status["installed"] is False
    assert status["version"] is None and status["interface_version"] is None
    assert status["modules"] == {name: False for name in fortran.MODULES}
    assert fortran.load_extension("_rrtmg_lw") is None
    message = str(fortran.unavailable_error("RRTMGLongwave", "CorkLongwaveRadiation"))
    assert "RRTMGLongwave" in message and "pip install climt-fortran" in message
    assert message.rstrip().endswith("Use CorkLongwaveRadiation instead.")


def test_available(monkeypatch, tmp_path):
    install_fake(monkeypatch, tmp_path)
    status = fortran.fortran_status()
    assert status == {
        "state": "available", "installed": True, "version": "9.9.9",
        "interface_version": 1, "detail": "",
        "modules": {name: True for name in fortran.MODULES},
    }
    assert fortran.load_extension("_emanuel_convection").LOADED is True


@pytest.mark.parametrize("wrappers_only", [(), fortran.MODULES])
def test_stub(monkeypatch, tmp_path, wrappers_only):
    # A stub of climt-fortran 2 ships the Python wrappers without their
    # compiled _cffi modules: still a stub, not broken.
    install_fake(monkeypatch, tmp_path, present=(), wrappers_only=wrappers_only)
    status = fortran.fortran_status()
    assert status["state"] == "stub"
    assert status["detail"] == "climt_fortran contains no compiled extensions"
    message = str(fortran.unavailable_error("SimplePhysics"))
    assert "pip install --force-reinstall --no-deps --no-binary climt-fortran climt-fortran" in message
    assert "gfortran" in message and "instead" not in message


@pytest.mark.parametrize("interface", [2, 0, None])
def test_interface_mismatch(monkeypatch, tmp_path, interface):
    install_fake(monkeypatch, tmp_path, interface=interface)
    status = fortran.fortran_status()
    assert status["state"] == "interface_mismatch"
    assert not any(status["modules"].values())
    assert fortran.load_extension("_rrtmg_sw") is None
    assert "pip install -U climt climt-fortran" in str(fortran.unavailable_error("RRTMGShortwave"))


def test_broken_binary_is_reported_not_raised(monkeypatch, tmp_path):
    present = tuple(n for n in fortran.MODULES if n != "_rrtmg_lw")
    install_fake(monkeypatch, tmp_path, present=present, failing=("_rrtmg_lw",))
    status = fortran.fortran_status()
    assert status["state"] == "broken"
    assert status["modules"]["_rrtmg_lw"] is False and status["modules"]["_emanuel_convection"] is True
    assert "_rrtmg_lw" in status["detail"] and ABI_MESSAGE in status["detail"]
    assert fortran.load_extension("_rrtmg_lw") is None
    assert fortran.load_extension("_emanuel_convection") is not None
    message = str(fortran.unavailable_error("RRTMGLongwave"))
    assert ABI_MESSAGE in message and "pip install --force-reinstall climt-fortran" in message


def test_partly_built_companion_is_broken(monkeypatch, tmp_path):
    install_fake(monkeypatch, tmp_path, present=("_emanuel_convection",))
    status = fortran.fortran_status()
    assert status["state"] == "broken"
    assert "_rrtmg_sw" in status["detail"]


def test_wrapper_without_compiled_module_is_not_loaded(monkeypatch, tmp_path):
    present = tuple(n for n in fortran.MODULES if n != "_rrtmg_sw")
    install_fake(monkeypatch, tmp_path, present=present, wrappers_only=("_rrtmg_sw",))
    status = fortran.fortran_status()
    assert status["state"] == "broken"
    assert status["modules"]["_rrtmg_sw"] is False
    assert "_rrtmg_sw: compiled module not found" in status["detail"]
    assert fortran.load_extension("_rrtmg_sw") is None


def test_extension_module_counts_as_compiled(monkeypatch, tmp_path):
    # climt-fortran 1.x (and 2.x builds part-way through the move to cffi)
    # has Cython extension modules _<name>.<suffix> and no _cffi modules.
    install_fake(monkeypatch, tmp_path, present=())
    package = tmp_path / "climt_fortran"
    for name in fortran.MODULES:
        (package / (name + importlib.machinery.EXTENSION_SUFFIXES[0])).write_bytes(b"not a binary")
    status = fortran.fortran_status()
    # Found, so imported; the fake binary then fails to load.
    assert status["state"] == "broken"
    assert "compiled module not found" not in status["detail"]
    assert all(name in status["detail"] for name in fortran.MODULES)


def test_error_type_and_unknown_name(monkeypatch):
    monkeypatch.setitem(sys.modules, "climt_fortran", None)
    assert isinstance(fortran.unavailable_error("X"), ImportError)
    with pytest.raises(ValueError):
        fortran.load_extension("_not_a_module")


def test_results_are_cached(monkeypatch, tmp_path):
    install_fake(monkeypatch, tmp_path)
    assert fortran.fortran_status()["state"] == "available"
    monkeypatch.setitem(sys.modules, "climt_fortran", None)
    assert fortran.fortran_status()["state"] == "available"      # cached
    fortran._reset_cache()
    assert fortran.fortran_status()["state"] == "missing"


def test_non_import_error_in_extension_is_broken(monkeypatch, tmp_path):
    install_fake(monkeypatch, tmp_path)
    (tmp_path / "climt_fortran" / "_rrtmg_sw.py").write_text('raise ValueError("numpy.dtype size changed")\n')
    status = fortran.fortran_status()
    assert status["state"] == "broken"
    assert "numpy.dtype size changed" in status["detail"]
    assert fortran.load_extension("_rrtmg_sw") is None
    assert "numpy.dtype size changed" in str(fortran.unavailable_error("RRTMGShortwave"))


def test_non_import_error_in_package_is_broken(monkeypatch, tmp_path):
    install_fake(monkeypatch, tmp_path)
    (tmp_path / "climt_fortran" / "__init__.py").write_text('raise RuntimeError("bad dylib")\n')
    status = fortran.fortran_status()
    assert status["state"] == "broken" and status["installed"] is True
    assert not any(status["modules"].values())
    assert "RuntimeError: bad dylib" in status["detail"]
    assert fortran.load_extension("_emanuel_convection") is None
    assert "pip install --force-reinstall climt-fortran" in str(fortran.unavailable_error("X"))


def test_extensions_come_only_from_the_imported_package(monkeypatch, tmp_path):
    # setuptools' editable-install finder sits on sys.meta_path and maps any
    # climt_fortran.<child> to the checkout, whatever the parent's __path__ is.
    # Modules the imported package lacks must not be filled in from elsewhere,
    # or they would bypass its INTERFACE_VERSION check.
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    for name in fortran.MODULES:
        (elsewhere / (name + ".py")).write_text("LOADED = True\n")
        (elsewhere / (name + "_cffi.py")).write_text("")

    class ElsewhereFinder:
        @classmethod
        def find_spec(cls, fullname, path=None, target=None):
            if fullname.rpartition(".")[0] == "climt_fortran":
                return importlib.machinery.PathFinder.find_spec(fullname, [str(elsewhere)])
            return None

    monkeypatch.setattr(sys, "meta_path", sys.meta_path + [ElsewhereFinder])
    package_root = tmp_path / "package"
    package_root.mkdir()
    install_fake(monkeypatch, package_root, present=("_emanuel_convection",))
    status = fortran.fortran_status()
    assert status["state"] == "broken"
    assert status["modules"] == {name: name == "_emanuel_convection" for name in fortran.MODULES}


def test_ported_modules_are_not_loaded():
    assert "_simple_physics" not in fortran.MODULES and "_dcmip" not in fortran.MODULES
