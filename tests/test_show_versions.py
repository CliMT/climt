import io
import subprocess
import sys

import climt
from climt._core.versions import REPORTED_PACKAGES


def test_show_versions_returns_core_keys_and_every_reported_package():
    out = io.StringIO()
    info = climt.show_versions(file=out)
    assert info["climt"] == climt.__version__
    assert info["python"].count(".") == 2
    assert isinstance(info["fortran_extensions"], bool)
    assert isinstance(info["platform"], str) and info["platform"]
    for name in REPORTED_PACKAGES:
        assert name in info
    assert info["numpy"] is not None
    assert info["sympl"] is not None


def test_show_versions_prints_one_line_per_key_in_order():
    out = io.StringIO()
    info = climt.show_versions(file=out)
    lines = out.getvalue().splitlines()
    assert [line.split(":")[0] for line in lines] == list(info.keys())
    assert lines[0] == "climt: {}".format(climt.__version__)


def test_missing_package_is_none_and_printed_as_not_installed(monkeypatch):
    import climt._core.versions as versions
    monkeypatch.setattr(
        versions, "REPORTED_PACKAGES", ("numpy", "surely_not_a_real_package_xyz"))
    out = io.StringIO()
    info = versions.show_versions(file=out)
    assert info["surely_not_a_real_package_xyz"] is None
    assert "surely_not_a_real_package_xyz: not installed" in out.getvalue()


def test_show_versions_does_not_import_reported_packages():
    code = (
        "import sys, io, climt\n"
        "before = set(sys.modules)\n"
        "climt.show_versions(file=io.StringIO())\n"
        "new = set(sys.modules) - before\n"
        "bad = [m for m in ('jax', 'jaxlib', 'zarr', 'zmq', 'panel', 'holoviews',"
        " 'netCDF4', 'gfs_dynamical_core') if m in new]\n"
        "print('IMPORTED:' + ','.join(bad))\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert result.stdout.strip().splitlines()[-1] == "IMPORTED:"


def test_show_versions_survives_a_broken_compiled_extension(monkeypatch):
    # climt issue #174: a numpy ABI mismatch surfaces as ValueError on first
    # import of a compiled extension, not ImportError.
    def broken():
        raise ValueError("numpy.dtype size changed, may indicate binary incompatibility")

    monkeypatch.setattr(climt, "has_fortran_extensions", broken)
    out = io.StringIO()
    info = climt.show_versions(file=out)
    assert info["fortran_extensions"] is False
    line = next(text for text in out.getvalue().splitlines()
                if text.startswith("fortran_extensions:"))
    assert "ValueError" in line and "dtype size changed" in line


def test_show_versions_reports_the_fortran_companion():
    assert "climt-fortran" in REPORTED_PACKAGES
    info = climt.show_versions(file=io.StringIO())
    assert "climt-fortran" in info
