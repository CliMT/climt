import subprocess
import sys

import climt


def test_has_fortran_extensions_is_bool():
    assert isinstance(climt.has_fortran_extensions(), bool)


def test_no_import_warning():
    """Importing climt in a fresh interpreter must be silent.

    Regression guard for import-time noise (print/logging.warning/warnings.warn)
    about missing compiled Fortran extensions. A fresh subprocess is required
    because ``importlib.reload`` does not re-execute the module-level
    try/except blocks of already-imported components, and ``recwarn`` cannot
    see ``print``/``logging.warning`` output anyway.
    """
    result = subprocess.run(
        [sys.executable, "-c", "import climt"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    combined_output = (result.stdout + result.stderr).lower()
    noise_keywords = ("compiled", "fortran", "extension")
    offending_lines = [
        line
        for line in combined_output.splitlines()
        if any(keyword in line for keyword in noise_keywords)
    ]
    assert not offending_lines, (
        "import climt printed extension-related noise:\n"
        f"{offending_lines}\nstdout={result.stdout!r}\nstderr={result.stderr!r}"
    )


def test_climt_works_without_the_companion():
    code = (
        "import sys; sys.modules['climt_fortran'] = None\n"
        "import climt\n"
        "assert climt.has_fortran_extensions() is False\n"
        "try:\n"
        "    climt.RRTMGLongwave()\n"
        "except climt.FortranUnavailableError as error:\n"
        "    assert 'pip install climt-fortran' in str(error)\n"
        "else:\n"
        "    raise SystemExit('RRTMGLongwave did not raise')\n"
        "climt.GrayLongwaveRadiation()\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout == "" and result.stderr == ""


def test_components_never_import_extensions_from_the_old_location():
    import pathlib
    root = pathlib.Path(climt.__file__).parent
    offenders = [
        str(path) for path in root.rglob("*.py")
        if path.name != "fortran.py" and "climt_fortran" in path.read_text()
    ]
    offenders += [
        str(path) for path in (root / "_components").rglob("component.py")
        if "from . import _" in path.read_text()
    ]
    assert offenders == []
