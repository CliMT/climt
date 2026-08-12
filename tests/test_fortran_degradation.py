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
