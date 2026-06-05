# tests/test_build_experiments.py
"""Hash-bookkeeping tests for scripts/build_experiments.py.

The regenerator walks docs/**/sources.yml, content-hashes each artifact's deps,
re-runs only stale artifacts, and (in --check mode) reports staleness without
executing. These tests build a throwaway repo tree in tmp_path.
"""
import json
import os
import subprocess
import sys
import textwrap

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
DRIVER = os.path.join(ROOT, "scripts", "build_experiments.py")
PY = sys.executable


def _make_tree(tmp_path):
    """A minimal repo: one experiment whose single artifact copies a dep file."""
    exp = tmp_path / "docs" / "experiments" / "demo"
    (exp / "_artifacts").mkdir(parents=True)
    dep = tmp_path / "src" / "input.txt"
    dep.parent.mkdir(parents=True)
    dep.write_text("v1")
    # cmd reads src/input.txt and writes the artifact png (here just text bytes).
    out_rel = "_artifacts/out.png"
    cmd = (f'{PY} -c "open(\'{exp}/_artifacts/out.png\',\'w\')'
           f'.write(open(\'{dep}\').read())"')
    (exp / "sources.yml").write_text(textwrap.dedent(f"""\
        artifacts:
          - out: {out_rel}
            out_txt: _artifacts/out.txt
            cmd: {cmd}
            deps:
              - src/input.txt
        """))
    return tmp_path, dep, exp / out_rel


def _run(root, *args):
    return subprocess.run([PY, DRIVER, "--root", str(root), *args],
                          capture_output=True, text=True)


def test_first_run_regenerates_then_clean(tmp_path):
    root, dep, out = _make_tree(tmp_path)
    r1 = _run(root)
    assert r1.returncode == 0
    assert out.exists() and out.read_text() == "v1"
    # The captured-stdout sidecar is written when out_txt is declared.
    assert (out.parent / "out.txt").exists()
    # Second run: nothing stale, exits 0, --check passes.
    r2 = _run(root, "--check")
    assert r2.returncode == 0


def test_check_flags_stale_after_dep_change(tmp_path):
    root, dep, out = _make_tree(tmp_path)
    _run(root)                      # populate
    dep.write_text("v2")            # change a dep -> stale
    rc = _run(root, "--check")
    assert rc.returncode == 1
    assert "STALE" in rc.stdout
    # Regenerate, then it's clean and the output reflects the new dep.
    _run(root)
    assert out.read_text() == "v2"
    assert _run(root, "--check").returncode == 0


def test_missing_output_is_stale(tmp_path):
    root, dep, out = _make_tree(tmp_path)
    _run(root)
    out.unlink()
    assert _run(root, "--check").returncode == 1
