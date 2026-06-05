# Sub-project B Implementation Plan — Picket-Fence LW Discovery Post (Quarto)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship the picket-fence longwave "fast vs faithful" discovery story as the first Quarto Experiments post (`docs/experiments/2026-06-05-picket-fence-co2-bands/`), with a runnable companion notebook, a robust tropopause helper, and the minimal forward-compatible Quarto/artifact scaffolding to render and regenerate it.

**Architecture:** B conforms to the Quarto website spec (`2026-05-18-climt-website-design.md`) but builds only the *minimal* scaffolding to author/render/regenerate one Experiments post; the full site migration is a separate project. New code is two pure-Python, independently-tested units — `scripts/experiments/tropopause.py` (θ-curvature + cold-point locator) and `scripts/build_experiments.py` (hash-gated artifact regenerator) — plus a single figure generator and small `--save` hooks added to existing experiment scripts. The post prose lives in `index.qmd`; figures are regenerated into `_artifacts/` via `sources.yml` + `make experiments`; the notebook is included as "Try it yourself."

**Tech Stack:** Python (numpy, matplotlib, PyYAML), pytest, Quarto (manual prereq), the existing `climt` conda env (`conda run -n climt …`), Jupyter (`nbclient`/`nbconvert`, present in the env).

**Conventions for every task:** run Python via the climt env (`conda run --no-capture-output -n climt python …`). Commit at the end of each task. Do not touch the unrelated working-tree changes already present on this branch.

---

## File Structure

```
scripts/
├── experiments/
│   ├── tropopause.py                     # NEW: θ-curvature + cold-point locator (pure numpy)
│   └── make_subproject_B_figures.py      # NEW: single figure generator (--figure {kg,lbl,rce,throughput,all})
│   └── rce_moist_pf_vs_rrtmg.py          # MODIFY: add --save <npz>
│   └── rce_dry_pf_vs_rrtmg.py            # MODIFY: add --save <npz>
│   └── bench_pf_vs_rrtmg.py              # MODIFY: add --save <npz>
│   └── plot_kg_curves.py                 # MODIFY: add --out <png>
└── build_experiments.py                  # NEW: hash-gated artifact regenerator (--check)

tests/
├── test_tropopause.py                    # NEW: clean / runaway / degenerate profiles
└── test_build_experiments.py             # NEW: hash bookkeeping + --check

docs/
├── _quarto.yml                           # NEW: minimal website project config + Experiments listing
└── experiments/
    └── 2026-06-05-picket-fence-co2-bands/
        ├── index.qmd                     # NEW: the post (narrative + callouts + figure/notebook includes)
        ├── sources.yml                   # NEW: artifact manifest
        ├── picket_fence_co2_bands.ipynb  # NEW: companion notebook (committed pre-executed)
        └── _artifacts/                   # NEW: regenerated figures + captured stdout + copied npz
references.bib                            # NEW: only the entries the post cites (each with url=)
Makefile                                  # MODIFY: add `experiments` target
```

Each new module has one responsibility: `tropopause.py` knows potential temperature + curvature, nothing about climt state; `build_experiments.py` knows hashing + subprocess, nothing about figures; `make_subproject_B_figures.py` knows plotting, reading the existing `.npz`/scripts. They share no globals.

---

## Task 1: Tropopause helper (θ-curvature + cold-point) — TDD

**Files:**
- Create: `scripts/experiments/tropopause.py`
- Test: `tests/test_tropopause.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_tropopause.py
"""Unit tests for the robust tropopause locator (scripts/experiments/tropopause.py).

Profiles are ordered surface -> top (index 0 = surface), matching the RCE
scripts' state arrays. The old constant-theta+10K marker mislocated the
tropopause for runaway-warm columns (the 583 hPa artifact); these tests pin the
replacement's behaviour.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "experiments"))
from tropopause import find_tropopause, potential_temperature  # noqa: E402


def _synthetic_column(p_tropopause_hpa, nz=40):
    """RCE-like column: constant-theta troposphere, then a sharp stratospheric
    theta ramp above the prescribed tropopause pressure. Returns (T, p_hPa)."""
    p = np.linspace(1000.0, 50.0, nz)              # surface -> top, decreasing
    theta_trop = 300.0
    # theta constant below the tropopause, ramping up steeply above it.
    theta = np.where(
        p >= p_tropopause_hpa,
        theta_trop,
        theta_trop + 2.0 * (p_tropopause_hpa - p),  # steep rise into stratosphere
    )
    T = theta * (p / 1000.0) ** (287.0 / 1004.0)
    return T, p


def test_potential_temperature_surface_identity():
    T = np.array([300.0])
    p = np.array([1000.0])
    assert potential_temperature(T, p)[0] == pytest.approx(300.0)


def test_curvature_finds_prescribed_kink():
    T, p = _synthetic_column(p_tropopause_hpa=200.0, nz=40)
    out = find_tropopause(T, p)
    # The kink is at 200 hPa; the located level must be within one grid step.
    grid_step = abs(p[1] - p[0])
    assert abs(out["p_curvature"] - 200.0) <= grid_step


def test_runaway_warm_column_is_not_mislocated_low():
    """A deep, warm, moist-adiabat-ish troposphere whose theta drifts upward
    through the column -- the case where +theta0+10K fired too low (~583 hPa).
    The curvature locator must land high (low pressure), not in the mid-troposphere."""
    nz = 40
    p = np.linspace(1010.0, 40.0, nz)
    # theta drifts +30 K across the troposphere (moist adiabat), then a hard
    # stratospheric cap above 150 hPa.
    theta_trop = 290.0 + 30.0 * (1010.0 - np.clip(p, 150.0, 1010.0)) / (1010.0 - 150.0)
    theta_strat_extra = np.where(p < 150.0, 3.0 * (150.0 - p), 0.0)
    theta = theta_trop + theta_strat_extra
    T = theta * (p / 1000.0) ** (287.0 / 1004.0)
    out = find_tropopause(T, p)
    assert out["p_curvature"] < 250.0      # high up, near the 150 hPa cap
    assert out["p_curvature"] > 80.0       # but not pinned at the model top


def test_degenerate_short_column_returns_sentinel():
    T = np.array([280.0, 270.0, 260.0])    # < 4 levels
    p = np.array([1000.0, 700.0, 400.0])
    out = find_tropopause(T, p)
    assert out["p_curvature"] == pytest.approx(p[-1])   # model-top sentinel
    assert out["k_coldpoint"] == 2                       # coldest = top here


def test_coldpoint_reported_alongside_curvature():
    T, p = _synthetic_column(p_tropopause_hpa=200.0, nz=40)
    out = find_tropopause(T, p)
    assert set(out) == {"p_curvature", "p_coldpoint", "k_curvature", "k_coldpoint"}
    assert out["p_coldpoint"] == pytest.approx(p[int(np.argmin(T))])
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `conda run --no-capture-output -n climt python -m pytest tests/test_tropopause.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'tropopause'`.

- [ ] **Step 3: Write the implementation**

```python
# scripts/experiments/tropopause.py
"""Robust cross-climate tropopause locator for the picket-fence RCE columns.

Replaces the constant-theta+10K marker (formerly in rce_moist_pf_vs_rrtmg.py)
that mislocated the tropopause for runaway-warm columns -- the 583 hPa artifact.
Pure numpy, so it is Pyodide-safe and importable from the experiment notebook.

Primary criterion: theta-curvature. The tropopause is the kink where potential
temperature transitions from the convective, near-constant-theta troposphere to
the rapidly increasing stratosphere -- the interior level of maximum |d^2 theta /
d(ln p)^2|. Secondary cross-check: the cold-point (temperature minimum), which
rides high (low pressure) in the moist tropics.
"""
import numpy as np

_RD_OVER_CP = 287.0 / 1004.0   # dry-air R/cp, matches the RCE scripts
_P0_HPA = 1000.0


def potential_temperature(temperature, pressure):
    """theta = T (p0/p)^(Rd/cp); pressure in hPa, same shape as temperature."""
    temperature = np.asarray(temperature, dtype=float)
    pressure = np.asarray(pressure, dtype=float)
    return temperature * (_P0_HPA / pressure) ** _RD_OVER_CP


def find_tropopause(temperature, pressure):
    """Locate the tropopause from one column's T(p) profile.

    Parameters
    ----------
    temperature : (nz,) array of K, ordered surface -> top (index 0 = surface).
    pressure : (nz,) array of hPa, same ordering, strictly decreasing.

    Returns
    -------
    dict:
        p_curvature : tropopause pressure (hPa), the headline value (max upward
            curvature of theta vs ln p).
        p_coldpoint : pressure (hPa) of the temperature minimum (cross-check).
        k_curvature, k_coldpoint : the corresponding integer level indices.

    Degenerate columns (fewer than 4 levels, or no usable interior curvature)
    return p_curvature = pressure[-1] (model top) as a documented sentinel.
    """
    temperature = np.asarray(temperature, dtype=float)
    pressure = np.asarray(pressure, dtype=float)
    nz = pressure.size

    k_cold = int(np.argmin(temperature))

    if nz < 4:
        return {
            "p_curvature": float(pressure[-1]),
            "p_coldpoint": float(pressure[k_cold]),
            "k_curvature": nz - 1,
            "k_coldpoint": k_cold,
        }

    theta = potential_temperature(temperature, pressure)
    lnp = np.log(pressure)

    # Second derivative on the (non-uniform) ln-p grid; np.gradient twice handles
    # unequal spacing. The stratospheric kink is the interior extremum of largest
    # magnitude; exclude the two end levels where np.gradient is one-sided.
    dtheta = np.gradient(theta, lnp)
    d2theta = np.gradient(dtheta, lnp)
    interior = np.arange(1, nz - 1)
    k_curv = int(interior[np.argmax(np.abs(d2theta[interior]))])

    return {
        "p_curvature": float(pressure[k_curv]),
        "p_coldpoint": float(pressure[k_cold]),
        "k_curvature": k_curv,
        "k_coldpoint": k_cold,
    }
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `conda run --no-capture-output -n climt python -m pytest tests/test_tropopause.py -v`
Expected: PASS (5 passed). If `test_runaway_warm_column_is_not_mislocated_low` fails, the curvature peak is being captured by a near-surface boundary-layer kink — restrict `interior` to levels above the lowest 20% of the column (`interior = np.arange(max(1, nz // 5), nz - 1)`) and re-run; this is a legitimate TDD refinement, keep the algorithm note in the docstring in sync.

- [ ] **Step 5: Commit**

```bash
git add scripts/experiments/tropopause.py tests/test_tropopause.py
git commit -m "feat(pf): robust theta-curvature tropopause locator + tests"
```

---

## Task 2: Artifact regenerator `build_experiments.py` — TDD

**Files:**
- Create: `scripts/build_experiments.py`
- Test: `tests/test_build_experiments.py`

Depends on PyYAML. Verify first: `conda run --no-capture-output -n climt python -c "import yaml; print(yaml.__version__)"`. If missing, add `pyyaml` to the env (`conda run -n climt pip install pyyaml`) and note it in the post's Contributing aside.

- [ ] **Step 1: Write the failing test**

```python
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
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `conda run --no-capture-output -n climt python -m pytest tests/test_build_experiments.py -v`
Expected: FAIL — driver does not exist (`can't open file '.../scripts/build_experiments.py'`), tests error.

- [ ] **Step 3: Write the implementation**

```python
# scripts/build_experiments.py
"""Regenerate experiment/chapter artifacts declared in docs/**/sources.yml.

Each sources.yml lists artifacts; each artifact maps an output figure (and
optional captured-stdout file) to the command that produces it and the source
files whose content determines its validity. Deps are content-hashed (sha256,
glob-expanded) and compared to a stored _artifacts/.hashes.json; only stale
artifacts (output missing or dep-hash changed) are re-run. Hashes, not mtimes,
so branch switches and network filesystems don't trigger spurious rebuilds.

    python scripts/build_experiments.py            # regenerate stale artifacts
    python scripts/build_experiments.py --check    # CI: report + exit 1 if stale
"""
import argparse
import glob
import hashlib
import json
import os
import subprocess
import sys

import yaml


def _expand_deps(repo_root, deps):
    paths = []
    for pattern in deps:
        matches = sorted(glob.glob(os.path.join(repo_root, pattern), recursive=True))
        if not matches:
            raise FileNotFoundError(f"dep pattern matched no files: {pattern}")
        paths.extend(matches)
    return paths


def _hash_files(paths):
    h = hashlib.sha256()
    for p in sorted(set(paths)):
        h.update(p.encode())
        with open(p, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 16), b""):
                h.update(chunk)
    return h.hexdigest()


def _load_hashes(path):
    return json.load(open(path)) if os.path.exists(path) else {}


def _stale(repo_root, manifest_dir, artifact, stored):
    out = os.path.join(manifest_dir, artifact["out"])
    current = _hash_files(_expand_deps(repo_root, artifact["deps"]))
    if not os.path.exists(out):
        return True, current
    return (stored.get(artifact["out"]) != current), current


def _regenerate_manifest(repo_root, manifest_path, check_only):
    manifest_dir = os.path.dirname(manifest_path)
    manifest = yaml.safe_load(open(manifest_path))
    hashes_path = os.path.join(manifest_dir, "_artifacts", ".hashes.json")
    os.makedirs(os.path.dirname(hashes_path), exist_ok=True)
    stored = _load_hashes(hashes_path)
    any_stale = False
    for art in manifest["artifacts"]:
        stale, current = _stale(repo_root, manifest_dir, art, stored)
        if not stale:
            continue
        any_stale = True
        if check_only:
            print(f"STALE: {art['out']} (dep changed or output missing)")
            continue
        print(f"regenerating {art['out']} ...")
        done = subprocess.run(art["cmd"], shell=True, cwd=repo_root,
                              capture_output=True, text=True)
        if done.returncode != 0:
            sys.stderr.write(done.stdout + done.stderr)
            raise RuntimeError(f"cmd failed for {art['out']}: {art['cmd']}")
        if "out_txt" in art:
            with open(os.path.join(manifest_dir, art["out_txt"]), "w") as f:
                f.write(done.stdout)
        if not os.path.exists(os.path.join(manifest_dir, art["out"])):
            raise RuntimeError(
                f"cmd for {art['out']} did not write the output file")
        stored[art["out"]] = current
    if not check_only:
        with open(hashes_path, "w") as f:
            json.dump(stored, f, indent=2, sort_keys=True)
    return any_stale


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="report staleness and exit 1 if any; no execution")
    ap.add_argument("--root", default=os.getcwd(), help="repo root (default cwd)")
    args = ap.parse_args()
    manifests = sorted(glob.glob(
        os.path.join(args.root, "docs", "**", "sources.yml"), recursive=True))
    any_stale = False
    for m in manifests:
        any_stale |= _regenerate_manifest(args.root, m, args.check)
    if args.check and any_stale:
        sys.stderr.write("Stale artifacts. Run `make experiments` and commit.\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `conda run --no-capture-output -n climt python -m pytest tests/test_build_experiments.py -v`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/build_experiments.py tests/test_build_experiments.py
git commit -m "feat(docs): hash-gated experiment-artifact regenerator + tests"
```

---

## Task 3: `make experiments` target + minimal Quarto scaffolding

**Files:**
- Modify: `Makefile` (append a target; do not restructure the cookiecutter file)
- Create: `docs/_quarto.yml`
- Create: `references.bib`

- [ ] **Step 1: Add the `experiments` target to the Makefile**

Append (the file uses tab-indented recipes — use a real tab):

```makefile
experiments: ## regenerate stale docs experiment artifacts (hash-gated)
	conda run --no-capture-output -n climt python scripts/build_experiments.py

experiments-check: ## fail if any experiment artifact is stale (CI bookkeeping)
	conda run --no-capture-output -n climt python scripts/build_experiments.py --check
```

- [ ] **Step 2: Verify the target is registered**

Run: `make help | grep experiments`
Expected: both `experiments` and `experiments-check` lines print with their `##` help text.

- [ ] **Step 3: Write the minimal Quarto project config**

```yaml
# docs/_quarto.yml
# Minimal, forward-compatible seed of the full climt-website config
# (docs/superpowers/specs/2026-05-18-climt-website-design.md). The site
# migration extends this with nav, theme, and quartodoc; nothing here may
# contradict that spec.
project:
  type: website
  output-dir: _site

website:
  title: "climt"
  navbar:
    left:
      - text: "Experiments"
        href: experiments/index.qmd

format:
  html:
    theme: cosmo
    toc: true
    code-fold: false

# Experiments listing is provided by docs/experiments/index.qmd (Task 9).
```

- [ ] **Step 4: Write the citations file (only entries the post cites)**

```bibtex
% references.bib -- only entries cited by the picket-fence post. Each carries a
% url=; the website migration appends further entries.
@article{parmentier2014,
  author  = {Parmentier, V. and Guillot, T.},
  title   = {A non-grey analytical model for irradiated atmospheres},
  journal = {Astronomy \& Astrophysics},
  volume  = {562},
  pages   = {A133},
  year    = {2014},
  doi     = {10.1051/0004-6361/201322342},
  url     = {https://doi.org/10.1051/0004-6361/201322342}
}
@article{mlawer1997,
  author  = {Mlawer, E. J. and Taubman, S. J. and Brown, P. D. and Iacono, M. J.
             and Clough, S. A.},
  title   = {Radiative transfer for inhomogeneous atmospheres: RRTM, a validated
             correlated-k model for the longwave},
  journal = {Journal of Geophysical Research},
  volume  = {102},
  number  = {D14},
  pages   = {16663--16682},
  year    = {1997},
  doi     = {10.1029/97JD00237},
  url     = {https://doi.org/10.1029/97JD00237}
}
@article{mlawer2012,
  author  = {Mlawer, E. J. and Payne, V. H. and Moncet, J.-L. and Delamere, J. S.
             and Alvarado, M. J. and Tobin, D. C.},
  title   = {Development and recent evaluation of the MT\_CKD model of continuum
             absorption},
  journal = {Philosophical Transactions of the Royal Society A},
  volume  = {370},
  pages   = {2520--2556},
  year    = {2012},
  doi     = {10.1098/rsta.2011.0295},
  url     = {https://doi.org/10.1098/rsta.2011.0295}
}
```

Add further `@key` entries only when the prose actually cites them in Task 9 (e.g. a WMO tropopause reference); never pre-populate unused entries.

- [ ] **Step 5: Commit**

```bash
git add Makefile docs/_quarto.yml references.bib
git commit -m "build(docs): minimal Quarto scaffolding + make experiments target"
```

---

## Task 4: Add `--save`/`--out` hooks to existing experiment scripts

The figure generator (Task 5) consumes clean `.npz`/`.png` from the existing scripts. Add small, focused output hooks. Each sub-step is a minimal edit + a smoke check.

**Files:**
- Modify: `scripts/experiments/rce_moist_pf_vs_rrtmg.py`
- Modify: `scripts/experiments/rce_dry_pf_vs_rrtmg.py`
- Modify: `scripts/experiments/bench_pf_vs_rrtmg.py`
- Modify: `scripts/experiments/plot_kg_curves.py`

- [ ] **Step 1: `rce_moist` — add `--save` writing final-state arrays per column**

In the argparse block (near the existing `--table` / `--co2` args), add:

```python
ap.add_argument("--save", default=None,
                help="if given, write final-state arrays to this .npz "
                     "(keys: <label>__{p_hpa,T,q,hr_lw,T_sfc,olr})")
```

After the existing `np.savez(npz, ...)` report block (the script already builds `history`), add — replacing the **fixed** `npz` path logic with an `--save`-aware write (keep the existing default-path save for back-compat):

```python
if args.save:
    payload = {}
    for label in COLUMNS:
        s = history[label][-1]
        safe = label.replace(" ", "_")
        payload[f"{safe}__p_hpa"] = np.asarray(s["p_hpa"])
        payload[f"{safe}__T"] = np.asarray(s["T"])
        payload[f"{safe}__q"] = np.asarray(s["q"])
        payload[f"{safe}__hr_lw"] = np.asarray(s["hr_lw"])
        payload[f"{safe}__T_sfc"] = np.asarray(s["T_sfc"])
        payload[f"{safe}__olr"] = np.asarray(s["olr"])
    np.savez(args.save, **payload)
    print(f"saved final state -> {args.save}")
```

- [ ] **Step 2: `rce_dry` — same `--save` addition**

Apply the identical `--save` argparse line and the identical save block to `rce_dry_pf_vs_rrtmg.py` (it has the same `COLUMNS`/`history`/state-dict shape). Repeat the code rather than importing across scripts — they are standalone experiment scripts.

- [ ] **Step 3: `bench` — add `--save` writing the two throughput scalars**

In `bench_pf_vs_rrtmg.py`'s `main()`, after the two timing prints (`mean_r`, `mean_p` known), add an argparse `--save` and:

```python
if args.save:
    np.savez(args.save,
             rrtmg_us_per_col=np.asarray(mean_r / NCOL * 1e3),
             pf_us_per_col=np.asarray(mean_p / NCOL * 1e3),
             ncol=np.asarray(NCOL))
    print(f"saved throughput -> {args.save}")
```

(Add `import argparse` + an `ap` in `main()` if not present; `import numpy as np` is already there.)

- [ ] **Step 4: `plot_kg_curves` — add `--out` to override the figure path**

Add an argparse `--out` defaulting to the script's current hard-coded
`debug_data/kg_curves_window_band.png`, and use it in the `savefig` call, so the
generator can write straight into `_artifacts/`.

- [ ] **Step 5: Smoke-check each hook is wired (no full runs yet)**

Run: `conda run --no-capture-output -n climt python scripts/experiments/rce_moist_pf_vs_rrtmg.py --help`
Expected: `--save` appears in the help. Repeat `--help` for `rce_dry`, `bench`, and confirm `--out` for `plot_kg_curves`.

- [ ] **Step 6: Commit**

```bash
git add scripts/experiments/rce_moist_pf_vs_rrtmg.py scripts/experiments/rce_dry_pf_vs_rrtmg.py scripts/experiments/bench_pf_vs_rrtmg.py scripts/experiments/plot_kg_curves.py
git commit -m "feat(pf): add --save/--out hooks to RCE/bench/kg scripts for figure pipeline"
```

---

## Task 5: Figure generator `make_subproject_B_figures.py`

**Files:**
- Create: `scripts/experiments/make_subproject_B_figures.py`

One command, `--figure {kg,lbl,rce,throughput,all}` + `--out <png>` (+ `--npz` where a figure reads a produced `.npz`). Each figure is one function. The `kg`, `lbl`, and `throughput` figures are self-contained; the `rce` figure reads two `.npz` produced by `rce_moist --save` (driven from `sources.yml`, Task 6).

- [ ] **Step 1: Write the generator**

```python
# scripts/experiments/make_subproject_B_figures.py
"""Single entry point for the picket-fence Experiments-post hero figures.

    --figure kg          window-band k(g) shelf+cliff (dry vs moist)
    --figure lbl         PF vs RRTMG vs line-by-line OLR overlay
    --figure rce         before/after moist RCE temperature + tropopause markers
    --figure throughput  PF vs RRTMG throughput bar (us/col)

Reads existing committed data (debug_data/*.npz) and the .npz produced by the
RCE/bench --save hooks; writes one PNG to --out. Wired into the experiment's
sources.yml so `make experiments` regenerates stale figures only.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DEBUG = os.path.join(REPO, "debug_data")


def fig_kg(out):
    """Delegate to plot_kg_curves.py --out (it owns the k(g) construction)."""
    subprocess.run(
        [sys.executable, os.path.join(HERE, "plot_kg_curves.py"), "--out", out],
        cwd=REPO, check=True)


def fig_lbl(out):
    """PF/RRTMG scalar OLR vs line-by-line spectrum, moist + dry at 400 ppm."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, kind in zip(axes, ("moist", "dry")):
        d = np.load(os.path.join(DEBUG, f"lbl_olr_spec_{kind}_co2_400ppm.npz"))
        ax.plot(d["nu"], d["olr_spec"], lw=0.4, color="#444",
                label="line-by-line")
        ax.set_title(f"{kind}: LBL={float(d['total']):.1f}  "
                     f"PF={float(d['olr_pf']):.1f}  "
                     f"RRTMG={float(d['olr_rrtmg']):.1f} W/m$^2$")
        ax.set_xlabel("wavenumber (cm$^{-1}$)")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("spectral OLR (W/m$^2$/cm$^{-1}$)")
    fig.suptitle("Line-by-line truth vs PF and RRTMG (CO$_2$=400 ppm)")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out, dpi=140)


def _load_rce(npz_path):
    d = np.load(npz_path)
    cols = sorted({k.split("__")[0] for k in d.files})
    return {c: {k.split("__")[1]: d[k] for k in d.files if k.startswith(c + "__")}
            for c in cols}


def fig_rce(out, before_npz, after_npz):
    """Before/after moist RCE temperature profiles with tropopause markers."""
    sys.path.insert(0, HERE)
    from tropopause import find_tropopause
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    for ax, npz, title in (
        (axes[0], before_npz, "before: 4-band default"),
        (axes[1], after_npz, "after: 14-band CO$_2$ table"),
    ):
        cols = _load_rce(npz)
        for label, s in cols.items():
            ax.plot(s["T"], s["p_hpa"], "-o", ms=3, label=label)
            tp = find_tropopause(s["T"], s["p_hpa"])
            ax.plot(np.interp(tp["p_curvature"], s["p_hpa"][::-1], s["T"][::-1]),
                    tp["p_curvature"], "*", ms=16, mec="k")
        ax.set_yscale("log"); ax.invert_yaxis(); ax.grid(alpha=0.3)
        ax.set_title(title); ax.set_xlabel("Temperature (K)")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("Pressure (hPa)")
    fig.suptitle("Moist RCE: PF vs RRTMG (★ = θ-curvature tropopause)")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out, dpi=140)


def fig_throughput(out, npz):
    d = np.load(npz)
    fig, ax = plt.subplots(figsize=(5, 5))
    vals = [float(d["rrtmg_us_per_col"]), float(d["pf_us_per_col"])]
    ax.bar(["RRTMG-LW", "PF-LW\n(14b×8g, njit)"], vals,
           color=["#888", "#1f77b4"])
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:.1f}", ha="center", va="bottom")
    ax.set_ylabel("µs / column"); ax.set_title("Longwave throughput (lower is better)")
    fig.tight_layout()
    fig.savefig(out, dpi=140)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", required=True,
                    choices=["kg", "lbl", "rce", "throughput", "all"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--before-npz")
    ap.add_argument("--after-npz")
    ap.add_argument("--npz")
    a = ap.parse_args()
    if a.figure == "kg":
        fig_kg(a.out)
    elif a.figure == "lbl":
        fig_lbl(a.out)
    elif a.figure == "rce":
        fig_rce(a.out, a.before_npz, a.after_npz)
    elif a.figure == "throughput":
        fig_throughput(a.out, a.npz)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke-test the two self-contained figures**

Run:
```bash
conda run --no-capture-output -n climt python scripts/experiments/make_subproject_B_figures.py --figure lbl --out /tmp/lbl.png
conda run --no-capture-output -n climt python scripts/experiments/make_subproject_B_figures.py --figure kg --out /tmp/kg.png
```
Expected: both print `wrote …` and the PNGs exist (`ls -l /tmp/lbl.png /tmp/kg.png`). `rce`/`throughput` are exercised in Task 6 once their `.npz` exist.

- [ ] **Step 3: Commit**

```bash
git add scripts/experiments/make_subproject_B_figures.py
git commit -m "feat(pf): single figure generator for the discovery post"
```

---

## Task 6: Experiment `sources.yml` + first regeneration

**Files:**
- Create: `docs/experiments/2026-06-05-picket-fence-co2-bands/sources.yml`
- Generated: `docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/*` (figures, captured stdout, `.hashes.json`)

The manifest lists the two RCE `.npz` (slow, deps = climt source + tables) **before** the figures that consume them, and the throughput `.npz`, so the in-order driver produces inputs before dependents.

- [ ] **Step 1: Write the manifest**

```yaml
# docs/experiments/2026-06-05-picket-fence-co2-bands/sources.yml
# Paths in `out`/`out_txt` are relative to this file's dir; cmd/deps paths are
# relative to the repo root. The driver runs entries in order.
artifacts:
  # --- data: two moist-RCE final states (slow; rerun only if PF source/tables change)
  - out: _artifacts/rce_moist_before.npz
    cmd: conda run --no-capture-output -n climt python scripts/experiments/rce_moist_pf_vs_rrtmg.py --table earth_low_res_lw_4band_ngpt2_before.nc --days 200 --save docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/rce_moist_before.npz
    deps:
      - scripts/experiments/rce_moist_pf_vs_rrtmg.py
      - climt/_components/picket_fence/**/*.py
      - climt/_data/picket_fence/correlated_k/earth_low_res_lw_4band_ngpt2_before.nc
  - out: _artifacts/rce_moist_after.npz
    cmd: conda run --no-capture-output -n climt python scripts/experiments/rce_moist_pf_vs_rrtmg.py --table earth_low_res_lw.nc --days 200 --save docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/rce_moist_after.npz
    deps:
      - scripts/experiments/rce_moist_pf_vs_rrtmg.py
      - climt/_components/picket_fence/**/*.py
      - climt/_data/picket_fence/correlated_k/earth_low_res_lw.nc
  # --- data: throughput scalars
  - out: _artifacts/throughput.npz
    out_txt: _artifacts/throughput.txt
    cmd: conda run --no-capture-output -n climt python scripts/experiments/bench_pf_vs_rrtmg.py --save docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/throughput.npz
    deps:
      - scripts/experiments/bench_pf_vs_rrtmg.py
      - climt/_components/picket_fence/**/*.py
  # --- figures
  - out: _artifacts/kg_window_band.png
    cmd: conda run --no-capture-output -n climt python scripts/experiments/make_subproject_B_figures.py --figure kg --out docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/kg_window_band.png
    deps:
      - scripts/experiments/make_subproject_B_figures.py
      - scripts/experiments/plot_kg_curves.py
  - out: _artifacts/lbl_overlay.png
    cmd: conda run --no-capture-output -n climt python scripts/experiments/make_subproject_B_figures.py --figure lbl --out docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/lbl_overlay.png
    deps:
      - scripts/experiments/make_subproject_B_figures.py
      - debug_data/lbl_olr_spec_moist_co2_400ppm.npz
      - debug_data/lbl_olr_spec_dry_co2_400ppm.npz
  - out: _artifacts/rce_before_after.png
    cmd: conda run --no-capture-output -n climt python scripts/experiments/make_subproject_B_figures.py --figure rce --before-npz docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/rce_moist_before.npz --after-npz docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/rce_moist_after.npz --out docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/rce_before_after.png
    deps:
      - scripts/experiments/make_subproject_B_figures.py
      - scripts/experiments/tropopause.py
      - docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/rce_moist_before.npz
      - docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/rce_moist_after.npz
  - out: _artifacts/throughput.png
    cmd: conda run --no-capture-output -n climt python scripts/experiments/make_subproject_B_figures.py --figure throughput --npz docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/throughput.npz --out docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/throughput.png
    deps:
      - scripts/experiments/make_subproject_B_figures.py
      - docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/throughput.npz
```

- [ ] **Step 2: Regenerate all artifacts (slow: the two RCE runs dominate)**

Run: `make experiments`
Expected: each `regenerating …` line prints; ends without error. The two RCE entries take the longest (200-day integrations). Confirm outputs:
`ls docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/` shows `kg_window_band.png lbl_overlay.png rce_before_after.png throughput.png throughput.npz throughput.txt rce_moist_before.npz rce_moist_after.npz .hashes.json`.

- [ ] **Step 3: Verify idempotence (`--check` is clean immediately after)**

Run: `make experiments-check`
Expected: exit 0, no `STALE` lines.

- [ ] **Step 4: Copy the LBL `.npz` into `_artifacts/` for the notebook's offline path**

```bash
cp debug_data/lbl_olr_spec_moist_co2_400ppm.npz debug_data/lbl_olr_spec_dry_co2_400ppm.npz \
   docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/
```
(The notebook's RRTMG/LBL cells read these, so the post is self-contained for a future browser port.)

- [ ] **Step 5: Commit**

```bash
git add docs/experiments/2026-06-05-picket-fence-co2-bands/sources.yml docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/
git commit -m "feat(pf): experiment sources.yml + regenerated hero artifacts"
```

---

## Task 7: Companion notebook (authored, executed, committed pre-executed)

**Files:**
- Create: `docs/experiments/2026-06-05-picket-fence-co2-bands/picket_fence_co2_bands.ipynb`

Author the notebook cell-by-cell to the spec below, then execute it in place. Each code cell mirrors the existing scripts / figure generator (reuse, don't reinvent). Markdown cells carry one-line framing only — the full narrative lives in `index.qmd`.

- [ ] **Step 1: Author the cells**

Cell list (md = markdown, py = code):

1. **md** — Title + "Try it yourself: reproduce every figure in the post."
2. **py — setup.** Imports; `sys.path.insert` for `tropopause`; define `ARTIF = "_artifacts"`; construct `PicketFenceLongwave(table="earth_low_res_lw")` and the "before" `PicketFenceLongwave(table="earth_low_res_lw_4band_ngpt2_before.nc")`; `RRTMG_AVAILABLE` try/except around `from climt import RRTMGLongwave`.
3. **md** — "The symptom."
4. **py — symptom.** Load `_artifacts/rce_moist_before.npz` via the `_load_rce` pattern; print T_sfc, q_sfc, and `find_tropopause(...)['p_curvature']`; assert the headline numbers (`T_sfc` warm, tropopause ≈ 583 hPa) so the cell self-checks.
5. **md** — "k-distribution sandbox."
6. **py — k(g).** Build and plot a k(g) curve for the window band (reuse `plot_kg_curves` logic or read its arrays); show the dry-vs-moist shelf+cliff. Pure-Python (Pyodide-live).
7. **md** — "Isolating the culprit."
8. **py — LBL overlay.** Load `_artifacts/lbl_olr_spec_moist_co2_400ppm.npz`; overlay PF vs RRTMG scalar OLR on the LBL spectrum; print the three OLR numbers.
9. **md** — "The dead-ends."
10. **py — dead-ends.** Show (from committed forward-flux artifacts or a short inline run) that 8→16→32 g-points barely move the bias; comment that finer vertical resolution behaves the same.
11. **md** — "Before → after."
12. **py — before/after.** Load both `rce_moist_*.npz`; plot side by side with `find_tropopause` markers; print the +11.9 K → +1.9 K and 583 → 196 hPa contrast.
13. **md** — "Which diagnostic to trust."
14. **py — subtlety.** Print fixed-profile single-column LBL−PF (+11 for both) next to the converged RCE ΔT (~2 K).
15. **md** — "Two more things: fast & the free knob."
16. **py — coda.** Load `_artifacts/throughput.npz` (PF vs RRTMG µs/col); a CO₂-sweep cell calling PF at a few CO₂ values (pure-Python, Pyodide-live).
17. **md** — "Pyodide note": which cells above are pure-Python-live (6, 16 CO₂ sweep, PF construction) vs which need pre-baked `.npz` (4, 8, 12, 14). This is the cell-boundary map the live-cells follow-up consumes.

- [ ] **Step 2: Execute the notebook end-to-end in the climt env**

Run:
```bash
cd docs/experiments/2026-06-05-picket-fence-co2-bands
conda run --no-capture-output -n climt jupyter nbconvert --to notebook --execute --inplace picket_fence_co2_bands.ipynb
```
Expected: exits 0, no cell raises. (The asserts in cells 4/12 guard the headline numbers.)

- [ ] **Step 3: Commit the pre-executed notebook**

```bash
git add docs/experiments/2026-06-05-picket-fence-co2-bands/picket_fence_co2_bands.ipynb
git commit -m "feat(pf): companion notebook (reproduces every post figure)"
```

---

## Task 8: The post `index.qmd` + Experiments listing

**Files:**
- Create: `docs/experiments/2026-06-05-picket-fence-co2-bands/index.qmd`
- Create: `docs/experiments/index.qmd` (the listing page)

Discovery-driven narrative, advanced-undergrad voice, layered via callouts. Every figure include points at `_artifacts/`; every hard number matches the source table below. Code excerpts use Quarto `include` from live source (not retyped). Citations use `@parmentier2014` / `@mlawer1997` / `@mlawer2012` against `references.bib`.

- [ ] **Step 1: Write the front-matter + the Experiments listing page**

`docs/experiments/index.qmd`:

```markdown
---
title: "Experiments"
listing:
  contents: "*/index.qmd"
  type: default
  sort: "date desc"
  fields: [date, title, description]
---
```

`index.qmd` front-matter:

```markdown
---
title: "What it takes to build a longwave scheme that rivals RRTMG"
description: "A picket-fence correlated-k warm bias, root-caused and fixed."
date: 2026-06-05
author: "Joy Monteiro"
bibliography: ../../../references.bib
categories: [radiation, picket-fence, correlated-k]
---
```

- [ ] **Step 2: Write the body to the section spec**

Sections (each beat from the spec's narrative spine; numbers from the source table):

1. **The symptom** — moist RCE before-table: **+11.9 K** T_sfc vs RRTMG, q_sfc **6.49 vs 2.02 g/kg** (3.2× runaway), tropopause **583 hPa (broken)**. Embed `_artifacts/rce_before_after.png` (left panel referenced here, full contrast revisited in §7).
2. **Primer** — `callout-note`: band-mean transmission can't use a mean k; the k-distribution trick; g-points; the "correlated" assumption. Adapted from `docs/notes/corr-k-explainer-seed.md`. Forward-link placeholder comment: once RT-walkthrough Ch 3/4 exist, replace deep math with `@`-links.
3. **Isolating the culprit** — linepyline matches RRTMG; PF over-traps ~17 (dry)/~36 (moist) W/m². Embed `_artifacts/lbl_overlay.png`.
4. **The dead-ends** — `callout-warning`: 8→16→32 g-points don't help; finer vertical res doesn't help; wing-boost rejected. The spectral-correlation floor.
5. **Root cause #1: the continuum** — MT_CKD folded in + linearly interpolated over-traps ~17 W/m²; fix = decouple (line-only k-dist + band-grey continuum, log-X), as RRTMG's `tauself`/`taufor` [@mlawer2012; @mlawer1997]. `callout-tip` for the engineering detail.
6. **Root cause #2: band structure** — the window band lumps two populations; embed `_artifacts/kg_window_band.png` (shelf+cliff); fix = 14-band partition.
7. **The payoff + diagnostic-trust lesson** — before→after table (T_sfc +11.9→+1.9–2.1 K, dry +0.8 K, q 2.2 vs 2.0, tropopause 583→196 hPa). Embed `_artifacts/rce_before_after.png` (full). `callout-note`: fixed-profile LBL−PF reads +11 for both prototype and shipped table yet RCE converges to ~2 K — the fixed-profile metric overstates a self-adjusting column's error.
8. **Two more things** — fast-vs-faithful (`@njit` → ~540×, **37.9 vs 45.2 µs/col**, embed `_artifacts/throughput.png`) and the free CO₂ knob (log-k, fidelity-free; LOO 5.5% vs 30.9%).
9. **Answer key** — the 14-band recipe (edges, 8 g-pts, decoupled continuum, runtime CO₂ axis, interp scheme).
10. **Try it yourself** — Quarto include of the notebook:

```markdown
{{< embed picket_fence_co2_bands.ipynb echo=true >}}
```
plus a download link to the raw `.ipynb`.

Add any new citation (e.g. a WMO tropopause reference if §2 cites one) to `references.bib` with a `url=`.

- [ ] **Step 3: Prose self-check**

Re-read `index.qmd`: every `_artifacts/…png` referenced exists (Task 6); every number matches the source table; every `@key` resolves in `references.bib`; no "TODO"/placeholder text remains. Fix inline.

- [ ] **Step 4: Commit**

```bash
git add docs/experiments/index.qmd docs/experiments/2026-06-05-picket-fence-co2-bands/index.qmd references.bib
git commit -m "docs(pf): the picket-fence discovery post (first Experiments post)"
```

---

## Task 9: Render verification, cross-references, and final wiring

**Files:**
- Modify (one line): `docs/superpowers/specs/2026-03-29-pyodide-pure-python-support-design.md`

- [ ] **Step 1: Full regeneration + check from a clean state (hard gate)**

Run: `make experiments && make experiments-check`
Expected: regeneration completes; `--check` exits 0. This proves the pipeline reproduces the committed artifacts.

- [ ] **Step 2: Run the full test suite for the new code (hard gate)**

Run: `conda run --no-capture-output -n climt python -m pytest tests/test_tropopause.py tests/test_build_experiments.py -v`
Expected: all pass.

- [ ] **Step 3: Quarto render (soft gate — manual prereq)**

If Quarto is installed (`quarto --version`): `quarto render docs/` and confirm the post + listing render with figures, the embedded notebook, and resolved citations (no `???` cross-refs). If Quarto is **not** installed, record this as a known manual step (it ships with the full website migration); do not block B on it. Note the result either way.

- [ ] **Step 4: One-line forward pointer in the Pyodide spec**

Add a single cross-reference line under that spec's deferred/notes area:

```markdown
- The picket-fence Experiments post (`docs/experiments/2026-06-05-picket-fence-co2-bands/`)
  is a target for the website spec's deferred live-cells follow-up; its notebook
  already maps pure-Python-live vs pre-baked-`.npz` cells.
```

- [ ] **Step 5: Update graphify + commit**

```bash
graphify update . && conda run --no-capture-output -n climt python scripts/augment_graph.py
git add -A docs/superpowers/specs/2026-03-29-pyodide-pure-python-support-design.md graphify-out
git commit -m "docs(pf): cross-link pyodide spec + refresh knowledge graph"
```

- [ ] **Step 6: Finish the branch**

Invoke `superpowers:finishing-a-development-branch` to choose merge / PR / cleanup.

---

## Self-Review (author checklist, completed at plan-writing time)

**Spec coverage:**
- Quarto Experiments post → Task 8. Companion notebook (included) → Task 7. Tropopause helper → Task 1. Figure/artifact pipeline (`sources.yml`, `build_experiments.py`, `make experiments`) → Tasks 2, 3, 6. Minimal Quarto scaffolding (`_quarto.yml`, `references.bib`) → Task 3. Single figure generator → Task 5. `--save`/`--out` hooks → Task 4. Pyodide deferral via website mechanism → Task 9 Step 4. Non-goals (no full migration, no walkthrough chapters, no CI) → respected (nothing in the plan builds them).
- All four hero figures (k(g), LBL overlay, before/after RCE, throughput) → Task 5 + Task 6 manifest.
- Hard gates (tropopause tests, build_experiments tests, notebook execution, `make experiments`) → Tasks 1, 2, 7, 9. Soft gate (quarto render) → Task 9 Step 3.

**Placeholder scan:** prose/notebook tasks specify exact sections, callout types, figure files, and numbers (not "write the content"); the only intentional forward-comment is the RT-walkthrough `@`-link note, which is conditional on chapters that are out of scope.

**Type/name consistency:** `find_tropopause` returns `{p_curvature, p_coldpoint, k_curvature, k_coldpoint}` — used identically in Task 1 tests, Task 5 `fig_rce`, and Task 7 cells. `make_subproject_B_figures.py --figure {kg,lbl,rce,throughput}` flags match the `sources.yml` cmds in Task 6. `--save`/`--out` hook names (Task 4) match their `sources.yml` invocations (Task 6). Artifact filenames in `sources.yml` match the `index.qmd` includes (Task 8).
