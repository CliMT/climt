# Modelling Tour of the Climate System — Radiation Tranche — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a new `docs/modelling-tour/` section — six in-browser pages that supply the computational half of the EC2213 course notes, using CORK's per-band diagnostics on prescribed profiles, with no time integration anywhere.

**Architecture:** Each page is a `format: live-html` Quarto page whose computational core lives in importable, natively-tested Python under `docs/modelling-tour/_tour/`. Page cells stay thin: build state → call a helper → plot. Two spectral tables serve distinct jobs — the shipped 14×8 table for anything quantitative or swept, and a new ~100-band table, served as a same-origin site asset, for single-call spectra.

**Tech Stack:** Python, numpy, matplotlib, sympl (`UnytBackend`), climt, Quarto, `quarto-live` (Pyodide + micropip), pytest.

**Spec:** `docs/superpowers/specs/2026-08-12-modelling-tour-radiation-design.md` is the source of truth. Where this plan and the spec disagree, stop and reconcile.

## Global Constraints

- **Branch:** all work continues on `feature/modelling-tour-radiation`. Do not branch off `develop`.
- **Conda env:** run every Python/pytest command in the `climt` conda env, e.g. `conda run -n climt python -m pytest ...`.
- **No time integration.** No `AdamsBashforth`, no stepping loops, no `Stepper.__call__` in a loop, on any page in this tranche. Single component calls and parameter sweeps only.
- **Backend:** `sympl.set_backend(climt.UnytBackend())` on every page and in every test.
- **Browser components only** on pages: `CorkLongwaveRadiation`, `CorkShortwaveRadiation`. `RRTMGLongwave`, `RRTMGShortwave`, `EmanuelConvection`, `SimplePhysics`, `BergerSolarInsolation` and `DcmipInitialConditions` are compiled and absent under Pyodide — never import them in a page.
- **Test marker:** anything over ~5 s gets `@pytest.mark.slow` (see `tox.ini`); CI runs `-m "not slow"`.
- **Commit after every task.** End commit messages with `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`.
- **Diffusivity convention:** the course notes' `τ_∞` is the **diffusivity-scaled** column optical depth, i.e. `τ_notes = D · τ_normal`. Any gray table calibration must use the same `D` the component runs with. Getting this wrong silently produces a wrong OLR and a non-zero heating rate.
- **Units:** `specific_humidity` in kg/kg; all other gases as mole fractions (mol/mol), e.g. `mole_fraction_of_carbon_dioxide_in_air`.
- **Experiment-artifact gate:** editing any file under `climt/_components/cork/` invalidates content-hashed docs artifacts and turns the docs workflow red. Only Task 1 does this. Regenerate with `make experiments ONLY='<glob>'` before committing — see Task 1 Step 9, which covers which artifacts are affected and what must not change. Verify with `conda run -n climt python scripts/build_experiments.py --check`.

## External prerequisite — LANDED (NOT a task here)

`load_k_table` in `climt/_components/cork/optics/correlated_k.py` used to return numpy's
lazy `NpzFile`, which re-decompresses the array from its zip on every `__getitem__`.
It now materialises every array into a plain dict at load time (`_load_npz_table`),
so both the `.npz` and `.nc` backends return the same type and the table can be indexed
in a hot loop for free.

Merged in PR #224 as `85ac334` on `develop`, and merged into this branch. **Do not
re-implement it**, and do not treat the performance figures in this plan as provisional:
they assumed this fix, and it is in place. Measured on the reference machine at nz=18,
single column: `CorkLongwaveRadiation` 57.4 → 3.6 ms/call with numba, and 89.5 → 36.8
with `NUMBA_DISABLE_JIT=1`, which is the Pyodide path the browser cost model is built on.

Sanity-check before starting, in case you are working from a stale base:

```bash
conda run -n climt python -c "
from climt._components.cork.optics.correlated_k import load_k_table
t = load_k_table('earth_low_res_lw')
print('materialised' if isinstance(t, dict) else f'NOT YET: {type(t).__name__}')"
```

Expected: `materialised`. If it says `NOT YET`, your branch predates the merge — run
`git merge origin/develop` rather than proceeding, because the quoted browser costs will
be ~2.5× off.

## File structure

| File | Responsibility |
|---|---|
| `climt/_components/cork/lw/kernels.py` | *modify* — thread `diffusivity_factor` through the transport kernel |
| `climt/_components/cork/lw/component.py` | *modify* — `diffusivity_factor` constructor arg |
| `scripts/generate_gray_default_table.py` | *modify* — `--diffusivity` flag; resolve own output by path |
| `scripts/generate_tour_spectrum_table.py` | *create* — wrapper generating the ~100-band table |
| `docs/modelling-tour/_tour/soundings.py` | prescribed profile builders (lapse-rate, analytic gray) |
| `docs/modelling-tour/_tour/spectra.py` | Planck, brightness temperature, band aggregation, weighting function |
| `docs/modelling-tour/_tour/tables.py` | locate/fetch the hi-res table, native and browser, with fallback |
| `docs/modelling-tour/index.qmd` | section landing page, course-chapter map |
| `docs/modelling-tour/0{1..6}-*.qmd` | the six pages |
| `docs/_quarto.yml` | *modify* — navbar, sidebar, render list |
| `tests/test_modelling_tour.py` | the physics each page claims, checked natively |
| `tests/test_spectrum_table.py` | validation gate for the new table |

Tasks 1–3 are library/tooling changes. Tasks 4–5 build the tested helper layer. Task 6 scaffolds the section. Tasks 7–12 are the pages, each shippable alone. Task 13 adds the hi-res table and switches pages 1–3 onto it. Task 14 closes out.

---

## Task 1: `diffusivity_factor` constructor argument on `CorkLongwaveRadiation`

The course notes derive the two-stream equations with a diffusivity factor of **2**; CORK hardcodes the Elsasser value **1.66** in two places. Page 4 needs `D = 2` to reproduce the notes' algebra exactly.

**Files:**
- Modify: `climt/_components/cork/lw/kernels.py`
- Modify: `climt/_components/cork/lw/component.py`
- Test: `tests/test_cork_diffusivity.py` (create)

**Interfaces:**
- Consumes: nothing.
- Produces: `CorkLongwaveRadiation(optics=..., table=..., diffusivity_factor=1.66)`; the diagnostic `longwave_transmittance_per_band` equals `exp(-diffusivity_factor * longwave_optical_depth_per_band)`.

- [x] **Step 1: Write the failing test**

Create `tests/test_cork_diffusivity.py`:

```python
"""The longwave diffusivity factor is settable per component.

The EC2213 course notes derive the two-stream equations with D = 2 (mu = 0.5);
cork's default is the Elsasser value 1.66. The modelling-tour pages need to
reproduce the notes' algebra exactly, so D must be a constructor argument.
"""
import numpy as np
import pytest
import sympl

import climt
from climt import CorkLongwaveRadiation, get_default_state, get_grid


@pytest.fixture(autouse=True)
def _unyt_backend():
    sympl.set_backend(climt.UnytBackend())


def _run(diffusivity_factor=None):
    kwargs = {"optics": "correlated_k", "table": "earth_low_res_lw"}
    if diffusivity_factor is not None:
        kwargs["diffusivity_factor"] = diffusivity_factor
    lw = CorkLongwaveRadiation(**kwargs)
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=18))
    return lw(state)


def test_default_is_elsasser_166():
    _, diag = _run()
    tau = diag["longwave_optical_depth_per_band"].values
    trans = diag["longwave_transmittance_per_band"].values
    np.testing.assert_allclose(trans, np.exp(-1.66 * tau), rtol=1e-6)


def test_transmittance_follows_the_requested_factor():
    _, diag = _run(diffusivity_factor=2.0)
    tau = diag["longwave_optical_depth_per_band"].values
    trans = diag["longwave_transmittance_per_band"].values
    np.testing.assert_allclose(trans, np.exp(-2.0 * tau), rtol=1e-6)


def test_larger_diffusivity_lowers_olr():
    """A longer effective path means a higher, colder emission level."""
    _, diag_166 = _run(diffusivity_factor=1.66)
    _, diag_2 = _run(diffusivity_factor=2.0)
    olr_166 = float(diag_166["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    olr_2 = float(diag_2["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    assert olr_2 < olr_166
```

- [x] **Step 2: Run the tests to verify they fail**

Run: `conda run -n climt python -m pytest tests/test_cork_diffusivity.py -v`
Expected: `test_default_is_elsasser_166` PASSES (current behaviour); the other two FAIL with `TypeError: __init__() got an unexpected keyword argument 'diffusivity_factor'`.

- [x] **Step 3: Thread the factor through the transport kernel**

In `climt/_components/cork/lw/kernels.py`, add the parameter to the jitted kernel signature:

```python
@njit(parallel=True)
def _lw_transport_kernel(
    tau, planck_source, surface_source, emissivity, weights,
    up_band, down_band, up_broad, down_broad,
    diag_trans, diag_up_gpt, diag_dn_gpt, want_diag, diffusivity_factor,
):
```

Replace **both** occurrences of `np.exp(-DIFFUSIVITY_FACTOR * tau[b, g, k, i])` (the
upward sweep and the downward sweep) with:

```python
                    trans = np.exp(-diffusivity_factor * tau[b, g, k, i])
```

Leave `DIFFUSIVITY_FACTOR = 1.66` defined at module level — it stays the default.

- [x] **Step 4: Thread it through the `lw_transport` wrapper**

Still in `kernels.py`, change the wrapper signature:

```python
def lw_transport(
    T, T_surface, tau, planck_source, surface_source, emissivity, weights, sigma,
    diagnostics_level=0, diffusivity_factor=DIFFUSIVITY_FACTOR,
):
```

Add to its docstring `Args:` block, after `diagnostics_level`:

```
        diffusivity_factor: two-stream diffusivity D in trans = exp(-D * tau).
            Defaults to the Elsasser value 1.66. The EC2213 course notes use 2.
```

and pass it at the call site:

```python
    _lw_transport_kernel(
        tau, planck_source, surface_source, emissivity, weights,
        up_band, down_band, up_broad, down_broad,
        diag_trans, diag_up_gpt, diag_dn_gpt, want_diag, diffusivity_factor,
    )
```

- [x] **Step 5: Add the constructor argument**

In `climt/_components/cork/lw/component.py`, import the default and accept the argument.
Change the import on line 17 to:

```python
from .kernels import DIFFUSIVITY_FACTOR, lw_transport, planck_sources_kernel
```

Add the parameter to `__init__`, after `rosseland_mean_fit`:

```python
        rosseland_mean_fit="freedman2014",
        diffusivity_factor=DIFFUSIVITY_FACTOR,
```

and store it as the first line of the body, before `self._optics_mode` is set:

```python
        self._diffusivity_factor = diffusivity_factor
```

- [x] **Step 6: Use it at both call sites in `array_call`**

In the `lw_transport(...)` call, add the argument after `diagnostics_level`:

```python
            diagnostics_level=self._diagnostics_level,
            diffusivity_factor=self._diffusivity_factor,
```

and replace the hardcoded constant (currently `D = 1.66  # diffusivity factor`):

```python
        D = self._diffusivity_factor
```

- [x] **Step 7: Run the tests to verify they pass**

Run: `conda run -n climt python -m pytest tests/test_cork_diffusivity.py -v`
Expected: 3 passed.

- [x] **Step 8: Verify no regression in the existing cork suite**

Run: `conda run -n climt python -m pytest tests/ -k "cork or grey or gray" -v -m "not slow"`
Expected: all pass. The default is unchanged, so `tests/test_grey_limit.py` and
`tests/test_gray_default_table.py` must be unaffected.

- [x] **Step 9: Regenerate the experiment artifacts this task invalidates**

Editing anything under `climt/_components/cork/` trips the docs workflow's staleness
gate (`.github/workflows/docs.yml`, "Verify experiment artifacts are not stale"), because
the two RCE artifacts in `docs/experiments/2026-06-05-cork-co2-bands/sources.yml`
content-hash the glob `climt/_components/cork/**/*.py`. Skipping this turns the docs job
red on the PR with `STALE: _artifacts/rce_moist_before.npz` and `_after`. Confirm first:

```bash
conda run -n climt python scripts/build_experiments.py --check
```

Regenerate **only** that manifest — a full sweep also re-runs the radiative-transfer
figures for nothing, and the two RCE integrations alone take over an hour:

```bash
make experiments ONLY='docs/experiments/*cork-co2-bands/**'
```

Expected: `rce_moist_before.npz` and `rce_moist_after.npz` come back **bitwise
identical** — `git status` should not list them. Changing `D`'s *plumbing* must not
change results at the default `D = 1.66`; if those files show as modified, the
default leaked and Step 8's suite is not catching it. Stop and fix rather than
committing the diff.

The run also prints `skipping _artifacts/throughput.npz (manual; pass --manual to
run it)` and the same for `throughput.png`. That is correct, not a failure. Those
are wall-clock benchmarks carrying `manual: true`, so they are outside the gate —
their numbers track the machine as much as the code, and re-recording them on a
dev laptop would commit hardware noise as a diff. Leave them alone; do not pass
`--manual` to make the message go away.

Two artifacts in `docs/radiative-transfer/sources.yml` also depend on cork sources —
`06_picket_fence_opacity.png` on `optics/parmentier.py` and `07_two_stream_phases.png`
on `sw/kernels.py`. This task touches neither, so they stay green. Widen the `--only`
glob if a later task edits those files.

- [x] **Step 10: Commit**

```bash
git add climt/_components/cork/lw/kernels.py climt/_components/cork/lw/component.py tests/test_cork_diffusivity.py
git add docs/experiments/2026-06-05-cork-co2-bands/_artifacts/
git commit -m "feat(cork): diffusivity_factor as a constructor argument

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

Then confirm the gate is clean before moving on: `conda run -n climt python
scripts/build_experiments.py --check` must exit 0.

---

## Task 2: Gray table generator — `--diffusivity` flag and path resolution

Two defects block the tour's gray table. The generator resolves its own output **by basename**, so it only works writing into the package data directory. And its calibration hardcodes D = 1.66, which is wrong when the component will run at D = 2.

**Files:**
- Modify: `scripts/generate_gray_default_table.py`
- Test: `tests/test_gray_default_table.py` (extend)

**Interfaces:**
- Consumes: `CorkLongwaveRadiation(diffusivity_factor=...)` from Task 1.
- Produces: `python scripts/generate_gray_default_table.py --output PATH --total-optical-depth TAU --diffusivity D`, writing a table for which `D * sum_layers(k * m_layer) == TAU`.

- [x] **Step 1: Write the failing test**

Append to `tests/test_gray_default_table.py`:

```python
def test_generator_honours_diffusivity_and_absolute_paths(tmp_path):
    """Calibration must use the same D the component will run with.

    The EC2213 notes' tau_inf is the diffusivity-SCALED column optical depth,
    so a table built for D=2 must satisfy 2 * sum(k * m) == tau_inf.
    """
    import subprocess
    import sys

    import numpy as np
    import sympl

    import climt
    from climt import CorkLongwaveRadiation, get_default_state, get_grid

    out = tmp_path / "tour_gray_lw.nc"
    tau_inf, D = 4.0, 2.0
    subprocess.run(
        [sys.executable, "scripts/generate_gray_default_table.py",
         "--output", str(out), "--total-optical-depth", str(tau_inf),
         "--diffusivity", str(D)],
        check=True,
    )

    sympl.set_backend(climt.UnytBackend())
    lw = CorkLongwaveRadiation(optics="correlated_k", table=str(out),
                               diffusivity_factor=D)
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=28))
    _, diag = lw(state)
    tau_col = float(diag["longwave_optical_depth_per_band"].values[:, 0, 0, 0].sum())
    np.testing.assert_allclose(D * tau_col, tau_inf, rtol=1e-3)
```

- [x] **Step 2: Run it to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_gray_default_table.py::test_generator_honours_diffusivity_and_absolute_paths -v`
Expected: FAIL — `generate_gray_default_table.py: error: unrecognized arguments: --diffusivity`.

- [x] **Step 3: Resolve the output by full path**

In `scripts/generate_gray_default_table.py`, `_column_mass_path` currently passes
`os.path.splitext(os.path.basename(nc_path))[0]`. `load_k_table` accepts filesystem
paths, so pass the path itself:

```python
    cork = CorkLongwaveRadiation(
        optics="correlated_k",
        table=nc_path,
    )
```

- [x] **Step 4: Add the `--diffusivity` flag**

Add the argument in `main()`:

```python
    ap.add_argument("--diffusivity", type=float, default=DIFFUSIVITY_FACTOR,
                    help="Diffusivity factor D the component will run with. "
                         "The calibration target is D * sum(k * m) == "
                         "--total-optical-depth, so this MUST match the D "
                         "passed to CorkLongwaveRadiation. The EC2213 notes "
                         "use D=2; cork's default is 1.66.")
```

and use it in the calibration and the report:

```python
    k_value = args.total_optical_depth / (args.diffusivity * sum_m)
    write_single_band_table(args.output, k_value=k_value)
    print(f"wrote {args.output}")
    print(f"  column mass path sum_m = {sum_m:.6f}")
    print(f"  calibrated k = {k_value:.6g}  "
          f"(raw total tau = {k_value*sum_m:.6f}, "
          f"D-scaled at D={args.diffusivity} = "
          f"{args.diffusivity*k_value*sum_m:.6f})")
```

- [x] **Step 5: Run the tests to verify they pass**

Run: `conda run -n climt python -m pytest tests/test_gray_default_table.py -v`
Expected: all pass, including the pre-existing tests (the default `--diffusivity` is
1.66, so the shipped table's calibration is unchanged).

- [x] **Step 6: Commit**

```bash
git add scripts/generate_gray_default_table.py tests/test_gray_default_table.py
git commit -m "fix(scripts): gray table generator honours --diffusivity and full paths

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 3: Generate and ship the tour's gray table

Page 4 needs a gray table whose D-scaled column optical depth is a known `τ_∞`, at the notes' D = 2.

**Files:**
- Create: `climt/_data/cork/correlated_k/tour_gray_lw.nc`
- Create: `climt/_data/cork/correlated_k/tour_gray_lw.npz`
- Modify: `climt/_data/cork/correlated_k/MANIFEST.md`

**Interfaces:**
- Produces: table name `"tour_gray_lw"`, resolvable by `load_k_table("tour_gray_lw")`, for which `2.0 * sum_layers(k * m_layer) == 4.0`.

- [x] **Step 1: Generate the `.nc` and convert to `.npz`**

```bash
conda run -n climt python scripts/generate_gray_default_table.py \
    --output climt/_data/cork/correlated_k/tour_gray_lw.nc \
    --total-optical-depth 4.0 --diffusivity 2.0
conda run -n climt python scripts/convert_ck_table_to_npz.py \
    climt/_data/cork/correlated_k/tour_gray_lw.nc
```

Expected: the script reports `D-scaled at D=2.0 = 4.000000`.

- [x] **Step 2: Record the checksum**

```bash
shasum -a 256 climt/_data/cork/correlated_k/tour_gray_lw.nc
```

- [x] **Step 3: Document it in the manifest**

Add a row to the "File listing" table in `climt/_data/cork/correlated_k/MANIFEST.md`
(paste the real sha256 from Step 2), and a paragraph after the `single_band_gray_lw`
paragraph:

```markdown
`tour_gray_lw.nc` is a synthetic single-band, constant-k fixture for the
Modelling Tour's gray-equilibrium page. Unlike `single_band_gray_lw` (calibrated
at cork's default D = 1.66), it is calibrated at the **EC2213 course notes'
D = 2** to a D-scaled column optical depth of `tau_inf = 4`. The notes' `tau` is
the diffusivity-scaled depth, so the page must construct the component with
`diffusivity_factor=2.0` for the analytic solution to hold. Regenerate with

```sh
python scripts/generate_gray_default_table.py \
    --output climt/_data/cork/correlated_k/tour_gray_lw.nc \
    --total-optical-depth 4.0 --diffusivity 2.0
python scripts/convert_ck_table_to_npz.py \
    climt/_data/cork/correlated_k/tour_gray_lw.nc
```
```

- [x] **Step 4: Verify it loads by name and has the right depth**

```bash
conda run -n climt python -c "
import sympl, climt
from climt import CorkLongwaveRadiation, get_default_state, get_grid
sympl.set_backend(climt.UnytBackend())
lw = CorkLongwaveRadiation(optics='correlated_k', table='tour_gray_lw', diffusivity_factor=2.0)
st = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=28))
_, d = lw(st)
print('D-scaled column tau =', 2.0 * d['longwave_optical_depth_per_band'].values[:,0,0,0].sum())"
```

Expected: `D-scaled column tau = 4.000...`

- [x] **Step 5: Commit**

```bash
git add climt/_data/cork/correlated_k/tour_gray_lw.nc climt/_data/cork/correlated_k/tour_gray_lw.npz climt/_data/cork/correlated_k/MANIFEST.md
git commit -m "feat(data): tour_gray_lw table calibrated at the notes' D=2

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 4: `_tour/soundings.py` — prescribed profiles

**Files:**
- Create: `docs/modelling-tour/_tour/soundings.py`
- Create: `docs/modelling-tour/_tour/__init__.py` (empty)
- Test: `tests/test_modelling_tour.py` (create)

**Interfaces:**
- Produces:
  - `saturation_vapour_pressure(T) -> ndarray` (Pa)
  - `lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8, gamma=6.5e-3, T_strat=200.0) -> (T, q)`
  - `analytic_gray_equilibrium(p, ps, tau_inf=4.0, T_e=255.0) -> (T, T_surf, tau_star)`
  - `apply_sounding(state, T, q=None, T_surf=None) -> None` (mutates in place)

- [x] **Step 1: Write the failing test**

Create `tests/test_modelling_tour.py`:

```python
"""The science the Modelling Tour pages claim, checked natively.

Pyodide cells cannot run in CI, so each page's computational core lives in
``docs/modelling-tour/_tour/`` as importable Python and is exercised here. This
mirrors the pattern of ``tests/test_live_rce_demo.py``.

``docs/modelling-tour`` contains a hyphen and is not a valid package path, so
the helpers are loaded by file path.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest
import sympl

import climt
from climt import CorkLongwaveRadiation, get_default_state, get_grid

REPO_ROOT = Path(__file__).resolve().parent.parent
TOUR = REPO_ROOT / "docs/modelling-tour/_tour"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, TOUR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(autouse=True)
def _unyt_backend():
    sympl.set_backend(climt.UnytBackend())


@pytest.fixture
def soundings():
    return _load("soundings")


def test_lapse_rate_sounding_shape_and_bounds(soundings):
    p = np.linspace(1e5, 1e3, 40)
    T, q = soundings.lapse_rate_sounding(p, 1e5, T_surf=288.0)
    assert T.shape == p.shape and q.shape == p.shape
    assert T[0] == pytest.approx(288.0, abs=0.5)     # bottom ~ surface
    assert T.min() == pytest.approx(200.0, abs=1e-6)  # stratosphere floor
    assert np.all(np.diff(T) <= 1e-9)                 # monotone with height
    assert np.all(q > 0) and q.max() < 0.05           # sane humidity


def test_saturation_vapour_pressure_matches_known_values(soundings):
    # Bolton (1980): ~611.2 Pa at 0 C, ~2339 Pa at 20 C
    assert soundings.saturation_vapour_pressure(273.15) == pytest.approx(611.2, rel=1e-3)
    assert soundings.saturation_vapour_pressure(293.15) == pytest.approx(2339.0, rel=2e-2)


def test_analytic_gray_equilibrium_reproduces_the_notes(soundings):
    """T_skin = 2^-0.25 T_e and T_g = T_e (1 + tau_inf/2)^0.25."""
    p = np.linspace(1e5, 1.0, 60)
    T, T_surf, tau_star = soundings.analytic_gray_equilibrium(p, 1e5, tau_inf=4.0,
                                                              T_e=255.0)
    assert T[-1] == pytest.approx(255.0 * 0.5 ** 0.25, rel=1e-3)
    assert T_surf == pytest.approx(255.0 * (1 + 4.0 / 2) ** 0.25, rel=1e-9)
    assert tau_star[0] == pytest.approx(4.0, rel=1e-9)   # deepest at the surface
    assert tau_star[-1] == pytest.approx(0.0, abs=1e-3)  # zero at TOA
```

- [x] **Step 2: Run it to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -v`
Expected: FAIL — `FileNotFoundError` / `spec_from_file_location` returns None, because
`docs/modelling-tour/_tour/soundings.py` does not exist.

- [x] **Step 3: Write the implementation**

Create `docs/modelling-tour/_tour/__init__.py` as an empty file, and
`docs/modelling-tour/_tour/soundings.py`:

```python
"""Prescribed vertical profiles for the Modelling Tour pages.

Nothing here integrates in time. Every tour page states a profile outright and
asks the radiation code what it makes of it, which keeps each cell to a single
component call and keeps the causal structure visible.

Kept importable and natively testable: the pages import these functions rather
than defining profiles inline, and ``tests/test_modelling_tour.py`` guards them.
"""
import numpy as np

RD = 287.0        # gas constant for dry air, J/(kg K)
G = 9.81          # gravitational acceleration, m/s^2
EPSILON = 0.622   # ratio of molar masses, water vapour to dry air


def saturation_vapour_pressure(T):
    """Saturation vapour pressure over liquid water, in Pa (Bolton 1980)."""
    Tc = np.asarray(T) - 273.15
    return 611.2 * np.exp(17.67 * Tc / (Tc + 243.5))


def lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8, gamma=6.5e-3,
                        T_strat=200.0, q_floor=1e-7):
    """A troposphere at a constant lapse rate under an isothermal stratosphere.

    Hydrostatic balance with a constant lapse rate gives T(p) = T_surf *
    (p/ps)**(RD*gamma/G); the profile is then clipped at ``T_strat``. Humidity
    is set at fixed relative humidity ``rh`` and floored so the stratosphere
    stays inside the k-table's H2O axis.

    Args:
        p: (nz,) mid-level pressure, Pa, surface first.
        ps: surface pressure, Pa.
        T_surf: surface temperature, K.
        rh: relative humidity, dimensionless (0-1).
        gamma: lapse rate, K/m.
        T_strat: isothermal stratosphere temperature, K.
        q_floor: minimum specific humidity, kg/kg.

    Returns:
        (T, q) each (nz,) — air temperature in K, specific humidity in kg/kg.
    """
    p = np.asarray(p, dtype=float)
    T = T_surf * (p / ps) ** (RD * gamma / G)
    T = np.maximum(T, T_strat)
    e = rh * saturation_vapour_pressure(T)
    q = EPSILON * e / np.maximum(p - e, 1.0)
    return T, np.maximum(q, q_floor)


def analytic_gray_equilibrium(p, ps, tau_inf=4.0, T_e=255.0):
    """The course notes' closed-form gray radiative-equilibrium profile.

    Chapter 8 of *Principles of Planetary Climate* derives

        T(tau)  = T_e * [(1 + tau_inf - tau) / 2]**0.25
        T_ground = T_e * (1 + tau_inf / 2)**0.25

    with ``tau`` measured UP from the surface, so ``tau_star = tau_inf - tau``
    is measured down from space. **``tau_inf`` is the diffusivity-SCALED column
    optical depth** — a component evaluating this profile must be built with the
    same ``diffusivity_factor`` the table was calibrated against, or the heating
    rate will not vanish.

    Optical depth is taken linear in pressure (a well-mixed absorber):
    ``tau_star = tau_inf * p / ps``.

    Args:
        p: (nz,) mid-level pressure, Pa, surface first.
        ps: surface pressure, Pa.
        tau_inf: diffusivity-scaled column optical depth.
        T_e: emission temperature, K.

    Returns:
        (T, T_surf, tau_star) — (nz,) air temperature in K, scalar surface
        temperature in K, and (nz,) optical depth measured down from space.
    """
    p = np.asarray(p, dtype=float)
    tau_star = tau_inf * p / ps
    T = T_e * ((1.0 + tau_star) / 2.0) ** 0.25
    T_surf = T_e * (1.0 + tau_inf / 2.0) ** 0.25
    return T, T_surf, tau_star


def apply_sounding(state, T, q=None, T_surf=None):
    """Write a prescribed profile into a climt state, in place.

    Args:
        state: a climt state dict from ``get_default_state``.
        T: (nz,) air temperature, K.
        q: (nz,) specific humidity in kg/kg, or None to leave it alone.
        T_surf: scalar surface temperature in K, or None to leave it alone.
    """
    state["air_temperature"].values[:, 0, 0] = T
    if T_surf is not None:
        state["surface_temperature"].values[:] = T_surf
    if q is not None and "specific_humidity" in state:
        state["specific_humidity"].values[:, 0, 0] = q
```

- [x] **Step 4: Run the tests to verify they pass**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -v`
Expected: 3 passed.

- [x] **Step 5: Commit**

```bash
git add docs/modelling-tour/_tour/__init__.py docs/modelling-tour/_tour/soundings.py tests/test_modelling_tour.py
git commit -m "feat(tour): prescribed sounding helpers

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 5: `_tour/spectra.py` — Planck, brightness temperature, weighting function

**Files:**
- Create: `docs/modelling-tour/_tour/spectra.py`
- Test: `tests/test_modelling_tour.py` (extend)

**Interfaces:**
- Produces:
  - `planck_flux(nu, T) -> ndarray` — π·B in W m⁻² / cm⁻¹, `nu` in cm⁻¹
  - `brightness_temperature(flux_density, nu) -> ndarray` — K
  - `band_centres(band_limits) -> ndarray`, `band_widths(band_limits) -> ndarray`
  - `band_limits_of(component) -> ndarray` of shape `(nband, 2)`
  - `spectral_olr(diagnostics, band_limits) -> (olr_per_band, flux_density)`
  - `emission_weight(tau_layer, diffusivity_factor) -> ndarray` — per-layer contribution
  - `tau_star_cumulative(tau_layer, diffusivity_factor) -> ndarray`

- [x] **Step 1: Write the failing test**

Append to `tests/test_modelling_tour.py`:

```python
@pytest.fixture
def spectra():
    return _load("spectra")


def test_brightness_temperature_inverts_planck(spectra):
    nu = np.array([300.0, 667.0, 900.0, 1400.0])
    for T in (200.0, 255.0, 288.0, 320.0):
        flux = spectra.planck_flux(nu, T)
        np.testing.assert_allclose(
            spectra.brightness_temperature(flux, nu), T, rtol=1e-8)


def test_per_band_olr_sums_to_the_broadband_diagnostic(spectra):
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=30))
    _, diag = lw(state)
    limits = spectra.band_limits_of(lw)
    olr_band, _ = spectra.spectral_olr(diag, limits)
    broadband = float(diag["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    assert olr_band.sum() == pytest.approx(broadband, rel=1e-9)


def test_band_limits_cover_the_expected_features(spectra):
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    limits = spectra.band_limits_of(lw)
    assert limits.shape == (14, 2)
    # the 15 um CO2 core and the window edges must be exact band boundaries
    edges = set(np.round(limits.ravel(), 1))
    for edge in (630.0, 700.0, 800.0, 980.0, 1080.0, 1180.0):
        assert edge in edges


def test_emission_weight_peaks_where_tau_star_is_one(spectra):
    """The notes' weighting function: emission * transmission peaks at tau*=1."""
    tau_layer = np.full(200, 0.05)      # uniform layers, column tau = 10
    w = spectra.emission_weight(tau_layer, diffusivity_factor=1.0)
    tau_star = spectra.tau_star_cumulative(tau_layer, diffusivity_factor=1.0)
    assert tau_star[-1] == pytest.approx(0.0, abs=0.06)   # ~0 at TOA
    assert tau_star[0] == pytest.approx(10.0, rel=0.02)   # deepest at surface
    np.testing.assert_allclose(tau_star[np.argmax(w)], 1.0, atol=0.1)
```

- [x] **Step 2: Run it to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -v -k spectra or "planck or olr or band or emission"`
Expected: FAIL — `spectra.py` does not exist.

- [x] **Step 3: Write the implementation**

Create `docs/modelling-tour/_tour/spectra.py`:

```python
"""Spectral diagnostics built on cork's per-band output.

The EC2213 notes stop at "the OLR is this integrated over all lambda and over
all theta". Cork performs that integral; these helpers un-do it, turning
``*_per_band`` diagnostics into the spectra and weighting functions the notes
describe analytically.

Wavenumbers are in cm^-1 throughout, so ``c`` is carried in cm/s and spectral
flux densities come out in W m^-2 / cm^-1 — directly comparable to published
IRIS/Nimbus spectra.
"""
import numpy as np

H_PLANCK = 6.62607015e-34    # J s
C_LIGHT = 2.99792458e10      # cm/s  (note: cm, to pair with nu in cm^-1)
K_BOLTZ = 1.380649e-23       # J/K

# 2 h c^2 nu^3 comes out in W cm^-2 sr^-1 / cm^-1; 1e4 converts to per m^-2,
# and pi integrates isotropic intensity into a flux.
_PLANCK_PREFACTOR = np.pi * 1e4 * 2.0 * H_PLANCK * C_LIGHT ** 2


def planck_flux(nu, T):
    """Blackbody spectral flux density pi*B, in W m^-2 / cm^-1.

    Args:
        nu: wavenumber(s) in cm^-1.
        T: temperature(s) in K.
    """
    nu = np.asarray(nu, dtype=float)
    return _PLANCK_PREFACTOR * nu ** 3 / np.expm1(
        H_PLANCK * C_LIGHT * nu / (K_BOLTZ * np.asarray(T, dtype=float)))


def brightness_temperature(flux_density, nu):
    """Temperature of a blackbody radiating ``flux_density`` at ``nu``.

    The exact inverse of :func:`planck_flux`.

    Args:
        flux_density: spectral flux density in W m^-2 / cm^-1.
        nu: wavenumber(s) in cm^-1.
    """
    nu = np.asarray(nu, dtype=float)
    flux_density = np.asarray(flux_density, dtype=float)
    return (H_PLANCK * C_LIGHT * nu / K_BOLTZ) / np.log1p(
        _PLANCK_PREFACTOR * nu ** 3 / flux_density)


def band_limits_of(component):
    """(nband, 2) band edges in cm^-1 for a correlated-k cork component."""
    return np.asarray(component._table["band_wavenumber_limits"], dtype=float)


def band_centres(band_limits):
    """(nband,) band-centre wavenumbers in cm^-1."""
    return 0.5 * (band_limits[:, 0] + band_limits[:, 1])


def band_widths(band_limits):
    """(nband,) band widths in cm^-1."""
    return band_limits[:, 1] - band_limits[:, 0]


def spectral_olr(diagnostics, band_limits, column=0):
    """Per-band OLR and spectral flux density at the top of the atmosphere.

    Args:
        diagnostics: the dict returned by a ``CorkLongwaveRadiation`` call.
        band_limits: (nband, 2) band edges, from :func:`band_limits_of`.
        column: horizontal index to extract.

    Returns:
        (olr_per_band, flux_density) — W m^-2 per band, and W m^-2 / cm^-1.
    """
    up = diagnostics["upwelling_longwave_flux_in_air_per_band"].values
    olr_band = np.asarray(up[-1, column, 0, :] if up.ndim == 4
                          else up[-1, column, :], dtype=float)
    return olr_band, olr_band / band_widths(band_limits)


def tau_star_cumulative(tau_layer, diffusivity_factor):
    """Diffusivity-scaled optical depth measured DOWN from space.

    Args:
        tau_layer: (nz,) per-layer optical depth, surface first.
        diffusivity_factor: D, so the result is D * sum of layers above.

    Returns:
        (nz,) tau* at each mid-level — largest at the surface, ~0 at the top.
    """
    tau_layer = np.asarray(tau_layer, dtype=float)
    # Cumulative sum of everything strictly above each level.
    above = np.cumsum(tau_layer[::-1])[::-1] - tau_layer
    return diffusivity_factor * (above + 0.5 * tau_layer)


def emission_weight(tau_layer, diffusivity_factor):
    """How much of the OLR each layer actually emits.

    The notes derive W ∝ tau* exp(-tau*), peaking at tau* = 1, by assuming
    thin layers and an exponential density profile. This is the exact discrete
    equivalent for the model's own grid: a layer emits ``1 - exp(-D tau_k)`` of
    a blackbody and that emission is attenuated by ``exp(-tau*_above)`` on the
    way out. It reduces to the notes' form as the layers get thin.

    Args:
        tau_layer: (nz,) per-layer optical depth, surface first.
        diffusivity_factor: D.

    Returns:
        (nz,) dimensionless weight per layer, surface first.
    """
    tau_layer = np.asarray(tau_layer, dtype=float)
    above = np.cumsum(tau_layer[::-1])[::-1] - tau_layer
    emitted = 1.0 - np.exp(-diffusivity_factor * tau_layer)
    transmitted = np.exp(-diffusivity_factor * above)
    return emitted * transmitted
```

- [x] **Step 4: Run the tests to verify they pass**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -v`
Expected: 7 passed. **As built: 10 passed** — see the log below; the weighting-function
test was rewritten and parametrised over D.

- [x] **Step 5: Commit**

```bash
git add docs/modelling-tour/_tour/spectra.py tests/test_modelling_tour.py
git commit -m "feat(tour): spectral diagnostic helpers

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

### Log — the weighting-function test as written in Step 1 was wrong

`test_emission_weight_peaks_where_tau_star_is_one` fed `emission_weight` **200 layers of
equal optical depth** (`np.full(200, 0.05)`) and expected the peak at `tau* = 1`. It came
out at `tau* = 0.025` — the topmost layer — and that is the correct answer for that input,
not a bug in the helper. With every layer emitting the same `1 - exp(-D tau_k)`, the only
term that varies is the attenuation `exp(-D tau*_above)`, which is largest at the top, so
the weight is monotone by construction and has no interior peak.

The notes' `W ∝ tau* exp(-tau*)` assumes layers of equal **thickness** over an
exponentially stratified absorber, where `dtau/dz ∝ exp(-z/H)` makes the per-layer optical
depth grow towards the surface — and that is what the pages will actually feed it, since
cork's per-layer `tau` on a real column behaves that way.

Fixed by replacing the fixture, not the physics: a `_exponential_absorber()` helper in the
test builds `tau_layer` from `tau_inf * exp(-z/H)` on a uniform-z edge grid. The peak then
lands at `tau* = 0.999` (D=1), `0.984` (D=1.66), `0.980` (D=2), so the test is now
parametrised over all three. `tau_star[0]` is `9.80`, not `10.0`, for `tau_inf = 10` —
the half-layer offset at the optically thick bottom layer — hence `rel=0.03`.

The equal-tau case is kept as its own test
(`test_emission_weight_is_monotone_for_equal_optical_depth_layers`) so the trap stays
documented: **a radiating level needs a grid that resolves height, not optical depth.**
Task 9's page cannot demonstrate the peak on an equal-tau grid.

---

## Task 6: Section scaffolding and the landing page

**Files:**
- Create: `docs/modelling-tour/index.qmd`
- Modify: `docs/_quarto.yml`

**Interfaces:**
- Produces: the `Modelling Tour` navbar entry and sidebar; `docs/modelling-tour/**/*.qmd` in the render list.

- [x] **Step 1: Add the section to the Quarto config**

In `docs/_quarto.yml`, add to `project: render:` after the `radiative-transfer` line:

```yaml
    - "modelling-tour/**/*.qmd"
```

Add to `website: navbar: left:` after the "Radiative Transfer" entry:

```yaml
      - text: "Modelling Tour"
        href: modelling-tour/index.qmd
```

Add a sidebar block after the `Radiative Transfer` sidebar block:

```yaml
    - title: "Modelling Tour"
      style: "docked"
      contents:
        - modelling-tour/index.qmd
        - modelling-tour/01-emissivity-spectrum.qmd
        - modelling-tour/02-window-measured.qmd
        - modelling-tour/03-radiating-level.qmd
        - modelling-tour/04-gray-equilibrium-tested.qmd
        - modelling-tour/05-co2-knob.qmd
        - modelling-tour/06-water-vapour-limit.qmd
```

- [x] **Step 2: Write the landing page**

Create `docs/modelling-tour/index.qmd`. It must contain, in order:

1. Title: `A Modelling Tour of the Climate System`.
2. An opening paragraph stating what this is: the computational half of
   *Principles of Planetary Climate* (EC2213), whose notes derive analytically and
   carry no figures. Link to <https://joymonteiro.github.io/principles_planetary_climate/>.
3. A `.callout-note` titled "Runs in your browser" saying every page runs climt itself
   via Pyodide — no install — and that the first cell on each page downloads a wheel and
   takes a few seconds.
4. A paragraph stating the two threads carried by every page: the science, and *using
   climt*. Say plainly that the craft thread is cumulative — page N assumes page N−1.
5. This chapter map as a markdown table:

| Page | Course chapter | What you compute |
|---|---|---|
| [ε is not a number](01-emissivity-spectrum.qmd) | 5 — The Shell Model | the OLR spectrum your single ε averages over |
| [The window, measured](02-window-measured.qmd) | 6, 7 | per-band optical depth; the window's real width |
| [Where photons come from](03-radiating-level.qmd) | 7 — Basics of Radiative Transfer | the weighting function, and a radiating level per band |
| [Gray equilibrium, tested](04-gray-equilibrium-tested.qmd) | 8 — Gray Gas radiative equilibrium | your analytic profile, checked by a radiation code |
| [The CO₂ knob](05-co2-knob.qmd) | 4 — A deeper look at the zero-dimensional model | the forcing per doubling, and which bands supply it |
| [Water vapour, and the limit](06-water-vapour-limit.qmd) | 4 and beyond | the feedback, drawn; and the runaway limit |

6. A short "What this tranche does not cover" section naming time integration, convection
   and surface fluxes, and saying they belong to a later tranche built on chapters 10–12.
7. A "Going deeper" line linking to `../radiative-transfer/index.qmd` for the
   research-grade treatment.

- [x] **Step 3: Render and check**

```bash
cd docs && QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render modelling-tour/index.qmd
```

Expected: renders with no error; the sidebar shows the six page entries (they 404 until
Tasks 7–12 land, which is expected at this point).

- [x] **Step 4: Commit**

```bash
git add docs/_quarto.yml docs/modelling-tour/index.qmd
git commit -m "docs(tour): section scaffolding and landing page

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 7: Page 1 — ε is not a number

**Files:**
- Create: `docs/modelling-tour/01-emissivity-spectrum.qmd`
- Test: `tests/test_modelling_tour.py` (extend)

**Interfaces:**
- Consumes: `soundings.lapse_rate_sounding`, `soundings.apply_sounding`,
  `spectra.planck_flux`, `spectra.brightness_temperature`, `spectra.band_limits_of`,
  `spectra.spectral_olr`.

- [x] **Step 1: Write the failing test**

Append to `tests/test_modelling_tour.py`:

```python
@pytest.mark.slow
def test_page1_window_is_warmer_than_the_co2_core(soundings, spectra):
    """The page's central claim: brightness temperature swings ~200 K to ~285 K."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, q = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)
    soundings.apply_sounding(state, T, q, T_surf=288.0)

    _, diag = lw(state)
    limits = spectra.band_limits_of(lw)
    _, flux_density = spectra.spectral_olr(diag, limits)
    tb = spectra.brightness_temperature(flux_density, spectra.band_centres(limits))

    centres = spectra.band_centres(limits)
    co2_core = tb[(centres > 630) & (centres < 700)]
    window = tb[(centres > 800) & (centres < 1180)]

    assert co2_core.max() < 230.0          # radiating from the stratosphere
    assert window.min() > 270.0            # radiating from near the surface
    assert window.min() < 288.0            # but not from the surface itself
    assert window.mean() - co2_core.mean() > 60.0
```

- [x] **Step 2: Run it to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -v -m slow`
Expected: FAIL — the fixtures exist, so this should actually PASS on the physics. If it
fails, the sounding or the spectral helpers are wrong; fix them before writing the page.
(This test guards the page's claim; it is written first so the page cannot assert
something untrue.)

- [x] **Step 3: Write the page**

Create `docs/modelling-tour/01-emissivity-spectrum.qmd`.

Front matter — identical on all six pages:

```yaml
---
title: "ε is not a number"
format: live-html
engine: jupyter
pyodide:
  packages:
    - "climt==0.30.0"
---
```

Then, in order:

1. `{{< include ../_includes/climt-live-setup.qmd >}}`
2. A `.callout-tip` titled "Where this fits" saying: chapter 5 of the course notes builds
   a shell of emissivity ε and gets a surface temperature out of it; this page asks what
   that single number is standing in for.
3. Prose recapping the shell result: one ε, one radiating level, one emission temperature.
   State the question the page answers: *is the atmosphere's emissivity one number?*
4. The setup cell:

```{pyodide}
#| autorun: true
import numpy as np
import matplotlib.pyplot as plt
import sympl, climt
from climt import CorkLongwaveRadiation, get_default_state, get_grid
sympl.set_backend(climt.UnytBackend())

T_SURF, RH = 288.0, 0.8          # <- change these and re-run

lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
p = state["air_pressure"].values[:, 0, 0]
ps = float(state["surface_air_pressure"].values.ravel()[0])

T, q = lapse_rate_sounding(p, ps, T_surf=T_SURF, rh=RH)
apply_sounding(state, T, q, T_surf=T_SURF)
print(f"prescribed a {T_SURF:.0f} K surface, RH = {RH:.0%}, {len(p)} levels")
```

5. Prose introducing the single call and what comes back — name
   `upwelling_longwave_flux_in_air_per_band` explicitly and say it has a band axis.
6. The figure cell:

```{pyodide}
tendencies, diagnostics = lw(state)
limits = band_limits_of(lw)
nu, width = band_centres(limits), band_widths(limits)
olr_band, flux_density = spectral_olr(diagnostics, limits)
tb = brightness_temperature(flux_density, nu)

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
nu_fine = np.linspace(20, 1800, 400)
for T_ref in (300, 280, 260, 240, 220, 200):
    ax.plot(nu_fine, planck_flux(nu_fine, T_ref), color="0.8", lw=0.8)
    ax.annotate(f"{T_ref} K", (1800, planck_flux(1800.0, T_ref)),
                fontsize=7, color="0.5", va="center")
ax.step(limits[:, 1], flux_density, where="pre", color="tab:blue", lw=2)
ax.set(xlabel="wavenumber (cm$^{-1}$)",
       ylabel="spectral OLR (W m$^{-2}$ / cm$^{-1}$)", xlim=(0, 1800))
ax.set_title(f"OLR spectrum — total {olr_band.sum():.1f} W m$^{{-2}}$")

ax2.step(limits[:, 1], tb, where="pre", color="tab:blue", lw=2)
ax2.axhline(T_SURF, color="k", ls=":", lw=0.8)
ax2.annotate(f"surface, {T_SURF:.0f} K", (1450, T_SURF + 2), fontsize=8)
for lo, hi, label in [(630, 700, "CO$_2$"), (800, 1180, "window")]:
    ax2.axvspan(lo, hi, color="0.92", zorder=0)
    ax2.annotate(label, (0.5 * (lo + hi), 196), fontsize=8, ha="center")
ax2.set(xlabel="wavenumber (cm$^{-1}$)", ylabel="brightness temperature (K)",
        xlim=(0, 1800), ylim=(190, 300))
ax2.set_title("Where each band radiates from")
fig.tight_layout()

sigma_T4 = 5.670374419e-8 * T_SURF ** 4
print(f"surface emission  sigma T^4 = {sigma_T4:.1f} W/m2")
print(f"OLR                         = {olr_band.sum():.1f} W/m2")
print(f"greenhouse effect           = {sigma_T4 - olr_band.sum():.1f} W/m2")
print(f"brightness temperature ranges {tb.min():.0f} K to {tb.max():.0f} K")
```

7. Prose reading the figure. State the numbers the test pins: the CO₂ core radiates near
   200 K — from the stratosphere — while the window radiates near 285 K, just below the
   surface. Then the payoff sentence: a single ε is an average over a curve that swings
   ~85 K, and the shell model's radiating level is the average of levels that differ by
   more than 10 km.
8. A `.callout-note` titled "climt craft" covering: components are objects you construct;
   `get_grid` and `get_default_state` build a state that satisfies the component's
   requirements; calling a component returns `(tendencies, diagnostics)`; and
   `lw.diagnostic_properties` lists what you can ask for. Include:

```{pyodide}
for name, spec in sorted(lw.diagnostic_properties.items()):
    print(f"{name:52s} {spec['units']:12s} {spec['dims']}")
```

9. `## Exercises`, with `### Physics` and `### Code` subsections:
   - *Physics 1:* raise `T_SURF` by 10 K. Does the window brightness temperature rise by
     10 K? Does the CO₂ core? Explain the difference.
   - *Physics 2:* the greenhouse effect printed above is σT⁴ − OLR. Which bands
     contribute most to it? Compute `planck_flux(nu, T_SURF) * width - olr_band`.
   - *Physics 3:* set `RH = 0.1`. Which part of the spectrum changes most, and why?
   - *Code 1:* print the band limits alongside each band's OLR, sorted by contribution.
   - *Code 2:* using `lw.input_properties`, list every quantity this component needs and
     its units.
10. `## Going deeper` linking to `../radiative-transfer/01-why-nongrey.qmd`.

- [x] **Step 4: Render the page**

```bash
cd docs && QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render modelling-tour/01-emissivity-spectrum.qmd
```

Expected: renders with no error.

- [x] **Step 5: Run the guard test**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py::test_page1_window_is_warmer_than_the_co2_core -v -m slow`
Expected: PASS.

- [x] **Step 6: Commit**

```bash
git add docs/modelling-tour/01-emissivity-spectrum.qmd tests/test_modelling_tour.py
git commit -m "docs(tour): page 1, the spectrum behind a single emissivity

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

### Log — how the `_tour` helpers reach the browser (applies to Tasks 8–12)

The plan's page cells call `lapse_rate_sounding(...)`, `band_limits_of(...)` and friends
by bare name, but `docs/_includes/climt-live-setup.qmd` defines none of them and
`docs/modelling-tour/_tour/` is not importable under Pyodide. Rather than duplicate the
helper bodies into an include (the `integrate_to_equilibrium` precedent), page 1 hands
the modules to quarto-live's **`pyodide: resources:`**, which fetches each listed path
relative to the page URL at document setup and writes it into the Pyodide filesystem
under the cwd. So the browser imports the *same* files `tests/test_modelling_tour.py`
imports natively — no second copy to keep in sync.

Two pieces make that work, and **pages 2–6 need both**:

1. Front matter gains, alongside `packages:`:

```yaml
  resources:
    - _tour/soundings.py
    - _tour/spectra.py
```

2. `docs/_quarto.yml` gained a `project: resources:` entry for
   `modelling-tour/_tour/*.py`. Quarto omits underscore-prefixed paths from `_site` by
   default, so without it the fetch 404s. Already landed with this task — do not add it
   again; extend the glob if a later task adds a helper module (Task 13's `tables.py`
   is covered by the existing `*.py` glob).

The page's setup cell then does `sys.path.insert(0, "_tour")` before importing.

Verified: `_site/modelling-tour/_tour/{soundings,spectra}.py` are served (HTTP 200), and
the rendered page's base64 `<script type="vfs-file">` manifest lists exactly
`["_tour/soundings.py", "_tour/spectra.py"]`. **Not** verified end-to-end in a browser —
a local run of the rendered page never finished booting Pyodide in the available session
(the tab stopped responding to the automation extension while downloading), so the first
real in-browser check of this mechanism is still outstanding. Native execution of the
page's cells is guarded by `test_page1_window_is_warmer_than_the_co2_core`.

### Log — `plt.show()` is required in the browser (applies to Tasks 8–12)

A figure cell that ends without `plt.show()` prints its text output and draws **nothing**
in quarto-live: the runtime paints the figure a cell explicitly shows, not whatever
matplotlib has open when the cell finishes. Native checking cannot catch this — a script
that ends in `fig.savefig(...)` is perfectly happy — so every page cell that produces a
figure must end on `plt.show()`, and the only way to confirm it is to load the rendered
page. (This is why the RCE include's `_draw_evolution` ends that way.) Caught on page 1
by a human loading the page, after the cell shipped printing its numbers and no figure.

### Log — figure and prose deviations from the step's literal code

- `ax.step(limits[:, 1], ..., where="pre")` as written drops the first band entirely
  (the line starts at the first upper edge, so band 10–250 cm⁻¹ never gets a segment).
  Replaced with `ax.stairs(values, edges)` on `edges = np.append(limits[0, 0],
  limits[:, 1])`, which draws all 14.
- The Planck reference curves are labelled at their own peaks, not at `nu = 1800`, where
  all six labels collapse into each other inside the bottom 10% of the axis.
- The last band (1800–3250 cm⁻¹) returns a brightness temperature of **300 K**, above the
  288 K surface — the band is 1450 cm⁻¹ wide and Planck is nowhere near linear across it,
  so the band-mean flux density is not the flux density at the band centre (πB at 288 K
  is ~240× larger at 1800 than at 3250 cm⁻¹, so the band mean is set by its low-wavenumber
  edge while the inversion is asked about the centre, 2525 cm⁻¹). The printed range is
  restricted to `nu < 1800`, and the same for the figure's x-axis.

  This is **promoted to a teaching beat, not hidden**: the cell prints the last band's
  T_b explicitly, a short section poses the "nothing here is 300 K — what could do that?"
  question, and the explanation sits in a collapsed callout to open after the class has
  had a go. Both panels' `xlim` therefore run the full width of the table, to 3250 cm⁻¹,
  rather than stopping at 1800 — an artifact you tell the class to look at must be *on*
  the figure. Sharing the range across the two panels is what makes the point land: the
  OLR panel is flat on the floor exactly where the temperature panel goes hot. Cost is
  ~45% horizontal compression of the 0–1800 detail, which the CO₂ core and window
  survive; the alternative (different ranges per panel) invites reading the same
  wavenumber at two different x positions and was rejected. It closes by noting the band carries 4.8 of the 246 W m⁻², i.e. the diagnostic
  is nonsense exactly where the energy is negligible. Physics exercise 2 reuses it (the
  same band is the one with a negative per-band greenhouse contribution, −1.9 W m⁻²).

Numbers as built, at `T_SURF = 288 K`, `RH = 0.8`, nz = 40: OLR 245.9 W m⁻², σT⁴ 390.1,
greenhouse effect 144.2; brightness temperature 200.3 K in the CO₂ core (630–700) and
282.7 / 284.7 / 283.5 K across the three window bands.


---

## Task 8: Page 2 — The window, measured

**Files:**
- Create: `docs/modelling-tour/02-window-measured.qmd`
- Test: `tests/test_modelling_tour.py` (extend)

**Interfaces:**
- Consumes: Task 4 and Task 5 helpers.

- [x] **Step 1: Write the failing test**

Append to `tests/test_modelling_tour.py`:

```python
@pytest.mark.slow
def test_page2_window_is_transparent_and_co2_core_is_opaque(soundings, spectra):
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, q = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)
    soundings.apply_sounding(state, T, q, T_surf=288.0)
    _, diag = lw(state)

    tau = diag["longwave_optical_depth_per_band"].values[:, 0, 0, :]  # (nz, nband)
    column_tau = tau.sum(axis=0)
    centres = spectra.band_centres(spectra.band_limits_of(lw))

    co2_core = column_tau[(centres > 630) & (centres < 700)]
    window = column_tau[(centres > 800) & (centres < 1180)]
    assert co2_core.min() > 1.0                 # opaque
    assert window.max() < co2_core.min()        # window is the transparent part


@pytest.mark.slow
def test_page2_removing_absorbers_opens_the_window(soundings):
    """A near-N2 atmosphere: OLR should approach the surface blackbody flux."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, q = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)

    soundings.apply_sounding(state, T, q, T_surf=288.0)
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 280e-6
    _, diag_earth = lw(state)
    olr_earth = float(diag_earth["upwelling_longwave_flux_in_air"].values[-1, 0, 0])

    soundings.apply_sounding(state, T, np.full_like(q, 1e-6), T_surf=288.0)
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 10e-6
    _, diag_thin = lw(state)
    olr_thin = float(diag_thin["upwelling_longwave_flux_in_air"].values[-1, 0, 0])

    sigma_T4 = 5.670374419e-8 * 288.0 ** 4
    assert olr_thin > olr_earth + 100.0        # the window opened
    assert olr_thin > 0.85 * sigma_T4          # approaching a bare surface
```

- [x] **Step 2: Run to verify**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -v -m slow -k page2`
Expected: PASS (helpers already exist). If either fails, the claim is wrong — reconcile
with the spec before writing the page. **As run: 2 passed.**

- [x] **Step 3: Write the page**

Create `docs/modelling-tour/02-window-measured.qmd` with the standard front matter
(title: `The window, measured`) and the setup include.

Content, in order:

1. `.callout-tip` "Where this fits": chapter 7 introduces `dτ = κ_λ ρ dz / cos θ` and
   notes that κ depends on wavelength. Chapter 6's lecture asked directly: *what would the
   window be if you had a pure nitrogen atmosphere?* This page answers it by turning the
   absorbers down.
2. Prose: optical depth is per-wavelength; the "atmospheric window" is not an assumption
   but a measurable range of wavenumbers where column optical depth is small.
3. A cell computing column optical depth per band on the standard sounding and plotting it
   against wavenumber on a log y-axis, with the window and CO₂ core shaded, plus a second
   panel showing `longwave_transmittance_per_band` at the surface layer. Reuse the sounding
   construction from page 1 verbatim (the reader may arrive here first).
4. Prose defining the window quantitatively from the plot: the bands where column τ < 1.
5. The pure-nitrogen experiment cell — a knob block:

```{pyodide}
CO2_PPM = 280.0      # <- try 10 (near-zero) or 10000
HUMIDITY_SCALE = 1.0 # <- try 0.001 for a dry atmosphere

state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = CO2_PPM * 1e-6
apply_sounding(state, T, np.maximum(q * HUMIDITY_SCALE, 1e-7), T_surf=T_SURF)
_, diag = lw(state)
olr_band, flux_density = spectral_olr(diag, limits)
sigma_T4 = 5.670374419e-8 * T_SURF ** 4
print(f"OLR = {olr_band.sum():.1f} W/m2   "
      f"({100 * olr_band.sum() / sigma_T4:.0f}% of the surface's {sigma_T4:.0f})")
```

6. Prose on the result: with both absorbers near zero the window widens until the
   atmosphere is nearly transparent and OLR approaches σT⁴ — there is no greenhouse effect
   left. This is the answer to the lecture's question.
7. `.callout-note` "climt craft": the state *contract*. Show `lw.input_properties` and
   explain that each entry fixes a name, units and dims, that climt converts units for you
   but will not guess names, and that `specific_humidity` is kg/kg while other gases are
   mole fractions.
8. Exercises — *Physics:* find the CO₂ concentration at which the 15 µm core's column τ
   first exceeds 1; explain why the window's τ barely responds to CO₂. *Code:* write a
   loop printing column τ per band for three CO₂ values; use `lw.input_properties` to
   discover the correct name for the CO₂ field rather than guessing.
9. `## Going deeper` → `../radiative-transfer/05-gas-overlap.qmd`.

- [x] **Step 4: Render**

```bash
cd docs && QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render modelling-tour/02-window-measured.qmd
```

- [x] **Step 5: Commit**

```bash
git add docs/modelling-tour/02-window-measured.qmd tests/test_modelling_tour.py
git commit -m "docs(tour): page 2, the window measured

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

### Log — the `resources:` mechanism is now verified in a real browser

Page 1's log left this open: the `pyodide: resources:` fetch of `_tour/*.py` had only
been checked at the level of the served files and the base64 vfs manifest, never by
actually booting the page. It has now been checked end to end, on this page, serving
`docs/_site` over `python3 -m http.server` and driving Chrome:

- the wheel installs (`climt 0.30.0 ready — has_fortran: False`);
- the autorun setup cell's `from soundings import ...` / `from spectra import ...`
  succeed after `sys.path.insert(0, "_tour")`, and it prints its line;
- both figure cells draw, so the `plt.show()` rule from page 1 holds in practice;
- editing the knob cell in place (`CO2_PPM = 10`, `HUMIDITY_SCALE = 0.001`) and
  pressing *Run Code* reproduces the native numbers to the last digit — OLR 354.4
  W m⁻² (91% of σT⁴), greenhouse 35.7, 980 cm⁻¹ of transparent bands.

Pages 3–6 can use the mechanism without re-verifying it.

### Log — the shared setup include renders ~200 lines of RCE code on every tour page

`{{< include ../_includes/climt-live-setup.qmd >}}` carries its cell with
`#| edit: false`, which makes quarto-live render it **read-only, not hidden**. So every
tour page opens with the full `integrate_with_snapshots` / `_draw_evolution` body —
roughly four screens of scrolling before the reader reaches "Where this fits". None of
the six tour pages calls any of it; they only want the imports and the wheel.

Not fixed here, because the include is shared with the radiative-transfer pages where the
code *is* the point. Options for **Task 14**: add `#| output: false`-style hiding to the
include's cell, or split a minimal `climt-live-boot.qmd` (imports only) for the tour and
leave the RCE helpers in the existing include. Page 1 has the same problem.

### Log — deviations from the step's literal content

- **Physics exercise 1 as specified is unanswerable.** The step asks students to "find the
  CO₂ concentration at which the 15 µm core's column τ first exceeds 1". There is none: at
  the *bottom* of the table's CO₂ axis the 630–700 band is already at τ = 33 on the
  standard sounding, and at τ = 20 with water driven to its own floor — so it is CO₂ at
  10 ppm doing this, not water vapour. Above the floor τ is very nearly linear in
  concentration (τ ≈ 2.02 × ppm^0.999, dry), which extrapolates to τ = 1 near **0.5 ppm**,
  a factor of 20 below anything the table resolves. Rewritten as a sweep over 10 / 100 / 1000 /
  10 000 ppm contrasting the core's hundred-fold response with the window's ~10% one, and
  a third exercise now poses the real puzzle the original was reaching for — *if the band
  is already opaque, how can more CO₂ change the OLR at all?* — pointing at page 3.
- **The 10 ppm floor is promoted to a teaching beat.** `co2_vmr_grid` runs 1e-5 … 1e-2, and
  cork silently clamps below it, so `CO2_PPM = 1` and `CO2_PPM = 10` give identical answers.
  The page says so explicitly rather than letting a student discover a "bug".
- **The transmittance panel plots two curves, not one.** The step asks for
  `longwave_transmittance_per_band` at the surface layer; on its own it reads as if the
  atmosphere were transparent everywhere, because it is a *single layer* on a 40-level
  grid. Plotted against the column transmittance `exp(−D τ_col)` it becomes a callout on
  reading `dims` before using a diagnostic — which feeds the craft note's third point.
- **The diffusivity factor is introduced here**, not left to Task 10. Converting τ to a
  transmittance needs it, chapter 7's `1 / cos θ` is the natural hook, and it forward-refs
  page 4 where the notes' `D = 2` starts to matter.

Numbers as built, `T_SURF = 288 K`, RH 80%, 280 ppm, nz = 40: OLR 246.7 W m⁻² (63% of
σT⁴ = 390.1); column τ = 0.43 / 0.22 / 0.68 in the three window bands (transmittance
49 / 70 / 32%), 580 in the CO₂ core, 1.6 × 10⁴ in the rotational band; window 800–1180
cm⁻¹, 380 cm⁻¹ wide. At 10 ppm and humidity × 0.001: OLR 354.4 (91%), greenhouse 35.7,
transparent bands 350–1400 cm⁻¹ spanning 980 cm⁻¹, and the CO₂ core still at τ = 20.

---

## Task 9: Page 3 — Where photons come from

**Files:**
- Create: `docs/modelling-tour/03-radiating-level.qmd`
- Test: `tests/test_modelling_tour.py` (extend)

**Interfaces:**
- Consumes: `spectra.emission_weight`, `spectra.tau_star_cumulative`,
  `spectra.brightness_temperature`.

- [x] **Step 1: Write the failing test**

Append to `tests/test_modelling_tour.py`:

```python
@pytest.mark.slow
def test_page3_brightness_temperature_matches_the_emission_level(soundings, spectra):
    """T_b(band) should be close to the air temperature where tau* = 1."""
    D = 1.66
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=60))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, q = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)
    soundings.apply_sounding(state, T, q, T_surf=288.0)
    _, diag = lw(state)

    limits = spectra.band_limits_of(lw)
    _, flux_density = spectra.spectral_olr(diag, limits)
    tb = spectra.brightness_temperature(flux_density, spectra.band_centres(limits))
    tau = diag["longwave_optical_depth_per_band"].values[:, 0, 0, :]

    # For each optically thick band, the emission-weighted temperature should
    # track the brightness temperature.
    for band in range(tau.shape[1]):
        if tau[:, band].sum() * D < 1.5:
            continue                       # optically thin: surface dominates
        w = spectra.emission_weight(tau[:, band], D)
        t_weighted = float((w * T).sum() / w.sum())
        assert abs(t_weighted - tb[band]) < 25.0
```

- [x] **Step 2: Run to verify**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -v -m slow -k page3`
Expected: PASS. If the tolerance is too tight for some band, widen it to 30 K and note
which band — do not weaken the test below 30 K without reconciling with the spec.

- [x] **Step 3: Write the page**

Create `docs/modelling-tour/03-radiating-level.qmd`, title `Where photons come from`.

Content:

1. `.callout-tip` "Where this fits": chapter 7 derives the weighting function
   `W ∝ τ* e^{−τ*}` and shows it peaks at τ* = 1 — the radiating level. It then notes
   that for a gray shell with ε ≠ 1 the radiating level "is not defined at all". This page
   shows why: there isn't one.
2. Prose restating the competition the lecture describes — go down and there are more
   molecules but more absorption above; go up and there is less absorption but almost
   nothing to emit. The product peaks in between.
3. A cell computing, per band, `emission_weight` and `tau_star_cumulative`, then plotting:
   left panel, the weighting function against pressure for four chosen bands (the CO₂
   core at 630–700, a window band at 800–980, the H₂O rotational band at 10–250, and the
   6.3 µm band at 1400–1600), with a log pressure axis; right panel, the temperature
   sounding with each band's τ* = 1 level marked and its brightness temperature plotted
   as a point, so the two visibly coincide.
4. Prose: the τ* = 1 levels span from near the surface (window) to the stratosphere (CO₂
   core). The "radiating height" of the shell model is an average over that spread.
5. A knob cell scaling humidity and re-drawing, showing emission levels rise as the
   atmosphere moistens — the mechanism behind the water vapour feedback, which page 6
   quantifies.
6. `.callout-note` "climt craft": `mid_levels` vs `interface_levels`. Explain that
   temperature and optical depth live on mid-levels (`nz` values) while fluxes live on
   interfaces (`nz+1` values), that the last interface index is the top of the atmosphere,
   and that this is why OLR is `[...][-1]`. Show both shapes:

```{pyodide}
print("air_temperature                     ", state["air_temperature"].values.shape)
print("longwave_optical_depth_per_band     ",
      diagnostics["longwave_optical_depth_per_band"].values.shape)
print("upwelling_longwave_flux_in_air      ",
      diagnostics["upwelling_longwave_flux_in_air"].values.shape)
```

7. Exercises — *Physics:* which band's radiating level moves most when humidity doubles,
   and why is it not the window? Estimate the radiating height in km using a 7.5 km scale
   height. *Code:* write a function returning the τ* = 1 pressure for a given band by
   interpolating `tau_star_cumulative`; apply it to all 14 bands and print a table.
8. `## Going deeper` → `../radiative-transfer/07-two-stream.qmd`.

- [x] **Step 4: Render**

```bash
cd docs && QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render modelling-tour/03-radiating-level.qmd
```

- [x] **Step 5: Commit**

```bash
git add docs/modelling-tour/03-radiating-level.qmd tests/test_modelling_tour.py
git commit -m "docs(tour): page 3, the weighting function and per-band radiating levels

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

### Log — Step 1's test was wrong, and being wrong is the page's best material

`test_page3_brightness_temperature_matches_the_emission_level` asserted
`|T_weighted − T_b| < 25 K` for every optically thick band. It fails, by 42 K at
700–800 cm⁻¹ and by 100 K at 1800–3250. Step 2 permits widening to 30 K; neither
failure fits under 30 K, so this is the reconciliation the step asks for.

**The error is real, one-signed, and explicable.** For every band, `T_b` comes out
*warmer* than the band-mean weighting function predicts — 0.7 K in the CO₂ core,
20 K in the rotational band, 36 K at 700–800. The cause is that
`longwave_optical_depth_per_band` is the **g-weighted mean** τ
(`component.py:313-316`), and within one band the eight g-points' k values span
2.7 × 10³ to 2.6 × 10⁹ (measured at 230 K, 100 hPa). The mean is dominated by the
strong lines, so it describes a band that is opaque throughout, while the OLR
actually escapes through the weak g-points from much lower, warmer air. The bands
with the smallest error are the genuinely saturated ones — the CO₂ core, where
every g-point is opaque at the same stratospheric temperature.

The 1800–3250 failure is a different bug and also worth keeping: `T_b = 300 K`
exceeds the 288 K surface, because `brightness_temperature` inverts a band-mean
flux density at the band *centre*, and Planck is nowhere near linear over 1450
cm⁻¹ of the shortwave tail.

Neither is a defect in the helpers, so the fix was to guard what is true. Four
tests replace the one:

- `test_page3_radiating_levels_are_a_distribution_not_a_height` — the spec's
  actual claim, and it lands harder than planned: the CO₂ core's τ* = 1 level is
  at 8 hPa, the twelve bands that have one span 2.8 → 993 hPa (a factor of 350),
  and **800–980 and 980–1080 never reach τ* = 1 at all**, so for them no
  radiating level exists in any sense.
- `test_page3_brightness_temperature_tracks_the_emission_level` — correlation
  > 0.9 over an 85 K spread, bias strictly one-signed against `T_weighted`
  (0.7 … 40.7 K) and one-signed with a single documented exception against
  `T(τ* = 1)` (−3.3 … 35.8 K).
- `test_page3_band_centre_brightness_temperature_fails_on_the_widest_band` —
  pins the artifact to exactly the last band.
- `test_page3_moistening_raises_the_emission_levels` — the knob.

**The one negative gap.** 1080–1180 gives `T_b` 3.3 K *colder* than `T(τ* = 1)`.
Its column τ is 0.68, so τ* crosses 1 inside the bottom model layer, at 993 hPa;
the "level" has degenerated into the ground, and what space sees is the surface
plus a veil of colder air. The construction fails at both ends — for saturated
bands because a band is not a wavelength, for transparent ones because there is
nothing to find. The page says this; the test allows exactly one negative.

### Log — deviations from Step 3's literal content

- **Step 3.3 asks for the τ* = 1 levels and the brightness temperatures to
  "visibly coincide".** They do not, per the above, so the right panel plots
  *both* markers per band joined by a segment, and the gap is the section that
  follows. Hiding it would have made the figure a lie.
- **A fifth cell was added**, printing the k spread per band straight from
  `lw._table["k_coefficients"]`. Without it the 36 K gap is an assertion; with it
  the reader can see the three-to-nine orders of magnitude that cause it. It
  follows page 2's precedent of reading `lw._table` in the open.
- **Code exercise 1 was rewritten.** The step asks the student to write the
  τ* = 1 interpolator, but the page needs it for its own figure, so it is given.
  The exercise now asks them to rank the bands by the gap and check that ranking
  against the k-spread column.
- **Physics exercises 1 and 3 were rewritten** after checking their answers.
  Exercise 1 as specified steered towards "it is not the window", which is false
  — the window bands move most, from no level at all to 730 hPa. Exercise 3 as
  specified said the OLR falls "through a couple of bands"; measured, a 200 → 180 K
  stratosphere moves the OLR by only 1.7 W m⁻², of which −1.9 is the CO₂ core
  alone and several water bands go slightly *positive* (a colder stratosphere is
  a drier one). Both now pose the measured behaviour as the question.

Numbers as built, nz = 60, 288 K, RH 80%, 280 ppm, D = 1.66: OLR 246.7 W m⁻²;
τ* = 1 at 2.8 hPa (10–250), 8.1 hPa (630–700), 30.2 hPa (1400–1600), 993 hPa
(1080–1180), none for 800–980 or 980–1080. At humidity × 4: OLR 204.1; the
rotational band rises to 1.3 hPa, 6.3 µm to 6.3 hPa, and the two window bands
acquire levels at 731.8 and 819.8 hPa — the window closing, which is page 6's
subject.

---

## Task 10: Page 4 — Gray equilibrium, tested

The strongest page: the notes' closed-form solution validated by a radiation code, then broken by real spectroscopy. **Verified during design** — heating rate rms 0.0012 K/day for the gray case, 2.33 K/day for the same profile at 14 bands.

**Files:**
- Create: `docs/modelling-tour/04-gray-equilibrium-tested.qmd`
- Test: `tests/test_modelling_tour.py` (extend)

**Interfaces:**
- Consumes: `soundings.analytic_gray_equilibrium`; the `tour_gray_lw` table (Task 3); the
  `diffusivity_factor` argument (Task 1).

- [x] **Step 1: Write the failing test**

Append to `tests/test_modelling_tour.py`:

```python
TAU_INF, D_NOTES, T_E = 4.0, 2.0, 255.0


def _gray_equilibrium_state(soundings, nz=60):
    lw = CorkLongwaveRadiation(optics="correlated_k", table="tour_gray_lw",
                               diffusivity_factor=D_NOTES)
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=nz))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, T_surf, _ = soundings.analytic_gray_equilibrium(p, ps, TAU_INF, T_E)
    soundings.apply_sounding(state, T, T_surf=T_surf)
    return lw, state, T, T_surf


@pytest.mark.slow
def test_page4_analytic_profile_is_an_equilibrium(soundings):
    """The notes' closed form must give ~zero heating and OLR = sigma T_e^4."""
    lw, state, _, _ = _gray_equilibrium_state(soundings)
    tendencies, diag = lw(state)
    H = tendencies["air_temperature"].values[:, 0, 0] * 86400.0    # K/day
    olr = float(diag["upwelling_longwave_flux_in_air"].values[-1, 0, 0])

    assert np.abs(H).max() < 0.05
    assert olr == pytest.approx(5.670374419e-8 * T_E ** 4, rel=2e-3)


@pytest.mark.slow
def test_page4_same_profile_is_not_an_equilibrium_for_a_real_spectrum(soundings):
    """The whole motivation for non-grey radiation, in one comparison."""
    lw_gray, state_gray, T, T_surf = _gray_equilibrium_state(soundings)
    H_gray = lw_gray(state_gray)[0]["air_temperature"].values[:, 0, 0] * 86400.0

    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=60))
    soundings.apply_sounding(state, T, np.full(T.shape, 1e-6), T_surf=T_surf)
    H_real = lw(state)[0]["air_temperature"].values[:, 0, 0] * 86400.0

    rms = lambda x: float(np.sqrt((x ** 2).mean()))
    assert rms(H_real) > 100 * rms(H_gray)
    assert np.abs(H_real).max() > 5.0
```

- [x] **Step 2: Run to verify it fails, then passes**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -v -m slow -k page4`
Expected: PASS once Tasks 1–4 are in. If `test_page4_analytic_profile_is_an_equilibrium`
fails with a non-zero heating rate, the table's diffusivity calibration and the
component's `diffusivity_factor` disagree — re-read the Global Constraints note about
`τ_notes = D · τ_normal`.

- [x] **Step 3: Write the page**

Create `docs/modelling-tour/04-gray-equilibrium-tested.qmd`, title
`Gray equilibrium, tested`.

Content:

1. `.callout-tip` "Where this fits": chapter 8 derives
   `T(τ) = T_e[(1 + τ_∞ − τ)/2]^{1/4}`, the skin temperature `2^{−1/4} T_e` and the
   surface discontinuity — and then stops. This page checks it, then breaks it.
2. A `.callout-important` titled "Which τ?" stating plainly that the notes' `τ_∞` is the
   **diffusivity-scaled** optical depth, that the notes use `D = 2` while cork defaults to
   the Elsasser value 1.66, and that this page passes `diffusivity_factor=2.0` and uses a
   table calibrated to match. Say why it matters: get it wrong and the heating rate does
   not vanish, for reasons that have nothing to do with the physics.
3. A cell prescribing the analytic profile and printing the analytic landmarks:

```{pyodide}
#| autorun: true
import numpy as np, matplotlib.pyplot as plt, sympl, climt
from climt import CorkLongwaveRadiation, get_default_state, get_grid
sympl.set_backend(climt.UnytBackend())

TAU_INF, T_E = 4.0, 255.0        # <- change TAU_INF and re-run
SIGMA = 5.670374419e-8

lw = CorkLongwaveRadiation(optics="correlated_k", table="tour_gray_lw",
                           diffusivity_factor=2.0)
state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=60))
p = state["air_pressure"].values[:, 0, 0]
ps = float(state["surface_air_pressure"].values.ravel()[0])
T, T_surf, tau_star = analytic_gray_equilibrium(p, ps, TAU_INF, T_E)
apply_sounding(state, T, T_surf=T_surf)

print(f"analytic skin temperature = {T_E * 0.5**0.25:.2f} K   (top level {T[-1]:.2f} K)")
print(f"analytic ground           = {T_surf:.2f} K")
print(f"surface-air discontinuity = {T_surf - T[0]:.2f} K")
```

4. Prose on the discontinuity: the ground is hotter than the air touching it, which is a
   real prediction of pure radiative equilibrium and the reason convection must exist.
5. The verification cell:

```{pyodide}
tendencies, diagnostics = lw(state)
H = tendencies["air_temperature"].values[:, 0, 0] * 86400.0
olr = float(diagnostics["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
print(f"OLR = {olr:.3f} W/m2   (sigma T_e^4 = {SIGMA * T_E**4:.3f})")
print(f"heating rate: max |H| = {np.abs(H).max():.4f} K/day")
```

6. Prose: the heating rate is zero to within a thousandth of a K/day. The analytic
   solution is not merely plausible — it *is* the equilibrium of the radiative transfer
   the model solves. State that this is a genuine validation, in both directions.
7. The demolition cell:

```{pyodide}
lw_real = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
state_real = get_default_state([lw_real], grid_state=get_grid(nx=1, ny=1, nz=60))
apply_sounding(state_real, T, np.full(T.shape, 1e-6), T_surf=T_surf)
H_real = lw_real(state_real)[0]["air_temperature"].values[:, 0, 0] * 86400.0

fig, ax = plt.subplots(figsize=(6, 5))
ax.plot(H, p / 100, label="gray table, D = 2", lw=2)
ax.plot(H_real, p / 100, label="14 real bands", lw=2)
ax.axvline(0, color="k", lw=0.8, ls=":")
ax.invert_yaxis(); ax.set_yscale("log")
ax.set(xlabel="longwave heating rate (K/day)", ylabel="pressure (hPa)")
ax.legend(); ax.set_title("The same profile, judged by two spectra")
print(f"rms heating: gray {np.sqrt((H**2).mean()):.4f}, "
      f"real {np.sqrt((H_real**2).mean()):.4f} K/day")
```

8. Prose — the payoff. The identical profile that the gray atmosphere calls a perfect
   equilibrium, the real spectrum calls wrong by several K/day, strongly cooling in some
   layers and heating in others. Everything the gray model concluded about lapse rate and
   surface temperature inherits that error. This is why the rest of the tour is non-grey.
9. `.callout-note` "climt craft": tendencies vs diagnostics. `lw(state)` returns both; the
   heating rate is a *tendency* — what the profile would do next — not a state variable.
   Note that `TendencyComponent`s return rates while `Stepper`s return updated states, and
   that this tranche never steps: every number here is an instantaneous rate.
10. Exercises — *Physics:* vary `TAU_INF` from 0.5 to 8 and record ground temperature;
    compare with `T_e(1 + τ_∞/2)^{1/4}`. Check the notes' claim that the radiative lapse
    rate approaches `g/4R` for large τ_∞ and explain why `c_p ≈ 3.5R < 4R` means Earth's
    atmosphere must convect. *Code:* set `diffusivity_factor=1.66` while keeping the same
    table and explain the resulting non-zero heating rate.
11. `## Going deeper` → `../radiative-transfer/01-why-nongrey.qmd`.

- [x] **Step 4: Render**

```bash
cd docs && QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render modelling-tour/04-gray-equilibrium-tested.qmd
```

- [x] **Step 5: Commit**

```bash
git add docs/modelling-tour/04-gray-equilibrium-tested.qmd tests/test_modelling_tour.py
git commit -m "docs(tour): page 4, the analytic gray equilibrium validated then broken

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

### Log — both tests passed first time; the numbers as built

Unlike Tasks 5 and 9, Step 1's tests needed no revision. nz = 60, `tour_gray_lw`,
`D = 2`, `τ_∞ = 4`, `T_e = 255 K`:

| quantity | value |
|---|---|
| D-scaled column τ | 4.0000 (target 4) |
| analytic skin / top model level | 214.43 / 214.70 K (top level is at 1.31 hPa, τ* = 0.005) |
| ground / bottom air / jump | 335.60 / 320.57 / 15.03 K |
| OLR vs σT_e⁴ | 239.854 vs 239.758 W m⁻² (+0.04%) |
| max \|H\| / rms H | 0.00193 / 0.00116 K/day |
| net flux constancy | 239.797 … 239.881 W m⁻² (spread 0.08) |
| downwelling at TOA | exactly 0 |

Same profile, `earth_low_res_lw`, `q = 10⁻⁶`: OLR **600.07** W m⁻², rms H
**2.327** K/day (a factor of 2006 on the gray rms), max \|H\| 16.74 K/day in the
bottom layer (heating — a 335.6 K ground radiating into near-transparent cold
air), 55 of 60 layers cooling, 0.4–0.5 K/day through the troposphere rising to
1.76 K/day at the lid.

`diffusivity_factor=1.66` on the `D = 2` table (Code exercise 1): OLR 263.5,
max \|H\| 0.653, cooling at the bottom and heating at the top. The τ_∞ the
profile *would* be an equilibrium for is 1.66 × 2 = 3.32.

### Log — Physics exercise 2 as specified reaches a false conclusion

Step 3.10 asks the student to "check the notes' claim that the radiative lapse
rate approaches `g/4R` for large τ_∞ and explain why `c_p ≈ 3.5R < 4R` means
Earth's atmosphere must convect." The premise is right and the conclusion is
backwards.

With τ linear in pressure — which is what `analytic_gray_equilibrium` assumes and
what `tour_gray_lw` delivers — `dT/dz = −(g/4R)·τ*/(1 + τ*)`, so the steepest
point is the surface and the limit is `g/4R = 8.545 K km⁻¹`. The dry adiabat is
`g/c_p = 9.771 K km⁻¹`. Because `c_p ≈ 3.5R` is *less* than `4R`, `g/4R` is
*less* than `g/c_p`, so the interior is statically **stable** — measured at
τ_∞ = 1, 4 and 20, the unstable fraction of the column is 0.00 in every case.
The instability is real but lives elsewhere: in the 15 K surface–air
discontinuity, and in pressure broadening.

Rewritten as two exercises. Exercise 2 now asks the student to derive `g/4R`,
compare it with `g/c_p`, discover that the standard argument does not go through
in this form, and find where the instability actually is. Exercise 3 is new and
supplies the resolution: with `τ* ∝ p²` the exponent doubles, the limit becomes
`g/2R = 17.09 K km⁻¹`, and at τ_∞ = 4 the surface rate is 13.67 K km⁻¹ with the
bottom 42% of the column unstable. Verified numerically at τ_∞ = 1, 4, 20 for
both exponents.

### Log — `unyt` named explicitly in the pyodide package list

Added `- "unyt"` to the `pyodide: packages:` block of this page **and of pages
1–3**, keeping the section's front matter uniform. unyt is already a hard
dependency of climt (`setup.py`), so micropip resolved it transitively before;
naming it makes the fast state backend a stated requirement rather than a
transitive accident, and it is pure Python so there is no wheel to build.

Measured natively with `NUMBA_DISABLE_JIT=1` (the Pyodide code path), 60 levels,
single column: `tour_gray_lw` 2.26 ms/call on `UnytBackend` against 2.80 ms/call
on `DataArrayBackend`; `earth_low_res_lw` 99.3 vs 99.8 ms/call. The saving is
state-wrapping overhead, so it is ~20% where the radiative transfer is cheap
(this page's single-band table, re-run by every exercise sweep) and in the noise
where it is not. Later tasks should keep the `- "unyt"` line when they copy the
front matter.

### Log — deviations from Step 3's literal content

- **Step 3.6 asks for prose saying the heating rate is "zero to within a
  thousandth of a K/day".** It is 0.0019, so the page says 0.002 and converts it
  to a timescale (one kelvin in eighteen months) that the demolition section then
  contrasts against the real spectrum's one kelvin in ninety minutes.
- **A `τ_column` print was added to the verification cell**, echoing the D-scaled
  column optical depth against `TAU_INF`. Without it the "Which τ?" callout is a
  warning the reader cannot check; with it, the calibration is visible in the
  output of the cell that depends on it.
- **A third Code exercise was added** on the flux profiles, because chapter 8's
  third analytic result — net flux constant with height — is the one the page
  never uses, and it is checkable in four lines (verified: 239.797 … 239.881).
- **The craft callout gained a paragraph on
  `air_temperature_tendency_from_longwave`**, which appears in the printed
  diagnostics list holding the same numbers as the tendency. A reader who sees
  both and sums them double-counts the radiation, so the page says why both
  exist.

---

## Task 11: Page 5 — The CO₂ knob

**Files:**
- Create: `docs/modelling-tour/05-co2-knob.qmd`
- Test: `tests/test_modelling_tour.py` (extend)

**Interfaces:**
- Consumes: Task 4 and Task 5 helpers. Uses `earth_low_res_lw` (14×8) — **not** the
  hi-res table: this page sweeps, and sweeps must stay cheap.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_modelling_tour.py`:

```python
def _olr_at_co2(lw, state, soundings, T, q, T_surf, co2_ppm):
    soundings.apply_sounding(state, T, q, T_surf=T_surf)
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = co2_ppm * 1e-6
    _, diag = lw(state)
    return float(diag["upwelling_longwave_flux_in_air"].values[-1, 0, 0])


@pytest.mark.slow
def test_page5_co2_doubling_forcing_is_canonical(soundings):
    """Students measure ~3.7 W/m2 per doubling themselves."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, q = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)

    olr = {c: _olr_at_co2(lw, state, soundings, T, q, 288.0, c)
           for c in (280, 560, 1120)}
    first, second = olr[280] - olr[560], olr[560] - olr[1120]

    assert 3.0 < first < 4.5           # canonical ~3.7 W/m2
    assert 3.0 < second < 4.5
    assert abs(first - second) < 1.0   # logarithmic, not linear


@pytest.mark.slow
def test_page5_forcing_comes_from_the_wings_not_the_core(soundings, spectra):
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, q = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)
    limits = spectra.band_limits_of(lw)
    centres = spectra.band_centres(limits)

    per_band = {}
    for co2 in (280, 1120):
        soundings.apply_sounding(state, T, q, T_surf=288.0)
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = co2 * 1e-6
        per_band[co2], _ = spectra.spectral_olr(lw(state)[1], limits)

    delta = per_band[280] - per_band[1120]
    core = delta[(centres > 630) & (centres < 700)].sum()
    wings = delta[((centres > 500) & (centres < 630))
                  | ((centres > 700) & (centres < 800))].sum()
    assert wings > core          # the saturated core cannot contribute much
```

- [ ] **Step 2: Run to verify**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -v -m slow -k page5`
Expected: PASS (design-time values were 3.59 and 3.86 W/m²).

- [ ] **Step 3: Write the page**

Create `docs/modelling-tour/05-co2-knob.qmd`, title `The CO₂ knob`.

Content:

1. `.callout-tip` "Where this fits": chapter 4 introduces radiative forcing, the Planck
   response and climate sensitivity, and takes the forcing for a doubling of CO₂ as given.
   This page computes it.
2. A cell sweeping CO₂ across the table's runtime axis, plotting OLR against log₁₀(CO₂),
   with the 280→560 and 560→1120 forcings annotated. Use 12 points, spaced
   `np.logspace(np.log10(10), np.log10(10000), 12)` — at ~0.12 s/call that is under 2 s.
   State the cost in the prose so readers know why 12 and not 200.
3. Prose: the curve is a straight line against log CO₂. State the measured forcing per
   doubling and compare with the canonical 3.7 W/m². Note that they just measured, from
   line-by-line-derived absorption data, one of the most consequential numbers in climate
   science.
4. A per-band decomposition cell — bar chart of `ΔOLR` per band between 280 and 1120 ppm,
   with the core and wings shaded differently.
5. Prose on band saturation: the 15 µm core is already opaque, so adding CO₂ there changes
   nothing; the forcing comes from the wings, where τ is near 1 and more absorber still
   moves the emission level. This is *why* forcing is logarithmic rather than linear.

   **Open this section on the contrast page 2 set up, using page 2's own numbers.** Column
   optical depth in the core is very nearly *linear* in CO₂ — measured there as
   τ ≈ 2.02 × ppm^0.999 in a dry column, so 10 → 10 000 ppm multiplies τ by roughly a
   thousand, from 20 to 20 000. The forcing over that same range is a handful of W m⁻² per
   doubling, near-constant. A thousand-fold change in the absorber buys a *logarithmic*
   change in the energy budget, and the gap between those two statements is the whole
   physical content of the page. Resolve it with page 3's weighting function: the OLR
   depends on the *temperature at which* τ\* = 1 sits, and τ\* = 1 moves logarithmically
   in height when τ scales linearly with concentration. Forcing tracks where the emission
   level is, not how opaque the column is.

   Add a guard test for the linear half of that claim, so the page cannot assert it
   loosely — fit `log(tau_core)` against `log(co2_ppm)` over 30–10 000 ppm and assert the
   exponent is within a few percent of 1, alongside the existing constant-per-doubling
   assertion in `test_page5_forcing_is_logarithmic`. The two together are the beat.
6. `.callout-note` "climt craft": sweeping a parameter. Emphasise that the state is
   mutated in place and re-passed, that diagnostics must be collected inside the loop
   because the next call overwrites them, and that the sounding is re-applied each
   iteration so nothing accumulates.
7. Exercises — *Physics:* combine the measured forcing with a Planck feedback of
   −3.3 W m⁻² K⁻¹ to estimate the no-feedback climate sensitivity; compare with the ~1.2 K
   in the notes. Repeat the sweep at RH = 0 and explain why the forcing changes.
   *Code:* rewrite the sweep to also record window and core OLR, and plot all three.
8. `## Going deeper` → `../radiative-transfer/04-correlated-k.qmd`.

- [ ] **Step 4: Render, then commit**

```bash
cd docs && QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render modelling-tour/05-co2-knob.qmd
```

```bash
git add docs/modelling-tour/05-co2-knob.qmd tests/test_modelling_tour.py
git commit -m "docs(tour): page 5, measuring the CO2 forcing band by band

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 12: Page 6 — Water vapour, and the limit

**Files:**
- Create: `docs/modelling-tour/06-water-vapour-limit.qmd`
- Test: `tests/test_modelling_tour.py` (extend)

**Interfaces:**
- Consumes: Task 4 and Task 5 helpers; adds `CorkShortwaveRadiation`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_modelling_tour.py`:

```python
@pytest.mark.slow
def test_page6_fixed_rh_suppresses_olr_relative_to_fixed_q(soundings):
    """The water vapour feedback: OLR rises more slowly when moisture responds."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    _, q_ref = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)

    def olr(T_surf, fixed_q):
        T, q = soundings.lapse_rate_sounding(p, ps, T_surf=T_surf, rh=0.8)
        soundings.apply_sounding(state, T, q_ref if fixed_q else q, T_surf=T_surf)
        return float(lw(state)[1]["upwelling_longwave_flux_in_air"].values[-1, 0, 0])

    d_fixed_q = olr(298.0, True) - olr(288.0, True)
    d_fixed_rh = olr(298.0, False) - olr(288.0, False)
    assert d_fixed_q > d_fixed_rh > 0        # moisture damps the OLR response


@pytest.mark.slow
def test_page6_olr_saturates_at_high_surface_temperature(soundings):
    """Approach to the Simpson-Nakajima limit on saturated soundings."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])

    def olr(T_surf):
        T, q = soundings.lapse_rate_sounding(p, ps, T_surf=T_surf, rh=1.0)
        soundings.apply_sounding(state, T, q, T_surf=T_surf)
        return float(lw(state)[1]["upwelling_longwave_flux_in_air"].values[-1, 0, 0])

    warm = olr(310.0) - olr(300.0)
    hot = olr(340.0) - olr(330.0)
    assert hot < warm            # the OLR response flattens
    assert olr(340.0) < 400.0    # nowhere near sigma T^4 = 757 W/m2
```

- [ ] **Step 2: Run to verify**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -v -m slow -k page6`
Expected: PASS. If `test_page6_olr_saturates` fails, check the sounding's `q_floor` and
that H₂O VMR stays inside the table's axis (max 1.0, i.e. q ≈ 0.38 kg/kg).

- [ ] **Step 3: Write the page**

Create `docs/modelling-tour/06-water-vapour-limit.qmd`, title
`Water vapour, and the limit`.

Content:

1. `.callout-tip` "Where this fits": chapter 4 lists water vapour as the strongest
   positive feedback. This page draws it, then follows it to the point where it stops
   being a feedback and becomes a limit.
2. A cell computing OLR against surface temperature twice — fixed RH and fixed q — over
   280–310 K in 8 steps each, plotted together with the gap shaded.
3. Prose: the gap between the curves *is* the feedback. Quantify it as the difference in
   dOLR/dTₛ and connect it to the feedback parameter in the notes.
4. Prose on the window's role: OLR stays nearly linear in Tₛ over the Earth-like range
   because the window keeps radiating from near the surface even as the rest of the
   spectrum saturates.
5. A cell pushing saturated soundings from 300 K to 350 K and plotting OLR, showing the
   flattening, with a per-band panel showing the window closing.
6. Prose: this is the Simpson–Nakajima limit. Above it, absorbed sunlight cannot be
   balanced by longwave emission at any surface temperature, and the ocean evaporates
   entirely — the runaway greenhouse. Note the model is being pushed to the edge of its
   table and say so honestly.
7. A shortwave cell introducing `CorkShortwaveRadiation` with `earth_low_res_sw`, showing
   absorbed shortwave against the OLR curve so the crossing point — the equilibrium — is
   visible, and the absence of a crossing above the limit is visible too.
8. `.callout-note` "climt craft": running two components on one state. Explain that
   `get_default_state([lw, sw], ...)` builds a state satisfying both, that each returns
   its own tendencies and diagnostics, and that summing tendencies is what a time
   integrator would do — which is exactly where the next tranche begins.
9. Exercises — *Physics:* find the absorbed shortwave flux above which no equilibrium
   exists. Explain why raising CO₂ does not produce a runaway but raising insolation does.
   *Code:* write a function returning equilibrium Tₛ by bisection on `OLR(Tₛ) − ASR`, and
   show it fails to converge above the limit.
10. A closing `.callout-note` "Where next" saying the tour continues with chapters 10–12 —
    turbulent heat exchange, dry convection and radiative-convective equilibrium — where
    time integration finally enters.
11. `## Going deeper` → `../radiative-transfer/08-multiplanet.qmd`.

- [ ] **Step 4: Render, then commit**

```bash
cd docs && QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render modelling-tour/06-water-vapour-limit.qmd
```

```bash
git add docs/modelling-tour/06-water-vapour-limit.qmd tests/test_modelling_tour.py
git commit -m "docs(tour): page 6, the water vapour feedback and the runaway limit

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 13: The hi-res spectrum table, its validation gate, and wiring pages 1–3

Everything so far runs on the shipped 14×8 table. This task adds the ~100-band table and switches pages 1–3 onto it **with fallback**, so a generation failure degrades rather than breaks.

**Files:**
- Create: `scripts/generate_tour_spectrum_table.py`
- Create: `docs/modelling-tour/_data/earth_spectrum_lw.npz`
- Create: `docs/modelling-tour/_tour/tables.py`
- Create: `tests/test_spectrum_table.py`
- Modify: `docs/modelling-tour/01-emissivity-spectrum.qmd`, `02-window-measured.qmd`,
  `03-radiating-level.qmd`

**Interfaces:**
- Produces: `tables.spectrum_table(prefer_hires=True) -> str` returning a table name or
  path usable as `CorkLongwaveRadiation(table=...)`, falling back to
  `"earth_low_res_lw"` when the hi-res asset is unavailable.

- [ ] **Step 1: Write the generation script**

Create `scripts/generate_tour_spectrum_table.py`. It must:

- Accept `--output`, `--ngpt` (default 4), `--nsub` (default 7, sub-divisions per existing
  band) and `--co2-nodes` (default 5).
- Build the band edge list by subdividing each of the 14 shipped bands into `--nsub`
  equal-wavenumber pieces, so **every existing edge is preserved** (10, 250, 350, 500,
  630, 700, 800, 980, 1080, 1180, 1250, 1400, 1600, 1800, 3250). 14 × 7 = 98 bands.
- Shell out to `scripts/generate_cork_tables_linepyline.py --scenario earth --kind lw`
  with `--bands` set to that edge list and `--ngpt` as given.
- Thin the CO₂ axis to `--co2-nodes` log-spaced nodes over 10–10 000 ppm.
- Print the resulting file size and sha256.

- [ ] **Step 2: Generate the table**

```bash
conda run -n linepyline python scripts/generate_tour_spectrum_table.py \
    --output /tmp/earth_spectrum_lw.nc --ngpt 4 --nsub 7 --co2-nodes 5
conda run -n climt python scripts/convert_ck_table_to_npz.py /tmp/earth_spectrum_lw.nc
mkdir -p docs/modelling-tour/_data
cp /tmp/earth_spectrum_lw.npz docs/modelling-tour/_data/
ls -lh docs/modelling-tour/_data/earth_spectrum_lw.npz
```

Expected: ~5–6 MB. If it exceeds 8 MB, reduce `--ngpt` to 2 and re-run before proceeding.

- [ ] **Step 3: Write the validation gate**

Create `tests/test_spectrum_table.py`:

```python
"""Validation gate for the Modelling Tour's high-resolution spectrum table.

The hi-res table exists to make pages 1-3 look like a real IRIS spectrum. It is
only allowed to do that if it agrees with the shipped 14-band reference on the
quantities the tour states numerically. Band edges are a strict refinement of
the 14-band grid, so aggregation is exact summation with no interpolation error.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest
import sympl

import climt
from climt import CorkLongwaveRadiation, get_default_state, get_grid

REPO_ROOT = Path(__file__).resolve().parent.parent
HIRES = REPO_ROOT / "docs/modelling-tour/_data/earth_spectrum_lw.npz"
TOUR = REPO_ROOT / "docs/modelling-tour/_tour"

pytestmark = pytest.mark.skipif(not HIRES.exists(),
                                reason="hi-res table not generated")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, TOUR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(autouse=True)
def _unyt_backend():
    sympl.set_backend(climt.UnytBackend())


def _column(table, nz=40, co2_ppm=280.0):
    soundings = _load("soundings")
    lw = CorkLongwaveRadiation(optics="correlated_k", table=table)
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=nz))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, q = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)
    soundings.apply_sounding(state, T, q, T_surf=288.0)
    if "mole_fraction_of_carbon_dioxide_in_air" in state:
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = co2_ppm * 1e-6
    return lw, lw(state)


@pytest.mark.slow
def test_band_edges_refine_the_reference_grid():
    spectra = _load("spectra")
    ref, _ = _column("earth_low_res_lw")
    hi, _ = _column(str(HIRES))
    ref_edges = set(np.round(spectra.band_limits_of(ref).ravel(), 3))
    hi_edges = set(np.round(spectra.band_limits_of(hi).ravel(), 3))
    assert ref_edges <= hi_edges, "every 14-band edge must survive in the hi-res grid"


@pytest.mark.slow
def test_broadband_olr_agrees_within_2_watts():
    _, (_, diag_ref) = _column("earth_low_res_lw")
    _, (_, diag_hi) = _column(str(HIRES))
    olr_ref = float(diag_ref["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    olr_hi = float(diag_hi["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    assert abs(olr_hi - olr_ref) < 2.0


@pytest.mark.slow
def test_aggregated_per_band_olr_agrees_within_5_percent():
    spectra = _load("spectra")
    ref, (_, diag_ref) = _column("earth_low_res_lw")
    hi, (_, diag_hi) = _column(str(HIRES))
    ref_limits = spectra.band_limits_of(ref)
    hi_limits = spectra.band_limits_of(hi)
    olr_ref, _ = spectra.spectral_olr(diag_ref, ref_limits)
    olr_hi, _ = spectra.spectral_olr(diag_hi, hi_limits)

    hi_centres = spectra.band_centres(hi_limits)
    for band, (lo, high) in enumerate(ref_limits):
        mask = (hi_centres > lo) & (hi_centres < high)
        aggregated = olr_hi[mask].sum()
        assert aggregated == pytest.approx(olr_ref[band], rel=0.05), (
            f"band {band} ({lo:.0f}-{high:.0f} cm^-1) disagrees")


@pytest.mark.slow
def test_heating_rate_profile_agrees_within_half_a_kelvin_per_day():
    _, (tend_ref, _) = _column("earth_low_res_lw")
    _, (tend_hi, _) = _column(str(HIRES))
    H_ref = tend_ref["air_temperature"].values[:, 0, 0] * 86400.0
    H_hi = tend_hi["air_temperature"].values[:, 0, 0] * 86400.0
    assert np.abs(H_hi - H_ref).max() < 0.5


@pytest.mark.slow
def test_co2_forcing_survives_the_thinned_axis():
    """Thinning the CO2 axis 10 -> 5 nodes must not distort the forcing."""
    def olr(table, ppm):
        _, (_, diag) = _column(table, co2_ppm=ppm)
        return float(diag["upwelling_longwave_flux_in_air"].values[-1, 0, 0])

    f_ref = olr("earth_low_res_lw", 280) - olr("earth_low_res_lw", 560)
    f_hi = olr(str(HIRES), 280) - olr(str(HIRES), 560)
    assert abs(f_hi - f_ref) < 0.3
```

- [ ] **Step 4: Run the gate**

Run: `conda run -n climt python -m pytest tests/test_spectrum_table.py -v -m slow`
Expected: all pass. **If `test_aggregated_per_band_olr_agrees_within_5_percent` or
`test_heating_rate_profile_agrees_within_half_a_kelvin_per_day` fails, 4 g-points is too
coarse.** Regenerate at `--ngpt 8 --nsub 4` (56 bands, same file size) and re-run. Do not
relax the tolerances.

- [ ] **Step 5: Write the table locator**

Create `docs/modelling-tour/_tour/tables.py`:

```python
"""Locate the spectrum table, in the browser and natively.

The high-resolution table is too large for the wheel, so it ships as a static
asset beside these pages. The docs site is same-origin with the page, so Pyodide
can fetch it directly -- unlike a GitHub release asset, which sends no CORS
headers.

Pages call :func:`spectrum_table` and pass the result straight to
``CorkLongwaveRadiation(table=...)``. If the asset is missing the shipped
14-band table is returned instead, so a page degrades to a coarser spectrum
rather than failing.
"""
import os

FALLBACK = "earth_low_res_lw"
ASSET = "earth_spectrum_lw.npz"


def spectrum_table(prefer_hires=True, base_url="_data"):
    """Return a table name or path for CorkLongwaveRadiation.

    Args:
        prefer_hires: if False, always return the shipped 14-band table.
        base_url: directory (native) or URL prefix (browser) holding the asset.

    Returns:
        A table name or filesystem path.
    """
    if not prefer_hires:
        return FALLBACK

    local = os.path.join(base_url, ASSET)
    if os.path.isfile(local):
        return local

    try:                                    # browser: fetch into the Pyodide FS
        import pyodide_http                 # noqa: F401
    except ImportError:
        pass
    try:
        from pyodide.http import open_url
        data = open_url(f"{base_url}/{ASSET}")
        with open(ASSET, "wb") as handle:
            handle.write(data.getvalue().encode("latin-1"))
        return ASSET
    except Exception:
        return FALLBACK
```

- [ ] **Step 6: Switch pages 1–3 onto the hi-res table**

In each of `01-emissivity-spectrum.qmd`, `02-window-measured.qmd` and
`03-radiating-level.qmd`, replace the component construction

```python
lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
```

with

```python
TABLE = spectrum_table()          # ~100 bands if available, 14 otherwise
lw = CorkLongwaveRadiation(optics="correlated_k", table=TABLE)
print(f"using {lw.num_longwave_bands} bands")
```

Add a `.callout-note` to page 1 only, titled "A bigger download", stating that this page
fetches a ~5 MB spectral table once (cached afterwards) to draw a finely-resolved
spectrum, and that pages 4–6 do not need it.

Add `spectrum_table` to the helper imports in `docs/_includes/climt-live-setup.qmd`'s
boot cell, or import it in each page's setup cell — whichever matches how `soundings` and
`spectra` were wired in Tasks 7–9. Keep it consistent across the three pages.

- [ ] **Step 7: Re-render pages 1–3 and re-run the page tests**

```bash
cd docs && QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render modelling-tour/
```

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py tests/test_spectrum_table.py -v -m slow`
Expected: all pass. Page tests still pin `earth_low_res_lw` explicitly, so they are
unaffected by the switch.

- [ ] **Step 8: Commit**

```bash
git add scripts/generate_tour_spectrum_table.py docs/modelling-tour/_data/earth_spectrum_lw.npz docs/modelling-tour/_tour/tables.py tests/test_spectrum_table.py docs/modelling-tour/0[123]-*.qmd
git commit -m "feat(tour): high-resolution spectrum table with a validation gate

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 14: Static fallbacks, full render, and link check

**Files:**
- Create: `docs/modelling-tour/_artifacts/generate.py`
- Create: `docs/modelling-tour/_artifacts/*.png`
- Modify: the six page files (add fallback images)

- [ ] **Step 1: Write the artifact generator**

Create `docs/modelling-tour/_artifacts/generate.py` producing one PNG per page, reusing
`_tour/soundings.py` and `_tour/spectra.py` so the figures cannot drift from the pages'
physics. Name them `01-spectrum.png` … `06-runaway.png`.

- [ ] **Step 2: Generate and inspect**

```bash
conda run -n climt python docs/modelling-tour/_artifacts/generate.py
ls -la docs/modelling-tour/_artifacts/
```

Expected: six PNGs. Open at least two and confirm they are not blank.

- [ ] **Step 3: Add fallbacks to the pages**

On each page, immediately after the main figure cell, add a collapsed callout:

```markdown
::: {.callout-note collapse="true"}
## If the live cell did not run

![](_artifacts/01-spectrum.png)

This is the figure the cell above produces. Live cells need WebAssembly; if your
browser blocks it, the static version is here.
:::
```

- [ ] **Step 4: Full render**

```bash
cd docs && QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render
```

Expected: the whole site renders with no error, and `_site/modelling-tour/` contains
`index.html` plus six page HTML files.

- [ ] **Step 5: Check every cross-link resolves**

```bash
cd docs && grep -ohE '\]\((\.\./)?[a-z0-9./-]+\.qmd\)' modelling-tour/*.qmd | sort -u | \
  sed -E 's/^\]\(//; s/\)$//' | while read -r link; do
    target=$(cd modelling-tour && realpath -m "$link")
    [ -f "$target" ] || echo "BROKEN: $link"
  done
echo "link check done"
```

Expected: `link check done` with no `BROKEN:` lines.

- [ ] **Step 6: Run the full test suite**

Run: `conda run -n climt python -m pytest tests/ -v -m "not slow"`
Then: `conda run -n climt python -m pytest tests/test_modelling_tour.py tests/test_spectrum_table.py tests/test_cork_diffusivity.py -v -m slow`
Expected: all pass.

Also run the artifact staleness gate, which pytest does not cover and which the docs
workflow fails the branch on:

```bash
conda run -n climt python scripts/build_experiments.py --check
```

Expected: exit 0, no output. If it reports anything stale, a cork source edit somewhere
in Tasks 1–13 skipped its regeneration step — go back and run `make experiments` for the
manifest concerned rather than fixing it up here.

- [ ] **Step 7: Commit**

```bash
git add docs/modelling-tour/_artifacts docs/modelling-tour/*.qmd
git commit -m "docs(tour): static figure fallbacks and full-site render

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Self-review

**Spec coverage.** Positioning and site architecture → Task 6. Six pages → Tasks 7–12,
each anchored to the chapter the spec names. Two tables with distinct jobs → Tasks 3 and
13. `.npz` requirement and same-origin hosting → Task 13. Prerequisites → Task 1 in this
plan; the `NpzFile` fix is external and explicitly flagged. Notation reconciliation →
Tasks 1–3 plus page 4's `.callout-important`. Page anatomy → applied uniformly in Tasks
7–12. Testing strategy, including the validation gate → Tasks 4, 5, 13, 14. Craft thread →
one `.callout-note` per page, cumulative in the order the spec tabulates. Out-of-scope
items are named in the landing page and page 6's closing callout.

**Type consistency.** `spectral_olr` returns `(olr_per_band, flux_density)` and is called
that way in Tasks 7, 11, 13. `band_limits_of(component)` is used consistently rather than
reaching into `_table` from page code. `analytic_gray_equilibrium` returns
`(T, T_surf, tau_star)` and is unpacked as three values in Task 10 and the Task 4 test.
`emission_weight` and `tau_star_cumulative` both take `(tau_layer, diffusivity_factor)`
positionally in Tasks 5 and 9. `spectrum_table()` returns a single string in Task 13.
`diffusivity_factor` is the argument name throughout — not `D`, which is used only as a
local variable and in prose.

**Known soft spots, deliberately left to the implementer:** Task 13 Step 6 says to wire
`spectrum_table` "whichever matches how `soundings` and `spectra` were wired in Tasks
7–9", because that choice is made in Task 7 and depends on whether the shared include or
per-page imports proves cleaner under Pyodide. Task 14 Step 1 describes the artifact
generator's outputs and constraints rather than its full source, since it is mechanical
reuse of figure code written in Tasks 7–12.

---

## Experiment log — 2026-08-15: the live-RCE stratosphere

Triggered by the live-RCE page (`docs/radiative-transfer/09-live-rce.qmd`), not by a
task here, but the findings constrain how the tour's pages may talk about skin
temperature, so they are recorded with the rest of the CORK work.

**Symptom.** The isothermal stratosphere was never visible on the live-RCE figure, and
the non-grey column showed no isothermal layer at all — cooling monotonically to 120 K
at 3.9 hPa and still falling after 2000 simulated days.

**Two distinct causes, one real and one not.**

1. *The axis.* `_draw_evolution` in `docs/_includes/climt-live-setup.qmd` used
   `invert_yaxis()` on a linear pressure scale. The gray column's isothermal layer
   spans 19.9 → 3.9 hPa: **1.6% of a linear axis, ~19–29% of a log axis.** Fixed by
   switching to log-p over the interface range (0.2 hPa top, so no flux point clips).

2. *The non-grey profile is not a bug.* Verified three ways:
   - `CorkLongwaveRadiation` with the **gray** table settles at **214.71 K**, against
     the analytic skin temperature `(OLR/2σ)^(1/4) = 214.48 K`, isothermal to 0.85 K
     across the top two levels. The transport honours the skin-temperature limit.
   - **RRTMG-LW in the identical dry configuration** also cools monotonically to a very
     cold lid (147 K) with a *steeper* top gradient than CORK (15.6 K vs 7.5 K across
     the top two levels). Surface temperature and OLR agree to 0.3 K / 0.4 W m⁻².
   - The committed moist-RCE artifact (`docs/experiments/2026-06-05-cork-co2-bands`,
     "after: 14-band" panel) shows the same structure in both codes.

**Mechanism.** The skin-temperature cancellation (`eps_b = k_b * m` cancels between
absorption and emission, so T is independent of layer mass) is correct, but T is
independent of *height* only if `F_up,b` is. That needs the atmosphere thin **in the
absorbing bands**. Band 5 (630–700 cm⁻¹, CO₂ 15 µm core) carries **96% of the topmost
layer's opacity at τ = 1.29 in that layer alone**, so `F_up,5` is emitted locally by
ever-colder layers, falls with height, and drags the skin temperature down with it.
Optically thin in the windows, optically thick where it counts. Dropping CO₂ to
1e-5 mol/mol (the table's `co2_vmr` floor) makes the top thin and restores a finite
equilibrium temperature for the top two layers, which at default CO₂ have **none** at
any temperature down to 100 K.

**Consequence for the tour.** Do not assert an isothermal stratosphere on any non-grey
page. It is a property of the *gray* pages only. The live-RCE annotation is now
conditional on `|T[-1] - T[-2]| < 2 K`, which separates the two cases cleanly at both
33 days and equilibrium.

**Open follow-up (not chased).** CORK runs **27–36 K colder than RRTMG through the
stratosphere** while matching it in the troposphere. Both show the same qualitative
structure, so this does not affect the conclusion above, but the gap is real. Most
likely stratospheric **ozone** longwave absorption, which RRTMG takes from the default
state and `earth_low_res_lw` has no term for — see
`docs/superpowers/plans/2026-06-01-cork-ozone-lw.md` — plus coarser band structure.

**Probes** (kept, all under `scripts/experiments/`): `live_rce_axis_check.py`,
`render_live_rce_figure.py` (execs the shipped `{pyodide}` cell so the native render is
the page's own code), `skin_temperature_probe.py` (`--budget`, `--sweep`, `--co2`),
`live_rce_cork_vs_rrtmg.py`.

**Rejected approach.** A block sweep that holds levels `k_lo..top` at one temperature is
confounded when the block is optically thick: the flux arriving at its top is then
emitted by the block itself and tracks `B(T_top)`, so no zero crossing can exist and the
result says nothing about the physics. Only meaningful once the block is genuinely thin.
