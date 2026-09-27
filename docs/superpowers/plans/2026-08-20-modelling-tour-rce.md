# Modelling Tour of the Climate System — Radiative-Convective Equilibrium Tranche — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add six in-browser pages `07`–`12` to `docs/modelling-tour/`, which introduce time integration, a slab surface, a turbulent boundary layer, dry convective adjustment and moist convection, and arrive at radiative-convective equilibrium — and retire `docs/radiative-transfer/09-live-rce.qmd`, whose demo page 07 replaces properly.

**Architecture:** Same shape as tranche 1. Each page is a `format: live-html` Quarto page whose computation lives in importable, natively-tested Python under `docs/modelling-tour/_tour/`; page cells stay thin. Four new helper modules carry this tranche's machinery: `stepping.py` (the time loops and the evolution figure), `states.py` (save/load a sympl state as `.npz` with provenance), `budgets.py` (TOA and surface energy budgets), and `assets.py` (browser/native asset resolution, lifted out of `tables.py`). Pages 11 and 12 load shipped equilibrium states rather than spinning one up; a residual test in CI guards those states against staleness.

**Tech Stack:** Python, numpy, matplotlib, sympl (`UnytBackend`, `UnytTimeDelta`, `AdamsBashforth`), climt, Quarto, `quarto-live` (Pyodide + micropip), pytest.

**Spec:** `docs/superpowers/specs/2026-08-20-modelling-tour-rce-design.md` is the source of truth. Where this plan and the spec disagree, stop and reconcile. Tranche 1's spec and plan (`2026-08-12-modelling-tour-radiation-design.md`, `2026-08-12-modelling-tour-radiation.md`) are the pattern this one continues; read tranche 1's plan Tasks 7–12 before writing any page here.

## Global Constraints

- **Branch:** all work on `feature/modelling-tour-rce`, branched off `develop`. Create it in Task 1.
- **Conda env:** run every Python/pytest command in the `climt` conda env, e.g. `conda run -n climt python -m pytest ...`.
- **Backend:** `sympl.set_backend(climt.UnytBackend())` on every page and in every test.
- **Timestep type:** `climt.UnytTimeDelta`, never `datetime.timedelta`. Under `UnytBackend` a plain `timedelta` raises `unyt.exceptions.UnitOperationError: The <ufunc 'add'> operator for unyt_arrays with units 'degK' ... and 'degK/s' ... is not well defined`, because `timedelta.total_seconds()` returns a bare float that does not cancel the tendency's `/s`. `UnytTimeDelta.total_seconds()` returns a `unyt_quantity` in seconds, which does.
- **Browser components only** on pages: `CorkLongwaveRadiation`, `SlabSurface`, `SimpleBoundaryLayer`, `DryConvectiveAdjustment`, `GridScaleCondensation`, `EmanuelConvectionPython`. `RRTMGLongwave`, `RRTMGShortwave`, the Fortran `EmanuelConvection`, `SimplePhysics`, `BergerSolarInsolation` and `DcmipInitialConditions` are compiled and absent under Pyodide — never import them in a page.
- **`CorkShortwaveRadiation` is not used anywhere in this tranche.** Every page prescribes the absorbed shortwave flux at the surface: `SOLAR = 240.0` W m⁻², written into `downwelling_shortwave_flux_in_air` at index 0 with `upwelling_shortwave_flux_in_air` zeroed. This is deliberate (spec, *Deferred*); do not "improve" it by adding a shortwave component.
- **Pyodide has no numba.** `NUMBA_DISABLE_JIT=1` is the correct native proxy for the browser, and every performance number in this plan was measured that way unless it says otherwise. Any code path that only works when numba compiles it away is broken in the browser — Task 1 exists because of exactly that.
- **Wind:** every configuration with a `SimpleBoundaryLayer` carries `stepping.wind_relaxation` (pages 08–12, and both shipped equilibria). A column has no dynamics, so nothing drives a wind, and the boundary layer's own drag removes any wind prescribed as an initial condition — measured, three initialisations spanning 0–10 m s⁻¹ all end at 0.000 m s⁻¹. Default target 5 m s⁻¹ over a 24 h timescale, `z0 = 1e-3`. A page that omits it is reporting the surface fluxes of a windless planet.
- **Grid:** `nz = 28`, `nx = ny = 1`, single column, on every page in this tranche. (Tranche 1 used varying `nz` per page; this tranche does not, because the shipped equilibrium states pin it.)
- **Test marker:** anything over ~5 s gets `@pytest.mark.slow` (see `tox.ini`); CI runs `-m "not slow"`.
- **Tests are parametrised over the table the page declares.** This is the tranche 1 lesson that cost the most: pages 1–3 ran on the 56-band table while every test ran on the 14-band one, so every number in those pages was guarded by nothing. A test for page N asserts on page N's table.
- **Commit after every task.** End commit messages with `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`.
- **Units:** `specific_humidity` in kg/kg; all other gases as mole fractions (mol/mol), e.g. `mole_fraction_of_carbon_dioxide_in_air`. `air_pressure` in Pa in state, but `EmanuelConvectionPython` declares its inputs in mbar — sympl converts, you do not.
- **Every number a page quotes comes from a cell the reader can run,** or from Task 8's measurement log with the configuration named beside it.

## What was measured before this plan was written

These are real measurements on the reference machine, not estimates, and they replace the spec's cost table where they disagree. Single column, `nz = 28`, `dt = 1 h`, `NUMBA_DISABLE_JIT=1` (the Pyodide proxy), 30 steps after a 5-step warm-up. Browser ≈ 3.5× native, per tranche 1's calibration.

| stack | native, no JIT | ≈ browser |
|---|---|---|
| page 07 gray (`single_band_gray_lw`): LW + slab | 2.4 ms | ~8 ms |
| page 07 non-grey (`earth_low_res_lw`): LW + slab | 48.5 ms | ~0.17 s |
| page 08 gray: + `SimpleBoundaryLayer` | 3.0 ms | ~11 ms |
| page 11 dry RCE: 14-band + SBL + dry adjustment | 52.8 ms | ~0.18 s |
| page 12 moist RCE: + Emanuel + condensation | 53.5 ms | ~0.19 s |

Radiation dominates by ~18×; every component this tranche adds is nearly free. **Cost is set almost entirely by the timestep and the step count**, which is what Task 8 measures.

With numba on (native, not the browser) the same five stacks run 2.5 / 2.5 / 2.9 / 4.5 / 6.1 ms — the 14-band radiation is ~20× cheaper JIT-compiled, which is why the `NUMBA_DISABLE_JIT=1` proxy is mandatory and a native timing is worthless here.

## File structure

| File | Responsibility |
|---|---|
| `climt/_components/simple_boundary_layer/component.py` | *modify* — coerce the timestep to a plain float (Task 1) |
| `climt/_components/cork/lw/component.py` | *modify* — raise on a non-finite longwave result instead of returning NaN (Task 2) |
| `docs/modelling-tour/_tour/assets.py` | *create* — browser/native asset path resolution, lifted out of `tables.py` |
| `docs/modelling-tour/_tour/tables.py` | *modify* — becomes a caller of `assets.py` |
| `docs/modelling-tour/_tour/stepping.py` | *create* — `integrate`, `integrate_with_snapshots`, the evolution figure |
| `docs/modelling-tour/_tour/budgets.py` | *create* — TOA and surface energy budget diagnostics |
| `docs/modelling-tour/_tour/states.py` | *create* — save/load a sympl state as `.npz`, with provenance |
| `docs/modelling-tour/_data/rce_dry_equilibrium.npz` | *create* — page 11's shipped equilibrium |
| `docs/modelling-tour/_data/rce_moist_equilibrium.npz` | *create* — page 12's shipped equilibrium |
| `scripts/generate_tour_equilibria.py` | *create* — regenerates both, defaults matching what shipped |
| `scripts/experiments/tour_rce_measurements.py` | *create* — the Task 8 measurement campaign |
| `scripts/serve_wheel.py` | *moved* from `docs/radiative-transfer/_live/` |
| `docs/modelling-tour/07-stepping-a-column.qmd` … `12-moist-rce.qmd` | the six pages |
| `docs/modelling-tour/index.qmd` | *modify* — chapter map gains six rows; the "nothing integrates in time" section is rewritten |
| `docs/_includes/climt-live-boot.qmd` | *modify* — absorbs the unreleased-wheel note; loses the reference to the deleted setup include |
| `docs/_includes/climt-live-setup.qmd` | *delete* |
| `docs/radiative-transfer/09-live-rce.qmd` | *delete* |
| `docs/radiative-transfer/_live/rce_helpers.py` | *delete* |
| `docs/_quarto.yml` | *modify* — sidebar, the `live-html` comment |
| `tests/test_modelling_tour.py` | the tranche 2 helpers and each page's physics claim |
| `tests/test_live_rce_demo.py` | *delete* — folded into the above |
| `tests/test_simple_boundary_layer.py` | *modify* — a no-numba regression test |

Tasks 1–2 are library fixes. Tasks 3–7 build and retire infrastructure. Task 8 is the spec's Task 0 — the measurement campaign every page number descends from. Task 9 ships the equilibrium states. Tasks 10–15 are the pages, each shippable alone. Task 16 closes out.

---

## Task 1: `SimpleBoundaryLayer` is broken in the browser — coerce the timestep

**This task is not in the spec, and it blocks four of the six pages.** It was found while measuring for this plan.

`SimpleBoundaryLayer.array_call` passes `timestep.total_seconds()` straight into `_boundary_layer_kernel`. Under `UnytBackend` that is a `unyt_quantity` carrying units of seconds, not a float. With numba compiling the kernel the units are stripped on the way in and nothing goes wrong, which is why every existing test passes. With `NUMBA_DISABLE_JIT=1` — that is, **in Pyodide, where there is no numba at all** — the quantity survives into `_diffuse_profile`, and:

```
File "climt/_components/simple_boundary_layer/component.py", line 82, in _diffuse_profile
    diag[0] = diag[0] + surface_exchange
unyt.exceptions.UnitOperationError: The <ufunc 'add'> operator for unyt_arrays with
units 'dimensionless' (dimensions '1') and 's' (dimensions '(time)') is not well defined.
```

Pages 08, 09, 11 and 12 all use this component. Without this fix they raise on their first cell in every reader's browser while passing every test on the developer's machine.

The precedent for the fix is already in the codebase: `climt/_components/sea_ice/component.py:333` does `dt = float(timestep.total_seconds())` with a comment explaining exactly this. `land_ice`, `second_best` and `bucket_hydrology` carry the same note. This component was missed.

**Files:**
- Modify: `climt/_components/simple_boundary_layer/component.py:448`
- Test: `tests/test_simple_boundary_layer.py` (modify)

**Interfaces:**
- Consumes: nothing.
- Produces: `SimpleBoundaryLayer()(state, timestep)` works under `climt.UnytBackend` with `NUMBA_DISABLE_JIT=1`. No signature changes.

- [ ] **Step 1: Create the branch**

```bash
git checkout develop && git pull
git checkout -b feature/modelling-tour-rce
```

- [ ] **Step 2: Write the failing test**

Append to `tests/test_simple_boundary_layer.py`. Two tests: one that fails on the
unfixed source with the JIT on (which is how CI runs), and one that exercises the
component end to end under the backend the pages use.

```python
# ------------------------------------------------- the no-numba / browser path

def test_boundary_layer_kernel_is_handed_a_plain_float_timestep():
    """Guards the Pyodide path: the kernel must never receive a unyt quantity.

    ``timestep.total_seconds()`` is a plain float for ``datetime.timedelta``
    but a ``unyt_quantity`` in seconds for ``climt.UnytTimeDelta``, which is
    the timestep type the UnytBackend requires. numba strips those units on
    the way into the kernel; the pure-Python path does not, and the units then
    collide inside ``_diffuse_profile``, where a dimensionless tridiagonal
    diagonal is added to a seconds-carrying exchange coefficient.

    ``NUMBA_DISABLE_JIT=1`` reproduces the real failure, but CI runs with the
    JIT on, where the bug is invisible. Asserting on the source is the only
    way to keep it from coming back in routine CI.
    """
    import inspect
    from climt._components.simple_boundary_layer import component as sbl

    source = inspect.getsource(sbl.SimpleBoundaryLayer.array_call)
    assert "float(timestep.total_seconds())" in source, (
        "array_call must coerce the timestep to a plain float before calling "
        "_boundary_layer_kernel: UnytTimeDelta.total_seconds() returns a "
        "unyt_quantity in seconds, and without numba to strip the units it "
        "collides with the dimensionless tridiagonal diagonal in "
        "_diffuse_profile. See sea_ice/component.py for the same coercion.")


def test_bulk_fluxes_run_under_the_unyt_backend():
    """The component runs with the timestep type the modelling-tour pages use.

    Restores the default backend afterwards: ``sympl.set_backend`` is global,
    and the rest of this file builds states expecting DataArrays.
    """
    import sympl

    sympl.set_backend(climt.UnytBackend())
    try:
        component = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
        state = climt.get_default_state(
            [component], grid_state=climt.get_grid(nx=1, ny=1, nz=28))
        diagnostics, new_state = component(state, climt.UnytTimeDelta(hours=1))
    finally:
        sympl.set_backend(sympl.DataArrayBackend())

    assert np.all(np.isfinite(new_state["air_temperature"].values))
    assert np.all(np.isfinite(
        diagnostics["surface_upward_sensible_heat_flux"].values))
```

- [ ] **Step 3: Run the tests to verify the first one fails**

Run: `conda run -n climt python -m pytest tests/test_simple_boundary_layer.py -k "plain_float_timestep or unyt_backend" -v`
Expected: `test_boundary_layer_kernel_is_handed_a_plain_float_timestep` FAILS with the assertion message above; `test_bulk_fluxes_run_under_the_unyt_backend` PASSES (the JIT is hiding the bug, which is the point).

- [ ] **Step 4: Reproduce the real failure**

```bash
NUMBA_DISABLE_JIT=1 conda run -n climt python -c "
import sympl, climt
sympl.set_backend(climt.UnytBackend())
bl = climt.SimpleBoundaryLayer()
s = climt.get_default_state([bl], grid_state=climt.get_grid(nx=1, ny=1, nz=28))
bl(s, climt.UnytTimeDelta(hours=1))"
```

Expected: `unyt.exceptions.UnitOperationError ... 'dimensionless' ... and 's'`.
This is the failure a reader would see in their browser.

- [ ] **Step 5: Apply the fix**

In `climt/_components/simple_boundary_layer/component.py`, in `array_call`'s call
to `_boundary_layer_kernel` (line 448), change:

```python
            timestep.total_seconds(),
```

to:

```python
            # A plain float, not the unyt_quantity UnytTimeDelta.total_seconds()
            # returns: numba strips the units on the way in, but Pyodide has no
            # numba, and the pure-Python kernel then adds seconds to the
            # dimensionless tridiagonal diagonal in _diffuse_profile. Same
            # coercion, same reason, as sea_ice/component.py.
            float(timestep.total_seconds()),
```

- [ ] **Step 6: Verify both the test and the real path**

```bash
conda run -n climt python -m pytest tests/test_simple_boundary_layer.py -v
NUMBA_DISABLE_JIT=1 conda run -n climt python -c "
import sympl, climt, numpy as np
sympl.set_backend(climt.UnytBackend())
bl = climt.SimpleBoundaryLayer()
s = climt.get_default_state([bl], grid_state=climt.get_grid(nx=1, ny=1, nz=28))
d, n = bl(s, climt.UnytTimeDelta(hours=1))
print('OK', float(d['surface_upward_sensible_heat_flux'].values.ravel()[0]))"
```

Expected: the whole file passes, and the second command prints `OK` and a finite flux.

- [ ] **Step 7: Check the sibling components for the same bug**

```bash
conda run -n climt python -m pytest tests/ -k "boundary_layer or slab or dry_conv or condensation" -q
NUMBA_DISABLE_JIT=1 conda run -n climt python -c "
import sympl, climt
sympl.set_backend(climt.UnytBackend())
dt = climt.UnytTimeDelta(hours=1)
for c in (climt.GridScaleCondensation(), climt.DryConvectiveAdjustment(),
          climt.EmanuelConvectionPython()):
    s = climt.get_default_state([c], grid_state=climt.get_grid(nx=1, ny=1, nz=28))
    s['specific_humidity'].values[:] = 5e-3
    c(s, dt); print(type(c).__name__, 'OK')"
```

Expected: all three print `OK`. They were checked while writing this plan and
were clean; this re-checks after the edit rather than trusting the note.

- [ ] **Step 8: Note the artifact gate**

This task edits `climt/_components/`, but **not** under `climt/_components/cork/`,
so the content-hashed docs experiment artifacts are untouched. Confirm:

```bash
conda run -n climt python scripts/build_experiments.py --check
```

Expected: clean. If it reports a mismatch, stop — something else in the tree
is dirty and needs sorting out before this branch grows.

- [ ] **Step 9: Commit**

```bash
git add climt/_components/simple_boundary_layer/component.py tests/test_simple_boundary_layer.py
git commit -m "fix(bl): SimpleBoundaryLayer works without numba, so it works in the browser

array_call forwarded UnytTimeDelta.total_seconds() -- a unyt_quantity in
seconds -- into the kernel. numba stripped the units, so every test passed;
Pyodide has no numba, and _diffuse_profile then added seconds to a
dimensionless tridiagonal diagonal and raised UnitOperationError. Coerced to
float, matching sea_ice/component.py.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 2: `CorkLongwaveRadiation` raises on a non-finite result instead of returning NaN

Six pages invite the experiment that triggers this: "what happens if I make the timestep bigger?" Past the stability limit the column temperature goes non-finite, `net_band = up_band[b] - down_band[b]` emits `RuntimeWarning: invalid value encountered in subtract`, and the component returns NaN tendencies. In the browser that surfaces as a blank figure — no error, no message, nothing to search for.

The component should detect it and raise, naming the timestep as the likely cause.

**Files:**
- Modify: `climt/_components/cork/lw/component.py` (in `array_call`, after `lw_transport` returns)
- Test: `tests/test_cork_lw.py` (modify)

**Interfaces:**
- Consumes: nothing.
- Produces: `CorkLongwaveRadiation(...)(state)` raises `ValueError` when the computed longwave fluxes are not all finite. Message names the offending quantity and the timestep. No signature changes; no new constructor argument (a silent-NaN escape hatch would defeat the point).

- [ ] **Step 1: Write the failing test**

Append to `tests/test_cork_lw.py`:

```python
def test_nonfinite_input_raises_instead_of_returning_nan():
    """A blown-up column must produce a message, not a blank figure.

    The modelling-tour pages invite readers to raise the timestep until the
    integration goes unstable. When it does, the temperature goes non-finite,
    the two-stream solve propagates that into the fluxes, and the component
    used to return NaN tendencies with only a numpy RuntimeWarning. Under
    Pyodide a warning is invisible and the figure just comes out empty.
    """
    import sympl

    sympl.set_backend(climt.UnytBackend())
    try:
        longwave = climt.CorkLongwaveRadiation(
            optics="correlated_k", table="earth_low_res_lw")
        state = climt.get_default_state(
            [longwave], grid_state=climt.get_grid(nx=1, ny=1, nz=28))
        state["air_temperature"].values[10, 0, 0] = np.nan

        with pytest.raises(ValueError, match="non-finite"):
            longwave(state)
    finally:
        sympl.set_backend(sympl.DataArrayBackend())
```

- [ ] **Step 2: Run it to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_cork_lw.py::test_nonfinite_input_raises_instead_of_returning_nan -v`
Expected: FAIL — `DID NOT RAISE <class 'ValueError'>` (the call returns NaN tendencies, possibly with a `RuntimeWarning`).

- [ ] **Step 3: Add the guard**

In `climt/_components/cork/lw/component.py`, in `array_call`, immediately after
the `lw_result` unpacking (the `if self._diagnostics_level > 0:` / `else:` pair
that assigns `up_band, down_band, up_broad, down_broad`) and **before**
`net_flux = up_broad - down_broad`, insert:

```python
        # A non-finite result here is almost always an unstable time
        # integration upstream: the temperature profile blew up, and the
        # two-stream solve carried it into the fluxes. Returning NaN
        # tendencies makes that a blank figure in the browser and a silent
        # wrong answer everywhere else, so name it instead. The check is two
        # array reductions on already-computed arrays -- negligible against
        # the transport solve it follows.
        if not (np.all(np.isfinite(up_broad)) and np.all(np.isfinite(down_broad))):
            culprit = ("air_temperature"
                       if not np.all(np.isfinite(T_flat))
                       else "surface_temperature"
                       if not np.all(np.isfinite(T_surf_flat))
                       else "the longwave transport solve")
            raise ValueError(
                "CorkLongwaveRadiation produced non-finite longwave fluxes; "
                f"{culprit} is non-finite on input or became so. The usual "
                "cause is a time step past this configuration's stability "
                "limit -- halve it and re-run. (If you are calling the "
                "component directly on a prescribed profile, check that "
                "profile for NaN or inf.)")
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `conda run -n climt python -m pytest tests/test_cork_lw.py::test_nonfinite_input_raises_instead_of_returning_nan -v`
Expected: PASS.

- [ ] **Step 5: Verify the guard fires on the real failure mode**

The message claims a too-large timestep is the usual cause. Check that it is:

```bash
conda run -n climt python -c "
import sympl, climt
from sympl import AdamsBashforth
sympl.set_backend(climt.UnytBackend())
lw = climt.CorkLongwaveRadiation(optics='correlated_k', table='single_band_gray_lw')
surf = climt.SlabSurface()
s = climt.get_default_state([lw, surf], grid_state=climt.get_grid(nx=1, ny=1, nz=28))
s['ocean_mixed_layer_thickness'].values[:] = 2.0
s['downwelling_shortwave_flux_in_air'].values[:] = 0.0
s['downwelling_shortwave_flux_in_air'].values[0, ...] = 240.0
s['upwelling_shortwave_flux_in_air'].values[:] = 0.0
m = AdamsBashforth([lw, surf])
dt = climt.UnytTimeDelta(hours=24)
try:
    for i in range(200):
        d, n = m(s, dt); s.update(n); s.update(d); s['time'] += dt
    print('no blow-up in 200 steps at dt=24h')
except ValueError as e:
    print('raised at step', i, '->', str(e)[:90])"
```

Expected: it raises, and the message begins `CorkLongwaveRadiation produced non-finite longwave fluxes`. Record the step index it reached in the task log — page 07's "make the timestep bigger" exercise quotes it. If it does *not* blow up at 24 h with a 2 m slab, record that instead and adjust the exercise on page 07 accordingly; the spec's "NaN at 24 h" was measured on a 5 m slab at `nz = 18`.

- [ ] **Step 6: Run the full cork suite — this edits `cork/`**

```bash
conda run -n climt python -m pytest tests/test_cork_lw.py tests/test_cork_integration.py tests/test_cork_validation.py -q
```

Expected: all pass. A test that deliberately feeds a NaN profile and expects NaN back would now fail; there is none, but check the output rather than assuming.

- [ ] **Step 7: Regenerate the content-hashed docs artifacts**

Editing anything under `climt/_components/cork/` invalidates the dependency hashes that gate `docs/experiments/`, and leaves the docs workflow red. This is the bill tranche 1 documented (its plan, "a one-line no-op in `cork/` costs two 200-day integrations"), and it is unavoidable here.

```bash
conda run -n climt python scripts/build_experiments.py --check
```

Expected: reports stale artifacts. Regenerate exactly those:

```bash
make experiments ONLY='<the globs --check named>'
conda run -n climt python scripts/build_experiments.py --check
```

Expected on the second `--check`: clean. **The regenerated artifact *values* must not change** — this task adds a raise on a path that never fires in a healthy run. Inspect the diff: if any number moved, the guard is firing somewhere it should not, and that is a bug in Step 3, not a new result to accept.

- [ ] **Step 8: Commit**

```bash
git add climt/_components/cork/lw/component.py tests/test_cork_lw.py docs/experiments
git commit -m "fix(cork): a blown-up column says so, instead of drawing nothing

Past the stability limit the longwave solve returned NaN with only a numpy
RuntimeWarning, which in the browser is a blank figure and no error. Raise
ValueError naming the timestep as the likely cause. Docs artifacts
regenerated for the dependency hash; no values changed.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 3: `_tour/assets.py` — one place that finds a site asset

`_tour/tables.py` already owns browser/native path resolution for the 56-band spectrum table. This tranche adds two more assets (the equilibrium `.npz` states), and duplicating that logic three times is how it drifts. Lift it into `assets.py`; `tables.py` becomes a caller. **This is the only refactor this tranche makes to tranche 1 code**, and it is in service of the current goal, not tidiness.

The mechanism, unchanged: `docs/_quarto.yml` lists `modelling-tour/_data/*.npz` under `project: resources:` so Quarto publishes underscore-prefixed paths, and each page lists what it needs under `pyodide: resources:` so quarto-live stages it into the Pyodide filesystem at the same relative path it has on disk. So the lookup is a plain `os.path.isfile` in both environments.

**Files:**
- Create: `docs/modelling-tour/_tour/assets.py`
- Modify: `docs/modelling-tour/_tour/tables.py`
- Test: `tests/test_modelling_tour.py` (modify)

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `assets.resolve(name, base_url="_data") -> str | None` — filesystem path to a staged asset, or `None` if it is not there.
  - `assets.fetch_same_origin(name, base_url="_data") -> str | None` — last-resort synchronous browser fetch into the Pyodide filesystem; returns the local filename, or `None` outside Pyodide or on failure.
  - `assets.locate(name, base_url="_data") -> str | None` — `resolve` then `fetch_same_origin`.
  - `tables.spectrum_table(prefer_hires=True, base_url="_data")` — unchanged signature and unchanged behaviour, now implemented over `assets.locate`.
  - Modules under `_tour/` may `import assets` by bare name. Pages already do `sys.path.insert(0, "_tour")`; the test loader gains the same.

- [ ] **Step 1: Write the failing test**

In `tests/test_modelling_tour.py`, first make cross-module imports work the way
they do in the browser. Replace `_load` with:

```python
# The pages do `sys.path.insert(0, "_tour")` and then import helpers by bare
# name; _tour modules import each other the same way. Match that here so a
# module loaded by path can still `import assets`.
if str(TOUR) not in sys.path:
    sys.path.insert(0, str(TOUR))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, TOUR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
```

and add `import sys` to the imports at the top of the file.

Then add:

```python
@pytest.fixture
def assets():
    return _load("assets")


def test_assets_resolve_finds_a_staged_file_from_any_working_directory(
        assets, monkeypatch, tmp_path):
    """Resolution must not depend on the caller's working directory.

    In the browser the page's working directory *is* the page directory, so
    `_data/x.npz` resolves directly. Natively — a test run, a static-figure
    render, someone poking at it from the repo root — it does not, and the
    fallback is the path relative to this module's own location.
    """
    monkeypatch.chdir(tmp_path)
    found = assets.resolve("earth_spectrum_lw.npz")
    assert found is not None and os.path.isfile(found)


def test_assets_resolve_returns_none_for_a_missing_asset(assets, tmp_path,
                                                        monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert assets.resolve("no_such_asset.npz") is None


def test_assets_locate_falls_back_to_resolve_outside_pyodide(assets, tmp_path,
                                                             monkeypatch):
    """`locate` is `resolve` plus a browser-only fetch. Natively they agree."""
    monkeypatch.chdir(tmp_path)
    assert assets.locate("earth_spectrum_lw.npz") == \
        assets.resolve("earth_spectrum_lw.npz")
    assert assets.locate("no_such_asset.npz") is None


def test_assets_honours_a_non_default_base_url(assets, tmp_path, monkeypatch):
    """Page 11 and 12 pass the same `_data`, but the argument must work."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "elsewhere").mkdir()
    (tmp_path / "elsewhere" / "probe.npz").write_bytes(b"x")
    assert assets.resolve("probe.npz", base_url="elsewhere") == \
        os.path.join("elsewhere", "probe.npz")
```

Add `import os` at the top of the test file if it is not already there.

- [ ] **Step 2: Run it to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k assets -v`
Expected: FAIL — `FileNotFoundError` / `ModuleNotFoundError` for `_tour/assets.py`, which does not exist yet.

- [ ] **Step 3: Write `_tour/assets.py`**

```python
"""Locate a Modelling Tour site asset, in the browser and natively.

Some of what these pages need is too large or too configuration-specific to
ship inside the climt wheel: the 56-band spectrum table (~5 MB) and the two
radiative-convective equilibrium states. Those live beside the pages as static
site assets under ``_data/``.

``docs/_quarto.yml`` lists ``modelling-tour/_data/*.npz`` under
``project: resources:`` so Quarto publishes them at all (it skips
underscore-prefixed paths otherwise), and each page lists the ones it uses
under ``pyodide: resources:`` so quarto-live stages them into the Pyodide
filesystem at document setup, at the same relative path they have on disk.

That is why the lookup below is a plain ``os.path.isfile`` in both
environments: the browser case needs no network code of its own and no CORS
negotiation, because the docs site is same-origin with the page. (A GitHub
release asset would not work -- it sends no CORS headers.)

This module is the single place that knows any of the above.
``tables.py`` and ``states.py`` are its callers.
"""
import os

DEFAULT_BASE_URL = "_data"


def _candidates(name, base_url):
    """Paths to try, in order, for a staged asset."""
    yield os.path.join(base_url, name)
    # `_tour/` is one level below the page, so the asset sits next to our own
    # parent directory. This is what makes resolution work from a native
    # render or a test run, where the working directory is not the page's.
    here = os.path.dirname(os.path.abspath(__file__))
    yield os.path.join(here, os.pardir, base_url, name)


def resolve(name, base_url=DEFAULT_BASE_URL):
    """Return a filesystem path to ``name``, or None if it is not staged.

    Args:
        name: the asset's file name, e.g. ``"rce_dry_equilibrium.npz"``.
        base_url: directory holding it, relative to the page.

    Returns:
        A path that ``open()`` will accept, or ``None``.
    """
    for path in _candidates(name, base_url):
        if os.path.isfile(path):
            return path
    return None


def fetch_same_origin(name, base_url=DEFAULT_BASE_URL):
    """Browser-only last resort: fetch a same-origin asset into the FS.

    Used when a page forgot to declare the asset under ``pyodide: resources:``.
    Synchronous by necessity -- quarto-live does not await a cell's async work
    before starting the next one, so an ``await``-based fetch would race.

    Returns the local file name on success, or ``None`` outside Pyodide or if
    the fetch fails for any reason.
    """
    try:
        from pyodide.http import open_url

        data = open_url("{}/{}".format(base_url, name))
        with open(name, "wb") as handle:
            handle.write(data.getvalue().encode("latin-1"))
        return name
    except Exception:
        return None


def locate(name, base_url=DEFAULT_BASE_URL):
    """``resolve``, then the browser fetch. ``None`` if neither works."""
    return resolve(name, base_url) or fetch_same_origin(name, base_url)
```

- [ ] **Step 4: Rewrite `_tour/tables.py` over it**

Replace the whole file with:

```python
"""Choose the longwave absorption table a Modelling Tour page runs on.

The high-resolution 56-band table is a site asset rather than wheel data --
``_tour/assets.py`` explains the mechanism and owns the path resolution. This
module is the policy on top of it: which table a page gets, and what happens
when the asset is absent.

Pages call :func:`spectrum_table` and pass the result straight to
``CorkLongwaveRadiation(table=...)``. If the asset cannot be found, the shipped
14-band table is returned instead, so a page degrades to a coarser spectrum
rather than failing -- which is also what happens on the pages that never ask
for it.
"""
import assets

FALLBACK = "earth_low_res_lw"
ASSET = "earth_spectrum_lw.npz"


def spectrum_table(prefer_hires=True, base_url=assets.DEFAULT_BASE_URL):
    """Return a table name or path for ``CorkLongwaveRadiation``.

    Args:
        prefer_hires: if False, always return the shipped 14-band table.
        base_url: directory (or URL prefix) holding the asset, relative to the
            page.

    Returns:
        A table name or filesystem path, always usable as ``table=``.
    """
    if not prefer_hires:
        return FALLBACK
    return assets.locate(ASSET, base_url) or FALLBACK
```

- [ ] **Step 5: Run the whole tour test file**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -m "not slow" -v`
Expected: the four new `assets` tests pass, and **every pre-existing `spectrum_table` test still passes unchanged** — `test_spectrum_table_finds_the_asset_from_any_working_directory`, `..._honours_prefer_hires`, `..._falls_back_when_the_asset_is_missing`, `..._result_is_usable_as_a_table_argument`. Those four are the contract; if any needed editing to pass, the refactor changed behaviour and is wrong.

- [ ] **Step 6: Verify the browser import path still works**

Pages 1–3 declare `_tour/tables.py` under `pyodide: resources:` but not `assets.py`. `tables.py` now imports it, so every page that lists `tables.py` must also list `assets.py`. Add the line to the `resources:` block of each:

```bash
grep -ln "_tour/tables.py" docs/modelling-tour/*.qmd
```

Expected: `01-emissivity-spectrum.qmd`, `02-window-measured.qmd`, `03-radiating-level.qmd`. In each, directly above `    - _tour/tables.py`, insert:

```yaml
    # tables.py imports assets.py, which owns the path resolution both it and
    # the equilibrium-state loader need. Staged the same way, or the import
    # fails in the browser while passing every native test.
    - _tour/assets.py
```

- [ ] **Step 7: Confirm the pages still run natively**

Run: `conda run -n climt python docs/modelling-tour/_artifacts/generate.py 01 02 03`
Expected: three PNGs written, no traceback. This execs the pages' own cells, so it is the cheapest check that the import rearrangement did not break them.

`git diff --stat docs/modelling-tour/_artifacts/` should show the three PNGs either unchanged or differing only in encoding noise. If a figure's *content* changed, stop: the refactor altered which table a page resolved.

- [ ] **Step 8: Commit**

```bash
git add docs/modelling-tour/_tour/assets.py docs/modelling-tour/_tour/tables.py \
        docs/modelling-tour/0{1,2,3}-*.qmd tests/test_modelling_tour.py
git commit -m "refactor(tour): one place that finds a site asset

tables.py owned browser/native path resolution; this tranche adds two more
assets that need the same logic. Lifted into _tour/assets.py, with tables.py
as a caller and its four existing tests unchanged as the contract.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 4: `_tour/stepping.py` — the time loops, the forcings that keep them honest, and the evolution figure

The live-RCE page's `integrate_to_equilibrium` and `integrate_with_snapshots` move here, corrected. They are the same loop the deleted page ran, with four changes: `UnytTimeDelta` instead of `datetime.timedelta`, stepper components allowed in the snapshot loop (page 08 onward has them), recording split from drawing so the recording can be tested without matplotlib, and no duplicated copy anywhere — the keep-two-copies-in-sync contract the old docstring described dies with this task.

**Files:**
- Create: `docs/modelling-tour/_tour/stepping.py`
- Test: `tests/test_modelling_tour.py` (modify)

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `stepping.integrate(tendency_components, stepper_components, state, timestep, n_steps) -> state`
  - `stepping.integrate_with_history(tendency_components, stepper_components, state, timestep, n_steps, n_snapshots=8) -> (state, history)`
  - `history` is a plain dict: `{"days": (n,) float array, "olr": (n,), "surface_temperature": (n,), "snapshots": [ {"day", "T", "H", "U", "D"} , ... ]}`
  - `stepping.draw_evolution(history, state, solar, title="") -> None` (calls `plt.show()`)
  - `stepping.DAY_SECONDS = 86400.0`
  - `stepping.UnytRelaxation(quantity_name, units, tendency_units)` — a `sympl.RelaxationTendencyComponent` subclass whose tendency units `unyt` can parse (Step 8)
  - `stepping.wind_relaxation(state, speed, timescale_hours=24.0, quantity_name="eastward_wind", initialise=True) -> UnytRelaxation` — builds the component **and** writes the two target fields it needs into `state`. `initialise=False` when attaching to a state loaded from disk, which already carries a spun-up wind profile

- [ ] **Step 1: Write the failing test**

Add to `tests/test_modelling_tour.py`:

```python
@pytest.fixture
def stepping():
    return _load("stepping")


def _gray_column(nz=28, slab_depth=2.0, table="single_band_gray_lw"):
    """Page 07's column: gray longwave over a slab, prescribed absorbed solar."""
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k", table=table)
    surface = climt.SlabSurface()
    state = climt.get_default_state([longwave, surface],
                                    grid_state=get_grid(nx=1, ny=1, nz=nz))
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    return [longwave, surface], state


def test_integrate_advances_the_clock_and_the_state(stepping):
    components, state = _gray_column()
    start = state["time"]
    T0 = state["air_temperature"].values.copy()

    timestep = climt.UnytTimeDelta(hours=12)
    out = stepping.integrate(components, [], state, timestep, 10)

    assert out is state                       # stepped in place, and returned
    assert (out["time"] - start).total_seconds() == pytest.approx(10 * 12 * 3600)
    assert not np.allclose(out["air_temperature"].values, T0)


def test_integrate_keeps_the_flux_diagnostics(stepping):
    """The update-order gotcha that broke the old demo once already.

    The prognostic state is applied *before* the diagnostics. Reverse it and
    the stepper carries the longwave fluxes forward at their pre-step value,
    SlabSurface consumes stale fluxes, and the surface heats without bound.
    """
    components, state = _gray_column()
    stepping.integrate(components, [], state, climt.UnytTimeDelta(hours=12), 20)

    for flux in ("upwelling_longwave_flux_in_air",
                 "downwelling_longwave_flux_in_air"):
        values = state[flux].values
        assert np.all(np.isfinite(values)) and np.all(values >= 0.0)
    surface = float(state["surface_temperature"].values.ravel()[0])
    assert 150.0 < surface < 400.0


def test_integrate_runs_stepper_components_too(stepping):
    """Page 08 onward puts a Stepper in the same loop as the tendencies."""
    components, state = _gray_column()
    boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="bulk")
    state = climt.get_default_state(components + [boundary_layer],
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    state["ocean_mixed_layer_thickness"].values[:] = 2.0
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0

    stepping.integrate(components, [boundary_layer], state,
                       climt.UnytTimeDelta(hours=1), 10)

    # The stepper's own diagnostics survive into the state.
    assert "surface_upward_sensible_heat_flux" in state
    assert np.all(np.isfinite(
        state["surface_upward_sensible_heat_flux"].values))
    assert np.all(np.isfinite(state["air_temperature"].values))


def test_integrate_rejects_a_plain_timedelta(stepping):
    """The UnytTimeDelta gotcha, asserted rather than left as prose.

    Page 07 teaches this; if sympl ever stops raising, the page's craft note is
    wrong and this test is where we find out.
    """
    from datetime import timedelta

    components, state = _gray_column()
    with pytest.raises(Exception) as caught:
        stepping.integrate(components, [], state, timedelta(hours=12), 1)
    assert "degK" in str(caught.value) or "unit" in str(caught.value).lower()


def test_history_records_every_step_and_n_snapshots(stepping):
    components, state = _gray_column()
    state, history = stepping.integrate_with_history(
        components, [], state, climt.UnytTimeDelta(hours=12), 12,
        n_snapshots=4)

    assert history["days"].shape == (12,)
    assert history["olr"].shape == (12,)
    assert history["surface_temperature"].shape == (12,)
    assert history["days"][0] == pytest.approx(0.5)     # 12 h, in days
    assert history["days"][-1] == pytest.approx(6.0)
    assert len(history["snapshots"]) == 4
    first, last = history["snapshots"][0], history["snapshots"][-1]
    assert first["day"] < last["day"]
    for key in ("T", "H", "U", "D"):
        assert np.all(np.isfinite(first[key]))
    assert first["T"].shape == (28,)          # mid levels
    assert first["U"].shape == (29,)          # interface levels


def test_history_olr_matches_the_final_state(stepping):
    """The recorded series is the state's own diagnostic, not a recomputation."""
    components, state = _gray_column()
    state, history = stepping.integrate_with_history(
        components, [], state, climt.UnytTimeDelta(hours=12), 8)

    assert history["olr"][-1] == pytest.approx(
        float(state["upwelling_longwave_flux_in_air"].values[-1, 0, 0]))
    assert history["surface_temperature"][-1] == pytest.approx(
        float(state["surface_temperature"].values.ravel()[0]))
```

Add near the top of the test file, beside the other module-level constants:

```python
from climt import get_grid

SOLAR = 240.0     # prescribed absorbed shortwave at the surface, W m^-2, on
                  # every tranche 2 page. There is no shortwave component.
```

- [ ] **Step 2: Run it to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k "stepping or integrate or history" -v`
Expected: FAIL — `_tour/stepping.py` does not exist.

- [ ] **Step 3: Write `_tour/stepping.py`**

```python
"""Step a column forward in time, and draw what happened.

Every page from 07 on runs one of these two loops. They are the single copy:
the deleted ``docs/radiative-transfer/09-live-rce.qmd`` kept one version in a
Quarto include and a testable twin in a ``.py`` beside it, with a docstring
contract asking the two to be edited together and a test enforcing it. Pages
here import this module through ``pyodide: resources:`` instead, so there is
nothing to keep in sync.

Two component categories, and only two:

* **tendency components** (``CorkLongwaveRadiation``, ``SlabSurface``,
  ``EmanuelConvectionPython``) return rates of change, and are wrapped in a
  ``sympl.AdamsBashforth`` time integrator that turns them into a new state;
* **stepper components** (``SimpleBoundaryLayer``,
  ``DryConvectiveAdjustment``, ``GridScaleCondensation``) return the new state
  themselves, and are called directly.

``EmanuelConvectionPython`` is an ``ImplicitTendencyComponent``, a third sympl
protocol, and putting one inside a ``TendencyStepper`` makes sympl warn that
results "may lead to scientifically invalid results". It goes in the tendency
list anyway -- that is how climt has always run it, and page 12 explains the
warning rather than restructuring the loop around it. So: two categories.

**Update order is load-bearing.** The new prognostic state is applied *before*
the diagnostics. Reverse it and the integrator carries flux diagnostics forward
at their pre-step value, ``SlabSurface`` consumes stale longwave fluxes, and
the surface heats without bound. This broke the demo once already.

**Timestep type is load-bearing.** Pass ``climt.UnytTimeDelta``, not
``datetime.timedelta``. Under ``UnytBackend`` a bare ``timedelta`` makes
``total_seconds()`` return a unitless float, which will not cancel a tendency's
``/s``, and sympl raises ``UnitOperationError``.
"""
import numpy as np
from sympl import AdamsBashforth

DAY_SECONDS = 86400.0


def integrate(tendency_components, stepper_components, state, timestep,
              n_steps):
    """Step ``state`` forward ``n_steps`` times, in place.

    Args:
        tendency_components: components returning tendencies; wrapped in one
            ``AdamsBashforth``.
        stepper_components: components returning a new state; called in the
            order given, after the tendency update, within each step.
        state: a climt state dict. Modified in place.
        timestep: a ``climt.UnytTimeDelta``.
        n_steps: number of steps.

    Returns:
        The same ``state`` object, for convenience.
    """
    model = AdamsBashforth(list(tendency_components))
    for _ in range(n_steps):
        diagnostics, new_state = model(state, timestep)
        state.update(new_state)
        state.update(diagnostics)
        for stepper in stepper_components:
            stepper_diagnostics, stepper_state = stepper(state, timestep)
            state.update(stepper_state)
            state.update(stepper_diagnostics)
        state["time"] += timestep
    return state


def integrate_with_history(tendency_components, stepper_components, state,
                           timestep, n_steps, n_snapshots=8):
    """``integrate``, recording the evolution as it goes.

    Two scalars are recorded every step -- they cost nothing -- and the full
    profiles at ``n_snapshots`` evenly spaced steps, which do not. The first
    snapshot is taken after step 1 rather than before step 1, so the radiation
    diagnostics it reads are populated.

    Returns:
        ``(state, history)``. ``history`` is a plain dict of numpy arrays and
        a list of snapshot dicts -- no sympl objects, so it survives being
        passed around a page and pickles cleanly:

        * ``days`` -- (n_steps,) elapsed time in days
        * ``olr`` -- (n_steps,) upwelling longwave flux at the top interface
        * ``surface_temperature`` -- (n_steps,) K
        * ``snapshots`` -- list of dicts with ``day``, ``T`` (mid levels),
          ``H`` (longwave heating rate, K/day), ``U`` and ``D`` (up and
          downwelling longwave flux on interface levels)
    """
    model = AdamsBashforth(list(tendency_components))
    day = timestep.total_seconds() / DAY_SECONDS
    # total_seconds() returns a unyt quantity under UnytBackend; the history is
    # plain numbers, so strip the units here rather than in six page cells.
    day = float(day)

    snap_at = set(np.unique(
        np.linspace(1, n_steps, n_snapshots).round().astype(int)).tolist())
    snapshots = []
    days, olr, surface_temperature = [], [], []

    for step in range(1, n_steps + 1):
        diagnostics, new_state = model(state, timestep)
        state.update(new_state)
        state.update(diagnostics)
        for stepper in stepper_components:
            stepper_diagnostics, stepper_state = stepper(state, timestep)
            state.update(stepper_state)
            state.update(stepper_diagnostics)
        state["time"] += timestep

        days.append(step * day)
        olr.append(float(
            state["upwelling_longwave_flux_in_air"].values[-1, 0, 0]))
        surface_temperature.append(float(
            state["surface_temperature"].values.ravel()[0]))
        if step in snap_at:
            snapshots.append(dict(
                day=step * day,
                T=state["air_temperature"].values[:, 0, 0].copy(),
                H=state["air_temperature_tendency_from_longwave"]
                .values[:, 0, 0].copy(),
                U=state["upwelling_longwave_flux_in_air"]
                .values[:, 0, 0].copy(),
                D=state["downwelling_longwave_flux_in_air"]
                .values[:, 0, 0].copy(),
            ))

    history = dict(
        days=np.array(days),
        olr=np.array(olr),
        surface_temperature=np.array(surface_temperature),
        snapshots=snapshots,
    )
    return state, history
```

- [ ] **Step 4: Run the recording tests**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k "stepping or integrate or history" -v`
Expected: all pass. `draw_evolution` is not tested here — it is added in Step 5 and checked by rendering, not by assertion.

- [ ] **Step 5: Add `draw_evolution` to `_tour/stepping.py`**

Append. This is the deleted page's `_draw_evolution`, unchanged in its layout decisions — including the log-pressure axis, which tranche 1's experiment log established is not optional (on a linear axis the gray column's isothermal layer is 1.6% of the panel and invisible) — with the signature changed to take a `history` dict and the state it came from.

```python
def draw_evolution(history, state, solar, title=""):
    """Render one annotated 2x2 figure of the integration's evolution.

    quarto-live runs each Pyodide cell to completion in a Web Worker and paints
    only the cell's *final* figure, so a frame-by-frame animation is impossible
    here -- an earlier version tried and only ever showed the last frame.
    Instead the whole run happens quietly behind quarto-live's spinner and this
    draws the evolution at the end: temperature, longwave heating rate and
    up/downwelling flux profiles at several times, plus OLR and surface
    temperature as time series approaching energy balance.

    Args:
        history: the dict returned by :func:`integrate_with_history`.
        state: the state it was recorded from -- read only for the pressure
            coordinates.
        solar: prescribed absorbed shortwave at the surface, W m^-2, drawn as
            the equilibrium target on the OLR panel.
        title: figure suptitle.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    snaps = history["snapshots"]
    days = history["days"]
    olr = history["olr"]
    tsurf = history["surface_temperature"]
    p = state["air_pressure"].values[:, 0, 0] / 100.0
    p_int = state["air_pressure_on_interface_levels"].values[:, 0, 0] / 100.0

    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(vmin=snaps[0]["day"], vmax=snaps[-1]["day"])

    fig = plt.figure(figsize=(11.5, 8.4))
    # right=0.88 reserves a margin for the colorbar. Letting colorbar(ax=[...])
    # steal the space instead does not work: its length comes from the default
    # aspect, so the bar is centred at half height and the panels it did not
    # shrink -- the OLR panel and its twinned temperature axis -- run their
    # labels straight through it.
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.30, right=0.88)
    ax_T = fig.add_subplot(gs[0, 0])
    ax_H = fig.add_subplot(gs[0, 1], sharey=ax_T)
    ax_F = fig.add_subplot(gs[1, 0], sharey=ax_T)
    ax_O = fig.add_subplot(gs[1, 1])

    for k, snap in enumerate(snaps):
        colour = cmap(norm(snap["day"]))
        final = (k == len(snaps) - 1)
        width = 2.6 if final else 1.3
        order = 5 if final else 2
        ax_T.plot(snap["T"], p, color=colour, lw=width, zorder=order)
        ax_H.plot(snap["H"], p, color=colour, lw=width, zorder=order)
        ax_F.plot(snap["U"], p_int, color=colour, lw=width, zorder=order)
        ax_F.plot(snap["D"], p_int, color=colour, lw=width, ls="--",
                  zorder=order)

    ax_T.set_ylabel("pressure (hPa)")
    ax_F.set_ylabel("pressure (hPa)")
    ax_T.set_xlabel("temperature (K)")
    ax_T.set_title("Temperature profile")
    ax_H.set_xlabel("LW heating rate (K day$^{-1}$)")
    ax_H.set_title("Longwave heating rate")
    ax_F.set_xlabel("LW flux (W m$^{-2}$)")
    ax_F.set_title("LW flux: up (solid) / down (dashed)")
    # Log-pressure y axis, always. On a linear axis the entire upper atmosphere
    # is crushed into the top sliver of the panel: the gray column's isothermal
    # layer spans ~20 -> 4 hPa, which is 1.6% of a linear axis and ~19% of this
    # one. Limits run over the *interface* levels so the flux panel's topmost
    # point is not clipped; high -> low inverts, and the shared y applies that
    # to all three profile panels.
    ax_T.set_yscale("log")
    ax_T.set_ylim(p_int.max(), p_int.min())
    for ax in (ax_T, ax_H, ax_F):
        ax.grid(alpha=0.3, which="both")
        ax.yaxis.set_major_formatter(ScalarFormatter())   # 1000, not 10^3
    ax_H.axvline(0.0, color="0.4", lw=0.8)

    T_final = snaps[-1]["T"]
    ax_T.annotate("warm surface", xy=(T_final[0], p[0]),
                  xytext=(0.42, 0.14), textcoords="axes fraction", fontsize=8,
                  arrowprops=dict(arrowstyle="->", color="0.3"))
    # Only claim a skin temperature when the profile has one. The gray column
    # relaxes to a uniform top; the non-grey column keeps cooling to the model
    # lid, because its topmost layers absorb only in the strong bands, where
    # the upwelling flux was already emitted by cold levels below. Labelling
    # that "isothermal" would be a caption contradicting its own figure.
    isothermal_top = abs(T_final[-1] - T_final[-2]) < 2.0
    ax_T.annotate("isothermal top\n(skin temperature)" if isothermal_top
                  else "no isothermal top:\nstill cooling at the model lid",
                  xy=(T_final[-1], p[-1]),
                  xytext=(0.30, 0.72), textcoords="axes fraction", fontsize=8,
                  arrowprops=dict(arrowstyle="->", color="0.3"))
    ax_H.text(0.04, 0.06,
              "cooling (< 0) aloft\nrelaxes toward 0\nat equilibrium",
              transform=ax_H.transAxes, fontsize=8, va="bottom")
    ax_F.text(0.55, 0.9, "surface emits up;\nspace sends nothing down",
              transform=ax_F.transAxes, fontsize=8, va="top")

    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    cax = fig.add_axes([0.945, 0.30, 0.013, 0.42])
    fig.colorbar(mappable, cax=cax).set_label(
        "time (days)  —  early (dark) → late (yellow)")

    ax_O.plot(days, olr, color="#c92a2a", lw=1.9,
              label="OLR (top of atmosphere)")
    ax_O.axhline(solar, color="0.35", ls="--", lw=1.2,
                 label=f"absorbed solar = {solar:.0f} W m$^{{-2}}$")
    ax_O.set_xlabel("time (days)")
    ax_O.set_ylabel("longwave flux (W m$^{-2}$)", color="#c92a2a")
    ax_O.tick_params(axis="y", labelcolor="#c92a2a")
    ax_O.set_title("Approach to energy balance")
    ax_O.grid(alpha=0.3)

    ax_Ts = ax_O.twinx()
    ax_Ts.plot(days, tsurf, color="#1c7ed6", lw=1.6)
    ax_Ts.set_ylabel("surface temperature (K)", color="#1c7ed6")
    ax_Ts.tick_params(axis="y", labelcolor="#1c7ed6")

    gap = solar - olr[-1]
    message = (f"still {gap:+.0f} W m$^{{-2}}$ from balance\n"
               "(raise n_steps to close it)"
               if abs(gap) > 2 else
               "OLR ≈ absorbed solar:\nnear equilibrium")
    ax_O.annotate(message, xy=(days[-1], olr[-1]), xytext=(0.12, 0.30),
                  textcoords="axes fraction", fontsize=8,
                  arrowprops=dict(arrowstyle="->", color="0.3"))
    ax_O.legend(loc="lower right", fontsize=7)

    if title:
        fig.suptitle(title, fontsize=13, y=0.98)
    # Required in the browser: quarto-live paints the figure a cell shows, and
    # a cell that only creates one paints nothing.
    plt.show()
```

- [ ] **Step 6: Check the figure renders**

```bash
conda run -n climt python -c "
import sys, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, sympl, climt
sys.path.insert(0, 'docs/modelling-tour/_tour')
import stepping
sympl.set_backend(climt.UnytBackend())
lw = climt.CorkLongwaveRadiation(optics='correlated_k', table='single_band_gray_lw')
surf = climt.SlabSurface()
s = climt.get_default_state([lw, surf], grid_state=climt.get_grid(nx=1, ny=1, nz=28))
s['ocean_mixed_layer_thickness'].values[:] = 2.0
s['downwelling_shortwave_flux_in_air'].values[:] = 0.0
s['downwelling_shortwave_flux_in_air'].values[0, ...] = 240.0
s['upwelling_shortwave_flux_in_air'].values[:] = 0.0
s, h = stepping.integrate_with_history([lw, surf], [], s, climt.UnytTimeDelta(hours=12), 200)
stepping.draw_evolution(h, s, 240.0, title='smoke test')
plt.savefig('/tmp/stepping_smoke.png', dpi=100, bbox_inches='tight')
print('wrote /tmp/stepping_smoke.png; Tsurf', float(s['surface_temperature'].values.ravel()[0]))"
```

Expected: a PNG is written with four populated panels and no matplotlib warnings. Open it and check the log-pressure axis reads `1000`, not `10^3`.

- [ ] **Step 8: Add wind relaxation — the momentum source a column does not have**

**Why this is here and not on page 08 alone.** A single column has no dynamics driving a wind, and `SimpleBoundaryLayer`'s own surface drag removes any wind you prescribe within days. Measured: three initialisations spanning 0–10 m s⁻¹ all converge to a lowest-level wind of 0.000 m s⁻¹ and the same equilibrium within 0.05 K. A column that has spun itself down is not a column with realistically weak winds — it is a column with *no* winds, and its surface fluxes are correspondingly wrong.

climt's own examples already handle this by brute force: `examples/column_code_with_slab.py:120` and `examples/gmd_radiative_convective.py:118` both end every timestep with `state['eastward_wind'].values[:] = 3.` A Newtonian relaxation is the cleaner form of the same idea — it is a *component*, so it composes, it carries units, it appears in the component list where a reader can see it, and its timescale is a physical parameter rather than a hidden assignment.

**Two obstacles, both real, both found by measurement:**

1. `sympl.RelaxationTendencyComponent` builds its tendency units with pint: `str(ureg('m s^-1') / ureg('s'))` gives `'1.0 meter / second ** 2'`, which **`unyt` cannot parse** — it raises `ValueError: Failed to parse unit '1.0 meter / second ** 2'`. So the component is unusable under `UnytBackend` as shipped. A two-line subclass overriding `tendency_properties` fixes it.
2. It requires `equilibrium_<quantity>` and `<quantity>_relaxation_timescale` in the state, which `get_default_state` does not provide. They have to be created with `sympl.get_backend().create_quantity(...)`.

Both are hidden in the helper, so the page teaches *Newtonian relaxation* rather than a units bug. The craft note on page 08 points at the shim and says why it exists.

Append to `_tour/stepping.py`:

```python
class UnytRelaxation(sympl.RelaxationTendencyComponent):
    """``sympl.RelaxationTendencyComponent`` with unyt-parseable tendency units.

    The stock component derives its tendency units with pint --
    ``str(ureg("m s^-1") / ureg("s"))`` is ``"1.0 meter / second ** 2"`` -- and
    unyt cannot parse that string, so constructing a state around it under
    ``climt.UnytBackend`` raises ``ValueError: Failed to parse unit``. The
    relaxation arithmetic itself is fine; only the units label is wrong for
    this backend. So: take the tendency units as an argument and hand unyt
    something it understands.
    """

    def __init__(self, quantity_name, units, tendency_units, **kwargs):
        self._tendency_units = tendency_units
        super(UnytRelaxation, self).__init__(quantity_name, units, **kwargs)

    @property
    def tendency_properties(self):
        return {self._quantity_name: {"dims": ["*"],
                                      "units": self._tendency_units}}


def wind_relaxation(state, speed, timescale_hours=24.0,
                    quantity_name="eastward_wind", initialise=True):
    """Give a column a wind, and keep it -- returns a tendency component.

    A single column has no dynamics, so nothing drives a wind, and
    ``SimpleBoundaryLayer``'s surface drag removes any wind you prescribe
    within days. Every configuration from page 8 on therefore needs a momentum
    source, or its surface fluxes are those of a dead-calm column.

    climt's examples do this by assigning the wind back at the end of every
    timestep (``examples/column_code_with_slab.py``). This is the same idea as
    a component: ``dx/dt = -(x - x_eq)/tau``, which composes with everything
    else, carries its units, and puts the assumption in the component list
    where a reader can see it.

    Writes ``equilibrium_<quantity>`` and ``<quantity>_relaxation_timescale``
    into ``state`` -- ``get_default_state`` does not know about them, because
    they belong to the relaxation rather than to any physics component.

    Args:
        state: the state to add the target fields to. Modified in place.
        speed: the wind the column relaxes toward, m/s.
        timescale_hours: relaxation timescale. Long compared with a timestep
            and short compared with the run, or it is either a hard reset or
            no forcing at all.
        quantity_name: the wind component to relax.
        initialise: also set the wind itself to ``speed``, which is what you
            want when building a fresh state. Pass ``False`` when attaching
            this to a state loaded from disk -- that state already carries a
            spun-up wind profile, sheared by the drag, and flattening it to a
            uniform ``speed`` would throw away part of the equilibrium.

    Returns:
        A tendency component to put in the tendency list.
    """
    backend = sympl.get_backend()
    template = state[quantity_name]

    def _field(name, value, units):
        return backend.create_quantity(
            np.full(template.values.shape, value), name, units, template.dims)

    state["equilibrium_" + quantity_name] = _field(
        "equilibrium_" + quantity_name, speed, "m s^-1")
    state[quantity_name + "_relaxation_timescale"] = _field(
        quantity_name + "_relaxation_timescale", timescale_hours * 3600.0, "s")
    if initialise:
        state[quantity_name].values[:] = speed

    return UnytRelaxation(quantity_name, "m s^-1", "m s^-2")
```

Add `import sympl` to the module's imports (it currently imports only `numpy` and `sympl.AdamsBashforth`).

Add the tests:

```python
def test_wind_relaxation_holds_the_wind_up(stepping):
    """Without a momentum source the column spins down; with one it does not.

    Measured without relaxation: the lowest-level wind reaches 0.000 m/s and
    the equilibrium is independent of how fast it started.
    """
    longwave = climt.CorkLongwaveRadiation(
        optics="correlated_k", table=PAGE7_GRAY_TABLE,
        diffusivity_factor=PAGE7_DIFFUSIVITY)
    surface = climt.SlabSurface()
    boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                               roughness_length=1e-3)
    state = climt.get_default_state([longwave, surface, boundary_layer],
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    state["ocean_mixed_layer_thickness"].values[:] = 2.0
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0

    relaxation = stepping.wind_relaxation(state, 8.0, timescale_hours=24.0)
    assert "equilibrium_eastward_wind" in state
    assert state["eastward_wind_relaxation_timescale"].attrs["units"] == "s"

    stepping.integrate([longwave, surface, relaxation], [boundary_layer],
                       state, climt.UnytTimeDelta(hours=1), 500)

    lowest = float(state["eastward_wind"].values[0, 0, 0])
    assert lowest > 1.0, (
        f"lowest-level wind {lowest:.3f} m/s — the relaxation is not holding "
        "the wind up against the surface drag")
    assert lowest < 8.0, (
        "the lowest level should sit below the target: drag is still acting, "
        "which is the point")


def test_unyt_relaxation_units_are_parseable(stepping):
    """Pins the reason the subclass exists.

    The stock component's tendency units come from pint and read
    '1.0 meter / second ** 2', which unyt refuses. If sympl ever emits
    something unyt can parse, this shim can go -- and this is where we notice.
    """
    import sympl as _sympl

    stock = _sympl.RelaxationTendencyComponent("eastward_wind", "m s^-1")
    assert "**" in stock.tendency_properties["eastward_wind"]["units"], (
        "sympl's tendency units are no longer pint-formatted — re-check "
        "whether stepping.UnytRelaxation is still needed")

    ours = stepping.UnytRelaxation("eastward_wind", "m s^-1", "m s^-2")
    assert ours.tendency_properties["eastward_wind"]["units"] == "m s^-2"
```

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k "relaxation" -v`
Expected: both pass. If `test_unyt_relaxation_units_are_parseable` fails on its first assertion, sympl changed and the shim may be removable — check before deleting it, because the state fields still need creating either way.

- [ ] **Step 9: Commit**

```bash
git add docs/modelling-tour/_tour/stepping.py tests/test_modelling_tour.py
git commit -m "feat(tour): _tour/stepping.py, the tranche's single time loop

Lifts integrate_to_equilibrium and the evolution figure out of the live-RCE
page's Quarto include and its .py twin into one importable module: UnytTimeDelta
instead of timedelta, steppers allowed in the snapshot loop, and recording
split from drawing so the recording is testable without matplotlib. Adds
wind_relaxation: a column has no dynamics, so its wind decays to zero under
the boundary layer's own drag unless something holds it up.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 5: Retire `radiative-transfer/09-live-rce.qmd` and consolidate the includes

The live-RCE page's demo is what page 07 will teach properly, with corrected numbers. It goes now, before the pages are written, so that nothing is written against a page that is about to disappear. Its loop and figure already survive as `_tour/stepping.py` (Task 4); its two physics claims survive as tests, moved here.

**The include consolidation.** `docs/_includes/climt-live-setup.qmd` (227 lines) had exactly one consumer: this page. It is deleted, leaving `climt-live-boot.qmd` as the site's single boot include for all twelve live pages. This deliberately undoes fix #2 of tranche 1's whole-branch review, which duplicated a few import lines across two includes so that two page families could not break each other. There is now one family.

Every inbound reference was checked. There are five, and none is a public link — see the spec for the full accounting.

**Files:**
- Delete: `docs/radiative-transfer/09-live-rce.qmd`
- Delete: `docs/_includes/climt-live-setup.qmd`
- Delete: `docs/radiative-transfer/_live/rce_helpers.py` (and the now-empty `_live/`)
- Delete: `tests/test_live_rce_demo.py`
- Move: `docs/radiative-transfer/_live/serve_wheel.py` → `scripts/serve_wheel.py`
- Modify: `docs/_includes/climt-live-boot.qmd`
- Modify: `docs/_quarto.yml`
- Modify: `docs/modelling-tour/04-gray-equilibrium-tested.qmd:305`
- Modify: `docs/modelling-tour/01-emissivity-spectrum.qmd:9`
- Move + rewrite: `scripts/experiments/render_live_rce_figure.py` → `scripts/experiments/render_page7_figure.py`
- Delete: `scripts/experiments/live_rce_axis_check.py`
- Modify: `tests/test_modelling_tour.py`

**Interfaces:**
- Consumes: `stepping.integrate` (Task 4).
- Produces: nothing new. `scripts/serve_wheel.py` keeps its CLI: `python scripts/serve_wheel.py <dir-with-the-wheel> [port]`.

- [ ] **Step 1: Move the two surviving physics tests into the tour suite**

The old file has four tests. `test_helper_twin_matches_the_page_include` dies with the twin contract it enforced — there is one copy now, so there is nothing to compare. `test_flux_diagnostics_survive_the_update_order` is already covered by `test_integrate_keeps_the_flux_diagnostics` from Task 4, but only for the gray table; the parametrised version moves. The two claims about spectral structure move as page 07's claims.

Add to `tests/test_modelling_tour.py`:

```python
# ------------------------------------------------------------------ page 7
#
# Inherited from the deleted tests/test_live_rce_demo.py, on the deleted
# radiative-transfer/09-live-rce.qmd. Same two claims, re-measured on page 7's
# configuration: UnytBackend rather than DataArrayBackend, nz=28 rather than 18,
# dt=12 h rather than 2 h, a 2 m slab rather than 5 m.

PAGE7_STEPS = 300      # far enough to separate the two columns unambiguously;
                       # the page itself runs longer. See the log below.

# The tables page 7 actually declares, and the diffusivity that goes with the
# gray one. NOT `single_band_gray_lw`, which is what the deleted page used:
# page 7's gray half must reproduce page 4's analytic profile, and that only
# works against the table page 4 calibrated (D = 2, D * sum(tau) = 4).
# Per the global constraint, a test asserts on the table its page declares.
PAGE7_GRAY_TABLE = "tour_gray_lw"
PAGE7_DIFFUSIVITY = 2.0
PAGE7_NONGREY_TABLE = "earth_low_res_lw"


def _page7_gray_column(nz=28, slab_depth=2.0):
    """Page 7's gray column: page 4's table and diffusivity, over a slab."""
    longwave = climt.CorkLongwaveRadiation(
        optics="correlated_k", table=PAGE7_GRAY_TABLE,
        diffusivity_factor=PAGE7_DIFFUSIVITY)
    surface = climt.SlabSurface()
    state = climt.get_default_state([longwave, surface],
                                    grid_state=get_grid(nx=1, ny=1, nz=nz))
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    return [longwave, surface], state


@pytest.fixture(scope="module")
def page7_columns():
    """Both of page 7's integrations, run once and shared."""
    import sympl as _sympl

    _sympl.set_backend(climt.UnytBackend())
    stepping_module = _load("stepping")
    timestep = climt.UnytTimeDelta(hours=12)

    gray_components, gray_state = _page7_gray_column()
    nongrey_components, nongrey_state = _gray_column(
        table=PAGE7_NONGREY_TABLE)
    return {
        PAGE7_GRAY_TABLE: stepping_module.integrate(
            gray_components, [], gray_state, timestep, PAGE7_STEPS),
        PAGE7_NONGREY_TABLE: stepping_module.integrate(
            nongrey_components, [], nongrey_state, timestep, PAGE7_STEPS),
    }


def _olr(state):
    return float(state["upwelling_longwave_flux_in_air"].values[-1, 0, 0])


def _surface_temperature(state):
    return float(state["surface_temperature"].values.ravel()[0])


@pytest.mark.slow
def test_page7_nongrey_column_is_the_more_efficient_radiator(page7_columns):
    """Page 7's comparison: spectral windows let the column radiate better.

    The non-grey table resolves atmospheric windows through which surface
    emission escapes more or less directly to space. The single-band gray
    column has no such windows, so at matched conditions the non-grey column
    emits more to space and runs cooler at the surface.
    """
    gray = page7_columns[PAGE7_GRAY_TABLE]
    nongrey = page7_columns[PAGE7_NONGREY_TABLE]

    assert _olr(nongrey) > _olr(gray) + 10.0, (
        f"non-grey OLR {_olr(nongrey):.1f} should exceed gray "
        f"{_olr(gray):.1f} W/m^2 — windows radiate to space more efficiently")
    assert _surface_temperature(nongrey) < _surface_temperature(gray) - 1.0, (
        f"non-grey surface {_surface_temperature(nongrey):.1f} K should be "
        f"cooler than gray {_surface_temperature(gray):.1f} K")


@pytest.mark.slow
def test_page7_nongrey_column_cools_faster_aloft(page7_columns):
    """Page 7's second comparison: a steeper temperature drop-off aloft.

    Band-resolved absorption concentrates cooling in the strongly absorbing
    bands high in the column, which the gray column smears out.
    """
    top_gradient = {}
    for name, state in page7_columns.items():
        temperature = state["air_temperature"].values[:, 0, 0]
        # Level index increases upward; mean of the top three level-to-level
        # differences, in K per level.
        top_gradient[name] = float(np.mean(np.diff(temperature)[-3:]))

    assert top_gradient[PAGE7_NONGREY_TABLE] < \
        top_gradient[PAGE7_GRAY_TABLE] - 1.0, (
        f"non-grey top-of-column gradient {top_gradient[PAGE7_NONGREY_TABLE]:.2f}"
        f" K/level should be markedly steeper than gray "
        f"{top_gradient[PAGE7_GRAY_TABLE]:.2f} K/level")


@pytest.mark.slow
@pytest.mark.parametrize("table", [PAGE7_GRAY_TABLE, PAGE7_NONGREY_TABLE])
def test_page7_flux_diagnostics_survive_the_update_order(page7_columns, table):
    """The update-order gotcha, on both of page 7's tables.

    ``stepping.integrate`` applies the prognostic state *before* the
    diagnostics. Reversing that clobbers the longwave fluxes SlabSurface
    consumes, and the surface heats without bound.
    """
    state = page7_columns[table]

    for flux in ("upwelling_longwave_flux_in_air",
                 "downwelling_longwave_flux_in_air"):
        assert flux in state, f"{flux} missing — diagnostics were clobbered"
        values = state[flux].values
        assert np.all(np.isfinite(values)), f"{flux} has non-finite values"
        assert np.all(values >= 0.0), f"{flux} must be non-negative"

    temperature = state["air_temperature"].values
    assert np.all(np.isfinite(temperature))
    assert temperature.min() > 150.0 and temperature.max() < 400.0, (
        f"air temperature ran to [{temperature.min():.1f}, "
        f"{temperature.max():.1f}] K — suspect a flux-coupling regression")
    assert 150.0 < _surface_temperature(state) < 400.0
```

- [ ] **Step 2: Measure the thresholds, then tighten them**

The three thresholds above (`+10.0`, `-1.0`, `-1.0`) are deliberately loose placeholders on a configuration nobody has run. Measure the real separation and record it:

```bash
conda run -n climt python -m pytest tests/test_modelling_tour.py -k page7 -v -m slow -s
conda run -n climt python -c "
import sys, sympl, climt, numpy as np
sys.path.insert(0, 'docs/modelling-tour/_tour'); import stepping
sympl.set_backend(climt.UnytBackend())
for table, D in (('tour_gray_lw', 2.0), ('earth_low_res_lw', 1.66)):
    lw = climt.CorkLongwaveRadiation(optics='correlated_k', table=table,
                                     diffusivity_factor=D)
    surf = climt.SlabSurface()
    s = climt.get_default_state([lw, surf], grid_state=climt.get_grid(nx=1, ny=1, nz=28))
    s['ocean_mixed_layer_thickness'].values[:] = 2.0
    s['downwelling_shortwave_flux_in_air'].values[:] = 0.0
    s['downwelling_shortwave_flux_in_air'].values[0, ...] = 240.0
    s['upwelling_shortwave_flux_in_air'].values[:] = 0.0
    stepping.integrate([lw, surf], [], s, climt.UnytTimeDelta(hours=12), 300)
    T = s['air_temperature'].values[:, 0, 0]
    print(f'{table:22s} OLR {float(s[\"upwelling_longwave_flux_in_air\"].values[-1,0,0]):7.2f}'
          f'  Tsurf {float(s[\"surface_temperature\"].values.ravel()[0]):7.2f}'
          f'  top grad {float(np.mean(np.diff(T)[-3:])):6.2f} K/level')"
```

Write the two lines of output into a `### Log — page 7's two columns as measured` section at the end of this task, then set each threshold to roughly **half the measured separation**, so the tests check the physics rather than fitting the numbers. If a separation is smaller than the placeholder, that is a finding, not a threshold to relax — page 07's whole reveal rests on it, so stop and reconcile with the spec.

- [ ] **Step 3: Delete the page, the include, the twin and the old test**

```bash
git rm docs/radiative-transfer/09-live-rce.qmd
git rm docs/_includes/climt-live-setup.qmd
git rm docs/radiative-transfer/_live/rce_helpers.py
git rm tests/test_live_rce_demo.py
git mv docs/radiative-transfer/_live/serve_wheel.py scripts/serve_wheel.py
rm -rf docs/radiative-transfer/_live
```

- [ ] **Step 4: Fix `serve_wheel.py`'s own usage lines**

`scripts/serve_wheel.py`'s docstring names its old path twice (lines 8 and 13). Update both to `python scripts/serve_wheel.py <dir-with-the-wheel> [port]` and `python scripts/serve_wheel.py /tmp/climt_wh 8912`.

Run: `conda run -n climt python scripts/serve_wheel.py --help 2>&1 | head -5` (or just read the file back — it is 36 lines).
Expected: no stale `docs/radiative-transfer/_live` path remains: `grep -rn "_live" scripts/serve_wheel.py` returns nothing.

- [ ] **Step 5: Move the unreleased-wheel instructions into the boot include**

That instruction — how to preview an unreleased wheel — lived on the deleted page, and `01-emissivity-spectrum.qmd:9` pointed at it. It moves into `docs/_includes/climt-live-boot.qmd`'s header comment, which is where anyone looking for it now is.

Replace the whole HTML comment block at the top of `docs/_includes/climt-live-boot.qmd` with:

```
<!--
The site's single quarto-live boot include. Every `format: live-html` page
includes it once, near the top, with a Quarto include shortcode pointing here
(see modelling-tour/01-emissivity-spectrum.qmd for the exact line). Do NOT
repeat the literal include shortcode in this comment: Quarto resolves include
shortcodes at the text level, before HTML comments are stripped, so a literal
one here would make this file include itself.

This include is deliberately minimal: imports and nothing else. The time
loops and figures the RCE pages need live in `modelling-tour/_tour/stepping.py`
and arrive by the `resources:` mechanism, declared only by the pages that use
them — so pages 01-06 download nothing they do not call. An earlier design put
those ~200 lines in a second include, whose cell carried `#| edit: false`,
which quarto-live renders **read-only, not hidden**: every page including it
opened with four screens of code it never ran.

The climt wheel itself is installed by micropip at document setup via the
including page's `pyodide: packages:` front-matter block — resolved from PyPI,
pinned to a released version.

PREVIEWING UNRELEASED CHANGES. To run a page against your working copy rather
than the released wheel, build the pure wheel and serve it *with CORS headers*
(plain `python -m http.server` will NOT work, and GitHub release assets do not
send CORS either), then point the page's front matter at the local URL:

    CLIMT_PURE_PYTHON=1 python -m pip wheel . --no-deps -w /tmp/climt_wh
    python scripts/serve_wheel.py /tmp/climt_wh 8912
-->
```

Leave the `{pyodide}` cell below it untouched.

- [ ] **Step 6: Re-point the two cross-references**

`docs/modelling-tour/01-emissivity-spectrum.qmd:9` — replace

```
    # radiative-transfer/09-live-rce.qmd for how to preview an unreleased wheel.
```

with

```
    # ../_includes/climt-live-boot.qmd for how to preview an unreleased wheel.
```

`docs/modelling-tour/04-gray-equilibrium-tested.qmd:305` — replace

```
[The live RCE demo](../radiative-transfer/09-live-rce.qmd) does it the other way,
if you want to watch.
```

with

```
[Page 7](07-stepping-a-column.qmd) does it the other way, if you want to watch —
and finds this same profile without being told it.
```

That link is forward-pointing within the tour, so it dangles until Task 10 creates page 07. Task 16's link check is where that is caught if Task 10 slips.

- [ ] **Step 7: Update `docs/_quarto.yml`**

Two edits. Remove the sidebar entry:

```yaml
        - radiative-transfer/09-live-rce.qmd
```

And rewrite the comment above the top-level `theme:` block (it names the deleted page as one of two `live-html` families):

```yaml
# HTML appearance, declared at the *top level* rather than under `format: html:`
# on purpose. The Modelling Tour pages set `format: live-html` -- a separate
# format contributed by the r-wasm/live extension -- and Quarto merges a
# `format: html:` block only into documents whose format is literally `html`.
# Those pages therefore rendered against bare `cosmo` with no climt.scss,
# visibly out of step with the rest of the site. Top-level keys are merged into
# every format the project produces, `live-html` included.
```

Do **not** add the six new pages to the sidebar here — each page's own task adds its own line, so a half-finished branch never advertises a page that does not exist.

- [ ] **Step 8: Fix the two experiment scripts that read the deleted include**

`scripts/experiments/render_live_rce_figure.py` extracts the plotting cell out of `climt-live-setup.qmd`; `scripts/experiments/live_rce_axis_check.py` does something similar. Neither can work now, and both are superseded: the figure code is importable.

Rename the first to match what it renders and rewrite it to import rather than extract:

```bash
git mv scripts/experiments/render_live_rce_figure.py scripts/experiments/render_page7_figure.py
```

```python
"""Render page 7's evolution figure natively, from the shipped code itself.

The plotting used to live inside a ```{pyodide} cell in a Quarto include and
had to be extracted and exec'd. It is now `_tour/stepping.py`, so this just
imports it -- what lands in debug_data/ is what a reader's browser draws.

    conda run -n climt python scripts/experiments/render_page7_figure.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402  (must follow the backend choice)
import sympl                       # noqa: E402

import climt                       # noqa: E402

sys.path.insert(0, os.path.join("docs", "modelling-tour", "_tour"))
import stepping                    # noqa: E402

SOLAR = 240.0
NZ, DT_HOURS, N_STEPS, SLAB_DEPTH = 28, 12, 675, 2.0
OUT_DIR = "debug_data"


def column(table, diffusivity):
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k", table=table,
                                           diffusivity_factor=diffusivity)
    surface = climt.SlabSurface()
    state = climt.get_default_state(
        [longwave, surface], grid_state=climt.get_grid(nx=1, ny=1, nz=NZ))
    state["ocean_mixed_layer_thickness"].values[:] = SLAB_DEPTH
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    return [longwave, surface], state


def main():
    sympl.set_backend(climt.UnytBackend())
    os.makedirs(OUT_DIR, exist_ok=True)
    # Page 7's own tables and diffusivities: `tour_gray_lw` with D = 2 is what
    # reproduces page 4's analytic profile, and is what the page runs.
    for table, diffusivity, slug, title in (
            ("tour_gray_lw", 2.0, "gray", "Gray radiative equilibrium"),
            ("earth_low_res_lw", 1.66, "nongrey",
             "Non-grey radiative equilibrium (14 bands)")):
        components, state = column(table, diffusivity)
        state, history = stepping.integrate_with_history(
            components, [], state, climt.UnytTimeDelta(hours=DT_HOURS),
            N_STEPS)
        stepping.draw_evolution(history, state, SOLAR, title=title)
        path = os.path.join(OUT_DIR, f"page7_{slug}.png")
        plt.savefig(path, dpi=110, bbox_inches="tight")
        plt.close("all")
        print(f"wrote {path}  Tsurf="
              f"{float(state['surface_temperature'].values.ravel()[0]):.2f} K")


if __name__ == "__main__":
    main()
```

`live_rce_axis_check.py` was written to settle tranche 1's log-pressure-axis question, which is settled and encoded in `draw_evolution`'s comment. Delete it: `git rm scripts/experiments/live_rce_axis_check.py`.

- [ ] **Step 9: Verify nothing still points at what was deleted**

```bash
grep -rn "09-live-rce\|climt-live-setup\|rce_helpers\|_live/" \
     docs/ scripts/ tests/ Makefile 2>/dev/null \
  | grep -v "docs/_site/\|docs/.quarto/\|superpowers/plans/\|superpowers/specs/\|superpowers/HANDOFF"
```

Expected: **no output.** Hits under `superpowers/` are historical records and are left alone, exactly as plans are; hits under `_site/` and `.quarto/` are build output cleared by the next render.

- [ ] **Step 10: Render the site and run the suite**

```bash
cd docs && quarto render && cd ..
conda run -n climt python -m pytest tests/test_modelling_tour.py -v
conda run -n climt python scripts/experiments/render_page7_figure.py
```

Expected: the render completes with no "file not found" for the deleted include; `docs/_site/radiative-transfer/09-live-rce.html` is gone; the sidebar shows no live-RCE entry; the tests pass (page 7's three are `slow`, so add `-m slow` to see them); and two PNGs land in `debug_data/`.

- [ ] **Step 11: Commit**

```bash
git add -A docs scripts tests
git commit -m "refactor(docs): the live-RCE page retires into the tour

Its loop and figure are _tour/stepping.py now, and its two physics claims are
tests in tests/test_modelling_tour.py, re-measured on page 7's configuration.
climt-live-setup.qmd had exactly one consumer and goes with it, leaving
climt-live-boot.qmd as the site's single boot include; the unreleased-wheel
instructions move into its header. serve_wheel.py moves to scripts/.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 6: `_tour/budgets.py` — energy budgets as the convergence test

Pages 11 and 12 do not decide they have reached equilibrium by looking at a curve. They read the top-of-atmosphere and surface energy budgets and check that both are near zero. That is also the residual test that guards the shipped states in CI (Task 9), and the acceptance criterion the generator script uses (Task 9). One module, three consumers.

**Files:**
- Create: `docs/modelling-tour/_tour/budgets.py`
- Test: `tests/test_modelling_tour.py` (modify)

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `budgets.pressure_thickness(state) -> (nz,) ndarray` in Pa
  - `budgets.column_enthalpy(state) -> float`, J m⁻², dry plus latent
  - `budgets.toa_imbalance(state) -> float`, W m⁻², positive = net energy into the planet
  - `budgets.surface_imbalance(state) -> float`, W m⁻², positive = net energy into the surface
  - `budgets.precipitation_rate(state, timestep) -> float`, mm day⁻¹, convective plus grid-scale
  - `budgets.evaporation_rate(state) -> float`, mm day⁻¹, implied by the surface latent heat flux
  - `budgets.summary(state) -> dict` with keys `surface_temperature`, `olr`, `absorbed_shortwave`, `toa_imbalance`, `surface_imbalance`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_modelling_tour.py`:

```python
@pytest.fixture
def budgets():
    return _load("budgets")


def test_toa_imbalance_is_absorbed_solar_minus_olr(budgets):
    """With no atmospheric shortwave absorption, TOA net = surface SW - OLR.

    Every page in this tranche prescribes the absorbed shortwave at the
    surface and runs no shortwave component, so the planet's absorbed
    shortwave *is* the surface's. The function must read that from the state
    rather than take it as an argument, or a page can quote a budget that
    disagrees with its own configuration.
    """
    components, state = _gray_column()
    stepping_module = _load("stepping")
    stepping_module.integrate(components, [], state,
                              climt.UnytTimeDelta(hours=12), 5)

    olr = float(state["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    assert budgets.toa_imbalance(state) == pytest.approx(SOLAR - olr, abs=1e-9)


def test_toa_imbalance_closes_on_a_converged_column(budgets):
    """The whole point: near equilibrium, the budget is near zero."""
    components, state = _gray_column()
    _load("stepping").integrate(components, [], state,
                                climt.UnytTimeDelta(hours=12), 700)
    assert abs(budgets.toa_imbalance(state)) < 1.0


def test_surface_imbalance_counts_every_term(budgets):
    """Radiation in, radiation out, and the two turbulent fluxes out."""
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table="single_band_gray_lw")
    surface = climt.SlabSurface()
    boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="bulk")
    state = climt.get_default_state([longwave, surface, boundary_layer],
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    state["ocean_mixed_layer_thickness"].values[:] = 2.0
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    _load("stepping").integrate([longwave, surface], [boundary_layer], state,
                                climt.UnytTimeDelta(hours=1), 10)

    expected = (
        SOLAR
        + float(state["downwelling_longwave_flux_in_air"].values[0, 0, 0])
        - float(state["upwelling_longwave_flux_in_air"].values[0, 0, 0])
        - float(state["surface_upward_sensible_heat_flux"].values.ravel()[0])
        - float(state["surface_upward_latent_heat_flux"].values.ravel()[0])
    )
    assert budgets.surface_imbalance(state) == pytest.approx(expected, abs=1e-9)


def test_surface_imbalance_works_without_turbulent_fluxes(budgets):
    """Pages 7 has no boundary layer; the two flux terms are then zero."""
    components, state = _gray_column()
    _load("stepping").integrate(components, [], state,
                                climt.UnytTimeDelta(hours=12), 5)
    expected = (
        SOLAR
        + float(state["downwelling_longwave_flux_in_air"].values[0, 0, 0])
        - float(state["upwelling_longwave_flux_in_air"].values[0, 0, 0])
    )
    assert budgets.surface_imbalance(state) == pytest.approx(expected, abs=1e-9)


def test_column_enthalpy_matches_the_conservation_suite_formula(budgets):
    """Same integral tests/test_conservation.py uses, so page 10 can quote it.

    Cp_dry * T + Lv * q, mass-weighted by dp/g. If these two ever disagree,
    page 10's enthalpy-conservation claim is guarded by a test measuring
    something else.
    """
    from sympl import get_constant

    components, state = _gray_column()
    state["specific_humidity"].values[:] = 4e-3

    Cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure",
                       "J/kg/degK")
    Lv = get_constant("latent_heat_of_condensation", "J/kg")
    g = get_constant("gravitational_acceleration", "m/s^2")
    p_int = state["air_pressure_on_interface_levels"].values[:, 0, 0]
    dp = p_int[:-1] - p_int[1:]
    T = state["air_temperature"].values[:, 0, 0]
    q = state["specific_humidity"].values[:, 0, 0]
    expected = float(np.sum((Cpd * T + Lv * q) * dp / g))

    assert budgets.column_enthalpy(state) == pytest.approx(expected, rel=1e-12)


def test_evaporation_rate_inverts_the_latent_heat_flux(budgets):
    """Page 12 closes the moisture budget against this."""
    from sympl import get_constant

    components, state = _gray_column()
    state["surface_upward_latent_heat_flux"].values[:] = 80.0
    Lv = get_constant("latent_heat_of_condensation", "J/kg")
    assert budgets.evaporation_rate(state) == pytest.approx(
        80.0 / float(Lv) * 86400.0, rel=1e-12)


def test_summary_reports_the_five_numbers_the_pages_quote(budgets):
    components, state = _gray_column()
    _load("stepping").integrate(components, [], state,
                                climt.UnytTimeDelta(hours=12), 5)
    report = budgets.summary(state)
    assert set(report) == {"surface_temperature", "olr", "absorbed_shortwave",
                           "toa_imbalance", "surface_imbalance"}
    assert all(np.isfinite(value) for value in report.values())
    assert report["absorbed_shortwave"] == pytest.approx(SOLAR)
```

- [ ] **Step 2: Run it to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k "budget or enthalpy or evaporation or summary or imbalance" -v`
Expected: FAIL — `_tour/budgets.py` does not exist.

- [ ] **Step 3: Write `_tour/budgets.py`**

```python
"""Read a column's energy and water budgets out of its state.

Pages 11 and 12 do not decide a column has equilibrated by squinting at a
curve; they read these numbers. The same functions are the acceptance criterion
in ``scripts/generate_tour_equilibria.py`` and the residual test that keeps the
shipped equilibrium states from going stale, so "converged" means one thing
across the pages, the generator and CI.

**Sign convention: positive is into the box.** ``toa_imbalance`` positive means
the planet is gaining energy and will warm; ``surface_imbalance`` positive
means the surface is gaining energy.

**Level indexing.** Index 0 is the bottom (the surface interface, or the
lowest model layer); index -1 is the top. climt states are ordered this way
throughout; a page that gets it backwards reports the OLR as the downwelling
longwave at the ground.

**The shortwave.** No page in this tranche runs a shortwave component -- the
absorbed shortwave is written into ``downwelling_shortwave_flux_in_air`` at the
surface interface and the upwelling shortwave is zero. So the atmosphere
absorbs no sunlight, and the planet's absorbed shortwave is exactly the
surface's. That is why :func:`toa_imbalance` reads it out of the state instead
of taking it as an argument: a page cannot then quote a budget computed against
a ``SOLAR`` it did not actually use.
"""
import numpy as np
from sympl import get_constant

DAY_SECONDS = 86400.0


def _column(state, name):
    """A quantity's single column as a plain 1-D float array."""
    return np.asarray(state[name].values, dtype=float).reshape(
        state[name].values.shape[0], -1)[:, 0]


def _scalar(state, name, default=0.0):
    """A surface scalar as a plain float, or ``default`` if absent."""
    if name not in state:
        return float(default)
    return float(np.asarray(state[name].values, dtype=float).ravel()[0])


def pressure_thickness(state):
    """Layer mass in pressure units, (nz,) in Pa, bottom first."""
    p_int = _column(state, "air_pressure_on_interface_levels")
    return p_int[:-1] - p_int[1:]


def column_enthalpy(state):
    """Column-integrated moist enthalpy, J m^-2.

    ``Cp_dry * T + Lv * q``, mass-weighted by ``dp / g`` -- the same integral
    ``tests/test_conservation.py`` uses, so page 10's claim that dry convective
    adjustment conserves it is guarded by a test that already exists.
    """
    Cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure",
                       "J/kg/degK")
    Lv = get_constant("latent_heat_of_condensation", "J/kg")
    g = get_constant("gravitational_acceleration", "m/s^2")

    T = _column(state, "air_temperature")
    q = (_column(state, "specific_humidity")
         if "specific_humidity" in state else np.zeros_like(T))
    return float(np.sum((float(Cpd) * T + float(Lv) * q)
                        * pressure_thickness(state) / float(g)))


def absorbed_shortwave(state):
    """Net downward shortwave at the surface, W m^-2 -- and at the TOA.

    Equal at both, because no page here puts a shortwave absorber in the air.
    """
    return (_column(state, "downwelling_shortwave_flux_in_air")[0]
            - _column(state, "upwelling_shortwave_flux_in_air")[0])


def olr(state):
    """Outgoing longwave radiation: upwelling longwave at the top, W m^-2."""
    return float(_column(state, "upwelling_longwave_flux_in_air")[-1])


def toa_imbalance(state):
    """Net downward energy flux at the top of the atmosphere, W m^-2.

    Positive means the planet is gaining energy. At equilibrium this is zero,
    and how close to zero is the convergence test.
    """
    return float(absorbed_shortwave(state) - olr(state))


def surface_imbalance(state):
    """Net downward energy flux into the surface, W m^-2.

    Shortwave in, longwave down in, longwave up out, and the two turbulent
    fluxes out. The turbulent terms are zero on pages that run no boundary
    layer, which is what makes those pages' surface--air discontinuity the
    thing page 8 then goes and erases.
    """
    return float(
        absorbed_shortwave(state)
        + _column(state, "downwelling_longwave_flux_in_air")[0]
        - _column(state, "upwelling_longwave_flux_in_air")[0]
        - _scalar(state, "surface_upward_sensible_heat_flux")
        - _scalar(state, "surface_upward_latent_heat_flux")
    )


def evaporation_rate(state):
    """Evaporation implied by the surface latent heat flux, mm day^-1.

    ``LHF / Lv`` is a mass flux in kg m^-2 s^-1, which over water is mm s^-1.
    Page 12 closes its moisture budget by comparing this with
    :func:`precipitation_rate`.
    """
    Lv = get_constant("latent_heat_of_condensation", "J/kg")
    return (_scalar(state, "surface_upward_latent_heat_flux")
            / float(Lv) * DAY_SECONDS)


def precipitation_rate(state, timestep):
    """Total precipitation, mm day^-1: convective plus grid-scale.

    ``EmanuelConvectionPython`` reports ``convective_precipitation_rate``
    already in mm day^-1. ``GridScaleCondensation`` reports
    ``precipitation_amount`` in kg m^-2 accumulated over the step it was
    called with, so it needs the timestep to become a rate. Either may be
    absent; a page with neither gets 0.0.

    Args:
        state: a climt state after at least one step.
        timestep: the ``UnytTimeDelta`` the step was taken with.
    """
    convective = _scalar(state, "convective_precipitation_rate")
    grid_scale = _scalar(state, "precipitation_amount")
    seconds = float(timestep.total_seconds())
    return float(convective + grid_scale * DAY_SECONDS / seconds)


def summary(state):
    """The five numbers pages 11 and 12 print beside every experiment."""
    return dict(
        surface_temperature=_scalar(state, "surface_temperature"),
        olr=olr(state),
        absorbed_shortwave=float(absorbed_shortwave(state)),
        toa_imbalance=toa_imbalance(state),
        surface_imbalance=surface_imbalance(state),
    )
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k "budget or enthalpy or evaporation or summary or imbalance" -v`
Expected: all pass. `test_toa_imbalance_closes_on_a_converged_column` takes ~2 s (700 gray steps at 2.4 ms) — under the `slow` threshold, so leave it unmarked.

- [ ] **Step 5: Commit**

```bash
git add docs/modelling-tour/_tour/budgets.py tests/test_modelling_tour.py
git commit -m "feat(tour): _tour/budgets.py, one definition of converged

TOA and surface energy budgets, column enthalpy and the moisture budget, read
out of the state. The pages, the equilibrium generator and the CI residual
test all use these, so 'converged' means the same thing in all three.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 7: `_tour/states.py` — save a model state, ship it, load it back

Pages 11 and 12 load an equilibrium rather than spinning one up. That needs a format, and the format needs provenance, because an equilibrium is only meaningful attached to the configuration that produced it.

**Reconstitution is deliberately strict.** `load` does not build a state out of the file. It calls `get_default_state` for the page's own component list and *overwrites* the arrays, checking each quantity's dims and units as it goes. So if a future sympl or climt renames a quantity, changes its dimension order, or changes its units, the page fails at a named quantity with a message — instead of quietly producing a state with the wrong shape and an equilibrium that is not one.

**Files:**
- Create: `docs/modelling-tour/_tour/states.py`
- Test: `tests/test_modelling_tour.py` (modify)

**Interfaces:**
- Consumes: `assets.locate` (Task 3).
- Produces:
  - `states.save(path, state, provenance) -> str` — writes `.npz`, returns the path
  - `states.load(name, components, grid_state=None, base_url="_data") -> (state, provenance)` — `name` is a file name resolved through `assets.locate`, or a path
  - `states.describe(provenance) -> str` — the block pages print above their first figure
  - `states.SKIPPED` — the tuple of state keys not saved (`"time"`, and coordinate quantities rebuilt by `get_default_state`)

- [ ] **Step 1: Write the failing test**

Add to `tests/test_modelling_tour.py`:

```python
@pytest.fixture
def states():
    return _load("states")


def _page11_components():
    """Page 11's exact component list, used by the round-trip tests."""
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table="earth_low_res_lw")
    surface = climt.SlabSurface()
    boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="bulk")
    adjustment = climt.DryConvectiveAdjustment()
    return [longwave, surface, boundary_layer, adjustment]


def test_state_round_trips_values_dims_and_units(states, tmp_path):
    """Save a state, load it back, and get the same state.

    Values, dimension names and units all have to survive. Values alone is not
    enough: a page that reloads a state with the right numbers under the wrong
    dims produces a plausible, wrong figure.
    """
    components = _page11_components()
    state = climt.get_default_state(components,
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    state["air_temperature"].values[:] = np.linspace(288.0, 200.0, 28
                                                     ).reshape(28, 1, 1)
    state["surface_temperature"].values[:] = 291.5
    state["specific_humidity"].values[:] = 3e-3

    path = states.save(str(tmp_path / "probe.npz"), state,
                       dict(note="round-trip probe"))
    reloaded, provenance = states.load(
        path, _page11_components(), grid_state=get_grid(nx=1, ny=1, nz=28))

    assert provenance["note"] == "round-trip probe"
    for name in ("air_temperature", "surface_temperature", "specific_humidity",
                 "air_pressure", "ocean_mixed_layer_thickness"):
        np.testing.assert_allclose(reloaded[name].values, state[name].values)
        assert tuple(reloaded[name].dims) == tuple(state[name].dims)
        assert reloaded[name].attrs["units"] == state[name].attrs["units"]


def test_load_records_the_climt_version_and_the_configuration(states, tmp_path):
    """Provenance travels with the state or the state means nothing."""
    components = _page11_components()
    state = climt.get_default_state(components,
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    path = states.save(str(tmp_path / "probe.npz"), state,
                       dict(table="earth_low_res_lw", nz=28, dt_hours=1.0,
                            n_steps=1700, slab_depth_m=2.0, solar=SOLAR,
                            co2_ppm=330.0,
                            components=["CorkLongwaveRadiation", "SlabSurface"],
                            toa_imbalance=0.02, surface_imbalance=-0.01))
    _, provenance = states.load(path, _page11_components(),
                                grid_state=get_grid(nx=1, ny=1, nz=28))

    assert provenance["climt_version"] == climt.__version__
    assert provenance["saved_at"]                      # ISO timestamp, non-empty
    assert provenance["table"] == "earth_low_res_lw"
    assert provenance["nz"] == 28
    assert provenance["components"][0] == "CorkLongwaveRadiation"


def test_load_fails_loudly_on_a_dimension_mismatch(states, tmp_path):
    """The whole reason load rebuilds instead of unpickling.

    A state saved at nz=28 must not quietly load into an nz=20 grid.
    """
    components = _page11_components()
    state = climt.get_default_state(components,
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    path = states.save(str(tmp_path / "probe.npz"), state, {})

    with pytest.raises(ValueError, match="air_temperature"):
        states.load(path, _page11_components(),
                    grid_state=get_grid(nx=1, ny=1, nz=20))


def test_load_fails_loudly_on_a_missing_quantity(states, tmp_path):
    """A component list that needs something the file does not carry."""
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table="earth_low_res_lw")
    thin_state = climt.get_default_state(
        [longwave], grid_state=get_grid(nx=1, ny=1, nz=28))
    path = states.save(str(tmp_path / "thin.npz"), thin_state, {})

    with pytest.raises(ValueError, match="not in the saved state"):
        states.load(path, _page11_components(),
                    grid_state=get_grid(nx=1, ny=1, nz=28))


def test_load_returns_none_for_a_missing_asset(states):
    """A page whose data asset did not stage degrades, it does not crash."""
    assert states.load("no_such_equilibrium.npz", _page11_components(),
                       grid_state=get_grid(nx=1, ny=1, nz=28)) == (None, None)


def test_describe_names_the_configuration(states, tmp_path):
    text = states.describe(dict(
        climt_version="0.31.0", saved_at="2026-08-20T10:00:00",
        table="earth_low_res_lw", nz=28, dt_hours=1.0, n_steps=1700,
        slab_depth_m=2.0, solar=240.0, co2_ppm=330.0,
        wind_m_s=5.0, wind_timescale_hours=24.0, roughness_length_m=1e-3,
        components=["CorkLongwaveRadiation", "SlabSurface"],
        toa_imbalance=0.02, surface_imbalance=-0.01))
    for token in ("0.31.0", "earth_low_res_lw", "28", "1700", "240", "330",
                  "5.0 m/s"):
        assert token in text
```

- [ ] **Step 2: Run it to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k "states or round_trip or provenance or describe" -v`
Expected: FAIL — `_tour/states.py` does not exist.

- [ ] **Step 3: Write `_tour/states.py`**

```python
"""Save a climt state to a small .npz, ship it, and load it back.

Pages 11 and 12 perturb an equilibrium rather than spending several thousand
in-browser steps finding one. The equilibrium therefore has to travel: as a
file in ``docs/modelling-tour/_data/``, published and staged by the same
mechanism as the spectrum table (see ``_tour/assets.py``).

**An equilibrium is only meaningful attached to its configuration**, so every
file carries provenance -- climt version, table, grid, timestep, step count,
slab depth, prescribed shortwave, CO2, the component list, and the TOA and
surface imbalances it finished at. :func:`describe` turns that into the block
the pages print above their first figure, so a reader can never be looking at
an equilibrium without also looking at what produced it.

**Loading rebuilds; it does not unpickle.** :func:`load` calls
``get_default_state`` for the caller's own component list and overwrites the
arrays, checking dims and units per quantity. A sympl or climt change that
renames a quantity, reorders its dimensions or changes its units then fails at
a named quantity with a message -- instead of yielding a state that looks fine
and is not an equilibrium.
"""
import datetime
import json

import numpy as np
from climt import get_default_state

import assets

# Rebuilt by get_default_state, or not an array at all. Saving them would only
# create a second place for the grid to be defined.
SKIPPED = ("time",)


def save(path, state, provenance):
    """Write ``state`` and its provenance to ``path`` as an .npz.

    Args:
        path: output file path.
        state: a climt state dict.
        provenance: dict of configuration facts. ``climt_version`` and
            ``saved_at`` are added here; anything else is the caller's.

    Returns:
        ``path``.
    """
    import climt

    meta = dict(provenance)
    meta["climt_version"] = climt.__version__
    meta["saved_at"] = datetime.datetime.utcnow().isoformat(timespec="seconds")

    arrays = {}
    layout = {}
    for name, quantity in state.items():
        if name in SKIPPED:
            continue
        values = np.asarray(quantity.values)
        arrays["value__" + name] = values
        layout[name] = dict(dims=list(quantity.dims),
                            units=quantity.attrs.get("units"))

    arrays["__meta__"] = np.array(json.dumps(meta))
    arrays["__layout__"] = np.array(json.dumps(layout))
    np.savez_compressed(path, **arrays)
    return path


def load(name, components, grid_state=None, base_url=assets.DEFAULT_BASE_URL):
    """Rebuild a saved state for ``components``.

    Args:
        name: the asset's file name (resolved through ``_tour/assets.py``) or
            a direct path.
        components: the component list the page will run. The state is built
            for *these*, not for whatever produced the file.
        grid_state: passed to ``get_default_state``; must match the saved grid.
        base_url: directory holding the asset, relative to the page.

    Returns:
        ``(state, provenance)``, or ``(None, None)`` if the asset is not
        available -- so a page whose resource failed to stage can say so
        rather than traceback.

    Raises:
        ValueError: if the file is there but does not fit ``components`` --
            a missing quantity, or a dims/units mismatch.
    """
    path = name if "/" in name or "\\" in name else None
    if path is None:
        path = assets.locate(name, base_url)
    elif not _exists(path):
        path = assets.locate(name, base_url)
    if path is None:
        return None, None

    with np.load(path, allow_pickle=False) as handle:
        meta = json.loads(str(handle["__meta__"]))
        layout = json.loads(str(handle["__layout__"]))
        saved = {key[len("value__"):]: handle[key][...]
                 for key in handle.files if key.startswith("value__")}

    state = get_default_state(components, grid_state=grid_state)

    for quantity_name, quantity in state.items():
        if quantity_name in SKIPPED:
            continue
        if quantity_name not in saved:
            raise ValueError(
                f"{quantity_name!r} is required by these components but is "
                f"not in the saved state {path!r}. The file was written for a "
                f"different component list ({meta.get('components')}); "
                "regenerate it with scripts/generate_tour_equilibria.py.")

        expected = layout[quantity_name]
        if tuple(expected["dims"]) != tuple(quantity.dims):
            raise ValueError(
                f"{quantity_name!r} was saved with dims "
                f"{tuple(expected['dims'])} but this state wants "
                f"{tuple(quantity.dims)}. Check grid_state (the file says "
                f"nz={meta.get('nz')}).")
        if saved[quantity_name].shape != quantity.values.shape:
            raise ValueError(
                f"{quantity_name!r} was saved with shape "
                f"{saved[quantity_name].shape} but this state wants "
                f"{quantity.values.shape}. Check grid_state (the file says "
                f"nz={meta.get('nz')}).")
        if expected["units"] != quantity.attrs.get("units"):
            raise ValueError(
                f"{quantity_name!r} was saved in {expected['units']!r} but "
                f"this state wants {quantity.attrs.get('units')!r}. A units "
                "convention changed under this file; regenerate it.")

        quantity.values[...] = saved[quantity_name]

    return state, meta


def describe(provenance):
    """A one-block human summary of what an equilibrium state is.

    Pages print this above their first figure. It is deliberately plain text
    rather than markdown: it goes through ``print()`` in a ``{pyodide}`` cell.
    """
    components = ", ".join(provenance.get("components", []))
    return "\n".join([
        "equilibrium state, as shipped",
        "  climt          {}".format(provenance.get("climt_version")),
        "  saved          {}".format(provenance.get("saved_at")),
        "  LW table       {}".format(provenance.get("table")),
        "  grid           nz = {}, single column".format(provenance.get("nz")),
        "  spin-up        {} steps at dt = {} h".format(
            provenance.get("n_steps"), provenance.get("dt_hours")),
        "  slab depth     {} m".format(provenance.get("slab_depth_m")),
        "  absorbed SW    {} W/m^2 (prescribed at the surface)".format(
            provenance.get("solar")),
        "  CO2            {} ppm".format(provenance.get("co2_ppm")),
        "  wind           relaxed to {} m/s over {} h, z0 = {} m".format(
            provenance.get("wind_m_s"), provenance.get("wind_timescale_hours"),
            provenance.get("roughness_length_m")),
        "  components     {}".format(components),
        "  residual       TOA {:+.3f}, surface {:+.3f} W/m^2".format(
            provenance.get("toa_imbalance", float("nan")),
            provenance.get("surface_imbalance", float("nan"))),
    ])


def _exists(path):
    import os

    return os.path.isfile(path)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k "states or round_trip or provenance or describe" -v`
Expected: all pass.

If `test_load_fails_loudly_on_a_dimension_mismatch` does not raise, check which guard should have fired — for a pure `nz` change the *shape* check catches it, since `dims` are names and stay the same. Both guards are needed; do not delete the one that did not fire.

- [ ] **Step 5: Check the file size a real state costs**

```bash
conda run -n climt python -c "
import sys, os, tempfile, sympl, climt
sys.path.insert(0, 'docs/modelling-tour/_tour'); import states
sympl.set_backend(climt.UnytBackend())
comps = [climt.CorkLongwaveRadiation(optics='correlated_k', table='earth_low_res_lw'),
         climt.SlabSurface(), climt.SimpleBoundaryLayer(surface_fluxes='bulk'),
         climt.DryConvectiveAdjustment()]
s = climt.get_default_state(comps, grid_state=climt.get_grid(nx=1, ny=1, nz=28))
p = os.path.join(tempfile.mkdtemp(), 'probe.npz')
states.save(p, s, dict(table='earth_low_res_lw', nz=28))
print(os.path.getsize(p), 'bytes')"
```

Expected: a few kB. The spec's budget is "four orders of magnitude below the table already committed" — the 56-band table is 5.6 MB, so anything under ~50 kB is comfortably within it. If a state comes out over 100 kB, find out which quantity is large (`longwave_optical_thickness_due_to_aerosol` is `(nband, nz, 1, 1)`) and add it to `SKIPPED` only if `get_default_state` reproduces it exactly — check by round-tripping before and after.

- [ ] **Step 6: Commit**

```bash
git add docs/modelling-tour/_tour/states.py tests/test_modelling_tour.py
git commit -m "feat(tour): _tour/states.py, an equilibrium that travels with its provenance

Save a state to .npz; load it by rebuilding get_default_state for the page's
own components and overwriting arrays, checking dims, shape and units per
quantity. A renamed or re-dimensioned quantity fails by name instead of
producing a plausible non-equilibrium.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 8: The measurement campaign — the spec's Task 0

**No `.qmd` in this tranche is written until this task is done and its log is filled in.** Every number the six pages quote descends from here. The spec is explicit that if the measurements come back bad, the levers are — in order — a thinner slab, a faster perturbation, then fewer levels. **Degrading the physics to make a cell finish is not a lever.** A run that takes minutes while the lecture continues is the accepted trade.

Priors going in, from tranche 1's experiment log and from the measurements at the head of this plan: gray radiative equilibrium converges in ~300 steps at `dt = 12 h`; non-grey needs ~1950 at `dt = 6 h`; `dt = 24 h` is unstable for gray. Emanuel is the unknown — the scheme is conventionally run at 10–20 minutes, and if that holds, page 12's perturbation is thousands of steps.

**Every stack with a boundary layer carries the wind relaxation from Task 4 Step 8.** These measurements are of the configuration the pages actually ship, and a column without a momentum source is not that configuration — its sensible heat flux is roughly a third of the relaxed column's. `column()` below adds it wherever `with_boundary_layer` is set.

**Files:**
- Create: `scripts/experiments/tour_rce_measurements.py`
- Modify: this plan (the log section at the end of this task)

**Interfaces:**
- Consumes: `_tour/stepping.py`, `_tour/budgets.py` (Tasks 4, 6).
- Produces: no code other pages import. It produces **numbers**, recorded in this plan, which Tasks 9–15 then use.

- [ ] **Step 1: Write the measurement script**

Create `scripts/experiments/tour_rce_measurements.py`:

```python
"""Measure everything the Modelling Tour's RCE pages need to quote.

The spec's Task 0: no page is written until these five are known. Run under
NUMBA_DISABLE_JIT=1 for anything cost-related -- Pyodide has no numba, so a
JIT-compiled timing is off by ~20x on the 14-band radiation and worthless as a
browser estimate.

    NUMBA_DISABLE_JIT=1 conda run -n climt python \\
        scripts/experiments/tour_rce_measurements.py all

Individual measurements:

    ... tour_rce_measurements.py stability     # 1: largest stable timestep
    ... tour_rce_measurements.py convergence   # 2: steps to equilibrium
    ... tour_rce_measurements.py cost          # 3: per-step cost per page
    ... tour_rce_measurements.py independence  # 4: start-independence
    ... tour_rce_measurements.py supersat      # 5: what condensation removes

Output goes to stdout as a table and to debug_data/tour_rce_measurements.json.
"""
import json
import os
import sys
import time

import numpy as np
import sympl

import climt

sys.path.insert(0, os.path.join("docs", "modelling-tour", "_tour"))
import budgets      # noqa: E402
import stepping     # noqa: E402

OUT_DIR = "debug_data"
SOLAR = 240.0
NZ = 28
SLAB_DEPTH = 2.0
WIND = 5.0           # m/s, the wind the relaxation holds the column at
Z0 = 1e-3            # surface roughness length, m
GRAY = "single_band_gray_lw"
BAND14 = "earth_low_res_lw"


def column(table, with_boundary_layer=False, with_dry_adjustment=False,
           moist=False, slab_depth=SLAB_DEPTH, nz=NZ, co2_ppm=None,
           wind=WIND, roughness_length=Z0):
    """Build one page's exact stack. Returns (tendencies, steppers, state)."""
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k", table=table)
    surface = climt.SlabSurface()
    tendencies = [longwave, surface]
    steppers = []
    if with_boundary_layer:
        steppers.append(climt.SimpleBoundaryLayer(
            surface_fluxes="bulk", roughness_length=roughness_length))
    if with_dry_adjustment:
        steppers.append(climt.DryConvectiveAdjustment())
    if moist:
        tendencies.append(climt.EmanuelConvectionPython())
        steppers.append(climt.GridScaleCondensation())

    state = climt.get_default_state(tendencies + steppers,
                                    grid_state=climt.get_grid(nx=1, ny=1, nz=nz))
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    if co2_ppm is not None:
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = co2_ppm * 1e-6
    if moist:
        # A saturated surface, so the latent flux is on. Page 9's knob.
        state["surface_specific_humidity"].values[:] = 0.015
    if with_boundary_layer:
        # Without this the column spins itself down to dead calm within days
        # and its surface fluxes are those of a windless planet. See
        # _tour/stepping.wind_relaxation.
        tendencies.append(stepping.wind_relaxation(state, wind))
    return tendencies, steppers, state


STACKS = {
    "07-gray":   dict(table=GRAY),
    "07-14band": dict(table=BAND14),
    "08-bl":     dict(table=GRAY, with_boundary_layer=True),
    "09-moist-bl": dict(table=BAND14, with_boundary_layer=True, moist=False),
    "11-dry-rce": dict(table=BAND14, with_boundary_layer=True,
                       with_dry_adjustment=True),
    "12-moist-rce": dict(table=BAND14, with_boundary_layer=True,
                         with_dry_adjustment=True, moist=True),
}


# ------------------------------------------------------------------ 1: stability

def stability(hours=(0.25, 0.5, 1, 2, 3, 6, 12, 24), n_steps=200):
    """Largest timestep each stack survives n_steps without going non-finite.

    "Survives" is: no exception (Task 2's guard raises on non-finite fluxes),
    no NaN anywhere, and a surface temperature still in [150, 400] K.
    """
    results = {}
    for name, kwargs in STACKS.items():
        largest = None
        for dt_hours in hours:
            tendencies, steppers, state = column(**kwargs)
            timestep = climt.UnytTimeDelta(hours=dt_hours)
            try:
                stepping.integrate(tendencies, steppers, state, timestep,
                                   n_steps)
                T = state["air_temperature"].values
                surface = float(state["surface_temperature"].values.ravel()[0])
                ok = (np.all(np.isfinite(T)) and 150.0 < surface < 400.0)
            except Exception as error:                   # noqa: BLE001
                ok = False
                print(f"    {name} dt={dt_hours}h -> {type(error).__name__}: "
                      f"{str(error)[:70]}")
            if ok:
                largest = dt_hours
            else:
                break
        results[name] = largest
        print(f"  {name:14s} largest stable dt = {largest} h "
              f"({n_steps} steps)")
    return results


# ---------------------------------------------------------------- 2: convergence

def steps_to_equilibrium(tendencies, steppers, state, timestep,
                         threshold=0.5, max_steps=6000, check_every=25):
    """Steps until |TOA imbalance| stays under `threshold` W/m^2.

    Returns (n_steps, final_toa, final_surface_temperature) or
    (None, ...) if it never gets there.
    """
    for step in range(0, max_steps, check_every):
        stepping.integrate(tendencies, steppers, state, timestep, check_every)
        imbalance = budgets.toa_imbalance(state)
        if abs(imbalance) < threshold:
            return (step + check_every, imbalance,
                    float(state["surface_temperature"].values.ravel()[0]))
    return (None, budgets.toa_imbalance(state),
            float(state["surface_temperature"].values.ravel()[0]))


def convergence(dt_hours=None):
    """Steps to equilibrium from a cold start, per stack, at its stable dt.

    `dt_hours` overrides the per-stack choice; otherwise each stack uses the
    largest dt `stability()` cleared, capped at 12 h (beyond which the
    snapshots on page 7 are too coarse to show an approach).
    """
    stable = stability()
    results = {}
    for name, kwargs in STACKS.items():
        dt = dt_hours or min(stable[name] or 1, 12)
        tendencies, steppers, state = column(**kwargs)
        started = time.time()
        n_steps, imbalance, surface = steps_to_equilibrium(
            tendencies, steppers, state, climt.UnytTimeDelta(hours=dt))
        elapsed = time.time() - started
        results[name] = dict(dt_hours=dt, n_steps=n_steps,
                             toa_imbalance=imbalance,
                             surface_temperature=surface,
                             wall_seconds=elapsed)
        days = None if n_steps is None else n_steps * dt / 24.0
        print(f"  {name:14s} dt={dt:5.2f}h  {str(n_steps):>6s} steps "
              f"({'' if days is None else f'{days:.0f} d'})  "
              f"TOA {imbalance:+7.3f}  Tsurf {surface:6.2f} K  "
              f"[{elapsed:.1f} s native]")
    return results


def co2_doubling_response(name="11-dry-rce", dt_hours=1.0, threshold=0.5):
    """From equilibrium, double CO2 and re-converge. The measured sensitivity."""
    kwargs = STACKS[name]
    tendencies, steppers, state = column(**kwargs)
    steps_to_equilibrium(tendencies, steppers, state,
                         climt.UnytTimeDelta(hours=dt_hours),
                         threshold=threshold)
    before = float(state["surface_temperature"].values.ravel()[0])
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] *= 2.0
    n_steps, imbalance, after = steps_to_equilibrium(
        tendencies, steppers, state, climt.UnytTimeDelta(hours=dt_hours),
        threshold=threshold)
    print(f"  {name}: 2xCO2  {before:.2f} -> {after:.2f} K "
          f"(dT = {after - before:+.2f} K) in {n_steps} steps at "
          f"dt = {dt_hours} h")
    return dict(before=before, after=after, warming=after - before,
                n_steps=n_steps, toa_imbalance=imbalance)


# ---------------------------------------------------------------------- 3: cost

def cost(n=30, warmup=5, dt_hours=1.0):
    """Per-step wall clock for each page's exact component list."""
    results = {}
    for name, kwargs in STACKS.items():
        tendencies, steppers, state = column(**kwargs)
        timestep = climt.UnytTimeDelta(hours=dt_hours)
        stepping.integrate(tendencies, steppers, state, timestep, warmup)
        started = time.time()
        stepping.integrate(tendencies, steppers, state, timestep, n)
        per_step = (time.time() - started) / n * 1000.0
        results[name] = per_step
        print(f"  {name:14s} {per_step:7.2f} ms/step native "
              f"(~{per_step * 3.5 / 1000:.3f} s/step browser)")
    return results


# -------------------------------------------------------------- 4: independence

def independence(name="11-dry-rce", dt_hours=1.0, threshold=0.5,
                 tolerance=0.5):
    """Two different initial conditions must reach the same equilibrium.

    Load-bearing: shipping an equilibrium state is only honest if the
    equilibrium does not depend on how you got there, and page 12 says so in
    prose. Start A is climt's default profile; start B is 40 K warmer at the
    surface with a steeper lapse rate.
    """
    kwargs = STACKS[name]
    finals = {}
    for label in ("default", "warm-steep"):
        tendencies, steppers, state = column(**kwargs)
        # Note: the wind is NOT part of the initial condition to vary. It is
        # relaxed to a fixed target, so it is a boundary condition here, the
        # same in both runs. Varying it would test something else.
        if label == "warm-steep":
            p = state["air_pressure"].values[:, 0, 0]
            surface_pressure = float(
                state["surface_air_pressure"].values.ravel()[0])
            T = 328.0 - 9.0e-3 * 7000.0 * np.log(surface_pressure / p)
            state["air_temperature"].values[:, 0, 0] = np.maximum(T, 170.0)
            state["surface_temperature"].values[:] = 330.0
        n_steps, imbalance, surface = steps_to_equilibrium(
            tendencies, steppers, state, climt.UnytTimeDelta(hours=dt_hours),
            threshold=threshold, max_steps=12000)
        finals[label] = dict(
            n_steps=n_steps, toa_imbalance=imbalance,
            surface_temperature=surface,
            profile=state["air_temperature"].values[:, 0, 0].copy())
        print(f"  {label:11s} {str(n_steps):>6s} steps  "
              f"Tsurf {surface:6.2f} K  TOA {imbalance:+7.3f}")

    difference = np.abs(finals["default"]["profile"]
                        - finals["warm-steep"]["profile"])
    surface_difference = abs(finals["default"]["surface_temperature"]
                             - finals["warm-steep"]["surface_temperature"])
    print(f"  max |dT| through the column: {difference.max():.3f} K; "
          f"surface: {surface_difference:.3f} K")
    verdict = ("START-INDEPENDENT" if difference.max() < tolerance
               else "PATH-DEPENDENT — DO NOT SHIP A STATE")
    print(f"  -> {verdict} (tolerance {tolerance} K)")
    return dict(max_profile_difference=float(difference.max()),
                surface_difference=float(surface_difference),
                verdict=verdict,
                default=finals["default"]["surface_temperature"],
                warm_steep=finals["warm-steep"]["surface_temperature"])


# ------------------------------------------------------------------ 5: supersat

def supersaturation(dt_hours=1.0, n_steps=200):
    """How much supersaturation GridScaleCondensation removes on page 9.

    Page 9's stack is the boundary layer over a saturated surface with 14-band
    radiation. Without a moisture sink the column supersaturates: the BL mixes
    moist air into colder air aloft. This runs the stack twice -- once with
    condensation, once without -- and reports the peak relative humidity each
    reaches, and the precipitation the sink produces.
    """
    from sympl import get_constant

    Rd = float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
    Rv = float(get_constant("gas_constant_of_vapor_phase", "J/kg/degK"))
    epsilon = Rd / Rv

    def saturation_specific_humidity(T, p):
        # Bolton (1980), the same form _tour/soundings.py uses.
        e_sat = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
        return epsilon * e_sat / np.maximum(p - (1.0 - epsilon) * e_sat, 1.0)

    results = {}
    for label, with_sink in (("with condensation", True),
                             ("without condensation", False)):
        longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                               table=BAND14)
        surface = climt.SlabSurface()
        boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                                   roughness_length=Z0)
        steppers = [boundary_layer]
        if with_sink:
            steppers.append(climt.GridScaleCondensation())
        state = climt.get_default_state([longwave, surface] + steppers,
                                        grid_state=climt.get_grid(nx=1, ny=1,
                                                                  nz=NZ))
        state["ocean_mixed_layer_thickness"].values[:] = SLAB_DEPTH
        state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
        state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
        state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
        state["surface_specific_humidity"].values[:] = 0.015
        relaxation = stepping.wind_relaxation(state, WIND)

        peak = 0.0
        timestep = climt.UnytTimeDelta(hours=dt_hours)
        for _ in range(n_steps):
            stepping.integrate([longwave, surface, relaxation], steppers,
                               state, timestep, 1)
            T = state["air_temperature"].values[:, 0, 0]
            q = state["specific_humidity"].values[:, 0, 0]
            p = state["air_pressure"].values[:, 0, 0]
            relative_humidity = q / saturation_specific_humidity(T, p)
            peak = max(peak, float(np.max(relative_humidity)))
        precipitation = budgets.precipitation_rate(state, timestep)
        results[label] = dict(peak_relative_humidity=peak,
                              precipitation_mm_day=precipitation)
        print(f"  {label:22s} peak RH {peak * 100:6.1f}%   "
              f"precip {precipitation:6.3f} mm/day")
    return results


# -------------------------------------------------------------------- driver

MEASUREMENTS = {
    "stability": stability,
    "convergence": convergence,
    "cost": cost,
    "independence": independence,
    "supersat": supersaturation,
    "co2": co2_doubling_response,
}


def main():
    sympl.set_backend(climt.UnytBackend())
    os.makedirs(OUT_DIR, exist_ok=True)
    if os.environ.get("NUMBA_DISABLE_JIT") != "1":
        print("WARNING: numba JIT is ON. Cost numbers will not reflect the "
              "browser. Re-run with NUMBA_DISABLE_JIT=1.\n")

    requested = sys.argv[1:] or ["all"]
    names = list(MEASUREMENTS) if requested == ["all"] else requested
    collected = {}
    for name in names:
        print(f"\n== {name} ==")
        collected[name] = MEASUREMENTS[name]()

    path = os.path.join(OUT_DIR, "tour_rce_measurements.json")
    with open(path, "w") as handle:
        json.dump(collected, handle, indent=2, default=float)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run measurement 3 first — it is the cheapest and gates the rest**

```bash
NUMBA_DISABLE_JIT=1 conda run -n climt python \
  scripts/experiments/tour_rce_measurements.py cost
```

Expected: six lines, roughly matching the table at the head of this plan (2–3 ms for the gray stacks, ~50 ms for the 14-band ones). If the 14-band stacks come out under 10 ms, `NUMBA_DISABLE_JIT=1` did not take effect and every number that follows is wrong — check the warning line printed at startup.

- [ ] **Step 3: Run measurement 1 — stability**

```bash
NUMBA_DISABLE_JIT=1 conda run -n climt python \
  scripts/experiments/tour_rce_measurements.py stability
```

Expected: a largest-stable-`dt` per stack. **The number that decides page 12 is `12-moist-rce`.** If it comes back at 0.25 h, page 12's perturbation at ~70 simulated days is ~6700 steps ≈ 21 minutes in-browser, and the levers open — in the spec's order.

- [ ] **Step 4: Run measurement 2 — convergence, and the CO₂ response**

```bash
NUMBA_DISABLE_JIT=1 conda run -n climt python \
  scripts/experiments/tour_rce_measurements.py convergence
NUMBA_DISABLE_JIT=1 conda run -n climt python \
  scripts/experiments/tour_rce_measurements.py co2
```

Expected: steps to equilibrium per stack, and a 2×CO₂ surface warming for the dry stack. `convergence` runs `stability` first, so it is the long one — allow tens of minutes.

- [ ] **Step 5: Run measurement 4 — start-independence. This one is a gate.**

```bash
NUMBA_DISABLE_JIT=1 conda run -n climt python \
  scripts/experiments/tour_rce_measurements.py independence
```

Expected: `START-INDEPENDENT`, with a max profile difference well under 0.5 K.

**If it prints `PATH-DEPENDENT`, stop.** Do not proceed to Task 9. Shipping an equilibrium state is only honest if the equilibrium does not depend on the path taken to it, and page 12 says so in prose. Investigate before writing anything: the usual causes are a convergence threshold too loose to have actually converged either run (tighten `threshold` and re-run) or a genuine multiple equilibrium (which is a finding worth a page of its own, and a reason to reconsider the tranche's shape with the spec's author).

Run it for `12-moist-rce` as well as the default `11-dry-rce` — the moist column is the one with a plausible physical route to path dependence:

```bash
NUMBA_DISABLE_JIT=1 conda run -n climt python -c "
import sys; sys.argv = ['x']
sys.path.insert(0, 'scripts/experiments')
import sympl, climt, tour_rce_measurements as m
sympl.set_backend(climt.UnytBackend())
m.independence(name='12-moist-rce', dt_hours=0.25)"
```

- [ ] **Step 6: Run measurement 5 — the supersaturation page 09 quotes**

```bash
NUMBA_DISABLE_JIT=1 conda run -n climt python \
  scripts/experiments/tour_rce_measurements.py supersat
```

Expected: two peak relative humidities. The one without condensation should exceed 100%; the difference is the number page 09 quotes for what `GridScaleCondensation` removes.

- [ ] **Step 7: Fill in the log below, in this file**

Paste the actual output into the `### Log — Task 0 measurements` section at the end of this task. Every one of the six pages will cite it. Then, for each page, write the one line it needs:

| Page | The number it needs from here |
|---|---|
| 07 | `dt`, step count and wall clock for the gray and 14-band columns; the top level's distance from `Te/2^(1/4)` |
| 08 | the surface–air discontinuity, radiative-equilibrium value and radiative-plus-turbulent value |
| 09 | the peak relative humidity with and without condensation; the Bowen ratio over a saturated surface |
| 10 | nothing — no time loop |
| 11 | steps and `dt` for the shipped state; the 2×CO₂ warming and how many steps it took |
| 12 | Emanuel's largest stable `dt`; steps for the moist state; the moist 2×CO₂ warming |

- [ ] **Step 8: Decide, and record the decision**

If page 12's browser cell exceeds ~10 minutes at the measured `dt`, apply the levers in order and record which was used and why:

1. **Thinner slab** — 1 m instead of 2 m halves the response time. Depth is already page 07's knob, and its lesson *is* response time, so this costs the tranche nothing.
2. **A faster perturbation** — a surface-temperature nudge instead of a CO₂ doubling reaches a new equilibrium sooner.
3. **Fewer levels** — `nz = 20` instead of 28. Costs the flat skin layer page 07 wants; take it last.

**Not a lever: degrading the physics.** No shortening the spin-up below convergence, no gray table where the page says 14-band, no dropping condensation.

- [ ] **Step 9: Commit**

```bash
git add scripts/experiments/tour_rce_measurements.py \
        docs/superpowers/plans/2026-08-20-modelling-tour-rce.md
git commit -m "docs(plan): Task 0 measured — the numbers the six pages will quote

Adds scripts/experiments/tour_rce_measurements.py and records its output in
the plan: stable timesteps, steps to equilibrium, per-step cost per page,
start-independence, and the supersaturation page 9 quotes.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

### Log — Task 0 measurements

**Machine:** Apple Silicon, 5 performance + 6 efficiency cores (`hw.perflevel0.logicalcpu`=5, `hw.perflevel1.logicalcpu`=6). **Date:** 2026-09-03. Same reference machine as tranche 1. Cost/stability/convergence measured with `NUMBA_DISABLE_JIT=1` (the browser proxy); the moist-stack numbers (see below) measured with numba **on** — see the note there.

**1 — Per-step cost (`cost`, numba off, the browser proxy):**

```
  07-gray         2.61 ms/step native (~0.009 s/step browser)
  07-14band      49.41 ms/step native (~0.173 s/step browser)
  08-bl           3.70 ms/step native (~0.013 s/step browser)
  09-moist-bl    51.37 ms/step native (~0.180 s/step browser)
  11-dry-rce     53.92 ms/step native (~0.189 s/step browser)
  12-moist-rce   56.78 ms/step native (~0.199 s/step browser)
```
Radiation dominates ~18×; every non-radiative component this tranche adds is nearly free. Cost is set by dt and step count.

**2 — Largest stable dt (`stability`, 200-step blow-up test):**

```
  07-gray        24 h      08-bl         24 h      11-dry-rce    24 h
  07-14band      24 h      09-moist-bl    6 h      12-moist-rce   6 h
```
**CAVEAT — the 200-step test is unreliable for the moist stacks.** It reports `12-moist-rce` "stable at 6 h", but over a real (thousands-of-steps) convergence run the moist column goes non-finite at 6 h (a `condensibles.py` log-of-negative precedes the blow-up). The repo's own Emanuel examples all run this column at **5 min** (`examples/gmd_radiative_convective_python_emanuel.ipynb`, `gmd_radiative_convective_unyt.py`, `gmd_radiative_convective.py`: `UnytTimeDelta(minutes=5)`); `column_code_with_slab.py` uses 10 min. Page 12's dt is **5 min**, not 6 h.

**3 — Steps to equilibrium (`convergence`, cold start, numba off; `|TOA|<0.5 W/m²`):**

```
  stack        dt     n_steps  (days)  TOA      Tsurf     dT(sfc-air)  Ttop     Ttop/Te   wall
  07-gray      12 h    600      300 d   +0.419   335.43 K   +15.16 K    214.96   0.8428    1.5 s
  07-14band    12 h   2300     1150 d   -0.499   267.85 K    +7.08 K    113.08   0.4433  108.6 s
  08-bl        12 h    625      312 d   +0.458   332.40 K    +6.40 K    214.95   0.8427    2.3 s
  09-moist-bl   6 h   2450      612 d   -0.177   289.21 K    +4.24 K    124.12   0.4866  122.3 s
  11-dry-rce   12 h   1650      825 d   -0.495   266.60 K    +3.88 K    119.81   0.4697   85.0 s
```
(`Te = (240/σ)^¼ = 255.06 K`; skin temp `Te/2^¼ = 214.48 K` — the gray tops sit right at it.) The `convergence` run CRASHED on the 6th stack (`12-moist-rce` at 6 h — the blow-up above), so its shard holds these five; the moist stack was measured separately at 5 min:

**3b — Moist stack at dt = 5 min, numba ON (`scripts/experiments/tour_rce_moist_probe.py`):** *(TOA-only gate: all three runs stopped ~7 K warm of the true equilibrium. Superseded by §7.)* **Superseded 2026-09-27 (Task 9 log):** every moist number here was made with `surface_specific_humidity` fixed at 0.015 kg kg⁻¹, which is 150–250 % of saturation at these surface temperatures. The shipped moist states are now made over a surface held saturated at its own temperature (`SurfaceHumidity(1.0)`).

```
  default start   n_steps=123400 (428 d)  TOA=+0.262  Tsurf=286.67 K  [856 s native, numba on]
  2xCO2           n_steps= 20000 ( 69 d)  TOA=+0.487  Tsurf=288.85 K  [136 s]
  warm-steep      n_steps=118800 (412 d)  TOA=-0.479  Tsurf=286.74 K  [802 s]
```
Numba is ON here on purpose: step count is a physics quantity (numba-independent), and numba makes the run ~10× faster (~7 ms/step vs the 56.8 ms/step browser cost recorded in §1). The offline state generator (Task 9) runs natively, so it gets numba; only the browser cell is bound by §1's cost.

**4 — Start-independence (`independence`; two starts: climt default vs +40 K/steeper):**

```
  11-dry-rce  @12 h: default 1650 st (266.60 K), warm-steep 1525 st (266.52 K)
              -> surface |dT| 0.078 K; max column |dT| 2.01 K
  12-moist-rce @5 min: default 123400 st (286.67 K), warm-steep 118800 st (286.74 K)
              -> surface |dT| 0.064 K; max column |dT| 3.84 K
```
**Reading — the spec author's decision (this is the plan's Step 5 gate, resolved 2026-09-03).** The script prints `PATH-DEPENDENT` because both stacks exceed its 0.5 K whole-column tolerance, but **the spec author has ruled this is not a path dependence and these are accepted as converged runs.** The evidence supports that: in both stacks the **surface** is reproducible to <0.08 K, both starts met the `|TOA|<0.5` gate, and the residual spread is confined to the slow-relaxing upper column (a whole-column energy gate does not pin the stratosphere, whose radiative relaxation time is far longer than the surface's — a benign feature, not a second attractor). So **both equilibria are shippable as-is; Task 9 does not need a stricter convergence gate.** The `PATH-DEPENDENT` string is a tolerance artifact of the measurement script, not a finding. Page 12 may claim start-independence. **But see §7:** these two moist runs were stopped by the TOA-only gate, both about 7 K warm of the shipped state. The claim was re-measured on converged runs and holds. The earlier `independence-moist` = 21 K was separately pure under-resolution (dt=0.25 h × 12000 steps = 125 d, vs the 428 d needed) and is superseded by the 5-min run above. All of §4's moist runs used the fixed 0.015 kg kg⁻¹ surface, superseded 2026-09-27 (Task 9 log).

**5 — 2×CO₂ response (`co2`, from equilibrium, double CO₂, re-converge):**

```
  11-dry-rce  @12 h:  266.60 -> 267.51 K  (+0.913 K) in   175 steps  TOA +0.409
  12-moist-rce @5 min: 286.67 -> 288.85 K  (+2.179 K) in 20000 steps  TOA +0.487
```
**Superseded by §7: both warmings were measured from runs the TOA-only gate stopped early. Do not quote them.** The moist one is also a fixed-0.015 surface; see the Task 9 log (2026-09-27). The moist response is ~2.4× the dry — water-vapour feedback via Emanuel + condensation. (Each pair of end-states sits on opposite sides of the 0.5 W/m² gate, so both warmings carry a few-tenths-K convergence-slop; fine for a teaching page, quote as ≈.)

**6 — Supersaturation page 09 quotes (`supersat`, page-09 stack, dt=1 h, 200 steps):**

```
  with condensation      peak RH 100.0%   precip 2.693 mm/day   SH 48.54  LH 75.42  Bowen 0.644
  without condensation   peak RH 456.7%   precip 0.000 mm/day   SH 79.50  LH 32.33  Bowen 2.459
```
Without a moisture sink the boundary layer drives the column to **457 %** relative humidity; `GridScaleCondensation` holds it at 100 % and rains out 2.69 mm/day. Bowen ratio over the saturated surface is 0.64. *(This "saturated surface" was the fixed 0.015 kg kg⁻¹. Page 09 re-measured it with `SurfaceHumidity`; see Task 12's log, "Supersaturation (replaces Task 0 §6)".)*

**7 — Re-measured against the shipped states (2026-09-26).** §3–§5 stopped every run at the first time |TOA| < 0.5 W m⁻². The shipped equilibria (Task 9) also require the surface temperature to stop moving, and they sit elsewhere: dry 266.48 K rather than 266.60 K, moist **279.97 K rather than 286.67 K**. The moist warm-start run below shows why. At step 120 000 it crossed TOA +0.49 W m⁻² at 286.1 K while still cooling 0.2 K per window. The TOA-only gate would have stopped it there, which is Task 0's 286.67 K after 123 400 steps. The TOA imbalance oscillates through zero on its way down, and the old gate stopped on a crossing. Every moist number in §3b–§5 was taken at such a crossing.

> **Superseded 2026-09-27: the moist half of this section.** The moist rows below, and the moist reading and table entries after them, describe states made over the fixed 0.015 kg kg⁻¹ surface (246 % RH at 279.97 K). Both moist states were regenerated over a saturated surface with a moist-specific trend gate. The column now settles at 286.0 K, TOA ≈ −1.1 W m⁻² (not +0.3), and its 2×CO₂ warming is +2.24 K (not +1.82 K). See the Task 9 log, 2026-09-27. The dry rows are unaffected.

Re-measured with `scripts/experiments/tour_rce_shipped_remeasure.py`, which takes its configuration, gate and loop from `scripts/generate_tour_equilibria.py`. The moist 2×CO₂ came from `generate_tour_equilibria.py --moist-2xco2`, which ships its result (Task 9 log). Linux, 4 cores, numba on, `NUMBA_NUM_THREADS=1`:

```
  dry 2xCO2 (from rce_dry_equilibrium.npz, dt 12 h):
      266.478 -> 267.678 K  (+1.200 K)  strict gate at 2350 steps; +1.197 K by step 1000 (TOA -0.028)
  moist 2xCO2 (from rce_moist_equilibrium.npz, dt 5 min):
      279.965 -> 281.613 K  (+1.647 K)  73 600 steps (256 d), TOA +0.500
  moist, both shipped states stepped on 60 000 steps (settle-moist; mean over the last 30 000):
      base   settles at 279.862 K, TOA +0.306 (std 0.008)   (file: 279.965 K)
      2xCO2  settles at 281.683 K, TOA +0.347 (std 0.010)   (file: 281.613 K; still +0.016 K/100 d)
      settled 2xCO2 warming: +1.82 K
  dry start-independence (warm-steep start, strict gate):
      3350 steps, Tsurf 266.478 vs 266.478 K (7e-6 K); column max |dT| 0.024 K at 3 hPa
  moist start-independence (warm-steep start, strict gate):
      195 600 steps, Tsurf 279.967 vs 279.965 K (0.002 K); below 200 hPa 0.008 K; max 1.72 K at 3 hPa
```

**Reading.**
- **Start-independence holds, and now on converged runs.** Two starts ~50 K apart agree at the surface to 0.002 K (moist) and 7 × 10⁻⁶ K (dry), and through the troposphere to under 0.01 K. The only visible spread is at the top model level (3 hPa), where the radiative relaxation time is longest. §4's reading was right about that, just on the wrong runs.
- **Dry sensitivity: +1.20 K, not +0.91 K.** Page 11 runs it live for 1000 steps (Task 14 Step 1).
- **Moist sensitivity: +1.82 K**, about 1.5× the dry response. It is not §5's +2.18 K, and not the +1.647 K between the shipped files. **The moist column never reaches TOA = 0.** `EmanuelConvectionPython` is not fully energy-conserving (a known property of the scheme, confirmed with the spec author 2026-09-26), so the column settles with a steady TOA imbalance of about +0.3 W m⁻². The cold start with the gate at 0.1 W m⁻² shows it: that run swept through ±0.1 near step 205 000 and parked at +0.30, 279.863 K. So the ±0.5 gate stops a spin-up wherever TOA first enters the band on its way to +0.3, not where the column settles. The shipped base stopped 0.10 K warm of its settled state and the 2×CO₂ state 0.07 K cool, and both files are kept as shipped. The sensitivity is the difference between the *settled* states, measured by `settle-moist`. (An earlier draft of this section quoted ≈ 2.1 K from a Gregory regression to TOA = 0. That was wrong: it assumed an equilibrium this column does not have.)
- **A tighter gate does not help.** Neither 0.1 nor any threshold below +0.3 is a convergence criterion for this column. A threshold can only be passed transiently, while TOA sweeps through it. If the moist states are ever regenerated, stop on a flat trend in TOA (for example a 30-day running mean), not on its magnitude.
- **Why the moist gate is loose and the dry one is not.** The drift window is 1000 *steps*: 500 days at 12 h, but only 3.5 days at 5 min. For the moist column the drift gate barely constrains anything, and |TOA| < 0.5 does all the work. The dry states are unaffected (both imbalances < 0.15 W m⁻²).

**The number each page takes from here (Step 7):**

| Page | Number(s) |
|---|---|
| 07 | gray dt 12 h / 600 steps / 1.5 s, top 214.96 K = skin `Te/2^¼`=214.48 K; 14-band dt 12 h / 2300 steps / 108.6 s, top 113.08 K (0.443 `Te`) |
| 08 | surface–air discontinuity: gray radiative-only +15.16 K → gray radiative+turbulent +6.40 K |
| 09 | peak RH 100 % (with) vs 457 % (without) condensation; Bowen 0.64 over a saturated surface |
| 10 | — (no time loop) |
| 11 | dry state dt 12 h / **3500 steps / 266.48 K** (shipped); 2×CO₂ **+1.20 K**, run live for 1000 steps (§7) |
| 12 | *Superseded 2026-09-27 (Task 9 log): moist state 295 800 steps / 286.01 K, TOA ≈ −1.1, 2×CO₂ +2.24 K, start-independent to ≈ 0.05 K at the surface.* Was: Emanuel dt 5 min; moist state **200 550 steps (≈ 696 d) / 279.97 K** (shipped); 2×CO₂ **+1.82 K** between the settled states (not the +1.65 K between the files), run offline (§7); settles at TOA ≈ +0.3, not 0 (Emanuel is not fully conservative); start-independent to 0.002 K at the surface (§7) |

**Decisions taken from these numbers:**

- **Slab depth for pages 11 and 12:** **2 m** (unchanged, the shipped value). No lever applied — see below.
- **`dt` for each page:** 07 gray **12 h**, 07 14-band **12 h**; 08 **12 h**; 09 illustrative loop **1 h** (equilibrium at 6 h); 11 dry **12 h**; 12 moist **5 min** (the Emanuel convention — **not** the 6 h the 200-step stability test wrongly cleared).
- **Step counts for the two shipped equilibria:** dry (page 11) **3500 steps @ 12 h** (1750 d); moist (page 12) **200 550 steps @ 5 min** (≈ 696 d), superseded 2026-09-27 by **295 800** (≈ 1027 d, saturated surface, moist trend gate). These are the shipped values (Task 9 log). This line originally said ≈1650 and ≈123 400, the TOA-only-gate numbers §7 supersedes.
- **Levers applied, if any, and why:** **None.** Page 12 is never spun up live in the browser — 200 550 steps × 56.8 ms ≈ **3.2 h** in-browser is impossible (this read 123 400 steps ≈ 117 min before §7) — so it **loads a precomputed equilibrium** (Task 9, already the plan's design), generated natively-with-numba in ~14 min (measured, §3b). Because the browser never runs the spin-up, the browser-cell-time levers (thinner slab / faster perturbation / fewer levels) do not apply, and the shipped-state config stays 2 m / 5 min / nz = 28. ~~Task 9 ships both states directly from this convergence, so no stricter gate or re-run is needed.~~ **Superseded:** Task 9 added the surface-stationarity gate and re-ran both states (Task 9 log), and §7 re-measured everything quoted from them. Page 12's 2×CO₂ re-equilibration is 73 600 steps, about 3.9 h in the browser, so it ships as a third state, `rce_moist_2xco2_equilibrium.npz`. (Since 2026-09-27 it is 180 550 steps, ≈ 9.5 h.)

---

## Task 9: Ship the two equilibrium states, and guard them with a residual test

Pages 11 and 12 load an equilibrium. This task produces the two files, the generator that made them, and the CI test that notices when the physics moves under them.

**The guard is a residual test, deliberately not dependency hashing.** `build_experiments.py --check` hashes `cork/**/*.py` and re-runs anything downstream when the hash moves. That machinery is what made a class-attribute one-liner in `cork/lw/component.py` re-run 57 600 five-minute steps twice, and a content hash cannot tell a no-op from a physics change. Instead: load the shipped state, step it ten steps under the page's own component list, and assert the TOA imbalance stays under threshold and the surface temperature drifts less than ~0.05 K. Seconds in CI, and it fails exactly when the physics genuinely moved — at which point the generator is re-run deliberately.

**Files:**
- Create: `scripts/generate_tour_equilibria.py`
- Create: `docs/modelling-tour/_data/rce_dry_equilibrium.npz`
- Create: `docs/modelling-tour/_data/rce_moist_equilibrium.npz`
- Modify: `docs/modelling-tour/_data/README.md`
- Modify: `tests/test_modelling_tour.py`

**Interfaces:**
- Consumes: `_tour/stepping.py`, `_tour/budgets.py`, `_tour/states.py`.
- Produces:
  - `generate_tour_equilibria.dry_components() -> list`, `moist_components() -> list` — **the single definition of pages 11 and 12's component lists.** The pages import nothing from here (it is a script, not a `_tour` module), but the residual test does, so the test and the generator can never disagree about what the state was made of.
  - Two `.npz` files loadable by `states.load`.
  - CLI: `python scripts/generate_tour_equilibria.py [--dry] [--moist] [--out DIR]`, whose **defaults reproduce exactly what shipped** — the lesson `generate_tour_spectrum_table.py` taught by not doing so.

- [x] **Step 1: Write the residual test first** — done, verbatim from below; appended after the existing `states` tests.

Add to `tests/test_modelling_tour.py`:

```python
# --------------------------------------------------- the shipped equilibria
#
# These states are not regenerated by build_experiments.py's dependency
# hashing. A content hash over cork/**/*.py cannot tell a no-op from a physics
# change, and tranche 1 paid for that twice. The guard is instead: load the
# shipped state, step it a little under the page's own components, and check
# it does not move. Seconds to run, and it fails exactly when the physics did.

DATA = REPO_ROOT / "docs/modelling-tour/_data"


def _generator():
    """The generator script, imported for its component-list definitions."""
    spec = importlib.util.spec_from_file_location(
        "generate_tour_equilibria",
        REPO_ROOT / "scripts/generate_tour_equilibria.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("asset, build", [
    ("rce_dry_equilibrium.npz", "dry_components"),
    ("rce_moist_equilibrium.npz", "moist_components"),
])
def test_shipped_equilibrium_is_still_an_equilibrium(states, asset, build):
    """The residual test. Ten steps must not move the shipped state.

    If this fails, the physics changed under a state that was computed before
    the change. Do not relax the thresholds: re-run
    `scripts/generate_tour_equilibria.py` (its defaults are what shipped) and
    commit the new file, then re-check every number the page quotes.
    """
    generator = _generator()
    components = getattr(generator, build)()
    tendencies, steppers = generator.split(components)

    state, provenance = states.load(
        str(DATA / asset), components,
        grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))
    assert state is not None, f"{asset} did not load"

    # The wind relaxation is rebuilt on the loaded state, exactly as the pages
    # do it and as the generator did. Omitting it would let the drag start
    # spinning the column down, and the "shipped state does not move" check
    # would be measuring that instead of the physics it is guarding.
    stepping_module = _load("stepping")
    tendencies = tendencies + [stepping_module.wind_relaxation(
        state, provenance["wind_m_s"], provenance["wind_timescale_hours"],
        initialise=False)]

    before = float(state["surface_temperature"].values.ravel()[0])
    timestep = climt.UnytTimeDelta(hours=provenance["dt_hours"])
    stepping_module.integrate(tendencies, steppers, state, timestep, 10)
    after = float(state["surface_temperature"].values.ravel()[0])

    budgets_module = _load("budgets")
    imbalance = budgets_module.toa_imbalance(state)

    assert abs(imbalance) < 1.0, (
        f"{asset}: TOA imbalance {imbalance:+.3f} W/m^2 after 10 steps — the "
        "shipped equilibrium is stale. Re-run "
        "scripts/generate_tour_equilibria.py.")
    assert abs(after - before) < 0.05, (
        f"{asset}: surface temperature drifted {after - before:+.4f} K in 10 "
        "steps — the shipped equilibrium is stale. Re-run "
        "scripts/generate_tour_equilibria.py.")


@pytest.mark.parametrize("asset", ["rce_dry_equilibrium.npz",
                                   "rce_moist_equilibrium.npz"])
def test_shipped_equilibrium_provenance_is_complete(states, asset):
    """Every field `states.describe` prints must actually be there.

    A page prints this block above its first figure. A `None` in it is a page
    telling a reader nothing while looking like it told them something.
    """
    generator = _generator()
    build = ("dry_components" if "dry" in asset else "moist_components")
    state, provenance = states.load(
        str(DATA / asset), getattr(generator, build)(),
        grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))

    for field in ("climt_version", "saved_at", "table", "nz", "dt_hours",
                  "n_steps", "slab_depth_m", "solar", "co2_ppm", "components",
                  "wind_m_s", "wind_timescale_hours", "roughness_length_m",
                  "toa_imbalance", "surface_imbalance"):
        assert provenance.get(field) is not None, f"{asset} lacks {field}"
    assert "None" not in states.describe(provenance)
```

- [x] **Step 2: Run it to verify it fails** — 4 failed (FileNotFoundError on the generator), as expected.

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k "shipped_equilibrium" -v`
Expected: FAIL — the generator script and both `.npz` files do not exist.

- [x] **Step 3: Write the generator** — written, with two deviations from the draft below, both recorded in the log.

Create `scripts/generate_tour_equilibria.py`:

```python
"""Generate the two radiative-convective equilibrium states pages 11 and 12 ship.

Pages 11 and 12 perturb an equilibrium rather than spending thousands of
in-browser steps finding one. This produces those states.

**The defaults here are exactly what shipped.** Running this with no arguments
must reproduce the committed files to within the convergence threshold. That is
not a nicety: `scripts/generate_tour_spectrum_table.py` shipped with defaults
that differed from what it had been run with, and reconstructing the real
invocation cost an afternoon.

    conda run -n climt python scripts/generate_tour_equilibria.py
    conda run -n climt python scripts/generate_tour_equilibria.py --moist
    conda run -n climt python scripts/generate_tour_equilibria.py --out /tmp

Regenerate deliberately, when `tests/test_modelling_tour.py`'s residual test
says the physics moved -- never on a schedule, and never from a dependency
hash. See that test for why.
"""
import argparse
import os
import sys

import numpy as np
import sympl

import climt

sys.path.insert(0, os.path.join("docs", "modelling-tour", "_tour"))
import budgets      # noqa: E402
import states       # noqa: E402
import stepping     # noqa: E402

DATA_DIR = os.path.join("docs", "modelling-tour", "_data")

# ---- the configuration the shipped states were made with. Changing any of
# ---- these invalidates both files and every number on pages 11 and 12.
NZ = 28
SOLAR = 240.0                 # prescribed absorbed shortwave at the surface
SLAB_DEPTH_M = 2.0
CO2_PPM = 330.0
TABLE = "earth_low_res_lw"
DRY_DT_HOURS = 1.0
MOIST_DT_HOURS = 0.25         # replace with Task 8's measured value
CONVERGENCE_W_M2 = 0.1        # |TOA imbalance| accepted as equilibrium
MAX_STEPS = 40000
SURFACE_SPECIFIC_HUMIDITY = 0.015    # saturated surface, page 9's default
WIND_M_S = 5.0                # the wind the relaxation holds the column at
WIND_TIMESCALE_HOURS = 24.0
ROUGHNESS_LENGTH_M = 1e-3


# The wind relaxation is NOT in these lists. It is built per-state by
# `stepping.wind_relaxation`, because it has to write its target fields into
# the state it will act on -- so it cannot exist before the state does. It is
# appended to the tendency list in `build_state`, and `split` never sees it.
# The residual test rebuilds it the same way; see `components_for`.


def dry_components():
    """Page 11's physics components. The single definition; the test imports it."""
    return [
        climt.CorkLongwaveRadiation(optics="correlated_k", table=TABLE),
        climt.SlabSurface(),
        climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                  roughness_length=ROUGHNESS_LENGTH_M),
        climt.DryConvectiveAdjustment(),
    ]


def moist_components():
    """Page 12's physics components. The single definition; the test imports it."""
    return [
        climt.CorkLongwaveRadiation(optics="correlated_k", table=TABLE),
        climt.SlabSurface(),
        climt.EmanuelConvectionPython(),
        climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                  roughness_length=ROUGHNESS_LENGTH_M),
        climt.DryConvectiveAdjustment(),
        climt.GridScaleCondensation(),
    ]


def split(components):
    """(tendency components, stepper components), in the order given.

    A sympl Stepper declares ``output_properties``; a TendencyComponent and an
    ImplicitTendencyComponent do not. ``EmanuelConvectionPython`` is the latter
    and lands in the tendency list, which is where climt has always run it --
    see _tour/stepping.py and page 12 for the warning that provokes.
    """
    tendencies = [c for c in components if not hasattr(c, "output_properties")]
    steppers = [c for c in components if hasattr(c, "output_properties")]
    return tendencies, steppers


def build_state(components, moist):
    """The state, and the wind relaxation that has to be built alongside it.

    Returns ``(state, relaxation)``. The relaxation is a tendency component,
    but it cannot be constructed before the state exists, because it writes
    its equilibrium and timescale fields into that state.
    """
    state = climt.get_default_state(
        components, grid_state=climt.get_grid(nx=1, ny=1, nz=NZ))
    state["ocean_mixed_layer_thickness"].values[:] = SLAB_DEPTH_M
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = CO2_PPM * 1e-6
    if moist:
        state["surface_specific_humidity"].values[:] = \
            SURFACE_SPECIFIC_HUMIDITY
    else:
        # Page 11's column is genuinely dry: no water vapour anywhere, so its
        # 14-band radiation sees CO2 alone. The page says so, and it is why
        # page 11 equilibrates well below Earth's surface temperature.
        state["specific_humidity"].values[:] = 0.0
        state["surface_specific_humidity"].values[:] = 0.0
    relaxation = stepping.wind_relaxation(
        state, WIND_M_S, timescale_hours=WIND_TIMESCALE_HOURS)
    return state, relaxation


def run_to_equilibrium(components, moist, dt_hours, check_every=50):
    tendencies, steppers = split(components)
    state, relaxation = build_state(components, moist)
    tendencies = tendencies + [relaxation]
    timestep = climt.UnytTimeDelta(hours=dt_hours)
    for step in range(0, MAX_STEPS, check_every):
        stepping.integrate(tendencies, steppers, state, timestep, check_every)
        imbalance = budgets.toa_imbalance(state)
        surface = float(state["surface_temperature"].values.ravel()[0])
        print(f"  step {step + check_every:6d}  "
              f"({(step + check_every) * dt_hours / 24.0:7.1f} d)  "
              f"TOA {imbalance:+8.4f} W/m^2  Tsurf {surface:7.3f} K")
        if abs(imbalance) < CONVERGENCE_W_M2:
            return state, step + check_every
    raise SystemExit(
        f"did not converge to {CONVERGENCE_W_M2} W/m^2 in {MAX_STEPS} steps "
        f"at dt = {dt_hours} h; last TOA imbalance "
        f"{budgets.toa_imbalance(state):+.4f}")


def generate(kind, out_dir):
    moist = (kind == "moist")
    components = moist_components() if moist else dry_components()
    dt_hours = MOIST_DT_HOURS if moist else DRY_DT_HOURS
    print(f"{kind} equilibrium: dt = {dt_hours} h, nz = {NZ}, "
          f"slab = {SLAB_DEPTH_M} m, CO2 = {CO2_PPM} ppm")

    state, n_steps = run_to_equilibrium(components, moist, dt_hours)
    provenance = dict(
        table=TABLE, nz=NZ, dt_hours=dt_hours, n_steps=n_steps,
        slab_depth_m=SLAB_DEPTH_M, solar=SOLAR, co2_ppm=CO2_PPM,
        wind_m_s=WIND_M_S, wind_timescale_hours=WIND_TIMESCALE_HOURS,
        roughness_length_m=ROUGHNESS_LENGTH_M,
        components=[type(c).__name__ for c in components] + ["UnytRelaxation"],
        toa_imbalance=budgets.toa_imbalance(state),
        surface_imbalance=budgets.surface_imbalance(state),
        convergence_threshold_w_m2=CONVERGENCE_W_M2,
    )
    path = os.path.join(out_dir, f"rce_{kind}_equilibrium.npz")
    states.save(path, state, provenance)
    print(states.describe({**provenance,
                           "climt_version": climt.__version__,
                           "saved_at": "(just now)"}))
    print(f"wrote {path} ({os.path.getsize(path) / 1024:.1f} kB)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry", action="store_true",
                        help="generate only the dry equilibrium")
    parser.add_argument("--moist", action="store_true",
                        help="generate only the moist equilibrium")
    parser.add_argument("--out", default=DATA_DIR,
                        help=f"output directory (default: {DATA_DIR})")
    args = parser.parse_args()

    sympl.set_backend(climt.UnytBackend())
    os.makedirs(args.out, exist_ok=True)
    wanted = [k for k, on in (("dry", args.dry), ("moist", args.moist)) if on]
    for kind in wanted or ["dry", "moist"]:
        generate(kind, args.out)


if __name__ == "__main__":
    main()
```

- [x] **Step 4: Set `MOIST_DT_HOURS` from Task 8's measurement** — set to `5.0 / 60.0` h (5 min), not the 6 h the `stability` probe reported. Task 8's CAVEAT and §Decisions both say 5 min: the 200-step stability test clears 6 h but the moist column goes non-finite over a real run there, and the repo's Emanuel examples all use 5 min. `DRY_DT_HOURS` was likewise set to 12 h (not the draft's 1.0) per §Decisions and page 11's quoted numbers — confirmed with the spec author. See the log.

`MOIST_DT_HOURS = 0.25` above is a placeholder marked as one. Replace it with the largest stable timestep Task 8's `stability` measured for `12-moist-rce`, and delete the `# replace with` comment. If Task 8 has not run, stop — this task depends on it.

- [ ] **Step 5: Generate both states**

Run with the JIT **on**: this is a native generation, not a browser cost model, and the 14-band radiation is ~20× faster compiled.

```bash
conda run -n climt python scripts/generate_tour_equilibria.py 2>&1 | tail -40
ls -la docs/modelling-tour/_data/*.npz
```

Expected: both files written, each a few kB, each preceded by a converging trace and a provenance block. Record the final TOA imbalance, surface temperature and step count for both in the log below.

If the dry run does not converge, it is almost certainly not stuck — a 2 m slab with a ~2 W m⁻² K⁻¹ feedback has a response time of about three weeks, so `MAX_STEPS = 40000` at `dt = 1 h` is 4.5 years of simulated time and ample. Read the trace: a TOA imbalance that is decreasing but slowly needs more steps; one that is oscillating needs a smaller `dt`; one that is flat and non-zero means a component is not seeing what you think it is.

- [ ] **Step 6: Run the residual test against the real files**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k "shipped_equilibrium" -v`
Expected: all four pass (two assets × two tests), in seconds.

- [ ] **Step 7: Prove the generator's defaults reproduce what shipped**

The generator's contract is that running it with no arguments reproduces the committed files. Check it, rather than trusting it:

```bash
conda run -n climt python scripts/generate_tour_equilibria.py --out /tmp/eq_check
conda run -n climt python -c "
import sys, numpy as np
sys.path.insert(0, 'docs/modelling-tour/_tour')
for kind in ('dry', 'moist'):
    a = np.load(f'docs/modelling-tour/_data/rce_{kind}_equilibrium.npz')
    b = np.load(f'/tmp/eq_check/rce_{kind}_equilibrium.npz')
    T_a, T_b = a['value__air_temperature'], b['value__air_temperature']
    print(kind, 'max |dT| =', float(np.abs(T_a - T_b).max()), 'K')"
```

Expected: `max |dT|` below ~0.05 K for both — the convergence threshold's worth of slack, not more. A larger difference means the run is not deterministic at the level the pages quote, and the pages' precision has to come down to match.

- [x] **Step 8: Document both files in `_data/README.md`** — appended, following the `earth_spectrum_lw.npz` shape; added a timestep row and a "what converged means here" note explaining the surface-temperature-stationarity gate.

Append a section, following the shape of the existing `earth_spectrum_lw.npz` one:

```markdown
## `rce_dry_equilibrium.npz` and `rce_moist_equilibrium.npz`

Two single-column radiative-convective equilibrium states, loaded by pages 11
and 12 so those pages can run perturbation experiments instead of spending
thousands of in-browser steps spinning up. `_tour/states.py` loads them;
`_tour/assets.py` finds them.

| | dry (page 11) | moist (page 12) |
|---|---|---|
| components | Cork LW, SlabSurface, SimpleBoundaryLayer, DryConvectiveAdjustment | + EmanuelConvectionPython, GridScaleCondensation |
| water vapour | none — a CO₂-only atmosphere | evolves, supplied by the surface |
| LW table | `earth_low_res_lw` (14 bands) | `earth_low_res_lw` |
| grid | nz = 28, one column | nz = 28, one column |
| slab depth | 2 m | 2 m |
| absorbed SW | 240 W m⁻², prescribed at the surface | same |
| CO₂ | 330 ppm | 330 ppm |
| size | a few kB | a few kB |

Each file carries its own provenance — climt version, timestep, step count,
and the TOA and surface imbalances it finished at. `_tour/states.describe()`
prints it, and pages 11 and 12 print it above their first figure, so a reader
is never looking at an equilibrium without also seeing what produced it.

**Page 11's column is genuinely dry.** Its 14-band radiation sees CO₂ alone, so
it equilibrates well below Earth's surface temperature and no number on that
page is comparable with tranche 1's. The difference between page 11 and page 12
is, to first order, the water vapour greenhouse plus its feedback.

### Regenerating

Only when `tests/test_modelling_tour.py::test_shipped_equilibrium_is_still_an_equilibrium`
fails — that test is the staleness guard, deliberately in place of the
dependency-hash machinery in `scripts/build_experiments.py`, which cannot tell
a no-op in `cork/` from a physics change and has charged for the difference.

```sh
conda run -n climt python scripts/generate_tour_equilibria.py
conda run -n climt python -m pytest tests/test_modelling_tour.py -k shipped_equilibrium
```

The script's defaults are exactly what shipped. After regenerating, re-check
every number pages 11 and 12 quote — the state moved, so they may have too.
```

- [ ] **Step 9: Confirm the states are published and stageable**

`docs/_quarto.yml` already lists `modelling-tour/_data/*.npz` under `project: resources:`, so both new files are published without further edits. Verify:

```bash
cd docs && quarto render modelling-tour/index.qmd && cd ..
ls docs/_site/modelling-tour/_data/
```

Expected: both `.npz` files present alongside `earth_spectrum_lw.npz`. Each page that loads one declares it under its own `pyodide: resources:` — that is Tasks 14 and 15.

- [ ] **Step 10: Commit**

```bash
git add scripts/generate_tour_equilibria.py \
        docs/modelling-tour/_data/rce_dry_equilibrium.npz \
        docs/modelling-tour/_data/rce_moist_equilibrium.npz \
        docs/modelling-tour/_data/README.md tests/test_modelling_tour.py
git commit -m "feat(tour): two shipped RCE equilibria, and the residual test that guards them

Pages 11 and 12 perturb an equilibrium rather than spinning one up in the
browser. Staleness is caught by stepping the shipped state ten steps and
checking it does not move -- not by dependency hashing, which cannot tell a
no-op in cork/ from a physics change.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

### Log — the two equilibria as generated

**Run 2026-09-04, climt 0.31.0, numba ON, reference machine (Apple Silicon).**

| | dry | moist |
|---|---|---|
| `dt` | 12 h | 5 min |
| steps to converge | 3500 | 200 550 |
| simulated days | 1750 d | ≈ 696 d |
| final TOA imbalance (W m⁻²) | −0.130 | −0.485 |
| final surface imbalance (W m⁻²) | +1.384 (instantaneous — see note) | −1.565 (instantaneous — see note) |
| surface temperature (K) | 266.478 | 279.965 |
| file size | 28.0 kB | 31.2 kB |
| native wall clock | ~35 s (after numba compile) | ~25 min (inferred from the two files' `saved_at` stamps, 04:34:58 → 04:59:51 UTC; not timed directly) |

The moist figures are read back from `rce_moist_equilibrium.npz`'s own
provenance (`__meta__`), which the run wrote on 2026-09-04. Two things to
note. 200 550 steps is well over Task 0's cold-start ~123 000; that is the
same surface-temperature-stationarity gate that took the dry run from 1650 to
3500 (deviation 2 below). And the final TOA imbalance, −0.485 W m⁻², only just
clears the 0.5 gate, so pages 12 onward should not quote it to better than
±0.5 W m⁻².

**Two deviations from the draft generator in Step 3, both forced by Task 8 and
by the residual test written first:**

1. **The four pre-Task-8 constants were reconciled to the measurements.**
   `DRY_DT_HOURS` 1.0 → **12** (§Decisions / page 11, confirmed with the spec
   author); `MOIST_DT_HOURS` 0.25 → **5/60 h** (Task 8 CAVEAT: the 200-step
   stability test wrongly clears 6 h); `CONVERGENCE_W_M2` 0.1 → **0.5** (the
   gate the campaign measured with, per §4 "ships directly from this
   convergence"); `MAX_STEPS` 40000 → **400000** (the moist run needs
   ~123 000). `build_state` sets up the same initial conditions the Task 8
   measurement did: CO₂ 330 ppm, `specific_humidity` 0 to start, and
   `surface_specific_humidity` **0.015** for the moist column (a saturated
   surface, page 9's default; both `build_state` and
   `tour_rce_measurements.py` set it) or 0 for the dry one. So the two runs
   differ only in where they were stopped.

2. **Convergence is TOA balance AND surface-temperature stationarity, not the
   draft's TOA-only gate.** The draft (TOA-only, |TOA|<0.5) stopped the dry run
   at 1450 steps on a *transient dip* of the oscillating TOA, while the surface
   was still +3.5 W m⁻² out of balance and descending — the residual test
   (written first) caught it: the surface drifted 0.0513 K > 0.05 in ten steps.
   Adding a surface **flux** gate fails too: the one-step flux lag
   `_tour/stepping.py` documents makes `surface_imbalance` oscillate ±several
   W m⁻² step to step forever (measured: a clean +3.6 / −8.0 / +1.6 / −2.0
   limit cycle), so it never sits under any small threshold. The surface
   **temperature**, though, settles into a ±0.05 K limit cycle. So the second
   gate is that the mean surface temperature over one 1000-step window equals
   the previous window's mean to within 0.02 K (means cancel the oscillation
   whatever its period). This is why the shipped dry state is 3500 steps, not
   Task 0's cold-start 1650 (measured numba-OFF, TOA-only): those numbers stop
   mid-transient. Page 11 should quote the shipped 3500 / 266.48 K, not 1650 /
   266.60 K.


**A third state, added 2026-09-26: `rce_moist_2xco2_equilibrium.npz`.** Page 12's 2×CO₂ experiment is 73 600 steps from the moist state, about 3.9 h in the browser, so it ships instead of running (Task 8 log §7, Task 15 cell 4). `generate_tour_equilibria.py --moist-2xco2` makes it. It loads the shipped moist state, rebuilds the wind relaxation with `initialise=False`, sets CO₂ to 660 ppm and re-converges under the same two gates. The generator's `run_to_equilibrium` now takes the state to step rather than building one, so the same loop serves the spin-ups, this perturbation and the start-independence re-measurement. The file records `perturbed_from`, `perturbed_from_saved_at`, `base_surface_temperature_k` and `warming_k`. `test_shipped_2xco2_state_is_the_shipped_moist_state_perturbed` pins the stamp to the shipped base, and the residual and provenance tests cover the new file.

| | moist 2×CO₂ |
|---|---|
| steps after doubling | 73 600 (≈ 256 d) |
| surface temperature | 281.613 K (+1.647 K on the base file; +1.82 K once both settle, see Task 8 log §7) |
| final TOA imbalance | +0.500 W m⁻² (just inside the gate) |
| file size | 31.3 kB |
| native wall clock | ~9 min, Linux, `NUMBA_NUM_THREADS=1` (≈ 7 ms/step) |

**Both moist states regenerated, 2026-09-27: a saturated surface, and a gate the moist column can pass.** The tables above are history. They describe files that no longer ship.

*Why.* `build_state` wrote `surface_specific_humidity = 0.015` once, and nothing updated it. At the 279.97 K the old base settled at, q_sat(Ts, ps) is 6.09 g/kg, so the surface sat at **246 % relative humidity** (220 % in the 2×CO₂ state). It evaporated from something wetter than water: LH 140.6 W m⁻² against SH 2.6, a Bowen ratio of 0.02. Page 09 teaches that a fixed surface q is wrong for exactly this reason, and introduces `stepping.SurfaceHumidity(rh)` to fix it. Page 12 cannot load a state built the way page 09 calls wrong.

*What changed in `scripts/generate_tour_equilibria.py`.*
1. `moist_components()` now includes `stepping.SurfaceHumidity(1.0)`. `split` keeps list order, so it is the first stepper, ahead of the boundary layer, which is page 09's order. The fixed value and its constant are gone. The residual test, `split` and the provenance all read this one list. Provenance records `surface_relative_humidity = 1.0`, and `test_shipped_moist_states_have_a_saturated_surface` checks both files for it. The dry path is unchanged (surface q = 0).
2. **The moist column has its own gate, `MOIST_GATE`.** The first regeneration used the old gate and was rejected. That gate was an instantaneous |TOA| < 0.5 plus a 1000-step drift window. It stopped at 194 400 steps (286.221 K, file TOA −0.18), but the 30-day mean TOA was still −1.7 W m⁻² and falling, and ten steps on the file read −1.32, failing the residual test. The saturated column is noisy: episodic convection moves the surface flux ±40 W m⁻² and the instantaneous TOA by ~0.5 W m⁻² (1 sd). Its late drift is slow, ~0.05 K per 30 days, and at 5 min a 1000-step window is only 3.5 days. Stepping that state on a further 150 000 steps showed where the column settles: **Ts 286.00 K, TOA −1.1 W m⁻², steady.** The stratosphere is in radiative balance (−0.02 W m⁻² of LW heating above 150 hPa), and the slab is not storing anything measurable. So −1.1 W m⁻² is the Emanuel scheme's non-conservation, as the +0.3 was for the old column (Task 8 §7), now with the other sign. *(Superseded 2026-09-27, Task 15 log: the −1.1 W m⁻² is not Emanuel's. Metering column enthalpy per component, `DryConvectiveAdjustment` adds +0.97 W m⁻² and Emanuel +0.03. The adjustment conserves a moist-c_p enthalpy that the rest of the stack, counting dry c_p, does not. The gate reasoning is unaffected.)* No |TOA| threshold below that can be a convergence test. The new gate therefore works on *trends* over 60-day (17 280-step) windows. The mean surface temperature and the mean TOA must each match the previous window's, to 0.01 K and 0.1 W m⁻², and the window-mean TOA must be inside ±2 W m⁻². That last condition rules out the turning point near step 28 000, where the surface peaks while TOA is still −25. Replayed on the measured trajectory, the 30-day windows tried first stopped 0.1 K short. The 60-day gate stops 0.01 K from the settled state. The dry gate is unchanged. Provenance records `toa_gate`, `toa_steady_w_m2`, `window_mean_toa_w_m2`, `window_mean_surface_temperature_k`, and, for 2×CO₂, `window_mean_warming_k`.
3. **The residual test compares against the file's own settled TOA.** It now requires |TOA after 10 steps − `window_mean_toa_w_m2`| < 1.0. Dry files lack the key and are compared against 0, as before. The 0.05 K surface-drift check is unchanged.

*As generated* (Linux, 4 cores, numba on, `NUMBA_NUM_THREADS=1`, ≈ 7 ms/step, generator defaults: `--moist --moist-2xco2`):

| | moist | moist 2×CO₂ |
|---|---|---|
| steps | 295 800 (≈ 1027 d) | 180 550 after doubling (≈ 627 d) |
| native wall clock | 34 min | 21 min |
| surface temperature (file / last-window mean) | 285.989 / 286.014 K | 288.251 / 288.237 K |
| TOA imbalance (file / last-window mean) | −1.270 / −1.095 W m⁻² | −1.197 / −1.029 W m⁻² |
| surface RH | 100 % | 100 % |
| SH / LH (30-day mean after loading) | 27.9 / 76.3 W m⁻² | 25.3 / 83.2 W m⁻² |
| Bowen ratio | 0.37 | 0.30 |
| precipitation / evaporation | 2.64 / 2.64 mm day⁻¹ | 2.88 / 2.88 mm day⁻¹ |
| lowest-level q | 3.2–3.3 g kg⁻¹ | 3.8–4.1 g kg⁻¹ |
| warming | | +2.262 K file to file; +2.223 K window mean to window mean |
| residual test (10 steps) | TOA −0.645, ΔTs −0.0013 K | TOA −1.101, ΔTs −0.0040 K |
| file size | 31.3 kB | 31.4 kB |

For comparison, the old files measured over 1000 steps under their own fixed-q physics: base 279.965 K, 246 % RH, SH 2.6 / LH 140.6 W m⁻², Bowen 0.018, P 4.82 mm day⁻¹, TOA −0.09, lowest-level q 5.7 g kg⁻¹. The 2×CO₂ state: 281.613 K, 220 %, SH 6.7 / LH 133.1, Bowen 0.050, P 4.60, TOA +0.50, q 6.0 g kg⁻¹.

*Re-measured* (`tour_rce_shipped_remeasure.py`, 2026-09-27): `settle-moist` steps each shipped state on 60 000 steps and averages the last 30 000. The base settles at **285.990 K**, TOA −1.045 (std 0.48), trend −0.026 K per 100 d, which is 0.001 K from its file. The 2×CO₂ state settles at **288.228 K**, TOA −1.066 (std 0.54), trend −0.001 K per 100 d, 0.02 K below its file. **Settled 2×CO₂ warming +2.24 K.** `indep-moist` starts warm and steep (328 K air at 9 K/km, 330 K surface) under the new gate. It stopped at 307 700 steps with Tsurf 286.037 K, against the shipped 285.989 K. The column's maximum |ΔT| is 0.66 K at 3 hPa, and 0.50 K below 200 hPa. Both are single-instant comparisons in a column whose surface wanders ~0.05 K step to step and whose convection is episodic, so the 0.05 K surface difference is at the noise level. **Start-independent to ≈ 0.05 K at the surface.** Task 8 §7's 0.002 K belonged to the smoother fixed-q column.

The saturated column is 6 K warmer than the fixed-q one. The old surface sent about 140 W m⁻² into latent heat, far more than the column could radiate away, and it ended cold and wet near the ground. Its 2×CO₂ response grows from +1.82 K to +2.24 K. That is about 1.9× the dry +1.20 K rather than 1.5×, *(superseded 2026-09-27: page 11 measures and quotes the dry sensitivity as +1.15 K, 30-day mean minus 30-day mean, Task 14 log; against that the moist +2.24 K is also about 1.9×)* because the surface now dries and moistens with its own temperature, so a warmer surface evaporates more and the water-vapour feedback is no longer muted. The browser cost of the 2×CO₂ experiment is now 180 550 × ~0.19 s ≈ 9.5 h, so it still ships.

---

## Task 10: Page 07 — Stepping a column

*Course chapters 9–10, and the bridge from tranche 1.*

The first page in the tour where anything moves. A single column over a slab with prescribed absorbed shortwave at the surface, stepped with `AdamsBashforth`, drawing the evolution.

**The reveal, and it was verified before this plan was written.** Started from climt's arbitrary default profile, the gray column *finds* the analytic equilibrium page 04 prescribed and verified. Page 04 checked that ℋ ≈ 0 on a profile handed to it; page 07 produces that profile from nothing. Measured, at `table="tour_gray_lw"`, `diffusivity_factor=2.0`, `nz=28`, 2 m slab, `SOLAR=240`, `dt=12 h`, 900 steps (450 simulated days):

| | model | analytic (chapter 8) | difference |
|---|---|---|---|
| surface temperature | 335.53 K | 335.60 K | −0.07 K |
| max \|T − T(τ)\| through the column | — | — | **0.103 K** |
| rms | — | — | 0.077 K |
| top level | 215.06 K | `Te/2^(1/4)` = 214.43 K | **+0.63 K** |
| OLR | 239.97 W m⁻² | 240 | TOA imbalance +0.03 W m⁻² |

Those are the numbers the page quotes, and the tolerances the tests use. **Note the table:** `tour_gray_lw` with `D = 2`, not `single_band_gray_lw`. The reveal only lands against the table and diffusivity page 04 calibrated — `single_band_gray_lw` reproduces climt's *default* gray scheme, which is a different τ. Page 07 then switches to `earth_low_res_lw` for the second half, where the whole point is that a different table finds a different equilibrium.

**Files:**
- Create: `docs/modelling-tour/07-stepping-a-column.qmd`
- Modify: `docs/_quarto.yml` (sidebar), `docs/modelling-tour/index.qmd` (chapter map), `docs/modelling-tour/_artifacts/generate.py` (fallback figure)
- Test: `tests/test_modelling_tour.py`

**Interfaces:**
- Consumes: `_tour/stepping.py`, `_tour/budgets.py`, `_tour/soundings.py` (`analytic_gray_equilibrium`).
- Produces: nothing other pages import. Pages 08–12 assume its craft content: `TendencyComponent` vs `Stepper`, `AdamsBashforth`, `UnytTimeDelta`, update order, `state["time"]`.

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_modelling_tour.py`:

```python
# `PAGE7_GRAY_TABLE`, `PAGE7_DIFFUSIVITY` and `_page7_gray_column` were added
# in Task 5, when the inherited comparison tests moved onto page 7's real
# tables. Add only the two analytic constants here.
PAGE7_TAU_INF = 4.0
PAGE7_TE = 255.0


@pytest.fixture(scope="module")
def page7_gray_equilibrium():
    """The gray column, stepped to equilibrium once and shared."""
    import sympl as _sympl

    _sympl.set_backend(climt.UnytBackend())
    components, state = _page7_gray_column()
    return _load("stepping").integrate(
        components, [], state, climt.UnytTimeDelta(hours=12), 900)


@pytest.mark.slow
def test_page7_gray_column_finds_page4s_analytic_profile(
        page7_gray_equilibrium, soundings):
    """Page 7's reveal, and the strongest claim in the tranche.

    Page 4 checked that the heating rate vanishes on a profile it was handed.
    Page 7 hands the model nothing and gets the same profile back. Measured
    max |dT| is 0.103 K; the threshold is 0.5 K, comfortable but not vacuous.
    """
    state = page7_gray_equilibrium
    p = state["air_pressure"].values[:, 0, 0]
    surface_pressure = float(state["surface_air_pressure"].values.ravel()[0])

    T_analytic, T_ground_analytic, _ = soundings.analytic_gray_equilibrium(
        p, surface_pressure, tau_inf=PAGE7_TAU_INF, T_e=PAGE7_TE)
    T = state["air_temperature"].values[:, 0, 0]

    assert np.max(np.abs(T - T_analytic)) < 0.5, (
        f"max |T - T(tau)| = {np.max(np.abs(T - T_analytic)):.3f} K; the "
        "column did not find chapter 8's profile")
    surface = float(state["surface_temperature"].values.ravel()[0])
    assert abs(surface - T_ground_analytic) < 0.5, (
        f"surface {surface:.2f} K vs analytic {T_ground_analytic:.2f} K")


@pytest.mark.slow
def test_page7_top_level_sits_near_the_skin_temperature(
        page7_gray_equilibrium):
    """The isothermal top the figure labels is the skin temperature.

    Measured +0.63 K above Te/2^(1/4); the threshold is 1.0 K. The residual is
    real and physical -- the top level is a finite layer, not the tau -> 0
    limit -- so this is not a tolerance to tighten toward zero.
    """
    T_top = float(page7_gray_equilibrium["air_temperature"].values[-1, 0, 0])
    skin = PAGE7_TE / 2.0 ** 0.25
    assert abs(T_top - skin) < 1.0, (
        f"top level {T_top:.2f} K vs skin temperature {skin:.2f} K")


@pytest.mark.slow
def test_page7_gray_column_reaches_energy_balance(page7_gray_equilibrium,
                                                  budgets):
    """The convergence claim, read off the budget rather than a curve."""
    assert abs(budgets.toa_imbalance(page7_gray_equilibrium)) < 0.2


@pytest.mark.slow
def test_page7_mixed_layer_depth_changes_speed_not_equilibrium(soundings):
    """Page 7's knob, and the cleanest lesson in the tranche.

    1 m and 5 m slabs reach the *same* equilibrium at very different speeds.
    Heat capacity sets response time; it does not set where you end up.
    """
    stepping_module = _load("stepping")
    finals = {}
    for depth in (1.0, 5.0):
        components, state = _page7_gray_column(slab_depth=depth)
        stepping_module.integrate(components, [], state,
                                  climt.UnytTimeDelta(hours=12), 1400)
        finals[depth] = float(state["surface_temperature"].values.ravel()[0])

    assert abs(finals[1.0] - finals[5.0]) < 0.5, (
        f"1 m settled at {finals[1.0]:.2f} K and 5 m at {finals[5.0]:.2f} K — "
        "slab depth must not change the equilibrium")
```

- [ ] **Step 2: Run them to verify they fail, then pass**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k page7 -m slow -v`
Expected: the three new equilibrium tests should **pass immediately** — they test physics, not the page, and the physics is already there. That is fine and expected: they are the page's contract, written before the page so the page cannot quote a number they do not cover.

`test_page7_mixed_layer_depth_changes_speed_not_equilibrium` is the one that may not: 1400 steps at `dt = 12 h` is 700 days, and a 5 m slab is 2.5× slower than the 2 m one measured above. If it fails with the two temperatures still converging toward each other, raise the step count; if they settle at genuinely different values, stop — the knob's whole lesson is wrong and the spec needs revisiting.

- [ ] **Step 3: Write the page**

Create `docs/modelling-tour/07-stepping-a-column.qmd`. Front matter — copy the shape from `04-gray-equilibrium-tested.qmd` and change the `resources:` list:

```yaml
---
title: "Stepping a column"
format: live-html
engine: jupyter
pyodide:
  packages:
    # Installed by micropip at document setup, before any cell runs; pinned to
    # the released version this page was written against. See
    # ../_includes/climt-live-boot.qmd for how to preview an unreleased wheel.
    - "climt==0.31.0"
    - "unyt"
  resources:
    # Fetched into the Pyodide filesystem at document setup — the same files
    # tests/test_modelling_tour.py imports natively. See page 1 for the
    # mechanism and docs/_quarto.yml `project: resources:` for why they ship.
    - _tour/soundings.py
    - _tour/stepping.py
    - _tour/budgets.py
---

{{< include ../_includes/climt-live-boot.qmd >}}
```

The page's anatomy is tranche 1's, unchanged: chapter anchor callout, visible prescribed state, the calls, a figure, one clearly marked knob, Physics and Code exercises, a "going deeper" link.

**Cell 0 — setup and the column.** This is the cell `_artifacts/generate.py` will exec, so it must stand alone.

```python
import sys
sys.path.insert(0, "_tour")   # helper modules, served alongside this page

import numpy as np
import sympl
import climt
from climt import (get_default_state, get_grid, CorkLongwaveRadiation,
                   SlabSurface, UnytTimeDelta)

import budgets
import soundings
import stepping

sympl.set_backend(climt.UnytBackend())

SOLAR = 240.0        # absorbed shortwave at the surface, W/m^2 — prescribed,
                     # because this tour runs no shortwave component.
NZ = 28
SLAB_DEPTH = 2.0     # metres of water. The knob.


def gray_column(table="tour_gray_lw", diffusivity=2.0, slab_depth=SLAB_DEPTH):
    """One air column over a slab ocean, warmed by a fixed absorbed sunlight."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table=table,
                               diffusivity_factor=diffusivity)
    surface = SlabSurface()
    state = get_default_state([lw, surface], grid_state=get_grid(nx=1, ny=1, nz=NZ))
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    return [lw, surface], state


components, state = gray_column()
print("starting from climt's default profile:")
print("  surface", float(state["surface_temperature"].values.ravel()[0]), "K")
print("  air     ", state["air_temperature"].values[:, 0, 0].round(1)[:5], "... K")
```

**Cell 1 — the gray run and the figure.** This is the headline cell; `_artifacts/generate.py` captures its figure.

```python
DT_HOURS, N_STEPS = 12, 900     # 450 simulated days

components, state = gray_column()
state, history = stepping.integrate_with_history(
    components, [], state, UnytTimeDelta(hours=DT_HOURS), N_STEPS)
stepping.draw_evolution(history, state, SOLAR,
                        title="Gray radiative equilibrium, found rather than assumed")
print(budgets.summary(state))
```

**Cell 2 — the reveal, against page 04.**

```python
p = state["air_pressure"].values[:, 0, 0]
ps = float(state["surface_air_pressure"].values.ravel()[0])
T_analytic, T_ground, _ = soundings.analytic_gray_equilibrium(
    p, ps, tau_inf=4.0, T_e=255.0)
T = state["air_temperature"].values[:, 0, 0]

print(f"max |T_model - T(tau)| = {np.abs(T - T_analytic).max():.3f} K")
print(f"surface: model {float(state['surface_temperature'].values.ravel()[0]):.2f} K"
      f"  analytic {T_ground:.2f} K")
print(f"top level: {T[-1]:.2f} K   skin temperature Te/2^(1/4): {255.0 / 2**0.25:.2f} K")
```

**Cell 3 — change one string.**

```python
components, state = gray_column(table="earth_low_res_lw", diffusivity=1.66)
state, history = stepping.integrate_with_history(
    components, [], state, UnytTimeDelta(hours=DT_HOURS), N_STEPS)
stepping.draw_evolution(history, state, SOLAR,
                        title="Non-grey radiative equilibrium (14 bands)")
print(budgets.summary(state))
```

**Cell 4 — the knob.** Mixed-layer depth: 1 m and 5 m, same equilibrium, different speed. Plot the two surface-temperature time series on one axis.

**Prose the page owes the reader, each in its own callout:**

1. *Where this fits* — chapters 9–10, and the callback to page 04: "Page 4 checked that ℋ ≈ 0 on a profile we wrote down. This page writes nothing down."
2. *climt craft: three kinds of component* — `TendencyComponent` returns rates and needs an integrator; `Stepper` returns a new state; the third protocol, `ImplicitTendencyComponent`, exists and page 12 meets it. `AdamsBashforth` wraps the tendency components.
3. *climt craft: the timestep is not a `timedelta`* — show the `UnitOperationError` and explain it. `total_seconds()` on a plain `timedelta` returns a bare float, which will not cancel a tendency's `/s`; `UnytTimeDelta.total_seconds()` returns seconds *with units*, which does.
4. *climt craft: update order* — the prognostic state is applied before the diagnostics, or the stepper carries flux diagnostics forward at their pre-step value and `SlabSurface` consumes stale longwave. This broke the previous version of this demo once already.
5. *This column is dry, and it does not matter* — `CorkLongwaveRadiation` with a gray table does not take `specific_humidity` as an input at all; its optical depth is a function of pressure alone. So the gray column reproduces page 04's profile whatever the reader does to humidity. Say it here rather than leaving them to wonder. **And climt's default `specific_humidity` is 0.0 everywhere**, so the 14-band column in cell 3 is a CO₂-only atmosphere — which is exactly what makes page 12's comparison worth drawing. (The page this one replaces described "spectral windows letting surface radiation escape" as though water vapour were present. It was not.)
6. *Knob* — mixed-layer depth, clearly marked.
7. *Physics exercises* — (i) predict, before running, whether a 10 m slab changes the equilibrium; (ii) from the OLR time series, estimate the e-folding time and check it against `ρ c_p h / (4σT³)`.
8. *Code exercises* — (i) raise `DT_HOURS` until it breaks, and read the error (Task 2's guard names the timestep); (ii) print `state["time"]` before and after and confirm the clock advanced by `N_STEPS × DT_HOURS`.
9. *Going deeper* — link to `../radiative-transfer/07-two-stream.qmd`.

Every number in the prose comes from a cell above it. Do not write a number the reader cannot regenerate by pressing run.

- [ ] **Step 4: Add the page to the sidebar and the chapter map**

In `docs/_quarto.yml`, after `modelling-tour/06-water-vapour-limit.qmd`:

```yaml
        - modelling-tour/07-stepping-a-column.qmd
```

In `docs/modelling-tour/index.qmd`, add a row to the chapter map table:

```markdown
| [Stepping a column](07-stepping-a-column.qmd) | 9, 10 | radiative equilibrium, found rather than assumed |
```

Leave the "What this tranche does not cover" section alone for now — Task 16 rewrites it once all six pages exist, so a half-finished branch does not contradict itself in two places.

- [ ] **Step 5: Add the static fallback figure**

In `docs/modelling-tour/_artifacts/generate.py`, add to `FIGURES`:

```python
    "07-stepping-a-column.qmd": [(1, "07-stepping.png")],
```

Run: `conda run -n climt python docs/modelling-tour/_artifacts/generate.py 07`
Expected: `07-stepping.png` written. This execs the page's own cells, so it is also the check that the page runs at all. Note the wall clock: 900 gray steps at 2.5 ms JIT-on is ~2 s, so if this takes minutes something is wrong.

Then add the collapsed fallback callout to the page, following the pattern in `04-gray-equilibrium-tested.qmd`.

- [ ] **Step 6: Render and check in a browser**

```bash
cd docs && quarto render modelling-tour/07-stepping-a-column.qmd && quarto preview
```

Expected: the page loads; the first cell prints the climt version; the gray cell runs in ~20 s (900 steps × ~8 ms browser) and draws four populated panels; cell 3 takes ~2.5 minutes (900 × ~0.17 s) and the page says so before the cell, so the spinner is expected rather than alarming. Check the log-pressure axis and that the gray figure's top annotation reads "isothermal top (skin temperature)" while the 14-band one reads "no isothermal top".

**If the 14-band cell is over ~3 minutes**, reduce `N_STEPS` for that cell only and say in the prose that it stops short of equilibrium and by how much — do not switch tables or drop levels to make it finish.

- [ ] **Step 7: Commit**

```bash
git add docs/modelling-tour/07-stepping-a-column.qmd docs/modelling-tour/index.qmd \
        docs/modelling-tour/_artifacts/generate.py docs/modelling-tour/_artifacts/07-stepping.png \
        docs/_quarto.yml tests/test_modelling_tour.py
git commit -m "docs(tour): page 7, where the column finds page 4's profile by itself

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 11: Page 08 — Turbulent heat exchange

*Course chapter 10.*

Adds `SimpleBoundaryLayer(surface_fluxes='bulk')` to page 07's gray column. The scheme computes Frierson (2006) bulk fluxes from its own exchange coefficient, applies them implicitly, and reports the applied fluxes as diagnostics — so the flux the reader plots is the flux the model used.

**The reveal.** Page 04 derived radiative equilibrium's surface–air temperature discontinuity as an analytic feature and left it there. This is the process that erases it — 15.15 K in radiative equilibrium, down to 6.07 K with a boundary layer over a wind-forced column.

**And the page has a second subject**, which the measurements below forced: a column with no dynamics stops moving, and you have to give it a wind. That turns out to be the better half of the page.

**Measured before this plan was written**, on page 07's gray column (`tour_gray_lw`, `D = 2`, `nz = 28`, 2 m slab, `SOLAR = 240`), each run to a TOA imbalance under 0.02 W m⁻² — 12 000 steps at `dt = 1 h`, 500 simulated days:

| configuration | surface − lowest air (K) | sensible heat flux (W m⁻²) | surface T (K) |
|---|---|---|---|
| radiative equilibrium (no boundary layer) | **15.15** | 0 | 335.53 |
| + boundary layer, `z0 = 3.21e-5` (default) | **11.37** | 9.7 | 334.88 |
| + boundary layer, `z0 = 1e-3` | 10.04 | 14.9 | 334.50 |
| + boundary layer, `z0 = 0.1` | **6.93** | 30.7 | 333.25 |

**The wind problem, and the page's best material.** Setting `eastward_wind` to 0, 5 or 10 m s⁻¹ and integrating changes the converged jump by 0.05 K, and the lowest-level wind ends at 0.000 m s⁻¹ in all three. Nothing drives a wind in a column with no dynamics, and the scheme's own surface drag removes the one you prescribe within days. **So the page must supply a momentum source**, and doing so is the second half of what it teaches.

climt's own examples do this by brute force — `examples/column_code_with_slab.py:120` and `examples/gmd_radiative_convective.py:118` both end every timestep with `state['eastward_wind'].values[:] = 3.` The page shows that idiom, names it as what climt has always done, and then shows the cleaner form: a **Newtonian relaxation**, `dx/dt = −(x − x_eq)/τ`, as a `sympl` component. It composes, it carries units, it appears in the component list where a reader can see the assumption, and its timescale is a physical parameter rather than a hidden assignment. `_tour/stepping.wind_relaxation` (Task 4, Step 8) builds it.

**With the wind held up, wind speed becomes a real knob.** Measured, `z0 = 1e-3`, 12 000 steps at `dt = 1 h`, every run converged to |TOA| < 0.02 W m⁻²:

| configuration | lowest-level wind (m s⁻¹) | surface − air (K) | sensible heat flux (W m⁻²) | surface T (K) |
|---|---|---|---|---|
| no momentum source | **0.00** | 10.04 | 14.9 | 334.50 |
| relax to 2 m s⁻¹ | 1.34 | 9.31 | 19.9 | 334.13 |
| relax to 5 m s⁻¹ | 2.88 | 7.50 | 34.9 | 332.97 |
| relax to 10 m s⁻¹ | 5.07 | **6.07** | **49.9** | 331.56 |
| hard reset to 5 m s⁻¹ | 5.00 | 5.84 | 42.9 | 332.12 |

**Read the last two rows together — that comparison is a lesson on its own.** Relaxing toward 5 m s⁻¹ leaves the lowest level at 2.88, because the surface drag is still acting on it and the relaxation only pulls; the hard reset pins it at exactly 5.00 and so overrides the drag entirely at the one level where the exchange happens. The reset therefore delivers more flux (42.9 vs 34.9 W m⁻²) from the same nominal wind. The relaxation is the more honest scheme: it lets the drag win near the ground, which is what drag does. Say this on the page — it is the difference between a forcing and an assignment, and a student who sees it will never confuse the two again.

**And `z0` remains a knob too**, measured at a dead-calm column (so these are lower bounds now that the wind is held up):

| `z0` (m) | surface − air (K) | sensible heat flux (W m⁻²) |
|---|---|---|
| 3.21e-5 (default) | 11.37 | 9.7 |
| 1e-3 | 10.04 | 14.9 |
| 1e-1 | 6.93 | 30.7 |

**Re-measure the `z0` row at the page's own wind** before quoting it — the table above was taken before the relaxation existed, and the two knobs are not independent.

**Files:**
- Create: `docs/modelling-tour/08-turbulent-heat-exchange.qmd`
- Modify: `docs/_quarto.yml`, `docs/modelling-tour/index.qmd`, `docs/modelling-tour/_artifacts/generate.py`
- Test: `tests/test_modelling_tour.py`

**Interfaces:**
- Consumes: `_tour/stepping.py`, `_tour/budgets.py`.
- Produces: nothing other pages import. Pages 09, 11, 12 assume its craft content: a `Stepper` in the same loop as tendency components; flux diagnostics; the `'bulk'`/`'external'`/`None` modes.

- [ ] **Step 1: Write the failing tests**

```python
PAGE8_Z0 = 1e-3
PAGE8_WIND = 5.0
PAGE8_STEPS = 12000        # dt = 1 h; every configuration converged by here


def _page8_column(z0=PAGE8_Z0, wind=PAGE8_WIND, slab_depth=2.0, nz=28,
                  boundary_layer=True):
    """Page 8's column: page 7's gray column, a boundary layer, and a wind.

    ``wind=None`` omits the momentum source, which is how the page shows what
    happens without one -- the column spins itself down to dead calm.
    """
    stepping_module = _load("stepping")
    longwave = climt.CorkLongwaveRadiation(
        optics="correlated_k", table=PAGE7_GRAY_TABLE,
        diffusivity_factor=PAGE7_DIFFUSIVITY)
    surface = climt.SlabSurface()
    steppers = ([climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                           roughness_length=z0)]
                if boundary_layer else [])
    state = climt.get_default_state([longwave, surface] + steppers,
                                    grid_state=get_grid(nx=1, ny=1, nz=nz))
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0

    tendencies = [longwave, surface]
    if boundary_layer and wind is not None:
        tendencies.append(stepping_module.wind_relaxation(state, wind))
    return tendencies, steppers, state


def _surface_air_jump(state):
    return (float(state["surface_temperature"].values.ravel()[0])
            - float(state["air_temperature"].values[0, 0, 0]))


@pytest.mark.slow
def test_page8_turbulence_shrinks_the_surface_air_discontinuity():
    """Page 8's reveal. Measured 15.15 K -> 7.50 K at z0=1e-3, 5 m/s.

    Thresholds sit either side of the measurement with room, so this checks
    that turbulence erodes the discontinuity, not that it erodes it by exactly
    7.65 K.
    """
    stepping_module = _load("stepping")
    timestep = climt.UnytTimeDelta(hours=1)

    jumps = {}
    for label, boundary_layer in (("radiative", False), ("turbulent", True)):
        tendencies, steppers, state = _page8_column(
            boundary_layer=boundary_layer)
        stepping_module.integrate(tendencies, steppers, state, timestep,
                                  PAGE8_STEPS)
        jumps[label] = _surface_air_jump(state)

    assert jumps["radiative"] > 13.0, (
        f"radiative-equilibrium jump {jumps['radiative']:.2f} K — page 4's "
        "discontinuity should be there before anything erodes it")
    assert jumps["turbulent"] < jumps["radiative"] - 4.0, (
        f"turbulent jump {jumps['turbulent']:.2f} K vs radiative "
        f"{jumps['radiative']:.2f} K — the boundary layer should erode it")


@pytest.mark.slow
def test_page8_a_column_with_no_momentum_source_spins_itself_down():
    """Why the page needs a wind forcing at all.

    Nothing drives a wind in a column with no dynamics, and the boundary
    layer's own surface drag removes any wind prescribed as an initial
    condition. Measured: three initialisations spanning 0-10 m/s all end at a
    lowest-level wind of 0.000 and the same equilibrium within 0.05 K.
    """
    stepping_module = _load("stepping")
    finals = {}
    for wind in (0.0, 5.0, 10.0):
        tendencies, steppers, state = _page8_column(wind=None)
        state["eastward_wind"].values[:] = wind
        stepping_module.integrate(tendencies, steppers, state,
                                  climt.UnytTimeDelta(hours=1), PAGE8_STEPS)
        finals[wind] = (_surface_air_jump(state),
                        float(state["eastward_wind"].values[0, 0, 0]))

    jumps = [jump for jump, _ in finals.values()]
    assert max(jumps) - min(jumps) < 0.5, (
        f"initial wind changed the equilibrium jump: {finals} — without a "
        "momentum source it must not, and page 8's argument for adding one "
        "depends on that")
    for _, lowest_wind in finals.values():
        assert abs(lowest_wind) < 0.1, (
            f"lowest-level wind {lowest_wind:.3f} m/s — the drag should have "
            "removed it entirely")


@pytest.mark.slow
@pytest.mark.parametrize("wind, min_flux, max_jump", [
    (2.0, 15.0, 10.0),
    (5.0, 28.0, 8.5),
    (10.0, 42.0, 7.0),
])
def test_page8_wind_speed_is_a_knob_once_it_is_held_up(wind, min_flux,
                                                       max_jump):
    """Page 8's knob, over the range the page tells readers it was tested on.

    Measured sensible heat fluxes: 19.9, 34.9, 49.9 W/m^2; jumps 9.31, 7.50,
    6.07 K. The thresholds allow ~25% either way — this checks the monotone
    response, not a fit.
    """
    tendencies, steppers, state = _page8_column(wind=wind)
    _load("stepping").integrate(tendencies, steppers, state,
                                climt.UnytTimeDelta(hours=1), PAGE8_STEPS)
    flux = float(state["surface_upward_sensible_heat_flux"].values.ravel()[0])
    assert flux > min_flux, f"wind={wind}: sensible heat flux {flux:.2f} W/m^2"
    assert _surface_air_jump(state) < max_jump


@pytest.mark.slow
@pytest.mark.parametrize("z0", [3.21e-5, 1e-3, 1e-1])
def test_page8_rougher_surface_moves_more_heat(z0):
    """The second knob, over the range the page tells readers it was tested on.

    Measured at a dead-calm column (before the wind relaxation existed): 9.7,
    14.9 and 30.7 W/m^2, jumps 11.37, 10.04 and 6.93 K. Those are lower bounds
    now that the wind is held up, and the two knobs are NOT independent — a
    rougher surface both exchanges more heat and exerts more drag. Re-measure
    at the page's own wind and set the thresholds from that; this asserts only
    the monotone ordering, which holds either way.
    """
    stepping_module = _load("stepping")
    fluxes = {}
    for roughness in (3.21e-5, z0):
        tendencies, steppers, state = _page8_column(z0=roughness)
        stepping_module.integrate(tendencies, steppers, state,
                                  climt.UnytTimeDelta(hours=1), PAGE8_STEPS)
        fluxes[roughness] = float(
            state["surface_upward_sensible_heat_flux"].values.ravel()[0])

    if z0 > 3.21e-5:
        assert fluxes[z0] > fluxes[3.21e-5], (
            f"z0={z0}: {fluxes[z0]:.2f} W/m^2 vs {fluxes[3.21e-5]:.2f} at the "
            "default — a rougher surface must move more heat")


@pytest.mark.slow
def test_page8_relaxation_lets_the_drag_win_near_the_ground():
    """The comparison that distinguishes a forcing from an assignment.

    Relaxing toward 5 m/s leaves the lowest level at ~2.9, because the surface
    drag is still acting and the relaxation only pulls. Assigning the wind back
    every step -- the idiom in examples/column_code_with_slab.py -- pins it at
    exactly 5.0 and so overrides the drag at the one level where the exchange
    happens, delivering more flux (42.9 vs 34.9 W/m^2) from the same nominal
    wind. The relaxation is the more honest scheme; the page says why.
    """
    stepping_module = _load("stepping")

    tendencies, steppers, state = _page8_column(wind=5.0)
    stepping_module.integrate(tendencies, steppers, state,
                              climt.UnytTimeDelta(hours=1), PAGE8_STEPS)
    relaxed_wind = float(state["eastward_wind"].values[0, 0, 0])
    relaxed_flux = float(
        state["surface_upward_sensible_heat_flux"].values.ravel()[0])

    tendencies, steppers, reset_state = _page8_column(wind=None)
    timestep = climt.UnytTimeDelta(hours=1)
    for _ in range(PAGE8_STEPS):
        stepping_module.integrate(tendencies, steppers, reset_state,
                                  timestep, 1)
        reset_state["eastward_wind"].values[:] = 5.0
    reset_flux = float(
        reset_state["surface_upward_sensible_heat_flux"].values.ravel()[0])

    assert 1.0 < relaxed_wind < 4.5, (
        f"relaxed lowest-level wind {relaxed_wind:.2f} m/s — it should sit "
        "well below the 5 m/s target, because the drag is still acting")
    assert reset_flux > relaxed_flux, (
        f"hard reset {reset_flux:.1f} vs relaxation {relaxed_flux:.1f} W/m^2 "
        "— pinning the lowest level overrides the drag where the exchange "
        "happens, so it must move more heat")


def test_page8_no_flux_mode_conserves_the_column():
    """The third mode: with surface_fluxes=None the diffusion conserves.

    Cheap, so unmarked: ten steps is enough, because conservation is exact
    rather than asymptotic.
    """
    budgets_module = _load("budgets")
    longwave = climt.CorkLongwaveRadiation(
        optics="correlated_k", table=PAGE7_GRAY_TABLE,
        diffusivity_factor=PAGE7_DIFFUSIVITY)
    boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes=None)
    state = climt.get_default_state([longwave, boundary_layer],
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    state["specific_humidity"].values[:] = 4e-3
    before = budgets_module.column_enthalpy(state)

    timestep = climt.UnytTimeDelta(hours=1)
    for _ in range(10):
        diagnostics, new_state = boundary_layer(state, timestep)
        state.update(new_state)
        state.update(diagnostics)

    after = budgets_module.column_enthalpy(state)
    assert after == pytest.approx(before, rel=1e-9), (
        "surface_fluxes=None must conserve every column integral exactly")
```

- [ ] **Step 2: Run them**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k page8 -v`
Expected: the unmarked conservation test passes at once; the four `slow` ones take ~4 minutes together (12 000 steps × ~2.9 ms × several runs) and pass against the measurements above. If `test_page8_no_flux_mode_conserves_the_column` fails, that is a library finding, not a test to loosen — the component's docstring promises exact conservation in that mode.

- [ ] **Step 3: Write the page**

`docs/modelling-tour/08-turbulent-heat-exchange.qmd`. Front matter as page 07's, with `resources:` listing `_tour/stepping.py` and `_tour/budgets.py` (no `soundings.py` — this page prescribes nothing).

**Cells:**

0. Setup, and the column builder taking `z0` and `wind`, as `_page8_column` above but written for a reader — with the `stepping.wind_relaxation(state, wind)` line visible, not tucked away.
1. **Headline figure**: the two converged profiles near the surface, side by side — radiative equilibrium and radiative-plus-turbulent — with the sensible heat flux time series beneath. Two panels, the upper one zoomed to the lowest ~200 hPa so the discontinuity is visible; label the jump in K on both. This is the cell `generate.py` captures.
2. **The column that stops moving.** Build it with *no* momentum source, initialise `eastward_wind` at 10 m s⁻¹, run it, and plot the wind profile at several times decaying to zero — plus the lowest-level wind as a time series. Print the final value: `0.000`. This is a cell and a figure, not a claim.
3. **Giving it a wind.** Introduce `stepping.wind_relaxation`, show `dx/dt = −(x − x_eq)/τ` beside it, print the component's `input_properties` so the reader sees `equilibrium_eastward_wind` and `eastward_wind_relaxation_timescale`, and print the two fields the helper wrote into the state. Re-run and plot the wind profile that now persists.
4. **The knob**: wind at 2, 5 and 10 m s⁻¹, plotting the jump and the sensible heat flux against it. Quote the measured table. Say the range was tested and that outside it is untested.
5. **Relaxation versus assignment.** Run the same column two ways: relaxed to 5 m s⁻¹, and hard-reset to 5 m s⁻¹ at the end of every step, the way `examples/column_code_with_slab.py` does. Print both lowest-level winds (2.88 vs 5.00) and both fluxes (34.9 vs 42.9 W m⁻²).
6. **`z0`, the second knob**, at 3.21e-5, 1e-3 and 1e-1, at the page's own wind.

**Prose:**

1. *Where this fits* — chapter 10, and the callback to page 04's discontinuity.
2. *climt craft: a Stepper in the loop* — `SimpleBoundaryLayer` returns `(diagnostics, new_state)`, not tendencies, so it is called directly inside the same loop rather than handed to `AdamsBashforth`. Point at `_tour/stepping.py`'s two lists.
3. *climt craft: diagnostics that are fluxes* — the scheme reports the flux it *applied*, evaluated post-solve, so the column budgets close to round-off. That is why the plotted flux is the model's flux and not a re-derivation.
4. *climt craft: the three surface-flux modes, and the double-counting trap* — `'bulk'` computes and applies its own; `'external'` consumes prescribed fluxes; `None` applies none and conserves exactly. Pairing `'bulk'` with a component that already applies surface fluxes applies them twice. Name `SimplePhysics(surface_fluxes=True)` as the classic case, noting it is a compiled component and absent here — the trap is real elsewhere.
5. **The page's second subject, and it needs three callouts of its own:**

   a. *Your column stops moving, and it takes its heat flux with it.* Cell 2's result. A single column has no dynamics, so nothing drives a wind, and the boundary layer's own surface drag removes the one you prescribed. Whatever wind you set is an initial condition, and initial conditions do not survive. A GCM's boundary layer is being run in a configuration that starves half of it — say plainly that this is what a column model is, and that noticing it is the skill.

   b. *climt craft: Newtonian relaxation, and what a forcing is.* `dx/dt = −(x − x_eq)/τ`. Two parameters, both physical: what you relax toward, and how fast. Show the climt idiom first — `state['eastward_wind'].values[:] = 3.` at the bottom of the loop, citing `examples/column_code_with_slab.py:120` — then `sympl.RelaxationTendencyComponent` as the same idea done as a component: it composes with everything else, carries units, and puts the assumption in the component list where a reader can see it rather than in a line at the bottom of a loop where they cannot. Note that τ has to be long compared with the timestep and short compared with the run, or it is a hard reset or no forcing at all.

   c. *A forcing is not an assignment.* Cell 5's comparison, which is the sharpest thing on this page. Relaxing toward 5 m s⁻¹ leaves the lowest level at 2.88, because the drag is still acting and the relaxation only pulls; the hard reset pins it at exactly 5.00 and overrides the drag at the one level where the exchange happens, so it moves more heat (42.9 vs 34.9 W m⁻²) from the same nominal wind. Neither is *wrong*; they are different physical statements, and the reader should be able to say which one they meant.

   Mention, in one sentence, that `_tour/stepping.py` subclasses `sympl.RelaxationTendencyComponent` because the stock version's tendency units come out of pint as `'1.0 meter / second ** 2'`, which `unyt` cannot parse — a small, real friction between two libraries, worth seeing once. Do not dwell on it.

6. *Knobs* — wind speed (2–10 m s⁻¹) and `z0` (3.21e-5 to 1e-1), both clearly marked with their tested ranges, and a note that they are not independent: a rougher surface both exchanges more heat and exerts more drag.
7. *Physics exercises* — (i) predict the sign of the surface temperature change when the boundary layer is switched on, and explain it from the surface energy budget before running; (ii) at 10 m s⁻¹ the jump is 6.1 K rather than 0 — what would have to change for it to close entirely?; (iii) the relaxation timescale is 24 h. Predict what happens at 1 h and at 100 days, then check.
8. *Code exercises* — (i) switch to `surface_fluxes=None` and verify column enthalpy is conserved with `budgets.column_enthalpy`; (ii) print `SimpleBoundaryLayer.diagnostic_properties` in each of the three modes and explain why the two flux entries move between inputs and diagnostics; (iii) relax `northward_wind` as well and check whether the flux responds to the wind *speed* or to the eastward component alone.
9. *Going deeper* — `../user-guide/` boundary layer material if it exists; otherwise the Frierson (2006) reference.

- [ ] **Step 4: Sidebar, chapter map, static fallback, render, commit**

Same six mechanical steps as Task 10, Steps 4–7, with:

```yaml
        - modelling-tour/08-turbulent-heat-exchange.qmd
```

```python
    "08-turbulent-heat-exchange.qmd": [(1, "08-turbulent.png"),
                                       (2, "08-spindown.png")],
```

Two figures: the headline comparison, and the column spinning itself down. The second earns its place — it is the evidence for the page's whole second half.

```markdown
| [Turbulent heat exchange](08-turbulent-heat-exchange.qmd) | 10 | bulk fluxes, and the discontinuity they erode |
```

The headline cell is 2 × 12 000 gray steps: ~70 s native JIT-on for `generate.py`, and ~4 minutes in the browser at ~11 ms/step. **State the expected browser time in the prose above the cell.** If that is judged too long, the lever is a thinner slab (which shortens the spin-up proportionally), not fewer steps to an unconverged state — an unconverged discontinuity is not the number the page is about.

```bash
git commit -m "docs(tour): page 8, the process that erodes page 4's discontinuity

A column with no dynamics spins its own wind down within days, so the page
also teaches the momentum source it needs: climt's assign-it-back idiom, then
sympl's RelaxationTendencyComponent as the cleaner form, and the difference
between a forcing and an assignment.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 12: Page 09 — Moisture changes the buoyancy

*Course chapter 10, continued.*

Turns on the latent flux by giving the slab a saturated surface humidity, and introduces the virtual potential temperature θ_v = θ(1 + 0.61q) as the quantity buoyancy actually depends on. Radiation switches to the 14-band table, so the water vapour the boundary layer supplies is radiatively active — the greenhouse effect tranche 1 measured, now supplied by the surface rather than typed in. `GridScaleCondensation` enters here as the minimum moisture sink: the boundary layer mixes moist air into colder air aloft, and without a sink the column supersaturates.

**Numbers this page needs come from Task 8's `supersat` measurement.** Do not write the page before that log is filled in.

**This page carries page 08's wind relaxation, and so does every page after it.** A column with no momentum source spins itself down within days, and the latent flux this page is about is a bulk flux — it scales with the wind. Without the forcing, page 09's Bowen ratio is that of a dead-calm planet. `_page9_column` and the page's own cells build it the same way page 08 did, and the page says in one sentence that it is inheriting the fix rather than re-teaching it.

**A spec assumption to check rather than inherit.** The spec's first claim is *"the latent flux exceeds the sensible flux over a saturated ~288 K surface"*. This page's column is not at 288 K — page 07's gray column equilibrates near 335 K, and the 14-band column with the water vapour this page supplies will land somewhere else again. Measure the converged surface temperature first, then write the claim at that temperature. If the Bowen ratio does not come out below 1, that is a finding to record and reconcile, not a claim to soften.

**Files:**
- Create: `docs/modelling-tour/09-moisture-and-buoyancy.qmd`
- Modify: `docs/_quarto.yml`, `docs/modelling-tour/index.qmd`, `docs/modelling-tour/_artifacts/generate.py`
- Test: `tests/test_modelling_tour.py`

**Interfaces:**
- Consumes: `_tour/stepping.py`, `_tour/budgets.py`.
- Produces: nothing other pages import. Pages 11 and 12 assume its content: θ_v, the Bowen ratio, and that `GridScaleCondensation` is in the stack from here on.

- [ ] **Step 1: Write the failing tests**

```python
PAGE9_DT_HOURS = 1.0
PAGE9_STEPS = 4000       # replace with Task 8's converged step count


def _page9_column(surface_relative_humidity=1.0, nz=28, slab_depth=2.0,
                  condensation=True, wind=PAGE8_WIND):
    """Page 9's column: 14-band radiation, boundary layer, moist surface.

    Carries page 8's wind relaxation. Every configuration from page 8 on does:
    a column with no momentum source spins itself down to dead calm, and its
    surface fluxes — including the latent flux this page is about — are then
    those of a windless planet.
    """
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table="earth_low_res_lw")
    surface = climt.SlabSurface()
    steppers = [climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                          roughness_length=PAGE8_Z0)]
    if condensation:
        steppers.append(climt.GridScaleCondensation())
    state = climt.get_default_state([longwave, surface] + steppers,
                                    grid_state=get_grid(nx=1, ny=1, nz=nz))
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["surface_specific_humidity"].values[:] = (
        surface_relative_humidity * 0.015)
    relaxation = _load("stepping").wind_relaxation(state, wind)
    return [longwave, surface, relaxation], steppers, state


def _virtual_potential_temperature(state):
    """theta_v = theta * (1 + 0.61 q), from state quantities."""
    from sympl import get_constant

    Rd = float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
    Cp = float(get_constant("heat_capacity_of_dry_air_at_constant_pressure",
                            "J/kg/degK"))
    p = state["air_pressure"].values[:, 0, 0]
    T = state["air_temperature"].values[:, 0, 0]
    q = state["specific_humidity"].values[:, 0, 0]
    theta = T * (1.0e5 / p) ** (Rd / Cp)
    return theta, theta * (1.0 + 0.61 * q)


def test_page9_virtual_potential_temperature_identity():
    """theta_v - theta = 0.61 q theta, to machine precision.

    Cheap and exact: this is the definition, and the page derives it in a cell.
    If a reader's arithmetic disagrees with the page, it is the page that is
    wrong, so pin it.
    """
    tendencies, steppers, state = _page9_column()
    state["specific_humidity"].values[:] = 8e-3
    theta, theta_v = _virtual_potential_temperature(state)
    np.testing.assert_allclose(theta_v - theta, 0.61 * 8e-3 * theta, rtol=1e-12)


@pytest.mark.slow
def test_page9_latent_flux_exceeds_sensible_over_a_saturated_surface():
    """Page 9's headline claim, at this page's own converged temperature.

    The spec phrases it "over a saturated ~288 K surface"; this column is not
    at 288 K. The claim is about the partition, so it is asserted on the
    Bowen ratio at whatever temperature the column reaches, and the page
    quotes that temperature beside it.
    """
    tendencies, steppers, state = _page9_column()
    _load("stepping").integrate(tendencies, steppers, state,
                                climt.UnytTimeDelta(hours=PAGE9_DT_HOURS),
                                PAGE9_STEPS)
    sensible = float(
        state["surface_upward_sensible_heat_flux"].values.ravel()[0])
    latent = float(state["surface_upward_latent_heat_flux"].values.ravel()[0])
    surface = float(state["surface_temperature"].values.ravel()[0])

    assert latent > sensible, (
        f"Bowen ratio {sensible / latent:.2f} at surface {surface:.1f} K — "
        "over a saturated surface the latent flux should dominate")


@pytest.mark.slow
def test_page9_condensation_removes_the_supersaturation():
    """The reason GridScaleCondensation is in the stack from this page on.

    Without a sink the boundary layer mixes moist air into colder air aloft
    and the column supersaturates. Fill the two peak relative humidities in
    from Task 8's `supersat` measurement before tightening this.
    """
    stepping_module = _load("stepping")
    peaks = {}
    for label, condensation in (("with", True), ("without", False)):
        tendencies, steppers, state = _page9_column(condensation=condensation)
        stepping_module.integrate(tendencies, steppers, state,
                                  climt.UnytTimeDelta(hours=PAGE9_DT_HOURS),
                                  PAGE9_STEPS)
        peaks[label] = _peak_relative_humidity(state)

    assert peaks["without"] > peaks["with"], (
        f"peak RH with condensation {peaks['with']:.3f}, without "
        f"{peaks['without']:.3f} — the sink must reduce it")
    assert peaks["with"] < 1.05, (
        f"peak RH {peaks['with']:.3f} with condensation active — the sink is "
        "not keeping up")


@pytest.mark.slow
def test_page9_surface_relative_humidity_moves_the_bowen_ratio():
    """Page 9's knob, over the range the page tells readers it was tested on."""
    stepping_module = _load("stepping")
    ratios = {}
    for relative_humidity in (0.4, 1.0):
        tendencies, steppers, state = _page9_column(
            surface_relative_humidity=relative_humidity)
        stepping_module.integrate(tendencies, steppers, state,
                                  climt.UnytTimeDelta(hours=PAGE9_DT_HOURS),
                                  PAGE9_STEPS)
        sensible = float(
            state["surface_upward_sensible_heat_flux"].values.ravel()[0])
        latent = float(
            state["surface_upward_latent_heat_flux"].values.ravel()[0])
        ratios[relative_humidity] = sensible / max(latent, 1e-9)

    assert ratios[0.4] > ratios[1.0], (
        f"Bowen ratio {ratios[0.4]:.2f} at RH 0.4 should exceed "
        f"{ratios[1.0]:.2f} at RH 1.0 — a drier surface partitions more into "
        "sensible heat")
```

Add the helper `_peak_relative_humidity(state)` beside the others, using the same Bolton form `_tour/soundings.saturation_vapour_pressure` uses — import it from `soundings` rather than restating it, so the page, the test and the measurement script all agree:

```python
def _peak_relative_humidity(state):
    soundings_module = _load("soundings")
    from sympl import get_constant

    epsilon = (float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
               / float(get_constant("gas_constant_of_vapor_phase",
                                    "J/kg/degK")))
    T = state["air_temperature"].values[:, 0, 0]
    q = state["specific_humidity"].values[:, 0, 0]
    p = state["air_pressure"].values[:, 0, 0]
    e_sat = soundings_module.saturation_vapour_pressure(T)
    q_sat = epsilon * e_sat / np.maximum(p - (1.0 - epsilon) * e_sat, 1.0)
    return float(np.max(q / q_sat))
```

- [ ] **Step 2: Run them, and record what converged**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k page9 -v`

Before trusting the `slow` ones, check `PAGE9_STEPS` is actually enough:

```bash
NUMBA_DISABLE_JIT=1 conda run -n climt python -c "
import sys; sys.path.insert(0, 'scripts/experiments')
import sympl, climt, tour_rce_measurements as m
sympl.set_backend(climt.UnytBackend())
m.convergence(dt_hours=1.0)" 2>&1 | grep 09-
```

Set `PAGE9_STEPS` from that, delete the `# replace with` comment, and record the converged surface temperature and Bowen ratio in the log below — the page quotes both.

- [ ] **Step 3: Write the page**

`docs/modelling-tour/09-moisture-and-buoyancy.qmd`. Front matter as page 08's.

**Cells:**

0. Setup and the column builder, taking `surface_relative_humidity`.
1. **Headline figure**: three panels — θ and θ_v side by side (the gap between them *is* the moisture's contribution to buoyancy), the specific humidity profile, and the sensible/latent flux time series. Captured by `generate.py`.
2. **θ_v derived from state quantities**, in code, with units handled explicitly. Print `θ_v − θ` and `0.61 q θ` next to each other so the identity is visible rather than asserted.
3. **The knob**: surface relative humidity at 0.4 and 1.0, plotting the Bowen ratio.
4. **What condensation removes**: run with and without `GridScaleCondensation`, print both peak relative humidities and the precipitation rate.

**Prose:**

1. *Where this fits* — chapter 10 continued; the callback to tranche 1's page 06, which measured the water vapour greenhouse on a profile the reader typed in. Here the surface supplies it.
2. *Why θ_v and not θ* — buoyancy depends on density, density depends on the mixture, and 0.61 is `Rv/Rd − 1`. A moist parcel is lighter than a dry one at the same temperature.
3. *climt craft: components that write the same variable* — the boundary layer and condensation both write `specific_humidity`. The order in the stepper list is a physical choice: mixing first then condensing means the sink sees the mixed profile. Say which order this page uses and why, and note that the reverse is defensible and gives a different answer.
3b. *One sentence, not a section: the wind relaxation is still here.* Page 08 introduced it; every page from here on carries it, because the latent flux is a bulk flux and scales with the wind. Point back at page 08 and move on.
4. *climt craft: deriving a quantity from state with correct units* — show the `(1e5/p)**(Rd/Cp)` line and where `Rd`, `Cp` come from (`sympl.get_constant`), not hardcoded.
5. *Why condensation is in the stack from here on* — the boundary layer mixes moist air into colder air aloft; without a sink the column supersaturates, and quote the two peak relative humidities from cell 4.
6. *Knob* — surface relative humidity, tested range 0.4–1.0.
7. *Physics exercises* — (i) at what q does θ_v − θ reach 1 K, and where in the column is that?; (ii) predict the Bowen ratio's direction as the surface dries, then run it.
8. *Code exercises* — (i) swap the order of the two moisture-writing components and quantify the difference; (ii) remove `GridScaleCondensation` and plot relative humidity against height.
9. *Going deeper* — tranche 1's `06-water-vapour-limit.qmd`.

- [ ] **Step 4: Sidebar, chapter map, static fallback, render, commit**

```yaml
        - modelling-tour/09-moisture-and-buoyancy.qmd
```
```python
    "09-moisture-and-buoyancy.qmd": [(1, "09-buoyancy.png")],
```
```markdown
| [Moisture changes the buoyancy](09-moisture-and-buoyancy.qmd) | 10 | θ_v, the Bowen ratio, and why a sink is needed |
```

```bash
git commit -m "docs(tour): page 9, the surface supplies its own greenhouse gas

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

### Log — page 9's column as converged

Measured 2026-09-26 with `scripts/experiments/tour_page9_measurements.py` (Linux, 4 cores, numba on for physics, `NUMBA_NUM_THREADS=1`; `cost` with `NUMBA_DISABLE_JIT=1`). Column exactly as the page builds it: `earth_low_res_lw` (default diffusivity), nz 28, 2 m slab, SOLAR 240, `SimpleBoundaryLayer(bulk, z0=1e-3)`, `wind_relaxation` 5 m/s / 24 h, `GridScaleCondensation` last. Every flux is a 30-day mean (converged) or a 10-day mean (the page's 30-day cells).

**Deviation from Step 1's `_page9_column`: the surface humidity.** The plan wrote `surface_specific_humidity = RH * 0.015`, a fixed number that nothing updates. At the 289.4 K this column settles at, 0.015 kg/kg is **131 % of saturation**. The spec's "saturated surface" and the "surface relative humidity" knob are then neither saturated nor a relative humidity. Measured with the fixed 0.015: Ts 289.39 K, SH 31.7, LH 92.8, Bowen **0.34** (surface RH 1.32). The page and tests instead use a new `stepping.SurfaceHumidity(rh)` stepper, first in the stepper list, that sets `q_s = rh * q_sat(Ts, ps)` every step with `soundings.saturation_specific_humidity`, the formula `GridScaleCondensation` condenses to. That is also what `SimplePhysics` does internally. All numbers below use it.

**Converged cold start, dt 1 h, 24 000 steps (1000 d), surface RH 1.0:**

```
  day 1000   Ts 289.37 K  SH 44.7  LH 63.5  Bowen 0.704  jump 7.46 K  BL depth median 723 m (mean 1703)
             q0 6.98 g/kg  theta_v - theta (lowest level) 1.20 K  peak RH 100.0 %  P 2.19 mm/d  TOA -0.14
  day  900   Ts 289.40 K  Bowen 0.702  TOA -0.22      day 800  Ts 289.47  TOA -0.37
  day  700   Ts 289.59 K  Bowen 0.694  TOA -0.63      day 30   Ts 293.65 (30-d mean)  TOA -117
```
Evaporation 63.5 W m⁻² / Lv = 2.19 mm/day = precipitation. Converges in ~900 days (21 600 steps at 1 h): over an hour in a browser at ~0.18 s/step. So the page runs **30-day cold starts live** (spec: "yes, live, short") and quotes the converged numbers from the script. `PAGE9_STEPS = 21600`.

**The spec's claim, checked at the column's own temperature:** latent exceeds sensible over a saturated surface at **289.37 K**, Bowen **0.70**. Holds, and the column is near the spec's ~288 K. (Well above Earth's ocean-mean ~0.1–0.2: this column's lowest air is saturated.)

**Timestep (converged, RH 1.0):**

```
  dt 30 min  Ts 289.33  SH 44.1  LH 63.3  Bowen 0.696
  dt  1 h    Ts 289.37  SH 44.7  LH 63.5  Bowen 0.704
  dt  3 h    Ts 289.42  SH 46.8  LH 64.1  Bowen 0.730
  dt  6 h    Ts 289.41  SH 49.3  LH 64.5  Bowen 0.764
```
Surface temperature is dt-independent; the fluxes are not (as on page 8). dt = 1 h.

**Knob and comparisons (converged, 1000 d):**

```
  RH 0.4      Ts 288.32 K  SH 67.1  LH 16.2  Bowen 4.15   BL median 993 m
  RH 0 (dry)  Ts 266.58 K  SH 17.2  LH  0.0               BL median 445 m
  condense first (order swapped)  Ts 289.48  Bowen 0.709  end-of-step peak RH 103.9 % (max 127.6 %)
```
The drier surface is 1.05 K *colder* converged (less vapour, less greenhouse), though 1.9 K warmer at day 30. Moist minus dry surface: **22.8 K**, the water-vapour greenhouse the surface supplied. Boundary-layer top: 723 m moist vs 445 m dry (medians); page 8's gray column sat near 550 m.

**The page's cells (30-day cold starts at 1 h; 10-day means):**

```
  headline RH 1.0   Ts 291.14  SH 46.6  LH 72.6  Bowen 0.64  TOA -93.8  P 2.44   (lowest level at the last step: q 7.30 g/kg, theta_v - theta 1.25 K)
  knob RH 0.4       Ts 293.06  SH 69.8  LH 18.6  Bowen 3.75
  condensation cell, from the headline state, +10 days:
    with       peak RH 100.0 %                          P 1.89 mm/d  LH 64.8  Bowen 0.69
    without    peak RH 219 % after 1 d, 707 % after 10 d, at 836 hPa   LH 36.0  Bowen 2.12
```

**Supersaturation (replaces Task 0 §6).** Task 0's `supersat` (457 % vs 100 %, Bowen 0.64) was 200 steps from a cold start with the fixed 0.015 surface and a loop of one-step `integrate` calls, which is forward Euler. Its 457 % is superseded by the 707 % above. With no sink at all from a cold start, the column has no equilibrium: lowest-level q 16.8 g/kg at day 30, 24.4 at day 60, 44.0 at day 100, 266 at day 360, surface 298 → 305 → 315 → 348 K, and the radiation returns non-finite fluxes after **8881 steps (370 d)**.

**Cost:** 34.7 ms/step native, JIT off, on this machine (Task 0's reference machine: 51.4 ms). Page cells, JIT off, here: headline 25.2 s, knob 25.2 s, condensation 19.0 s. Browser at 0.18 s/step: about 2 min, 2 min, 1.5 min.

---

## Task 13: Page 10 — Dry convection

*Course chapter 11.*

**The only page in the tranche with no time loop.** Prescribe a superadiabatic profile, call `DryConvectiveAdjustment` once, show before and after in both T and θ. Static stability and enthalpy-conserving mixing seen on their own, with no convergence noise in the way — the same single-call discipline that made tranche 1 legible.

**Verified before this plan was written.** A 28-level column with Γ = 14 K/km through the lowest 8 layers over a 6.5 K/km profile above, `specific_humidity = 0`, one call at `dt = 1 h`:

| | value |
|---|---|
| levels adjusted | 8 |
| θ spread across the adjusted layers, after | **0.0 K** (exactly) |
| column enthalpy, relative change | **1.8 × 10⁻¹⁶** |
| max \|ΔT\| | 1.94 K |
| θ before, layers 0–7 | 298.83 → 295.02 K (decreasing upward: unstable) |
| θ after, layers 0–7 | 296.89 K, every one |

Those are page 10's numbers, and they are as clean as this tranche gets.

**Files:**
- Create: `docs/modelling-tour/10-dry-convection.qmd`
- Modify: `docs/_quarto.yml`, `docs/modelling-tour/index.qmd`, `docs/modelling-tour/_artifacts/generate.py`
- Test: `tests/test_modelling_tour.py`

**Interfaces:**
- Consumes: `_tour/budgets.py` (`column_enthalpy`), `_tour/soundings.py`.
- Produces: nothing other pages import. Pages 11 and 12 assume its content: static stability, θ, and what adjustment does.

- [ ] **Step 1: Write the failing tests**

```python
def _superadiabatic_column(nz=28, unstable_levels=8, gamma_unstable=14e-3,
                           gamma=6.5e-3, T_surf=300.0):
    """Page 10's prescribed profile: superadiabatic near the ground."""
    from sympl import get_constant

    adjustment = climt.DryConvectiveAdjustment()
    state = climt.get_default_state([adjustment],
                                    grid_state=get_grid(nx=1, ny=1, nz=nz))
    Rd = float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
    g = float(get_constant("gravitational_acceleration", "m/s^2"))
    p = state["air_pressure"].values[:, 0, 0]
    surface_pressure = float(state["surface_air_pressure"].values.ravel()[0])
    z = -(Rd * 250.0 / g) * np.log(p / surface_pressure)

    T = T_surf - gamma * z
    T[:unstable_levels] = T_surf - gamma_unstable * z[:unstable_levels]
    state["air_temperature"].values[:, 0, 0] = np.maximum(T, 200.0)
    state["specific_humidity"].values[:] = 0.0
    return adjustment, state


def _potential_temperature(state):
    from sympl import get_constant

    Rd = float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
    Cp = float(get_constant("heat_capacity_of_dry_air_at_constant_pressure",
                            "J/kg/degK"))
    p = state["air_pressure"].values[:, 0, 0]
    T = state["air_temperature"].values[:, 0, 0]
    return T * (1.0e5 / p) ** (Rd / Cp)


def test_page10_adjustment_makes_theta_uniform(budgets):
    """Page 10's first claim. Measured: exactly 0.0 K spread."""
    adjustment, state = _superadiabatic_column()
    before = state["air_temperature"].values[:, 0, 0].copy()

    diagnostics, new_state = adjustment(state, climt.UnytTimeDelta(hours=1))
    state.update(new_state)

    changed = np.abs(state["air_temperature"].values[:, 0, 0] - before) > 1e-6
    assert changed.sum() >= 6, (
        f"only {changed.sum()} levels adjusted — the prescribed instability "
        "should reach at least six")
    theta = _potential_temperature(state)
    spread = float(theta[changed].max() - theta[changed].min())
    assert spread < 1e-6, (
        f"theta spread {spread:.2e} K across the adjusted layers — an "
        "adjusted layer is by definition neutrally stratified")


def test_page10_adjustment_conserves_column_enthalpy(budgets):
    """Page 10's second claim, and the test of whether you understood it.

    Measured relative change: 1.8e-16, which is machine precision on a sum of
    28 terms. tests/test_conservation.py::TestDryConvectionConservation asserts
    the same thing through a different route; this one asserts it on page 10's
    own profile, with budgets.column_enthalpy — the function the page calls.
    """
    adjustment, state = _superadiabatic_column()
    before = budgets.column_enthalpy(state)

    diagnostics, new_state = adjustment(state, climt.UnytTimeDelta(hours=1))
    state.update(new_state)

    after = budgets.column_enthalpy(state)
    assert after == pytest.approx(before, rel=1e-12), (
        f"enthalpy changed by {abs(after - before) / before:.2e} relative — "
        "dry adjustment mixes, it does not heat")


def test_page10_a_stable_column_is_left_alone():
    """The control. Nothing to adjust means nothing adjusted."""
    adjustment, state = _superadiabatic_column(gamma_unstable=6.5e-3)
    before = state["air_temperature"].values.copy()
    diagnostics, new_state = adjustment(state, climt.UnytTimeDelta(hours=1))
    np.testing.assert_allclose(new_state["air_temperature"].values, before)


def test_page10_adjustment_returns_a_state_not_tendencies():
    """The craft claim: a Stepper's signature is different, and visibly so."""
    adjustment, state = _superadiabatic_column()
    result = adjustment(state, climt.UnytTimeDelta(hours=1))
    assert isinstance(result, tuple) and len(result) == 2
    diagnostics, new_state = result
    assert "air_temperature" in new_state
    assert not any("tendency" in key for key in new_state)
```

- [ ] **Step 2: Run them**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k page10 -v`
Expected: all four pass in under a second — there is no integration here. None is marked `slow`.

- [ ] **Step 3: Write the page**

`docs/modelling-tour/10-dry-convection.qmd`. Front matter with `resources:` listing `_tour/budgets.py` only — no stepping, because there is no time loop, and saying so in the front matter is itself part of the page's point.

**Cells:**

0. Setup, and the prescribed superadiabatic profile, built visibly.
1. **Headline figure**: four panels — T before/after and θ before/after, on log pressure. The θ panel is where the lesson is: a decreasing-with-height θ becomes a vertical line. Captured by `generate.py`.
2. **The conservation check**, written as code: `budgets.column_enthalpy` before and after, printed with the relative difference. The page frames this as *the test of whether you understood the scheme*, not as a formality.
3. **The knob**: depth and strength of the initial instability — sweep `unstable_levels` in (4, 8, 12) and `gamma_unstable` in (10, 14, 20) K/km, plotting how far the adjusted layer reaches.

**Prose:**

1. *Where this fits* — chapter 11. And a note on why this page has no time loop: every other page in the tranche integrates, and this one deliberately does not, because static stability is a property of a profile and putting a convergence history in front of it would only add noise.
2. *climt craft: calling a Stepper directly* — it takes a timestep and returns `(diagnostics, new_state)`. No tendencies anywhere, and no `AdamsBashforth`. Contrast with page 07's `TendencyComponent`. Note that the timestep is required by the signature but that adjustment is instantaneous — the scheme does not use it.
3. *Why θ and not T* — a column is stable when θ increases with height, whatever T is doing. This is the whole content of the figure.
4. *Enthalpy is conserved to machine precision* — quote the 1.8 × 10⁻¹⁶, and say what it means: the scheme moves heat, it does not create it. Point at `tests/test_conservation.py::TestDryConvectionConservation`, which asserts this independently, so the page's claim is guarded by a test that existed before the page did.
5. *Knob* — depth and strength of the instability.
6. *Physics exercises* — (i) predict the adjusted θ from the initial profile alone, by conserving enthalpy, before running the scheme; (ii) make the instability deeper than the troposphere and explain what stops the adjustment.
7. *Code exercises* — (i) write the conservation check yourself rather than calling `budgets.column_enthalpy`, and compare; (ii) set `specific_humidity` non-zero and re-run — the scheme uses a moisture-dependent heat capacity, so what changes and by how much?
8. *Going deeper* — page 11, where this scheme joins radiation.

- [ ] **Step 4: Sidebar, chapter map, static fallback, render, commit**

```yaml
        - modelling-tour/10-dry-convection.qmd
```
```python
    "10-dry-convection.qmd": [(1, "10-dry-convection.png")],
```
```markdown
| [Dry convection](10-dry-convection.qmd) | 11 | enthalpy-conserving mixing to a dry adiabat |
```

This page's cells are milliseconds; the browser cost is the wheel download and nothing else.

```bash
git commit -m "docs(tour): page 10, one call and a conservation law

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 14: Page 11 — Dry RCE

*Course chapter 12.*

Everything dry, together, to equilibrium: longwave radiation, slab surface, boundary layer, dry convective adjustment. The page **loads the shipped equilibrium state** from Task 9 rather than spinning one up, then runs perturbation experiments from it.

**The reveal.** The convecting layer sits on a ~9.8 K/km dry adiabat — not the 6.5 K/km every tranche 1 page assumed — and the tropopause **emerges** as the level where convection stops, instead of being prescribed as a 200 K cap. Compared against page 07's radiative equilibrium: a much cooler surface, a similar upper atmosphere.

**Two things this page must say plainly.** Its column is genuinely dry, so its 14-band radiation sees CO₂ alone and it equilibrates well below Earth's surface temperature; **no number on this page is comparable with tranche 1's**. And it is perturbing an equilibrium computed for one specific configuration — `states.describe(provenance)` printed above the first figure is how the reader knows which.

**Files:**
- Create: `docs/modelling-tour/11-dry-rce.qmd`
- Modify: `docs/_quarto.yml`, `docs/modelling-tour/index.qmd`, `docs/modelling-tour/_artifacts/generate.py`
- Test: `tests/test_modelling_tour.py`

**Interfaces:**
- Consumes: `_tour/states.py`, `_tour/stepping.py`, `_tour/budgets.py`, `_tour/assets.py`, `docs/modelling-tour/_data/rce_dry_equilibrium.npz`, and `scripts/generate_tour_equilibria.dry_components` (via the test only — the page rebuilds its own list, visibly, because a reader has to see it).
- Produces: nothing other pages import. Page 12 assumes its content: component ordering, the budgets as convergence test, loading and saving a state.

- [ ] **Step 1: Write the failing tests**

```python
DRY_ADIABAT_K_PER_KM = 9.8


def _lapse_rate_profile(state):
    """dT/dz in K/km at mid levels, from the hydrostatic thickness."""
    from sympl import get_constant

    Rd = float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
    g = float(get_constant("gravitational_acceleration", "m/s^2"))
    T = state["air_temperature"].values[:, 0, 0]
    p = state["air_pressure"].values[:, 0, 0]
    # dz between mid levels from the hypsometric relation on the layer mean.
    T_mean = 0.5 * (T[:-1] + T[1:])
    dz = (Rd * T_mean / g) * np.log(p[:-1] / p[1:])
    return -np.diff(T) / dz * 1000.0


@pytest.mark.slow
def test_page11_convecting_layer_is_near_dry_adiabatic(states):
    """Page 11's reveal, and the number tranche 1 assumed away.

    Only the convecting layer is checked: above the tropopause the profile is
    radiative and steeply stable, and averaging that in would hide the claim.
    The convecting layer is identified as the levels dry adjustment actually
    touches, by calling it once on the loaded state.
    """
    generator = _generator()
    components = generator.dry_components()
    state, provenance = states.load(
        str(DATA / "rce_dry_equilibrium.npz"), components,
        grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))

    adjustment = [c for c in components
                  if isinstance(c, climt.DryConvectiveAdjustment)][0]
    before = state["air_temperature"].values[:, 0, 0].copy()
    _, adjusted = adjustment(state, climt.UnytTimeDelta(hours=1))
    touched = np.abs(adjusted["air_temperature"].values[:, 0, 0] - before) > 1e-8

    lapse = _lapse_rate_profile(state)
    # `touched` is per mid level; the lapse rate lives between them.
    convecting = touched[:-1] & touched[1:]
    assert convecting.sum() >= 3, (
        "dry adjustment touched fewer than four levels in the shipped "
        "equilibrium — there is no convecting layer to make a claim about")

    mean_lapse = float(np.mean(lapse[convecting]))
    assert abs(mean_lapse - DRY_ADIABAT_K_PER_KM) < 1.0, (
        f"convecting-layer lapse rate {mean_lapse:.2f} K/km is not within "
        f"1 K/km of dry adiabatic ({DRY_ADIABAT_K_PER_KM})")
    assert mean_lapse > 8.0, (
        f"lapse rate {mean_lapse:.2f} K/km — page 11's whole point is that it "
        "is NOT the 6.5 K/km tranche 1 prescribed")


@pytest.mark.slow
def test_page11_tropopause_emerges_rather_than_being_prescribed(states):
    """Above the convecting layer the profile is stable, and nothing capped it.

    Tranche 1 prescribed a 200 K isothermal stratosphere. Here the top of
    convection is wherever the scheme stops, and the air above it is whatever
    radiation makes it.
    """
    generator = _generator()
    components = generator.dry_components()
    state, _ = states.load(str(DATA / "rce_dry_equilibrium.npz"), components,
                           grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))
    T = state["air_temperature"].values[:, 0, 0]

    assert not np.any(np.isclose(T, 200.0, atol=0.05)), (
        "a level sitting at exactly 200.0 K suggests a prescribed cap leaked "
        "in from tranche 1's soundings")
    lapse = _lapse_rate_profile(state)
    assert lapse[-1] < lapse[0], (
        "the upper column should be less steeply lapsing than the convecting "
        "layer below it")


@pytest.mark.slow
def test_page11_co2_doubling_warms_the_surface_by_a_measured_amount(states):
    """Page 11's knob: a *measured* dry climate sensitivity.

    Tranche 1's page 5 could only compute the forcing. This integrates to the
    new equilibrium and reads the warming off. Fill the expected range in from
    Task 8's `co2` measurement before tightening.
    """
    generator = _generator()
    components = generator.dry_components()
    tendencies, steppers = generator.split(components)
    state, provenance = states.load(
        str(DATA / "rce_dry_equilibrium.npz"), components,
        grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))

    stepping_module = _load("stepping")
    tendencies = tendencies + [stepping_module.wind_relaxation(
        state, provenance["wind_m_s"], provenance["wind_timescale_hours"],
        initialise=False)]

    before = float(state["surface_temperature"].values.ravel()[0])
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] *= 2.0
    stepping_module.integrate(
        tendencies, steppers, state,
        climt.UnytTimeDelta(hours=provenance["dt_hours"]), PAGE11_2XCO2_STEPS)
    after = float(state["surface_temperature"].values.ravel()[0])

    warming = after - before
    assert 0.3 < warming < 6.0, (
        f"2xCO2 dry surface warming {warming:+.2f} K is outside any "
        "defensible range — check that the run reached equilibrium")
    imbalance = _load("budgets").toa_imbalance(state)
    assert abs(imbalance) < 0.5, (
        f"TOA imbalance {imbalance:+.3f} W/m^2 — the perturbed run has not "
        f"equilibrated in {PAGE11_2XCO2_STEPS} steps, so the warming is a "
        "lower bound, not a sensitivity")
```

Define `PAGE11_2XCO2_STEPS = 1000` beside the other page constants, and tighten the `0.3 < warming < 6.0` bracket to **`0.72 < warming < 1.68`** (±40% of 1.20 K).

**Use the re-measurement, not Task 8's `co2` number.** Task 8 measured +0.91 K in 175 steps from its own 1650-step dry run, which the Task 9 log shows stopped partway through settling. Re-measured from the shipped `rce_dry_equilibrium.npz` under the generator's own gate (`scripts/experiments/tour_rce_shipped_remeasure.py co2-dry`, 2026-09-26):

```
  step    500  TOA +0.183  Tsurf 267.632 K
  step   1000  TOA -0.028  Tsurf 267.676 K
  step   2350  TOA -0.042  Tsurf 267.678 K   <- the strict gate passes
  dry 2xCO2: 266.478 -> 267.678 K (+1.200 K)
```

The strict gate cannot pass before 2000 steps (it compares two 1000-step windows), but the warming has arrived by step 1000: +1.197 K, 99.8% of the final value, with |TOA| < 0.05. So the page runs **1000 steps live** (~3 min at ~0.18 s/step) and reads convergence off `budgets.summary`, as cell 3 already says. Unlike page 12's, this experiment does not need to ship.

- [ ] **Step 2: Run them**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k page11 -m slow -v`
Expected: the first two pass in seconds (they only load and inspect); the third is the long one.

If `test_page11_convecting_layer_is_near_dry_adiabatic` fails with a lapse rate near 6.5, check whether the shipped state was generated with `specific_humidity` non-zero — page 11's column is supposed to be dry, and moisture in it would relax the lapse rate toward moist adiabatic, which is page 12's result arriving a page early.

- [ ] **Step 3: Write the page**

`docs/modelling-tour/11-dry-rce.qmd`. Front matter's `resources:` block lists `_tour/assets.py`, `_tour/states.py`, `_tour/stepping.py`, `_tour/budgets.py` **and the data asset**:

```yaml
  resources:
    - _tour/assets.py
    - _tour/states.py
    - _tour/stepping.py
    - _tour/budgets.py
    # The shipped dry equilibrium. Staged into the Pyodide filesystem at the
    # same relative path it has on disk, which is why _tour/assets.py can find
    # it with a plain os.path.isfile in both environments.
    - _data/rce_dry_equilibrium.npz
```

**Cells:**

0. Setup; build the component list **visibly**, in the order the page defends, and load the state:

```python
components = [
    CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw"),
    SlabSurface(),
    SimpleBoundaryLayer(surface_fluxes="bulk", roughness_length=1e-3),
    DryConvectiveAdjustment(),
]
tendency = [c for c in components if not hasattr(c, "output_properties")]
stepper = [c for c in components if hasattr(c, "output_properties")]

state, provenance = states.load("rce_dry_equilibrium.npz", components,
                                grid_state=get_grid(nx=1, ny=1, nz=28))

# Page 8's momentum source, rebuilt on the loaded state. It has to come after
# the state exists, because it writes its target fields into it -- which is
# also why it is not in `components` above and not in the shipped file's
# component list. The provenance records the wind it was spun up at; use the
# same one, or you are perturbing a different model.
# initialise=False: the loaded state already carries a spun-up, drag-sheared
# wind profile. Flattening it to a uniform 5 m/s would discard part of the
# equilibrium we just loaded.
tendency.append(stepping.wind_relaxation(state, provenance["wind_m_s"],
                                         provenance["wind_timescale_hours"],
                                         initialise=False))

print(states.describe(provenance))
print(budgets.summary(state))
```

1. **Headline figure**: the equilibrium profile with the dry adiabat drawn over it, the convecting layer shaded, page 07's radiative-equilibrium profile as a dashed comparison, and the lapse-rate profile in a second panel with 9.8 and 6.5 K/km marked. Captured by `generate.py`.
2. **The lapse rate, measured**: print the mean over the convecting layer and compare it with 9.8 and with the 6.5 tranche 1 assumed.
3. **The knob**: double CO₂, integrate to the new equilibrium, print the warming and the number of steps. Print `budgets.summary` before and after, so convergence is read rather than asserted.
4. **Saving your own**: `states.save` the perturbed equilibrium to a local file and reload it — the reader has now done what the page's own data asset did.

**Prose:**

1. *Where this fits* — chapter 12, and the closing of the loop opened on page 01.
2. *This column is dry, and its numbers are not Earth's* — the paragraph the spec requires. 14-band radiation over a CO₂-only atmosphere; it equilibrates well below Earth's surface temperature; nothing here is comparable with tranche 1's numbers, and page 12 is where the missing gas goes back in.
3. *You are perturbing a shipped equilibrium* — say it, and show `states.describe` output as the evidence. Say why: spinning this up in the browser is thousands of steps, and watching that is page 07's job, not this one's.
4. *climt craft: ordering components defensibly* — radiation, then surface, then boundary layer, then adjustment. Give the reason for each adjacency and say which orderings are defensible alternatives. This is a physical argument, not a formatting one. Note where the wind relaxation sits and why it has to be constructed after the state rather than declared with the rest.
5. *climt craft: the budgets as a convergence test* — `budgets.toa_imbalance` and `budgets.surface_imbalance`, and what value counts as converged and why. Explicitly: eyeballing a flattening curve is not a convergence test.
6. *climt craft: saving and loading a state* — cell 4, and the provenance that has to travel with it.
7. *Knob* — CO₂, with the measured sensitivity and the step count it took.
8. *Physics exercises* — (i) why is the convecting layer's lapse rate 9.8 and not 6.5, given that Earth's is close to 6.5?; (ii) compare the upper column here with page 07's radiative equilibrium and explain why they are similar while the surfaces are not.
9. *Code exercises* — (i) halve CO₂ instead and check the response is roughly symmetric in log CO₂; (ii) remove `DryConvectiveAdjustment` from the stepper list, re-run, and recover page 07's profile.
10. *Going deeper* — page 12.

- [ ] **Step 4: Sidebar, chapter map, static fallback, render, commit**

```yaml
        - modelling-tour/11-dry-rce.qmd
```
```python
    "11-dry-rce.qmd": [(1, "11-dry-rce.png")],
```
```markdown
| [Dry RCE](11-dry-rce.qmd) | 12 | a 9.8 K/km adiabat, an emergent tropopause, a measured sensitivity |
```

**Browser cost check.** Cells 0–2 are instant (loading a state). Cell 3 is the CO₂ perturbation — `PAGE11_2XCO2_STEPS` = 1000 × ~0.18 s ≈ **3 minutes**. **Print the expected wall time in the prose above that cell.** If it exceeds ~10 minutes, apply Task 8 Step 8's levers in order.

```bash
git commit -m "docs(tour): page 11, where the lapse rate stops being an assumption

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

### Log — page 11's equilibrium and its sensitivity

**Run 2026-09-26, climt 0.31.0, Linux.** Every number from `scripts/experiments/tour_page11_measurements.py` (named per row) or the page's own cells; timings with `NUMBA_DISABLE_JIT=1`.

| | value | where |
|---|---|---|
| surface temperature | 266.48 K as shipped; 266.487 K 30-day mean (limit cycle 266.43–266.53) | cell 0; `convection` |
| TOA imbalance | −0.130 W m⁻² as shipped; +0.011 30-day mean | cell 0; `convection` |
| convecting layer (θ uniform) | 968–534 hPa in the shipped state, 11 levels, θ = 264.85 K | cell 2 |
| its lapse rate | **9.76 K/km** (g/c_p), every step | cell 2; `convection` |
| top of convection, step to step | 589–370 hPa on 56 of 60 steps, 836 hPa on 4; median ≈ 480 hPa | `convection` |
| tropopause | no temperature minimum: T falls to 111.73 K at the lid; lapse < 9 K/km above ~340 hPa, ~1.3 K/km at the top | cell 2 |
| lowest three layers | 4.9, 2.2, 8.1 K/km — stable; `SimpleBoundaryLayer` diffuses T, not θ | cell 2 |
| page 7's 14-band column at equilibrium | 10 000 steps at 12 h: T_s 267.78 K, TOA 0.000, top 101.2 K | `radiative` |
| RCE − page 7 | surface **−1.30 K**; air **+25.4 K** at 695 hPa, **9–10 K** at the top (cell 2 prints +10.5; the shipped top level still cools 1.17 K in 10 000 more steps) | cell 2; `settle` |
| OLR split, one LW call | air alone (surface at 1 K) 27.69 W m⁻²; surface shining through 212.27 = 74 % of σT_s⁴; surface +1 K, air fixed: +3.37 W m⁻² K⁻¹ (79 % of 4σT³) | `radiative` |
| 2×CO₂ instantaneous forcing | 4.08 W m⁻² | cell 3 |
| 2×CO₂ warming | **+1.15 K** (30-day mean perturbed − 30-day mean unperturbed) in **1000 steps**; within 0.1 K from day 164; 3.54 W m⁻² K⁻¹ | cell 3 |
| 0.5×CO₂ | −1.10 K | `co2` |
| cell 3 wall time | 1060 steps (60 unperturbed for the baseline): 36.6–39.7 s native → **≈ 2.2 min browser** (0.13 s/step) | page harness |
| page 7's column cost | 32.5 ms/step native → 10 000 steps ≈ 19 min browser | — |

**Three findings that differ from this task's text.**

1. *"A much cooler surface, a similar upper atmosphere"* (the reveal) does not hold. Against page 7's converged radiative equilibrium the surface is only 1.30 K cooler, and the air is 25 K warmer mid-troposphere and 10 K warmer at the top. The column is optically thin (74 % of σT_s⁴ reaches space), so warming the air adds only 4.43 W m⁻² of OLR, and a surface 1.30 K cooler at 3.37 W m⁻² K⁻¹ removes about the same. (A first draft zeroed the *air* to isolate the surface and got 85 %; that changes the temperature-dependent absorption and gives a surface term, 244.06, larger than the whole OLR. Chill the surface instead.) The page says so, and physics exercise 1 is built around it.
2. *The 2×CO₂ warming is +1.15 K, not +1.20 K.* The re-measurement ran `run_to_equilibrium`, which calls `integrate` in 50-step chunks and so restarts `AdamsBashforth` every 50 steps, and quoted single samples of a ±0.05 K limit cycle. The page runs one 1000-step `integrate` and compares 30-day means of the perturbed and unperturbed columns. The test bracket is ±40 % of 1.15 (0.69–1.61).
3. *The convecting layer cannot be found by calling the adjustment on the loaded state.* The file was saved straight after an adjustment, so the call moves nothing by more than ~1e-13 K. The test steps once without the adjustment, then calls it. Code exercise 2 (drop the adjustment, "recover page 7's profile") does not recover it: the boundary layer mixes up to 7 km in an unstable column, and the air stays 3.5–28.5 K warmer than page 7's after 1000 steps.

---

## Task 15: Page 12 — Moist RCE

*Course chapter 12 and beyond.*

Adds `EmanuelConvectionPython`. Loads the shipped moist equilibrium and perturbs it.

**The reveal.** Latent heating relaxes the lapse rate toward the moist adiabat — near the 6.5 K/km tranche 1 assumed, so the assumption is finally *explained* rather than replaced. Precipitation appears as a diagnostic, the moisture budget closes against the surface latent heat flux, and the climate sensitivity differs from page 11's dry one.

**The warning this page will print, and must explain rather than obey.** `EmanuelConvectionPython` is an `ImplicitTendencyComponent`, so putting it inside `AdamsBashforth` makes sympl emit, verbatim:

```
UserWarning: Using an ImplicitTendencyComponent in sympl TendencyStepper objects
may lead to scientifically invalid results. Make sure the component follows the
same numerical assumptions as the TendencyStepper used.
```

**Emanuel goes inside `AdamsBashforth` anyway** — that is how climt has always run it, and it works. The warning is generic to the component protocol, not a finding about this scheme at these timesteps, and restructuring the loop around it would buy nothing and cost the reader a special case to remember. What the page owes is the sentence explaining it, because the warning *will* appear in the reader's browser output on the first run and an unexplained warning in a teaching page reads as a bug.

There is a **second** warning from the same construction, which the page should also account for rather than let a reader puzzle over:

```
UserWarning: TimeSteppers should be given individual Prognostics rather than a
list, and will not accept lists in a later version.
```

That one is about `AdamsBashforth([a, b, c])` versus `AdamsBashforth(a, b, c)`. Decide once, in `_tour/stepping.py`, whether to pass a list or splat it — splatting silences it and is what sympl wants — and if the warning is silenced there, say nothing about it on the page. **Check this while writing the page**, because the two warnings appear together and the prose must match what the reader actually sees.

**Files:**
- Create: `docs/modelling-tour/12-moist-rce.qmd`
- Modify: `docs/_quarto.yml`, `docs/modelling-tour/index.qmd`, `docs/modelling-tour/_artifacts/generate.py`, possibly `docs/modelling-tour/_tour/stepping.py` (the splat)
- Test: `tests/test_modelling_tour.py`

**Interfaces:**
- Consumes: `_tour/states.py`, `_tour/stepping.py`, `_tour/budgets.py`, `_tour/assets.py`, `docs/modelling-tour/_data/rce_moist_equilibrium.npz`, `docs/modelling-tour/_data/rce_moist_2xco2_equilibrium.npz`.
- Produces: nothing. This is the last page.

- [ ] **Step 1: Decide the `AdamsBashforth` list-versus-splat question**

```bash
conda run -n climt python -W error::UserWarning -c "
import sympl, climt
from sympl import AdamsBashforth
sympl.set_backend(climt.UnytBackend())
lw = climt.CorkLongwaveRadiation(optics='correlated_k', table='earth_low_res_lw')
try:
    AdamsBashforth(lw, climt.SlabSurface())
    print('splat: no warning')
except UserWarning as w:
    print('splat warns:', w)"
```

If splatting is clean, change `_tour/stepping.py`'s two `AdamsBashforth(list(tendency_components))` calls to `AdamsBashforth(*tendency_components)`, re-run the whole tour suite, and commit that as part of this task. One warning left for page 12 to explain is better than two.

- [ ] **Step 2: Write the failing tests**

```python
MOIST_ADIABAT_K_PER_KM = 6.5     # nominal; the real one varies with T


@pytest.mark.slow
def test_page12_lapse_rate_lies_between_dry_and_moist_adiabatic(states):
    """Page 12's reveal: latent heating relaxes the lapse rate.

    Through the lower troposphere the moist column should lapse markedly less
    steeply than page 11's dry 9.8 K/km, and land near the 6.5 K/km every
    tranche 1 page assumed -- which is the point: the assumption is explained,
    not replaced.
    """
    generator = _generator()
    components = generator.moist_components()
    state, _ = states.load(str(DATA / "rce_moist_equilibrium.npz"), components,
                           grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))

    lapse = _lapse_rate_profile(state)
    p = state["air_pressure"].values[:, 0, 0]
    lower = (p[:-1] > 5.0e4)          # lower troposphere, below 500 hPa
    mean_lapse = float(np.mean(lapse[lower]))

    assert mean_lapse < DRY_ADIABAT_K_PER_KM - 1.0, (
        f"lower-tropospheric lapse rate {mean_lapse:.2f} K/km is not visibly "
        "less than dry adiabatic — latent heating should have relaxed it")
    assert 4.0 < mean_lapse < 8.5, (
        f"lapse rate {mean_lapse:.2f} K/km is outside the moist-adiabatic "
        "neighbourhood")


@pytest.mark.slow
def test_page12_moist_column_is_warmer_than_the_dry_one(states):
    """The water vapour greenhouse plus its feedback, as a difference.

    Tranche 1's page 6 drew this from prescribed profiles. Here two models
    that supply their own water (or do not) produce it.
    """
    generator = _generator()
    dry, _ = states.load(str(DATA / "rce_dry_equilibrium.npz"),
                         generator.dry_components(),
                         grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))
    moist, _ = states.load(str(DATA / "rce_moist_equilibrium.npz"),
                           generator.moist_components(),
                           grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))

    dry_surface = float(dry["surface_temperature"].values.ravel()[0])
    moist_surface = float(moist["surface_temperature"].values.ravel()[0])
    assert moist_surface > dry_surface + 5.0, (
        f"moist {moist_surface:.2f} K vs dry {dry_surface:.2f} K — adding a "
        "greenhouse gas the surface supplies should warm the column")


@pytest.mark.slow
def test_page12_precipitation_balances_evaporation_at_equilibrium(states):
    """The moisture budget closes against the surface latent heat flux.

    P ~= LHF / Lv. At equilibrium the column cannot be accumulating water.
    """
    generator = _generator()
    components = generator.moist_components()
    tendencies, steppers = generator.split(components)
    state, provenance = states.load(
        str(DATA / "rce_moist_equilibrium.npz"), components,
        grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))

    stepping_module = _load("stepping")
    tendencies = tendencies + [stepping_module.wind_relaxation(
        state, provenance["wind_m_s"], provenance["wind_timescale_hours"],
        initialise=False)]
    timestep = climt.UnytTimeDelta(hours=provenance["dt_hours"])
    stepping_module.integrate(tendencies, steppers, state, timestep, 50)

    budgets_module = _load("budgets")
    precipitation = budgets_module.precipitation_rate(state, timestep)
    evaporation = budgets_module.evaporation_rate(state)

    assert precipitation > 0.0, "an equilibrium moist column must precipitate"
    assert abs(precipitation - evaporation) < 0.25 * evaporation, (
        f"precipitation {precipitation:.3f} mm/day vs evaporation "
        f"{evaporation:.3f} mm/day — the moisture budget does not close")


@pytest.mark.slow
def test_page12_emanuel_inside_adamsbashforth_warns_and_still_works():
    """Pins the warning the page explains, and that the results are usable.

    If sympl ever stops emitting this, page 12's craft note describes
    something the reader will not see, and this is where we find out.
    """
    import warnings

    from sympl import AdamsBashforth

    convection = climt.EmanuelConvectionPython()
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table="earth_low_res_lw")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        AdamsBashforth(longwave, climt.SlabSurface(), convection)

    messages = " ".join(str(w.message) for w in caught)
    assert "ImplicitTendencyComponent" in messages, (
        "sympl no longer warns about an ImplicitTendencyComponent inside a "
        "TendencyStepper — page 12's craft note needs rewriting")
```

- [ ] **Step 3: Run them**

Run: `conda run -n climt python -m pytest tests/test_modelling_tour.py -k page12 -m slow -v`

If `test_page12_precipitation_balances_evaporation_at_equilibrium` fails with precipitation well below evaporation, the shipped moist state has not equilibrated in *moisture* even though its energy budget closed. That is a real finding: add a moisture-budget criterion to `scripts/generate_tour_equilibria.py`'s convergence test alongside the TOA one, regenerate, and record it in Task 9's log.

- [ ] **Step 4: Write the page**

`docs/modelling-tour/12-moist-rce.qmd`. Front matter's `resources:` as page 11's, with `_data/rce_moist_equilibrium.npz` **and `_data/rce_moist_2xco2_equilibrium.npz`**.

**Cells:**

0. Setup; the component list built visibly, Emanuel in the tendency list; load the state; rebuild page 08's wind relaxation on it with `initialise=False`, exactly as page 11's cell 0 does; print `states.describe(provenance)`, whose wind line is now part of what the reader checks. **This is the cell that emits the warning** — the prose immediately above it says so, and says what it means, so the reader meets the explanation before the warning.
1. **Headline figure**: the moist equilibrium profile against page 11's dry one, with the dry adiabat and a moist adiabat drawn over both, and the specific humidity profile in a second panel. Captured by `generate.py`.
2. **The lapse rate, measured**, in the lower troposphere, against 9.8 and 6.5.
3. **The moisture budget**: precipitation and evaporation printed side by side, with the residual.
4. **The knob, run offline**: load `rce_moist_2xco2_equilibrium.npz` beside the base state, print its `states.describe` block (660 ppm, 180 550 steps), and compare the two states: the surface warming, the two temperature profiles and the two humidity profiles. Compare the sensitivity with page 11's dry +1.20 K and say why they differ. *(Superseded 2026-09-27: page 11's number is +1.15 K.)* **The cell does not integrate.** At dt = 5 min the re-equilibration is 180 550 steps, about nine and a half hours in the browser, so it ships (see `_data/README.md`, and `generate_tour_equilibria.py --moist-2xco2`). The prose above the cell says so, and says that page 11 ran its own doubling live, so the reader knows what they are being handed and why. The short live run a reader *can* afford (say 500 steps from the doubled base, ~1.5 min) shows the warming starting, and is an optional code exercise, not the cell.

   **Which number to quote: +2.24 K**, about 1.9× page 11's dry +1.20 K. *(Superseded 2026-09-27: page 11 quotes +1.15 K, and page 12 compares against that.)* It is the difference between the two states *as they settle* when stepped on (`tour_rce_shipped_remeasure.py settle-moist`, Task 9 log 2026-09-27). The files agree with it: +2.262 K file to file, and +2.223 K between the gate's last-window means, which is what the provenance's `window_mean_warming_k` records. The moist gate now stops on flat 60-day trends, so the files sit within about 0.02 K of where the column settles. Quote the settled number and give the file difference as a cross-check. *(Until 2026-09-27 this said +1.82 K. That number came from states made over a fixed, 246 %-saturated surface; see the Task 9 log.)*

   **The page must not claim the moist column balances at the top.** `EmanuelConvectionPython` is not fully energy-conserving (known). This column settles with TOA ≈ −1.1 W m⁻² and stays there, and the file's own `window_mean_toa_w_m2` records it. `budgets.summary` will print that number, so the page says what it is before a reader finds it: the scheme's non-conservation, visible as a residual that does not decay. The stratosphere is in radiative balance and the slab is not storing anything, so this is not a slow adjustment. Contrast it with page 11's dry column, whose residual does decay. *(Superseded 2026-09-27, Task 15 log: the residual is `DryConvectiveAdjustment`'s moist-c_p bookkeeping, +0.97 of it, not Emanuel's, +0.03. Page 12 says so.)* The *moisture* budget (cell 3) is a separate claim, and still closes (P 2.64 vs E 2.64 mm day⁻¹ over 30 days).
5. **Timestep sensitivity, demonstrated rather than asserted**: run a short perturbation at two timesteps and compare. This is where the page cashes the cheque the warning note writes.

**Prose:**

1. *Where this fits* — chapter 12 and beyond; the arc closing.
2. *The warning you are about to see* — the paragraph the spec specifies, in one clear sentence: this is sympl saying it cannot verify that an implicitly-formulated component is safe inside a multi-stage stepper; here it is, climt has always run it this way, and cell 5 is where you check that claim yourself.
3. *climt craft: `ImplicitTendencyComponent`* — the third and last component protocol. `_tour/stepping.py` still has exactly two lists, tendency and stepper, and Emanuel goes in the first. Say why that is a decision and not an oversight.
4. *The lapse rate, and where 6.5 came from* — the tranche's closing move. Latent heating relaxes the dry adiabat toward the moist one; the number tranche 1 assumed on every page is now produced.
5. *climt craft: precipitation and convective diagnostics* — `convective_precipitation_rate` in mm/day from Emanuel, `precipitation_amount` in kg m⁻² per step from `GridScaleCondensation`, and why `budgets.precipitation_rate` needs the timestep to combine them.
6. *The moisture budget closes* — quote P and E and the residual.
7. *Knob* — CO₂ and surface temperature, with both sensitivities and the reason they differ from page 11's. Say that this experiment was run offline and why (step count × browser cost), and that the shipped file's provenance says exactly what was run.
8. *Physics exercises* — (i) why is the moist sensitivity different from the dry one, and which sign would you have predicted?; (ii) the moist adiabat is not a single number — at what level does the lapse rate here cross 6.5, and why there?
9. *Code exercises* — (i) halve the timestep and quantify the change in the equilibrium, which is cell 5 done properly; (ii) drop `GridScaleCondensation` and see what Emanuel alone does to the humidity profile.
10. *Going deeper, and where the tour ends* — say what is deliberately absent: no clouds, no sea ice, no shortwave component, no latitude. Name `CorkShortwaveRadiation` as deferred and why (it has not been stress-tested), so a reader who goes looking knows it exists and knows why it is not here.

- [ ] **Step 5: Sidebar, chapter map, static fallback, render, commit**

```yaml
        - modelling-tour/12-moist-rce.qmd
```
```python
    "12-moist-rce.qmd": [(1, "12-moist-rce.png")],
```
```markdown
| [Moist RCE](12-moist-rce.qmd) | 12 and beyond | the lapse rate tranche 1 assumed, finally explained |
```

**Browser cost check.** Cell 4 no longer integrates. It loads a second shipped state, so it is instant. That decision was taken because the live version would be 73 600 steps × ~0.19 s ≈ 3.9 h (since 2026-09-27, 180 550 steps ≈ 9.5 h). The levers in Task 8 Step 8 would not have saved it: they change the column, so the page would no longer be perturbing the equilibrium it loaded. Cell 5 (timestep sensitivity) is now the expensive one. Keep it to a short perturbation, and **print its expected wall time above the cell**.

```bash
git commit -m "docs(tour): page 12, and the 6.5 K/km that page 1 assumed

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

### Log — page 12's equilibrium, budget and sensitivity

*Fill in: surface temperature, lower-tropospheric lapse rate, precipitation and evaporation and their residual, the timestep-sensitivity result, browser wall time for cell 5.*

**Page 12 as written, 2026-09-27** (`scripts/experiments/tour_page12_measurements.py`; cells timed natively with `NUMBA_DISABLE_JIT=1`, browser ≈ 3.5×).

*Step 1 was already done:* `_tour/stepping.py` splats (`AdamsBashforth(*tendency_components)`), so only the `ImplicitTendencyComponent` warning appears; `test_page12_emanuel_inside_adamsbashforth_warns_and_still_works` pins both facts.

*The spec's reveal is false in this configuration.* The lapse rate does not relax to near 6.5 K/km:

| | value |
|---|---|
| surface temperature | 285.99 K (file); 285.996 K, 30-day mean; +19.51 K on page 11's dry 266.48 |
| adjusted (θ-uniform) layer | 968–792 hPa, 9.76 K/km, RH 32–83 % |
| 792–500 hPa lapse rate | 9.48 K/km (loaded state and 30-day mean); page 11's dry column 9.75 over the same levels |
| below 500 hPa, mean lapse | 9.07 K/km (30-day), 2×CO₂ 9.02; dominated by the boundary-layer layers, so the page quotes the 792–500 band |
| moist adiabat on the column's own T | 5.57 K/km at the lowest level, 6.5 at 876 hPa (270.6 K), 9.2 at 534 hPa |
| saturated parcel, mean 1000–500 hPa | 6.52 K/km from 284 K; 6.19 from 286 K (6.24 from this column's surface) |
| P / E, 30 days | 2.645 / 2.640 mm/day, quoted as 2.65 / 2.64 (the table below rounds P to 2.64); Emanuel supplies 0.030 of it |
| P / E, cell 3 (2 days) | 2.62 / 2.59, residual +0.03; 2-day means of P − E scatter by 0.17 (sd) |
| SH / LH, Bowen | 27.96 / 76.40 W/m², 0.37 (30 days); cell 3: 27.8 / 74.9 |
| grid-scale condensation heating | 745–370 hPa, up to 2.5 K/day |

*Correction, 2026-09-27 (review of 36c8062).* This log and the first version of the page said Emanuel's "CAPE is never positive (mean −4 J/kg)". Wrong: `atmosphere_convective_available_potential_energy` from `EmanuelConvectionPython` (`pure_python_v3`) is the scheme's internal CAPEM + BYP, a cloud-top sum in K·Δp/p, not CAPE. Traced through `_convect_functional_np` (one day at 5 min, JIT off): Emanuel's own lowest-level parcel is up to 0.89 K *cooler* (virtual) than the air between 968 and 911 hPa (last call of the day), then buoyant from 836 hPa, +2.4 K at 745 and +6.1 K at 370, up to about 180 hPa. Its positive area is about 1700 J/kg (mean over the day). What limits the scheme is its closure, which looks only at cloud base: DTMA = TVPPLCL − TVAPLCL + DTMAX(0.9) + DTPBL averages −0.27 K (range −0.63 to +0.07), DTPBL alone −0.70, because the parcel crosses the stable 1002–968 hPa layer and arrives in the dry-adjusted layer ≈ 0.9 K cooler than it, right at DTMAX. (`tour_page12_measurements.py closure`.) So the mass flux, updated as (1 − 0.1·dt/300 s)·CBMF + 0.01·DTMA each call, stays near zero: mean 1.1e-3, max 4.8e-3 kg m⁻² s⁻¹. (The "> 0 on 63 % of steps" in the first version came from the post-step state over 5 days; the closure trace over 1 day gives 33 % of calls. Neither is quoted now.) At dt = 10 min (`closure(steps=144, minutes=10.0)`) DTMA averages −0.44 and the mean mass flux 1e-4. The page carried the same claim; the follow-up commit fixes it.

*Why only the CO₂ knob.* The spec lists "CO₂, and surface temperature". Here the surface temperature is prognostic (the slab under 240 W m⁻²), not a parameter, so there is nothing honest to turn; the page's knob is CO₂ alone.

Why: the grid-scale condensation does ~99 % of the raining, at the levels that saturate above the dry-adjusted layer. `nogsc` (condensation removed): within ~10 days Emanuel carries all the rain (2.56 mm/day) and 792–589 hPa lapses at 7.2–8.5 K/km, within 0.8 K/km of its own moist adiabat, with 157 % RH at 792 hPa. `nodca` (adjustment removed): Emanuel carries the rain but the column cools steadily (278.9 K after 300 d, not converged).

*The TOA residual is not Emanuel's.* Metering column enthalpy (c_p,dry T + L_v q) per component over 5 days: `DryConvectiveAdjustment` +0.967 W/m², Emanuel +0.026, everything else 0 (boundary layer = SH + LH exactly); TOA −1.104. The adjustment conserves enthalpy with a moist heat capacity (its moist-c_p change is 0.0000), while every other component counts dry c_p. The page says so; `_data/README.md`, the generator comment and the residual-test comment were corrected.

*Sensitivity.* File difference +2.26 K; 30-day mean minus 30-day mean +2.233 K; settled +2.24 K (quoted); page 11's dry +1.15 K. Radiation-only decomposition (file states): forcing +4.69 W/m²; temperatures +8.93 W/m² (3.95 per K); vapour −4.27 (−1.89 per K); with the vapour fixed the warming would be 1.19 K. So the water-vapour feedback is the whole difference.

*Timestep* (60 days from the shipped state, last-30-day means, dt = 2.5 / 5 / 10 / 20 min): Emanuel's rain 0.056 / 0.029 / 0.007 / 0.000 mm/day; Ts 285.978 / 286.007 / 286.043 / 286.065 K (286.0649, so 286.06 on the page); lapse below 500 hPa 9.07 / 9.07 / 9.07 / 9.08; TOA −1.05 / −1.11 / −1.16 / −1.96. Cell 5 (2 days, 10 min vs cell 3's 5 min): Emanuel 0.002 vs 0.029, P 2.75 vs 2.62, Ts 286.02 vs 285.99, largest 2-day-mean ΔT −0.26 K at 268 hPa. Emanuel's cloud-base mass flux does *not* scale as 1/dt (1.34e-3 / 1.22e-3 / 1.02e-3 / 3.6e-4), so the page does not claim a mechanism.

*Cell cost* (43.7 ms/step native, JIT off; ~0.15 s/step browser): cell 0 1.5 s native (~5 s browser); cell 3, 576 steps, 25.3 s (~1.5 min); cell 4, no stepping, 1.4 s (~5 s); cell 5, 288 steps, 12.3 s (~45 s). Code exercises: 1152 steps ~3 min; 2880 steps ~7 min.

**Current, 2026-09-27** (saturated surface, moist trend gate; Task 9 log). Generated with `generate_tour_equilibria.py --moist --moist-2xco2`; the 2×CO₂ state is stamped with the base's `saved_at` 2026-09-27T02:22:56:

| | base | 2×CO₂ |
|---|---|---|
| CO₂ | 330 ppm | 660 ppm |
| steps (dt = 5 min) | 295 800 (≈ 1027 d) | 180 550 after doubling (≈ 627 d) |
| surface temperature (file / last-window mean) | 285.989 / 286.014 K | 288.251 / 288.237 K |
| TOA imbalance (file / last-window mean) | −1.270 / −1.095 W m⁻² | −1.197 / −1.029 W m⁻² |
| surface RH | 100 % | 100 % |
| SH / LH, Bowen ratio | 27.9 / 76.3 W m⁻², 0.37 | 25.3 / 83.2 W m⁻², 0.30 |
| P / E | 2.64 / 2.64 mm day⁻¹ | 2.88 / 2.88 mm day⁻¹ |
| file difference | | +2.262 K (window means +2.223 K) |
| settled, 60 000 steps on (mean of last 30 000) | 285.990 K, TOA −1.045 | 288.228 K, TOA −1.066 |
| **settled 2×CO₂ warming** | | **+2.24 K** |

*Superseded: the states made over a fixed 0.015 kg kg⁻¹ surface (246 % / 220 % RH).* Already measured, 2026-09-26 (`generate_tour_equilibria.py --moist-2xco2`, from the shipped base saved 2026-09-04T04:59:51):

| | base | 2×CO₂ |
|---|---|---|
| CO₂ | 330 ppm | 660 ppm |
| steps (dt = 5 min) | 200 550 | 73 600 after doubling (≈ 256 d) |
| surface temperature | 279.965 K | 281.613 K |
| TOA imbalance | −0.485 W m⁻² | +0.500 W m⁻² |
| file difference | | +1.647 K (not the number to quote) |
| settled, 60 000 steps on (mean of last 30 000) | 279.862 K, TOA +0.306 | 281.683 K, TOA +0.347 |
| **settled 2×CO₂ warming** | | **+1.82 K** |

---

## Task 16: Close out — the index, the full render, and the link check

Six pages exist. The section still describes itself as a tranche in which nothing integrates in time.

**Files:**
- Modify: `docs/modelling-tour/index.qmd`
- Modify: `docs/modelling-tour/_data/README.md` (if anything shifted)
- Modify: `HISTORY.rst`
- Test: the whole suite, and a real browser

**Interfaces:**
- Consumes: everything.
- Produces: a shippable branch.

- [ ] **Step 1: Rewrite the index's "what this does not cover" section**

`docs/modelling-tour/index.qmd` currently says:

```markdown
## What this tranche does not cover

Nothing here integrates in time. There are no time-stepping loops, no convection, and
no surface fluxes: each page makes a single radiation call, or a sweep of them, on a
profile you prescribe. [...]

Radiative-convective equilibrium, moist convection and surface energy balance build on
chapters 10–12 of the notes, and belong to a later tranche.
```

Replace with a two-part structure: what the first six pages do (single calls on prescribed profiles — and *why*, which is still a good reason worth keeping), what the second six add (time integration, a surface, turbulence, convection), and then the genuine gaps:

```markdown
## Two halves, and why

**Pages 1–6 never integrate in time.** Each makes a single radiation call, or a sweep
of them, on a profile you prescribe. That is deliberate: it keeps the causal chain
short enough to see, so that when a number changes you know precisely which line
changed it.

**Pages 7–12 do nothing else.** They step a column forward, give it a surface, a
boundary layer and convection, and arrive at radiative-convective equilibrium. Page 7
is where the two halves meet: the profile page 4 wrote down by hand is the profile
page 7 produces from nothing.

The arc closes on a number. Every page in the first half prescribed Γ = 6.5 K/km.
Page 11 shows that radiation plus dry convection gives 9.8. Page 12 shows that latent
heating relaxes it back toward 6.5 — which is where the assumption came from.

## What the tour still does not cover

- **Sunlight in the atmosphere.** Every page prescribes the absorbed shortwave flux at
  the surface. There is no shortwave component anywhere, so surface albedo and zenith
  angle are not knobs, and the surface energy balance is one-sided in the shortwave.
  `CorkShortwaveRadiation` exists and is deliberately deferred until it has been
  exercised properly.
- **Clouds, sea ice, and land.** No cloud scheme, no `IceSheet`, no `SecondBEST` or
  `BucketHydrology`.
- **Anything with a horizontal dimension.** Every page is one column. There is no
  latitude, no dynamics, and no transport — which is why, on page 8, wind speed turns
  out not to be a knob.
```

Then extend the chapter-map table with the six rows added in Tasks 10–15, if any were missed.

- [ ] **Step 2: Regenerate every static fallback figure**

```bash
conda run -n climt python docs/modelling-tour/_artifacts/generate.py
```

Expected: all thirteen PNGs (tranche 1's seven plus tranche 2's six) rewritten with no traceback. This execs every page's cells in order, so it is the single cheapest check that all twelve pages still run.

`git diff --stat docs/modelling-tour/_artifacts/` — tranche 1's seven should be unchanged in content. If any moved, find out why before committing: the only tranche 1 code this branch touched is `tables.py` (Task 3), and that refactor was contract-preserving.

- [ ] **Step 3: Run the full suite, both markers**

```bash
conda run -n climt python -m pytest tests/ -m "not slow" -q
conda run -n climt python -m pytest tests/test_modelling_tour.py -m slow -q
conda run -n climt python -m pytest tests/test_simple_boundary_layer.py tests/test_cork_lw.py \
                                    tests/test_conservation.py -q
```

Expected: green. The `slow` set is the long one — allow half an hour.

- [ ] **Step 4: Verify the browser path for real**

The tests prove the physics; they cannot prove the pages run under Pyodide, and this branch changed the two things most likely to break there (Task 1's numba-free path, and Task 3's cross-module import).

```bash
CLIMT_PURE_PYTHON=1 conda run -n climt python -m pip wheel . --no-deps -w /tmp/climt_wh
conda run -n climt python scripts/serve_wheel.py /tmp/climt_wh 8912    # leave running
cd docs && quarto preview
```

Point each of pages 07–12's `pyodide: packages:` at the local wheel URL, then open all six. For each, confirm:

1. the boot cell prints `climt <version> ready — has_fortran: False`;
2. every cell runs to a figure with no traceback in the browser console;
3. page 12 emits the `ImplicitTendencyComponent` warning **and** the page's prose above that cell explains it;
4. no cell exceeds the wall time its own prose promises.

**Revert the front-matter wheel pins to `climt==0.31.0` before committing.** A page pinned at `http://localhost:8912` is broken for every reader.

- [ ] **Step 5: Check every link**

```bash
cd docs && quarto render && cd ..
grep -rn "](\.\./\|](0[1-9]\|](1[0-2]" docs/modelling-tour/*.qmd | grep -o "](\([^)]*\))" | sort -u
```

Walk the list. In particular: `04-gray-equilibrium-tested.qmd`'s re-pointed link now resolves to `07-stepping-a-column.qmd`, and nothing anywhere resolves to `radiative-transfer/09-live-rce.qmd`.

```bash
grep -rn "09-live-rce" docs/_site/ | head
```

Expected: no output.

- [ ] **Step 6: Re-check the experiment-artifact gate**

Task 2 edited `cork/` and regenerated the artifacts. Confirm nothing has drifted since:

```bash
conda run -n climt python scripts/build_experiments.py --check
```

Expected: clean.

- [ ] **Step 7: HISTORY.rst**

Add an entry naming, in this order: the `SimpleBoundaryLayer` no-numba fix (a **user-visible bug fix**, not a docs change — anyone running climt under Pyodide or with the JIT disabled hit it), the CORK non-finite guard, the six new pages, and the removal of `radiative-transfer/09-live-rce.qmd`.

- [ ] **Step 8: Request a whole-branch review**

Use `superpowers:requesting-code-review` against `develop`. Point the reviewer at the two things this branch does that a diff will not make obvious:

- **page 08's second subject** — a column with no dynamics spins its own wind down, so every page from 08 on carries a Newtonian relaxation as a momentum source, and both shipped equilibria were spun up with one. Found by measurement; it changes what page 08 teaches and what the shipped states contain;
- **the two library fixes**, one of which (`SimpleBoundaryLayer`) is a real bug in shipped code that every existing test passed over.

- [ ] **Step 9: Commit and open the PR**

```bash
git add -A
git commit -m "docs(tour): close out the RCE tranche — index, fallbacks, history

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
git push -u origin feature/modelling-tour-rce
gh pr create --base develop --title "Modelling Tour tranche 2: radiative-convective equilibrium" --body "$(cat <<'BODY'
Six new pages, 07-12, in `docs/modelling-tour/`, continuing PR #226. Where tranche 1
made single radiation calls on prescribed profiles, this tranche integrates in time
and arrives at radiative-convective equilibrium.

It also absorbs and deletes `docs/radiative-transfer/09-live-rce.qmd`, whose demo
page 07 now teaches properly, and consolidates the two Pyodide boot includes into one.

Two library fixes ride along, both found while measuring for this work:

- `SimpleBoundaryLayer` forwarded a unit-carrying timestep into its kernel. numba
  stripped the units, so every test passed; Pyodide has no numba, and the component
  raised `UnitOperationError` in the browser. This is a user-visible bug.
- `CorkLongwaveRadiation` returned NaN with only a numpy `RuntimeWarning` when the
  integration went unstable, which in the browser is a blank figure and no error. It
  now raises, naming the timestep.

Spec: `docs/superpowers/specs/2026-08-20-modelling-tour-rce-design.md`
Plan: `docs/superpowers/plans/2026-08-20-modelling-tour-rce.md`

🤖 Generated with [Claude Code](https://claude.com/claude-code)
BODY
)"
```

---

## Self-review

Run against the spec, section by section.

**Spec coverage.**

| Spec section | Task |
|---|---|
| Six pages 07–12 | 10, 11, 12, 13, 14, 15 |
| `_tour/assets.py`, `tables.py` becomes a caller | 3 |
| `_tour/stepping.py` | 4 |
| `_tour/states.py` | 7 |
| `_tour/budgets.py` | 6 |
| `_data/rce_dry_equilibrium.npz`, `rce_moist_equilibrium.npz` (and `rce_moist_2xco2_equilibrium.npz`, added after §7) | 9 |
| `scripts/generate_tour_equilibria.py` | 9 |
| `serve_wheel.py` moved | 5 |
| Delete `09-live-rce.qmd` and its five inbound references | 5 |
| Delete `climt-live-setup.qmd`, `rce_helpers.py` | 5 |
| Fold `test_live_rce_demo.py` into `test_modelling_tour.py` | 5 |
| Corrections inherited not copied (`dt`, `nz`, `UnytBackend`) | 4, 8, 10 |
| Library change: NaN guard | 2 |
| Task 0 measurement campaign, all five items | 8 |
| Residual test, not dependency hashing | 9 |
| Static fallback figures | 10–15, 16 |
| Tests parametrised over the table the page declares | Global Constraints; enforced per page task |
| `CorkShortwaveRadiation` deferred | Global Constraints; stated on page 12 and the index |

**Three additions the spec does not contain, all from measurement:**

1. **Task 1** — `SimpleBoundaryLayer` is broken without numba, which means broken in the browser. Four of the six pages use it. This is a blocker, not an improvement, and it is sequenced first.
2. **Wind relaxation (Task 4 Step 8), and page 08's second subject.** A column with no dynamics spins its own wind down: three initialisations spanning 0–10 m s⁻¹ all end at a lowest-level wind of 0.000 m s⁻¹ and the same equilibrium within 0.05 K. So the spec's "wind speed" knob does not work *as specified* — but rather than dropping it, every configuration from page 08 on gets a momentum source, and page 08 teaches it: climt's own assign-it-back idiom (`examples/column_code_with_slab.py:120`), then `sympl.RelaxationTendencyComponent` as the cleaner form, then the difference between the two. With the wind held up, wind speed *is* a knob — sensible heat flux 14.9 → 49.9 W m⁻², jump 10.04 → 6.07 K over 0–10 m s⁻¹.
3. **`sympl.RelaxationTendencyComponent` is unusable under `UnytBackend` as shipped** — its tendency units come out of pint as `'1.0 meter / second ** 2'`, which `unyt` refuses to parse. A two-line subclass in `_tour/stepping.py` fixes it, with a test pinned to the reason so the shim can be removed if sympl changes.

**A consequence worth stating, because it propagates:** the wind forcing is not confined to page 08. Pages 09, 11 and 12 carry it, and so do both shipped equilibrium states and `scripts/generate_tour_equilibria.py`. Page 08 teaching the fix and later pages silently reverting to a dead-calm column would be incoherent, and their surface fluxes would be wrong by a factor of three.

**One spec assumption flagged for checking rather than inherited:** page 09's "latent flux exceeds sensible over a saturated ~288 K surface". This tranche's columns do not sit at 288 K. Task 12 measures the converged temperature first and writes the claim there.

**Numbers verified before this plan was written**, so the pages have real targets rather than hopes: page 07's gray column converges to 0.103 K of page 04's analytic profile and 0.63 K of the skin temperature; page 08's surface–air discontinuity goes 15.15 K (radiative) → 10.04 (dead calm) → 7.50 (5 m s⁻¹) → 6.07 (10 m s⁻¹), and 11.37 → 6.93 K across the `z0` range; relaxation to 5 m s⁻¹ leaves the lowest level at 2.88 m s⁻¹ against a hard reset's 5.00, and moves 34.9 W m⁻² against 42.9; page 10's adjustment gives exactly uniform θ and conserves enthalpy to 1.8 × 10⁻¹⁶; the per-step cost of all five stacks, native and JIT-free.

**Type and name consistency.** `stepping.wind_relaxation(state, speed, timescale_hours=24.0, quantity_name="eastward_wind", initialise=True)` is called with `initialise=False` on every path that loads a shipped state (Tasks 9, 14, 15) and with the default when building fresh (Tasks 8, 9's generator, 11, 12) — flattening a loaded, drag-sheared wind profile to a uniform value would discard part of the equilibrium. `stepping.integrate(tendency, stepper, state, timestep, n_steps)` has the same signature in Tasks 4, 5, 6, 8, 9, 10–15. `budgets.toa_imbalance(state)` takes no `solar` argument anywhere. `states.load(name, components, grid_state=...)` returns `(state, provenance)` in Tasks 7, 9, 14, 15. `assets.locate(name, base_url)` is used by `tables.py` and `states.py` only. `generate_tour_equilibria.split(components)` is the single definition of which component goes in which list, used by the generator, the residual test and pages 11 and 12.

**Placeholders.** Three constants are deliberately marked as placeholders and each has a step that replaces it from a measurement: `MOIST_DT_HOURS` (Task 9 Step 4), `PAGE9_STEPS` (Task 12 Step 2), `PAGE11_2XCO2_STEPS` (Task 14 Step 1). Three test thresholds are loose pending Task 8 and each has a step that tightens it from the measurement. No step says "add error handling", "write tests for the above", or "similar to Task N".

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-08-20-modelling-tour-rce.md`.

**Sequencing note that is not optional:** Task 1 before anything that runs `SimpleBoundaryLayer`, Task 4 before Task 5 (which deletes the code Task 4 lifts), Task 8 before Task 9, and Task 8 before any `.qmd` is written. Within Tasks 10–15 the pages are independent and could be parallelised, but each depends on Task 8's log being filled in.
