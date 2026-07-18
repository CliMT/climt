# In-Browser Non-Grey Radiative Equilibrium Demo

**Date**: 2026-07-19
**Status**: Draft — design approved, pending spec review

## Context

The climt website (Quarto, `docs/_quarto.yml`) has a "Radiative Transfer"
theory track — eight chapters (`docs/radiative-transfer/01-why-nongrey.qmd`
… `08-multiplanet.qmd`) that argue, with static figures, why non-grey
radiation matters. It is good theory, but it is *only* theory: the reader
cannot run anything.

Two things now make it possible to make that argument runnable in the browser:

1. **CORK** (correlated-k non-grey radiation, `CorkLongwaveRadiation` /
   `CorkShortwaveRadiation`) is 100% pure Python with a numba-*optional*
   kernel layer.
2. The **pure-Python Emanuel** port is available at the top level.

The pure-Python *science* for a real non-grey radiative equilibrium already
exists and runs natively. What is missing is (a) the packaging that lets climt
install and run under Pyodide, and (b) the front-end that puts editable,
runnable climt code inside the theory chapters.

This deliverable is **self-contained**: it finishes the minimal pure-Python
packaging *and* builds the live demo on top. It builds directly on the
packaging analysis in
[`2026-03-29-pyodide-pure-python-support-design.md`](./2026-03-29-pyodide-pure-python-support-design.md)
(Revision 2026-07-19 — CORK-aware Phase 2); that document's task list is the
authoritative source for the packaging items, which are reproduced here in the
scope this demo requires.

### Wider goal: a reusable template

This is not a one-off. The intent is to eventually host *lots* of runnable
climt code — e.g. the course notebooks in `~/github/model_tour_climate`, and
to let third parties embed climt in **their own** websites. So the wheel is
hosted as a **GitHub release asset** (a stable, versioned URL anyone can
`micropip.install`), and the Pyodide boot/install mechanism is written as a
**documented, copyable template**, not glue buried in one page.

## Decisions (from brainstorming)

| Question | Decision |
|----------|----------|
| Reader experience | **Editable Pyodide code cells** (quarto-live), inline in the chapter prose. Not a Panel dashboard, not a JupyterLite iframe. |
| Packaging scope | **In-scope.** One deliverable: finish the minimal pure-Python wheel, host it, then build the demo. |
| Placement | **Retrofit key chapters** (Ch.1, Ch.6) **plus one flagship RCE page.** |
| Flagship physics | **Pure radiative equilibrium** (radiation + `SlabSurface`, no convection). Cleanest isolation of the radiation difference; fastest to run. |
| Gray vs non-grey mechanism | **Same `CorkLongwaveRadiation` component for both**, swapping only `table=`: a single-band constant-k table (gray) vs the multi-band Earth table (non-grey). **No `GrayLongwaveRadiation` in the flagship.** The reader sees that grey-vs-non-grey lives entirely in the absorption *data*, not the code. This is already proven bit-exact by `tests/test_grey_limit.py`. |
| Ch.8 multiplanet | **Deferred.** Out of scope for this deliverable. |
| Wheel hosting | **GitHub release asset** (versioned URL, reusable as a template by third-party sites). |
| IceSheet TDMA | **In-scope.** Complete full scipy removal now rather than leave a loose end. |

## Architecture & delivery stack

### Runtime: quarto-live

Use the [`quarto-live`](https://github.com/r-wasm/quarto-live) extension. It
embeds a Pyodide runtime in the rendered HTML and renders editable code cells
with live output (including matplotlib). The reader edits the *actual* climt
code and re-runs in place; no separate notebook environment.

**Validate this extension first** (see Task Order step 0): confirm, on a
throwaway page, that quarto-live can (a) `micropip.install` a wheel from an
absolute URL, (b) import climt, (c) render a matplotlib figure. If quarto-live
proves immature for any of these, the fallback is a hand-rolled
Pyodide + CodeMirror include (`docs/_includes/`), same wheel and boot contract.

### The wheel: GitHub release asset

Build `climt-<version>-py3-none-any.whl` with `CLIMT_PURE_PYTHON=1` and attach
it to a **GitHub release**. Pages install it by pinned URL:

```python
import micropip
await micropip.install(
    "https://github.com/CliMT/climt/releases/download/<tag>/climt-<version>-py3-none-any.whl"
)
```

A pinned, versioned URL means (a) the demo is reproducible, (b) any external
site can copy the same one-liner. The tag/version is defined **once** in a
shared include so bumping it updates every page.

### Data: gray + non-grey CORK tables as `.npz`, no scipy

Both flagship cases run through `CorkLongwaveRadiation`; only the k-table
differs. **Both** tables must therefore load in-browser, so the `.npz`/numpy
reader is on the critical path (not optional). The k-tables ship **inside** the
wheel as `.npz`, read by a small numpy loader — no `scipy.io.netcdf_file`.

Tables shipped for the demo (all already exist as `.nc`, to be converted to
`.npz`):

- `single_band_unit_lw` — the **gray** table (1 band, 1 g-point, constant k;
  956 bytes). Built by `scripts/generate_single_band_unit_table.py`.
- `earth_low_res_lw` (3 MB `.nc`) + `earth_low_res_sw` (17 KB) — the
  **non-grey** Earth correlated-k tables. `.npz` (compressed) should be
  comparable or smaller; this dominates the page-weight budget.

Earth-only keeps the download bounded. (Mars/Venus tables are not needed here;
Ch.8 is deferred.)

### Boot boilerplate hidden, science visible

Each live page begins with **one collapsed setup cell** (`code-fold: true`)
that boots Pyodide, installs the wheel, and defines a tiny helper so the
*visible* cells stay short and show real climt calls:

```python
def integrate_to_equilibrium(components, state, timestep, n_steps):
    """Step a list of components forward to (near) radiative equilibrium."""
    ...
```

This setup cell **is the reusable template** — documented in an appendix page
so others can paste it into their own Quarto site (see "Reusable template").

## Content plan

### Flagship: `docs/radiative-transfer/09-live-rce.qmd` (new)

Title: "Radiative equilibrium: gray vs non-grey — live." Added to the section
sidebar in `_quarto.yml` after Ch.8 (or repositioned as the section's live
capstone).

- **Setup cell** (folded): boot + install + `integrate_to_equilibrium`.
- **Cell A — gray:** `CorkLongwaveRadiation(optics="correlated_k",
  table="single_band_unit_lw")` + `SlabSurface`; build a single column via
  `get_default_state`, integrate to RE, plot T(p). (Shortwave held identical
  to Cell B so only the LW spectral treatment differs.)
- **Cell B — non-grey:** the **same** `CorkLongwaveRadiation` call with
  `table="earth_low_res_lw"` (+ `CorkShortwaveRadiation` /`Instellation`) +
  `SlabSurface`; integrate to RE; overlay T(p) on Cell A. The diff between
  Cell A and Cell B is *one string* — the pedagogical punchline.
- **Cell C — the payoff:** OLR / per-band flux diagnostics showing the
  **stratosphere and stratospheric cooling** the non-grey table produces and
  the gray table cannot; call out gray's convectively-unstable troposphere as
  a teachable deficiency (motivates the deferred RCE follow-up).
- **Cell D — CO₂ knob:** editable ppm; re-run to watch the forcing move the
  effective emission level. All cells are editable.

Default step counts are tuned for an acceptable in-browser wait (see
Performance); a comment shows how to increase them for a tighter equilibrium.

### Ch.1 retrofit: `docs/radiative-transfer/01-why-nongrey.qmd`

- Make the "mean of an exponential is not the exponential of the mean" toy
  (Fig 1.1) **live** — pure numpy, trivial under Pyodide.
- Add a minimal live RE contrast (a stripped-down teaser of the flagship).
- Replace the "Try it yourself → open `examples/…ipynb`" callout with an
  in-page runnable cell.

### Ch.6 retrofit: `docs/radiative-transfer/06-picket-fence.qmd`

- Live `CorkLongwaveRadiation` per-band flux diagnostic, replacing the callout
  that currently points at a static notebook.

### Deferred

- Ch.8 multiplanet live "switch planet" (needs Mars/Venus `.npz` tables).
- Panel/dashboard UX.
- Radiative-*convective* equilibrium (convection on). Explicitly seeded as the
  natural next chapter by Cell C's "unstable troposphere" callout.

## Packaging work (in-scope prerequisite)

Minimal set required for a browser-runnable pure wheel, drawn from the Phase 2
revision. File paths per that document.

1. **`find_packages()` + `package_data`.** `setup.py` currently lists
   `packages=[...]` explicitly, so a built wheel omits subpackages
   (`climt._components`, `climt._core`, `climt._data`, …). Switch to
   `find_packages()`. Verify with a **real `pip wheel` + fresh-venv install**,
   not the editable checkout that masks this today. (`package_data` CORK path
   already fixed.)
2. **Silent Fortran degradation** + `has_fortran_extensions()`. Replace
   import-time warnings in the five Fortran component files with silent
   availability flags; raise a clear `ImportError` only at instantiation. Add
   `has_fortran_extensions()` to `climt/__init__.py`.
3. **scipy removal (complete):**
   - `climt/_core/initialization.py`: `scipy.interpolate.CubicSpline` →
     `np.interp` (top-level import, so required for a clean scipy-free
     `import climt`).
   - CORK k-table reader (`cork/optics/correlated_k.py`): the `scipy.io`
     import is already function-local (`_load_netcdf_table`), so it does not
     block `import climt`. But **both** flagship tables load at demo runtime,
     so this is critical-path: convert the shipped `.nc` tables the demo uses
     (`single_band_unit_lw`, `earth_low_res_lw`, `earth_low_res_sw`) to `.npz`,
     make `load_k_table` resolve `.npz` without scipy, and give a clear error
     if a `.nc` table is requested without scipy. (Conversion can reuse
     `scripts/generate_single_band_unit_table.py` for the gray table.)
   - IceSheet: new `climt/_core/tridiagonal.py` (Thomas/TDMA solver,
     numba-optional) replacing `scipy.sparse.spdiags + spsolve` in
     `climt/_components/surface_ice.py`. **Validate TDMA against scipy before
     removal.**
   - Drop `scipy` from `install_requires` in `setup.py`.
4. **Re-export CORK** in `climt/__init__.py` `__all__`
   (`CorkLongwaveRadiation`, `CorkShortwaveRadiation`) — currently only
   reachable via `climt._components`.
5. **`CLIMT_PURE_PYTHON=1` build path** → `py3-none-any` wheel; add the CI job
   in `.github/workflows/release_climt.yml` and attach the wheel to the GitHub
   release.
6. **Host + boot:** publish the wheel as a release asset; the shared setup cell
   installs it by pinned URL.

## Reusable template

Add a short appendix page (e.g. `docs/radiative-transfer/live-code-template.qmd`
or a `docs/get-started/` how-to) that documents the copy-paste recipe for
embedding runnable climt in *any* Quarto site: the quarto-live dependency, the
`micropip.install(<release-url>)` snippet, and the `integrate_to_equilibrium`
helper. This is what makes the deliverable a template for `model_tour_climate`
and third-party sites, not a one-off.

## Performance / UX

From the Phase 2 benchmark: a single-column RCE step is ~7 ms native; apply a
~3–10× WASM factor and a full equilibrium (~thousands of steps) is tens of
seconds to a few minutes in-browser. Mitigations:

- Small level count and a **loose default step budget** for the flagship, with
  a commented knob to tighten it.
- The folded setup cell installs the wheel **while the reader reads the prose**,
  so the wait overlaps with reading.
- A visible progress indication during the integration loop.
- The gray single-band table is far cheaper (1 band × 1 g-point) than the
  multi-band Earth table, so Cell A is the fast path and Cell B the payoff.

## Verification

1. **Real wheel, fresh venv:** `pip wheel` the pure wheel, install in a clean
   venv (no editable checkout, no scipy), and run an RE loop headless — the
   gray-table and non-grey-table CORK columns both reach a stable profile, and
   the non-grey column develops a stratosphere the gray one lacks.
2. **Clean import:** `import climt` triggers no scipy import and no Fortran
   warnings; `has_fortran_extensions()` returns False on the pure wheel;
   Fortran components raise a clear `ImportError` on instantiation.
3. **TDMA correctness:** `solve_tridiagonal` matches scipy on the IceSheet test
   cases before scipy is removed.
4. **Quarto build:** `quarto render` succeeds for the new/retrofitted pages.
5. **Real-browser end-to-end** (Playwright or manual): the flagship page boots
   Pyodide, installs the wheel from the release URL, runs the gray + CORK RE
   cells, and renders both T(p) plots. This is the definition of done.
6. **Page-weight budget** recorded: Pyodide + numpy + wheel + Earth tables.

## Risks

- **k-table download size** in the browser — mitigate with Earth-only `.npz`.
- **quarto-live Python maturity** (micropip-from-URL, matplotlib) — validate at
  step 0 before committing; hand-rolled Pyodide include is the fallback.
- **In-browser runtime** for a full equilibrium — mitigate with a loose default
  step budget and overlapped install.

## Task order

0. **Spike:** validate quarto-live (micropip-from-URL + climt import +
   matplotlib) on a throwaway page. Decide extension vs. hand-rolled fallback.
1. Packaging: `find_packages()` + verify via real wheel/fresh-venv install.
2. Silent Fortran degradation + `has_fortran_extensions()`.
3. scipy removal: `np.interp` → CORK `.npz`/numpy reader → TDMA (`tridiagonal.py`,
   validate vs scipy) → drop `scipy` from `install_requires`.
4. Re-export CORK at top-level `__all__`.
5. `CLIMT_PURE_PYTHON=1` build path + `py3-none-any` wheel + CI job; attach to a
   GitHub release.
6. Flagship page `09-live-rce.qmd` (gray + CORK RE, per-band flux, CO₂ knob).
7. Ch.1 and Ch.6 retrofits.
8. Reusable-template appendix page.
9. Prove RCE end-to-end in a real browser; record page-weight budget.
