# Handoff — in-browser non-grey RCE demo

**Branch:** `feature/pyodide-cork-prep` (PR #217 open; local commits **unpushed**).
**Env:** run all Python/pytest in the `climt` conda env.
**Plan/spec:** `docs/superpowers/plans/2026-07-19-in-browser-nongrey-rce-demo.md`.

## Where things stand

Phase A (pure-Python packaging, Tasks 1–8) was already committed before this
stretch. This session built **Task 10/11**: the flagship in-browser page
`docs/radiative-transfer/09-live-rce.qmd` (quarto-live + Pyodide) — a gray
column (`single_band_gray_lw`) vs a non-grey column (`earth_low_res_lw`),
each driven to radiative equilibrium and now showing **live-updating 3-panel
diagnostics** (temperature, LW heating rate, up/down LW flux).

New this session (all committed):
- `single_band_gray_lw` table + generator + test — mimics climt's default gray
  LW optical depth (τ0=1 linear); verified CORK≈GrayLongwaveRadiation.
- Shared boot include `docs/_includes/climt-live-setup.qmd` with
  `integrate_to_equilibrium` (sync, native-tested twin in
  `docs/radiative-transfer/_live/rce_helpers.py`) and `animate_to_equilibrium`
  (async, browser-only, the live multi-panel plot).
- `docs/radiative-transfer/_live/serve_wheel.py` — CORS static server for local
  preview.
- CI: prerelease guard on `release_climt.yml` merged to develop (PR #218).

## IMMEDIATE next step (was mid-verification when interrupted)

Just fixed a **`TypeError: unhashable type: 'pyodide.ffi.JsProxy'`**: the gray/
non-grey cells ended on `state = await animate_to_equilibrium(...)`, and
quarto-live hashes each cell's final value — a sympl state marshals to an
unhashable JsProxy. Fixed & committed (`8454656`): both cells now end on a
`print(f"surface temperature: … K")` (returns None). Page re-rendered. All demo
work is committed; branch is clean apart from this untracked handoff.

**Verify:** reload `http://127.0.0.1:7655/radiative-transfer/09-live-rce.html`
in a browser and confirm (a) no JsProxy error, (b) the gray cell auto-runs and
shows the 3-panel figure, (c) the `surface temperature: … K` print appears.
This was the exact step in progress when the session ended.

## Open question (the reason for the reload)

Does the animation actually update **frame-by-frame live**, or only show the
final frame? The 3-panel figure renders (verified: canvas 949×362), but I could
NOT confirm live frames via automation — the WebGL canvas is opaque to pixel
inspection and synthetic click events don't re-trigger quarto-live cells. **Ask
the user to eyeball it.** If it only shows the final frame, switch to a
guaranteed-live approach: render each frame to a PNG and swap an `<img>` src
(independent of the interactive matplotlib backend).

## How to run locally

```sh
conda activate climt
CLIMT_PURE_PYTHON=1 python -m pip wheel . --no-deps -w /tmp/climt_wh
python docs/radiative-transfer/_live/serve_wheel.py /tmp/climt_wh 8912   # leave running
cd docs && quarto preview radiative-transfer/09-live-rce.qmd
```
(Both a CORS wheel server on :8912 and a site server were running this session.)

## Key gotchas (hard-won)

- **GitHub release assets are CORS-blocked** → micropip can't fetch them in the
  browser. Wheel must be on a CORS host (PyPI) or bundled same-origin. Current
  demo uses the localhost CORS server (`pyodide: packages:` → `127.0.0.1:8912`).
  A `web-demo-preview` prerelease was created and **deleted** (unusable + to keep
  it off PyPI). See memory `project_pyodide_demo_hosting`.
- **quarto-live autorun race:** don't `await micropip.install(...)` in a cell;
  install via the page front-matter `pyodide: packages:` (runs at setup).
- **Cell-final-value JsProxy:** never end a `{pyodide}` cell on a `def` or on a
  Python object (e.g. a sympl state); end on `print(...)`/None.
- **Flux-coupling update order:** in the RCE loop apply `state.update(new_state)`
  **before** `state.update(diagnostics)`, else the LW fluxes SlabSurface needs
  get clobbered → surface heats without bound.

## Physics decisions (from the user)

- Gray baseline = `single_band_gray_lw` (Frierson/default τ0=1), NOT the τ≈28
  `single_band_unit_lw`.
- Prescribe a fixed surface SW (240 W/m²), **no** shortwave/ozone component;
  5 m slab; `sympl.AdamsBashforth`. → an **isothermal** stratosphere (skin temp),
  not Earth's ozone inversion (explained on the page).
- Non-grey cell is **click-to-run** (its 14-band integration is ~2 min in
  Pyodide). Gray auto-runs.
- Next real climt release will be ~0.70.0 (not 0.20.0); PyPI/deployed hosting
  deferred to then.

## Deferred (plan Tasks 12–15)

Ch.1/Ch.6 live-cell retrofits, reusable-template appendix, formal native science
test (`tests/test_live_rce_demo.py`), end-to-end weight-budget. Also: decide
whether to push the branch (updates PR #217).
