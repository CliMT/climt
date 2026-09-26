# Modelling Tour data assets

Files here are **site assets**, not package data. `docs/_quarto.yml` lists
`modelling-tour/_data/*.npz` under `project: resources:` so Quarto publishes them
(it skips underscore-prefixed paths otherwise), and pages 1–3 list the table in
their `pyodide: resources:` front matter so quarto-live stages it into the
Pyodide filesystem at the same relative path. `_tour/tables.py` is what the pages
call to find it.

## `earth_spectrum_lw.npz`

A 56-band longwave correlated-k table for Earth: the same physics as the shipped
`earth_low_res_lw`, at four times the spectral resolution and half the CO₂ axis.

| | shipped `earth_low_res_lw` | `earth_spectrum_lw` |
|---|---|---|
| bands | 14 | 56 (each shipped band split in 4) |
| g-points per band | 8 | 8 |
| CO₂ axis | 10 nodes, 10–10 000 ppm | 5 nodes, same range |
| (T, p, X_H₂O) grid | 12 × 8 × 7 | identical |
| H₂O continuum | decoupled, band-grey | identical |
| size | 2.7 MB (in the wheel) | 5.6 MB (site asset) |

It exists because 14 bars do not look like a spectrum and a 180 cm⁻¹ window band
cannot show the window narrowing. It is **not** a better table for quantitative
work in general — it is a finer one, built from the same line data, and the
places where the two disagree are documented in `tests/test_spectrum_table.py`.

The band edges are a strict **refinement** of the shipped grid, so aggregating
the 56 bands back to 14 is exact summation. That is what makes the validation
gate a real comparison rather than an interpolation exercise.

- sha256 (`.npz`): `7cede8082e1b54e418946f22263bc4531099b7582990631d5d4c09c0a61d58fa`
- sha256 (`.nc`, not committed): `ff8e201b2b7ae840c9f3f18e297041b920a63051e08b44e448cbe62be7c5000f`
- source: `linepyline:earth_hifi`, HITRAN 2024 + MT_CKD 4.3, pseudovoigt,
  dnu = 0.1 cm⁻¹

### `rce_moist_2xco2_equilibrium.npz` — page 12's 2×CO₂, run offline

`rce_moist_equilibrium.npz` with CO₂ doubled to 660 ppm and re-converged under
the same two gates, same components, same 5 min timestep. Page 12 loads it
beside the base state rather than running the experiment: it took 73 600 steps
(≈ 256 simulated days), which at the ~0.19 s a step the moist column costs in
the browser is nearly four hours. Page 11's dry 2×CO₂ still runs live — 1000
steps at 12 h, about three minutes.

Its provenance records what it was perturbed from (`perturbed_from`,
`perturbed_from_saved_at`, the base's surface temperature) and the warming
(`warming_k`). `test_shipped_2xco2_state_is_the_shipped_moist_state_perturbed`
checks that stamp against the shipped base, so regenerating the base without
this file fails a test instead of quietly changing the number page 12 quotes.

**What the warming does and does not mean.** The moist column never
reaches TOA = 0. `EmanuelConvectionPython` is not fully energy-conserving
(a known property of the scheme), so the column settles with a steady TOA
imbalance of about +0.3 W m⁻² that no amount of stepping removes. The
|TOA| < 0.5 gate stops a spin-up where TOA first dips inside that band, which
is not where the column settles. Stepped on for 60 000 steps, the base settles
0.10 K below its file (279.86 K, TOA +0.31) and the 2×CO₂ state 0.07 K above
its own (281.68 K, TOA +0.35). **Page 12 quotes the settled difference,
+1.82 K**, not the +1.647 K between the files.
`scripts/experiments/tour_rce_shipped_remeasure.py settle-moist` reproduces it.

### Regenerating

The generation step needs the `linepyline` conda env (it owns the HITRAN line
data); the conversion needs `climt` (it goes through climt's netCDF reader).
Takes about ten minutes on the reference machine.

```sh
conda run -n linepyline python scripts/generate_tour_spectrum_table.py \
    --output /tmp/earth_spectrum_lw.nc --ngpt 8 --nsub 4 --co2-nodes 5
conda run -n climt python scripts/convert_ck_table_to_npz.py /tmp/earth_spectrum_lw.nc
cp /tmp/earth_spectrum_lw.npz docs/modelling-tour/_data/
conda run -n climt python -m pytest tests/test_spectrum_table.py -m slow
```

`--ngpt 4 --nsub 7` (98 bands, 5.2 MB) was tried first and is equally valid on
every test except the per-band OLR one; see the log in
`docs/superpowers/plans/2026-08-12-modelling-tour-radiation.md`, Task 13.

## `rce_dry_equilibrium.npz`, `rce_moist_equilibrium.npz` and `rce_moist_2xco2_equilibrium.npz`

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
| timestep | 12 h | 5 min (the Emanuel convention) |
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

**What "converged" means here.** The generator runs until the top-of-atmosphere
imbalance is under 0.5 W m⁻² *and* the surface temperature has stopped moving.
The surface *flux* imbalance is not usable as a gate — the one-step flux lag
`_tour/stepping.py` documents makes it oscillate several W m⁻² step to step even
at equilibrium — so the surface criterion is that the mean surface temperature
over one 1000-step window equals the mean over the previous one. That is a
strict form of the drift the residual test below checks, and it is why the
shipped step counts (dry ≈ 3500 at 12 h) are larger than Task 0's cold-start
convergence measurement: Task 0 stopped at the first TOA crossing, which lands
mid-transient while the surface is still settling.

### Regenerating

Only when `tests/test_modelling_tour.py::test_shipped_equilibrium_is_still_an_equilibrium`
fails — that test is the staleness guard, deliberately in place of the
dependency-hash machinery in `scripts/build_experiments.py`, which cannot tell
a no-op in `cork/` from a physics change and has charged for the difference.

```sh
conda run -n climt python scripts/generate_tour_equilibria.py
conda run -n climt python -m pytest tests/test_modelling_tour.py -k shipped
```

With no arguments it writes all three files, the 2×CO₂ state last because it
starts from the moist state just written. `--moist-2xco2` alone rebuilds only
the 2×CO₂ state from the moist file already there. Set `NUMBA_NUM_THREADS=1`
when running several of these at once: a single column gains nothing from
numba's threads, and oversubscribed cores made each step ~10× slower.

The script's defaults are exactly what shipped. After regenerating, re-check
every number pages 11 and 12 quote — the state moved, so they may have too.
