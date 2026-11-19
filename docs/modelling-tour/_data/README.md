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
| components | Cork LW, SlabSurface, SimpleBoundaryLayer, DryConvectiveAdjustment | + EmanuelConvectionPython, SurfaceHumidity, GridScaleCondensation |
| water vapour | none — a CO₂-only atmosphere | evolves, supplied by the surface |
| surface humidity | 0 | saturated at the surface's own temperature, every step (`stepping.SurfaceHumidity(1.0)`, page 9) |
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
shipped step counts (dry ≈ 3650 at 12 h) are larger than Task 0's cold-start
convergence measurement: Task 0 stopped at the first TOA crossing, which lands
mid-transient while the surface is still settling.

**The moist column has a gate of its own** (`MOIST_GATE` in the generator).
Two things are different about it. First, it is noisy: episodic convection
swings the instantaneous TOA by about 0.5 W m⁻², while its late drift is only
about 0.05 K a month, and 1000 steps at 5 min is 3.5 days. Second, it never
reaches TOA = 0. It settles with a steady imbalance of about −0.45 W m⁻² that
no amount of stepping removes. About half of it is `DryConvectiveAdjustment`,
which conserves enthalpy with a moist heat capacity while the rest of the stack
counts dry air's: +0.23 W m⁻² of the source, against +0.15 from
`EmanuelConvectionPython` (`scripts/experiments/tour_page12_measurements.py
energy`; page 12 explains it). So the moist gate works on trends over 60-day
(17 280-step) windows. The mean surface temperature and the mean TOA must each
match the previous window's, to 0.01 K and 0.1 W m⁻², and the mean TOA must be
inside ±2 W m⁻². The file records the residual it settled at
(`window_mean_toa_w_m2`), and the residual test checks the state against that
rather than against zero.

**Until 2026-09-27 the moist files were made over a fixed surface specific
humidity of 0.015 kg kg⁻¹.** At the 280 K they settled at, that is 246 % of
saturation: a surface evaporating 140 W m⁻² from something wetter than water.
Both were regenerated over a saturated surface. The plan's Task 9 log has the
before and after.

**Since 2026-09-28 all three files are made with `SimpleBoundaryLayer`
diffusing dry static energy** (its new default; `HISTORY.rst`). Diffusing
temperature had left a stable layer near the ground that held Emanuel's
cloud-base closure negative, so grid-scale condensation did the raining and
the moist column's lapse rate stayed near dry adiabatic. With dry static energy
Emanuel does all the raining and the column above its boundary layer lies on
its moist adiabat. The moist column moved from 285.99 K to 279.79 K, and the
dry one from 266.48 K to 266.62 K.

### `rce_moist_2xco2_equilibrium.npz` — page 12's 2×CO₂, run offline

`rce_moist_equilibrium.npz` with CO₂ doubled to 660 ppm and re-converged under
the same gate, same components, same 5 min timestep. Page 12 loads it beside
the base state rather than running the experiment. It took 157 800 steps
(≈ 548 simulated days), and at the ~0.16 s a step the moist column costs in
the browser that is about seven hours. Page 11's dry 2×CO₂ still runs live:
1000 steps at 12 h, about two minutes.

Its provenance records what it was perturbed from (`perturbed_from`,
`perturbed_from_saved_at`, the base's surface temperature) and the warming,
both file to file (`warming_k`) and between the two gates' last-window means
(`window_mean_warming_k`).
`test_shipped_2xco2_state_is_the_shipped_moist_state_perturbed` checks that
stamp against the shipped base. Regenerating the base without this file then
fails a test instead of quietly changing the number page 12 quotes.

**The warming.** Stepped on for 60 000 steps, the base settles at
279.78 K (TOA −0.42) and the 2×CO₂ state at 281.90 K (TOA −0.28), within 0.02 K
of their files. **Page 12 quotes the settled difference, +2.12 K.** The file
difference, +2.10 K, agrees with it; the old fixed-humidity files did not,
because the old gate stopped them wherever TOA first crossed ±0.5.
`scripts/experiments/tour_rce_shipped_remeasure.py settle-moist` reproduces it.

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
