# Picket Fence LW: CO2 15 μm Band Refinement Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refine the Earth LW picket-fence band partition inside the 500–800 cm⁻¹ CO2 15 μm region so the line-core opacity is no longer diluted by band-averaging with the wings, and validate against RRTMG cool-to-space and radiative-equilibrium diagnostics.

**Architecture:** The PF component already reads `num_bands` from the table's `k_coefficients` shape (`climt/_components/picket_fence/lw/component.py:42`), so no component code needs to change. The fix is entirely in table generation + a new shipped `.nc` + diagnostic validation. Two-stretch quadrature (already wired into `scripts/generate_picket_fence_tables.py`) stays on and composes with the new band partition. All changes are additive — new files alongside existing ones — so rollback is delete-only.

**Tech Stack:** numpy, exo_k (in `radiation` conda env), climt (in `climt` conda env), matplotlib.

**Background (read before starting):**
- Prior diagnostic work in this session established: at isothermal 250 K with q=0, 376 ppm CO2, the shipped PF tables (`earth_low_res_lw_8gpt.nc` with band edges `10,500,800,1800,3250` and Gauss-Legendre quadrature) reproduce OLR correctly but collapse cool-to-space at p < ~10 hPa. Two-stretch quadrature with g_split sweep (0.90 → 0.99) only partially closes the gap, with the same kink shape persisting at all g_split values. Conclusion: the 500–800 cm⁻¹ band is too wide — one k-distribution can't represent both the 667 cm⁻¹ line cores and the wings without diluting the cores.
- Existing diagnostic scripts: `examples/cooling_to_space_T_sweep.py`, `examples/cooling_to_space_gsplit_sweep.py`, `examples/radiative_equilibrium_pf_convergence.py`.
- Chaverot k-table source: `/Users/joymonteiro/github/chaverot/data/correlated-k_tables/Earth_376ppmCO2_R500_10-30k.ktable.SI.h5` (see memory `reference_chaverot_ktables.md`).
- Conda envs: use `climt` for everything except table generation; use `radiation` (has `exo_k`) for `generate_picket_fence_tables.py`. See memory `feedback_conda_env.md`.
- RRTMG-LW reference: splits the CO2 15 μm region at 630 and 700 cm⁻¹ (band 3 = 500–630, band 4 = 630–700, band 5 = 700–820). This is the canonical answer for how fine the split needs to be.

---

### Task 1: Pick the candidate band partition and generate the test table

**Files:**
- Modify: none yet — we use existing CLI flags on `scripts/generate_picket_fence_tables.py`.
- Create: `climt/_data/picket_fence/correlated_k/earth_low_res_lw_co2refined.nc`

**Candidate band edges:** `10, 500, 630, 700, 800, 1800, 3250` (6 bands).
- Splits at 630 and 700 isolate the CO2 fundamental Q/R-branch core (630–700) from the P-branch wing (500–630) and the long-wavelength window edge (700–800).
- This matches RRTMG's split points around the band center, the minimum refinement that physically captures the issue.

**Quadrature:** Reuse the two-stretch quadrature with `g_split=0.97` and `ngpt=8`. Composition test: if cool-to-space now reaches RRTMG-like at low p, the band split was the missing ingredient.

- [x] **Step 1: Generate the refined-band test table**

Run from repo root `/Users/joymonteiro/github/climt`:
```bash
conda run -n radiation python scripts/generate_picket_fence_tables.py \
  --input /Users/joymonteiro/github/chaverot/data/correlated-k_tables/Earth_376ppmCO2_R500_10-30k.ktable.SI.h5 \
  --output climt/_data/picket_fence/correlated_k/earth_low_res_lw_co2refined.nc \
  --kind lw \
  --bands 10,500,630,700,800,1800,3250 \
  --ngpt 8 \
  --quadrature two-stretch \
  --g-split 0.97 \
  --mixture-molar-mass 0.028964
```

Expected output (last line): `wrote .../earth_low_res_lw_co2refined.nc  shape=(6, 8, 12, 8, 7) has_h2o_axis=True (k x2.079e+25 m^2/kg per m^2/molec)` — note `shape[0] == 6` for the 6 bands.

- [x] **Step 2: Verify the file exists and has the expected dimensions**

Run:
```bash
conda run -n climt python -c "
import numpy as np
from scipy.io import netcdf_file
nc = netcdf_file('climt/_data/picket_fence/correlated_k/earth_low_res_lw_co2refined.nc', 'r')
k = nc.variables['k_coefficients'].data
edges = nc.variables['band_edges'].data if 'band_edges' in nc.variables else None
print('k shape:', k.shape)
print('band edges:', edges)
print('gpt weights shape:', nc.variables['gpoint_weights'].data.shape)
nc.close()
"
```

Expected: `k shape: (1, 6, 8, 12, 8, 7)` (ngas, nband, ngpt, nT, nP, nX); 6 band edges shown; `gpoint_weights shape: (6, 8)`. If the band-edges variable has a different name, list `nc.variables.keys()` and adapt — do NOT change the table file.

- [x] **Step 3: Smoke-test that `PicketFenceLongwave` loads the new table**

Run:
```bash
conda run -n climt python -c "
from climt._components.picket_fence import PicketFenceLongwave
pf = PicketFenceLongwave(optics='correlated_k', table='earth_low_res_lw_co2refined')
print('num_bands:', pf._num_bands)
print('num_gpts:', pf._num_gpts)
"
```

Expected: `num_bands: 6`, `num_gpts: 8`. If the constructor errors, the loader has a band-count assumption — STOP and report; do not attempt a fix without going back to brainstorming.

- [x] **Step 4: Commit the new table**

```bash
cd /Users/joymonteiro/github/climt
git add climt/_data/picket_fence/correlated_k/earth_low_res_lw_co2refined.nc
git commit -m "feat(pf): add CO2-refined Earth LW k-table (6 bands, 2-stretch ngpt=8)"
```

---

### Task 2: Validate against the cool-to-space diagnostic

**Files:**
- Modify: `examples/cooling_to_space_gsplit_sweep.py` — replace the g_split sweep with a comparison: RRTMG vs shipped GL vs best two-stretch (g=0.97) vs refined-band-plus-two-stretch.

- [x] **Step 1: Edit the diagnostic script**

In `examples/cooling_to_space_gsplit_sweep.py`, replace the `LABELS`, `TABLES`, `COLORS`, `STYLES` block with:

```python
LABELS = [
    "RRTMG (16 bands)",
    "PF 4-band GL (shipped)",
    "PF 4-band 2s g=0.97",
    "PF 6-band 2s g=0.97 (refined)",
]
TABLES = {
    "PF 4-band GL (shipped)":      "earth_low_res_lw_8gpt",
    "PF 4-band 2s g=0.97":         "earth_low_res_lw_8gpt_2s097",
    "PF 6-band 2s g=0.97 (refined)": "earth_low_res_lw_co2refined",
}
COLORS = ["black", "#999999", "#ff7f0e", "#d62728"]
STYLES = ["-", "--", "-", "-"]
```

Also rename the output path on the `fig.savefig` line to `cooling_to_space_band_refinement.png` so we don't overwrite the gsplit plot:

```python
out_path = os.path.join(os.path.dirname(__file__), "cooling_to_space_band_refinement.png")
```

And update the title:

```python
ax.set_title(
    f"Cool-to-space: 4-band vs CO2-refined 6-band PF — T_iso = {T_ISO:.0f} K"
)
```

- [x] **Step 2: Run the diagnostic**

```bash
cd /Users/joymonteiro/github/climt
conda run -n climt python examples/cooling_to_space_gsplit_sweep.py
```

Expected stdout: 4 rows of OLR (~221.4 W/m²) and max|HR| values; OLR for all 3 PF variants should match within ~0.1 W/m² (table dimension change doesn't break total opacity).

Expected plot (`examples/cooling_to_space_band_refinement.png`): the refined-band red curve should be visibly closer to RRTMG above ~10 hPa than the 4-band orange curve — specifically, the "kink at 5 hPa" that persisted across all g_split values should be reduced or shifted higher. If the kink is *gone*, the band split was sufficient. If it's just moved, deeper refinement is needed (see Task 4 spike).

- [x] **Step 3: Eyeball acceptance**

Open `examples/cooling_to_space_band_refinement.png`. Three possible outcomes; record which one you see in your task report:

- **A. Strong improvement:** refined curve tracks RRTMG within ~20% across all pressures, no kink. → Proceed to Task 3.
- **B. Partial improvement:** kink shrinks but doesn't disappear. → Proceed to Task 3 anyway; Task 4 (spike for finer split) is recommended afterward.
- **C. No improvement:** refined curve ≈ 4-band curve. → STOP and report. Means our hypothesis is wrong; this needs reinvestigation, not more iteration.

Do not edit the script further to "fix" the plot. The diagnostic is honest evidence.

- [x] **Step 4: Commit**

```bash
git add examples/cooling_to_space_gsplit_sweep.py
git commit -m "test(pf): switch cool-to-space diag to compare 4-band vs CO2-refined 6-band"
```

---

### Task 3: Validate against radiative equilibrium

**Files:**
- Modify: `examples/radiative_equilibrium_pf_convergence.py` — replace the g-point-convergence comparison with a band-structure comparison.

- [x] **Step 1: Edit the script header constants**

In `examples/radiative_equilibrium_pf_convergence.py`, replace the `LABELS`, `COLORS`, `STYLES`, `HR_KEY` block with:

```python
LABELS = ["RRTMG (16 bands)", "PF 4-band GL", "PF 4-band 2s", "PF 6-band 2s (refined)"]
COLORS = ["black", "#999999", "#ff7f0e", "#d62728"]
STYLES = ["-", "--", "--", "-"]
HR_KEY = {
    "RRTMG (16 bands)":          "air_temperature_tendency_from_longwave",
    "PF 4-band GL":              "longwave_heating_rate",
    "PF 4-band 2s":              "longwave_heating_rate",
    "PF 6-band 2s (refined)":    "longwave_heating_rate",
}
```

- [x] **Step 2: Edit the component instantiation block**

Find the three `pf_2gpt = ...`, `pf_4gpt = ...`, `pf_8gpt = ...` lines (around `examples/radiative_equilibrium_pf_convergence.py:67-69`) and replace with:

```python
pf_gl   = PicketFenceLongwave(optics="correlated_k", table="earth_low_res_lw_8gpt")
pf_2s   = PicketFenceLongwave(optics="correlated_k", table="earth_low_res_lw_8gpt_2s097")
pf_ref  = PicketFenceLongwave(optics="correlated_k", table="earth_low_res_lw_co2refined")
```

Replace the `surface_pf_2/4/8` and `pf_state_2/4/8` blocks with three corresponding renames (`surface_pf_gl`, `surface_pf_2s`, `surface_pf_ref` and `pf_state_gl`, `pf_state_2s`, `pf_state_ref`).

Update the `for st in [...]` loops and the `steppers`/`states` dicts to use these new names with the new `LABELS`.

- [x] **Step 3: Run radiative equilibrium**

```bash
cd /Users/joymonteiro/github/climt
conda run -n climt python examples/radiative_equilibrium_pf_convergence.py
```

This runs `N_STEPS = 20000` steps and takes a few minutes. Expected: final TOA balance shows OLR ≈ 240 W/m² for all four runs (already verified for shipped tables). The temperature panel should show: refined-band T(p) closer to RRTMG than 4-band variants, especially in the stratosphere (above ~30 hPa). Surface T differences between PF variants and RRTMG should also shrink.

- [x] **Step 4: Record the result**

Note the printed `T_sfc` values for all 4 runs. The shipped 4-band 8-gpt got `T_sfc = 272.12 K` vs RRTMG `268.15 K` — a 4 K bias. Refined-band should reduce this bias. Write the numbers to the task report.

Also inspect `radiative_equilibrium_pf_convergence.png`: stratospheric T for refined-band should sit between RRTMG (~180 K at 1 hPa) and shipped 4-band (~120 K). If refined-band's strat T is now within ~10–20 K of RRTMG, that's success.

- [x] **Step 5: Commit**

```bash
git add examples/radiative_equilibrium_pf_convergence.py
git commit -m "test(pf): RCE comparison — 4-band vs CO2-refined 6-band Earth LW"
```

---

### Task 4 (CONDITIONAL — only run if Task 2 result was B or C): Deeper-split spike

If the 630/700 split helped but didn't close the gap, try one finer split before recommending broader changes. This is a research spike, not a deliverable — generate, test, and either keep or revert based on the result.

**Files:**
- Create: `climt/_data/picket_fence/correlated_k/earth_low_res_lw_co2refined7b.nc`

- [x] **Step 1: Generate a 7-band table**

```bash
cd /Users/joymonteiro/github/climt
conda run -n radiation python scripts/generate_picket_fence_tables.py \
  --input /Users/joymonteiro/github/chaverot/data/correlated-k_tables/Earth_376ppmCO2_R500_10-30k.ktable.SI.h5 \
  --output climt/_data/picket_fence/correlated_k/earth_low_res_lw_co2refined7b.nc \
  --kind lw \
  --bands 10,500,600,660,720,820,1800,3250 \
  --ngpt 8 \
  --quadrature two-stretch \
  --g-split 0.97 \
  --mixture-molar-mass 0.028964
```

Splits at 600/660/720/820 isolate a tight 660–720 cm⁻¹ "line core" band. Expect `shape[1] == 7`.

- [x] **Step 2: Add it to the cool-to-space diagnostic**

In `examples/cooling_to_space_gsplit_sweep.py`, append to `LABELS` and `TABLES`:

```python
LABELS.append("PF 7-band 2s g=0.97 (deeper)")
TABLES["PF 7-band 2s g=0.97 (deeper)"] = "earth_low_res_lw_co2refined7b"
COLORS.append("#9467bd")
STYLES.append(":")
```

- [x] **Step 3: Re-run the cool-to-space diagnostic**

```bash
conda run -n climt python examples/cooling_to_space_gsplit_sweep.py
```

Compare the 7-band purple curve to the 6-band red curve. If the gap to RRTMG is materially smaller (e.g., max|HR| within 0.2 K/day of RRTMG at T_iso=250 K), the deeper split is worth promoting. If 6-band and 7-band look essentially identical, the residual gap is not a band-count problem — it's likely the underlying Chaverot k-table's spectral resolution or the picket-fence integrator itself, and we should stop here and reconvene.

- [x] **Step 4: Record finding; do not commit unless 7-band is clearly better**

If 7-band is clearly better: commit both the new .nc and the diagnostic update with `git commit -m "feat(pf): add 7-band CO2-refined Earth LW table"`. If not: delete the .nc with `rm climt/_data/picket_fence/correlated_k/earth_low_res_lw_co2refined7b.nc` and revert the diagnostic edit. Either way, write the conclusion to the task report so we know what we tried.

---

### Task 5: Cleanup decision

Once Tasks 1–3 are done (and Task 4 if it ran), summarize what to keep vs. throw away.

- [ ] **Step 1: Inventory all the test tables**

```bash
ls -la /Users/joymonteiro/github/climt/climt/_data/picket_fence/correlated_k/earth_low_res_lw_8gpt_2s*.nc /Users/joymonteiro/github/climt/climt/_data/picket_fence/correlated_k/earth_low_res_lw_8gpt_twostretch.nc
```

These are leftovers from the g_split exploration (`_2s090`, `_2s095`, `_2s097`, `_2s099`, `_twostretch`). Five files. They're useful for re-running the g_split-sweep diagnostic but otherwise dead.

- [ ] **Step 2: Decide and commit**

Either:
- **Keep:** if we want future work to be able to re-verify the g_split sweep conclusion, leave them in place — they're <100 kB each.
- **Discard:** delete all five and commit the deletion. The g_split-sweep script will still exist as documentation of the procedure even without the data.

Whichever you choose, do not also delete `earth_low_res_lw_co2refined.nc` — that one (and possibly `_co2refined7b.nc`) is the deliverable.

Write the recommendation in the task report and let the user make the call before committing the cleanup.

---

## Self-Review Notes

- All 5 tasks have concrete files, commands, and expected outputs. No placeholders.
- Task 4 is correctly marked conditional on Task 2's outcome — explicit branching with stop conditions.
- The plan does not attempt to modify `PicketFenceLongwave` or any component code, because the component already reads band count from the table shape (confirmed at `climt/_components/picket_fence/lw/component.py:42`). If that assumption breaks during Task 1 Step 3, the plan instructs the worker to STOP and report, not to fix it inline.
- Quadrature change from prior session is preserved (already merged into `scripts/generate_picket_fence_tables.py`); this plan composes with it and does not re-touch it.
- Rollback path is simple: delete the new `.nc` file(s) and `git revert` the diagnostic/example edits.

---

## Experiment Log

> **Policy:** every experiment we run (diagnostic, opacity tweak, table refit, comparison) is logged below — what we tried, what we learned, what we missed, what we rejected. Negative results matter. Append, never overwrite.

### 2026-05-18 — Pure-LW radiative-equilibrium intercomparison (cold-stratosphere diagnostic)

**Hypothesis.** RRTMG settles to T_strat ≈ 190 K while all PF variants sit at 130–160 K in pure-LW radiative equilibrium (no SW heating). Initial guess: the picket-fence tables have insufficient LW opacity in the upper atmosphere.

**Method.**
- Ran `examples/radiative_equilibrium_pf_convergence.py` (modified to N_STEPS = 8000 ≈ 1000 days @ 3 h timestep, slab surface, SW = 240 W m⁻², CO₂ = 376 ppm, dry).
- Models: RRTMG (16 bands), PF 4-band shipped, PF 6-band 2s g=0.97 (CO₂-refined), PF 6-band GL (CO₂-refined), PF 10-band GL.
- Added `upwelling_longwave_flux_in_air_per_band` capture and a per-band brightness-temperature diagnostic at p ≈ 100 hPa using band-integrated Planck inversion. New plot: `examples/radiative_equilibrium_pf_per_band.png`.

**Result (per-band BT at 105 hPa vs simulated T_local):**

| Table | T_local | CO₂ core (band, BT) | Wing bands BT |
|---|---|---|---|
| PF 4-band shipped | 154 K | 500–800 cm⁻¹ lumped, BT = 233 K | n/a (lumped) |
| PF 6-band 2s g=0.97 | 157 K | 630–700, BT = 164 K | 500–630: 252 K, 700–800: 242 K |
| PF 6-band GL | 155 K | 630–700, BT = 164 K | 500–630: 251 K, 700–800: 243 K |
| PF 10-band GL | 149 K | 630–700, BT = 158 K | 500–630: 249 K, 700–820: 246 K |

All non-CO₂-core bands at 100 hPa carry brightness temperatures ≥ 230 K → essentially window-like (κ ≈ 0 aloft).

**What we learned.**
1. Picket-fence physics is internally consistent: T_local ≈ CO₂-core BT for all 6/10-band tables, exactly as the κ-weighted Planck-balance argument predicts. Cold stratosphere is *not* a bug in the integrator, time-stepping, or convergence — it's the equilibrium the spectral distribution wants.
2. The CO₂ core (630–700 cm⁻¹) is the *only* band constraining T_strat in the PF tables. Everything else is window-like at 100 hPa.
3. RRTMG's warmer stratosphere (~190 K) must come from the **CO₂ wings (500–630 and 700–820 cm⁻¹) carrying nontrivial opacity** in the stratosphere, raising their BT and pulling T_eq up via the κ-weighted balance. PF tables apparently have these wings as windows — opacity concentrated too tightly in the 630–700 core.
4. The 4-band shipped table lumps the CO₂ feature into a single 500–800 cm⁻¹ band; its BT is dragged up to 233 K by mixing with window flux, but the band κ is too weak to actually constrain T_strat — the layer ends up at 154 K (similar to the more refined runs).

**What we rejected.**
- "PF runs haven't converged yet." Rejected — Top-N-level T panel is flat, OLR plateaus at 240 W m⁻², TOA imbalance ≈ 0.
- "Missing stratospheric SW heating (O₃, NIR H₂O)." User confirmed: pure-LW comparison for *all* models including RRTMG, so this can't be the discriminator.
- "Broadband upwelling LW differs aloft." User noted: broadband upwelling profiles are similar across models; spectral distribution is what matters.

**Open question / next experiment.**
- Hypothesis to test: artificially boost CO₂-wing opacity (500–630 and 700–820 cm⁻¹) in the 6-band table by ~5–10× and rerun. If T_strat moves from ~155 K toward ~190 K, the wing-opacity diagnosis is confirmed and the fix path is "regenerate CO₂ k-tables with better wing representation" (e.g., more g-points distributed into the wings, or refit against LBL).

**Files touched.**
- `examples/radiative_equilibrium_pf_convergence.py` — N_STEPS to 8000, added per-band capture + diagnostic figure + brightness-T inversion helper.
- New diagnostic figure: `examples/radiative_equilibrium_pf_per_band.png`.

### 2026-05-18 — Wing-opacity boost experiment (REJECTED hypothesis)

**Hypothesis.** Multiplicatively boosting CO₂-wing k-coefficients (bands 1: 500–630, b3: 700–820 cm⁻¹) in `earth_low_res_lw_co2refined_gl` will warm T_strat toward RRTMG's ~190 K by activating the wings at 100 hPa.

**Method.** `examples/experiment_co2_wing_boost.py` copies the baseline table, multiplies `k_coefficients` in bands 1 and 3 by 1× / 3× / 10× / 30×, writes `*_wingboost{N}x.nc` (NETCDF3_CLASSIC, required by PF loader), and runs 1000 days of pure-LW RE for each.

**Result (1000 days, p ≈ 100 hPa):**

| variant | T_sfc | T_strat | b1 BT | b3 BT | b2 (core) BT |
|---|---|---|---|---|---|
| baseline | 272.0 | 155.0 | 250.9 | 242.6 | 163.8 |
| wing×3 | 273.9 | 153.6 | 245.7 | 235.3 | 162.2 |
| wing×10 | 276.0 | 152.5 | 238.6 | 226.2 | 161.3 |
| wing×30 | 278.1 | 152.2 | 230.8 | 216.5 | 161.1 |

OLR balanced at 240 W m⁻² in every case (TOA imbalance < 0.1).

**What we learned.**
1. Wing-band BTs at 100 hPa do drop monotonically with boost (b3: 243 → 217 K at 30×), confirming added opacity is reaching upper levels of the column.
2. But T_strat barely moves and actually *cools slightly* (155 → 152 K).
3. T_sfc warms ~6 K from baseline to 30× — the boost is operating as expected as a tropospheric greenhouse enhancement.
4. Mechanism: multiplicative k boost preserves the vertical *structure* of k(p). CO₂ wings have k concentrated at high pressure (collision broadening). Boosting by 30× still leaves κ_wing(100 hPa) ≪ κ_core(100 hPa), so the stratospheric κ-weighted balance remains pinned by the core. The wings' τ=1 surface migrates from ~500 hPa up to ~200 hPa — into the upper troposphere, not the stratosphere.
5. Slight stratospheric cooling explanation: warmer troposphere with more opacity reduces downwelling reaching 100 hPa in the now-more-absorbing wings; layer loses a previously-present heat source.

**What we rejected.**
- "Multiplicative wing boost will warm T_strat." Rejected: T_strat goes the *wrong* way.
- "Total wing-band column opacity is the controlling variable." Rejected: 30× more total wing opacity, virtually no T_strat change.

**Reframed hypothesis (to test next).**
The issue is the **vertical structure of k(p), not the total column k**. The wings need k that doesn't collapse at low pressure — i.e., line strength at stratospheric p that's a larger fraction of surface-level k than the current Chaverot rebinning provides.

**Candidate next experiments.**
1. **Pressure-selective boost.** Boost wing k only at low pressure (indices p ≲ 100 hPa). If T_strat now warms toward 190 K, the diagnosis is "k(p) collapses too fast aloft." If still flat, the wing story is fully dead.
2. **Inspect k(p) directly.** Plot k(p, T=fixed) in the wing bands vs core band, and compare to RRTMG band coefficients at the same pressures. Quantify the "wing-collapse" at low p.
3. **Core de-boost.** Try the opposite: divide core (b2) k by 3, 10, 30. If T_strat warms, the diagnosis flips — it's the core being too strong, not the wings being too weak.
4. **Look beyond CO₂.** RRTMG includes O₃ 9.6 μm and H₂O continuum LW; even at zero H₂O specific humidity, the O₃ band could be contributing stratospheric warming. The PF table at 376 ppm CO₂ + 0 ppm O₃ may simply lack a critical absorber. Worth checking what gases the Chaverot table has baked into the "premixed" background.

**Files touched.**
- `examples/experiment_co2_wing_boost.py` (new).
- New diagnostic figure: `examples/experiment_co2_wing_boost.png`.
- Throwaway `.nc` files in `climt/_data/picket_fence/correlated_k/`: `earth_low_res_lw_co2refined_gl_wingboost{3,10,30}x.nc` — do not commit; delete after experiment series concludes.

### 2026-05-18 — Direct opacity inspection (quantitative confirmation)

**Hypothesis (clarified).** The PF k-tables over-concentrate CO₂ opacity into the 630–700 cm⁻¹ core (τ_total ≫ 1000) while leaving the wings (500–630, 700–820 cm⁻¹) with too little total opacity to saturate anywhere except the lower troposphere.

**Method.** `examples/diagnose_pf_opacity_structure.py` equilibrates the column and then makes one forward call to capture `longwave_optical_depth_per_band` directly from the PF component diagnostic (no offline interp — guaranteed self-consistent with the production calc). Computes cumulative τ(TOA → p) per band and finds p(τ=1).

**Result (1000-day equilibrated, 6-band-GL):**

| Band | ν range (cm⁻¹) | τ(sfc) | p(τ=1) |
|---|---|---|---|
| b0 H₂O rot tail | 10–500 | 1.5 | 741 hPa |
| b1 CO₂ left wing | 500–630 | 3.1 | 469 hPa |
| **b2 CO₂ core** | **630–700** | **1011** | **4.4 hPa** |
| b3 CO₂ right wing | 700–800 | 6.9 | 358 hPa |
| b4 window | 800–1800 | 0.15 | never |
| b5 CO₂ 4.3 μm | 1800–3250 | 212 | 30 hPa |

10-band table: same picture; b3 (630–700) τ_total = 988, p(τ=1) = 4.4 hPa. Wings b2 (500–630) τ_total = 2.8, b4 (700–820) τ_total = 4.8.

**What we learned (quantitative).**
1. **Core overpacked**: τ_total ≈ 1000 in 630–700 cm⁻¹, saturating at 4 hPa. This is the *only* LW band reaching τ=1 in the upper stratosphere (other than b5 4.3 μm which carries negligible Planck flux at cold T).
2. **Wings underpacked**: τ_total only 3–7 in the 500–630 and 700–820 cm⁻¹ wings. They saturate at 360–470 hPa (mid-troposphere). They can't contribute to stratospheric energy balance because their κ stays small at low p AND their total column opacity is small.
3. **The earlier wing-boost experiment now makes sense in retrospect**: even at 30× boost, wing τ_total goes from 5–7 to ~150–200, and p(τ=1) only migrates from ~400 hPa up to ~50 hPa. The wings don't reach the upper stratosphere where they'd need to be to compete with the core at 4 hPa.
4. **Diagnosis of the underlying cause**: the Chaverot k-table rebinning to 8 g-points per band likely keeps the strongest absorbers (core lines) but loses the moderate-strength absorbers that populate the wings between line peaks. With only 8 g-points per band, you can't represent both a strong-line core and a moderate-line wing within the same band without losing one.

**Reframed hypothesis.**
The discrepancy with RRTMG is not "PF tables have wrong total CO₂ opacity" — they likely have *too much* in the core and *too little* in the wings. The fix is structural: regenerate the tables with either (a) finer g-point distribution that captures moderate-strength wing lines, or (b) different band edges that don't lump the strong core into a band where it dominates everything.

**Candidate next experiments.**
1. **Compare to LBL**: take the existing `examples/compare_lbl_pf_band_k.py` diagnostic and quantify how badly the Chaverot-rebinned τ underestimates wing opacity vs an LBL benchmark at, say, 100 hPa / 250 K.
2. **Split the core**: try a 7-band table that splits 630–700 into a tight core (665–675 cm⁻¹) plus shoulder regions (630–665, 675–700) — this gives the rebinner a fighting chance to put moderate lines in the shoulders without competing with the strong core. This is similar to the Task 4 suggestion already in the plan.
3. **Low-p-only wing boost**: as a final validation of the structural argument, boost wing k by, say, 100× only at p < 50 hPa (table grid indices). If T_strat now jumps from 155 K to 180 K, we've definitively isolated the vertical structure as the controlling factor and the fix is a g-point redistribution that puts more weight on the low-p line strengths.

**Files touched.**
- `examples/diagnose_pf_opacity_structure.py` (new; first version was buggy offline interp, rewrote to use the component's own per-band τ diagnostic).
- `examples/diagnose_pf_opacity_structure.png` — opacity & cumulative-τ plots for both 6-band and 10-band tables.

### 2026-05-18 — Proposed next experiment: more g-points per band, smaller band count

**Hypothesis.** The stratospheric cold bias is a **g-point resolution problem within the CO₂ band**, not a band-edge problem. With only 8 g-points per band the rebinner can only resolve a handful of points on the within-band k-distribution; for a band that contains both strong-line core absorbers and moderate-line wing absorbers, the 8 nodes get concentrated at the strong-line tail (high g) and effectively miss the moderate-strength contributors. Result: the band saturates at one altitude (driven by the strongest g-points) but has very little opacity *elsewhere* in the column — exactly what we observe (τ_total = 1011, p(τ=1) = 4 hPa, then nothing).

If we keep the band structure simple and lump the full CO₂ feature into a single 500–820 cm⁻¹ band but increase ngpt to 16 or 32, the rebinner has the budget to place g-points across the moderate-strength portion of the CDF. Predicted effect: b_CO₂ still has very large τ_total, but p(τ=1) migrates *down* from 4 hPa to ~50–100 hPa where moderate g-points start dominating; T_strat warms toward 190 K.

This proposal is preferred over band-splitting because it preserves interpretability (4 physically motivated bands, not 7+ tightly-cut shards) and provides a clean axis to vary (ngpt) that decouples band structure from quadrature resolution.

**Plan.**

1. **Generate three new tables** (requires `radiation` conda env per memory `feedback_conda_env.md`; source HDF5 per memory `reference_chaverot_ktables.md`):

   Source file: `/Users/joymonteiro/github/chaverot/data/correlated-k_tables/Earth_376ppmCO2_R500_10-30k.ktable.SI.h5`

   ```bash
   # 4-band, 16 g-points per band (64 total)
   conda run -n radiation python scripts/generate_picket_fence_tables.py \
     --input /Users/joymonteiro/github/chaverot/data/correlated-k_tables/Earth_376ppmCO2_R500_10-30k.ktable.SI.h5 \
     --output climt/_data/picket_fence/correlated_k/earth_low_res_lw_4band_16gpt.nc \
     --kind lw \
     --bands 10,500,820,1800,3250 \
     --ngpt 16 \
     --quadrature two-stretch \
     --g-split 0.97 \
     --mixture-molar-mass 0.028964

   # 4-band, 32 g-points per band (128 total)
   conda run -n radiation python scripts/generate_picket_fence_tables.py \
     --input /Users/joymonteiro/github/chaverot/data/correlated-k_tables/Earth_376ppmCO2_R500_10-30k.ktable.SI.h5 \
     --output climt/_data/picket_fence/correlated_k/earth_low_res_lw_4band_32gpt.nc \
     --kind lw \
     --bands 10,500,820,1800,3250 \
     --ngpt 32 \
     --quadrature two-stretch \
     --g-split 0.97 \
     --mixture-molar-mass 0.028964

   # 6-band, 16 g-points per band (same edges as co2refined_gl, to isolate ngpt effect)
   conda run -n radiation python scripts/generate_picket_fence_tables.py \
     --input /Users/joymonteiro/github/chaverot/data/correlated-k_tables/Earth_376ppmCO2_R500_10-30k.ktable.SI.h5 \
     --output climt/_data/picket_fence/correlated_k/earth_low_res_lw_6band_16gpt.nc \
     --kind lw \
     --bands 10,500,630,700,800,1800,3250 \
     --ngpt 16 \
     --quadrature gauss-legendre \
     --mixture-molar-mass 0.028964
   ```

   Band edges chosen: 4-band uses 10/500/820/1800/3250 (rotational / full CO₂ 15μm / window / NIR). 6-band keeps the existing co2refined_gl edges 10/500/630/700/800/1800/3250 so the ngpt comparison is pure (same bands, different g-point count).

2. **Add the three new tables to the diagnostic and the convergence script.**
   - `examples/diagnose_pf_opacity_structure.py` → add the three new tables to `TABLES`. Compare per-band τ(sfc) and p(τ=1).
   - `examples/radiative_equilibrium_pf_convergence.py` → add the three new tables as labels with distinct colors. Run for 1000 days.

3. **Decision criteria.**
   - **Confirmed (ngpt-resolution diagnosis correct):** For the 4-band 16-gpt and especially 32-gpt tables, the CO₂ band p(τ=1) drops from 4 hPa toward 50–100 hPa AND T_strat warms toward 180–190 K.
   - **6-band 16-gpt as ablation:** if this *also* warms T_strat substantially, the effect is purely from ngpt and band-structure doesn't matter. If only the 4-band variants warm, lumping helps because the 4-band CO₂ band sees a wider, more diverse k-distribution that benefits more from extra g-points.
   - **Rejected (ngpt doesn't help):** if 32-gpt looks identical to 8-gpt, the problem is upstream in the Chaverot source — the k-distribution we're sampling already lacks the moderate wing absorbers, regardless of how finely we sample it. Next path then is to look at the source CDF directly via Exo_k.

4. **What to log either way.**
   - Per-band τ table for all three new variants (analogous to the 2026-05-18 opacity table).
   - Final T_strat, OLR, T_sfc for each.
   - For the 4-band tables: report p(τ=1) for the CO₂ band and how it compares to the 4 hPa we got with the 6-band co2refined_gl.
   - For the 6-band 16-gpt table: report whether p(τ=1) for the b2 core differs from the 8-gpt version, and by how much.

5. **Cleanup.**
   - The three new `.nc` files are candidates for shipping if the experiment succeeds. If they don't help, delete and revert.
   - Throwaway wing-boost `.nc` files from the prior experiment can be deleted now — they served their purpose.

**Files to be touched (when this runs).**
- New: `climt/_data/picket_fence/correlated_k/earth_low_res_lw_{4band_16gpt,4band_32gpt,6band_16gpt}.nc`
- Modify: `examples/diagnose_pf_opacity_structure.py`, `examples/radiative_equilibrium_pf_convergence.py`

### 2026-05-18 — Result: ngpt boost REJECTED

**Method (as proposed).** Generated three new LW tables from the Chaverot HDF5 source using `scripts/generate_picket_fence_tables.py`:
- `earth_low_res_lw_4band_16gpt.nc` — 4 bands (10/500/820/1800/3250), 16 g-pts, two-stretch g=0.97
- `earth_low_res_lw_4band_32gpt.nc` — same bands, 32 g-pts, two-stretch g=0.97
- `earth_low_res_lw_6band_16gpt.nc` — 6 bands (10/500/630/700/800/1800/3250), 16 g-pts, GL

Ran the per-band brightness-T diagnostic from `examples/radiative_equilibrium_pf_convergence.py` (1000 days, pure LW, SW=240 W/m², 376 ppm CO₂, q=0). Compared baselines (4-band 8 gpt shipped, 6-band 8 gpt GL, 10-band 8 gpt GL) against the three new high-ngpt variants.

**Result (per-band BT at p ≈ 105 hPa, T_local in K).**

| Table | ngpt | T_local | CO₂-core band BT | Wing bands BT |
|---|---|---|---|---|
| PF 4-band GL (shipped) | 8 | 153.9 | 500–800 lumped, 233.4 | n/a |
| PF 4-band 2s g=0.97 | **16** | 156.3 | 500–820 lumped, 236.0 | n/a |
| PF 4-band 2s g=0.97 | **32** | 155.8 | 500–820 lumped, 236.3 | n/a |
| PF 6-band GL | 8 | 155.0 | 630–700, **163.8** | 500–630: 250.9, 700–800: 242.6 |
| PF 6-band GL | **16** | 156.4 | 630–700, **164.6** | 500–630: 250.9, 700–800: 242.5 |
| PF 10-band GL | 8 | 148.9 | 630–700, **157.7** | 500–630: 249.1, 700–820: 245.6 |

**What we learned.**
1. **ngpt is not the controlling variable.** Going from 8 → 16 → 32 g-points in the 4-band CO₂ band shifts T_local by < 3 K (153.9 → 156.3 → 155.8) and band BT by < 3 K (233.4 → 236.0 → 236.3). The 4-band 32-gpt result is statistically indistinguishable from 4-band 16-gpt — additional g-points are sampling the same coarse k-distribution and bring nothing new.
2. **Same conclusion for the 6-band ablation.** Doubling 6-band GL from 8 → 16 g-pts moves the CO₂-core BT by 0.8 K (163.8 → 164.6) and T_local by 1.4 K. Band structure already isolates the core; more samples within that core don't add resolved absorbers.
3. **The 10-band GL run remains the coldest** (T_local = 148.9 K, core BT = 157.7 K), because splitting at finer band edges *concentrates* the core's strongest opacity into an even narrower band, dropping core BT further. This confirms the 2026-05-18 opacity-structure diagnosis: finer band/g-point sampling makes the problem **worse**, not better, because each new degree of freedom sharpens the strong-line core without ever sourcing wing opacity that doesn't exist in the input Chaverot CDF.

**What we rejected.**
- "More g-points per CO₂ band will distribute k across moderate-strength wing absorbers and warm T_strat." **Rejected** — 4× more g-points (8 → 32) yields essentially no T_strat change.
- "4-band lumping benefits more from extra g-points than narrow bands do" (the alternative branch of the proposal). **Rejected** — 4-band 16/32-gpt sits within 1 K of 6-band 16-gpt and within 3 K of every other variant; lumping vs splitting doesn't change the answer.

**Result (opacity structure, from `examples/diagnose_pf_opacity_structure.py`).**

| Table | ngpt | Key CO₂ band | τ(sfc) | p(τ=1) hPa | T_strat (col min) |
|---|---|---|---|---|---|
| PF 6-band GL | 8 | b2 630–700 (core) | 1011 | 4.38 | 117.1 |
| PF 6-band GL | **16** | b2 630–700 (core) | **1309** | **2.40** | 113.0 |
| PF 10-band GL | 8 | b3 630–700 (core) | 988 | 4.38 | 112.4 |
| PF 4-band 2s | 8 (shipped uses GL) | b1 500–800 (lumped) | — | — | (153.9 in RCE) |
| PF 4-band 2s g=0.97 | **16** | b1 500–820 (lumped) | 253 | 12.0 | 111.6 |
| PF 4-band 2s g=0.97 | **32** | b1 500–820 (lumped) | 261 | 11.9 | 111.2 |
| PF 6-band 2s g=0.97 | 8 (co2refined) | (not in this run) | — | — | — |

Wing-band τ(sfc) for 6-band variants: 500–630 stays at ~3–5, 700–800 at ~7–8 across 8-gpt and 16-gpt — **wings get nothing from more g-points**.

**What we learned from τ-structure (sharpens the BT-only finding).**
1. **More g-points pushes the core's τ=1 surface UP, not down.** 6-band core p(τ=1) goes from 4.4 hPa (8-gpt) to 2.4 hPa (16-gpt) and τ_sfc from 1011 to 1309. Extra nodes are sampling the strong-line tail of the same CDF more precisely; they pile on opacity in the line core (already optically thick everywhere) rather than discovering moderate-strength absorbers in the wings.
2. **4-band lumping pushes the τ=1 surface from 4 hPa to ~12 hPa** for the CO₂ band — the lumped band's effective opacity at low p is much weaker because the strong core sits in a band whose Planck weighting is dominated by window-region flux. But T_strat is the same (~111 K) because the core's *radiating temperature* is still pinned by the few strong-g-pt nodes saturating aloft.
3. **The 32-gpt → 16-gpt comparison is the cleanest control**: identical band structure, doubled g-pt count, τ_sfc rises 3% and p(τ=1) shifts 1%. Quantitatively dead.

**Diagnosis (now load-bearing).** The cold-stratosphere bias is **upstream of the picket-fence integrator**. The Chaverot R=500 source k-distribution itself lacks the moderate-strength CO₂-wing line opacity at low pressure that warms RRTMG's stratosphere. No band edge choice or g-point count can recover absorbers that aren't in the input CDF. The fix path lies outside this plan — either (a) regenerate Chaverot's source k-tables with explicit line-by-line wing treatment, or (b) source a different LW k-table (e.g. RRTMG's own, or RFMIP-style LBL) for the Earth PF configuration.

**Files touched.**
- New (kept; useful as reference for future investigators verifying the ngpt-null result): `climt/_data/picket_fence/correlated_k/earth_low_res_lw_4band_16gpt.nc`, `_4band_32gpt.nc`, `_6band_16gpt.nc`.
- Modified: `examples/diagnose_pf_opacity_structure.py` (added 3 new tables to `TABLES`), `examples/radiative_equilibrium_pf_convergence.py` (added 3 new labels + states + steppers).

**Next steps (out of this plan's scope).**
- Inspect the Chaverot HDF5 CDF directly in g-space at p ≲ 100 hPa for 500–630 and 700–820 cm⁻¹ to quantify the wing-absorber deficit.
- Evaluate alternative source k-tables (RRTMG's, or an LBL-derived set) for the picket-fence Earth LW path.
- Document this conclusion in `climt/_components/picket_fence/lw/component.py` or a `docs/` note so future contributors don't re-walk this path.

### 2026-05-18 — Follow-up experiment #1: Chaverot CDF wing vs core direct inspection

**Why.** The 2026-05-18 result diagnosed the cold-stratosphere bias as upstream of the integrator — the Chaverot R=500 source CDF itself lacks moderate-strength wing absorbers in 500–630 and 700–820 cm⁻¹. Before pursuing the larger fixes (regenerating Chaverot's CDF via LBL, or swapping in RRTMG's own k-tables), the smallest possible experiment is to look at the source CDF directly and quantify the deficit.

**Method.** New script `scripts/experiments/inspect_chaverot_cdf_lw_wings.py`:
- Opens `Earth_376ppmCO2_R500_10-30k.ktable.SI.h5` with `exo_k.Ktable5d` (radiation env, not climt).
- `bin_down([500, 630, 700, 820])` to expose the three diagnostic bands while preserving the native g-grid (Ng=20).
- Slices at T=230 K and X_H2O=1e-6 (nearest grid points to dry low-stratosphere target).
- Plots k(g) for all three bands at p = 1, 10, 100, 1000, 10000 hPa.
- Computes g-weighted band-mean k and prints wing/core ratios.
- Output: `debug_data/chaverot_cdf_lw_wings_kg.png` and stdout table.

**Result.**

| p [hPa] | wing 500–630 [m²/mol] | core 630–700 [m²/mol] | wing 700–820 [m²/mol] | wing/core (left) | wing/core (right) |
|---:|---:|---:|---:|---:|---:|
| 1     | 1.76e-30 | 4.64e-28 | 2.54e-30 | 3.79e-03 | 5.48e-03 |
| 10    | 7.43e-30 | 1.98e-27 | 1.20e-29 | 3.75e-03 | 6.09e-03 |
| 100   | 3.84e-29 | 8.05e-27 | 4.59e-29 | 4.77e-03 | 5.70e-03 |
| 1000  | 1.87e-29 | 4.64e-27 | 3.61e-29 | 4.04e-03 | 7.77e-03 |
| 10000 | 1.54e-29 | 3.91e-27 | 3.03e-29 | 3.95e-03 | 7.75e-03 |

**What we learned.**
1. **Wing-absorber-deficit diagnosis confirmed quantitatively.** The wing bands (500–630, 700–820 cm⁻¹) carry mean k that is 0.4–0.8% of the 630–700 core across the full diagnostic pressure range (1 hPa → 10000 hPa). This is the smoking gun: it's not an artifact of climt's re-binning, it's not an integrator bug — the moderate-strength wing opacity simply isn't in the input CDF.
2. **The deficit is ~constant in pressure.** Wing/core ratios stay near 0.5% from 1 hPa to 10000 hPa with no clear pressure-broadening signal that would let extra g-points or band splits recover wing flux at low p (which is where it would warm T_strat). At all pressures the wings are an order of ~200× weaker absorber than the core — no integrator scheme can extract stratospheric warming from absorbers that aren't there.
3. **The R=500 sampling is likely the root cause.** Chaverot's R=500 CDF samples coarsely enough across the wing region (~1.4 cm⁻¹ between line cores in the strongest CO₂ wing regions) that moderate-strength lines either average into the strong core or are missed entirely. RRTMG's 16-band LW uses LBL-derived k-distributions per-band; that's where its wing absorbers come from.

**What we rejected.**
- "Wing opacity is being eaten by re-binning inside generate_picket_fence_tables.py." **Rejected** — the raw Chaverot CDF, before any climt-side processing, already shows the wing/core deficit at the same magnitude as the shipped tables suggest.
- "Wing pressure-broadening at high p would rescue stratospheric warming if we used the right band split." **Rejected** — wing/core ratio is essentially flat in p; no pressure regime has wing k within 1% of core k.

**Next steps (in priority order, all still out of this plan's scope).**
1. **Try a different source k-table for the same Earth LW config.** RRTMG ships its own k-tables that we know warm the stratosphere correctly. Extracting them into a form `generate_picket_fence_tables.py` can ingest tests the diagnosis end-to-end and, if it works, gives a shipping path.
2. **Compare against an LBL benchmark.** Run linepyline (or an existing LBL spectrum in chaverot/) at p=100 hPa / T=250 K, integrate to band-mean k for the three bands, compare to Chaverot. Quantifies the deficit against ground truth and tells us whether the fix has to go all the way to LBL or whether a higher-R Chaverot table would suffice.
3. **Document the conclusion in `climt/_components/picket_fence/lw/component.py`** so future contributors don't re-walk this diagnostic path.

**Files touched.** New: `scripts/experiments/inspect_chaverot_cdf_lw_wings.py`, `debug_data/chaverot_cdf_lw_wings_kg.png`. No code in `climt/` changed.

### 2026-05-18 — Follow-up experiment #2: Chaverot vs LBL ground truth → previous diagnosis INVERTED

**Why.** Experiment #1's "wing-absorber-deficit" framing was hand-wavy. R=500 (Δν≈1.4 cm⁻¹ at 700 cm⁻¹) is more than enough to resolve CO₂ v2 rotational structure (line spacing 0.4–1.5 cm⁻¹), so a 200× core/wing contrast might just be the genuine spectrum at low T. Needed LBL ground truth.

**Method.** New script `scripts/experiments/compare_chaverot_vs_lbl_wings.py`. Same three bands as #1. Two sources:
- **Chaverot at native R=500**, linearly interpolated in T to 250 K between the bracketing grid nodes (230 and 290 K); g-integrated per native sub-band; bandwidth-weighted across sub-bands in each diagnostic window; converted m²/molec → m²/kg-air via N_A/M_air.
- **linepyline LBL** from the pre-existing `lbl_k_iso.nc` (T=250 K isothermal, 60 p levels, dnu=0.1 cm⁻¹, CO₂-only at 376 ppm in air, pseudovoigt). Uniform-dnu mean over each window — already in m²/kg-air.

Reports both as m²/kg-air, plus wing/core ratios per source and LBL/Chaverot per band. Plot: k(ν) overlay at p=100 hPa across 480–840 cm⁻¹.

**Result (band-mean k, m²/kg-air, T=250 K, CO₂=376 ppm, dry).**

| p [hPa] | source | wing 500–630 | core 630–700 | wing 700–820 | w/c (left) | w/c (right) |
|---:|:--|---:|---:|---:|---:|---:|
| 10   | Chaverot | 4.93e-4 | 7.89e-2 | 8.09e-4 | 0.62% | 1.03% |
| 10   | LBL      | 7.21e-4 | 2.11e-1 | 3.26e-4 | 0.34% | 0.15% |
| 100  | Chaverot | 8.75e-4 | 1.48e-1 | 1.35e-3 | 0.59% | 0.92% |
| 100  | LBL      | 4.82e-4 | 1.19e-1 | 8.15e-4 | 0.40% | 0.68% |
| 1000 | Chaverot | 8.67e-4 | 1.23e-1 | 1.41e-3 | 0.70% | 1.15% |
| 1000 | LBL      | 4.86e-4 | 9.73e-2 | 1.08e-3 | 0.50% | 1.11% |

LBL/Chaverot absolute ratio: core = {2.68, 0.81, 0.79} at p = {10, 100, 1000} hPa; wings range 0.40–1.46.

**What we learned.**
1. **The 200× core/wing contrast is REAL CO₂ physics, not a CDF artifact.** LBL and Chaverot agree on wing/core ratios to within a factor of 2 across the full pressure range. Both put wing/core well under 1.5% everywhere. The previous "wing-absorber-deficit" framing in #1 was wrong.
2. **Chaverot R=500 is a faithful representation of LBL** for this band at this composition. Absolute core agreement within 20% at p ≥ 100 hPa. The 2.7× core overshoot at 10 hPa is the only large disagreement — most plausibly a line-shape divergence at low p (pseudovoigt vs Chaverot's source shape and far-wing/continuum treatment). Not large enough to be the root cause of the cold-strat bias.
3. **The cold-stratosphere bias is NOT caused by missing wing opacity in the source CDF.** That hypothesis is ruled out.

**What we rejected.**
- "Chaverot R=500 lacks moderate-strength CO₂ wing absorbers." **Rejected by LBL.** Wing opacity is correctly low because the wings are physically weak at 250 K, not because the table is missing them.
- "Bringing RRTMG's k-tables in is the path forward." **Deprioritized** — there's no longer a defensible reason to believe the LW k-distribution is the bottleneck.

**Revised diagnosis (open hypotheses).**
1. **Missing absorbers other than CO₂/H₂O.** The Chaverot Earth LW table bakes in CO₂+H₂O. RRTMG's 16 LW bands include O₃ (strong 9.6 µm absorption that warms the stratosphere), CH₄, N₂O, and the H₂O continuum. A missing O₃ band is the most likely candidate for a stratospheric warming deficit. **Next check:** inspect the Chaverot HDF5 attributes to see exactly which gases are baked in.
2. **Band-partition / Planck-weighting interaction.** Even with faithful k, slicing the band differently moves p(τ=1) and the radiating T. The earlier 10-band run was *coldest* — consistent with this. But this is a secondary effect; a 5–10 K shift, not the 30+ K bias.
3. **The radiative-equilibrium target may genuinely be this cold.** Observed stratospheric T is set by radiation + dynamics; PF/RRTMG-pure-radiative answers could both be cold and only RRTMG's tuning + bigger gas set masks it.

**Files touched.** New: `scripts/experiments/compare_chaverot_vs_lbl_wings.py`, `debug_data/chaverot_vs_lbl_wings_kspec.png`. No code in `climt/` changed.

**Note on experiment #1's plan entry.** The previous entry's "Diagnosis confirmed quantitatively" claim is **superseded by this result**. The wing/core ratio in Chaverot is real, not a deficiency. Future readers: rely on this entry, not on #1.

### 2026-05-18 — Correction: head-to-head IS the right framing; bias is real

**What I got wrong.** Last turn I speculated that the cold stratosphere might not be a real bias — that with O₃=0 the pure-radiative answer just *is* cold. That speculation was untethered from the actual data. `scripts/experiments/radiative_equilibrium_pf_convergence.py` already runs RRTMG alongside the PF variants under identical conditions (CO₂=376 ppm, q=0, no O₃, SW=240 W/m²) for 1000 days. The script's per-band diagnostic block excludes RRTMG (gated on `PF_TABLE[label] is not None`), so the printed `T_local` numbers in the experiment-#1 log were PF-only and I didn't pull RRTMG's T_strat from the history.

**What the head-to-head actually shows** (read from `radiative_equilibrium_pf_convergence.png`, same run that produced the experiment-#1 numbers):

| level | RRTMG | PF (any variant) | gap |
|---|---:|---:|---:|
| ~100 hPa | ~210–220 K | ~155–170 K | 40–60 K |
| 0.1–1 hPa | ~190–200 K | ~120–140 K | 60–70 K |

The LW heating rate panel shows the mechanism: RRTMG's stratospheric LW heating rate is ≈0 (radiative equilibrium / Schwarzschild balance), while every PF variant has substantial net cooling extending well into the upper atmosphere.

**Implications for prior diagnoses.**
- Experiments #1 and #2 (Chaverot wing/core inspection, Chaverot vs LBL) were looking at the *source k-table* and found it is faithful to LBL. That conclusion stands: **the k-table is not the bug.**
- The "cold-stratosphere bias is upstream of the integrator" claim from the 2026-05-18 ngpt-rejection entry was **wrong**. The bias is real and head-to-head with RRTMG under identical inputs, so it is *downstream* of the source CDF — in the band partition, the integrator, or in physics RRTMG includes that PF does not.

**Open hypotheses (revised).**
1. **PF integrator / two-stream divergence from RRTMG's source-function treatment.** RRTMG fits the in-band Planck source with a specific scheme (multiple source-function moments per layer). PF may use a single-Planck-fraction-per-band approximation that overestimates cooling-to-space when the band is very optically thick (the CO₂ core).
2. **CO₂ hot-band or 4.3 µm contribution dilution.** PF's band containing 2326 cm⁻¹ averages it with surrounding window, possibly weakening absorption that RRTMG resolves.
3. **Missing absorbers RRTMG silently includes even with q=0 and O₃=0.** Most likely candidates: N₂-N₂ CIA (rotational region), N₂-O₂ CIA, or some baseline CO₂ continuum / line-mixing treatment that RRTMG uses by default and the Chaverot table omits.

**Next concrete steps (in order).**
1. Add a print of `T_local` and column-min T for RRTMG to `radiative_equilibrium_pf_convergence.py`'s diagnostic block, plus a per-level RRTMG–PF temperature *and* heating-rate difference printout, to put exact numbers on the gap (no re-run needed if a `.npz` was saved; otherwise one re-run).
2. Inspect what gases/continua RRTMG's `RRTMGLongwave` component activates by default with our state inputs (climt wrapper level; or RRTMG fortran source comments). Specifically check whether N₂ CIA is on.
3. If RRTMG has a CIA/continuum source that Chaverot doesn't, that's the bias. If not, the bias is in band partition × integrator and needs an integrator-level diff.

### 2026-05-18 — Follow-up experiment #3: Quantify the RRTMG–PF gap per level

**Why.** The prior convergence run set up RRTMG and PF head-to-head under identical conditions but the diagnostic block printed `T_local` only for PF (gated on `PF_TABLE[label] is not None`). The previous plot makes the gap visible at eyeball precision (RRTMG ~210 K at 100 hPa, PF ~155–170 K) but we don't have exact numbers per level, per variant, and we don't have heating-rate diffs that would tell us *where in the column* the gap is generated.

**Method.** Modify `scripts/experiments/radiative_equilibrium_pf_convergence.py`:
- Add `--days N` argparse (default unchanged at 1000) so we can re-run on a fast budget.
- Save the final history as `.npz` so future analysis doesn't require re-running.
- Add a new diagnostic block AFTER the TOA balance print: for each PF label, print a per-level table of `(p, T_RRTMG, T_PF, ΔT, HR_RRTMG, HR_PF, ΔHR)` so we can localize where in the column PF cools more than RRTMG and by how much.

Then re-run at ~300 days (a compromise: the previous plot's top-3-level T(t) panel showed RRTMG levels off by day 200–300, while PF was still drifting in the upper levels through day 400–600 — so 300 days catches the bulk of the stratospheric pattern without the full 1000-day cost).

**Expected outputs.**
- Column-aligned table of ΔT(p) and ΔHR(p) per PF variant.
- Identification of the pressure layer where the heating-rate gap is concentrated. If the gap is broadband across the stratosphere, it points at integrator differences or a missing continuum. If concentrated in a narrow layer, it points at a specific band's saturation behavior.

### 2026-05-18 — Follow-up experiment #4: PF k vs LBL on dry-adiabat profile

**Why.** Experiment #2 spot-checked k-table fidelity at one (T, p) point (T=250 K, p=100 hPa, three bands) and found ~factor-of-2 agreement. Experiment #3 quantified the RRTMG-vs-PF bias as a column-internal redistribution: PF +4 to +7 K below 400 hPa, −25 to −27 K at 23–41 hPa, OLR matching at 241–243 W/m². To distinguish "k-table is fine, integrator is the bug" from "k-table itself diverges across the column," need a per-layer, per-band comparison along a realistic T(p) profile.

**Method.** New script `scripts/experiments/compare_pf_k_vs_lbl_profile.py`:
- LBL: `lbl_k_adiabat.nc` (dry adiabat T_surf=288 K, 60 levels, dnu=0.1, CO₂-only at 376 ppm in air, pseudovoigt).
- PF: load `k_coefficients(1, nband, ngpt, nT, nP, nX)` from each candidate PF table, g-integrate to band-mean k(b, T, P, X), bilinear-interp in (log p, T) at X_H2O=1e-6 (table minimum) to each LBL layer's (p, T).
- LBL band-mean κ over each PF band's wavenumber window (uniform-dnu mean).
- Compare CO₂-dominated bands only (LBL is CO₂-only — non-CO₂ bands like 10–500, 800–1800 give spurious 1e6+ ratios because PF has H₂O baseline and LBL doesn't; expected, not real).
- Restrict ratio summary to T ≥ 100 K to avoid double-sided extrapolation artifacts (PF T-grid bottoms at 50 K; LBL `T_FLOOR=80 K`).
- Spot-print at p ∈ {1000, 500, 200, 100, 50, 20} hPa with T from the adiabat.

**Result (PF/LBL ratio, `earth_low_res_lw_co2refined_gl` table).**

| p [hPa] | T [K] | 500–630 (wing) | 630–700 (core) | 700–800 (wing) | 1800–3250 |
|---:|---:|---:|---:|---:|---:|
| 1013 | 288.0 | **1.99** | **1.84** | 1.47 | 0.93 |
| 464  | 230.4 | 1.15 | 1.07 | 0.98 | 0.40 |
| 212  | 184.3 | 1.96 | 1.10 | 1.11 | 0.29 |
| 97   | 147.4 | **4.15** | 0.98 | 1.39 | 0.16 |
| 52   | 123.3 | **6.67** | 0.61 | 1.78 | 0.14 |
| 20   |  94.3 | **19.74** | 0.26 | 1.69 | 0.10 |

Same patterns for the 10-band and 4-band shipped tables (different band groupings, same trend within shared CO₂-dominated bands).

**What we learned.**
1. **Two systematic k-table biases align with the experiment-#3 RRTMG-vs-PF temperature pattern:**
   - At warm T (surface): PF is ~2× too opaque across all CO₂ bands. More absorption near surface → traps energy in the boundary layer → PF surface is warmer (matches +4 to +7 K tropospheric warm bias).
   - At cold T (strat): 500–630 cm⁻¹ wing PF/LBL grows from 1.15 at 230 K to **19.7 at 94 K**. More wing opacity at cold T → cool-to-space radiates from higher/colder layers → cold strat (matches −25 K bias at 23–41 hPa).
2. **The core 630–700 cm⁻¹ band agrees well (ratio 0.98–1.10) at T ≥ 147 K**, then under-shoots at colder layers (0.26 at 94 K). This partly offsets the wing overestimate but doesn't cancel it.
3. **1800–3250 cm⁻¹ (CO₂ v3) is uniformly under-opaque in PF (ratio 0.10–0.93)**, but Planck emission there is small at strat T, so radiative impact is limited.
4. **The integrator is no longer the leading suspect.** With these k-table errors of factor 2–20 at strat T, even a flawless integrator would produce a substantial cold bias. Whatever residual integrator effect exists is second-order.

**Mechanism candidates (in order of plausibility, none yet tested).**
1. `bin_down` distortion. Averaging native R=500 sub-bands into PF bands may bias wing-vs-core differently at low T.
2. GL rebin from 20 native g-points onto 8 PF g-points. GL nodes cluster toward (0, 1); at cold T where the strong-line tail dominates band-mean k, sparse sampling of the weak end may inflate the band-mean estimate.
3. Coarse T-grid linear interpolation between widely-spaced nodes (50, 110, 170, …). At T=147 K we interp between 110 and 170 K; line strengths can shift substantially over 40 K — linear-k may fail where log-k or T-aware interp would not.
4. Chaverot source CDF artifact at cold T. Experiment #2 confirmed source fidelity at T=250/p=100; never spot-checked at T < 200 K.

**Next experiment (proposed, not yet run): raw Chaverot R=500 vs LBL at (T=147 K, p=97 hPa) — the strat-bias hotspot.** Extends experiment #2 script to one more (T, p) point.
- If raw Chaverot agrees with LBL at the cold strat layer → bias is in our `bin_down` + rebin + T-interp pipeline (PF-side fix; mechanism candidates 1–3).
- If raw Chaverot also disagrees → bias is in the source CDF at low T (mechanism candidate 4; need a higher-T-resolution source or to regenerate Chaverot at finer T grid in the strat range).

**Files touched.** New: `scripts/experiments/compare_pf_k_vs_lbl_profile.py`, `debug_data/pf_vs_lbl_k_profile__earth_low_res_lw_co2refined_gl.png`, `..._10band_gl.png`, `..._8gpt_twostretch.png`. Modified: `scripts/experiments/radiative_equilibrium_pf_convergence.py` (added `--days` argparse, per-level RRTMG-vs-PF diagnostic block, `.npz` snapshot save). No code in `climt/` changed.

**Resume here after restart.** First action: write the (T=147 K, p=97 hPa) spot-check by extending `compare_chaverot_vs_lbl_wings.py` to add one (p, T) row, then decide between the pipeline-fix and source-fix paths.

---

### Experiment #5 — context update + T-grid resolution hypothesis (2026-05-29)

**New external evidence the bias is real.** Colleague running SOCRATES (`../other_libs/socrates_2511`) reports its radiative-equilibrium profile is *similar to RRTMG*. Two independent reference codes (RRTMG, SOCRATES) now agree against PF → the cold-strat bias is a PF defect, not "reality is cold." This rules out the possibility that RRTMG itself was the outlier.

**Question to test now:** is the bias "simply due to the low resolution of the tables we generate"? Three resolution axes exist: (a) number of bands, (b) number of g-points, (c) the **T-grid**. Axes (a) and (b) were already explored in experiments #1–4 (finer band splits and the 11-band tightcore table made things *worse*, not better; g-point count is not the driver). The **un-rejected axis is the T-grid.**

**Chaverot source T-grid (inspected via h5py):** `[50, 110, 170, 230, 290, 350, 410, 470, 530, 590, 650, 1000]` K. The strat layer (T≈147 K) is interpolated **linearly between the 110 K and 170 K nodes** — a 60 K gap, exactly where experiment #4 found the wing PF/LBL ratio exploding (4.15 at 147 K → 19.7 at 94 K). p-grid is decade-spaced `[1,10,100,…,1e7]` Pa.

**Generator cannot add T resolution.** `scripts/generate_picket_fence_tables.py:231` sets `"T_grid": np.asarray(k_binned.tgrid)` — the climt table inherits the Chaverot source T-grid verbatim. So if coarse T is the cause, the fix requires either (i) a finer-T source, or (ii) better interpolation in the loader.

**Loader uses LINEAR-in-k T interpolation.** `climt/_components/picket_fence/optics/correlated_k.py:96` `interpolate_k`: trilinear with linear-in-k in both T and P (`v000=base[iT,iP,iX]` …), log-linear only in the H2O (X) axis. CO₂ line strengths vary strongly and non-linearly over a 60 K gap at strat T; linear-in-k between 110 K and 170 K nodes may badly overestimate k(147 K). **If log-k interpolation is much closer to truth, the fix is cheap (loader-only, no table regen).**

**Experiment T1 (cheap, decisive — running now):** LBL κ(T) at fixed p≈100 hPa across a fine cold-T grid, band-mean over the diagnostic CO₂ bands, then compare `linear_interp[k(110),k(170)]` at 147 K and log-k interp vs true `k(147)`. Quantifies how much of the wing overestimate is *interpolation error* (loader-fixable) vs *intrinsic source error* (needs finer source).

**T1 RESULTS (2026-05-29).** Script: `scripts/experiments/lbl_kT_curvature.py` (linepyline env); data: `debug_data/lbl_kT_curvature.npz`. At p=1e4 Pa (~100 hPa), 376 ppm CO₂, target T=147 K interpolated between Chaverot nodes 110/170 K (w=0.617):

| band (cm⁻¹) | k_true(147) | linear_interp | log_interp | lin/true | log/true |
|---|---:|---:|---:|---:|---:|
| 500–630 (wing) | 3.41e-5 | 5.18e-5 | 2.56e-5 | **1.52** | 0.75 |
| 630–700 (core) | 1.14e-1 | 1.14e-1 | 1.14e-1 | 1.00 | 1.00 |
| 700–800 (wing) | 8.34e-5 | 1.21e-4 | 6.39e-5 | **1.45** | 0.77 |
| 1800–3250 (v3) | 4.65e-2 | 4.65e-2 | 4.65e-2 | 1.00 | 1.00 |

Full LBL k(T) wing 500–630: 9.6e-7 (94 K) → 4.0e-6 (110) → 9.7e-6 (123) → 3.4e-5 (147) → 8.2e-5 (170) → 1.2e-4 (184) → 3.5e-4 (230). **The wing k grows ~exponentially with T — a factor of ~20 across the single 110→170 K Chaverot interval (roughly doubling every ~13 K).**

**What T1 establishes.**
1. **Linear-in-k interpolation (what the loader does) overestimates wing k at the strat by ~1.5×**, purely from T-grid curvature — the 60 K-wide 110↔170 K bracket cannot follow the exponential rise. This is a *real* contributor to the cold-strat bias and is **directionally consistent**: too much wing opacity → cool-to-space from higher/colder layers → cold strat.
2. **The core (630–700) and v3 (1800–3250) bands are essentially flat in T and interpolate perfectly** (ratio ≈ 1.00). The bias lives entirely in the optically-thin wings, exactly as experiment #4 saw.
3. **But interpolation curvature alone (1.5×) does not explain experiment #4's 4.15× PF-table/LBL wing ratio at 147 K.** The remaining ~2.7× is *not* interpolation — it is intrinsic to the PF table at the nodes (bin_down + GL rebin, or the Chaverot source CDF at cold T). So "low resolution" is *one* contributor, not the whole story.
4. **log-k interpolation undershoots (0.75×) rather than overshooting** — directionally it would *help* the cold-strat bias (less wing opacity → warmer strat), and its geometric error (|ln 0.75|=0.29) is a bit smaller than linear's (|ln 1.52|=0.42). A loader switch to log-k is a cheap, directionally-correct partial fix, but neither scheme nails a curve this steep on a 60 K bracket — **the proper fix is a finer T-grid in the strat range (a node near ~140 K).**

**Decision for next step.** Two independent fixes are now on the table, in increasing cost:
- (cheap, loader-only, no regen) switch T-interpolation from linear-in-k to log-in-k in `correlated_k.py:interpolate_k`. Test impact on RCE strat T. Expected: partial warming.
- (expensive, regen) build a table with a finer strat T-grid. Generator inherits the Chaverot grid, so this needs either a finer-T Chaverot source or a linepyline-built Earth table on a custom T-grid. This also addresses the residual ~2.7× node error if it comes from rebin.

**T2 (next):** isolate node-fidelity from interp error — compare the *shipped PF table* band-mean k at the 110/170 K nodes against LBL at those same nodes, to confirm the ~2.7× residual is node/pipeline error rather than interpolation. Then decide regen vs loader-fix.

**T2 RESULTS (2026-05-29) — the bias is NOT the T-grid; it is baked into the cold table node.** Script: `scripts/experiments/pf_node_fidelity_vs_lbl.py` (radiation env). Compared shipped `earth_low_res_lw_co2refined_gl` band-mean k *directly at its T-nodes* 110/170 K, at the exact p-node 1e4 Pa (no p-interp, no T-interp), vs LBL. p=1e4 Pa and T=110/170 K are all exact Chaverot nodes; the table bands match the diagnostic windows exactly.

| band (cm⁻¹) | PF@110 | LBL@110 | **node110** | PF@170 | LBL@170 | node170 | PFlin(147) | LBL(147) | r147 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 500–630 (wing) | 5.77e-5 | 3.96e-6 | **14.55** | 1.99e-4 | 8.16e-5 | 2.44 | 1.45e-4 | 3.41e-5 | **4.24** |
| 630–700 (core) | 9.93e-2 | 1.10e-1 | 0.90 | 1.23e-1 | 1.16e-1 | 1.06 | 1.14e-1 | 1.14e-1 | 1.00 |
| 700–800 (wing) | 7.37e-6 | 1.12e-5 | 0.66 | 1.85e-4 | 1.89e-4 | 0.98 | 1.17e-4 | 8.34e-5 | 1.40 |
| 1800–3250 (v3) | 7.79e-3 | 4.74e-2 | 0.16 | 7.51e-3 | 4.60e-2 | 0.16 | 7.62e-3 | 4.65e-2 | 0.16 |

**Conclusions.**
1. **The cold-strat bias is a TABLE NODE error, not interpolation.** The 500–630 wing PF k is **14.5× too opaque at the 110 K node** and 2.4× at 170 K — *before any interpolation runs*. The r147=4.24 (reproducing experiment #4's 4.15) is just linear interp landing between an over-opaque cold node and a less-over-opaque warm node. T1's 1.5× interp-curvature is a minor effect riding on top of a node that is already 2.4–14.5× wrong.
2. **This rejects both the user's "coarse T-grid" hypothesis AND the cheap log-k loader fix as primary remedies.** A finer T-grid or log-k interp cannot rescue a node that is 14.5× too opaque. The error is in *how the table value at the cold node is produced*.
3. **The culprit band is uniquely 500–630 cm⁻¹** (the far wing of the 15 µm band). The core (630–700) is faithful (0.90–1.06) and the other wing (700–800) is faithful at the nodes. v3 (1800–3250) is uniformly 6× *under*-opaque but radiatively minor at strat T.
4. **Most likely mechanism: the bin_down + GL 20→8 g-point rebin inflating the weak/far-wing band-mean at cold T**, or a Chaverot source-CDF artifact at 110 K. At 110 K the 500–630 far wing is genuinely near-transparent (LBL k=4e-6); GL nodes cluster toward g=0,1 and sparsely sample the weak end, so a few mis-weighted strong sub-bins can dominate the band-mean. This is a *resolution* problem — but in the **g-point/CDF rebin, not the T-grid** (the user's instinct that "low resolution" is to blame may be right, just the wrong axis).

**Decisive next experiment (#6): raw Chaverot R=500 k vs LBL at (T=110 K, p=1e4 Pa) for the 500–630 wing.** This was the originally-proposed spot-check; T2 now makes the target temperature unambiguous — go to the *coldest* node (110 K), where the error is largest (14.5×), not 147 K.
- If raw Chaverot (pre-rebin) agrees with LBL at 110 K → the error is in *our* `bin_down` + GL rebin pipeline → fix on the climt-generator side (better quadrature, finer g-points, or log-space rebin).
- If raw Chaverot also reads ~14× high at 110 K → the Chaverot source CDF is wrong in the far wing at cold T → need a finer-T or higher-fidelity source (or switch this band to a linepyline-built table).

**Files added.** `scripts/experiments/lbl_kT_curvature.py`, `scripts/experiments/pf_node_fidelity_vs_lbl.py`, `debug_data/lbl_kT_curvature.npz`. No `climt/` code changed.

---

### Experiment #6 — raw Chaverot R=500 source vs LBL at the cold node (DECISIVE) (2026-05-29)

Script: `scripts/experiments/raw_chaverot_vs_lbl_coldnode.py` (radiation env, exo_k 1.3.0). Loads the native R=500 source `Earth_376ppmCO2_R500_10-30k.ktable.SI.h5` (4004 bins × 20 g-points), computes the uniform-ν band-mean k over the diagnostic bands at the exact (T,p) nodes **bypassing bin_down and the GL 20→8 rebin**, applies the generator's unit conversion (×N_A/M_air = 2.079e25, m²/molecule→m²/kg, per `generate_picket_fence_tables.py:359`), and compares to LBL.

| band (cm⁻¹) | T | native R=500 (m²/kg) | LBL (m²/kg) | **native/LBL** |
|---|---:|---:|---:|---:|
| 500–630 (wing) | 110 | 2.165e-4 | 3.96e-6 | **54.6** |
| 630–700 (core) | 110 | 1.232e-1 | 1.105e-1 | 1.12 |
| 700–800 (wing) | 110 | 1.339e-5 | 1.116e-5 | 1.20 |
| 500–630 (wing) | 170 | 3.514e-4 | 8.16e-5 | **4.31** |
| 630–700 (core) | 170 | 1.427e-1 | 1.161e-1 | 1.23 |
| 700–800 (wing) | 170 | 2.780e-4 | 1.890e-4 | 1.47 |

**DECISIVE CONCLUSION — the cold-strat bias is a Chaverot SOURCE defect in the 500–630 cm⁻¹ band, not a resolution or pipeline problem.**
1. **The raw R=500 source 500–630 wing is 54.6× too opaque at 110 K** (4.3× at 170 K) vs line-by-line — *before any of our processing*. The core (630–700) and the other wing (700–800) are faithful at both nodes (1.1–1.5×). So the bias is isolated to one band, and it lives in the source CDF.
2. **Our bin_down + GL 20→8 rebin does NOT create the bias — it partially MASKS it.** Native 54.6× → shipped table 14.5× (T2): the GL-8 rebin under-samples the steeply-skewed CDF's high-k tail and *reduces* the band mean. So a "better" (denser-g) rebin would make the shipped table *more* opaque, i.e. worse — which is exactly why experiments #1–4 (finer bands/g-points) made things worse. **This retroactively explains the #1–4 paradox.**
3. **This rejects every resolution hypothesis** — T-grid (exp T2), band count and g-point count (exps #1–4, and the masking mechanism above), and the cheap log-k loader fix. None of them touch the source-CDF error. The user's "low resolution tables" hypothesis is **answered: no.** The Chaverot Earth LW source is simply wrong in the 15 µm far-wing at strat temperatures.
4. **Likely spectroscopic cause:** 500–630 cm⁻¹ is the far red wing of the CO₂ ν₂ band. linepyline (pseudovoigt) applies a sub-Lorentzian far-wing/χ-factor treatment; the Chaverot source LBL evidently uses fuller Lorentzian wings (no sub-Lorentzian cutoff, or different continuum), over-counting far-wing absorption — and the relative over-count explodes as lines narrow at cold T (4×→55× from 170→110 K).

**RECOMMENDED FIX (for user decision): build the Earth LW table from linepyline, not Chaverot.** We already have the full linepyline table pipeline (`scripts/generate_pf_tables_linepyline.py` + `pf_table_builder/`) used for Titan/TRAPPIST-1e, and linepyline is the LBL reference that exonerates the core and convicts the source. Adding an `"earth"` CO₂=376 ppm LW scenario to that registry and regenerating `earth_low_res_lw` would replace the defective far-wing with the correct sub-Lorentzian opacity. This is the clean path; a Chaverot-side χ-factor patch is not available to us.

**Files added (exp #6).** `scripts/experiments/raw_chaverot_vs_lbl_coldnode.py`. No `climt/` code changed.

#### Addendum — RRTMG_LW resolution cross-check (rules resolution out for good)

Read from `climt/_lib/rrtmg_lw/{parrrtm,rrtmg_lw_init}.f90`. RRTMG_LW = **16 bands, 140 reduced g-points** over 10–3250 cm⁻¹. Band edges `wavenum1/2` and reduced g-points `ngc`:

| RRTMG band | range cm⁻¹ | g-points (ngc) |
|---|---|---|
| 3 | **500–630** | **16** |
| 4 | 630–700 | 14 |
| 5 | 700–820 | 16 |

So RRTMG covers the culprit 500–630 region with a *single band and 16 g-points*. The **Chaverot R=500 source has 116 native bins in 500–630, each a 20-point k-CDF** — i.e. ~116×20 spectral/quadrature sampling, *far finer than RRTMG's 16 g-points*. Yet RRTMG (and SOCRATES) get the strat right and the Chaverot-derived PF table is 54× too opaque there.

**Conclusion: resolution is definitively not the cause.** The accurate code (RRTMG, built from LBLRTM) has *lower* 500–630 resolution than the inaccurate one. What RRTMG has that the Chaverot source lacks is correct far-wing **spectroscopy** — LBLRTM's sub-Lorentzian CO₂ line-shape χ-factor and the MT_CKD/CKD CO₂ continuum, which cut the far-wing absorption that a plain-Lorentzian LBL over-counts (worst at cold T, narrow lines). This is the same mechanism linepyline (pseudovoigt + continuum) captures and that exonerates the core/700–800 bands. The fix remains: regenerate the Earth LW table from a line-shape-correct LBL (linepyline), not add resolution.

---

### Experiment #7 — linepyline-built Earth LW table: sanity check FAILED, root cause = `binning=True` (2026-05-29)

Generated `earth_low_res_lw_co2refined_linepyline.nc` (6 LW bands, grids matching the shipped table) via `scripts/generate_pf_tables_linepyline.py --scenario earth --kind lw`, then ran the node-fidelity sanity check (`pf_node_fidelity_vs_lbl.py <table>`) vs LBL.

**Sanity check FAILED.** The new table read ~5× too opaque in *every* band — including the core (630–700), which the LBL reference says should be ~1.0×:

| band | OLD Chaverot node110 | NEW linepyline node110 | expected |
|---|---:|---:|---:|
| 500–630 (wing) | 14.55 | 4.86 | ~1 |
| 630–700 (core) | 0.90 | **4.83** | ~1 |
| 700–800 (wing) | 0.66 | 8.14 | ~1 |
| 1800–3250 (v3) | 0.16 | 1.84 | ~1 |

**Root cause pinned by differential probe** (`scripts/experiments/probe_lpl_kappa_apis.py`, linepyline env). Computed band-mean CO₂ kappa at (110/170 K, p=1e4 Pa, dry 376 ppm) several ways and compared to `radiative_transfer` (the LBL ground truth, which matches Chaverot's core at 0.110). Key result — the builder's `get_kappa_hitran × mass_fraction` matches LBL **exactly (1.00×)** when `binning=False`, but is wildly off with `binning=True`:

| get_kappa_hitran setting | core@110 vs LBL | 700–800@170 vs LBL |
|---|---:|---:|
| `binning=True, dnu=0.5` ← **what the table used** | **9.5×** | **19.2×** |
| `binning=True, dnu=0.1` | 2.0× | 4.0× |
| `binning=False, dnu=0.1` | **1.00×** | **1.00×** |
| `binning=False, dnu=0.5` | 1.3× | 0.72× |

Findings:
1. **`mass_fraction` was never the problem.** `radiative_transfer` internally does the identical weighting `kappa += kappas[name] * q`, where `q = (M_species/mean_mol_weight)·f` = 5.711e-4 for CO₂ — bit-identical to the builder's `_mass_fraction`. The two code paths differ *only* in the kwargs passed to `get_kappa_hitran`.
2. **`line_shape` is irrelevant for this band** (lorentz ≈ pseudovoigt to <0.1%).
3. **`binning=True` is the bug.** linepyline's line-width binning approximation overestimates kappa by 2–4× at dnu=0.1 and 9.5–19× at dnu=0.5, worst in the cold/low-pressure strat and growing with the diagnostic band's sharpness — exactly the regime its own docstring warns about ("reduces accuracy, especially at low pressure/narrow linewidths in the upper atmosphere. Use with caution!"). The coarse `dnu=0.5` compounds it.
4. **The validated ground-truth recipe is `binning=False` + fine `dnu=0.1`** — this reproduces `radiative_transfer`/LBL exactly across all bands and temperatures.

**Latent-bug flag:** the shipped Titan and TRAPPIST-1e linepyline tables were generated with the old `binning=True` default and are therefore likely 2–19× too opaque in their cold upper layers. They were never validated against absolute LBL k. Should be regenerated with `binning=False` too.

**Code fix applied.** `scripts/generate_pf_tables_linepyline.py`: defaults changed to `dnu=0.1`, `line_shape="pseudovoigt"`, `binning=False`; added `--binning` opt-in flag (documented as low-accuracy preview only). No change to `pf_table_builder/`.

**Next:** regenerate the Earth LW table with the corrected defaults (`binning=False`, `dnu=0.1`) and re-run the sanity check — expect the 500–630 wing node110 to collapse from 14.5× toward ~1× and the core to return to ~1×.

**Files added (exp #7).** `scripts/experiments/probe_lpl_kappa_apis.py`. Code edited: `scripts/generate_pf_tables_linepyline.py`.

#### Experiment #7 — final result: corrected table PASSES (2026-05-29)

Two bugs were masking each other in the first linepyline build:
1. **`binning=True`** → kappa 9.5–19× too opaque (fixed: `binning=False`, `dnu=0.1`).
2. **`ngpt=2`** (CLI default) → GL-2 quadrature under-integrates the skewed k-CDF band-mean, collapsing it ~16× *low*. The shipped Chaverot table uses **`ngpt=8`** (GL-8 weights); a drop-in replacement must match. Fixed by regenerating with `--ngpt 8`.

Final node-fidelity check (`binning=False`, `dnu=0.1`, `ngpt=8`), band-mean k via GL-8 g-integration vs LBL:

| band (cm⁻¹) | OLD Chaverot @110 | NEW linepyline @110 | @170 | @147 (interp) |
|---|---:|---:|---:|---:|
| **500–630 (wing, the defect)** | **14.55** | **0.42** | 0.38 | 0.57 |
| 630–700 (core) | 0.90 | 0.87 | 0.72 | 0.78 |
| 700–800 (wing) | 0.66 | 0.55 | 0.63 | 0.91 |
| 1800–3250 (v3) | 0.16 | 0.08 | 0.08 | 0.08 |

**PASS.** The 500–630 over-opacity (14.5×, the sole cause of the cold-strat bias per exp #6) is eliminated → 0.42×. The core is faithful (0.87, matching the trusted Chaverot 0.90). No band is anomalously over-opaque.

The band-means sitting at 0.4–0.9× (not 1.0×) is the **GL-8 quadrature signature, not a kappa error**: the probe proved per-ν kappa is exactly 1.00× vs LBL, and the shipped Chaverot core shows the same 0.90× under-read. GL-8 under-samples the high-k line-center tail (g>0.98) that barely affects flux; the low-k g-points that dominate strat cooling-to-space are well sampled. So the table is radiatively sound. (A full RCE strat-T comparison vs RRTMG would be the definitive end-to-end confirmation; deferred.)

**Deliverable:** `climt/_data/picket_fence/correlated_k/earth_low_res_lw_co2refined_linepyline.nc` — validated drop-in for the defective Chaverot `earth_low_res_lw_co2refined_gl.nc`, identical grids/layout (ngpt=8), correct sub-Lorentzian far-wing opacity.

**Open follow-ups (not yet done):** (a) wire the linepyline table as the production Earth LW default if an RCE run confirms the strat warms to RRTMG levels; (b) build the matching Earth SW table; (c) regenerate Titan & TRAPPIST-1e linepyline tables with `binning=False` (they inherit the 2–19× binning over-opacity bug).

#### Experiment #8 — ngpt convergence of the GL band-mean (confirms exp #7's "metric artifact" claim) (2026-05-29)

Script: `scripts/experiments/ngpt_convergence.py` (linepyline env). Samples CO2 kappa once at the diagnostic node (T=110/170 K, p=1e4 Pa, dry 376 ppm) with the validated recipe (binning=False, dnu=0.1, pseudovoigt), then runs the REAL pipeline `kappa_to_k_coeffs` at ngpt = 2,4,8,16,32,64,128 and computes the g-integrated band-mean `sum(k_i*w_i)`. Convergence target = arithmetic mean of the sampled kappa (= LBL band-mean; the sampled arith_mean matched LBL to 4 digits at every band).

Band-mean / LBL ratio vs ngpt (T=110 K; 170 K nearly identical):

| band | g=2 | g=4 | g=8 | g=16 | g=32 | g=64 | g=128 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 500–630 | 0.04 | 0.13 | 0.41 | 0.74 | 1.03 | 0.98 | 0.99 |
| 630–700 (core) | 0.06 | 0.20 | 0.87 | 0.89 | 1.00 | 1.02 | 1.00 |
| 700–800 | 0.07 | 0.28 | 0.60 | 0.95 | 0.95 | 0.99 | 1.01 |
| 1800–3250 (v3) | 0.00 | 0.00 | 0.08 | 0.23 | 0.62 | 0.95 | 0.98 |

**CONFIRMED.** The band-mean climbs monotonically to 1.00× LBL as ngpt grows (within a few % by ngpt=32; the sharp v3 band needs ngpt≥64). So the 0.4–0.9× deficit seen at the shipped ngpt=8 (exp #7) is **pure Gauss-Legendre quadrature truncation, not a kappa error** — the sampled kappa is exactly LBL at every ngpt.

Two caveats this clarifies:
1. The **shipped Chaverot table is also ngpt=8** and suffers the identical under-read (its core sits at 0.90× for the same reason). Both tables are on equal footing at ngpt=8; the cold-strat fix (500–630: 14.5× → 0.42×) is a *kappa* fix, robust at any ngpt.
2. **Band-mean convergence ≠ flux convergence.** The under-sampled tail is the saturated line centers, which contribute little to flux; ngpt=8 may give accurate radiation despite the band-mean under-read. Only an RCE/flux test (follow-up (a)) can settle the production ngpt.

**Files added (exp #8).** `scripts/experiments/ngpt_convergence.py`. No `climt/` code changed.

---

### Follow-ups deferred to next session (2026-05-29)

The Earth LW linepyline table is built and validated (exps #7–#8). The following are explicitly deferred — pick these up next time:

1. **End-to-end RCE confirmation.** Run a radiative-equilibrium column with the new `earth_low_res_lw_co2refined_linepyline.nc` and confirm the PF stratosphere warms to RRTMG/SOCRATES levels (closing the original 40–70 K cold bias). This is the definitive test; the node-fidelity + convergence checks are necessary-but-not-sufficient proxies. Also use it to decide the production ngpt (band-mean convergence ≠ flux convergence — see exp #8).
2. **Wire as production Earth LW default** (only after (1) passes) — swap the loader/scenario default from the Chaverot `*_gl.nc` to the linepyline table.
3. **Build the matching Earth SW table** (`--scenario earth --kind sw`) with the corrected `binning=False` defaults, for a consistent LW+SW pair.
4. **Regenerate Titan & TRAPPIST-1e linepyline tables with `binning=False`.** They were built with the old `binning=True` default and inherit the 2–19× cold-layer over-opacity bug diagnosed in exp #7; they were never validated against absolute LBL k.

---

### Experiment #9 — (T, p) envelope of Chaverot R=500 vs LBL, for the Chaverot-group report (2026-05-29)

Built a (T, p) envelope comparison + figures to substantiate the source-defect finding for an outreach message to the Chaverot et al. group. Scripts (all in `scripts/experiments/`): `chaverot_vs_lbl_envelope_native.py` (radiation env, exo_k — native R=500 band-means + per-bin spectra at the table's own T/p nodes), `chaverot_vs_lbl_envelope_lbl.py` (linepyline env — LBL band-means + fine kappa(nu)), `plot_chaverot_vs_lbl_envelope.py` (matplotlib — heatmap + spectra). Outputs: `debug_data/chaverot_vs_lbl_{native,lbl}.npz`, `debug_data/chaverot_vs_lbl_{ratio,spectra}.png`.

**Grid-convergence caveat (important).** Initially sampled p = {100, 1e3, 1e4, 1e5} Pa, but a dnu-convergence spot check (dnu 0.1 → 0.005) showed the LBL band-means are only grid-converged at **p = 1e4 and 1e5 Pa**. At p ≤ 1e3 the CO2 lines are Doppler-narrow (~1e-3 cm⁻¹ at 110 K) and dnu=0.1 under-samples them (700–800 swings 20× at p=100), so those cells reflect OUR reference's grid limit, not source defects. The figures were therefore restricted to p = 1e4 (stratosphere; also the regime cross-validated against RRTMG/SOCRATES in exps #6 & addendum) and p = 1e5 (troposphere control).

**Converged result — native/LBL band-mean ratio (>1 = source too opaque):**

| band | p (Pa) | T=110 | T=170 | T=230 | T=290 |
|---|---|---:|---:|---:|---:|
| **500–630** | 1e4 | **54.7** | **4.3** | **2.5** | 1.1 |
| 500–630 | 1e5 | 2.3 | 1.2 | 1.1 | 2.3 |
| 630–700 (core) | 1e4 | 1.1 | 1.2 | 1.8 | 0.8 |
| 630–700 (core) | 1e5 | 1.0 | 1.1 | 1.1 | 1.8 |
| 700–800 | 1e4 | 1.2 | 1.5 | 1.9 | 1.3 |
| 700–800 | 1e5 | 0.6 | 1.0 | 1.0 | 1.5 |

**Findings for the report.** (1) The 500–630 cm⁻¹ far-wing over-opacity is a smooth, strongly T-dependent effect at strat pressure (54.7× → 1.1× from 110 → 290 K) — a clean far-wing line-shape signature, not a normalization error (core is faithful). (2) Secondary observation: the native 500–630 band-mean is non-monotonic in pressure at 110 K (2.17e-4 at 1e4 Pa vs 9.8e-6 at 1e3 and 1.1e-5 at 1e5), suggesting a possibly-localized table issue at that node on top of the broad trend — flagged as a question, not a hard claim. (3) Outreach message drafted (in conversation) with the two figures + `raw_chaverot_vs_lbl_coldnode.py` / envelope scripts offered as reproduction material.

**Files added (exp #9).** `chaverot_vs_lbl_envelope_native.py`, `chaverot_vs_lbl_envelope_lbl.py`, `plot_chaverot_vs_lbl_envelope.py`. No `climt/` code changed.

---

### Experiment #10 — RCE end-to-end test: wing fix is radiatively inert; "cold-strat bias" was mostly missing O₃ (2026-06-01)

The deferred follow-up (a) — the definitive end-to-end RCE confirmation of the linepyline table. Ran `scripts/experiments/radiative_equilibrium_pf_convergence.py` trimmed to a 4-way head-to-head: RRTMG (16 bands) vs PF 4-band shipped vs PF 6-band Chaverot (`co2refined_gl`, defective) vs PF 6-band linepyline (`co2refined_linepyline`, the exp #7 fix). Dry, CO₂=376 ppm, slab surface, SW=240 W/m², 1000 days.

**Finding 1 — linepyline ≈ Chaverot; the wing fix is radiatively inert for strat T.** Despite the 500–630 cm⁻¹ band-mean k differing **35× at the 110 K node** (file-level and in-run confirmed; the two tables loaded as distinct data), the equilibrium strat T differs by **≤2.9 K** between the two 6-band tables — and *colder* in the mid-strat (wrong sign vs exp #4's prediction). Root cause: strat T is set by the **630–700 cm⁻¹ core** (faithful, ~0.95× in both tables), not the 500–630 wing. The wing reaches τ=1 in the **mid-troposphere** (~469 hPa, per the exp opacity-structure table), so a 35× opacity change in a band that radiates from the warm lower atmosphere cannot move the stratospheric radiating level. **This falsifies exp #6's central claim that the 500–630 wing over-opacity *causes* the cold-strat bias.** The wing defect is real spectroscopy (and the linepyline table is a legitimate k-fidelity improvement) but it is NOT the driver of strat T.

**Finding 2 — the "cold-strat bias" was largely a missing-O₃ comparison artifact.** RRTMG's `get_default_state` populates a **climatological ozone profile (max ~3.9 ppm in the stratosphere)**; the script set CO₂ and zeroed `specific_humidity` but **never zeroed O₃**, while the PF CO₂(+H₂O) tables contain none. So every prior "head-to-head" compared **RRTMG-with-O₃ vs PF-without**. The 2026-05-18 "Correction" entry's assertion that the comparison ran with "no O₃" was never enforced in code and is **wrong**. Zeroing `mole_fraction_of_ozone_in_air` in the RRTMG state cooled RRTMG's middle atmosphere by **~30 K at 40–90 hPa** (90 hPa: 178.5→148.7 K; 40 hPa: 170.5→141.4 K; 10 hPa: 153.4→136.9 K; 2.4 hPa: 143.6→135.0 K). PF curves unchanged. The −26 K PF gap at 90 hPa collapsed to **+3 K**.

**Residual biases under matched inputs (dry, CO₂-only, no O₃).** Two separate, smaller effects remain:
- **Troposphere → ~50 hPa: PF is WARMER than RRTMG by +5 to +16 K** (peak ~+16 K at 600–700 hPa). Shared by all PF tables (linepyline did not change it). Consistent with exp #4's warm-T ~2× CO₂ over-opacity — never validated or fixed at warm nodes (230–290 K).
- **Above ~40 hPa: PF still COLDER by −7 K (23 hPa) → −16 K (model top, 2.4 hPa).** Candidate causes: PF integrator / two-stream source-function treatment vs RRTMG, or core-saturation behavior.

**Rejected.**
- "500–630 wing over-opacity causes the cold strat" (exp #6) — **falsified end-to-end.**
- "PF has a 20–70 K cold-strat bias vs RRTMG" — **mostly an O₃ artifact**; under matched no-O₃ inputs the 40–90 hPa strat gap is within ±3 K. The headline number that motivated exps #1–#9 was a comparison-setup error.

**Implication for the deliverable.** The linepyline table stands as a correct far-wing k-fidelity improvement, but its justification must be **reframed**: it is NOT a cold-stratosphere fix. Shipping it (deferred follow-up (b)) should be decided on k-fidelity grounds, not strat-T grounds.

**Files touched.** `scripts/experiments/radiative_equilibrium_pf_convergence.py` — trimmed to the 4-way comparison, added the linepyline table, and zeroed RRTMG O₃. No `climt/` code changed.

**Next.**
1. Decide whether the tropospheric warm bias is a real warm-T CO₂ over-opacity — validate linepyline/Chaverot core k vs LBL at warm nodes (230/290 K), the regime exp #7 never checked.
2. Investigate the upper-strat (>40 hPa) residual cold bias — integrator/source-function diff vs RRTMG, the hypothesis line the wing rabbit-hole displaced (see 2026-05-18 "Correction" open hypotheses).
3. Optionally bound the O₃ contribution precisely (e.g. re-run with O₃ matched on both sides, or as an RRTMG-only sensitivity).

---

### Experiment #11 — Warm-node k-fidelity: the tropospheric warm bias is NOT a CO₂ k-table over-opacity (2026-06-01)

Followed up exp #10's residual tropospheric warm bias (PF +5 to +16 K vs RRTMG, surface → ~50 hPa). exp #4 blamed a warm-T CO₂ over-opacity, but exps #7/#8 only validated k-fidelity at the COLD nodes (110/170 K). Script `scripts/experiments/pf_warm_node_fidelity.py` compares PF table band-mean k (Chaverot & linepyline, GL-8 g-integrated, read at exact T/p nodes, dry) vs LBL (exp #9 envelope, grid-converged at p=1e4/1e5 Pa) at the WARM nodes 230/290 K.

**PF/LBL core (630–700) band-mean:**

| node | Chaverot | linepyline |
|---|---:|---:|
| 290 K, 1000 hPa (surface) | **1.89** | **0.81** |
| 230 K, 1000 hPa | 0.89 | 0.83 |
| 290 K, 100 hPa | 0.67 | 0.60 |

**Findings.**
1. **linepyline is more faithful at warm T too.** Core 0.81× at the surface node vs Chaverot's 1.89×. exp #4's warm-T ~2× over-opacity is a Chaverot defect; linepyline largely fixes it. So linepyline is a better k-table across the full T range (cold far-wing AND warm core), not just the cold wing — strengthening its case as the k-fidelity deliverable.
2. **But this is radiatively inert in the lower troposphere.** The 630–700 core is optically thick at the surface in both tables (cumulative τ TOA→sfc is order 10²–10³; both ≫ 1 → saturated). A 2.3× band-mean difference in an already-saturated band cannot change surface fluxes — which is exactly why exp #10's RCE troposphere was within ~1 K between the two tables despite this k difference.
3. **Therefore the +15 K warm bias is NOT a CO₂ k-fidelity problem.** It is shared by all PF tables and is a PF-vs-RRTMG difference in the optically-thin bands / two-stream integrator / absorber set (RRTMG carries O₂=0.21 and possible CO₂ line-mixing/continuum the Chaverot/linepyline tables omit).
4. **Context — it is a PURE-radiative-equilibrium artifact.** `radiative_equilibrium_pf_convergence.py` has NO convective adjustment, so its troposphere is superadiabatic and unphysically warm for *both* models. The project goal is radiative-*convective* equilibrium with a reasonable tropopause; with convective adjustment the troposphere is reset to the adiabat and the pure-radiative warm overshoot should largely collapse.

**Rejected.** "The tropospheric warm bias is caused by PF CO₂ k-table over-opacity" — the over-opacity is real in Chaverot but (a) fixed in linepyline and (b) radiatively inert because the core is saturated near the surface. Not the lever.

**Tooling bug fixed.** First pass mis-converted `pressure_grid_log` (it is natural-log of Pa, not log10), sampling the wrong p-nodes; corrected to `np.exp(lpg)`.

**Next.** Add convective adjustment (`DryConvectiveAdjustment`) to the RCE script and re-check the tropopause for PF vs RRTMG — the direct test of the original goal. The CO₂ k-table warm bias is not worth chasing further (saturated/inert).

**Separately flagged (feature gap, not this plan).** PF has no ozone absorber. exp #10 showed O₃ accounts for ~30 K of RRTMG's stratospheric warmth; for PF to match a real atmosphere it needs an O₃ 9.6 µm LW band/absorber. Tracked as a future enhancement.

**Files added.** `scripts/experiments/pf_warm_node_fidelity.py`.

---

### Experiment #12 — Ascend the hierarchy: dry & moist RCE scripts (PF-LW vs RRTMG-LW tropopause) (2026-06-01)

Pure radiative equilibrium (exps #10–#11) was the simplest rung — used to characterise PF-LW in isolation. Per user direction, moving up the hierarchy to radiative-*convective* equilibrium, the configuration whose tropopause was the original target. Two separate single-column scripts were built so the two convection treatments can be compared independently. Both reuse the established PF-vs-RRTMG setup: prescribed fixed SW (240 W/m², identical on all columns), O₃ zeroed everywhere (apples-to-apples LW; PF has no O₃ absorber yet — see the spawned task), columns = RRTMG-LW / PF-LW linepyline / PF-LW shipped 4-band, slab surface.

- **Dry:** `scripts/experiments/rce_dry_pf_vs_rrtmg.py` — `SimplePhysics(boundary_layer=False)` + `DryConvectiveAdjustment`, dry column (q=0). Tendency stepper = `AdamsBashforth([lw, slab])`; SimplePhysics and dry adjustment applied by operator-splitting.
- **Moist:** `scripts/experiments/rce_moist_pf_vs_rrtmg.py` — `SimplePhysics(boundary_layer=True)` + `EmanuelConvection` (moisture spins up from surface evaporation). Tendency stepper = `AdamsBashforth([convection, lw, slab])`; `convection.current_time_step` set each step.

Both mirror the canonical `examples/column_code_with_slab.py` / `gmd_radiative_convective.py` coupling. Both **smoke-tested (14 steps) and run end-to-end without error**; the coupling, PF band-count ordering (RRTMG built first), and diagnostics are validated. Full equilibration runs (dry ~400 d, moist ~200 d) are the user's to run. Each writes a T(p)/heating(/moisture) figure to `debug_data/` with a cold-point tropopause marker, plus an `.npz`. One plotting bug fixed: the Emanuel `air_temperature_tendency_from_convection` diagnostic has dims `('*','mid_levels')` (reversed vs the radiation profiles) — flattened for the single column.

**Next.** Run both to equilibrium; compare the PF-LW vs RRTMG-LW tropopause pressure/temperature under dry vs moist convection. Expectation from exp #11: the pure-radiative +15 K tropospheric warm bias should largely collapse once convection sets the lapse rate, leaving the tropopause as the discriminating diagnostic.

---

### Caveat (2026-06-01) — PF moisture handling is mechanically correct but UNVALIDATED

Traced via graphify (`component.py` → `optics/correlated_k.py`) and confirmed numerically. For the Earth `premixed_bg` tables, PicketFenceLongwave handles water vapour correctly *mechanically*:
- declares `specific_humidity` as a live input (`lw/component.py:103`);
- converts kg/kg → H₂O mole fraction with the exact relation `x = q/(q+(1−q)·M_H₂O/M_dry)` (`component.py:247`);
- interpolates the mixture k along the tabulated H₂O VMR axis — trilinear, linear-in-k in T, log-p, **log-linear in log-X**, clamped to `[1e-6, 1.0]` (`correlated_k.py:96`);
- forms τ = k(T,p,X) × air-column-mass (dimensionally consistent).

The H₂O axis genuinely carries water opacity — at 290 K/1000 hPa, dry→wet (X 1e-6→0.1) the rotation band (10–500) and window (800–1800) band-mean k both jump ~10⁵×; CO₂ bands ×1.5–3.6. So the moist RCE script exercises real water vapour, not a no-op. **All three Earth tables (linepyline, shipped, Chaverot-gl) have a populated H₂O axis.**

**But the H₂O opacity has never been validated against LBL.** Every k-fidelity experiment (#1–#11) was at the DRY node (X=1e-6) and CO₂ bands only. Open risks for the *moist* comparison specifically:
1. **H₂O VMR axis vs LBL — RESOLVED (see exp #13 below): faithful within the known GL-8 factor.**
2. **H₂O continuum (MT_CKD) — RESOLVED: it IS included.** `pf_table_builder/kappa_sampling.py:sample_kappa_grid` defaults `include_mtckd_continuum=True`, and for H₂O it adds `rtm.get_kappa_mtckd(...)` while passing `remove_plinth=True` to the line call (avoiding double-counting). `generate_pf_tables_linepyline.py:build_table` calls `sample_kappa_grid` *without overriding* that flag, so the Earth linepyline table was built **with** the MT_CKD self/foreign continuum on the H₂O axis. (Confirmed by code path, not yet by LBL — the magnitude/correctness still rides on the next-session H₂O validation.)
3. **Coarse H₂O grid** (7 nodes / 6 decades, log-linear) — same interpolation-curvature risk as the T-grid (exp T1), partly mitigated by log-spacing.

**Implication:** the DRY RCE script rests on validated (CO₂-only, dry) physics; the MOIST script's PF-vs-RRTMG differences — especially in the window region — are provisional until the H₂O axis is validated.

---

### Experiment #13 — H₂O-axis fidelity vs LBL (lines + MT_CKD): PF moisture is validated (2026-06-01)

Script `scripts/experiments/pf_h2o_axis_fidelity.py` (linepyline env). LBL reference = `pf_table_builder.kappa_sampling.sample_kappa_grid` with the validated recipe (binning=False, dnu=0.1, pseudovoigt), sampled **with the MT_CKD continuum ON and OFF**, at tropospheric nodes (T=230/290 K, p=1e4/1e5 Pa) across all 7 H₂O VMR nodes. Compared to the PF table GL-8 band-mean k. Reads via `np.exp(pressure_grid_log)` (ln-Pa). Output log: `debug_data/pf_h2o_axis_fidelity.log`.

**linepyline table, T=290 K / 1000 hPa, PF/LBL across X=1e-6→1e-1 (physical range) | continuum× at X=0.1:**

| band | PF/LBL (dry→X=0.1) | cont× @X=0.1 |
|---|---|---|
| 10–500 (H₂O rotation) | 0.76 → 0.83 | 1.02 |
| 500–630 | 0.69 → 0.83 | 1.28 |
| 630–700 (CO₂ core) | 0.81 → 0.88 | 1.16 |
| 700–800 | 0.93 → 0.95 | 2.39 |
| **800–1800 (window)** | **0.84 → 0.89** | 1.03 |
| 1800–3250 | 0.47 → 0.60 | 1.02 |

**Findings.**
1. **The H₂O axis is faithful to LBL within the known GL-8 quadrature under-read.** Across the physical humidity range (X ≤ 0.1) the radiatively important bands sit at PF/LBL ≈ 0.76–1.0 — the *same* flat Gauss-Legendre ngpt=8 truncation signature (exp #8) that affects the already-validated CO₂ dry node (core 0.81–0.88 here vs 0.87 dry in exp #7). No H₂O-specific defect.
2. **Ratios are flat across the X nodes** (window 0.84→0.89, rotation 0.76→0.83) → the X-axis storage/interpolation introduces no drift; the ~0.85 offset is the table-wide quadrature factor, not a moisture error.
3. **MT_CKD continuum is confirmed present AND correctly matched.** `cont×` (LBL_on/LBL_off) > 1 exactly where expected — window and CO₂ wings under humidity — and PF k tracks the continuum-**ON** reference, not OFF. This upgrades caveat #2 from "confirmed by code path" to **confirmed against LBL**. PF is NOT missing window absorption.
4. **Only-large deviations are benign:** 1800–3250 (≈0.47, the sharp CO₂-v3 band that exp #8 showed needs ngpt≥64; radiatively minor at terrestrial T) and the X=1.0 corner (unphysical 100% H₂O), both pure GL-8 truncation.

**Conclusion.** PF handles moisture **properly** — mechanically (correct q→VMR, trilinear interp, τ=k·air-mass) *and* in fidelity (H₂O lines + continuum match LBL to the same ~0.85 GL-8 factor as CO₂). The MOIST RCE script therefore rests on validated physics, with the caveat that the whole table carries a uniform ~0.85 ngpt=8 band-mean under-read (radiatively benign per exp #8: the under-sampled tail is saturated line cores that contribute little flux). Caveats #1 and #2 RESOLVED; #3 (interpolation *between* X nodes) untested but low-risk given the flat node ratios.

**Open (minor).** The builder `sample_kappa_grid` passes array p/T to `get_kappa_hitran`, which the installed linepyline only accepts for length≥2 arrays (length-1 hits a `float()` path that raises) — not a problem for real table builds (full p_grid) but worth knowing for single-point regen.

**Files added.** `scripts/experiments/pf_h2o_axis_fidelity.py`, `debug_data/pf_h2o_axis_fidelity.log`.

---

### Experiment #14 — Dry RCE bug: the "dry" run was secretly moist (2026-06-02)

First dry-RCE run (400 d, 5 m slab) came out **way too warm and unconverged**: RRTMG OLR=239.6/T_sfc=281 K (≈balanced), but PF-linepyline OLR=213 (imbalance +27 → still warming) at T_sfc=306 K, PF-shipped 228.6/301 K. PF surface 25 K above RRTMG *and climbing* — impossible for dry convective adjustment, which can only *cool* the surface below the pure-radiative value (pure-RE PF equilibrium was ~265–272 K in exp #10).

**Root cause (verified by minimal test).** `SimplePhysics(boundary_layer=False)` keeps `surface_fluxes=True`, so over `area_type='sea'` it **evaporates** moisture: a q=0 init rose to ~0.8 g/kg within ~10 steps. The "dry" run was actually moist → H₂O greenhouse → runaway surface warmth. Confirmed: printing `specific_humidity` showed q climbing from 0.

**Fix.** `SimplePhysics(boundary_layer=False, use_external_surface_specific_humidity=True)` + hold `surface_specific_humidity=0` (set at init and re-zeroed each step). Verified: q stays exactly 0.0 while the sensible-heat flux still operates (~1 W/m²). Also: mixed-layer depth 5→1 m (5× faster surface equilibration; equilibrium is slab-independent), and the report now prints TOA imbalance + q_max as convergence/dryness checks.

**Also fixed the tropopause diagnostic.** Cold-point returned the model top (2.40 hPa) for every column because with O₃=0 and no stratospheric heating the temperature decreases monotonically upward — there is no cold-point tropopause. Replaced with the **convective-layer top** (top of the constant-potential-temperature well-mixed layer), the physically meaningful tropopause in dry RCE.

**Note.** The MOIST script's evaporative moistening is *intended* (that's the point of the moist run) — no change there. Only the dry script needed the dryness fix.

**Status.** Dry script fixed and smoke-tested (q=0 confirmed); awaiting a full converged run to read the PF-vs-RRTMG dry tropopause/surface comparison.

**Why did RRTMG NOT run away under the same accidental moistening? (verified)** A side-by-side moist trajectory (both columns, identical SimplePhysics evaporation, 1 m slab, from 260 K) shows: **both moisten essentially identically** (q_max ≈ 1 g/kg for both), so moisture asymmetry is NOT the cause. But PF's OLR is *lower* than RRTMG's at every step *despite PF being warmer* (day 40: PF 252 W/m² @ 284 K vs RRTMG 262 W/m² @ 276 K) → PF has a stronger total greenhouse. The T_sfc gap widens monotonically (+2 K day 5 → +6 K day 20 → +9 K day 40 → +25 K day 400). So the divergence is driven by **PF's intrinsic over-trapping (the exp #11 dry warm bias: thin bands / integrator, NOT CO₂/H₂O k which are LBL-faithful)**, with the Clausius-Clapeyron water-vapour feedback *amplifying* it as PF warms. RRTMG, lacking the bias, finds a moist balance; PF (seeded ~+15 K warm) climbs past where the feedback converges. **Implication: the MOIST (Emanuel) run will likely show the same PF over-warming — it is the warm bias under water-vapour feedback, an important PF finding, not a setup artifact.** Trajectory reproduced in conversation (not saved as a script).

---

### Experiment #15 — Moist RCE (Emanuel) confirms the warm-bias prediction (2026-06-02)

`rce_moist_pf_vs_rrtmg.py`, 200 d, 5 m slab, SW=240, O₃=0:

| column | T_sfc | q_sfc g/kg | OLR | imbalance (240−OLR) |
|---|---:|---:|---:|---:|
| RRTMG-LW | 283.5 | 2.02 | 254.1 | −14 (cooling) |
| PF-LW linepyline | 296.9 | 6.77 | 233.3 | +6.7 (warming) |
| PF-LW shipped 4band | 295.8 | 6.60 | 240.5 | −0.5 (≈balanced) |

**Findings.**
1. **Prediction confirmed.** PF runs ~13 K warmer than RRTMG with ~3× the surface water vapour — the exp #14 warm-bias-under-water-vapour-feedback mechanism, now in the moist (Emanuel) configuration. Not converged (RRTMG cooling, PF-lpl warming), so the converged gap will *exceed* 13 K.
2. **Table-independent.** PF-linepyline and PF-shipped both ≈296 K → the bias does not depend on the CO₂ table (consistent with exp #11/#13: CO₂ and H₂O k are LBL-faithful). The cause is the radiative scheme itself (thin bands / two-stream integrator).
3. **q diverged 3×** (6.8 vs 2.0 g/kg) vs the ~equal ~1 g/kg at day 40 of the dry trajectory — the C–C amplifier is now fully engaged as the surfaces separated.

**Caveat — tropopause diagnostic still broken in the moist script.** p_tp=2.40 hPa (model top) for all columns: the dry script's cold-point→convective-top fix was NOT applied to the moist script, and a θ-based convective top doesn't work for moist convection anyway (θ rises along the moist adiabat). A moist tropopause needs the convective-heating top or a WMO lapse-rate definition. Numbers above for T_sfc/q/OLR are valid; p_tp/T_tp are not.

**Conclusion / pivot.** Both dry and moist RCE are dominated by PF's intrinsic over-trapping. The CO₂ and H₂O k-tables are validated and NOT the cause. The "reasonable tropopause" goal cannot be met until the over-trapping is understood. **Highest-value next step: a fixed-profile forward comparison** — run RRTMG-LW and PF-LW once on one identical (T,p,q) profile and compare per-band up/down fluxes and per-level heating rates, to localize which bands/levels PF over-traps (isolates the radiative scheme difference from the convective/feedback dynamics). This finally targets the exp #11 open question (thin bands / two-stream / source-function treatment).

---

### Experiment #16 — Fixed-profile LW forward comparison: PF over-traps by 51 W/m² (2026-06-02)

`scripts/experiments/lw_forward_pf_vs_rrtmg.py`. ONE shared profile (288 K moist troposphere → ~200 K isothermal strat, q≈15 g/kg surface decaying upward), each scheme called once, NZ=40, CO₂=376, O₃=0, emissivity=1. No time-stepping, convection, or feedback.

| quantity | RRTMG | PF (linepyline) | PF−RRTMG |
|---|---:|---:|---:|
| OLR (TOA up) | 242.6 | **191.2** | **−51.4** |
| surface back-radiation | 337.1 | **383.2** | **+46.1** |
| surface up (control) | 390.10 | 390.07 | −0.03 ✓ |

**Findings.**
1. **PF over-traps by ~51 W/m² OLR on an IDENTICAL profile** (and radiates +46 W/m² more back to the surface). Surface emission matches to 0.03 W/m² (identical T_surf + emissivity=1) → the comparison is clean and the divergence is purely atmospheric. This single number fully explains the RCE warm bias: to restore OLR=240, PF must warm/moisten until it emits 51 W/m² more.
2. **Definitively a radiative-scheme problem, not k-tables or feedback.** No feedback here, and the k-tables are LBL-validated (exps #7/#13). PF's band-mean k is actually ~0.85× LBL (slightly *less* opaque) yet it emits *less* to space — so the over-trapping is in the SOLVER, not the absorption coefficients.
3. **Localized to the troposphere.** The up-flux deficit `d_up` grows from 0 at the surface to −51 W/m² by ~150 hPa, then stays constant through the (transparent) stratosphere — i.e. the differential absorption accumulates between the surface and ~150 hPa.
4. **PF per-band TOA OLR** (total 191.2): rotation 10–500 = 58.7; 500–630 = 30.8; CO₂ core 630–700 = 6.6; 700–800 = 20.6; **window 800–1800 = 69.8**; 1800–3250 = 4.8. (RRTMG per-band not exposed, so no direct spectral diff — but the wide 800–1800 "window" band lumping the transparent 8–12 µm window with H₂O-opaque regions into one 8-g-point k-distribution is a prime suspect, the same dilution problem that motivated the CO₂ band split.)

**Mechanism hypotheses (next).** Since it's the solver: (a) diffusivity factor (1.66?) wrong or mis-applied in the two-stream LW; (b) band-Planck-fraction weighting biasing emission; (c) the 800–1800 window band too wide (dilutes window transparency); (d) two-stream/Eddington closure or layer-transmittance formulation over-absorbing. **Next step: inspect the PF LW solver (`climt/_components/picket_fence/lw/kernels.py`) — diffusivity factor, layer transmittance, Planck source — to find the 51 W/m².**

**Files added.** `scripts/experiments/lw_forward_pf_vs_rrtmg.py`, `debug_data/lw_forward_pf_vs_rrtmg.png`.

**Dry-vs-moist split (`--dry` flag).** Same profile, q zeroed (CO₂-only):

| | OLR deficit (PF−RRTMG) | back-radiation excess |
|---|---:|---:|
| Moist | −51.4 | +46.1 |
| **Dry (CO₂-only)** | **−13.3** | **+20.9** |

Removing H₂O cuts the OLR over-trapping 51 → 13 W/m². So **~38 W/m² (≈¾) is H₂O/window-related** (the wide 800–1800 window band and/or H₂O handling — suspect #2), while a **~13 W/m² residual persists dry in CO₂-only** (a genuine solver-level discrepancy independent of H₂O — suspect #1). Dry PF per-band TOA OLR: rotation 100.0, 500–630 46.2, core 6.6, 700–800 26.0, **window 800–1800 = 141.5**, 1800–3250 7.1. Both components are real; the dry CO₂-only residual (−13 W/m² OLR, +21 back-radiation) is the **cleanest "PF Python two-stream vs RRTMG Fortran" signal** (the PF solver was originally built against RRTMG), so it is the place to start a kernel-level comparison. **Next: read `climt/_components/picket_fence/lw/kernels.py` — diffusivity factor, layer transmittance, Planck source — against the standard/RRTMG two-stream.**

**Kernel inspection (done 2026-06-02).** `climt/_components/picket_fence/lw/kernels.py` (`lw_transport` / `_lw_transport_single_gpt`) uses the standard Schwarzschild two-stream with **diffusivity factor D=1.66** (`DIFFUSIVITY_FACTOR = 1.66`, lines 34/78 — CORRECT, so a wrong/double-applied diffusivity is ruled out). Per g-point: `up[k+1]=up[k]·trans + B[k]·(1−trans)`, `down[k]=down[k+1]·trans + B[k]·(1−trans)`, `trans=exp(−1.66·tau[k])`. This is the **isothermal-layer** emissivity method (one Planck value `B(T_mid)` per layer, emitted equally up/down). It is correct in form but coarser than RRTMG's **linear-in-τ Planck source** (Clough/RRTMG source function, B varies within each layer) — the prime suspect for the ~13 W/m² dry residual, especially on thick layers.

---

### Experiment #17 — Task A: window-band split at 1250 cm⁻¹ recovers ~16 of the ~38 W/m² window penalty (2026-06-02)

Executed Task A from the handoff. Split the wide 800–1800 cm⁻¹ "window" band at 1250 cm⁻¹ (transparent 8–12 µm window 800–1250 vs H₂O ν₂ 1250–1800), generated a new 7-band table, and re-ran the fixed-profile forward comparison (exp #16's diagnostic).

**Build.** `scripts/generate_pf_tables_linepyline.py` `SCENARIOS["earth"]["lw_band_edges"]` → `[10,500,630,700,800,1250,1800,3250]`. Generated `earth_low_res_lw_co2refined_linepyline_winsplit.nc` (linepyline env, validated recipe: binning=False, dnu=0.1, pseudovoigt, ngpt=8). `k_coefficients.shape=(1, 7, 8, 12, 8, 7)` — 7 bands. Added a `--table` arg to `scripts/experiments/lw_forward_pf_vs_rrtmg.py` so tables can be swapped without overwriting figures.

**Result (fixed-profile forward, PF−RRTMG OLR deficit, W/m²):**

| case | exp #16 (lumped 800–1800) | exp #17 (split at 1250) | recovered |
|---|---:|---:|---:|
| Moist OLR deficit | −51.4 | **−31.7** | **+19.7** |
| Dry (CO₂-only) OLR deficit | −13.3 | −9.8 | +3.5 |
| Moist-only (H₂O/window) component | 38.1 | 21.9 | **+16.2** |

PF per-band TOA OLR (moist), total 211.0 (was 191.2): the new transparent **800–1250 band emits 82.5 W/m²** while 1250–1800 emits only 7.1 — vs the lumped 800–1800 band's 69.8 in exp #16. Surface up matched to 0.03 W/m² (clean comparison). Figures: `debug_data/lw_forward_pf_vs_rrtmg_{moist,dry}_earth_low_res_lw_co2refined_linepyline_winsplit.png`.

**What we learned.**
1. **The wide window band is confirmed as a real, dominant culprit.** Splitting at 1250 lets the transparent window photons escape and recovers ~16 of the ~38 W/m² moist window penalty (~42%) — directionally exactly as predicted, the "band too wide" pathology that started this plan, now in the window.
2. **But it is only a PARTIAL fix.** The moist OLR deficit collapsed from −51 to −32, NOT down to the dry −13 the decision rule hoped for. ~22 W/m² of moist-specific over-trapping survives. The 800–1250 band still lumps the truly transparent 8–10 µm window with the H₂O continuum / 10–12 µm region; a single split doesn't fully separate them.
3. **The split also slightly helped the dry case** (−13.3 → −9.8, +3.5 W/m²) — even with q=0 the 800–1800 band contains CO₂-free window that escapes better once isolated. The dry CO₂-only residual is now ~−10 W/m² (Task B territory: the isothermal-layer Planck source).

**Decision (per Task A's rule).** The deficit did NOT collapse all the way to −13, so a single 1250 split is necessary but not sufficient. Two follow-ups:
- **A2 (finer window split):** try splitting the window further (e.g. 800/980/1180/1800, isolating the cleanest 8–10 µm window) and re-run the forward test. If the moist deficit approaches the dry residual, adopt that partition as the shipped table.
- **B (dry residual):** the ~−10 W/m² dry CO₂-only deficit is now the cleanest remaining signal — run the NZ=40/80/160 resolution test; if it shrinks with resolution, implement the linear-in-τ Planck source in the kernel.

**Files touched.** `scripts/generate_pf_tables_linepyline.py` (earth lw_band_edges), `scripts/experiments/lw_forward_pf_vs_rrtmg.py` (`--table` arg + table-tagged figure name). New table: `earth_low_res_lw_co2refined_linepyline_winsplit.nc`. New figures in `debug_data/`.

---

### Experiment #18 — Task A2: finer window subdivision is WORSE; the transparent window wants to be WIDE (2026-06-02)

Followed up exp #17's partial window fix by subdividing the window further to isolate the cleanest 8–10 µm window. Generated `earth_low_res_lw_co2refined_linepyline_winsplit2.nc` (8 bands, edges `[10,500,630,700,800,980,1180,1800,3250]`) and re-ran the moist forward test.

**Result — finer split made it WORSE, not better:**

| table | window edges | moist OLR deficit (PF−RRTMG) | window-region (800–1800) PF emission |
|---|---|---:|---:|
| winsplit (exp #17) | …,800,**1250**,1800,… | **−31.7** | 89.5 (800–1250=82.5 + 1250–1800=7.1) |
| winsplit2 (this exp) | …,800,**980,1180**,1800,… | **−34.3** | 86.9 (800–980=41.5 + 980–1180=33.1 + 1180–1800=12.3) |

**What we learned (counterintuitive but important).**
1. **A transparent window must be kept WIDE, not subdivided.** Splitting the single 800–1250 transparent band (82.5 W/m² escape) into 800–980 + 980–1180 dropped total window emission ~3 W/m². This is the *opposite* of the CO₂-core lesson: for an absorption feature, splitting isolates the core and helps; for a transparent window, lumping is better because the low-k g-points of the clean window carry flux across the whole band's Planck weight. Subdividing strands each sub-band with only its own (less clean) k-distribution.
2. **The 1180–1250 transparent slice got lumped into the opaque 1180–1800 H₂O ν₂ band**, losing its transparency — the same dilution pathology, now self-inflicted. The boundary that matters is transparent-window↔H₂O-band (≈1250), not subdivisions within the window.
3. **Window subdivision is exhausted.** The best window partition is the single 1250 split (winsplit). The moist deficit plateaus at ~−32, far above the dry −10. **Per Task A2's decision rule, the residual ~22 W/m² moist over-trapping is NOT purely band width** → it is the solver / H₂O handling, i.e. Task B is the lever.

**Verdict.** `winsplit` (1250) stands as the window-band deliverable; `winsplit2` REJECTED (delete unless kept as a negative-result reference). Scenario `lw_band_edges` reverted to the 1250 partition.

**Files touched.** `scripts/generate_pf_tables_linepyline.py` (edges set to winsplit2, then reverted to 1250). New table: `earth_low_res_lw_co2refined_linepyline_winsplit2.nc` (rejected). Figure: `debug_data/lw_forward_pf_vs_rrtmg_moist_..._winsplit2.png`.

---

### Experiment #19 — Task B: dry solver residual is RESOLUTION-INDEPENDENT → a two-stream formulation difference (2026-06-02)

Ran the Task B resolution test: `lw_forward_pf_vs_rrtmg.py --dry` (CO₂-only, the cleanest solver signal) at NZ=40/80/160 with the winsplit table. Added an `--nz` arg.

**Result — flat with resolution:**

| NZ | dry OLR deficit (PF−RRTMG) | dry back-radiation excess |
|---:|---:|---:|
| 40 | −9.83 | +15.87 |
| 80 | −9.86 | +15.87 |
| 160 | −9.86 | +15.87 |

RRTMG OLR also converged (340.62/340.62/340.60). Both schemes are grid-converged; the gap is constant.

**What we learned.**
1. **The ~−10 W/m² dry residual does NOT shrink with resolution → it is NOT a layer-thickness/discretization error.** Per Task B's decision rule, this **rules out the isothermal-layer Planck source as a resolution-limited cause.** Strong evidence: the CO₂ core carries ~25 optical depths *per layer* at NZ=40, so an isothermal-vs-linear-in-τ source error would be large and would halve as layers halve — it didn't move at all.
2. **The residual is a fixed formulation difference in the two-stream solver**, not a numerical-resolution artifact. Candidates (term-by-term comparison needed): the constant D=1.66 diffusivity vs RRTMG's angle/path-dependent secant; band Planck-fraction weighting; the surface/TOA boundary terms; or the transmittance functional form.
3. **The true solver bias is larger than −10.** PF's band-mean k is ~0.85× LBL (less opaque, exp #8), which alone would *raise* PF's OLR. PF over-traps anyway (lower OLR + higher back-radiation, both directions) → the solver effect is partially masked by under-opacity and is genuinely the dominant remaining error.

**Verdict / next.** Both window-width (Tasks A/A2) and resolution (Task B) are exhausted as levers. The remaining over-trapping (~22 W/m² moist + ~10 W/m² dry) lives in the **two-stream solver formulation**. Next experiment is a **term-by-term comparison of `_lw_transport_single_gpt` against RRTMG's RTRN/two-stream** — starting with a diffusivity-factor sensitivity (vary D in the kernel, see if the dry −10 closes) since that's a one-line, resolution-independent knob consistent with a constant offset, and a check of the band Planck-fraction weighting.

**Files touched.** `scripts/experiments/lw_forward_pf_vs_rrtmg.py` (`--nz` arg). No `climt/` code changed.

---

### Experiment #20 — Task D: PF solver vs rte-rrtmgp term-by-term; source function is NOT the over-trapping (2026-06-02)

Per user direction, compared the PF LW kernel against the actual reference solver it was built from, `rte-rrtmgp` (`/Users/joymonteiro/github/other_libs/rte-rrtmgp`, `rte/kernels/mo_rte_solver_kernels.F90` + `rte/frontend/mo_rte_lw.F90`), BEFORE pursuing the heavyweight Task D sub-steps.

**Term-by-term comparison (`_lw_transport_single_gpt` vs `lw_solver_noscat`):**

| element | PF | rte-rrtmgp | difference |
|---|---|---|---|
| diffusivity | constant D=1.66 | default 1 angle, D=1/0.6097=**1.64** (Hogan 2023 Gauss-Jacobi-5) | negligible (~1%) |
| angular quadrature | single angle | default `n_quad_angs=1` | same |
| Planck source per band | `planck_frac[b,g,T]·σT⁴` lookup | Planck-fraction lookup | same approach ✓ |
| **source function** | **isothermal-layer**: `B(T_mid)·(1−trans)` up=down, midpoint only | **linear-in-τ** (Clough 1992 Eq 13): edge-anchored `(1−trans)·B(T_edge)+2·fact·(B_mid−B_edge)` | **the only real difference** |

This reordered the planned Task D: D1 (diffusivity) predicted null, D2 (Planck-fraction) a non-issue (PF already uses `planck_frac`), and D3 (linear-in-τ) identified as the sole formulation difference.

**D1 — confirmed null.** Setting `DIFFUSIVITY_FACTOR` 1.66→1.64 moved the dry deficit −9.83 → −9.65 (+0.18 W/m²). Diffusivity is not the lever. (rte-rrtmgp's secant 1.64 ≈ PF's 1.66.)

**D3 — implemented and TESTED; radiatively inert.** Implemented the Clough linear-in-τ source in the PF kernel (`kernels.py`): added interface-level Planck (`lev_src`) computed from log-p-interpolated interface temperatures (`_interface_temperatures` in `component.py`), and the edge-anchored source with the 3rd-order `fact` series for small τ. Backward-compatible (lev_source=None → exact isothermal fallback; all 20 PF kernel/validation tests pass; parmentier path unchanged).

Forward result (winsplit table): dry −9.83 → **−9.86**, moist −31.65 → **−32.02**. Essentially unchanged (~0.4 W/m²).

Resolution sweep (dry, PF OLR absolute) confirms it's inert even at GCM-coarse resolution:

| NZ | linear-in-τ | isothermal | RRTMG |
|---:|---:|---:|---:|
| 10 | 330.54 | 330.96 | 340.52 |
| 20 | 330.72 | 330.83 | 340.62 |
| 40 | 330.77 | 330.79 | 340.62 |

Both source treatments converge to ~330.8 from opposite sides (Δ ≤ 0.4 W/m²); the −10 W/m² deficit vs RRTMG persists at all resolutions. **This confirms exp #19: the source function is NOT the over-trapping mechanism.** The term-by-term comparison found the only solver-formulation difference and proved it inert.

**Decision pending (user):** keep linear-in-τ (more correct physics, matches RRTMG/RRTMGP, passes tests, ~0.4 W/m², but adds interface-T machinery) vs revert (simpler; inert here). Code currently retains it.

**Reframe — what the over-trapping is NOT (now exhaustively ruled out):** k-table band-means (exps #7/#11/#13), CO₂/window band width beyond the 1250 split (exps #17/#18), diffusivity (D1), source function (D3), vertical resolution (exp #19 + this). **What remains:** the comparison is NOT "same k, different solver" — RRTMG runs its own 16-band/140-g-point gas optics while PF runs a 7-band/8-g-point correlated-k table. The leading remaining hypothesis is the **correlated-k band STRUCTURE itself** (few wide bands + 8 g-points violate the spectral-correlation assumption and over-trap vs RRTMG's 16 bands) and/or genuine PF-table vs RRTMG-gas-optics differences. The window split helping +16 W/m² is direct evidence band structure matters.

**Decisive next experiment (#21): per-band OLR localization vs a line-faithful reference.** Run rte-rrtmgp (or RRTMG per-band) on the IDENTICAL fixed profile and compare per-band TOA up-flux to PF's per-band OLR (already emitted by `lw_forward_pf_vs_rrtmg.py`). This localizes WHICH bands PF over-traps — H₂O rotation (10–500), H₂O ν₂ (1250–1800), CO₂, or broadband — and distinguishes "band structure" from "table difference." That, not more solver work, is the way forward.

**Files touched.** `climt/_components/picket_fence/lw/kernels.py` (linear-in-τ source, `lev_planck_source` arg), `climt/_components/picket_fence/lw/component.py` (`_interface_temperatures`, interface Planck `lev_src`, wiring). `scripts/experiments/lw_forward_pf_vs_rrtmg.py` (`--nz` arg). No table changes.

---

### Experiment #21 — LBL OLR tiebreak: PF over-traps vs ground truth; it is the correlated-k BAND STRUCTURE (2026-06-02)

The user's decisive framing: the PF−RRTMG OLR gap (exp #16) could mean PF over-traps OR RRTMG under-traps. Add **linepyline (LBL) as a third reference** on the IDENTICAL fixed profile to settle which scheme the truth agrees with. Crucially, run the LBL two-stream at **D=1.66 (PF's diffusivity)** so the ONLY difference between LBL and PF is full spectral resolution vs PF's 7-band/8-g-point correlated-k — and the LBL is built by the SAME linepyline that made PF's k-table, so the k-data is identical too.

Scripts: `scripts/experiments/dump_forward_profile.py` (climt env — saves the exp #16 profile + PF/RRTMG OLR via one-RRTMG-per-process subprocess calls to the forward script; note a second RRTMG instance in one process returns garbage) and `scripts/experiments/lbl_olr_tiebreak.py` (linepyline env — `rtm.radiative_transfer` over 10–3250 cm⁻¹, dnu=0.1, pseudovoigt, binning=False, MT_CKD on, D=1.66; OLR binned onto PF band edges).

**Result — the LBL truth sits with RRTMG, not PF:**

| case | RRTMG | PF | **LBL** | LBL−PF | LBL−RRTMG | LBL closer to |
|---|---:|---:|---:|---:|---:|---|
| moist | 242.62 | 210.97 | **246.76** | **+35.8** | +4.1 | **RRTMG** |
| dry | 340.62 | 330.79 | **347.89** | **+17.1** | +7.3 | **RRTMG** |

**PF genuinely over-traps** (~36 W/m² moist, ~17 W/m² dry vs line-by-line). RRTMG is approximately correct (within ~4–7 of LBL). This rejects any lingering "maybe PF is right and RRTMG/O₃ was the artifact" — confirmed against ground truth.

**Per-band localization (LBL−PF OLR, W/m²; >0 = PF over-traps that band):**

| band (cm⁻¹) | moist Δ | dry Δ |
|---|---:|---:|
| 10–500 (H₂O rotation) | **+8.9** | **+10.6** |
| 500–630 | +3.7 | +0.9 |
| 630–700 (CO₂ core) | +0.5 | +0.5 |
| 700–800 | +3.4 | +1.0 |
| **800–1250 (window)** | **+16.4** | +0.3 |
| 1250–1800 (H₂O ν₂) | +4.3 | +3.9 |
| 1800–3250 | −1.3 | −0.2 |

**What we learned (the root cause, finally).**
1. **The over-trapping IS the correlated-k band structure.** Same k-source, same diffusivity, same profile, same two-stream family — the only difference is 7 wide bands / 8 g-points vs full spectral. That difference alone produces +17 to +36 W/m² of over-trapping. Not the k-tables (exps #7–13), not the solver source function (exp #20), not resolution (exp #19), not diffusivity (D1). **The discretization itself over-traps.**
2. **It concentrates in bands that mix transparent + absorbing spectra.** Moist: the 800–1250 window is the single worst (+16.4) — even after the exp #17 split, the band still lumps the transparent 8–12 µm window with the H₂O continuum, and correlated-k with 8 g-points can't let the transparent photons fully escape. The H₂O rotation (10–500, +9) and ν₂ (1250–1800, +4) are next. The CO₂ core (630–700) is nearly perfect (+0.5) — narrow and spectrally homogeneous.
3. **Even dry, the 10–500 band over-traps +10.6** — the wide H₂O-rotation band is spectrally inhomogeneous (CO₂ far-wing + window mix) and over-traps regardless of moisture.
4. **The exp #18 paradox is now fully explained.** Splitting the transparent window the WRONG way (lumping a transparent slice into the H₂O band) hurt; but the window genuinely over-traps (+16). The fix is to isolate a SPECTRALLY HOMOGENEOUS transparent window cleanly — not to subdivide arbitrarily.

**Implication / next.** The fundamental tension is the correlated-k spectral-correlation assumption breaking down for few wide, inhomogeneous bands. Paths: (a) more, spectrally-homogeneous bands (carefully — must separate transparent-vs-opaque, not split within a regime; the window and H₂O-rotation bands are the targets); (b) more g-points specifically in the inhomogeneous bands (exp #8 said band-MEAN k converges by ngpt~16, but FLUX trapping may need more — untested for flux); (c) accept ~10–35 W/m² as the price of a fast 7-band scheme and document it (the PLASIM lesson: tune/accept for GCM use). This is a design decision for the user, not a bug to chase.

**Files added.** `scripts/experiments/dump_forward_profile.py`, `scripts/experiments/lbl_olr_tiebreak.py`, `debug_data/forward_profile.npz`, `debug_data/lbl_olr_spec_{moist,dry}.npz`. No `climt/` code changed (linear-in-τ from exp #20 was reverted per user; the pre-existing diagnostic-rename refactor in `component.py` was preserved).

---

### Experiment #22 — Task F(b): more g-points does NOT reduce flux over-trapping → it's band width, not g-resolution (2026-06-02)

Tested whether the exp #21 over-trapping is g-point under-resolution (fixable by more g-points, same bands) or a fundamental wide-band spectral-correlation breakdown (needs narrower bands). Regenerated the winsplit table (same 7-band 1250 partition) at ngpt=16 and 32, ran the forward comparison, compared PF OLR to the fixed LBL truth (moist 246.76, dry 347.89).

| ngpt | moist PF OLR | dry PF OLR |
|---:|---:|---:|
| 8 | 210.97 | 330.79 |
| 16 | 210.98 | 330.55 |
| 32 | 210.96 | 330.50 |
| **LBL** | **246.76** | **347.89** |

**PF OLR is flat across ngpt 8→16→32 (≤0.3 W/m²).** The ~36 (moist)/~17 (dry) W/m² over-trapping vs LBL is unchanged.

**Conclusion.** g-point resolution is REJECTED as the lever (Task F option b is dead). The over-trapping is intrinsic to the **band width / spectral inhomogeneity**: within a wide inhomogeneous band, a single g-point's k maps to different wavenumbers at different layers (the correlated-k correlation assumption fails); adding g-points just samples the same broken within-band CDF more finely. Consistent with exp #8 (band-MEAN k converges by ngpt~16 while FLUX trapping does not improve at all). **Only narrower, spectrally-homogeneous bands (Task F option a) can reduce it** — or accept/tune it (option c).

**Files added (rejected, negative-result references).** `earth_low_res_lw_co2refined_linepyline_winsplit_{16,32}gpt.nc` — keep only as evidence the ngpt lever is dead; safe to delete. (Also `..._winsplit2.nc` from exp #18.)

---

### Experiment #23 — GOAL + band-structure exploration: splitting fixes DRY, not MOIST (2026-06-02)

**GOAL (user-set).** Design a PF LW band structure whose single-column TOA OLR matches line-by-line truth as well as RRTMG does: target **|LBL − PF| < 5 W/m²** for both moist and dry on the fixed exp#16 profile (RRTMG itself achieves 4.1 moist / 7.3 dry, exp #21). Downstream validation: the converged moist-RCE surface-T gap vs RRTMG should drop to within a few K. Rationale (user): poor single-column RCE → wrong tropical tropopause + amplified water-vapour feedback; "reasonable circulation despite bad radiation" is not assumed. Keep the band count as low as hits the target.

**Recent surface-T gap (answering "how far from RRTMG?").** Most recent RCE is exp #15 (moist, 200 d, 5 m slab, SW=240, O₃=0, the pre-winsplit 6-band linepyline table): RRTMG T_sfc = 283.5 K, PF = **296.9 K → +13.4 K and widening** (RRTMG cooling, PF warming, not converged). Table-independent (PF-shipped also ≈296). The ~+13 K maps to the ~30–36 W/m² single-column over-trapping (sensitivity ≈ 2.3 W/m²/K). No converged RCE for the winsplit table yet.

**Tooling.** Added `--bands` and `--no-continuum` to `scripts/generate_pf_tables_linepyline.py`; new `scripts/experiments/eval_band_structure.py` (climt env) computes PF per-band TOA OLR on the fixed profile and re-bins the saved LBL spectrum onto the candidate's band edges — no linepyline re-run per iteration (the LBL spectrum is profile-only).

**Candidates tested (ngpt=8, vs LBL):**

| table | bands | edges | moist LBL−PF | dry LBL−PF |
|---|---|---|---:|---:|
| winsplit | 7 | …,800,1250,1800,… | +35.8 | +17.1 |
| candA | 10 | 10,350,500,630,700,800,1050,1250,1500,1800,3250 | +30.3 | **+12.5** |
| candB | 12 | 10,250,500,630,700,800,980,1080,1180,1250,1500,1800,3250 | +30.6 | **+11.8** |

**Findings.**
1. **Band splitting fixes DRY (17→12) but barely moves MOIST (36→30).** The dry window 800–1250 splits to near-perfect (candB sub-bands +0.2/+0.1/+0.05/+0.02). But the **moist window stays +16.7 regardless of how finely it's split** (candA 2-way +11.4/+5.3; candB 4-way +8.8/+3.7/+2.7/+1.6 — same total). Every moist band over-traps ~+3–4 where dry over-traps ~0 → the moist excess is uniformly H₂O.
2. The CO₂ core (630–700) is perfect (+0.5) in every partition; the far-IR (10–500) responds to splitting (dry 10–250 −0.5, 250–500 +3.4 in candB).

**Conclusion.** The DRY over-trapping is band-structure-fixable (toward RRTMG's 7.3 with ~12 bands). The MOIST over-trapping is dominated by an H₂O effect that band width cannot touch → diagnosed decisively in exp #24.

### Experiment #24 — DECISIVE: the moist over-trapping is the H₂O MT_CKD CONTINUUM (2026-06-02)

Tested whether the immovable moist window over-trapping is the H₂O continuum or H₂O lines. Regenerated candA with `--no-continuum` and re-ran the LBL (`lbl_olr_tiebreak.py --no-continuum`) so both PF and reference omit the MT_CKD continuum; compared with `eval_band_structure.py <table> _nocont`.

| quantity | candA (continuum ON) | candA_nocont (continuum OFF) |
|---|---:|---:|
| moist total LBL−PF | +30.3 | **+12.85** |
| moist window 800–1050 LBL−PF | +11.4 | **+1.2** |
| moist window 1050–1250 LBL−PF | +5.3 | **+1.4** |
| dry total LBL−PF | +12.5 | +12.2 (unchanged) |

**Removing the continuum collapses the moist window over-trapping +16.7 → +2.5, and the moist total +30 → +13 (≈ dry).** So:
1. **~17 W/m² of the moist over-trapping is the H₂O MT_CKD continuum** being over-represented by the PF correlated-k. It is NOT band width (exp #23) and NOT g-points (exp #22) — it is *how the smooth continuum is folded into the line-sorted k-distribution*. Mechanism: the linepyline builder sums continuum into per-ν kappa BEFORE the correlated-k sort, so the line-driven g→ν mapping scrambles the continuum's smooth spectral structure and over-fills the transparent-window low-k g-points; and its quadratic-in-H₂O scaling is captured only by 7 X-VMR nodes + log-linear interpolation. **Contrast with RRTMG (verified in `rrtmg_lw_taumol.f90:331-342` + `rrtmg_lw_setcoef.f90:296-410`):** RRTMG keeps the continuum as a SEPARATE ADDITIVE term `taug(lay,ig) = colh2o*absa(...,ig) [lines] + tauself + taufor [continuum]`, where `tauself/taufor` use per-g-point reference coeffs `selfref(:,ig)/forref(:,ig)` scaled by `selffac ∝ H₂O²·p` / `forfac ∝ H₂O·p` — i.e. the continuum's quadratic-H₂O/pressure scaling is analytic and exact at every layer, never co-sorted with the lines. (NOTE: RRTMG's continuum is per-g-point, weakly g-dependent — NOT grey; an earlier note here mischaracterized it as grey.)
2. **The remaining ~12 W/m² (moist, continuum-off ≈ dry) is band-structure-fixable** (exp #23) — far-IR/ν₂ refinement.

**What it takes to hit the GOAL (|LBL−PF|<5, both cases):**
- (i) **~12 spectrally-homogeneous bands** (candB-style far-IR + window + ν₂ splits) → closes the dry/non-continuum part toward ~8–10 W/m².
- (ii) **Separate H₂O-continuum treatment** in the table builder: build the line-only k-distribution (sort lines WITHOUT continuum), and store the MT_CKD continuum as its own term added to every g-point's τ at runtime with explicit self(∝H₂O²·p)/foreign(∝H₂O·p) scaling — mirroring RRTMG's `tauself+taufor` structure (`rrtmg_lw_taumol.f90:331`). A band-grey-across-g continuum is an acceptable approximation (RRTMG's per-g coeffs are only weakly g-dependent); the essential change is decoupling it from the line CDF. → recovers the ~17 W/m² moist continuum penalty.
- Neither alone suffices for moist; both together should bring moist+dry to within a few W/m² of LBL (≈ RRTMG), i.e. the RCE surface-T gap from ~+13 K toward ~+2–3 K.

**Files.** `scripts/experiments/eval_band_structure.py` (new), `lbl_olr_tiebreak.py` (+`--no-continuum`), `generate_pf_tables_linepyline.py` (+`--bands`,`--no-continuum`). Experiment tables (negative/diagnostic, deletable): `earth_lw_candA.nc`, `earth_lw_candB.nc`, `earth_lw_candA_nocont.nc`, `..._winsplit_{16,32}gpt.nc`, `..._winsplit2.nc`. LBL spectra: `debug_data/lbl_olr_spec_{moist,dry}{,_nocont}.npz`.

**Next.** Implement the separate-continuum table-builder path (ii) and re-eval; in parallel adopt a candB-style ~12-band partition (i). Then run a converged moist RCE with the resulting table to confirm the surface-T gap closes. (The 7-band shipped winsplit table stays as the fast default per the user.)

---

### Experiment #25 — Decoupled + log-X-interpolated H₂O continuum: moist over-trapping +36→+15 (2026-06-02, autonomous)

Implemented the continuum decoupling (the exp #24 fix) and found the *real* mechanism is the X-axis interpolation, not the decoupling alone.

**Implementation (backward-compatible; absent `continuum_kappa` → unchanged):**
- `scripts/generate_pf_tables_linepyline.py` `--decouple-continuum`: builds the k-distribution from LINE-ONLY kappa (continuum off), computes the MT_CKD continuum as `full − lines`, band-grey-averages it over ν → `continuum_kappa(nband,nT,nP,nX)`, written by `write_lw_table(continuum_kappa=…)` (`pf_table_builder/netcdf_writer.py`). Also added `--bands` and `--no-continuum`.
- `climt/_components/picket_fence/optics/correlated_k.py`: loads `continuum_kappa`; new `interpolate_continuum` adds it grey-across-g to every g-point's τ in `_compute_ck_optical_depth_additive` with the air-column-mass scaling.
- 28 PF tests pass; existing tables unaffected.

**Two-step result on candC (= candB's 12 bands + decoupled continuum), vs LBL:**

| variant | moist LBL−PF | dry LBL−PF |
|---|---:|---:|
| winsplit 7-band (folded continuum) | +35.8 | +17.1 |
| candB 12-band (folded continuum) | +30.6 | +11.8 |
| candC decoupled, **grey/linear-X interp** | +35.8 (worse!) | +11.8 |
| candC decoupled, **log-X interp** | **+15.0** | +11.7 |

**Key finding — it's the X-INTERPOLATION, not the decoupling per se.** The H₂O self-continuum scales ~X². Stored on the widely log-spaced 7-node VMR grid and interpolated **linear-in-value** (what the line path does), a convex X² quantity is wildly over-estimated between nodes (~6× at surface X≈0.024, between the 1e-2 and 1e-1 nodes) → window over-traps. Decoupling alone with linear-X interp was *worse* (+35.8, grey smearing). Decoupling + **log-space (log-k vs log-X) interpolation** — exact for a power law — collapsed the moist window bands 800–1250 from +8.9 each to ~+0.5 and **moist total +35.8 → +15.0** (≈21 W/m² recovered). Dry unchanged (no continuum). This is PF's analogue of RRTMG's analytic `selffac ∝ H₂O²` scaling.

**Where we are vs the GOAL (|LBL−PF|<5):** moist +15.0, dry +11.7 (from +35.8/+17.1). The continuum piece is solved. Remaining residual is band-structure/H₂O-lines: moist offenders now 250–500 (+5.3, far-IR H₂O rotation), 1250–1500 (+3.6, ν₂); dry residual ~+12 is CO₂ correlated-k spread (10–250 +2.6, 250–500 +3.4, 1500–1800 +2.0) — partly more-bands-fixable, partly the intrinsic finite-band/8-gpt correlated-k floor (RRTMG uses 16 bands/140 g-pts to reach 7.3/4.1).

**Note — the line k-distribution likely has the SAME linear-X overshoot** (lines scale ~X¹, also convex on a log axis). Switching the line X-interpolation to log-k may further help moist; untested (separate from the continuum fix). Candidate next.

**Files.** `generate_pf_tables_linepyline.py` (+`--bands`,`--no-continuum`,`--decouple-continuum`), `pf_table_builder/netcdf_writer.py` (+`continuum_kappa`), `correlated_k.py` (+`continuum_kappa` load, `interpolate_continuum` log-space, additive τ add). Table: `earth_lw_candC_decoupcont.nc` (12-band, decoupled continuum). Validation deferred to RCE (in progress).

---

### Experiment #26 — Band refinement + ngpt floor: "what it takes" quantified (2026-06-02, autonomous)

With the continuum decoupled (exp #25), refined the band structure and re-probed g-points.

| table (all decoupled+log-X continuum unless noted) | bands | ngpt | moist LBL−PF | dry LBL−PF |
|---|---|---|---:|---:|
| winsplit (folded continuum) | 7 | 8 | +35.8 | +17.1 |
| candC | 12 | 8 | +15.0 | +11.7 |
| **candD** (+far-IR 250/350, ν₂ 1400/1600 splits) | 14 | 8 | **+11.5** | **+11.0** |
| candD | 14 | **16** | +11.5 | +11.4 |

**Findings.**
1. **More spectrally-homogeneous bands keep helping** (dry 17.1→11.7→11.0; moist 35.8→15.0→11.5). With 14 bands the moist and dry residuals converge (~+11) — the H₂O-specific over-trapping (continuum + extra lines) is essentially eliminated; what remains is a band-structure residual common to both.
2. **g-points do NOT close the floor (clean confirmation).** ngpt 8→16 on candD: +11.5/+11.4 ≈ unchanged. exp #22's null (folded continuum) was confounded; this is the clean repeat — the floor is NOT within-band CDF resolution.
3. **The ~+11 floor is the correlated-k spectral-correlation limit.** It concentrates in MODERATELY-opaque bands (350–500, 500–630, 700–800, 800–980 each ~+1–2.6) where the k-distribution's column-correlation breaks down; the saturated CO₂ core (+0.5) and the near-transparent 10–250 (~0) are fine. Only narrower/more bands reduce it (RRTMG reaches 4.1/7.3 with 16 bands/140 g-pts).

**"WHAT IT TAKES" (answer to the user's question), from +35.8/+17.1 (moist/dry) to the GOAL <5:**
- **Decoupled + log-X-interpolated H₂O continuum** — the single biggest lever, −21 W/m² moist (exp #25). Essential.
- **~14+ spectrally-homogeneous bands** — far-IR (10/250/350/500), full CO₂ split, window (800/980/1080/1180/1250), ν₂ (1250/1400/1600/1800). Gets to ~+11/+11.
- **g-points stay at 8** (16 doesn't help).
- **To reach +5** (RRTMG fidelity): push toward ~16–20 bands in the moderate-opacity regions; this is the residual correlated-k floor, diminishing returns. ~+11 may be an acceptable stopping point for a fast scheme (≈ +5 K RCE surface bias vs RRTMG, less than half the current +13 K).

**Files.** Tables `earth_lw_candD_decoupcont.nc` (14-band, ngpt8 — current best), `earth_lw_candD_ngpt16.nc` (ngpt null-test, deletable). RCE validation in progress.

---

### Experiment #27 — MOIST RCE validation: the continuum fix closes the surface-T gap to +2 K (2026-06-02, user-run)

Converged moist RCE (`rce_moist_pf_vs_rrtmg.py --days 200`, Emanuel convection, SW=240, O₃=0, 5 m slab), comparing RRTMG vs PF-candD (decoupled+log-X continuum, 14-band) vs PF-winsplit (7-band, folded continuum, current shipped).

| column | OLR | T_sfc | gap vs RRTMG | q_sfc (g/kg) |
|---|---:|---:|---:|---:|
| RRTMG-LW | 254.10 | 283.54 | — | 2.02 |
| **PF candD (decoupled continuum)** | 254.81 | **285.58** | **+2.04 K** | 2.38 |
| PF winsplit (folded continuum) | 249.81 | 290.52 | +6.98 K | 3.60 |

**The single-column continuum fix translates directly to RCE.** Surface-T gap vs RRTMG: **+13.4 K (original 6-band, exp #15) → +7.0 K (winsplit 7-band) → +2.0 K (candD)**. candD's OLR now matches RRTMG to 0.7 W/m² (winsplit still traps ~4.3 W/m² more → +7 K). The water-vapour feedback that amplified the old bias is tamed: candD q_sfc 2.38 vs RRTMG 2.02 (winsplit 3.60, ~1.8× RRTMG). This confirms exps #24–26 end-to-end: the decoupled+log-X continuum + 14-band partition brings PF moist RCE to within ~2 K of RRTMG.

**Caveat (200 d):** OLR≈254 > SW=240 → not fully equilibrated; the +2.0 K was a co-evolving snapshot.

**CONVERGED update (1500 d, user-run; OLR≈240 for all):**

| column | OLR | T_sfc | ΔT_sfc | p_tp | T_tp | q_sfc |
|---|---:|---:|---:|---:|---:|---:|
| RRTMG-LW | 240.06 | 276.29 | — | 121.9 | 153.7 | 0.95 |
| **PF candD** | 240.09 | 279.03 | **+2.74 K** | 121.9 | 156.3 | 1.33 |
| PF winsplit | 239.87 | 285.92 | +9.63 K | 121.9 | 161.2 | 2.35 |

Fully equilibrated, the surface gaps firm up to **candD +2.74 K, winsplit +9.63 K** (snapshot had +2.0/+7.0). candD's converged q_sfc 1.33 vs RRTMG 0.95 (winsplit 2.35 — runaway). The fix holds at equilibrium.

**Converged CURRENT-DEFAULT comparison (1500 d, RRTMG vs shipped 4-band `earth_low_res_lw_8gpt` vs candD):**

| column | OLR | T_sfc | ΔT_sfc | p_tp | T_tp | q_sfc |
|---|---:|---:|---:|---:|---:|---:|
| RRTMG-LW | 240.06 | 276.29 | — | 122 | 153.7 | 0.95 |
| **PF shipped default (4-band)** | 238.31 | **296.38** | **+20.1 K** | 583* | 259.8* | 6.94 |
| PF candD (improved) | 240.09 | 279.03 | +2.74 K | 122 | 156.3 | 1.33 |

**The shipped 4-band default runs +20.1 K warmer than RRTMG** at equilibrium (~7× surface q — full WV-feedback runaway; OLR 238.3 still slightly warming so +20 K is a floor). The improved candD table cuts this to +2.7 K. Strong motivation for shipping the improved table (opt-in at minimum). *Caveat: the default's p_tp=583 hPa / T_tp=259.8 K is NOT a real tropopause — the constant-θ+10 K criterion trips mid-troposphere because the runaway-warm moist column's θ rises 10 K along the moist adiabat by then. The θ-offset tropopause is not robust across very different climate states (fine for RRTMG/candD at 122 hPa); a WMO lapse-rate or θ-curvature criterion is needed for cross-climate tropopause comparison.

**Tropopause diagnostic fixed (constant-θ).** The moist script's `tropopause()` was a cold-point search (`argmin T` above 500 hPa) which always returned the model top (2.40 hPa) because T decreases monotonically to TOA with O₃=0. Replaced with the convective-θ-top criterion (first level where θ > θ_surface + 10 K; +10 vs the dry script's +2 because the moist adiabat drifts θ a few K through the troposphere). Recomputed from the saved `.npz`: tropopause **~122 hPa for all three** (same grid level at 40 layers — pressure is grid-limited), with **T_tp discriminating: RRTMG 153.7 < candD 156.3 (+2.6) < winsplit 161.2 (+7.5)** — mirroring the surface bias. (Code: `rce_moist_pf_vs_rrtmg.py:tropopause`.)

**Status of GOAL:** single-column |LBL−PF| candD = +11.5 moist / +11.0 dry (vs RRTMG's 4.1/7.3); RCE surface gap +2.0 K. The "few K" downstream target is met; the single-column <5 W/m² target needs more bands (residual correlated-k floor). **`earth_lw_candD_decoupcont.nc` is the recommended improved Earth-LW table** (7-band winsplit remains the fast default unless the +5 K reduction is wanted).

**Pending:** dry candB RCE (user running separately) — will confirm the dry/band-structure part. → DONE, exp #28.

### Experiment #28 — DRY RCE: candB surface within +1 K of RRTMG; tropopause ~50 hPa higher (2026-06-02, user-run)

Converged dry RCE (`rce_dry_pf_vs_rrtmg.py --days 400`, SimplePhysics no-BL held dry + DryConvectiveAdjustment, 1 m slab, SW=240, O₃=0), RRTMG vs PF-candB (12-band; continuum inert when dry, so this isolates the line k-distribution / band structure).

| column | OLR | imbal | T_sfc | gap | p_cnvtop | T_tp | q_max |
|---|---:|---:|---:|---:|---:|---:|---:|
| RRTMG-LW | 242.61 | −2.61 | 266.51 | — | 379 hPa | 201.1 | 0.000 |
| PF-LW candB (12-band) | 242.27 | −2.27 | 267.48 | **+0.97 K** | 330 hPa | 195.5 | 0.000 |

**Findings.**
1. **Dry surface-T gap is ~+1 K** (OLR matches RRTMG to 0.34 W/m²; both imbal ≈ −2.4, near-converged at 400 d). q_max=0.000 confirms the exp #14 truly-dry fix holds. The 12-band line k-distribution tracks RRTMG's surface energy balance almost exactly.
2. **Residual aloft:** candB's dry convective top is ~50 hPa higher (330 vs 379 hPa) and the tropopause ~5.6 K colder (195.5 vs 201.1 K) — the band-structure correlated-k residual (exp #26's ~+11 W/m² floor) manifesting in upper-tropospheric structure, relevant to the tropical-tropopause goal. Smaller than the surface gap but present; reducible only with more bands.

**Combined RCE picture (both fixes validated):** DRY candB +1.0 K surface; MOIST candD +2.0 K surface (exp #27) — both vs RRTMG, down from the original +13 K. The decoupled+log-X continuum (moist) + 14-band partition deliver a PF Earth-LW scheme that matches RRTMG RCE to ~1–2 K. Remaining work for the tropopause is more bands (diminishing returns).

---

## ════════ SESSION HANDOFF (2026-06-02) — READ THIS FIRST IN A CLEAN SESSION ════════

The next session has no memory of this work. Everything needed is below; experiment entries #1–#19 above hold the detail.

> **UPDATE (exps #17–#19, 2026-06-02): Tasks A, A2, B are all DONE. Both levers exhausted — the remaining error is the two-stream SOLVER.**
> - **Task A (exp #17):** window split at 1250 cm⁻¹ cut the moist OLR deficit −51 → −32 (recovered ~16 of ~38 W/m²) and dry −13 → −10. Partial. Table `earth_low_res_lw_co2refined_linepyline_winsplit.nc` (7 bands) — the window-band deliverable.
> - **Task A2 (exp #18):** finer window subdivision made it WORSE (−34). A transparent window must be kept WIDE, not split (opposite of the CO₂-core lesson). Window subdivision exhausted; `winsplit2` rejected.
> - **Task B (exp #19):** the dry residual is RESOLUTION-INDEPENDENT (−9.83/−9.86/−9.86 at NZ=40/80/160) → NOT the isothermal-layer Planck source / not a discretization error → a fixed two-stream **formulation** difference vs RRTMG.
> - **Task D (exp #20):** term-by-term vs rte-rrtmgp. Only solver-formulation difference = isothermal-layer vs linear-in-τ source. D1 (diffusivity) null; rte-rrtmgp's secant is 1.64≈PF's 1.66. **Linear-in-τ implemented (kernel + interface Planck) and TESTED — radiatively inert (~0.4 W/m²) at NZ=10/20/40.** Source function is NOT the cause (confirms exp #19). Decision pending: keep (correct physics, passes tests) vs revert (simpler, inert).
> - **Task E (exp #21) — DECISIVE, the root cause is found.** Added linepyline LBL as a 3rd reference on the identical profile at D=1.66 (PF-matched). **LBL OLR agrees with RRTMG, not PF** (moist: LBL 246.8 vs RRTMG 242.6 vs PF 211.0; dry: LBL 347.9 vs RRTMG 340.6 vs PF 330.8). PF genuinely over-traps ~36 (moist)/~17 (dry) W/m² vs ground truth. Same k-source, same D, same profile → **the cause is purely PF's correlated-k BAND STRUCTURE** (7 wide bands/8 g-pts). Localized: worst in spectrally-inhomogeneous bands mixing transparent+opaque — 800–1250 window (+16.4 moist), 10–500 H₂O rotation (+9–11), 1250–1800 ν₂ (+4); CO₂ core 630–700 near-perfect (+0.5).
> - **Task F(b) g-points — REJECTED (exp #22):** ngpt 8→16→32 leaves PF OLR flat; over-trapping is band width, not g-resolution.
> - **Band-structure exploration (exp #23) + GOAL set:** target |LBL−PF| < 5 W/m² (moist & dry) = RRTMG's own fidelity. Splitting (candA 10-band, candB 12-band) fixes DRY (17→12) but barely moves MOIST (36→30): the moist window over-traps +16.7 regardless of how finely split.
> - **Exp #24 DECISIVE — the moist excess is the H₂O MT_CKD CONTINUUM.** Continuum OFF: moist window +16.7→+2.5, moist total +30→+13 (≈dry). So moist over-trapping = ~17 W/m² continuum (folded into PF's line-sorted k-distribution; RRTMG instead adds it as a SEPARATE per-g-point term `tauself+taufor` with explicit H₂O²·p/H₂O·p scaling — `rrtmg_lw_taumol.f90:331`, NOT grey) + ~12 W/m² band-structure-fixable.
> - **WHAT IT TAKES (→ goal):** (i) ~12 spectrally-homogeneous bands (candB-style) for the non-continuum part + (ii) separate H₂O-continuum treatment in the table builder for the ~17 W/m² continuum part. Both together → moist+dry within a few W/m² of LBL → RCE surface-T gap ~+13 K → ~+2–3 K.
> - **Recent surface-T gap:** exp #15 moist RCE PF = 296.9 K vs RRTMG 283.5 K = **+13.4 K, widening** (pre-winsplit table; no converged winsplit RCE yet).
> - **Exps #25–26 (autonomous) — the continuum fix works; "what it takes" quantified.** Implemented decoupled H₂O continuum (line-only k-dist + separate band-grey `continuum_kappa`, added at runtime). Key: it must be **log-X interpolated** (self-continuum ~X²; linear-in-value over log-spaced X-nodes overshoots ~6× → grey/linear was WORSE). With log-X: **moist over-trapping +35.8→+15.0** (window 800–1250 +8.9→+0.5 each). Adding far-IR/ν₂ band splits (candD, 14 bands): **moist +11.5, dry +11.0** (H₂O-specific excess gone; moist≈dry). g-points 8→16 don't help (clean confirm of exp #22). Residual ~+11 = correlated-k spectral-correlation floor in moderate-opacity bands → needs ~16+ bands to approach the +5 goal. **Best table: `earth_lw_candD_decoupcont.nc` (14-band, decoupled continuum, ngpt8).**
> - **VALIDATED in RCE (exps #27–28, user-run):** MOIST candD surface-T gap vs RRTMG, **CONVERGED (1500 d, OLR≈240): +13.4 K (orig) → +9.6 K (winsplit) → +2.7 K (candD)**; q_sfc 1.33 vs 0.95 (WV feedback tamed; winsplit runaway 2.35). DRY candB **+1.0 K** surface (OLR match 0.34, q_max=0). Constant-θ tropopause (diagnostic fixed): ~122 hPa all columns (grid-limited at 40 lev), T_tp RRTMG 153.7 < candD 156.3 < winsplit 161.2 — tracks the surface bias. Residual: candB dry convective top ~50 hPa higher / ~5.6 K colder (band-structure floor, needs more bands). `earth_lw_candD_decoupcont.nc` = recommended improved table.
> - **Code changes (all backward-compatible):** `generate_pf_tables_linepyline.py` (+`--bands`,`--no-continuum`,`--decouple-continuum`), `pf_table_builder/netcdf_writer.py` (+`continuum_kappa`), `correlated_k.py` (+`continuum_kappa` load, `interpolate_continuum` log-space, additive τ-add), new `scripts/experiments/eval_band_structure.py`. 28 PF tests pass. linear-in-τ from exp #20 reverted; diagnostic-rename refactor preserved.
> - **CANDIDATE NEXT:** the LINE k-distribution likely shares the linear-X overshoot (lines ~X¹, also convex) — switching line X-interp to log-k may further cut moist (untested, riskier: log(0), affects all tables).

### Bottom line (where we are)
The picket-fence longwave scheme (`PicketFenceLongwave`) gives a radiative-convective equilibrium that is **far too warm** vs RRTMG (and SOCRATES): in moist RCE, PF surface ≈ +13 K and climbing, water-vapour-feedback amplified. **Root cause is now definitively identified (exp #21): PF's correlated-k BAND STRUCTURE over-traps.** On an identical (T,p,q) profile, line-by-line truth (linepyline at PF's D=1.66) gives OLR = 246.8 (moist)/347.9 (dry); RRTMG matches it (242.6/340.6) but **PF emits far too little (211.0/330.8) — over-trapping ~36/~17 W/m²**. Since the LBL uses the SAME k-source and diffusivity as PF, the only difference is PF's 7-wide-bands/8-g-points vs full spectral resolution. NOT the k-tables (LBL-validated, exps #7–13), NOT the solver source function (exp #20, linear-in-τ inert), NOT diffusivity (D1), NOT vertical resolution (exp #19), NOT convection/feedback. The fix is a band-structure design decision (Task F), not a bug. Worst bands: 800–1250 window (+16 moist), 10–500 H₂O rotation (+9–11), 1250–1800 ν₂ (+4); CO₂ core is near-perfect.

### Settled facts — DO NOT re-investigate (already proven above)
1. **CO₂ k-tables are LBL-faithful.** The Chaverot source had a 54× 500–630 far-wing defect (exp #6); the linepyline rebuild `earth_low_res_lw_co2refined_linepyline.nc` fixes it (exp #7). Validated at cold (exp #7) and warm (exp #11) nodes; core ≈0.85–0.9× LBL (a benign ngpt=8 quadrature under-read, exp #8).
2. **H₂O axis (incl. MT_CKD continuum) is LBL-faithful** within the same ~0.85 ngpt=8 factor across X≤0.1 (exp #13). PF reads `specific_humidity`, converts to VMR correctly, interpolates the tabulated mixture k.
3. **The original "cold stratosphere bias" was mostly a missing-O₃ comparison artifact** (exp #10): RRTMG ran with climatological O₃; PF has none. Under matched O₃=0, the strat gap at 40–90 hPa is within ±3 K.
4. **The 500–630 CO₂ far-wing fix (the linepyline table) is radiatively inert for strat T** (exp #10) — the strat is set by the faithful 630–700 core. The wing fix is a legitimate k improvement but NOT the cold-strat fix that exp #6 claimed.
5. **The warm bias is table-independent** (linepyline ≈ shipped, both ≈296 K in moist RCE, exp #15) → it is in the solver, not the absorption data.
6. **Diffusivity factor in the kernel is correct (1.66).** Ruled out.

### THE OPEN PROBLEM (quantified, exp #16 + dry split)
PF over-traps by 51 W/m² OLR on an identical moist profile. The dry (`--dry`, CO₂-only) test splits it into TWO independent causes:
- **~38 W/m² (≈¾, moist only) — the 800–1800 cm⁻¹ "window" band is TOO WIDE.** It lumps the transparent 8–12 µm window with the opaque H₂O ν₂ band into one 8-g-point k-distribution; correlated-k then can't let the window photons escape. *This is the same "band too wide" pathology that started this whole plan, now in the window.*
- **~13 W/m² (dry, CO₂-only) — the isothermal-layer Planck source.** The kernel emits one B(T_mid) per layer; RRTMG uses a linear-in-τ source. Resolution-dependent if that's the cause.

### NEXT TASKS (priority order) — each is self-contained

**TASK A — DONE (exp #17).** Window split at 1250 cm⁻¹ → `earth_low_res_lw_co2refined_linepyline_winsplit.nc` (7 bands). Moist −51 → −32, dry −13 → −10. Partial.

**TASK A2 — DONE (exp #18).** Finer subdivision (`…,800,980,1180,1800,…`) made it WORSE (−34). Transparent windows want to be wide. Window subdivision exhausted; `winsplit2` rejected. The 1250 partition is the window deliverable.

**TASK B — DONE (exp #19).** Dry residual resolution-independent (−9.86 at NZ=40/80/160) → NOT the isothermal-layer Planck source → a fixed two-stream formulation difference.

**TASK D — DONE (exp #20).** Term-by-term vs rte-rrtmgp. D1 (diffusivity) null. D2 (Planck-fraction) non-issue (PF already uses it). D3 (linear-in-τ source) — the only formulation difference — implemented (`kernels.py` + `_interface_temperatures`/`lev_src` in `component.py`) and TESTED: radiatively inert (~0.4 W/m²) at NZ=10/20/40. The solver source function is NOT the over-trapping. **Decision pending: keep linear-in-τ (correct physics, 20 tests pass) or revert (simpler).**

**TASK E — DONE (exp #21).** LBL tiebreak on the identical profile (linepyline at D=1.66): LBL OLR agrees with RRTMG, PF over-traps ~36 (moist)/~17 (dry) W/m². Root cause = correlated-k **band structure** (same k/D/profile, only spectral binning differs). Localized to spectrally-inhomogeneous bands: 800–1250 window (+16.4 moist), 10–500 H₂O rotation (+9–11), 1250–1800 ν₂ (+4); CO₂ core near-perfect.

**TASK F (DESIGN DECISION — user-driven, NOT a bug hunt).** The over-trapping is the fundamental correlated-k limitation for few wide inhomogeneous bands. Three paths, pick per project goal:
- **(a) more spectrally-homogeneous bands.** Target the offenders: isolate a CLEAN transparent window (the 8–12 µm region) separate from BOTH its edges and the H₂O ν₂; give the wide 10–500 H₂O-rotation band its own refinement. KEY lesson from exp #18: separate transparent-from-opaque regimes; do NOT subdivide within a single regime (that hurt). Re-run `lbl_olr_tiebreak.py` + `lw_forward_pf_vs_rrtmg.py` after each candidate partition; success = LBL−PF per-band → 0.
- **(b) more g-points — REJECTED (exp #22).** ngpt 8→16→32 leaves PF OLR flat (≤0.3 W/m²); the over-trapping is band width, not g-resolution. Do not revisit.
- **(c) accept & document.** ~10–35 W/m² is the price of a fast 7-band scheme; for GCM use the circulation is forgiving and the scheme can be tuned (PLASIM lesson). Document the bias and move on.

**TASK C (cleanup, after F): re-validate RCE** with whatever table/scheme F settles on.
- The MOIST RCE script's tropopause diagnostic is still broken (cold-point returns model top, 2.40 hPa). Replace with a convective-heating-top or WMO lapse-rate criterion (the dry script already uses a convective-θ-top — see `tropopause()` in `rce_dry_pf_vs_rrtmg.py`; θ won't work for moist, use the convective heating tendency instead). Then re-run dry+moist RCE with the fixed solver/table and compare the PF-vs-RRTMG tropopause — the original project goal.

### KEY FILES & ARTIFACTS
Experiment scripts (all in `scripts/experiments/`, run with the climt env unless noted):
- `lw_forward_pf_vs_rrtmg.py` — **the key diagnostic.** Fixed-profile single-call RRTMG vs PF; `--dry`, `--table <name>` (exp #17), and `--nz <N>` (exp #19) flags. Outputs `debug_data/lw_forward_pf_vs_rrtmg_{moist,dry}_<table>.png`.
- `rce_dry_pf_vs_rrtmg.py` — dry RCE (SimplePhysics no-BL + DryConvectiveAdjustment). FIXED to be truly dry (`use_external_surface_specific_humidity=True`, surface q=0); 1 m slab; convective-θ-top tropopause; prints imbalance + q_max.
- `rce_moist_pf_vs_rrtmg.py` — moist RCE (SimplePhysics BL + EmanuelConvection). Works; tropopause diagnostic still broken (Task C).
- `pf_warm_node_fidelity.py`, `pf_h2o_axis_fidelity.py` (linepyline env), `pf_node_fidelity_vs_lbl.py`, `raw_chaverot_vs_lbl_coldnode.py`, etc. — k-fidelity validation (all PASSED; see exps #11–#13).
- Note `pressure_grid_log` in the .nc tables is **natural-log of Pa** (use `np.exp`, not `10**`).

Tables (`climt/_data/picket_fence/correlated_k/`):
- `earth_low_res_lw_co2refined_linepyline.nc` — the validated 6-band LW deliverable (drop-in for the defective Chaverot `..._co2refined_gl.nc`). Built with MT_CKD continuum ON.
- `earth_low_res_lw_co2refined_linepyline_winsplit.nc` — 7-band variant (window split at 1250 cm⁻¹, exp #17). Recovers ~16–20 W/m² of OLR over-trapping; **the window-band deliverable** (Task A2 confirmed finer splits don't help).
- `earth_low_res_lw_co2refined_linepyline_winsplit2.nc` — 8-band finer-window variant (exp #18). **REJECTED** (made over-trapping worse); safe to delete, kept only as a negative-result reference.
- `earth_low_res_lw_8gpt.nc` — shipped 4-band Chaverot table (current production).

### ENVIRONMENT & CONVENTIONS
- climt Python: `/Users/joymonteiro/miniconda3/envs/climt/bin/python` (or `conda run --no-capture-output -n climt python …` — the `--no-capture-output` keeps tqdm live).
- linepyline env (table generation + LBL reference, has linepyline + MT_CKD + netCDF4): `/Users/joymonteiro/miniconda3/envs/linepyline/bin/python`.
- `graphify` CLI: `/Users/joymonteiro/miniconda3/envs/climt/bin/graphify`. After code changes run `graphify update .` then `python scripts/augment_graph.py` (done this session through exp #13's additions; the RCE/forward scripts added afterward are NOT yet re-indexed — run it next session).
- Experiment scripts → `scripts/experiments/`; figures → `debug_data/`. Use the `climt` conda env. Log every experiment to this plan.

### LOOSE ENDS
- **Spawned background task:** "Add O₃ 9.6 µm LW absorber to PicketFence" (PF has no ozone; needed for realistic stratosphere — exp #10). Independent of the warm-bias work.
- The graph needs re-indexing for the scripts added after exp #13 (see above).
- Convergence: RCE runs need long spin-up; check `imbal = SW−OLR` (|·|≲1 W/m²) and `q_max` before trusting RCE numbers. With the +51 W/m² solver bias unfixed, PF RCE will not converge to a sensible state — fix Task A/B first.
- **Docs seed:** `docs/notes/corr-k-explainer-seed.md` captures a worked correlated-k explainer (how it works + how the 800–1800 window-band lumping fails) for a future corr-k doc page, with the `k(g)` figure regenerated by `scripts/experiments/plot_kg_curves.py` → `debug_data/kg_curves_window_band.png`. The Task A window split is the live demonstration of that failure.

---

## 2026-06-04 — Sub-project A: CO₂-adjustable 14-band table (`earth_hifi_lw`) productionized

Implemented per design `docs/superpowers/specs/2026-06-03-pf-co2-adjustable-hifi-earth-lw-design.md` and plan `docs/superpowers/plans/2026-06-04-pf-co2-adjustable-hifi-earth-lw.md`. Code Tasks 1–8 landed (writer/loader/`interpolate_k`/optics/component/sampler/generation-driver + consolidated `@njit(parallel=True)` LW kernel, all TDD, bit-for-bit kernel regression). Table generated by user (Task 9): `earth_hifi_lw.nc`, `k_coefficients` shape `(1,14,8,12,8,7,10)`, CO₂ axis `logspace(-5,-2,10)` = 10–10000 ppm.

**Generation fidelity (re-cross-checked, not "hours" — minutes on a fast box, full LBL, not binned):**
- k dynamic range 5e-26 … 8e3 (29 orders); real pressure broadening (k rises 0.01→1.4 from 1 Pa→1e4 Pa) and Boltzmann T-structure (ν₂ peaks ~110 K).
- **earth_hifi_lw band-mean CO₂-band k @464 ppm = 0.628 vs validated candD = 0.608 (ratio 1.03).** earth_hifi == candD to 3%.

**Single-column LBL fidelity (`eval_band_structure.py`, fixed exp #16 profile, 376 ppm, D=1.66):**

| table | moist PF | moist LBL−PF | dry PF | dry LBL−PF |
|---|---|---|---|---|
| candD (validated ref) | 235.29 | **+11.47** | 336.90 | **+10.99** |
| earth_hifi_lw | 235.26 | **+11.50** | 336.90 | **+10.99** |

→ earth_hifi_lw is fidelity-**identical** to candD (Δ ≤ 0.03 W/m²). CO₂ interp at off-node 376 ppm (between the 215/464 nodes) is essentially exact. The single-column LBL−PF ≈ +11 is the *fixed-profile* number for **both** tables (the docstring's <5 W/m² target was never met by candD either); candD is nonetheless the table that converged to the design's +2.7 K moist-RCE result, i.e. the RCE metric — not fixed-profile OLR — is what gates the ship.

**A5 — CO₂-interpolation scheme confirmed (`eval_co2_interp_accuracy.py`, leave-one-out node test, no linepyline):**
- Empirically CO₂ band-mean k ∝ amount^≈1 (×2.154 per ⅓-decade node). For a power law, log-k-vs-log-C interpolation is exact; linear-k overshoots between log-spaced nodes.
- Pooled over all 8 interior nodes (2-node-apart **stress test** = upper bound; production adjacent-node error ~4× smaller):
  - **log-k**: median |rel err| **0.06%**, p95 **5.5%**, signed bias ~0.
  - **linear-k**: median 0.08%, p95 **30.9%**; CO₂ ν₂ band-mean error **+19.6%** @4641 ppm vs log-k **+1.3%**.
- **VERDICT: keep `_CO2_INTERP_LOGK = True` (geometric).** Note the design's stated *motivation* ("convex/saturating") was slightly off — k is ~linear, not convex — but log-k is exact for any power law, so the choice stands.
- *Correction to log:* the spec's code-changes line 58 ("linear-in-k vs log-X_CO₂") contradicts the Design section (line 46, "log-k"); implementation followed the Design section. A5 confirms it.

**New/changed validation tooling (committed):** `--co2`/`--table` flags on `rce_moist/dry_pf_vs_rrtmg.py`; `co2=` override on `eval_band_structure.py`; new `eval_co2_interp_accuracy.py`.

**STILL OPEN (needs user's linepyline env / long compute — Task 10c):**
- Off-node LBL OLR at 400/2000 ppm (regenerate `lbl_olr_spec_{moist,dry}.npz` at those CO₂ via the D=1.66 LBL path; the PF-side hook is in `eval_co2_interp_accuracy.py` Part 2).
- Converged RCE (moist+dry) at 280 and 1120 ppm vs RRTMG + the *true* ngpt=2 baseline (`--table earth_low_res_lw`), via the new `--co2`/`--table` flags. 2-day moist smoke with `--table earth_hifi_lw --co2 280` ran clean (576 steps, ~33 s).
- Then Task 11: flip default `earth_hifi_lw.nc` → `earth_low_res_lw.nc`, shipped-table integration test, MANIFEST, graphify re-index.
