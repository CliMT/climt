# PicketFence Longwave: Ozone (O₃) Absorption — Design

**Date:** 2026-06-01
**Status:** Design — awaiting review before implementation plan
**Author:** scoping session (climt `feature/picket-fence-radiation`)

## Goal

Give `PicketFenceLongwave` a longwave ozone absorber so it can reproduce a
realistic stratosphere. Experiment #10 of the CO₂ band-refinement plan
(`docs/superpowers/plans/2026-05-16-pf-co2-band-refinement.md`) showed that
RRTMG's climatological ozone profile (max ~3.9 ppm in the stratosphere)
accounts for **~30 K of stratospheric warming** in a dry CO₂-only radiative
column. PF currently has **no ozone absorber at all**, so every prior
PF-vs-RRTMG head-to-head was RRTMG-with-O₃ vs PF-without. The dominant LW
feature is the O₃ **9.6 µm band (~980–1080 cm⁻¹)**, which lies inside PF's
existing 800–1800 cm⁻¹ window band.

## Decision: O₃ as a dedicated additive gas (option 1)

O₃ is added as a **second gas** in the correlated-k table, alongside the
existing premixed air gas, with its own precomputed cross-section
`k_O3(T, p)` and a **runtime-variable O₃ amount**. Optical depth is linear in
amount: `τ_O3 = k_O3(T,p) · m_O3`, combined with the air optical depth by
**additive overlap**.

### Why this option (vs the alternatives considered)

| Option | What varies at runtime | k-table stores | Verdict |
|---|---|---|---|
| **1 — additive gas (chosen)** | O₃ amount (profile) | `k_O3(T,p)`, fixed | RRTMG-faithful, cheap, correct stratospheric peak |
| 2 — fixed bake-in | nothing (frozen) | premixed `k(T,p,X_H2O)` with O₃ mixed in at one VMR | Rejected: bakes a constant O₃ *mass* fraction with height (wrong shape — O₃ peaks in the stratosphere where air is thin); cannot vary at runtime |
| 3 — O₃ VMR axis | O₃ amount (profile) | `k_O3(T,p,X_O3)` — 7-D, quadrilinear | Rejected: models an O₃ *self*-concentration dependence of the cross-section that is physically negligible in the 9.6 µm band; high effort/risk, no RRTMG analog |

**RRTMG alignment.** RRTMG-LW is itself a correlated-k model. It ingests the
ozone profile as a runtime field and computes `k(T,p) · O₃_amount` (linear in
amount) — exactly option 1. RRTMG's "key-species η-binning" refines *overlap*
between two co-located major absorbers (e.g. H₂O:O₃ in the 9.6 µm band); it is
**not** an O₃-self-VMR axis (so it does not correspond to option 3). PF's
plain additive overlap is a coarser treatment of that overlap, consistent with
PF's existing coarse-band philosophy and adequate in a window band where CO₂ is
weak and O₃ is optically thin.

## Architecture fit

The PF correlated-k stack already has every primitive this needs; the work is
wiring them together in a combination that has **not been exercised before**.

### What already exists (no change)

- **`interpolate_k`** (`climt/_components/picket_fence/optics/correlated_k.py:96`)
  already loops over `ngas` inside the 6-D (H₂O-X-axis) branch. If `k_O3` is
  stored **broadcast (constant) across the `nX` H₂O nodes**, the X-interpolation
  returns `k_O3` regardless of H₂O VMR — correct and requires **no change**.
- **Additive overlap** (`_compute_ck_optical_depth_additive`) already sums
  `k_interp[ig] · gas_amounts[ig]` over gases on a shared g-point grid. Units
  work out: `gas0` is m²/kg-air × air-column-mass; `gas1` is m²/kg-O₃ ×
  O₃-column-mass. **No change.**
- **The netcdf writer** (`scripts/pf_table_builder/netcdf_writer.py`) already
  accepts `gas_names=(...)` and `overlap_method=...`. **No change.**
- **The ozone profile** — climt already ships `init_ozone`
  (`climt/_core/initialization.py:1070`, loading `ozone_profile.npy`) registered
  as `mole_fraction_of_ozone_in_air` (mol/mol, mid-levels). This is the **same
  climatological field RRTMG consumes**, so no new data is required and
  validation is apples-to-apples.
- **`sample_kappa_grid`** already supports O₃ as a HITRAN absorber.

### What changes

1. **Table builder** (`scripts/pf_table_builder/kappa_sampling.py` +
   `scripts/generate_pf_tables_linepyline.py`). Today non-H₂O absorbers are
   *premixed* into the air kappa (multiplied by mass fraction → m²/kg-air). For
   option 1 we need O₃'s cross-section as a **separate gas, per kg of O₃**
   (i.e. the raw `get_kappa_hitran` output, *not* scaled by mass fraction). The
   builder gains a notion of a "secondary additive gas":
   - Sample `k_O3_raw(nT, nP, nNu)` (no mass-fraction scaling, no X-axis).
   - Build its own k-distribution with the same `ngpt`/quadrature →
     `(nband, ngpt, nT, nP)`, then broadcast across `nX` → `(nband, ngpt, nT, nP, nX)`.
   - Stack with the premixed-air gas → `k_coefficients` shape
     `(2, nband, ngpt, nT, nP, nX)`.
   - Write `gas_names=("effective", "o3")`, `overlap_method="additive"`,
     `background_is_premixed="true"`.
   - g-point weights are quadrature-defined (identical per gas), so additive
     overlap's shared-weight assumption holds.

2. **Component** (`climt/_components/picket_fence/lw/component.py`). The
   `_premixed_bg` path currently assumes a single gas and fills only
   `gas_amounts[0]` (air mass, H₂O as the live X-axis). Extend it to also carry
   secondary additive gases:
   - In `input_properties`, when `"o3"` is in `gas_names`, declare
     `mole_fraction_of_ozone_in_air` (units `mole/mole`, dims `["mid_levels","*"]`,
     alias `o3`) **in addition to** `specific_humidity`.
   - In `array_call`, keep `gas_amounts[0]` = air column mass and `h2o_vmr` from
     specific humidity (unchanged), and additionally set
     `gas_amounts[i_o3]` = column mass of O₃ from `mole_fraction_of_ozone_in_air`
     (convert mole fraction → mass via `M_O3 / M_dry_air`, then
     `compute_column_amount`).
   - Overlap stays additive (the existing additive path already handles
     `ngas > 1`).

3. **New shipped table.** Regenerate the Earth LW table with the O₃ gas, e.g.
   `climt/_data/picket_fence/correlated_k/earth_low_res_lw_co2refined_o3.nc`
   (6-band edges `10,500,630,700,800,1800,3250`, `ngpt=8`, two-stretch,
   `g_split=0.97`, plus the O₃ secondary gas). Generated in the
   `radiation`/`linepyline` env per the plan's conventions. Existing tables are
   left untouched — additive, delete-only rollback.

### Data flow (per column, per layer)

```
mole_fraction_of_ozone_in_air ──► m_O3 (kg/m²)  ─┐
                                                 ├─► τ = k_air·m_air + k_O3·m_O3
specific_humidity ──► X_H2O ──► k_air(T,p,X_H2O) ┘     (additive, shared g-grid)
air mass ──► m_air (kg/m²) ─────────────────────┘
k_O3(T,p) [X-broadcast] ─────────────────────────► same g-points
                                                 ──► lw_transport → fluxes, heating
```

## Testing

- **Optics unit test (new combination).** No existing test combines the 6-D
  H₂O-X-axis with `ngas=2`. Add a test: build a tiny 6-D, 2-gas additive table
  (gas1 = X-broadcast constant); assert `interpolate_k` returns gas1's value
  independent of `h2o_vmr`, and that `_compute_ck_optical_depth_additive` sums
  the two contributions correctly.
- **Component test.** With the O₃ table loaded: (a) zero O₃ ⇒ identical fluxes
  to the no-O₃ table within tolerance; (b) nonzero O₃ adds optical depth only in
  the 800–1800 cm⁻¹ window band and increases stratospheric LW heating.
- **Builder test.** O₃ raw kappa is per-kg-O₃ (not mass-fraction-scaled);
  resulting `k_coefficients` has shape `(2, nband, ngpt, nT, nP, nX)` and gas1 is
  constant along the X axis.

## Validation (acceptance)

Extend `scripts/experiments/radiative_equilibrium_pf_convergence.py` to run
**PF-with-O₃ vs RRTMG-with-O₃**, both using climt's `init_ozone` profile (same
field). Acceptance: PF recovers most of the ~30 K stratospheric warming O₃
contributes in RRTMG (90 hPa, 40 hPa, 10 hPa nodes from exp #10), and the
matched-input strat gap stays within a few K — i.e. the missing-O₃ artifact
identified in exp #10 is closed. Figures to `debug_data/`.

## Effort & risk

**Effort: ~3–4 focused days.**

| Task | Est. |
|---|---|
| 1. Builder: O₃ as separate per-kg gas + stack + gas_names | 0.5–1 d |
| 2. Component: ingest `mole_fraction_of_ozone_in_air`, fill `gas_amounts[i_o3]` | 0.5 d |
| 3. Generate shipped O₃ table (radiation/linepyline env) + verify dims | 0.5 d |
| 4. Unit + component tests (6-D × ngas=2; O₃ on/off) | 0.5 d |
| 5. RRTMG RCE validation + figures | 0.5–1 d |
| 6. `graphify update .` + `python scripts/augment_graph.py` | minutes |

**Risk: medium.** Chief risks and mitigations:
- *Never-combined 6-D X-axis × `ngas=2` path.* Mitigated: `interpolate_k`'s 6-D
  branch and the additive-overlap loop already iterate `ngas`; the new optics
  unit test pins this exact combination.
- *Builder unit correctness* (per-kg-O₃ vs per-kg-air). Mitigated by the builder
  test and the O₃-zero ≡ no-O₃-table component check.
- *Additive-overlap g-point correlation approximation* (vs RRTMG's η-binning).
  Accepted as second-order in a window band; bounded by the RRTMG RCE validation.
- *Profile source.* De-risked — climt already ships the climatological profile
  RRTMG uses; no new data.

## Out of scope

- Shortwave O₃ (Hartley/Huggins UV bands) — separate feature.
- η-style key-species overlap binning — additive overlap is sufficient here.
- Replacing/retiring existing CO₂-only Earth tables — this ships alongside them.
```

