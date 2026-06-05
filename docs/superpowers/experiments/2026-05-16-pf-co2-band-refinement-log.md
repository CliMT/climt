# Experiment Log: Picket-Fence LW CO2 Band Refinement

**Date:** 2026-05-16
**Branch:** `feature/picket-fence-radiation`
**Related plan:** `docs/superpowers/plans/2026-05-16-pf-co2-band-refinement.md`

## Goal

Reduce the Earth-LW picket-fence cooling-to-space bias around the 5–10 hPa
"kink" by refining the 500–800 cm⁻¹ band (CO2 15 μm region) where one
correlated-k distribution can't represent both line cores and wings.

Reference truth was originally RRTMG, then upgraded to **line-by-line via
linepyline** (HITRAN2024 + MT_CKD continuum) after we recognized RRTMG is
itself band-averaged and may not be a clean ground truth for line-core
fidelity questions.

## Diagnostic Infrastructure (built during the session)

| Script | Env | Purpose |
|---|---|---|
| `examples/cooling_to_space_gsplit_sweep.py` | climt | Heating-rate vs pressure, isothermal 250 K column, multiple PF tables vs RRTMG |
| `examples/lbl_band_mean_k.py --profile {iso,adiabat}` | linepyline | Run LBL on a column, save `kappa(p, nu)` to netCDF |
| `examples/compare_lbl_pf_band_k.py` | climt | Single-pressure tabular comparison: PF vs LBL Planck-weighted band-mean k |
| `examples/compare_lbl_pf_heatmap.py` | climt | (band × pressure) heatmap of PF/LBL ratio, two profiles × all PF tables |

LBL was Planck-weighted band-averaged using each PF table's own band edges, so
the comparison is apples-to-apples per table.

## Tables Generated

All from `scripts/generate_picket_fence_tables.py` against the Chaverot Earth
k-table (376 ppm CO2, R=500).

| Filename | Bands | Edges (cm⁻¹) | ngpt | Quadrature |
|---|---|---|---|---|
| `earth_low_res_lw_8gpt.nc` (shipped) | 4 | 10,500,800,1800,3250 | 8 | Gauss-Legendre |
| `earth_low_res_lw_8gpt_2s097.nc` | 4 | 10,500,800,1800,3250 | 8 | two-stretch g=0.97 |
| `earth_low_res_lw_co2refined.nc` | 6 | 10,500,**630,700**,800,1800,3250 | 8 | two-stretch g=0.97 |
| `earth_low_res_lw_co2refined_gl.nc` | 6 | 10,500,630,700,800,1800,3250 | 8 | Gauss-Legendre |
| `earth_low_res_lw_co2refined_16gpt_2s097.nc` (deleted) | 6 | 10,500,630,700,800,1800,3250 | **16** | two-stretch g=0.97 |
| `earth_low_res_lw_10band_gl.nc` | 10 | 10,350,500,630,700,820,980,1180,1800,**2380**,3250 | 8 | Gauss-Legendre |
| `earth_low_res_lw_10band_2s097.nc` | 10 | same as 10band_gl | 8 | two-stretch g=0.97 |
| `earth_low_res_lw_11band_tightcore_2s097.nc` | 10 | 10,500,**600,660,720**,820,980,1180,1800,2380,3250 | 8 | two-stretch g=0.97 |

## Experiments and Findings

### Experiment 1 — Band partition (4 → 6) at 630/700 cm⁻¹

**Hypothesis:** Splitting 500–800 cm⁻¹ at 630/700 isolates the CO2 Q/R-branch
core (630–700) from the P-wing (500–630) and the long-λ window edge
(700–800), recovering line-core opacity that was averaged away.

**Method:** Generate `co2refined.nc` (6 bands, 2s, 8 g-pt), compare to shipped
4-band at p=10 hPa, T=250 K, q=0.

**Result:** PF/LBL ratio in the 630–700 line-core band = **0.450** (with 2s
quadrature). The 6-band table reduces the cooling-to-space kink visibly but
does not eliminate it.

**Interpretation:** Band split helps relative to the conflated 500–800 band,
but the residual deficit means a single correlated-k distribution within the
60–70 cm⁻¹ band still can't represent both the Lorentzian core and the
surrounding wings simultaneously.

---

### Experiment 2 — Quadrature: Gauss-Legendre vs two-stretch

**Hypothesis:** Two-stretch quadrature (g_split=0.97, 4+4 nodes) places more
g-points in the strong-line tail where the core lives, better resolving the
within-band k-distribution.

**Result (PF/LBL at p=10 hPa, iso 250 K):**

| Band (cm⁻¹) | 6-band GL | 6-band 2s | 10-band GL | 10-band 2s |
|---|---|---|---|---|
| 500–630   | 0.161 | 0.527 | 0.161 | 0.527 |
| 630–700   | 0.230 | 0.450 | 0.230 | 0.450 |
| 700–820   | 0.146 | 0.385 | 0.126 | 0.354 |
| 820–980   | —     | —     | 0.207 | 0.474 |
| 980–1180  | —     | —     | 0.130 | 0.408 |
| 1800–2380 | —     | —     | 0.138 | 0.650 |

**Interpretation:** Two-stretch is a **universal win**, roughly doubling the
recovered line-core opacity across every CO2-region band. Confirmed in the
heatmap: every iso 2s panel is consistently lighter blue than its GL
counterpart. **GL tables are now obsolete for Earth LW.**

The 1800–2380 cm⁻¹ (4.3 μm) gain is especially large (0.138 → 0.650, ~4.7×)
because that band has the same line-core-plus-wing structure as 15 μm, and
the same fix applies.

---

### Experiment 3 — Decoupling band count from line-core fidelity

**Hypothesis:** Maybe further outer-band refinement (10 bands vs 6 bands)
also helps the 15 μm core, by reducing g-point "competition" from unrelated
spectral features.

**Result:** 6-band 2s and 10-band 2s give **identical** 630–700 ratio (0.450).
GL versions also identical (0.230).

**Interpretation:** Line-core fidelity in the 630–700 band is fully decoupled
from outer-band refinement. The deficit lives entirely inside that one band's
k-distribution.

---

### Experiment 4 — More g-points (ngpt=8 → 16) in the refined 6-band table

**Hypothesis:** Doubling g-points lets the within-band CDF sample the strong-
absorption tail more finely.

**Result:** 630–700 ratio: 8gpt = **0.450** → 16gpt = **0.467**. A 3.8%
relative gain for 2× the table size and 2× the runtime cost.

**Interpretation:** The k(g) CDF is already well-sampled at 8 g-points; the
problem is not how finely we sample the CDF, it's what the CDF *contains*.
Table deleted as not worth shipping.

---

### Experiment 5 — Tighter line-core band (660–720 cm⁻¹)

**Hypothesis:** Even narrower band isolation (60 cm⁻¹ wide, centered on the
Q-branch) might finally separate the line core from any residual wing.

**Result:** PF/LBL ratio in 660–720 = **0.446**. Essentially identical to the
70 cm⁻¹-wide 630–700 result.

**Interpretation:** Decisive. Neither narrower bands nor more g-points moves
the line core past ~0.45. Both knobs are saturated.

---

### Cross-cutting Finding: The 0.45 Ceiling Is Source-Data-Limited

The factor-of-2 line-core deficit is **intrinsic to the Chaverot R=500
k-table** used as input. At ν=667 cm⁻¹ the source-table bin width is
~1.3 cm⁻¹; real CO2 rotational line cores at stratospheric (T=250 K,
p=10 hPa) conditions are pressure-broadened to ~0.01–0.05 cm⁻¹ wide and
spaced ~1.5 cm⁻¹ apart. A 1.3 cm⁻¹ bin averages one peak plus neighboring
wings into a single number *before* the correlated-k CDF is computed. The
PF tables can't recover what's been discarded upstream.

**To go below 0.45 you need a higher-R input k-table (R≥2000)** or direct
LBL-derived k-distributions. This is a separate, much larger project.

---

### Experiment 6 — Heatmap diagnostic: full (band × pressure) view, two profiles

**Method:** Repeat the band-mean comparison across the full p-grid for two
T(p) profiles (isothermal 250 K and dry adiabat T_s=288 K, T_floor=80 K).

**Findings from the heatmap:**

1. **Pressure broadening recovery is visible.** In every iso panel the CO2
   bands (565, 665, 760, …) get lighter toward higher p. At 10 hPa the 665
   cm⁻¹ cell is ~0.4–0.5; at the top of the visible range (~0.001 hPa) it's
   ~0.05. Textbook: pressure broadening fills in line cores so PF's smeared
   k becomes more faithful at higher p.

2. **The dry adiabat exposes cold-T table behavior.** Several bands go *deep
   red* (PF over-opaque by 10×+) at low p on the adiabat:
   - 1300/1490 cm⁻¹ window across all tables
   - 2090 cm⁻¹ (4.3 μm core) in 10-band variants and tightcore
   - 770 wing band on tightcore

   Both PF and LBL are evaluated at the same T(p), so this is the **PF
   table's k-coefficient grid behaving badly at very cold T (~100 K)**. The
   Chaverot source was sampled at warm/temperate-Earth T; its extrapolation
   to T<150 K is unreliable.

3. **Tightcore (11-band) does not beat 10-band 2s.** The 550/630/690/770
   cells are all dark blue with the same saturation. Once again: narrower
   bands don't recover what was lost upstream.

4. **Edge artifacts:** first/last bands (180 cm⁻¹, 2815 cm⁻¹) solid red
   across all panels — denominator-near-zero (LBL k ≈ 0 there). Ignore.

---

## Strengths and Weaknesses of Each Knob

| Knob | Strength | Weakness | Verdict |
|---|---|---|---|
| Band split at 630/700 (CO2 core/wing) | Doubles core opacity vs un-split 4-band; cheap | Saturates at 0.45 ratio | **Keep** as minimum |
| Two-stretch quadrature g=0.97 | Universal ~2× gain across all CO2 bands; same g-point count | None observed; strictly dominates GL | **Always on** |
| Outer-band refinement (4→6→10) | Captures 4.3 μm and window structure; lets each region have its own k-distribution | No effect on the 15 μm core itself | **Use 10-band** for completeness |
| More g-points (8→16) | Resolves CDF more finely | <4% gain at 2× cost; source data is the bottleneck | **Skip** |
| Narrower core band (60 cm⁻¹ vs 70) | None measurable | Same saturation as wider | **Skip** |
| Higher-R source table | Only path past 0.45 ceiling | Major upstream project (regenerate from HITRAN) | **Future work** |

## Reversal: Full-Column Heatmap Refutes the Initial 2s Recommendation

After the LBL output bug was fixed (linepyline returns p in hPa with doubled
NaN-filled levels; we now save kappa on our own explicit Pa grid), the
heatmap revealed the **troposphere**, which had been clipped out of all
earlier diagnostics. The new picture overturns the single-pressure (p=10 hPa)
recommendation.

**What the full-column iso row shows:**

| Table | Stratosphere (p ≲ 10 hPa) | Troposphere (p ≳ 100 hPa) | Whole column |
|---|---|---|---|
| `co2refined_gl` (6-band GL) | light-blue (under by ~3–5×) | pale (near 1.0) | **near-white almost everywhere** |
| `co2refined` (6-band 2s g=0.97) | pale to light-blue (under by ~2×) | **deep red (over by 10×+)** in 565, 750, 1300, 2525 | strong errors of opposite sign at top vs bottom |
| `10band_2s097` | similar to 6-band 2s | same tropospheric red, plus saturated 1490/2090 columns | strong errors at top and bottom |
| `10band_gl` | similar to 6-band GL | pale | balanced |

**Why this was missed earlier.** All single-pressure comparisons were
anchored at p=10 hPa because that's where the cool-to-space kink sits. At
10 hPa, 2s gives a line-core ratio of 0.45 vs GL's 0.23, so 2s "won." But
that strat-only win comes with a large **tropospheric over-opacity** cost
that the single-p table never showed.

**Physical interpretation.** g_split=0.97 places 4 of 8 g-points in the
upper 3% of the k-CDF — the strong-line tail. That's the right thing for
sharp Lorentzian cores at low p, where most band-mean opacity actually lives
in the tail. In the troposphere, pressure broadening flattens the lines:
the "strong tail" the 2s table was placing g-points in no longer carries
that much opacity, but the table still weights it heavily and reports
opacity that LBL says isn't there. GL's uniform g-sampling has no such
asymmetric bias.

**Corrected verdict per quadrature:**

| Regime | Best quadrature | Reason |
|---|---|---|
| Stratosphere (low p, sharp cores) | two-stretch | Resolves narrow lines in g-space |
| Troposphere (high p, broadened lines) | Gauss-Legendre | No asymmetric tail bias |
| Full column (Earth-LW default) | **Gauss-Legendre** | Fewer extreme errors at any altitude |

## Revised Shipping Recommendation

**`earth_low_res_lw_co2refined_gl.nc` is now the leading candidate for
shipping**, not `10band_2s097` as previously stated. The full-column
heatmap shows it is the most balanced table: near-perfect agreement
through the troposphere and lower stratosphere, with the line-core deficit
confined to the upper stratosphere/mesosphere (where any picket-fence
representation will struggle).

**The previous "shipping table" `earth_low_res_lw_10band_2s097.nc` is NOT
recommended.** Its low-p line-core gain is real, but it introduces large
(10×+) tropospheric over-opacity that will bias OLR and tropospheric HR.

**Important caveats this experiment did not resolve:**

1. **15 μm line-core opacity remains ~75–80% under-represented at low p** in
   the GL table (line-core ratio 0.23 at p=10 hPa). This is the source-data
   R=500 ceiling and persists across all PF configurations. Residual cool-
   to-space bias above ~1 hPa is unavoidable with the current input.
2. **Cold-T (T<150 K) behavior is unreliable across all tables.** Both GL
   and 2s tables on the dry-adiabat row show extreme PF over-opacity
   throughout the cold stratosphere — likely a Chaverot-table T-grid
   extrapolation issue, not a quadrature or band-partition issue. RCE
   diagnostics involving very cold stratospheric T should be validated
   against LBL before trusting upper-atmosphere fluxes.
3. **`g_split` was not re-tuned** on the refined band partition. The
   original g_split=0.97 was selected on the 4-band layout for stratospheric
   line-core fidelity. An intermediate value (0.93–0.95) on the refined
   6-band partition might capture most of the stratospheric 2s gain without
   the tropospheric overshoot — left for the follow-up session.

## Cleanup Recommendations

**Updated after the reversal** — most of the earlier "dead" calls were based
on the single-pressure (10 hPa) view, which we now know was misleading.
Re-evaluate before deletion:

- `earth_low_res_lw_co2refined_gl.nc` — **KEEP, now the leading shipping
  candidate.** Earlier called dead because beaten by 2s at p=10 hPa; full-
  column heatmap shows it is the most balanced overall.
- `earth_low_res_lw_co2refined.nc` (6-band 2s) — keep for now as the
  reference point in the 2s-vs-GL comparison; revisit after the follow-up
  g_split sweep.
- `earth_low_res_lw_10band_2s097.nc` — keep; tropospheric overshoot is a
  serious problem but the 4.3 μm split is the right structural choice. A
  10-band-GL or 10-band-2s with lower g_split would likely supersede it.
- `earth_low_res_lw_10band_gl.nc` — keep; counterpart to 10band_2s for the
  full-column comparison.
- `earth_low_res_lw_11band_tightcore_2s097.nc` — still no improvement over
  10band_2s; safe to delete unless used as a negative-result reference.
- `earth_low_res_lw_8gpt_2s090/095/099.nc`, `earth_low_res_lw_8gpt_twostretch.nc`
  — g_split-sweep leftovers on the 4-band layout. The conclusion that
  "g=0.97 is the right setting" was made on the 4-band cool-to-space
  diagnostic; the heatmap reversal suggests g_split itself needs re-tuning
  on the refined partitions. Keep until that sweep is re-run.

`earth_low_res_lw_16gpt_2s097.nc` was already deleted in-session.

## Future Work (in priority, revised after the reversal)

1. **Re-run cool-to-space and RCE diagnostics with `co2refined_gl` as the
   primary PF curve.** Test whether the more balanced full-column opacity
   gives a flatter HR profile despite the larger stratospheric line-core
   deficit, and whether the original cool-to-space "kink" actually returns
   or whether it was an artifact of the 4-band layout's even-worse line-
   core dilution.
2. **Re-tune g_split on the refined band partitions.** Try 0.90, 0.93, 0.95
   on the 6-band layout. Plot the full-column heatmap for each. The
   hypothesis is that a lower g_split keeps most of the stratospheric line-
   core gain without the tropospheric overshoot. If a clear winner emerges,
   it may displace `co2refined_gl` as the shipping recommendation.
3. **Investigate the cold-T behavior.** Check Chaverot table T-grid bounds
   (it goes 50–1000 K, so probably not a bounds issue) and whether the k(T)
   interpolation/extrapolation is the bug, or whether the source data
   itself is unreliable at very cold T. The dry-adiabat row of the heatmap
   is the diagnostic to keep watching.
4. **(Long-term) Build a higher-R source k-table** from HITRAN directly,
   with R≥2000 in the CO2 15 μm region, to break the 0.45 line-core
   ceiling. Independent of the quadrature/band-partition question, but
   sets the absolute upper bound on what any picket-fence representation
   can achieve.

## Diagnostic Bug Fixed Late in Session

The linepyline output's `p` coordinate is in **hPa** (despite the docs
saying "must be in Pa" — input is Pa, output is hPa), and for the 60-level
input it returned 120 levels by interleaving the original Pa values as
NaN-filled rows. The first heatmap pass dropped the NaN rows, which
inadvertently kept only the levels above ~10 hPa — the troposphere was
silently clipped from the visualization. After the fix (`lbl_band_mean_k.py`
now writes kappa on its own explicit Pa grid, asserting one finite
row per input level), the full column is included and the reversal above
became visible. **Lesson: always sanity-check the displayed range of a
plotted axis against the requested data range — the original heatmap
showed y from 10⁻³ to 10¹ hPa even though the input p-grid was 0.1 to
1013 hPa, and this was missed for one round of interpretation.**
