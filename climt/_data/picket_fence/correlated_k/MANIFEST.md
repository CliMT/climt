# Shipped correlated-k tables

All `.nc` tables below are derived from **G. Chaverot's Zenodo dataset 16795590**
(CC-BY 4.0, https://zenodo.org/records/16795590), re-binned to picket-fence
resolution by `scripts/generate_picket_fence_tables.py`.

## LW band structures

Two LW band structures are shipped:

- **Option B (CO2-isolating):** `10, 500, 800, 1800, 3250` cm⁻¹ — splits the
  CO₂ 15 µm band (500–800) from the atmospheric window (800–1800). Best for
  present-day climate sensitivity and CO₂ forcing studies.
- **Option D (hot-moist-optimised):** `10, 500, 1200, 2000, 3250` cm⁻¹ —
  resolves H₂O near-IR (2000–3250) and window closure. Use for runaway-
  greenhouse / high-T moist-atmosphere pedagogy. **Limitation:** molar mass is
  frozen at dry-air Earth (0.028964 kg/mol); at H₂O mole fractions above ~20%
  the bulk molar mass drops and optical depth becomes inaccurate (~20% error).

SW bands are unchanged: `3250, 8000, 14000, 30000` cm⁻¹.

## File listing

| File | Planet | Kind | Bands (cm⁻¹) | g-pts | sha256 |
|---|---|---|---|---|---|
| earth_hifi_lw.nc | Earth (CO₂ 10–10000 ppm runtime axis) | LW candD-class (14-band, decoupled+log-X H₂O continuum, runtime CO₂) | 10…3250 (14-band) | 8 | `3192c45f2c4bad5243023177841089568270f246280bc5bbbe0dc74bac74b760` |
| earth_low_res_lw.nc | Earth (CO₂ 10–10000 ppm runtime axis) | **DEFAULT** LW candD-class (14-band, decoupled+log-X H₂O continuum, runtime CO₂; byte-copy of earth_hifi_lw) | 10…3250 (14-band) | 8 | `3192c45f2c4bad5243023177841089568270f246280bc5bbbe0dc74bac74b760` |
| earth_low_res_lw_4band_ngpt2_before.nc | Earth (376 ppm CO₂) | LW option B — *"before" reference* (the pre-2026-06 4-band/ngpt2 default, kept for sub-project B) | 10, 500, 800, 1800, 3250 | 2 | `94570a87d18c0c1f022d633433a663a0940816dc34c25693eb9ffd4aed1bfb82` |
| earth_low_res_lw_4gpt.nc | Earth (376 ppm CO₂) | LW option B | 10, 500, 800, 1800, 3250 | 4 | `8c9f83e2b38cfcd577cbc67f02d750afca80d0e2c2b70443d20825b4ccd04fe7` |
| earth_low_res_lw_8gpt.nc | Earth (376 ppm CO₂) | LW option B | 10, 500, 800, 1800, 3250 | 8 | `15b06a4311eae5bc0fecd255c54e7a810c80d90634edd2a5d27ebe7c2ef49136` |
| earth_low_res_lw_hotmoist.nc | Earth (376 ppm CO₂) | LW option D | 10, 500, 1200, 2000, 3250 | 2 | `9ae21cd2edf1c54bf64c7896673f1e5f93838fb3aa9a870c3f29eaa6b54aaccd` |
| earth_low_res_sw.nc | Earth (376 ppm CO₂) | SW | 3250, 8000, 14000, 30000 | 2 | `8709799c08dfaaf015a81a931c1eec699649559803b8fd2831eeed68d6360bda` |
| mars_lw.nc | Mars (CO₂+N₂) | LW option B | 10, 500, 800, 1800, 3250 | 2 | `94e61290e54be8f3552b7762014bc09443bd5e7f0d94d39932c726c4198d78c7` |
| mars_sw.nc | Mars (CO₂+N₂) | SW | 3250, 8000, 14000, 30000 | 2 | `8e186ac238e1bf7f562348f3fc70bbc006358efeb4f0d26c4aab7756a73157bb` |
| venus_lw.nc | Venus (H₂O+CO₂) | LW option B | 10, 500, 800, 1800, 3250 | 2 | `f308088c1cc2d6f028e44b1be3a004d4c01000d6347a4bc5bfff21dd285fe6b8` |
| venus_sw.nc | Venus (H₂O+CO₂) | SW | 3250, 8000, 14000, 30000 | 2 | `606960889f56af3b7c49e93c728be04c1a80b5555df3643cff3af94c2adbf83f` |

## Format

These tables carry an `h2o_vmr_grid` axis — H2O is a **runtime variable**
driven from `specific_humidity` in state. The non-H2O background
(CO2, N2, continua, CIA) is pre-mixed into the k-coefficients at the
composition Chaverot built them for. A netCDF attribute
`background_is_premixed = "true"` flags this convention so the component
multiplies k by column mass of AIR (not of H2O) and uses trilinear
(T, logP, logX_H2O) interpolation.

**Supported experiments:** realistic moisture variability, standard
climate runs with fixed CO2 budget.

**Not supported:** CO2-doubling / gas-forcing experiments. CO2-CO2 CIA and
line absorption are frozen at Chaverot's fixed CO2 fraction; scaling a
CO2 channel linearly does not reproduce the ρ² CIA dependence. For that
use case, generate per-gas tables via the `linepyline` workflow **with
continuum and CIA contributions included**, and add a runtime CIA handler
to the component.

## Regenerating

```sh
# Requires exo_k; use the `radiation` conda env.
python scripts/regenerate_shipped_picket_fence_tables.py \
    --chaverot-dir /path/to/correlated-k_tables
shasum -a 256 climt/_data/picket_fence/correlated_k/*.nc
```

The `.npz` files in this directory are pre-existing unit-test fixtures and
are not subject to the Chaverot provenance above. `single_band_unit_lw.nc`
is a synthetic constant-k fixture for the grey-limit unit test.

## linepyline-derived tables

Built offline by `scripts/generate_pf_tables_linepyline.py` from
HITRAN 2024 line data (via the [linepyline](https://github.com/exoclime/linepyline)
package) plus MT_CKD 4.3 (H2O continuum) and HITRAN CIA 2021/2024 (Titan
N2-N2, N2-CH4).

| File | Atmosphere | LW bands (cm⁻¹) | SW bands (cm⁻¹) | g-pts | H2O axis | Notes |
|---|---|---|---|---|---|---|
| trappist1e_hab1_lw.nc | THAI Hab1: N2 + 400 ppm CO2 + H2O | 10, 500, 1250, 2500, 3250 | — | 2 | yes | Covers THAI Ben1 (radiation identical). |
| trappist1e_hab1_sw.nc | as above | — | 3250, 8000, 14000, 30000 | 2 | yes | Stellar source from `trappist1.npz`. |
| trappist1e_hab2_lw.nc | THAI Hab2: pure 1-bar CO2 | 10, 500, 1250, 2500, 3250 | — | 2 | no | Covers THAI Ben2. |
| trappist1e_hab2_sw.nc | as above | — | 3250, 8000, 14000, 30000 | 2 | no | |
| titan_lw.nc | 95% N2 + 5% CH4, incl. N2-N2 & N2-CH4 CIA | 10, 500, 1250, 2500, 3250 | — | 2 | no | CIA dominates rotational window (band 0). |
| titan_sw.nc | as above | — | 3250, 8000, 14000, 30000 | 2 | no | Stellar source from `sun.npz`. |

THAI protocol: Turbet et al. (2022), *ApJ* 924, 10; Sergeev et al. (2022), *ApJ* 924, 11.
Titan CIA: Karman et al. (2019), HITRAN CIA database (CC-BY).
