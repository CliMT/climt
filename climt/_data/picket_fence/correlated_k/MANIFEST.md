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
| earth_low_res_lw.nc | Earth (376 ppm CO₂) | LW option B | 10, 500, 800, 1800, 3250 | 2 | `94570a87d18c0c1f022d633433a663a0940816dc34c25693eb9ffd4aed1bfb82` |
| earth_low_res_lw_4gpt.nc | Earth (376 ppm CO₂) | LW option B | 10, 500, 800, 1800, 3250 | 4 | `8c9f83e2b38cfcd577cbc67f02d750afca80d0e2c2b70443d20825b4ccd04fe7` |
| earth_low_res_lw_8gpt.nc | Earth (376 ppm CO₂) | LW option B | 10, 500, 800, 1800, 3250 | 8 | `15b06a4311eae5bc0fecd255c54e7a810c80d90634edd2a5d27ebe7c2ef49136` |
| earth_low_res_lw_hotmoist.nc | Earth (376 ppm CO₂) | LW option D | 10, 500, 1200, 2000, 3250 | 2 | `9ae21cd2edf1c54bf64c7896673f1e5f93838fb3aa9a870c3f29eaa6b54aaccd` |
| earth_low_res_sw.nc | Earth (376 ppm CO₂) | SW | 3250, 8000, 14000, 30000 | 2 | `7e3678235f91cb3296d1b6e03e0873b31b59cf0947f8539d0a45e6811fd9a3e2` |
| mars_lw.nc | Mars (CO₂+N₂) | LW option B | 10, 500, 800, 1800, 3250 | 2 | `94e61290e54be8f3552b7762014bc09443bd5e7f0d94d39932c726c4198d78c7` |
| mars_sw.nc | Mars (CO₂+N₂) | SW | 3250, 8000, 14000, 30000 | 2 | `983ce8dc67460a54466f7a5925c2dcfde5ffd4d880ec614a50f107d0ed834064` |
| venus_lw.nc | Venus (H₂O+CO₂) | LW option B | 10, 500, 800, 1800, 3250 | 2 | `f308088c1cc2d6f028e44b1be3a004d4c01000d6347a4bc5bfff21dd285fe6b8` |
| venus_sw.nc | Venus (H₂O+CO₂) | SW | 3250, 8000, 14000, 30000 | 2 | `e2a040cc69d53058cc4300e07fdf9c493286a8e9802cfeec7c413dcf2cc3c2b7` |

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
