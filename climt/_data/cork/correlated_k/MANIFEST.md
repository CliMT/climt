# Shipped correlated-k tables

The Earth LW table (`earth_low_res_lw`) is built from HITRAN 2024 line data via
the [linepyline](https://github.com/exoclime/linepyline) pipeline (MT_CKD 4.3
H₂O continuum, decoupled log-X interpolation, runtime CO₂ axis). All other
planet tables are derived from **G. Chaverot's Zenodo dataset 16795590**
(CC-BY 4.0, https://zenodo.org/records/16795590), re-binned to cork
resolution by `scripts/generate_cork_tables.py`.

SW bands are unchanged across all planets: `3250, 8000, 14000, 30000` cm⁻¹.

## File listing

| File | Planet | Kind | Bands (cm⁻¹) | g-pts | sha256 |
|---|---|---|---|---|---|
| earth_low_res_lw.nc | Earth (CO₂ 10–10000 ppm runtime axis) | **DEFAULT** LW — 14-band candD-class, decoupled+log-X H₂O continuum, runtime CO₂ | 10…3250 (14-band) | 8 | `3192c45f2c4bad5243023177841089568270f246280bc5bbbe0dc74bac74b760` |
| earth_low_res_sw.nc | Earth (376 ppm CO₂) | SW | 3250, 8000, 14000, 30000 | 2 | `8709799c08dfaaf015a81a931c1eec699649559803b8fd2831eeed68d6360bda` |
| mars_lw.nc | Mars (CO₂+N₂) | LW option B | 10, 500, 800, 1800, 3250 | 2 | `94e61290e54be8f3552b7762014bc09443bd5e7f0d94d39932c726c4198d78c7` |
| mars_sw.nc | Mars (CO₂+N₂) | SW | 3250, 8000, 14000, 30000 | 2 | `8e186ac238e1bf7f562348f3fc70bbc006358efeb4f0d26c4aab7756a73157bb` |
| venus_lw.nc | Venus (H₂O+CO₂) | LW option B | 10, 500, 800, 1800, 3250 | 2 | `f308088c1cc2d6f028e44b1be3a004d4c01000d6347a4bc5bfff21dd285fe6b8` |
| venus_sw.nc | Venus (H₂O+CO₂) | SW | 3250, 8000, 14000, 30000 | 2 | `606960889f56af3b7c49e93c728be04c1a80b5555df3643cff3af94c2adbf83f` |

## Format

Tables with an `h2o_vmr_grid` axis treat H₂O as a **runtime variable**
interpolated from `specific_humidity`. A netCDF attribute
`background_is_premixed = "true"` means all other gases are baked into k and
the component scales by column mass of AIR. The Earth LW table additionally
has a `co2_vmr_grid` axis for runtime CO₂ (10–10000 ppm, 10 log nodes).

## Regenerating

```sh
# Earth LW: requires linepyline env.
python scripts/generate_cork_tables_linepyline.py --scenario earth --kind lw

# Planet tables (Mars, Venus, SW): requires exo_k; use the `radiation` env.
python scripts/regenerate_shipped_cork_tables.py \
    --chaverot-dir /path/to/correlated-k_tables
shasum -a 256 climt/_data/cork/correlated_k/*.nc
```

The `.npz` files in this directory are pre-existing unit-test fixtures and
are not subject to the Chaverot provenance above. `single_band_unit_lw.nc`
is a synthetic constant-k fixture for the grey-limit unit test.

`single_band_gray_lw.nc` is a synthetic single-band, constant-k fixture whose
diffusivity-scaled column optical depth is calibrated to climt's **default**
gray longwave optical depth (`tau0 * (1 - p/ps)`, `tau0 = 1`; see the
`longwave_optical_depth` entry in `climt/_core/initialization.py`). It lets the
CORK component stand in for `GrayLongwaveRadiation` as the "gray" baseline of
the in-browser radiative-equilibrium demo. Regenerate (scipy-free) with

```sh
python scripts/generate_gray_default_table.py \
    --output climt/_data/cork/correlated_k/single_band_gray_lw.nc
python scripts/convert_ck_table_to_npz.py \
    climt/_data/cork/correlated_k/single_band_gray_lw.nc
```

and verify with `tests/test_gray_default_table.py`.

## linepyline-derived tables

Built offline by `scripts/generate_cork_tables_linepyline.py` from
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
