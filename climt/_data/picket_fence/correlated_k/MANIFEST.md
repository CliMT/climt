# Shipped correlated-k tables

All `.nc` tables below are derived from **G. Chaverot's Zenodo dataset 16795590**
(CC-BY 4.0, https://zenodo.org/records/16795590), re-binned to picket-fence
resolution by `scripts/generate_picket_fence_tables.py`.

| File | Source mixture | Bands (cm^-1) | g-points | sha256 |
|---|---|---|---|---|
| earth_low_res_lw.nc | H2O+N2+CO2 (376 ppm) | 10, 500, 1250, 2500, 3250 | 2 | (fill after build) |
| earth_low_res_sw.nc | H2O+N2+CO2 (376 ppm) | 3250, 8000, 14000, 30000 | 2 | (fill after build) |
| mars_lw.nc | CO2+N2 (96% CO2) | 10, 500, 1250, 2500, 3250 | 2 | (fill after build) |
| mars_sw.nc | CO2+N2 (96% CO2) | 3250, 8000, 14000, 30000 | 2 | (fill after build) |
| venus_lw.nc | H2O+CO2 (96% CO2) | 10, 500, 1250, 2500, 3250 | 2 | (fill after build) |
| venus_sw.nc | H2O+CO2 (96% CO2) | 3250, 8000, 14000, 30000 | 2 | (fill after build) |

Tables are pre-mixed at the listed VMRs; gas-forcing experiments
(e.g., CO₂ doubling) are not supported by these tables. For per-gas tables,
use the `linepyline` workflow.

Fill the `sha256` column after generating files with:

```sh
cd climt/_data/picket_fence/correlated_k && shasum -a 256 *.nc
```

The `.npz` files in this directory are pre-existing unit-test fixtures and
are not subject to the Chaverot provenance above.
