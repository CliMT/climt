# Shipped correlated-k tables

All `.nc` tables below are derived from **G. Chaverot's Zenodo dataset 16795590**
(CC-BY 4.0, https://zenodo.org/records/16795590), re-binned to picket-fence
resolution by `scripts/generate_picket_fence_tables.py`.

| File | Source mixture | Bands (cm^-1) | g-points | sha256 |
|---|---|---|---|---|
| earth_low_res_lw.nc | H2O+N2+CO2 (376 ppm, H2O VMR=0.01) | 10, 500, 1250, 2500, 3250 | 2 | `6ea1378af06b54bd84c40c816dd74c0e7eeb7b56686608de81c924d226ad822c` |
| earth_low_res_sw.nc | H2O+N2+CO2 (376 ppm, H2O VMR=0.01) | 3250, 8000, 14000, 30000 | 2 | `6c27357fdc49c810e2dab55e5078abb1d7003e5d0a79290831427ad860e11eb1` |
| mars_lw.nc | CO2+N2 (H2O VMR=1e-5) | 10, 500, 1250, 2500, 3250 | 2 | `70399a0d68bfcaa2c63dc56e2bdd2a389b3052c5058e7faa66ce42c26abd6d4f` |
| mars_sw.nc | CO2+N2 (H2O VMR=1e-5) | 3250, 8000, 14000, 30000 | 2 | `6e53be02b9d1af2764a97053154389b2377d9e12dc1abe57f24dbb90c020aab9` |
| venus_lw.nc | H2O+CO2 (H2O VMR=3e-5) | 10, 500, 1250, 2500, 3250 | 2 | `9f277b6e22626ef3f94b276156e3d41040b42af8564f4a5254b5f624f0b620a4` |
| venus_sw.nc | H2O+CO2 (H2O VMR=3e-5) | 3250, 8000, 14000, 30000 | 2 | `2cd45198a9bec8bc99ca3fe597e36730291a23987e93f5027cdfd918e396bb31` |

Tables are pre-mixed at the listed VMRs; gas-forcing experiments
(e.g., CO₂ doubling) are not supported by these tables. For per-gas tables,
use the `linepyline` workflow.

Fill the `sha256` column after generating files with:

```sh
cd climt/_data/picket_fence/correlated_k && shasum -a 256 *.nc
```

The `.npz` files in this directory are pre-existing unit-test fixtures and
are not subject to the Chaverot provenance above.
