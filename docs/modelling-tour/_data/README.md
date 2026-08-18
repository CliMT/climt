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
