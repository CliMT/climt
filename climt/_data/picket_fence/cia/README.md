# HITRAN CIA flat files (offline build only)

These ASCII files are consumed by `scripts/pf_table_builder/cia.py` when
building the Titan opacity table. They are not loaded at climt runtime —
the CIA contribution is baked into the resulting `.nc` k-table.

Source: <https://hitran.org/cia/> (Karman et al., 2019; CC-BY).

Files expected (download once):
- `N2-N2_2018.cia`
- `N2-CH4_2018.cia`

If your HITRAN download supplies only `CH4-N2_…`, that is the symmetric
pair — the loader accepts either ordering.
