"""Regenerate every shipped picket-fence correlated-k table from the Chaverot
(Zenodo 16795590, CC-BY 4.0) source HDF5 files.

Run from the project root with an environment that has `exo_k` installed
(climt's runtime environments do **not** bundle exo_k; use the dedicated
`radiation` conda env or install exo_k manually):

    source ~/miniconda3/etc/profile.d/conda.sh && conda activate radiation
    python scripts/regenerate_shipped_picket_fence_tables.py \\
        --chaverot-dir /Users/joymonteiro/github/chaverot/data/correlated-k_tables

Each invocation writes five ``.nc`` files under
``climt/_data/picket_fence/correlated_k/``:

    earth_low_res_sw.nc
    mars_lw.nc   mars_sw.nc
    venus_lw.nc  venus_sw.nc

The Earth LW table (``earth_low_res_lw.nc``) is built from HITRAN 2024 line
data via ``scripts/generate_pf_tables_linepyline.py``; it is NOT regenerated
by this script.

This script is a thin wrapper around ``generate_picket_fence_tables.py`` — it
just fixes the per-planet source file, band edges, g-point count, and mixture
molar mass. The H2O VMR axis is preserved in the output (trilinear runtime
interpolation); we no longer collapse it at build time. See that script's
module docstring for the rationale.

CONTINUUM AND CIA
-----------------
Chaverot builds the LBL spectra with MT_CKD H2O continuum (self + foreign),
CO2 continuum, and CIA pairs (N2-N2, CO2-CO2, cross terms) folded in before
correlated-k compression. So continuum and CIA are already inside k at
Chaverot's background composition:
  * H2O self- and foreign-continuum: resolved along the X_H2O axis (correct).
  * N2-N2 / CO2-CO2 CIA and continuum: frozen at Chaverot's fixed CO2 / N2
    fraction. Gas-forcing experiments (e.g. 2×CO2) therefore get wrong CIA
    and wrong line absorption — these tables are not designed for that.

Future tables generated via the ``linepyline`` workflow MUST include
continuum and CIA terms in the LBL spectrum, or per-pair runtime CIA
handling must be added to the component. Dropping either reduces accuracy
by a large factor in CO2-rich atmospheres (Mars, Venus).

WHY EACH PARAMETER IS SET THE WAY IT IS
---------------------------------------
* ``bands`` — picket-fence resolution bands in wavenumber (cm^-1).
    LW bands cover 10–3250 cm^-1 (~3.08 µm to 1 mm), split at 500/1250/2500
      which roughly mark the H2O rotational / window / CO2 / H2O-vib boundaries.
    SW bands cover 3250–30000 cm^-1 (~0.33–3.08 µm), split at 8000/14000 to
      separate visible / near-IR / UV regions.

* ``ngpt`` — number of correlated-k g-points per band after re-binning. Two
    g-points is the low-resolution target for picket-fence use; it preserves
    the essential in-band line-strength variability while staying cheap.

* ``mixture-molar-mass`` — used to convert Chaverot k from m^2/molecule to
    m^2/kg-of-atmosphere via Avogadro / M_mix. Must match the *bulk* mixture
    composition the Chaverot table was computed for (not the trace gas):
        Earth: dry-air-like 0.028964 kg/mol
        Mars:  ~96% CO2 + 2.7% N2: 0.04355 kg/mol
        Venus: ~96.5% CO2 + 3.5% N2: 0.04344 kg/mol

SOURCE FILE MAPPING
-------------------
The Chaverot dataset ships multiple CO2 concentrations. We pick:
    Earth → Earth_376ppmCO2_R500_10-30k.ktable.SI.h5  (present-day)
    Mars  → CO2-N2_R500_1-30k.ktable.SI.h5            (Martian composition)
    Venus → H2O-CO2_R500_1-30k.ktable.SI.h5           (Venusian composition)
"""
import argparse
import os
import subprocess
import sys

CLIMT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(
    CLIMT_ROOT, "climt", "_data", "picket_fence", "correlated_k"
)
GENERATOR = os.path.join(CLIMT_ROOT, "scripts", "generate_picket_fence_tables.py")

# Per-planet, per-band-kind recipes. Every entry describes one output .nc.
#     source:  basename of the Chaverot HDF5 file (looked up under --chaverot-dir)
#     output:  filename under climt/_data/picket_fence/correlated_k/
#     kind:    "lw" | "sw" — controls which auxiliary arrays are written
#              (planck_fraction for LW, solar_source_per_gpoint + rayleigh for SW)
#     bands:   comma-separated N+1 wavenumber edges (cm^-1)
#     ngpt:    g-points per band after re-quadrature
#     molar_mass: kg/mol for the m^2/molecule → m^2/kg conversion
# The H2O VMR axis from the Chaverot Ktable5d is preserved in the output for
# all three planets (runtime trilinear interpolation).
RECIPES = [
    # --- Earth SW ---
    # Earth LW (earth_low_res_lw.nc) is built from linepyline; not in this script.
    {
        "source": "Earth_376ppmCO2_R500_10-30k.ktable.SI.h5",
        "output": "earth_low_res_sw.nc",
        "kind": "sw",
        "bands": "3250,8000,14000,30000",
        "ngpt": 2,
        "molar_mass": 0.028964,
    },
    # --- Mars LW: option B ---
    {
        "source": "CO2-N2_R500_1-30k.ktable.SI.h5",
        "output": "mars_lw.nc",
        "kind": "lw",
        "bands": "10,500,800,1800,3250",
        "ngpt": 2,
        "molar_mass": 0.04355,
    },
    {
        "source": "CO2-N2_R500_1-30k.ktable.SI.h5",
        "output": "mars_sw.nc",
        "kind": "sw",
        "bands": "3250,8000,14000,30000",
        "ngpt": 2,
        "molar_mass": 0.04355,
    },
    # --- Venus LW: option B ---
    {
        "source": "H2O-CO2_R500_1-30k.ktable.SI.h5",
        "output": "venus_lw.nc",
        "kind": "lw",
        "bands": "10,500,800,1800,3250",
        "ngpt": 2,
        "molar_mass": 0.04344,
    },
    {
        "source": "H2O-CO2_R500_1-30k.ktable.SI.h5",
        "output": "venus_sw.nc",
        "kind": "sw",
        "bands": "3250,8000,14000,30000",
        "ngpt": 2,
        "molar_mass": 0.04344,
    },
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chaverot-dir", required=True,
                    help="Directory containing the Chaverot *.ktable.SI.h5 files")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print the commands that would run, but do not execute")
    args = ap.parse_args()

    if not os.path.isdir(args.chaverot_dir):
        sys.exit(f"chaverot-dir not found: {args.chaverot_dir}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for r in RECIPES:
        source = os.path.join(args.chaverot_dir, r["source"])
        output = os.path.join(OUTPUT_DIR, r["output"])
        cmd = [
            sys.executable, GENERATOR,
            "--input", source,
            "--output", output,
            "--kind", r["kind"],
            "--bands", r["bands"],
            "--ngpt", str(r["ngpt"]),
            "--mixture-molar-mass", repr(r["molar_mass"]),
        ]
        print("$", " ".join(cmd))
        if not args.dry_run:
            subprocess.check_call(cmd)


if __name__ == "__main__":
    main()
