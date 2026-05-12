#!/usr/bin/env python
"""Generate a picket-fence correlated-k table for a named planetary scenario.

Usage:
    conda run -n linepyline python scripts/generate_pf_tables_linepyline.py \\
        --scenario trappist1e_hab2 --kind lw --output climt/_data/picket_fence/correlated_k/trappist1e_hab2_lw.nc

Scenarios:
    trappist1e_hab1  THAI Hab1: N2 bath + 400 ppm CO2 + variable H2O (axis)
    trappist1e_hab2  THAI Hab2: pure 1-bar CO2
    titan            Titan: 95% N2 + 5% CH4, N2-N2 and N2-CH4 CIA included

`--kind lw|sw` selects band edges and whether Planck or solar source is written.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_SCRIPT_DIR)
sys.path.insert(0, _SCRIPT_DIR)
from pf_table_builder.kappa_sampling import sample_kappa_grid
from pf_table_builder.k_distribution import kappa_to_k_coeffs
from pf_table_builder.planck_fraction import build_uniform_planck_fraction
from pf_table_builder.solar_source import (
    build_solar_source_per_gpoint, build_rayleigh_per_band,
)
from pf_table_builder.netcdf_writer import write_lw_table, write_sw_table
from pf_table_builder.cia import cia_kappa_on_grid


# Scenario configuration registry
SCENARIOS = {
    "trappist1e_hab1": dict(
        background_gas="N2",
        absorbers={"CO2": 4e-4},  # H2O handled via h2o_vmr_grid
        h2o_vmr_grid=np.array([1e-7, 1e-5, 1e-3, 1e-2, 1e-1]),
        T_grid=np.linspace(150.0, 350.0, 10),
        p_grid=np.logspace(0.0, 5.05, 15),     # 1 Pa to ~1.1e5 Pa
        stellar_spectrum="trappist1",
        mean_molar_mass_g=28.05,
        rayleigh_refractivity=2.97e-4,         # N2-dominated
        cia_files=[],
    ),
    "trappist1e_hab2": dict(
        background_gas=None,
        absorbers={"CO2": 1.0},
        h2o_vmr_grid=None,
        T_grid=np.linspace(150.0, 350.0, 10),
        p_grid=np.logspace(0.0, 5.05, 15),
        stellar_spectrum="trappist1",
        mean_molar_mass_g=44.01,
        rayleigh_refractivity=4.50e-4,         # CO2
        cia_files=[],
    ),
    "titan": dict(
        background_gas="N2",
        absorbers={"CH4": 0.05},
        h2o_vmr_grid=None,
        T_grid=np.linspace(70.0, 200.0, 10),
        p_grid=np.logspace(-1.0, 5.2, 15),
        stellar_spectrum="sun",
        mean_molar_mass_g=28.6,
        rayleigh_refractivity=2.97e-4,
        cia_files=[
            dict(path="climt/_data/picket_fence/cia/N2-N2_2018.cia",
                 pair=("N2", "N2"), vmr_a=0.95, vmr_b=0.95),
            dict(path="climt/_data/picket_fence/cia/N2-CH4_2018.cia",
                 pair=("N2", "CH4"), vmr_a=0.95, vmr_b=0.05),
        ],
    ),
}

_STELLAR_SPECTRA_DIR = os.path.join(
    _REPO_ROOT, "climt", "_data", "picket_fence", "stellar_spectra"
)


def _load_stellar_spectrum(name_or_path):
    """Load stellar spectrum from a built-in name or .npz path."""
    if os.path.isfile(name_or_path):
        path = name_or_path
    else:
        path = os.path.join(_STELLAR_SPECTRA_DIR, f"{name_or_path}.npz")
    data = np.load(path)
    return {"wavenumber": np.array(data["wavenumber"]),
            "irradiance": np.array(data["irradiance"])}


LW_BAND_EDGES = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
SW_BAND_EDGES = np.array([3250.0, 8000.0, 14000.0, 30000.0])


def build_table(scenario_name, kind, output_path, ngpt=2, dnu=0.5,
                line_shape="lorentz", binning=True):
    cfg = SCENARIOS[scenario_name]
    band_edges = LW_BAND_EDGES if kind == "lw" else SW_BAND_EDGES
    nu_grid = np.arange(band_edges[0], band_edges[-1] + dnu / 2, dnu)

    print(f"[{scenario_name}/{kind}] sampling kappa on "
          f"(nT={len(cfg['T_grid'])}, nP={len(cfg['p_grid'])}, nNu={len(nu_grid)})")
    kappa = sample_kappa_grid(
        background_gas=cfg["background_gas"],
        absorbers=cfg["absorbers"],
        h2o_vmr_grid=cfg["h2o_vmr_grid"],
        T_grid=cfg["T_grid"],
        p_grid=cfg["p_grid"],
        nu_grid=nu_grid,
        line_shape=line_shape,
        binning=binning,
    )

    # Add CIA contribution if any (CIA is independent of X_H2O)
    if cfg["cia_files"]:
        print(f"[{scenario_name}/{kind}] adding {len(cfg['cia_files'])} CIA pair(s)")
        for entry in cfg["cia_files"]:
            kappa_cia = cia_kappa_on_grid(
                entry["path"], pair=entry["pair"],
                vmr_a=entry["vmr_a"], vmr_b=entry["vmr_b"],
                background_gas=cfg["background_gas"],
                absorbers={**cfg["absorbers"],
                           cfg["background_gas"]: 1.0 - sum(cfg["absorbers"].values())},
                T_grid=cfg["T_grid"],
                p_grid=cfg["p_grid"],
                nu_grid=nu_grid,
            )
            # kappa_cia is (nT, nP, nNu); broadcast over X axis if present
            if kappa.ndim == 4:
                kappa = kappa + kappa_cia[:, :, None, :]
            else:
                kappa = kappa + kappa_cia

    print(f"[{scenario_name}/{kind}] building k-distribution (ngpt={ngpt})")
    k_coeffs, gpt_weights = kappa_to_k_coeffs(kappa, nu_grid, band_edges, ngpt)
    # Add the leading "ngas=1" axis expected by the climt loader.
    k_coeffs = k_coeffs[np.newaxis]

    # Rearrange axes to (ngas, nband, ngpt, nT, nP[, nX])
    if kappa.ndim == 4:
        # current k_coeffs shape: (1, nT, nP, nX, nband, ngpt)
        k_coeffs = np.moveaxis(k_coeffs, (1, 2, 3, 4, 5), (3, 4, 5, 1, 2))
    else:
        # current k_coeffs shape: (1, nT, nP, nband, ngpt)
        k_coeffs = np.moveaxis(k_coeffs, (1, 2, 3, 4), (3, 4, 1, 2))

    log_p_grid = np.log(cfg["p_grid"])

    if kind == "lw":
        planck_fraction = build_uniform_planck_fraction(
            len(band_edges) - 1, ngpt, len(cfg["T_grid"]),
        )
        write_lw_table(
            output_path,
            k_coefficients=k_coeffs,
            gpoint_weights=gpt_weights,
            T_grid=cfg["T_grid"],
            log_p_grid=log_p_grid,
            band_edges=band_edges,
            planck_fraction=planck_fraction,
            h2o_vmr_grid=cfg["h2o_vmr_grid"],
            source=f"linepyline:{scenario_name}",
        )
    else:
        # Load stellar spectrum and build solar source / Rayleigh
        spectrum = _load_stellar_spectrum(cfg["stellar_spectrum"])
        solar = build_solar_source_per_gpoint(spectrum, band_edges, ngpt)
        rayleigh = build_rayleigh_per_band(
            band_edges,
            mean_molar_mass_g=cfg["mean_molar_mass_g"],
            refractivity_298K=cfg["rayleigh_refractivity"],
        )
        write_sw_table(
            output_path,
            k_coefficients=k_coeffs,
            gpoint_weights=gpt_weights,
            T_grid=cfg["T_grid"],
            log_p_grid=log_p_grid,
            band_edges=band_edges,
            solar_source_per_gpoint=solar,
            rayleigh_coefficient=rayleigh,
            h2o_vmr_grid=cfg["h2o_vmr_grid"],
            source=f"linepyline:{scenario_name}",
        )
    print(f"[{scenario_name}/{kind}] wrote {output_path}  "
          f"k_coefficients.shape={k_coeffs.shape}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", choices=sorted(SCENARIOS), required=True)
    ap.add_argument("--kind", choices=("lw", "sw"), required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--ngpt", type=int, default=2)
    ap.add_argument("--dnu", type=float, default=0.5,
                    help="wavenumber resolution in cm^-1 (default 0.5)")
    ap.add_argument("--line-shape", default="lorentz",
                    choices=("lorentz", "voigt", "pseudovoigt"))
    args = ap.parse_args()
    build_table(args.scenario, args.kind, args.output,
                ngpt=args.ngpt, dnu=args.dnu, line_shape=args.line_shape)


if __name__ == "__main__":
    main()
