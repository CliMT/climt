import os
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
band_flux_dir = os.path.join(script_dir, "plots", "band_fluxes")
out_dir = os.path.join(script_dir, "plots")

COMPOSITION_CASES = ['Ben1_N2_CO2', 'Ben2_CO2_only']


def parse_last_block(path):
    """Read a band_fluxes_*.txt file and return the rows for the last
    (most recent / most converged) call_index, sorted by wl_short."""
    rows = []
    last_idx = None
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            idx = int(parts[0])
            if last_idx is None or idx > last_idx:
                last_idx = idx
                rows = []
            if idx == last_idx:
                rows.append({
                    'band': int(parts[1]),
                    'wl_short': float(parts[2]),
                    'wl_long': float(parts[3]),
                    'up_toa': float(parts[5]),
                    'col_abs': float(parts[10]),
                })
    rows.sort(key=lambda r: r['wl_short'])
    return rows


def band_edges_and_values(rows, metric):
    edges = [rows[0]['wl_short']] + [r['wl_long'] for r in rows]
    values = [r[metric] for r in rows]
    return edges, values


def plot_panel(ax, mars_rows, trap_rows, metric, xlim, title):
    edges_m, vals_m = band_edges_and_values(mars_rows, metric)
    edges_t, vals_t = band_edges_and_values(trap_rows, metric)

    ax.stairs(vals_t, edges_t, color='#0066cc', linewidth=2.0, label='TRAPPIST-1e', baseline=None)
    ax.stairs(vals_m, edges_m, color='#cc0066', linewidth=2.0, linestyle='--', label='DSA Mars', baseline=None)

    ax.axhline(0, color='0.6', linewidth=0.8, zorder=0)
    ax.set_xscale('log')
    ax.set_xlim(*xlim)
    ax.set_xlabel('Wavelength (μm)')
    ylabel = 'Outgoing flux at TOA (W/m²)' if metric == 'up_toa' else 'Net column absorption (W/m²)'
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11)
    ax.grid(True, which='both', ls=':', alpha=0.5)
    ax.legend(fontsize=9, loc='best')


def make_figure(case_name, case_title, metric):
    mars_lw = parse_last_block(os.path.join(band_flux_dir, f"band_fluxes_dsa_mars_{case_name}_lw.txt"))
    trap_lw = parse_last_block(os.path.join(band_flux_dir, f"band_fluxes_trappist1e_{case_name}_lw.txt"))
    mars_sw = parse_last_block(os.path.join(band_flux_dir, f"band_fluxes_dsa_mars_{case_name}_sw.txt"))
    trap_sw = parse_last_block(os.path.join(band_flux_dir, f"band_fluxes_trappist1e_{case_name}_sw.txt"))

    fig, (ax_lw, ax_sw) = plt.subplots(2, 1, figsize=(9, 9))

    plot_panel(ax_lw, mars_lw, trap_lw, metric, xlim=(3, 300),
               title='Longwave (band 1 extends to 10000 μm, clipped here)')
    plot_panel(ax_sw, mars_sw, trap_sw, metric, xlim=(0.2, 20),
               title='Shortwave')

    metric_label = 'Outgoing/reflected flux at TOA' if metric == 'up_toa' else 'Net column absorption'
    fig.suptitle(f'{metric_label} vs wavelength\n{case_title}', fontsize=13, fontweight='bold')
    fig.tight_layout()

    out_path = os.path.join(out_dir, f"flux_spectrum_{case_name}_{metric}.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


CASE_TITLES = {
    'Ben1_N2_CO2': 'Ben1/Hab1 (1 bar N2 + 400 ppm CO2)',
    'Ben2_CO2_only': 'Ben2/Hab2 (1 bar pure CO2)',
}

if __name__ == '__main__':
    for case_name in COMPOSITION_CASES:
        for metric in ('up_toa', 'col_abs'):
            make_figure(case_name, CASE_TITLES[case_name], metric)
