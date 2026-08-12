"""Convert a CORK correlated-k .nc table to .npz (numpy-native, scipy-free reader)."""
import sys
from pathlib import Path
import numpy as np
from climt._components.cork.optics.correlated_k import _load_netcdf_table


def convert(nc_path):
    nc_path = Path(nc_path)
    table = _load_netcdf_table(str(nc_path))  # uses scipy at build time (dev env only)
    out = nc_path.with_suffix(".npz")
    np.savez_compressed(out, **{k: np.asarray(v) for k, v in table.items()})
    print(f"wrote {out} ({out.stat().st_size} bytes)")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        convert(p)
