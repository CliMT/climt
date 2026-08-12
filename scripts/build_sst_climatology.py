"""Generate the bundled 12-month climatological SST for ``DataOcean``.

Build-time only: downloads the NOAA OI SST V2 1-degree monthly long-term mean
(1991-2020 base period) from NOAA PSL, reshapes it into the layout
:class:`climt.DataOcean` expects, and writes it as a netCDF-3 classic file
(``DataOcean`` opens boundary files with xarray's ``scipy`` engine, which reads
classic netCDF only -- the PSL source is netCDF-4/HDF5 and would not load).

The result is a single Jan..Dec climatology; ``DataOcean`` interpolates it by
day-of-year with a cyclic Dec->Jan wrap, so the model loops around it
indefinitely to build a long-term climate.

Run with the climt env:
    ~/miniconda3/envs/climt/bin/python scripts/build_sst_climatology.py
Output: climt/_data/data_ocean/earth_sst_climatology_1deg.nc

Source: NOAA OI SST V2, NOAA PSL, Boulder, Colorado, USA
        https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html
"""
import os
import urllib.request

import numpy as np
import xarray as xr

SRC_URL = (
    "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2/"
    "sst.ltm.1991-2020.nc"
)
CACHE = os.path.join(
    os.path.dirname(__file__), "..", "build_cache", "sst.ltm.1991-2020.nc",
)
OUT = os.path.join(
    os.path.dirname(__file__), "..", "climt", "_data", "data_ocean",
    "earth_sst_climatology_1deg.nc",
)


def _download(url, path):
    if os.path.exists(path):
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print("downloading", url)
    urllib.request.urlretrieve(url, path)


def main():
    _download(SRC_URL, CACHE)
    src = xr.open_dataset(CACHE, decode_times=False)

    sst = np.asarray(src["sst"].values, dtype=np.float32)   # (12, 180, 360)
    assert sst.shape[0] == 12, sst.shape
    lat = np.asarray(src["lat"].values, dtype=np.float64)   # ascending -89.5..89.5
    lon = np.asarray(src["lon"].values, dtype=np.float64)   # 0.5..359.5, 0-360

    # DataOcean reads ds["tos"], ds["lat"], ds["lon"]; it does not read the time
    # coordinate (it indexes months by array order and interpolates on the model
    # calendar). Store time as a plain 1..12 month index to avoid any calendar /
    # netCDF-3 time-encoding fragility.
    ds = xr.Dataset(
        {"tos": (("time", "lat", "lon"), sst)},
        coords={
            "time": np.arange(1, 13, dtype=np.int32),
            "lat": lat,
            "lon": lon,
        },
    )
    ds["tos"].attrs.update({
        "units": "degC",
        "long_name": "climatological sea surface temperature",
        "standard_name": "sea_surface_temperature",
    })
    ds["lat"].attrs.update({"units": "degrees_north", "axis": "Y"})
    ds["lon"].attrs.update({"units": "degrees_east", "axis": "X"})
    ds["time"].attrs.update({"units": "month", "long_name": "month of year"})
    ds.attrs.update({
        "title": "Climatological SST for climt.DataOcean (Jan-Dec, loops)",
        "source": "NOAA OI SST V2 monthly long-term mean, 1991-2020 base period",
        "source_url": SRC_URL,
        "acknowledgement": (
            "Data provided by the NOAA PSL, Boulder, Colorado, USA, "
            "from https://psl.noaa.gov"
        ),
        "note": (
            "Reshaped for climt: variable renamed sst->tos, written as "
            "netCDF-3 classic. See scripts/build_sst_climatology.py."
        ),
    })

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    ds.to_netcdf(OUT, engine="scipy", format="NETCDF3_CLASSIC")
    print("wrote", OUT, "shape", sst.shape)


if __name__ == "__main__":
    main()
