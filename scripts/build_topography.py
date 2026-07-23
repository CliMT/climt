"""Generate the bundled Earth topography + land-ice-thickness boundary file.

Build-time only. Derives **two** boundary fields from one ETOPO 2022 ice-surface /
bedrock DEM pair and writes them, on the same 2-degree grid as the bundled land
mask, as a netCDF-3 classic file (loaders open boundary files with xarray's
``scipy`` engine, which reads classic netCDF only):

- ``surface_geopotential`` = g * surface elevation, zeroed over sea cells.
  Read by the GFS dynamical core (``set_topography``) to force realistic
  orographically-driven circulation.
- ``land_ice_thickness`` = (ice-surface - bedrock) elevation, clamped >= 0 and
  zeroed over sea cells. Read by ``SeaIce``/``LandIce`` as the ice-sheet state.

The source is downsampled server-side over OPeNDAP (a ~0.25 deg strided read of
the 60 arc-second global tiles -- the full grids are ~900 MB each) and then
area-binned onto the 2 deg mask grid; sea cells are identified from the bundled
land mask so ocean bathymetry never leaks into the orography or ice thickness.

Run with the climt env:
    ~/miniconda3/envs/climt/bin/python scripts/build_topography.py
Output: climt/_data/topography/earth_topography_2deg.nc

Source: NOAA NCEI ETOPO 2022 Global Relief Model (ice-surface & bedrock),
        https://www.ncei.noaa.gov/products/etopo-global-relief-model
"""
import os

import numpy as np
import xarray as xr

G = 9.80665  # standard gravity, matches climt's gravitational_acceleration

_BASE = "https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s"
SURF_URL = _BASE + "/60s_surface_elev_netcdf/ETOPO_2022_v1_60s_N90W180_surface.nc"
BED_URL = _BASE + "/60s_bed_elev_netcdf/ETOPO_2022_v1_60s_N90W180_bed.nc"
STRIDE = 15  # 60 arc-sec * 15 = 0.25 deg source read

CACHE = os.path.join(
    os.path.dirname(__file__), "..", "build_cache", "etopo_quarter_deg.npz")
MASK = os.path.join(
    os.path.dirname(__file__), "..", "climt", "_data", "land_mask",
    "earth_landmask_2deg.nc")
OUT = os.path.join(
    os.path.dirname(__file__), "..", "climt", "_data", "topography",
    "earth_topography_2deg.nc")


def _read_strided():
    """Return (lat, lon[0-360], surface, bed) at ~0.25 deg, cached on disk."""
    if os.path.exists(CACHE):
        d = np.load(CACHE)
        return d["lat"], d["lon"], d["surface"], d["bed"]
    surf = xr.open_dataset(SURF_URL, decode_times=False)
    bed = xr.open_dataset(BED_URL, decode_times=False)
    lat = np.asarray(surf["lat"][::STRIDE].values, float)
    lon = np.asarray(surf["lon"][::STRIDE].values, float)
    print("reading ETOPO surface over OPeNDAP (~0.25 deg) ...")
    s = np.asarray(surf["z"][::STRIDE, ::STRIDE].values, float)
    print("reading ETOPO bedrock over OPeNDAP (~0.25 deg) ...")
    b = np.asarray(bed["z"][::STRIDE, ::STRIDE].values, float)
    lon = np.mod(lon, 360.0)  # -180..180 -> 0..360
    order = np.argsort(lon)
    lon, s, b = lon[order], s[:, order], b[:, order]
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    np.savez(CACHE, lat=lat, lon=lon, surface=s, bed=b)
    return lat, lon, s, b


def _bin_to_2deg(src_lat, src_lon, field, tgt_lat, tgt_lon):
    """Area-mean a fine (lat, lon) field onto the 2 deg target-cell centres."""
    lat_edges = np.concatenate(([tgt_lat[0] - 1.0], tgt_lat + 1.0))   # -90..90
    lon_edges = np.concatenate(([tgt_lon[0] - 1.0], tgt_lon + 1.0))   # 0..360
    ilat = np.clip(np.digitize(src_lat, lat_edges) - 1, 0, tgt_lat.size - 1)
    ilon = np.clip(np.digitize(src_lon, lon_edges) - 1, 0, tgt_lon.size - 1)
    acc = np.zeros((tgt_lat.size, tgt_lon.size))
    cnt = np.zeros((tgt_lat.size, tgt_lon.size))
    IL, JL = np.meshgrid(ilat, ilon, indexing="ij")
    np.add.at(acc, (IL, JL), field)
    np.add.at(cnt, (IL, JL), 1.0)
    return acc / np.maximum(cnt, 1.0)


def main():
    lat, lon, surface, bed = _read_strided()

    md = xr.open_dataset(MASK, engine="scipy")
    tgt_lat = np.asarray(md["lat"].values, float)   # -89..89 step 2
    tgt_lon = np.asarray(md["lon"].values, float)   # 1..359 step 2 (0-360)
    is_sea = np.asarray(md["area_type_code"].values) == 0

    surf2 = _bin_to_2deg(lat, lon, surface, tgt_lat, tgt_lon)
    bed2 = _bin_to_2deg(lat, lon, bed, tgt_lat, tgt_lon)

    # orography: elevation on land, zero over sea, no negative land cells
    surface_height = np.where(is_sea, 0.0, np.clip(surf2, 0.0, None))
    surface_geopotential = (G * surface_height).astype(np.float32)

    # ice thickness: ice-surface minus bedrock, clamped >= 0, zero over sea
    ice = np.clip(surf2 - bed2, 0.0, None)
    land_ice_thickness = np.where(is_sea, 0.0, ice).astype(np.float32)

    ds = xr.Dataset(
        {
            "surface_geopotential": (("lat", "lon"), surface_geopotential),
            "land_ice_thickness": (("lat", "lon"), land_ice_thickness),
        },
        coords={"lat": tgt_lat, "lon": tgt_lon},
    )
    ds["surface_geopotential"].attrs.update({
        "units": "m^2 s^-2", "long_name": "surface geopotential (g * orography)",
        "standard_name": "surface_geopotential"})
    ds["land_ice_thickness"].attrs.update({
        "units": "m", "long_name": "land ice-sheet thickness",
        "standard_name": "land_ice_thickness"})
    ds["lat"].attrs.update({"units": "degrees_north", "axis": "Y"})
    ds["lon"].attrs.update({"units": "degrees_east", "axis": "X"})
    ds.attrs.update({
        "title": "Earth orography + land-ice thickness for climt (2 deg)",
        "source": "NOAA NCEI ETOPO 2022 (ice-surface and bedrock), 60 arc-sec, "
                  "area-averaged to 2 deg",
        "source_url": "https://www.ncei.noaa.gov/products/etopo-global-relief-model",
        "gravity": G,
        "note": "sea cells (from earth_landmask_2deg.nc) zeroed in both fields; "
                "see scripts/build_topography.py",
    })

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    ds.to_netcdf(OUT, engine="scipy", format="NETCDF3_CLASSIC")
    print("wrote", OUT,
          "| max orog %.0f m" % surface_height.max(),
          "| max ice %.0f m" % land_ice_thickness.max())


if __name__ == "__main__":
    main()
