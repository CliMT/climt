# scripts/build_land_mask.py
"""Generate the bundled 2-degree Earth land/sea/land_ice mask.

Build-time only: requires `regionmask` and `numpy`. Run with the climt env:
    ~/miniconda3/envs/climt/bin/python scripts/build_land_mask.py
Output: climt/_data/land_mask/earth_landmask_2deg.nc
"""
import os
import numpy as np
import xarray as xr
import regionmask

OUT = os.path.join(
    os.path.dirname(__file__), "..", "climt", "_data", "land_mask",
    "earth_landmask_2deg.nc",
)

def main():
    lat = np.arange(-89.0, 90.0, 2.0)      # ascending, cell centres
    lon = np.arange(1.0, 360.0, 2.0)       # 0..360 cell centres
    land = regionmask.defined_regions.natural_earth_v5_0_0.land_110
    land_mask = land.mask(lon, lat)        # NaN over ocean, 0 over land
    code = np.where(np.isnan(land_mask.values), 0, 1).astype(np.int8)  # sea=0, land=1

    # Tag permanent ice sheets as land_ice (code 2): Antarctica (lat < -63)
    # and Greenland (58<lat<84, 285<lon<350) where currently land.
    LAT, LON = np.meshgrid(lat, lon, indexing="ij")
    antarctica = LAT < -63.0
    greenland = (LAT > 58.0) & (LAT < 84.0) & (LON > 285.0) & (LON < 350.0)
    code[(code == 1) & (antarctica | greenland)] = 2

    da = xr.DataArray(
        code, coords={"lat": lat, "lon": lon}, dims=("lat", "lon"),
        name="area_type_code",
    )
    da.attrs["flag_values"] = [0, 1, 2]
    da.attrs["flag_meanings"] = "sea land land_ice"
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    da.to_dataset().to_netcdf(OUT, engine="scipy")
    print("wrote", OUT, "shape", code.shape)

if __name__ == "__main__":
    main()
