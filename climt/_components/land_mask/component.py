import os
import numpy as np
import xarray as xr
from sympl import DiagnosticComponent

_DEFAULT_MASK = os.path.join(
    os.path.dirname(__file__), "..", "..", "_data", "land_mask",
    "earth_landmask_2deg.nc",
)
_CODE_TO_NAME = {0: b"sea", 1: b"land", 2: b"land_ice"}


class LandMask(DiagnosticComponent):
    """Set ``area_type`` to a present-day Earth land/sea/land_ice configuration.

    Nearest-neighbour maps a bundled (or user-supplied) categorical geography
    onto the model grid, at arbitrary resolution. ``sea_ice`` is not set here
    (it is dynamic and owned by :class:`SeaIce`).
    """

    input_properties = {
        "latitude": {"dims": ["*"], "units": "degrees_north"},
        "longitude": {"dims": ["*"], "units": "degrees_east"},
    }
    diagnostic_properties = {
        "area_type": {"dims": ["*"], "units": "dimensionless"},
    }

    def __init__(self, mask_dataset=None, include_land_ice=True, **kwargs):
        source = mask_dataset if mask_dataset is not None else _DEFAULT_MASK
        ds = source if isinstance(source, xr.Dataset) else xr.open_dataset(
            source, engine="scipy")
        self._src_lat = np.asarray(ds["lat"].values, dtype=float)
        self._src_lon = np.asarray(ds["lon"].values, dtype=float)
        self._src_code = np.asarray(ds["area_type_code"].values)
        self._include_land_ice = include_land_ice
        self._weights = None  # (lat_idx, lon_idx) cached on first call
        super(LandMask, self).__init__(**kwargs)

    def _build_weights(self, lat_flat, lon_flat):
        lon_wrapped = np.mod(lon_flat, 360.0)
        lat_idx = np.abs(self._src_lat[None, :] - lat_flat[:, None]).argmin(axis=1)
        lon_idx = np.abs(self._src_lon[None, :] - lon_wrapped[:, None]).argmin(axis=1)
        self._weights = (lat_idx, lon_idx)

    def array_call(self, state):
        lat = state["latitude"]
        lat_flat = np.reshape(np.asarray(lat), (-1,)).astype(float)
        lon_flat = np.reshape(np.asarray(state["longitude"]), (-1,)).astype(float)
        if self._weights is None:
            self._build_weights(lat_flat, lon_flat)
        lat_idx, lon_idx = self._weights
        codes = self._src_code[lat_idx, lon_idx]
        area = np.empty(codes.shape, dtype="S100")
        for c, name in _CODE_TO_NAME.items():
            sel = codes == c
            if c == 2 and not self._include_land_ice:
                name = b"land"
            area[sel] = name
        return {"area_type": np.reshape(area, np.asarray(lat).shape)}
