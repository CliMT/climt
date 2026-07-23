import os
import numpy as np
import xarray as xr
from sympl import DiagnosticComponent
# scipy is imported lazily inside the method that needs it so that
# ``import climt`` stays scipy-free for the Pyodide/pure-wheel target
# (see tests/test_no_scipy_import.py).

_DATA = os.path.join(os.path.dirname(__file__), "..", "..", "_data")
_DEFAULT_MASK = os.path.join(_DATA, "land_mask", "earth_landmask_2deg.nc")
_DEFAULT_TOPO = os.path.join(_DATA, "topography", "earth_topography_2deg.nc")
_CODE_TO_NAME = {0: b"sea", 1: b"land", 2: b"land_ice"}


class LandMask(DiagnosticComponent):
    """Set ``area_type`` and (optionally) the topographic boundary forcing.

    Nearest-neighbour maps a bundled (or user-supplied) categorical geography
    onto the model grid, at arbitrary resolution, to set ``area_type``.
    ``sea_ice`` is not set here (it is dynamic and owned by :class:`SeaIce`).

    When ``load_topography=True`` (the default), it *also* bilinearly maps the
    bundled Earth orography and ice-sheet thickness onto the model grid and
    emits ``surface_geopotential`` (read by the GFS dynamical core) and
    ``land_ice_thickness`` (read by :class:`SeaIce`/:class:`LandIce`). Both are
    non-negative and re-zeroed on ``sea`` cells using the ``area_type`` computed
    in the same call, so the categorical geography and the continuous
    topographic fields stay mutually consistent. (Loading a topographic forcing
    here is a deliberate separation-of-concerns compromise: it is a static
    boundary condition tied to the same geography, so one component owns the
    whole land/sea/ice/orography specification.)
    """

    input_properties = {
        "latitude": {"dims": ["*"], "units": "degrees_north"},
        "longitude": {"dims": ["*"], "units": "degrees_east"},
    }
    diagnostic_properties = {
        "area_type": {"dims": ["*"], "units": "dimensionless"},
    }

    def __init__(self, mask_dataset=None, include_land_ice=True,
                 topography_dataset=None, load_topography=True, **kwargs):
        source = mask_dataset if mask_dataset is not None else _DEFAULT_MASK
        ds = source if isinstance(source, xr.Dataset) else xr.open_dataset(
            source, engine="scipy")
        self._src_lat = np.asarray(ds["lat"].values, dtype=float)
        self._src_lon = np.asarray(ds["lon"].values, dtype=float)
        self._src_code = np.asarray(ds["area_type_code"].values)
        self._include_land_ice = include_land_ice
        self._weights = None  # (lat_idx, lon_idx) cached on first call

        self.diagnostic_properties = {
            "area_type": {"dims": ["*"], "units": "dimensionless"},
        }
        self._load_topography = load_topography
        if load_topography:
            tsrc = topography_dataset if topography_dataset is not None \
                else _DEFAULT_TOPO
            tds = tsrc if isinstance(tsrc, xr.Dataset) else xr.open_dataset(
                tsrc, engine="scipy")
            self._topo_lat = np.asarray(tds["lat"].values, dtype=float)
            self._topo_lon = np.asarray(tds["lon"].values, dtype=float)
            self._geopotential = np.asarray(
                tds["surface_geopotential"].values, dtype=float)
            self._ice_thickness = np.asarray(
                tds["land_ice_thickness"].values, dtype=float)
            self.diagnostic_properties.update({
                "surface_geopotential": {"dims": ["*"], "units": "m^2 s^-2"},
                "land_ice_thickness": {"dims": ["*"], "units": "m"},
            })
        super(LandMask, self).__init__(**kwargs)

    def _build_weights(self, lat_flat, lon_flat):
        lon_wrapped = np.mod(lon_flat, 360.0)
        lat_idx = np.abs(self._src_lat[None, :] - lat_flat[:, None]).argmin(axis=1)
        lon_idx = np.abs(self._src_lon[None, :] - lon_wrapped[:, None]).argmin(axis=1)
        self._weights = (lat_idx, lon_idx)

    def _interp_topo(self, field, lat_flat, lon_wrapped):
        from scipy.interpolate import RegularGridInterpolator
        interp = RegularGridInterpolator(
            (self._topo_lat, self._topo_lon), field,
            bounds_error=False, fill_value=None)   # finite edge extrapolation
        return np.clip(interp(np.c_[lat_flat, lon_wrapped]), 0.0, None)

    def array_call(self, state):
        lat = np.asarray(state["latitude"])
        shape = lat.shape
        lat_flat = np.reshape(lat, (-1,)).astype(float)
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

        diagnostics = {"area_type": np.reshape(area, shape)}

        if self._load_topography:
            lon_wrapped = np.mod(lon_flat, 360.0)
            geopotential = self._interp_topo(
                self._geopotential, lat_flat, lon_wrapped)
            ice = self._interp_topo(
                self._ice_thickness, lat_flat, lon_wrapped)
            # tie the continuous fields to the geography from this same call:
            # oceans are flat and ice-free.
            sea = area == b"sea"
            geopotential[sea] = 0.0
            ice[sea] = 0.0
            diagnostics["surface_geopotential"] = np.reshape(geopotential, shape)
            diagnostics["land_ice_thickness"] = np.reshape(ice, shape)

        return diagnostics
