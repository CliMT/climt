import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree
from sympl import DiagnosticComponent
from ._sst_interpolation import mid_month_values, interp_time
from ..._core.surface_fluxes import bulk_fluxes


class DataOcean(DiagnosticComponent):
    """Prescribe observed SSTs on ``area_type == 'sea'`` cells (AMIP-style).

    Loads a (time, lat, lon) sea-surface-temperature dataset (e.g. an AMIP
    boundary-condition file), fills any masked/NaN source points via
    nearest-valid lookup, optionally converts monthly means to mid-month
    values (Taylor/Williamson/Zwiers-style piecewise-linear reconstruction,
    see :mod:`._sst_interpolation`), and on every call time-interpolates and
    spatially interpolates the field onto the model grid. ``sea_surface_
    temperature`` is set everywhere (diagnostic, for coupling to sea-ice
    components); ``surface_temperature`` is overwritten only on
    ``area_type == "sea"`` cells, leaving land/ice cells untouched.
    """

    input_properties = {
        "latitude": {"dims": ["*"], "units": "degrees_north"},
        "longitude": {"dims": ["*"], "units": "degrees_east"},
        "area_type": {"dims": ["*"], "units": "dimensionless"},
        "surface_temperature": {"dims": ["*"], "units": "degK"},
    }
    diagnostic_properties = {
        "sea_surface_temperature": {"dims": ["*"], "units": "degK"},
        "surface_temperature": {"dims": ["*"], "units": "degK"},
    }

    def __init__(self, sst_dataset, sst_variable="tos",
                 time_interpolation="mid_month", relaxation_timescale=None,
                 compute_fluxes=False, **kwargs):
        self._compute_fluxes = compute_fluxes
        if compute_fluxes:
            self.input_properties = dict(self.input_properties)
            self.input_properties.update({
                "eastward_wind": {"dims": ["*"], "units": "m/s"},
                "northward_wind": {"dims": ["*"], "units": "m/s"},
                "air_temperature": {"dims": ["*"], "units": "degK"},
                "specific_humidity": {"dims": ["*"], "units": "kg/kg"},
                "surface_specific_humidity": {"dims": ["*"], "units": "kg/kg"},
                "air_density": {"dims": ["*"], "units": "kg/m^3"},
            })
            self.diagnostic_properties = dict(self.diagnostic_properties)
            self.diagnostic_properties.update({
                "surface_upward_sensible_heat_flux": {
                    "dims": ["*"], "units": "W m^-2"},
                "surface_upward_latent_heat_flux": {
                    "dims": ["*"], "units": "W m^-2"},
            })

        ds = sst_dataset if isinstance(sst_dataset, xr.Dataset) else \
            xr.open_dataset(sst_dataset, engine="scipy")
        sst = ds[sst_variable]
        units = sst.attrs.get("units", "K").strip().lower()
        if units in ("degc", "c", "celsius", "degrees_c", "degrees c",
                     "degrees celsius", "deg_c"):
            sst = sst + 273.15
        self._src_lat = np.asarray(ds["lat"].values, float)
        self._src_lon = np.asarray(ds["lon"].values, float)
        raw = np.asarray(sst.values, float)                 # (12, nlat, nlon)
        self._filled = self._fill_source(raw)               # no NaNs in source
        self._time_interp = time_interpolation
        self._mid = mid_month_values(self._filled) if time_interpolation == \
            "mid_month" else self._filled
        self._relaxation_timescale = relaxation_timescale
        self._weights = None
        super(DataOcean, self).__init__(**kwargs)

    def _fill_source(self, raw):
        # nearest-valid fill of masked/NaN source points, per month
        out = raw.copy()
        for t in range(out.shape[0]):
            layer = out[t]
            mask = ~np.isfinite(layer)
            if not mask.any():
                continue
            LON, LAT = np.meshgrid(self._src_lon, self._src_lat)
            good = ~mask
            tree = cKDTree(np.c_[LAT[good], LON[good]])
            _, idx = tree.query(np.c_[LAT[mask], LON[mask]])
            layer[mask] = layer[good][idx]
        return out

    def _build_weights(self, lat_flat, lon_flat):
        self._lat_flat = lat_flat
        self._lon_flat = np.mod(lon_flat, 360.0)

    def array_call(self, state):
        lat = np.asarray(state["latitude"])
        shape = lat.shape
        lat_flat = np.reshape(lat, (-1,)).astype(float)
        lon_flat = np.reshape(np.asarray(state["longitude"]), (-1,)).astype(float)
        if self._weights is None:
            self._build_weights(lat_flat, lon_flat)
        # time-interpolated source field (nlat, nlon)
        if self._time_interp == "mid_month":
            field = interp_time(self._mid, state["time"])
        else:
            field = self._mid[state["time"].month - 1]
        interp = RegularGridInterpolator(
            (self._src_lat, self._src_lon), field,
            bounds_error=False, fill_value=None)   # None -> linear extrapolation beyond
        # the source grid edges (finite, still bounded in practice since climt's
        # Gaussian model grids stay within ~87 deg of the poles, inside the usual
        # +/-88 deg extent of AMIP-style SST source grids).
        sst_flat = interp(np.c_[self._lat_flat, self._lon_flat])
        area = np.reshape(np.asarray(state["area_type"]).astype(str), (-1,))
        sea = area == "sea"
        surf = np.reshape(np.asarray(state["surface_temperature"]), (-1,)).astype(float).copy()
        surf[sea] = sst_flat[sea]
        assert np.all(np.isfinite(sst_flat[sea]))   # full-sea-coverage guarantee
        diagnostics = {
            "sea_surface_temperature": np.reshape(sst_flat, shape),
            "surface_temperature": np.reshape(surf, shape),
        }
        if self._compute_fluxes:
            fluxes = bulk_fluxes(
                wind_speed=np.sqrt(
                    np.asarray(state["eastward_wind"]) ** 2
                    + np.asarray(state["northward_wind"]) ** 2
                ),
                surface_temperature=np.reshape(surf, shape),
                air_temperature=np.asarray(state["air_temperature"]),
                surface_specific_humidity=np.asarray(state["surface_specific_humidity"]),
                air_specific_humidity=np.asarray(state["specific_humidity"]),
                air_density=np.asarray(state["air_density"]),
            )
            diagnostics["surface_upward_sensible_heat_flux"] = fluxes["sensible_heat_flux"]
            diagnostics["surface_upward_latent_heat_flux"] = fluxes["latent_heat_flux"]
        return diagnostics
