"""
Benchmark numba speedup for each optimized component.

Timings are wall-clock seconds for N iterations over a 8192-column,
30-level grid. The "Python" path calls the raw Python function via
.py_func on the @njit-decorated kernel (same code, no JIT).
"""
import time
import numpy as np
import sympl
import climt
from climt import get_grid, get_default_state
from datetime import timedelta

NCOL = 8192
NLEV = 30
TIMESTEP = timedelta(minutes=10)


def wall(fn, iters):
    fn()  # warm-up
    t0 = time.perf_counter()
    for _ in range(iters):
        fn()
    return time.perf_counter() - t0


def report(name, t_py, t_nb, iters, ncol):
    speedup = t_py / t_nb
    print(f"  Python : {t_py:.3f}s  ({t_py/iters*1000:.2f} ms/call)")
    print(f"  Numba  : {t_nb:.3f}s  ({t_nb/iters*1000:.2f} ms/call)")
    print(f"  Speedup: {speedup:.1f}x")


# ---------------------------------------------------------------------------
def bench_held_suarez():
    from climt import HeldSuarez
    from climt._components.held_suarez import _held_suarez_kernel_np

    ITERS = 50
    grid = get_grid(nx=NCOL, ny=1, nz=NLEV)
    sympl.set_backend(sympl.DataArrayBackend())
    comp = HeldSuarez()
    state = get_default_state([comp], grid_state=grid)
    state['eastward_wind'].values[:] = 10.0
    state['air_temperature'].values[:] = 300.0
    comp(state)  # trigger JIT compile

    print(f"\nHeldSuarez  ({NCOL} cols, {ITERS} iters)")
    t_nb = wall(lambda: comp(state), ITERS)
    t_py = wall(lambda: comp(state) if _held_suarez_kernel_np.py_func else None, ITERS)

    # Time pure Python kernel directly
    t = state['air_temperature'].values.reshape(NLEV, -1)
    u = state['eastward_wind'].values.reshape(NLEV, -1)
    v = state['northward_wind'].values.reshape(NLEV, -1)
    p = state['air_pressure'].values.reshape(NLEV, -1)
    ps = state['surface_air_pressure'].values.reshape(-1)
    lat = state['latitude'].values.reshape(-1)
    params = comp._params

    t_py = wall(lambda: _held_suarez_kernel_np.py_func(u, v, t, p, ps, lat, params), ITERS)
    t_nb2 = wall(lambda: _held_suarez_kernel_np(u, v, t, p, ps, lat, params), ITERS)
    report("HeldSuarez", t_py, t_nb2, ITERS, NCOL)


def bench_gray_radiation():
    from climt import GrayLongwaveRadiation, Frierson06LongwaveOpticalDepth
    from climt._components.radiation import _gray_lw_kernel_np, _frierson_tau_kernel_np

    ITERS = 50
    grid = get_grid(nx=NCOL, ny=1, nz=NLEV)
    sympl.set_backend(sympl.DataArrayBackend())
    tau_comp = Frierson06LongwaveOpticalDepth()
    lw_comp = GrayLongwaveRadiation()
    state = get_default_state([tau_comp, lw_comp], grid_state=grid)
    state.update(tau_comp(state))
    lw_comp(state)  # trigger JIT

    from sympl import get_constant
    sigma = float(get_constant("stefan_boltzmann_constant", "W/m^2/K^4"))
    g     = float(get_constant("gravitational_acceleration", "m/s^2"))
    cpd   = float(get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/K"))

    t_flat    = np.ascontiguousarray(state['air_temperature'].values.reshape(NLEV, -1))
    tau_flat  = np.ascontiguousarray(state['longwave_optical_depth_on_interface_levels'].values.reshape(NLEV+1, -1))
    ts_flat   = np.ascontiguousarray(state['surface_temperature'].values.reshape(-1))
    pint_flat = np.ascontiguousarray(state['air_pressure_on_interface_levels'].values.reshape(NLEV+1, -1))

    print(f"\nGrayLongwaveRadiation  ({NCOL} cols, {ITERS} iters)")
    t_py = wall(lambda: _gray_lw_kernel_np.py_func(t_flat, pint_flat, ts_flat, tau_flat, sigma), ITERS)
    t_nb = wall(lambda: _gray_lw_kernel_np(t_flat, pint_flat, ts_flat, tau_flat, sigma), ITERS)
    report("GrayLW", t_py, t_nb, ITERS, NCOL)

    lat_flat = np.ascontiguousarray(state['latitude'].values.reshape(-1))
    ps_flat  = np.ascontiguousarray(state['surface_air_pressure'].values.reshape(-1))
    tau0e, tau0p, fl = tau_comp._tau0e, tau_comp._tau0p, tau_comp._fl

    print(f"\nFrierson06LongwaveOpticalDepth  ({NCOL} cols, {ITERS} iters)")
    t_py = wall(lambda: _frierson_tau_kernel_np.py_func(lat_flat, pint_flat, ps_flat, tau0e, tau0p, fl), ITERS)
    t_nb = wall(lambda: _frierson_tau_kernel_np(lat_flat, pint_flat, ps_flat, tau0e, tau0p, fl), ITERS)
    report("Frierson tau", t_py, t_nb, ITERS, NCOL)


def bench_gsc():
    from climt import GridScaleCondensation
    from climt._components.grid_scale_condensation import _gsc_kernel_np

    ITERS = 20
    grid = get_grid(nx=NCOL, ny=1, nz=NLEV)
    sympl.set_backend(sympl.DataArrayBackend())
    comp = GridScaleCondensation()
    state = get_default_state([comp], grid_state=grid)
    state['specific_humidity'].values[:] = 0.02
    comp(state, TIMESTEP)  # trigger JIT

    t_flat    = np.ascontiguousarray(state['air_temperature'].values.reshape(NLEV, -1))
    q_flat    = np.ascontiguousarray(state['specific_humidity'].values.reshape(NLEV, -1))
    p_flat    = np.ascontiguousarray(state['air_pressure'].values.reshape(NLEV, -1))
    pint_flat = np.ascontiguousarray(state['air_pressure_on_interface_levels'].values.reshape(NLEV+1, -1))
    params = comp._params

    print(f"\nGridScaleCondensation  ({NCOL} cols, {ITERS} iters)")
    t_py = wall(lambda: _gsc_kernel_np.py_func(t_flat, q_flat, p_flat, pint_flat, params), ITERS)
    t_nb = wall(lambda: _gsc_kernel_np(t_flat, q_flat, p_flat, pint_flat, params), ITERS)
    report("GSC", t_py, t_nb, ITERS, NCOL)


def bench_dry_convection():
    from climt import DryConvectiveAdjustment
    from climt._components.dry_convection.component import _dry_adj_kernel_np, DryAdjParams
    from sympl import get_constant

    ITERS = 10
    grid = get_grid(nx=NCOL, ny=1, nz=NLEV)
    sympl.set_backend(sympl.DataArrayBackend())
    comp = DryConvectiveAdjustment()
    state = get_default_state([comp], grid_state=grid)
    state['air_temperature'].values[:5] += 20.0
    comp(state, TIMESTEP)  # trigger JIT

    t_flat    = np.ascontiguousarray(state['air_temperature'].values.reshape(NLEV, -1))
    q_flat    = np.ascontiguousarray(state['specific_humidity'].values.reshape(NLEV, -1))
    p_flat    = np.ascontiguousarray(state['air_pressure'].values.reshape(NLEV, -1))
    pint_flat = np.ascontiguousarray(state['air_pressure_on_interface_levels'].values.reshape(NLEV+1, -1))
    params = DryAdjParams(
        Cpd=float(get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK")),
        Cvap=float(get_constant("heat_capacity_of_vapor_phase", "J/kg/K")),
        Rdair=float(get_constant("gas_constant_of_dry_air", "J/kg/degK")),
        Pref=float(get_constant("reference_air_pressure", "Pa")),
        Rv=float(get_constant("gas_constant_of_vapor_phase", "J/kg/K"))
    )

    print(f"\nDryConvectiveAdjustment  ({NCOL} cols, {ITERS} iters)")
    t_py = wall(lambda: _dry_adj_kernel_np.py_func(t_flat, q_flat, p_flat, pint_flat, params), ITERS)
    t_nb = wall(lambda: _dry_adj_kernel_np(t_flat, q_flat, p_flat, pint_flat, params), ITERS)
    report("DryConv", t_py, t_nb, ITERS, NCOL)


def bench_berger():
    from climt import BergerSolarInsolation
    from climt._components.berger_solar_insolation import _get_solar_parameters_np
    from datetime import datetime

    ITERS = 100
    grid = get_grid(nx=NCOL, ny=1, nz=NLEV)
    sympl.set_backend(sympl.DataArrayBackend())
    comp = BergerSolarInsolation()
    state = get_default_state([comp], grid_state=grid)
    state['time'] = datetime(2000, 6, 21)
    comp(state)  # trigger JIT

    year = state['time'].year
    lambda_m0, ecc, omega, obl = comp._orbital_parameters[year]
    from climt._components.berger_solar_insolation import years_since_vernal_equinox
    yve = years_since_vernal_equinox(state['time'])
    frac_day = (state['time'].hour + state['time'].minute / 60.) / 24.
    lat_flat = np.ascontiguousarray(state['latitude'].values.reshape(-1).astype(np.float64))
    lon_flat = np.ascontiguousarray(state['longitude'].values.reshape(-1).astype(np.float64))
    solar_const = float(climt.get_default_state([comp], grid_state=grid)['solar_constant'].values.flat[0]) if 'solar_constant' in state else 1367.0

    from sympl import get_constant
    solar_const = float(get_constant('stellar_irradiance', 'W/m^2'))

    print(f"\nBergerSolarInsolation  ({NCOL} cols, {ITERS} iters)")
    t_py = wall(lambda: _get_solar_parameters_np.py_func(lambda_m0, ecc, omega, obl, yve, frac_day, lat_flat, lon_flat, solar_const), ITERS)
    t_nb = wall(lambda: _get_solar_parameters_np(lambda_m0, ecc, omega, obl, yve, frac_day, lat_flat, lon_flat, solar_const), ITERS)
    report("Berger", t_py, t_nb, ITERS, NCOL)


def bench_slab_surface():
    from climt import SlabSurface
    from climt._components.slab_surface import _slab_surface_kernel_np

    ITERS = 50
    grid = get_grid(nx=NCOL, ny=1, nz=NLEV)
    sympl.set_backend(sympl.DataArrayBackend())
    comp = SlabSurface()
    state = get_default_state([comp], grid_state=grid)
    comp(state)  # trigger JIT

    def _flat(key):
        return np.ascontiguousarray(state[key].values.reshape(-1))

    def _surf(key):
        v = state[key].values
        return np.ascontiguousarray((v[0] if v.ndim > 2 else v).reshape(-1))

    area_type_raw = state['area_type'].values.reshape(-1)
    area_map = {'land': 0, 'land_ice': 1, 'sea': 2, 'sea_ice': 3}
    area_code = np.array([area_map.get(str(a), 2) for a in area_type_raw], dtype=np.float64)

    sw_down  = _surf('downwelling_shortwave_flux_in_air')
    lw_down  = _surf('downwelling_longwave_flux_in_air')
    sw_up    = _surf('upwelling_shortwave_flux_in_air')
    lw_up    = _surf('upwelling_longwave_flux_in_air')
    lh       = _flat('surface_upward_latent_heat_flux')
    sh       = _flat('surface_upward_sensible_heat_flux')
    up_soil  = _flat('upward_heat_flux_at_ground_level_in_soil')
    hf_ice   = _flat('heat_flux_into_sea_water_due_to_sea_ice')
    sw_dens  = _flat('sea_water_density')
    sf_dens  = _flat('surface_material_density')
    hc_soil  = _flat('heat_capacity_of_soil')
    stc      = _flat('surface_thermal_capacity')
    omt      = _flat('ocean_mixed_layer_thickness')
    slt      = _flat('soil_layer_thickness')

    print(f"\nSlabSurface  ({NCOL} cols, {ITERS} iters)")
    t_py = wall(lambda: _slab_surface_kernel_np.py_func(sw_down, lw_down, sw_up, lw_up, lh, sh, area_code, up_soil, hf_ice, sw_dens, sf_dens, hc_soil, stc, omt, slt), ITERS)
    t_nb = wall(lambda: _slab_surface_kernel_np(sw_down, lw_down, sw_up, lw_up, lh, sh, area_code, up_soil, hf_ice, sw_dens, sf_dens, hc_soil, stc, omt, slt), ITERS)
    report("SlabSurface", t_py, t_nb, ITERS, NCOL)


def bench_instellation():
    from climt import Instellation
    from climt._components.instellation.component import (
        _instellation_kernel_np, days_from_2000, fractional_day
    )
    from datetime import datetime

    ITERS = 100
    grid = get_grid(nx=NCOL, ny=1, nz=NLEV)
    sympl.set_backend(sympl.DataArrayBackend())
    comp = Instellation()
    state = get_default_state([comp], grid_state=grid)
    state['time'] = datetime(2000, 6, 21, 12)
    comp(state)  # trigger JIT

    lat_flat = np.ascontiguousarray(state['latitude'].values.reshape(-1).astype(np.float64))
    lon_flat = np.ascontiguousarray(state['longitude'].values.reshape(-1).astype(np.float64))

    t = state['time']
    julian_centuries = days_from_2000(t) / 36525.0
    frac_day = fractional_day(t)

    print(f"\nInstellation  ({NCOL} cols, {ITERS} iters)")
    args = (lat_flat, lon_flat, julian_centuries, frac_day)
    t_py = wall(lambda: _instellation_kernel_np.py_func(*args), ITERS)
    t_nb = wall(lambda: _instellation_kernel_np(*args), ITERS)
    report("Instellation", t_py, t_nb, ITERS, NCOL)


if __name__ == "__main__":
    print(f"Numba speedup benchmark — {NCOL} columns, {NLEV} levels")
    print("=" * 55)
    bench_held_suarez()
    bench_gray_radiation()
    bench_gsc()
    bench_dry_convection()
    bench_berger()
    bench_slab_surface()
    bench_instellation()
    print()
