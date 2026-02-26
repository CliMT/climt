from sympl import TendencyComponent
import numpy as np
from .._core.backend import get_array_namespace, jit_compile, prange

# Area type mapping for JIT compatibility
# land: 0, land_ice: 1, sea: 2, sea_ice: 3
AREA_MAP = {"land": 0, "land_ice": 1, "sea": 2, "sea_ice": 3}

class SlabSurface(TendencyComponent):
    """
    Calculate the surface energy balance.

    This component assumes the surface is a slab of possibly
    varying heat capacity, and calculates the surface temperature.
    """

    input_properties = {
        "downwelling_longwave_flux_in_air": {
            "dims": ["*", "interface_levels"],
            "units": "W m^-2",
        },
        "downwelling_shortwave_flux_in_air": {
            "dims": ["*", "interface_levels"],
            "units": "W m^-2",
        },
        "upwelling_longwave_flux_in_air": {
            "dims": ["*", "interface_levels"],
            "units": "W m^-2",
        },
        "upwelling_shortwave_flux_in_air": {
            "dims": ["*", "interface_levels"],
            "units": "W m^-2",
        },
        "surface_upward_latent_heat_flux": {
            "dims": ["*"],
            "units": "W m^-2",
        },
        "surface_temperature": {
            "dims": ["*"],
            "units": "degK",
        },
        "surface_upward_sensible_heat_flux": {
            "dims": ["*"],
            "units": "W m^-2",
        },
        "surface_thermal_capacity": {
            "dims": ["*"],
            "units": "J kg^-1 degK^-1",
        },
        "surface_material_density": {
            "dims": ["*"],
            "units": "kg m^-3",
        },
        "upward_heat_flux_at_ground_level_in_soil": {
            "dims": ["*"],
            "units": "W m^-2",
        },
        "heat_flux_into_sea_water_due_to_sea_ice": {
            "dims": ["*"],
            "units": "W m^-2",
        },
        "area_type": {
            "dims": ["*"],
            "units": "dimensionless",
        },
        "soil_layer_thickness": {
            "dims": ["*"],
            "units": "m",
        },
        "ocean_mixed_layer_thickness": {
            "dims": ["*"],
            "units": "m",
        },
        "heat_capacity_of_soil": {
            "dims": ["*"],
            "units": "J kg^-1 degK^-1",
        },
        "sea_water_density": {
            "dims": ["*"],
            "units": "kg m^-3",
        },
    }

    tendency_properties = {
        "surface_temperature": {
            "dims": ["*"],
            "units": "degK s^-1",
        }
    }

    diagnostic_properties = {
        "depth_of_slab_surface": {
            "dims": ["*"],
            "units": "m",
        },
    }

    def array_call(self, state):
        sw_down_raw = state["downwelling_shortwave_flux_in_air"]
        xp = get_array_namespace(sw_down_raw)
        
        area_type_raw = state["area_type"]
        
        # Convert area_type to integer codes
        # land: 0, land_ice: 1, sea: 2, sea_ice: 3
        # Standard climt state uses strings/bytes for area_type.
        if hasattr(area_type_raw, 'astype'):
            area_type_str = area_type_raw.astype(str)
        else:
            # If it's a JAX array of objects/strings (unlikely but possible)
            area_type_str = np.asarray(area_type_raw).astype(str)
            
        area_type_code = np.zeros(area_type_str.shape, dtype=np.int32)
        for k, v in AREA_MAP.items():
            area_type_code[area_type_str == k] = v
            
        # For JAX, we need the codes as a JAX array
        if xp is not np:
            area_type_code = xp.array(area_type_code)

        # Flatten inputs for kernel
        def flat(x): return xp.reshape(x, (-1,))
        
        # Extract surface values (index -1 if vertical is first dim and TOA is 0)
        # We'll use -1 to be safe as it's the standard for surface in climt.
        lw_down_raw = state["downwelling_longwave_flux_in_air"]
        sw_up_raw = state["upwelling_shortwave_flux_in_air"]
        lw_up_raw = state["upwelling_longwave_flux_in_air"]
        
        def get_surface(arr):
            data = getattr(arr, 'data', arr)
            if xp is np:
                data = np.asarray(data)
            
            if data.ndim > 1:
                # sympl aligns ["*", "interface_levels"] -> vertical is last
                surf = data[..., 0]
            else:
                surf = data
                
            return flat(surf)
        
        sw_down = get_surface(sw_down_raw)
        lw_down = get_surface(lw_down_raw)
        sw_up = get_surface(sw_up_raw)
        lw_up = get_surface(lw_up_raw)
        
        lh = flat(getattr(state["surface_upward_latent_heat_flux"], 'data', state["surface_upward_latent_heat_flux"]))
        sh = flat(getattr(state["surface_upward_sensible_heat_flux"], 'data', state["surface_upward_sensible_heat_flux"]))
        
        up_heat_soil = flat(getattr(state["upward_heat_flux_at_ground_level_in_soil"], 'data', state["upward_heat_flux_at_ground_level_in_soil"]))
        heat_flux_sea_ice = flat(getattr(state["heat_flux_into_sea_water_due_to_sea_ice"], 'data', state["heat_flux_into_sea_water_due_to_sea_ice"]))
        
        sea_water_dens = flat(getattr(state["sea_water_density"], 'data', state["sea_water_density"]))
        surf_dens = flat(getattr(state["surface_material_density"], 'data', state["surface_material_density"]))
        
        heat_cap_soil = flat(getattr(state["heat_capacity_of_soil"], 'data', state["heat_capacity_of_soil"]))
        surf_therm_cap = flat(getattr(state["surface_thermal_capacity"], 'data', state["surface_thermal_capacity"]))
        
        ocean_mix_thick = flat(getattr(state["ocean_mixed_layer_thickness"], 'data', state["ocean_mixed_layer_thickness"]))
        soil_layer_thick = flat(getattr(state["soil_layer_thickness"], 'data', state["soil_layer_thickness"]))
        
        if xp is np:
            tend_ts, depth = _slab_surface_kernel_np(
                sw_down, lw_down, sw_up, lw_up, lh, sh,
                flat(area_type_code), up_heat_soil, heat_flux_sea_ice,
                sea_water_dens, surf_dens, heat_cap_soil, surf_therm_cap,
                ocean_mix_thick, soil_layer_thick
            )
        else:
            tend_ts, depth = _slab_surface_kernel_jax(
                sw_down, lw_down, sw_up, lw_up, lh, sh,
                flat(area_type_code), up_heat_soil, heat_flux_sea_ice,
                sea_water_dens, surf_dens, heat_cap_soil, surf_therm_cap,
                ocean_mix_thick, soil_layer_thick
            )
            
        return {"surface_temperature": xp.reshape(tend_ts, area_type_raw.shape)}, \
               {"depth_of_slab_surface": xp.reshape(depth, area_type_raw.shape)}

@jit_compile(backend=np)
def _slab_surface_kernel_np(sw_down, lw_down, sw_up, lw_up, lh, sh,
                           area_type, up_heat_soil, heat_flux_sea_ice,
                           sea_water_dens, surf_dens, heat_cap_soil, surf_therm_cap,
                           ocean_mix_thick, soil_layer_thick):
    
    ncol = area_type.size
    tend_ts = np.zeros(ncol)
    depth = np.zeros(ncol)
    
    for i in prange(ncol):
        net_heat_flux = sw_down[i] + lw_down[i] - sw_up[i] - lw_up[i] - sh[i] - lh[i]
        
        at = area_type[i]
        land_mask = (at == 0) or (at == 1)
        sea_mask = (at == 2) or (at == 3)
        land_ice_mask = (at == 1)
        sea_ice_mask = (at == 3)
        
        if land_ice_mask:
            net_heat_flux = -up_heat_soil[i]
        elif sea_ice_mask:
            net_heat_flux = heat_flux_sea_ice[i]
            
        if sea_mask:
            final_dens = sea_water_dens[i]
            d = ocean_mix_thick[i]
        else:
            final_dens = surf_dens[i]
            d = soil_layer_thick[i] if land_mask else 0.0
            
        if land_mask:
            final_therm_cap = heat_cap_soil[i]
        else:
            final_therm_cap = surf_therm_cap[i]
            
        depth[i] = d
        mass_slab = final_dens * d
        heat_cap_slab = mass_slab * final_therm_cap
        
        if heat_cap_slab != 0:
            val = net_heat_flux / heat_cap_slab
        else:
            val = 0.0
            
        if land_ice_mask or sea_ice_mask:
            val = 0.0
            
        tend_ts[i] = val
    
    return tend_ts, depth

def _slab_surface_kernel_jax(sw_down, lw_down, sw_up, lw_up, lh, sh,
                            area_type, up_heat_soil, heat_flux_sea_ice,
                            sea_water_dens, surf_dens, heat_cap_soil, surf_therm_cap,
                            ocean_mix_thick, soil_layer_thick):
    import jax.numpy as jnp
    
    net_heat_flux = sw_down + lw_down - sw_up - lw_up - sh - lh
    
    land_mask = (area_type == 0) | (area_type == 1)
    sea_mask = (area_type == 2) | (area_type == 3)
    land_ice_mask = (area_type == 1)
    sea_ice_mask = (area_type == 3)
    
    net_heat_flux = jnp.where(land_ice_mask, -up_heat_soil, net_heat_flux)
    net_heat_flux = jnp.where(sea_ice_mask, heat_flux_sea_ice, net_heat_flux)
    
    final_dens = jnp.where(sea_mask, sea_water_dens, surf_dens)
    final_therm_cap = jnp.where(land_mask, heat_cap_soil, surf_therm_cap)
    
    depth = jnp.where(sea_mask, ocean_mix_thick, jnp.where(land_mask, soil_layer_thick, 0.0))
    
    mass_slab = final_dens * depth
    heat_cap_slab = mass_slab * final_therm_cap
    
    tend_ts = net_heat_flux / heat_cap_slab
    tend_ts = jnp.where(land_ice_mask | sea_ice_mask, 0.0, tend_ts)
    
    return tend_ts, depth
