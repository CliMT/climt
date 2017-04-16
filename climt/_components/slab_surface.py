from climt._core.climt_components import ClimtPrognostic


class SlabSurface(ClimtPrognostic):
    """
    Calculate the surface energy balance.
    
    This component assumes the surface is a slab of possibly
    varying heat capacity, and calculates the surface temperature.
    """

    _climt_inputs = {
        'downwelling_longwave_flux_in_air': 'W m^-2',
        'downwelling_shortwave_flux_in_air': 'W m^-2',
        'upwelling_longwave_flux_in_air': 'W m^-2',
        'upwelling_shortwave_flux_in_air': 'W m^-2',
        'surface_upward_latent_heat_flux': 'W m^-2',
        'surface_upward_sensible_heat_flux': 'W m^-2',
        'surface_thermal_capacity': 'J kg^-1 degK^-1',
        'depth_slab_surface': 'm',
        'density_surface_material': 'kg m^-3'
    }

    _climt_tendencies = {
        'surface_temperature': 'degK s^-1',
    }

    _climt_diagnostics = {}

    quantity_descriptions = {
        'surface_thermal_capacity': {
            'dims': ['x', 'y'],
            'units': 'J kg^-1 degK^-1',
            'default_value': 4.1813e3
        },
        'depth_slab_surface': {
            'dims': ['x', 'y'],
            'units': 'm',
            'default_value': 50.
        },
        'density_surface_material': {
        'dims': ['x', 'y'],
        'units': 'kg m^-3',
        'default_value': 1000.
        }
    }
    """
    Quantities used for surface heat capacity calculation.
    """

    def __init__(self):
        """
        Initialise slab surface
        """
        return

    def __call__(self, state):
        """
        Calculate surface temperature.
        
        Args:
            state (dict):
                The state dictionary
                
        Returns:
            tendencies (dict), diagnostics(dict):
                * The surface temperature tendency
                * Any diagnostics

        """

        raw_arrays = self.get_numpy_arrays_from_state('_climt_inputs', state)

        net_heat_flux = (raw_arrays['downwelling_shortwave_flux_in_air'][:, :, 0]+
                         raw_arrays['downwelling_longwave_flux_in_air'][:, :, 0]  -
                         raw_arrays['upwelling_shortwave_flux_in_air'][:, :, 0] -
                         raw_arrays['upwelling_longwave_flux_in_air'][:, :, 0] -
                         raw_arrays['surface_upward_sensible_heat_flux'] -
                         raw_arrays['surface_upward_latent_heat_flux'])

        mass_surface_slab = raw_arrays['density_surface_material']*raw_arrays['depth_slab_surface']
        heat_capacity_surface = mass_surface_slab*raw_arrays['surface_thermal_capacity']

        tend_dict = self.create_state_dict_for('_climt_tendencies', state)
        diag_dict = self.create_state_dict_for('_climt_diagnostics', state)

        tend_dict['surface_temperature'].values[:] = net_heat_flux/heat_capacity_surface

        return tend_dict, diag_dict
