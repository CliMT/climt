from climt import ClimtImplicit
from sympl import replace_none_with_default

class IceSheet(ClimtImplicit):
    """
    1-d ice energy balance model.
    """

    _climt_inputs = {
        'downwelling_longwave_flux_in_air': 'W m^-2',
        'downwelling_shortwave_flux_in_air': 'W m^-2',
        'upwelling_longwave_flux_in_air': 'W m^-2',
        'upwelling_shortwave_flux_in_air': 'W m^-2',
        'surface_upward_latent_heat_flux': 'W m^-2',
        'surface_upward_sensible_heat_flux': 'W m^-2',
        'land_ice_thickness': 'm',
        'sea_ice_thickness': 'm',
        'surface_snow_thickness': 'm',
        'surface_albedo': 'dimensionless',
        'area_type': 'dimensionless'
    }

    _climt_outputs = {
        'land_ice_thickness': 'm',
        'sea_ice_thickness': 'm',
        'surface_snow_thickness': 'm'
    }

    _climt_diagnostics = {}

    def __init__(
        self,
        num_layers=50,
        thermal_conductivity_of_ice=2.22,
        thermal_conductivity_of_snow=0.2,
        density_of_ice=916.7,
        density_of_snow=100.,
        specific_heat_of_ice=2108.,
        specific_heat_of_snow=2108.,
        latent_heat_of_fusion=333550.,
        melting_point_of_ice=273.):
        """
        
        Args:
            num_layers(int):
                The number of layers to use to simulate the ice sheet.
                
            thermal_conductivity_of_ice (float, optional):
                The thermal conductivity of ice in :math:`W m^{-1} K^{-1}`.
            
            thermal_conductivity_of_snow (float, optional):
                The thermal conductivity of ice in :math:`W m^{-1} K^{-1}`.
            
            density_of_ice (float, optional):
                The density of ice in :math:`kg m^{-3}`.
            
            density_of_snow (float, optional):
                The density of snow in :math:`kg m^{-3}`.
            
            specific_heat_of_ice (float, optional):
                The specific heat of ice in :math:`J kg^{-1}`.
            
            specific_heat_of_snow (float, optional):
                The specific heat of snow in :math:`J kg^{-1}`.
            
            latent_heat_of_fusion (float, optional):
                The latent heat of fusion of ice in :math:`J kg^{-1}`
            
            melting_point_of_ice (float, optional):
                The melting point of ice in :math:`K`.
        """

        self._Kice = replace_none_with_default(
            'thermal_conductivity_of_ice', thermal_conductivity_of_ice)

        self._Ksnow = replace_none_with_default(
            'thermal_conductivity_of_snow', thermal_conductivity_of_snow)

        self._rho_ice = replace_none_with_default(
            'density_of_ice', density_of_ice)

        self._C_ice = replace_none_with_default(
            'specific_heat_of_ice', specific_heat_of_ice)

        self._rho_snow = replace_none_with_default(
            'density_of_snow', density_of_snow)

        self._C_snow = replace_none_with_default(
            'specific_heat_of_snow', specific_heat_of_snow)

        self._Lf = replace_none_with_default(
            'latent_heat_of_fusion_of_ice', latent_heat_of_fusion)

        self._temp_melt = replace_none_with_default(
            'melting_point_of_ice', melting_point_of_ice)

    def __call__(self, state, time_step):
        """
        Calculate new ice sheet height.
        
        Args:
            state (dict):
                The state dictionary.
                
            time_step (timedelta):
                The model timestep.

        Returns:
            new_state(dict), diagnostics(dict):
                * The new state calculated
                * Any diagnostics

        """

        raw_arrays = self.get_numpy_arrays_from_state('_climt_inputs', state)

        num_lons, num_lats = raw_arrays['area_type'].shape

        for lon in range(num_lons):
            for lat in range(num_lats):

                    area_type = raw_arrays['area_type'][lon, lat].astype(str)
                    if area_type == 'land_ice':
                        total_height = raw_arrays['land_ice_thickness'] + raw_arrays['surface_snow_thickness']
