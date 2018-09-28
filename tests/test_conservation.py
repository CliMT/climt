import climt
from sympl import (
    get_constant, AdamsBashforth,
    TendencyComponent, TendencyStepper)
from datetime import timedelta
import numpy as np


class Conservation(object):

    def get_model_state(self, component):

        state = climt.get_default_state([component])
        return self.modify_state(state)

    def modify_state(self, state):
        return state

    def get_component_instance(self):
        pass

    def get_quantity_forcing(self, state):
        pass

    def get_quantity_amount(self, state):
        pass

    def get_steppable_component(self):

        component = self.get_component_instance()

        if isinstance(component, TendencyComponent):
            return AdamsBashforth(component)
        else:
            return component

    def get_new_state_and_diagnostics(self, state,
                                      component, time_step):

        if isinstance(component, TendencyStepper):
            diag, state = component(state, time_step)
            state.update(diag)

            return state
        else:
            diag, new_state = component(state, time_step)
            state.update(new_state)
            state.update(diag)

            return state

    def test_quantity_is_conserved(self):

        component = self.get_steppable_component()
        state = self.get_model_state(component)
        time_step = timedelta(seconds=1)

        first_state = self.get_new_state_and_diagnostics(state,
                                                         component, time_step)

        old_amount = self.get_quantity_amount(first_state)

        new_state = self.get_new_state_and_diagnostics(first_state,
                                                       component, time_step)

        new_amount = self.get_quantity_amount(new_state)

        forcing_amount = self.get_quantity_forcing(new_state)*time_step.total_seconds()

        assert np.isclose(new_amount - old_amount, forcing_amount,
                          rtol=0, atol=1e-3)


def get_pressure_thickness(state):
    return (state['air_pressure_on_interface_levels'][:-1] -
            state['air_pressure_on_interface_levels'][1:]).rename(
                dict(interface_levels='mid_levels'))


def vertical_integral(state, quantity):
    g = get_constant('gravitational_acceleration', 'm/s^2')
    dp = get_pressure_thickness(state)

    return (quantity*dp/g).sum().values


def get_moist_enthalpy(state):
    Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/degK')
    Lv = get_constant('latent_heat_of_condensation', 'J/kg')

    dry_enthalpy = vertical_integral(state, Cpd*state['air_temperature'])

    moisture_enthalpy = 0.0
    if 'specific_humidity' in state:
        moisture_enthalpy = vertical_integral(state, Lv*state['specific_humidity'])

    return dry_enthalpy + moisture_enthalpy


def heat_capacity_including_condensible(q):

    Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/degK')
    Cvap = get_constant('heat_capacity_of_vapor_phase', 'J/kg/K')

    return Cpd*(1-q) + Cvap*q


def get_surface_forcing(state):

    surf_forcing = 0

    if 'upwelling_shortwave_flux_in_air' in state:
        surf_forcing += state['upwelling_shortwave_flux_in_air'][0].values

    if 'upwelling_longwave_flux_in_air' in state:
        surf_forcing += state['upwelling_longwave_flux_in_air'][0].values

    if 'downwelling_shortwave_flux_in_air' in state:
        surf_forcing -= state['downwelling_shortwave_flux_in_air'][0].values

    if 'downwelling_longwave_flux_in_air' in state:
        surf_forcing -= state['downwelling_longwave_flux_in_air'][0].values

    if 'surface_upward_sensible_heat_flux' in state:
        surf_forcing += state['surface_upward_sensible_heat_flux'].values

    if 'surface_upward_latent_heat_flux' in state:
        surf_forcing += state['surface_upward_latent_heat_flux'].values

    return surf_forcing


def get_top_of_atmosphere_forcing(state):

    toa_forcing = 0

    if 'downwelling_shortwave_flux_in_air' in state:
        toa_forcing += state['downwelling_shortwave_flux_in_air'][-1].values

    if 'upwelling_shortwave_flux_in_air' in state:
        toa_forcing -= state['upwelling_shortwave_flux_in_air'][-1].values

    if 'upwelling_longwave_flux_in_air' in state:
        toa_forcing -= state['upwelling_longwave_flux_in_air'][-1].values

    return toa_forcing


class AtmosphereEnergyConservation(Conservation):

    def get_quantity_amount(self, state):
        return get_moist_enthalpy(state)

    def get_quantity_forcing(self, state):
        return get_surface_forcing(state) + get_top_of_atmosphere_forcing(state)


class AtmosphereTracerConservation(Conservation):

    def get_quantity_amount(self, state):
        return vertical_integral(state, state[self.tracer_name])


class SurfaceEnergyConservation(Conservation):

    def get_quantity_forcing(self, state):
        return -get_surface_forcing(state)


class SlabSurfaceConservation(SurfaceEnergyConservation):

    def get_component_instance(self):
        return climt.SlabSurface()

    def get_quantity_amount(self, state):
        C = state['surface_thermal_capacity'].values
        T = state['surface_temperature'].values
        mass = state['sea_water_density'].values*state['ocean_mixed_layer_thickness'].values

        return mass*C*T


'''
class SeaIceConservation(SurfaceEnergyConservation):

    def get_component_instance(self):
        return climt.IceSheet()

    def get_quantity_amount(self, state):
        C_ice = get_constant('heat_capacity_of_solid_phase_as_ice', 'J/kg/degK')
        C_snow = get_constant('heat_capacity_of_solid_phase_as_snow', 'J/kg/degK')
        rho_ice = get_constant('density_of_solid_phase_as_ice', 'kg/m^3')
        rho_snow = get_constant('density_of_solid_phase_as_snow', 'kg/m^3')

        height_snow_ice = state['sea_ice_thickness'].values + state['surface_snow_thickness'].values
        snow_height_fraction = state['surface_snow_thickness'].values / height_snow_ice

        temp_profile = state['snow_and_ice_temperature'].values
        num_layers = temp_profile.shape[0]
        dz = float(height_snow_ice / num_layers)

        snow_level = int((1 - snow_height_fraction)*num_layers) - 1
        levels = np.arange(num_layers - 1)

        # Create vertically varying profiles
        rho_snow_ice = rho_ice*np.ones(num_layers - 1)
        rho_snow_ice[levels > snow_level] = rho_snow

        heat_capacity_snow_ice = C_ice*np.ones(num_layers - 1)
        heat_capacity_snow_ice[levels > snow_level] = C_snow

        mass_snow_ice = rho_snow_ice*dz

        temp_mean = (temp_profile[1:] + temp_profile[:-1])/2.

        heat_content = np.sum(mass_snow_ice*heat_capacity_snow_ice*temp_mean)

        return heat_content

    def get_quantity_forcing(self, state):

        return -state['heat_flux_into_sea_water_due_to_sea_ice'].values +\
            state['surface_downward_heat_flux_in_sea_ice'].values
'''

#####################
# Start Actual Tests
#####################


class TestRRTMGLongwaveConservation(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.RRTMGLongwave()


class TestRRTMGShortwaveConservation(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.RRTMGShortwave()


class TestRRTMGShortwaveConservationWithClouds(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.RRTMGShortwave()

    def modify_state(self, state):
        state['mass_content_of_cloud_liquid_water_in_atmosphere_layer'].loc[dict(mid_levels=slice(4, 8))] = 0.03
        state['cloud_area_fraction_in_atmosphere_layer'].loc[dict(mid_levels=slice(4, 8))] = 1.
        return state


class TestRRTMGShortwaveConservationWithCloudsAndMCICA(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.RRTMGShortwave(mcica=True)

    def modify_state(self, state):
        state['mass_content_of_cloud_liquid_water_in_atmosphere_layer'].loc[dict(mid_levels=slice(4, 8))] = 0.03
        state['cloud_area_fraction_in_atmosphere_layer'].loc[dict(mid_levels=slice(4, 8))] = 0.5
        return state


class TestRRTMGLongwaveConservationWithClouds(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.RRTMGLongwave()

    def modify_state(self, state):
        state['mass_content_of_cloud_liquid_water_in_atmosphere_layer'].loc[dict(mid_levels=slice(4, 8))] = 0.03
        state['cloud_area_fraction_in_atmosphere_layer'].loc[dict(mid_levels=slice(4, 8))] = 1.
        return state


class TestSimplePhysicsDryConservation(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.SimplePhysics(boundary_layer=False,
                                   use_external_surface_specific_humidity=True)

    def modify_state(self, state):
        state['eastward_wind'].values[:] = 3.
        return state


class TestSimplePhysicsConservation(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.SimplePhysics(boundary_layer=False)

    def modify_state(self, state):
        state['eastward_wind'].values[:] = 3.
        return state


class TestDryConvectionConservation(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.DryConvectiveAdjustment()

    def modify_state(self, state):
        unstable_level = 5
        state['air_temperature'][:unstable_level] += 10
        state['specific_humidity'][:unstable_level] = 0.05
        return state

    def get_quantity_amount(self, state):
        return vertical_integral(
            state,
            heat_capacity_including_condensible(state['specific_humidity']).values *
            state['air_temperature'].values)

    def get_quantity_forcing(self, state):
        return 0


class TestDryConvectionCondensibleConservation(AtmosphereTracerConservation):

    def get_component_instance(self):
        self.tracer_name = 'specific_humidity'
        return climt.DryConvectiveAdjustment()

    def modify_state(self, state):
        unstable_level = 5
        state['air_temperature'][:unstable_level] += 10
        state['specific_humidity'][:unstable_level] = 0.05
        return state

    def get_quantity_amount(self, state):
        return vertical_integral(
            state, state['specific_humidity'])

    def get_quantity_forcing(self, state):
        return 0


class TestSlabSurfaceOnlySensibleHeat(SlabSurfaceConservation):

    def modify_state(self, state):
        state['surface_upward_sensible_heat_flux'].values = 10.
        state['ocean_mixed_layer_thickness'].values = 1.

        return state


class TestSlabSurfaceOnlyLatentHeat(SlabSurfaceConservation):

    def modify_state(self, state):
        state['surface_upward_latent_heat_flux'].values = 40.
        state['ocean_mixed_layer_thickness'].values = 1.

        return state


class TestSlabSurfaceOnlyRadiative(SlabSurfaceConservation):

    def modify_state(self, state):
        state['upwelling_shortwave_flux_in_air'].values[:] = 40.
        state['upwelling_longwave_flux_in_air'].values[:] = 40.
        state['downwelling_shortwave_flux_in_air'].values[:] = 40.
        state['downwelling_longwave_flux_in_air'].values[:] = 40.
        state['ocean_mixed_layer_thickness'].values = 1.

        return state


'''
class TestSeaIceOnlySensibleHeatNoSnow(SeaIceConservation):

    def modify_state(self, state):
        state['surface_upward_sensible_heat_flux'].values = 10.
        state['snow_and_ice_temperature'].values[:] = 273.
        state['sea_ice_thickness'].values = 5.
        state['area_type'].values = 'sea_ice'
        return state
'''
