import climt
from sympl import (
    get_constant, AdamsBashforth,
    TendencyComponent, TendencyStepper)
from datetime import timedelta
import numpy as np
import logging

logger = logging.getLogger('Conservation')


class Conservation(object):

    def get_model_state(self, component):

        return climt.get_default_state([component])

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

        old_amount = self.get_quantity_amount(state)

        time_step = timedelta(seconds=1)
        new_state = self.get_new_state_and_diagnostics(state,
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


class TestRRTMGLongwaveConservation(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.RRTMGLongwave()


class TestRRTMGShortwaveConservation(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.RRTMGShortwave()


class TestSimplePhysicsDryConservation(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.SimplePhysics(boundary_layer=False,
                                   use_external_surface_specific_humidity=True)

    def get_model_state(self, component):
        state = climt.get_default_state([component])
        state['eastward_wind'].values[:] = 3.
        return state


class TestSimplePhysicsConservation(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.SimplePhysics(boundary_layer=False)

    def get_model_state(self, component):
        state = climt.get_default_state([component])
        state['eastward_wind'].values[:] = 3.
        return state


class TestDryConvectionConservation(AtmosphereEnergyConservation):

    def get_component_instance(self):
        return climt.DryConvectiveAdjustment()

    def get_model_state(self, component):
        unstable_level = 5

        state = climt.get_default_state([component], grid_state=climt.get_grid(nz=35))
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
        return climt.DryConvectiveAdjustment()

    def get_model_state(self, component):
        unstable_level = 5

        state = climt.get_default_state([component], grid_state=climt.get_grid(nz=35))
        state['air_temperature'][:unstable_level] += 10
        state['specific_humidity'][:unstable_level] = 0.05
        return state

    def get_quantity_amount(self, state):
        return vertical_integral(
            state, state['specific_humidity'])

    def get_quantity_forcing(self, state):
        return 0

