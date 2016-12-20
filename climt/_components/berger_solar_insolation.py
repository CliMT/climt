#from .._core.base_components import Diagnostic
#from ._berger_solar_insolation import get_parameters

class BergerSolarInsolation(object):#Diagnostic):
    inputs = ('time',)
    diagnostics = (
        'obliquity_relative_to_mean_ecliptic',
        'precessional_parameter',
        'eccentricity',)

    def __call__(self, state):
        """
        Gets diagnostics from the passed model state.

        Args:
            state (dict): A model state dictionary.

        Returns:
            diagnostics (dict): A dicitonary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the time of the input state.

        Raises:
            KeyError: If a required quantity is missing from the state.
            InvalidStateException: If state is not a valid input for the
                Prognostic instance.
        """
        t = years_since_jan_1_1950_ad(state['time'])

        pass

def years_since_jan_1_1950_ad(dt):
    pass

if __name__ == '__main__':
    import pyximport; pyximport.install()
    import _berger_solar_insolation
