import abc
from sympl import (Implicit, Diagnostic, Prognostic,
                   TimeStepper, PrognosticComposite,
                   UpdateFrequencyWrapper, ScalingWrapper)
from datetime import timedelta
from sympl import get_numpy_array
from .initialization import climt_quantity_descriptions, get_default_values
import numpy as np
import copy


class ArrayHandler(object):
    """
    Provide array handling functionality.

    This class provides methods to move back and forth from
    the :py:class:`~sympl.DataArray` and :py:class:`~numpy.array`
    representations of data in CliMT. The former is useful for human
    interaction, while the latter is required to pass onto cython or
    fortran extensions.

    Methods:
    """

    __metaclass__ = abc.ABCMeta

    def create_properties_dict(self, description):
        """

        Create properties dictionary using _climt_inputs/outputs/etc., dictionaries

        Args:

            description (dict):
                The climt component descriptions in _climt_inputs, etc.,

        Returns:
            properties_description (dict):
                An input/output/etc., properties dict conforming to Sympl requirements.
        """

        property_dict = {}

        for quantity in description.keys():

            # Use version in quantity descriptions if present
            if hasattr(self, 'quantity_descriptions') and quantity in self.quantity_descriptions.keys():
                property_dict[quantity] = copy.deepcopy(self.quantity_descriptions[quantity])
            else:
                # Use default version if not
                property_dict[quantity] = copy.deepcopy(climt_quantity_descriptions[quantity])

                # Update units
                property_dict[quantity]['units'] = copy.deepcopy(description[quantity])

        return property_dict

    def get_numpy_arrays_from_state(self, attribute, state, memory_layout='fortran'):
        """

        Extract inputs as numpy arrays from state.

        Returns arrays with dimensions (x,y,z) in the same order as specified in
        :code:`component` (first priority) or :code:`_quantity_descriptions` in
        :code:`initialization.py` (second priority).

        Args:

            attribute (string):
                The attribute (:code:`inputs`, :code:`tendencies`, :code:`outputs` or
                :code:`diagnostics`) for which to return numpy arrays.

            state (dict):
                The state dictionary.

            memory_layout (string, optional):
                String which is either :code:`'fortran'` or :code:`'c'`. This specifies
                the memory layout which the component expects the arrays in. If the arrays
                in :code:`state` are not in this memory layout and/or not memory aligned,
                a copy will be made.

        Returns:

            array_dict (dict):
                dictionary whose values are numpy arrays corresponding
                to the input quantities specified in :code:`component`. The returned arrays will
                be in the units specified in :code:`attribute`.

        Raises:

            NotImplementedError:
                If the component's :code:`inputs` attribute is not a dictionary.

            ValueError:
                If the :code:`memory_layout` argument is neither "fortran" nor "c".

        """

        quantities_to_extract = self.check_if_sane_and_return_attribute(attribute)

        if memory_layout not in ['fortran', 'c']:
            raise ValueError(
                'memory_layout can be either fortran or c')

        array_dict = {}

        for quantity in quantities_to_extract.keys():

            dims = self.get_dimensions_for(quantity)
            units = quantities_to_extract[quantity]

            new_array = get_array_from_state(
                state,
                quantity,
                units,
                dims)

            if memory_layout is 'fortran' and not new_array.flags['FARRAY']:
                new_array = np.asfortranarray(new_array)
            elif memory_layout is 'c' and not new_array.flags['CARRAY']:
                new_array = np.ascontiguousarray(new_array)

            array_dict[quantity] = new_array

        return array_dict

    def get_dimensions_for(self, quantity_name):

        if hasattr(self, 'quantity_descriptions'):
            if quantity_name in self.quantity_descriptions:
                return self.quantity_descriptions[quantity_name]['dims']

        if quantity_name in climt_quantity_descriptions:
            return climt_quantity_descriptions[quantity_name]['dims']

        # Should never come here.
        raise IndexError(
            '{} not described either by the component or by CliMT!'.format(quantity_name))

    def create_state_dict_for(self, attribute, state):
        """
        Create dictionaries to return to caller.

        Use quantities in :code:`component.attribute` to create a
        dictionary of DataArrays which is returned by the component
        to the caller.

        Args:

            self (Prognostic, Implicit, Diagnostic, ImplicitPrognostic, TimeStepper):
                component for which the output dictionary is required.

            attribute (basestring):
                The attribute of the component which should be used to create the
                dictionary. Typically, one of :code:`inputs`, :code:`tendencies`, :code:`outputs` or
                :code:`diagnostics`.

            state (dict):
                The state dictionary that was passed in to the component

        Returns:

            output_dict (dict):
                The dictionary whose keys are labels from :code:`component.attribute` and
                values are the appropriate DataArrays.

        """

        quantities_to_extract = self.check_if_sane_and_return_attribute(attribute)
        # quantities_to_extract is a dictionary whose keys are quantity names
        # and values are units produced by the code. Note that if this is a
        # tendency term, the final units returned to the caller must be in per second,
        # since the TimeStepper requires quantities in per seconds.

        output_state = {}
        for quantity in quantities_to_extract.keys():
            description = copy.deepcopy(climt_quantity_descriptions)

            if hasattr(self, 'quantity_descriptions'):
                if quantity in self.quantity_descriptions:
                    description[quantity] = self.quantity_descriptions[quantity]

            additional_dimensions = {}
            for dimension in description[quantity]['dims']:
                if dimension not in ['x', 'y', 'mid_levels', 'interface_levels']:
                    additional_dimensions[dimension] = state[dimension]

            # Set the units according to component's description, not global
            # description
            description[quantity]['units'] = quantities_to_extract[quantity]

            using_2d_coordinates = False
            x_coord = state['x']
            y_coord = state['y']
            mid_level_coord = state['mid_levels']
            interface_level_coord = state['interface_levels']

            if x_coord.ndim == 2:
                assert y_coord.ndim == 2
                using_2d_coordinates = True
                x_coord = state['logical_x_coordinate']
                y_coord = state['logical_y_coordinate']

            quantity_data_array = get_default_values(quantity,
                                                     x_coord, y_coord,
                                                     mid_level_coord,
                                                     interface_level_coord,
                                                     None,
                                                     description,
                                                     additional_dimensions)

            if using_2d_coordinates:
                physical_x = state['x']
                physical_y = state['y']

                quantity_data_array.coords[physical_x.label] = (
                    physical_x.dims, physical_x.values)

                quantity_data_array.coords[physical_y.label] = (
                    physical_y.dims, physical_y.values)

            output_state[quantity] = quantity_data_array

        return output_state

    def check_if_sane_and_return_attribute(self, attribute):
        """
        Check if attribute exists and is a dict

        Args:
            attribute (string):
                one of the attributes of the component.
        """

        if not hasattr(self, attribute):
            raise IndexError(
                'Component has no attribute called {}'.format(attribute))

        quantities_to_extract = getattr(self, attribute)

        if not isinstance(quantities_to_extract, dict):
            raise NotImplementedError(
                'This method will only work with components with a dict-like {} attribute'.format(attribute))

        return quantities_to_extract


def get_array_from_state(state,
                         quantity_name,
                         quantity_units,
                         quantity_dims):

    if quantity_name not in state:
        raise IndexError(
            'The input state does not contain {}'.format(quantity_name))

    return get_numpy_array(
        state[quantity_name].to_units(quantity_units), quantity_dims)


class ClimtPrognostic(ArrayHandler, Prognostic):
    """
    The base class to use for all CliMT Prognostics.

    It inherits from :py:class:`~sympl.Prognostic`, and adds a few methods that
    make array and state dictionary creation easier.

    Attributes:

        _climt_inputs (dict):
            The inputs expected by the component. The keys are the
            names of the quantities (preferably in CF convention),
            and the values are the units in which the component requires
            the quantity.

        _climt_tendencies (dict):
            The tendencies that are returned by the component. The units
            here are actually to tell CliMT about the units which are returned
            by the compiled extensions, if any. The developer **must** convert tendencies
            to (quantity)/second before returning them when the object is called.

        _climt_diagnostics (dict):
            The diagnostics that are returned by the component.

    """

    _climt_inputs = {}

    _climt_tendencies = {}

    _climt_diagnostics = {}

    @property
    def inputs(self):
        return tuple(self._climt_inputs.keys())

    @property
    def tendencies(self):
        return tuple(self._climt_tendencies.keys())

    @property
    def diagnostics(self):
        return tuple(self._climt_diagnostics.keys())

    @property
    def input_properties(self):
        return self.create_properties_dict(self._climt_inputs)

    @property
    def tendency_properties(self):
        return self.create_properties_dict(self._climt_tendencies)

    @property
    def diagnostic_properties(self):
        return self.create_properties_dict(self._climt_diagnostics)

    def piecewise_constant_version(self, update_time):
        """
        Returns component that updates once every :code:`update_time`.

        Args:
            update_time (timedelta):
                The time difference between updates.

        Returns:
            component (UpdateFrequencyWrapper):
                A "piecewise constant" component.

        """

        return UpdateFrequencyWrapper(self, update_time)

    def scaled_version(self,
                       input_scale_factors,
                       diagnostic_scale_factors,
                       tendency_scale_factors):
        """
        Returns component whose input/outputs/tendencies/diagnostics are scaled
        by the given scale factors.

        Args:
            input_scale_factors (dict):
                a dictionary whose keys are the inputs that will be scaled
                and values are floating point scaling factors.
            tendency_scale_factors (dict):
               a dictionary whose keys are the tendencies that will be scaled
               and values are floating point scaling factors.
            diagnostic_scale_factors (dict):
               a dictionary whose keys are the diagnostics that will be scaled
               and values are floating point scaling factors.

        """

        return ScalingWrapper(self,
                              input_scale_factors=input_scale_factors,
                              tendency_scale_factors=tendency_scale_factors,
                              diagnostic_scale_factors=diagnostic_scale_factors)


class ClimtDiagnostic(ArrayHandler, Diagnostic):
    """
    The base class to use for all CliMT Diagnostics.

    It inherits from :py:class:`~sympl.Diagnostic`, and adds a few methods that
    make array and state dictionary creation easier.

    Attributes:

        _climt_inputs (dict):
            The inputs expected by the component. The keys are the
            names of the quantities (preferably in CF convention),
            and the values are the units in which the component requires
            the quantity.

        _climt_diagnostics (dict):
            The diagnostics that are returned by the component.

    """

    _climt_inputs = {}

    _climt_tendencies = {}

    @property
    def inputs(self):
        return tuple(self._climt_inputs.keys())

    @property
    def diagnostics(self):
        return tuple(self._climt_diagnostics.keys())

    @property
    def input_properties(self):
        return self.create_properties_dict(self._climt_inputs)

    @property
    def diagnostic_properties(self):
        return self.create_properties_dict(self._climt_diagnostics)

    def scaled_version(self,
                       input_scale_factors,
                       diagnostic_scale_factors):
        """
        Returns component whose input/outputs/tendencies/diagnostics are scaled
        by the given scale factors.

        Args:
            input_scale_factors (dict):
                a dictionary whose keys are the inputs that will be scaled
                and values are floating point scaling factors.
            diagnostic_scale_factors (dict):
               a dictionary whose keys are the diagnostics that will be scaled
               and values are floating point scaling factors.

        """

        return ScalingWrapper(self,
                              input_scale_factors=input_scale_factors,
                              diagnostic_scale_factors=diagnostic_scale_factors)


class ClimtImplicit(ArrayHandler, Implicit):
    """
    The base class to use for all CliMT Implicits.

    It inherits from :py:class:`~sympl.Implicit`, and adds a few methods that
    make array and state dictionary creation easier.

    Attributes:

        _climt_inputs (dict):
            The inputs expected by the component. The keys are the
            names of the quantities (preferably in CF convention),
            and the values are the units in which the component requires
            the quantity.

        _climt_outputs (dict):
            The outputs returned by the component. These are new values
            of the state variables.

        _climt_diagnostics (dict):
            The diagnostics that are returned by the component.

    """

    _climt_inputs = {}

    _climt_outputs = {}

    _climt_diagnostics = {}

    @property
    def inputs(self):
        return tuple(self._climt_inputs.keys())

    @property
    def outputs(self):
        return tuple(self._climt_outputs.keys())

    @property
    def diagnostics(self):
        return tuple(self._climt_diagnostics.keys())

    @property
    def input_properties(self):
        return self.create_properties_dict(self._climt_inputs)

    @property
    def output_properties(self):
        return self.create_properties_dict(self._climt_outputs)

    @property
    def diagnostic_properties(self):
        return self.create_properties_dict(self._climt_diagnostics)

    def prognostic_version(self):
        """
        Returns a Prognostic component whose tendencies are the time differenced
        outputs of this Implicit component.

        """
        return ClimtTimeDifferenced(self)

    def scaled_version(self,
                       input_scale_factors,
                       diagnostic_scale_factors,
                       output_scale_factors):
        """
        Returns component whose input/outputs/tendencies/diagnostics are scaled
        by the given scale factors.

        Args:
            input_scale_factors (dict):
                a dictionary whose keys are the inputs that will be scaled
                and values are floating point scaling factors.
            output_scale_factors (dict):
               a dictionary whose keys are the tendencies that will be scaled
               and values are floating point scaling factors.
            diagnostic_scale_factors (dict):
               a dictionary whose keys are the diagnostics that will be scaled
               and values are floating point scaling factors.

        """

        return ScalingWrapper(self,
                              input_scale_factors=input_scale_factors,
                              output_scale_factors=output_scale_factors,
                              diagnostic_scale_factors=diagnostic_scale_factors)


class ClimtImplicitPrognostic(ClimtPrognostic):
    """
    The base class used mainly for convection schemes.

    Convection schemes return tendencies, but require the current model time step
    to work. This can be worked around by simply adding an attribute storing the
    time step in :code:`Prognostic`, but then we lose the distinction between these
    kind of schemes and pure Prognostics like radiative transfer codes. Depending on
    availability of time and the possibility of resolving this situation (by
    decomposing convection schemes into multiple components?), **this class might eventually
    be deprecated.**

    Attributes:

        current_time_step (timedelta): The current time step selected for the model. This
            should be updated on every iteration if the time step is not a constant.

    """

    current_time_step = timedelta(seconds=100.)


class ClimtSpectralDynamicalCore(ArrayHandler, TimeStepper):
    """
    The base class to use for spectral dynamical cores.

    Spectral Dynamical cores usually step the model in spectral space, and therefore
    the tendencies are required in spectral space as well. Therefore, spectral dycores
    cannot be simple :code:`Implicit` objects. :code:`SpectralDynamicalCore` objects should
    take in an initial state, convert it to spectral space, and provide future state in
    grid space. These objects need to perform the following steps:
        * take in a list of components which provide tendencies
        * create grid space quantities from the native spectral format
        * take in tendency terms in grid space
        * convert tendency terms to spectral space
        * step state forward

    Attributes:

        _climt_inputs (list): A list of quantities that the model steps forward
            by default -- like zonal_wind, etc.

        prognostics (list):
            A property that is set to a list of prognostics which
            constitute the "physics" for the dynamical core.

        _climt_inputs (dict):
            A list of quantities that the model steps forward
            by default -- like zonal_wind, etc.

        _climt_diagnostics (dict):
            The diagnostics that are returned by the dynamical core itself.

        inputs (tuple of str):
            The quantities required by the dynamical core when it is called.

        outputs (tuple of str):
            The quantities which are updated in the returned state.

        diagnostics (tuple of str):
            The quantities, calculated based on the input state (not new state),
            which can be used to diagnose behaviour of the model.

    """

    _climt_inputs = {}

    _climt_outputs = {}

    _climt_diagnostics = {}

    _prognostic = None

    @property
    def prognostics(self):
        return self._prognostic

    @prognostics.setter
    def prognostics(self, prognostic_list):
        self._prognostic = PrognosticComposite(*prognostic_list)

    @property
    def inputs(self):
        if self._prognostic is not None:
            return set(self._prognostic.inputs).union(set(self._climt_inputs.keys()))
        else:
            return set(self._climt_inputs.keys())

    @property
    def outputs(self):
        return set(self._climt_outputs.keys())

    @property
    def diagnostics(self):
        if self._prognostic is not None:
            return set(self._prognostic.diagnostics).union(set(self._climt_diagnostics.keys()))
        else:
            return set(self._climt_diagnostics.keys())

    @property
    def input_properties(self):
        return self.create_properties_dict(self._climt_inputs)

    @property
    def output_properties(self):
        return self.create_properties_dict(self._climt_outputs)

    @property
    def diagnostic_properties(self):
        return self.create_properties_dict(self._climt_diagnostics)

    def scaled_version(self,
                       input_scale_factors,
                       diagnostic_scale_factors,
                       output_scale_factors):
        """
        Returns component whose input/outputs/tendencies/diagnostics are scaled
        by the given scale factors.

        Args:
            input_scale_factors (dict):
                a dictionary whose keys are the inputs that will be scaled
                and values are floating point scaling factors.
            output_scale_factors (dict):
               a dictionary whose keys are the tendencies that will be scaled
               and values are floating point scaling factors.
            diagnostic_scale_factors (dict):
               a dictionary whose keys are the diagnostics that will be scaled
               and values are floating point scaling factors.

        """

        return ScalingWrapper(self,
                              input_scale_factors=input_scale_factors,
                              output_scale_factors=output_scale_factors,
                              diagnostic_scale_factors=diagnostic_scale_factors)


class ClimtTimeDifferenced(ClimtImplicitPrognostic):

    def __init__(self, implicit_component):
        """
        Create time differenced version of implicit_component.
        """

        self._implicit_component = implicit_component
        self._climt_inputs = copy.copy(implicit_component._climt_inputs)
        self._climt_diagnostics = copy.copy(implicit_component._climt_diagnostics)
        self._climt_tendencies = copy.copy(implicit_component._climt_outputs)

        for quantity in self._climt_tendencies:
            units = self._climt_tendencies[quantity]
            self._climt_tendencies[quantity] = units + ' s^-1'

    def __call__(self, state):
        """
        Create tendencies using first-order time differencing.
        """

        diagnostics, new_state = self._implicit_component(state,
                                                          self.current_time_step)

        tendencies = self.create_state_dict_for('_climt_tendencies', state)

        for quantity in tendencies.keys():
            tendencies[quantity].values = (
                new_state[quantity] - state[quantity])/self.current_time_step.total_seconds()

        return tendencies, diagnostics
