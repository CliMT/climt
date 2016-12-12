import pytest
from climt import (
    ConstantPrognostic, ConstantDiagnostic, RelaxationPrognostic, DataArray
)
import numpy as np


def test_constant_prognostic_empty_dicts():
    prog = ConstantPrognostic({}, {})
    tendencies, diagnostics = prog({})
    assert isinstance(tendencies, dict)
    assert isinstance(diagnostics, dict)
    assert len(tendencies) == 0
    assert len(diagnostics) == 0


def test_constant_prognostic_cannot_modify_through_input_dict():
    in_tendencies = {}
    in_diagnostics = {}
    prog = ConstantPrognostic(in_tendencies, in_diagnostics)
    in_tendencies['a'] = 'b'
    in_diagnostics['c'] = 'd'
    tendencies, diagnostics = prog({})
    assert len(tendencies) == 0
    assert len(diagnostics) == 0


def test_constant_prognostic_cannot_modify_through_output_dict():
    prog = ConstantPrognostic({}, {})
    tendencies, diagnostics = prog({})
    tendencies['a'] = 'b'
    diagnostics['c'] = 'd'
    tendencies, diagnostics = prog({})
    assert len(tendencies) == 0
    assert len(diagnostics) == 0


def test_constant_diagnostic_empty_dict():
    diag = ConstantDiagnostic({})
    diagnostics = diag({})
    assert isinstance(diagnostics, dict)
    assert len(diagnostics) == 0


def test_constant_diagnostic_cannot_modify_through_input_dict():
    in_diagnostics = {}
    diag = ConstantDiagnostic(in_diagnostics)
    in_diagnostics['a'] = 'b'
    diagnostics = diag({})
    assert isinstance(diagnostics, dict)
    assert len(diagnostics) == 0


def test_constant_diagnostic_cannot_modify_through_output_dict():
    diag = ConstantDiagnostic({})
    diagnostics = diag({})
    diagnostics['c'] = 'd'
    diagnostics = diag({})
    assert len(diagnostics) == 0


def test_relaxation_prognostic_at_equilibrium():
    prognostic = RelaxationPrognostic('quantity')
    state = {
        'quantity': DataArray(np.array([0., 1., 2.]), attrs={'units': 'degK'}),
        'quantity_relaxation_timescale': DataArray(
            np.array([1., 1., 1.]), attrs={'units': 's'}),
        'equilibrium_quantity': DataArray(
            np.array([0., 1., 2.]), attrs={'units': 'degK'}),
    }
    tendencies, diagnostics = prognostic(state)
    assert np.all(tendencies['quantity'].values == 0.)


def test_relaxation_prognostic_with_change():
    prognostic = RelaxationPrognostic('quantity')
    state = {
        'quantity': DataArray(np.array([0., 1., 2.]), attrs={'units': 'degK'}),
        'quantity_relaxation_timescale': DataArray(
            np.array([1., 1., 1.]), attrs={'units': 's'}),
        'equilibrium_quantity': DataArray(
            np.array([1., 3., 5.]), attrs={'units': 'degK'}),
    }
    tendencies, diagnostics = prognostic(state)
    assert np.all(tendencies['quantity'].values == np.array([1., 2., 3.]))


def test_relaxation_prognostic_with_change_different_timescale_units():
    prognostic = RelaxationPrognostic('quantity')
    state = {
        'quantity': DataArray(np.array([0., 1., 2.]), attrs={'units': 'degK'}),
        'quantity_relaxation_timescale': DataArray(
            np.array([1/60., 2/60., 3/60.]), attrs={'units': 'minutes'}),
        'equilibrium_quantity': DataArray(
            np.array([1., 3., 5.]), attrs={'units': 'degK'}),
    }
    tendencies, diagnostics = prognostic(state)
    assert np.all(tendencies['quantity'].values == np.array([1., 1., 1.]))


def test_relaxation_prognostic_with_change_different_equilibrium_units():
    prognostic = RelaxationPrognostic('quantity')
    state = {
        'quantity': DataArray(np.array([0., 1., 2.]), attrs={'units': 'm'}),
        'quantity_relaxation_timescale': DataArray(
            np.array([1., 2., 3.]), attrs={'units': 's'}),
        'equilibrium_quantity': DataArray(
            np.array([1., 3., 5.])*1e-3, attrs={'units': 'km'}),
    }
    tendencies, diagnostics = prognostic(state)
    assert np.all(tendencies['quantity'].values == np.array([1., 1., 1.]))


def test_relaxation_prognostic_caching_timescale_and_equilibrium():
    prognostic = RelaxationPrognostic(
        'quantity',
        relaxation_timescale=DataArray(
            np.array([1., 1., 1.]), attrs={'units': 's'}),
        equilibrium_value=DataArray(
            np.array([1., 3., 5.]), attrs={'units': 'degK'}))
    state = {
        'quantity': DataArray(np.array([0., 1., 2.]), attrs={'units': 'degK'}),
    }
    tendencies, diagnostics = prognostic(state)
    assert np.all(tendencies['quantity'].values == np.array([1., 2., 3.]))


def test_relaxation_prognostic_cannot_override_cached_values():
    prognostic = RelaxationPrognostic(
        'quantity',
        relaxation_timescale=DataArray(
            np.array([1., 1., 1.]), attrs={'units': 's'}),
        equilibrium_value=DataArray(
            np.array([1., 3., 5.]), attrs={'units': 'degK'}))
    state = {
        'quantity': DataArray(np.array([0., 1., 2.]), attrs={'units': 'degK'}),
        # random values which should not be used by the prognostic
        'quantity_relaxation_timescale': DataArray(
            np.array([5., 1., 18.]), attrs={'units': 's'}),
        'equilibrium_quantity': DataArray(
            np.array([4., 7., 2.])*1e-3, attrs={'units': 'km'}),
    }
    tendencies, diagnostics = prognostic(state)
    assert np.all(tendencies['quantity'].values == np.array([1., 2., 3.]))

if __name__ == '__main__':
    pytest.main([__file__])
