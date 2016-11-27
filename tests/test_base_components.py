import pytest
import unittest
import mock
import climt
from climt import (
    Prognostic, Diagnostic, Monitor, PrognosticCollection, DiagnosticCollection,
    MonitorCollection, SharedKeyException,
)


def test_empty_prognostic_collection():
    prognostic_collection = PrognosticCollection([])
    state = {'air_temperature': 273.15}
    prognostic_collection.ensure_state_is_valid_input(state)
    tendencies, diagnostics = prognostic_collection(state)
    assert len(tendencies) == 0
    assert len(diagnostics) == 0
    assert isinstance(tendencies, dict)
    assert isinstance(diagnostics, dict)


@mock.patch.object(climt.Prognostic, 'ensure_state_is_valid_input')
@mock.patch.object(climt.Prognostic, '__call__')
def test_prognostic_collection_calls_one_prognostic(mock_call, mock_ensure):
    mock_call.return_value = ({'air_temperature': 0.5}, {'foo': 50.})
    mock_ensure.return_value = None
    prognostic_collection = PrognosticCollection([Prognostic()])
    state = {'air_temperature': 273.15}
    prognostic_collection.ensure_state_is_valid_input(state)
    assert mock_ensure.called
    tendencies, diagnostics = prognostic_collection(state)
    assert mock_call.called
    assert tendencies == {'air_temperature': 0.5}
    assert diagnostics == {'foo': 50.}


@mock.patch.object(climt.Prognostic, 'ensure_state_is_valid_input')
@mock.patch.object(climt.Prognostic, '__call__')
def test_prognostic_collection_calls_two_prognostics(mock_call, mock_ensure):
    mock_call.return_value = ({'air_temperature': 0.5}, {})
    prognostic_collection = PrognosticCollection([Prognostic(), Prognostic()])
    state = {'air_temperature': 273.15}
    prognostic_collection.ensure_state_is_valid_input(state)
    assert mock_ensure.called
    assert mock_ensure.call_count == 2
    tendencies, diagnostics = prognostic_collection(state)
    assert mock_call.called
    assert mock_call.call_count == 2
    assert tendencies == {'air_temperature': 1.}
    assert diagnostics == {}


def test_empty_diagnostic_collection():
    diagnostic_collection = DiagnosticCollection([])
    state = {'air_temperature': 273.15}
    diagnostic_collection.ensure_state_is_valid_input(state)
    diagnostics = diagnostic_collection(state)
    assert len(diagnostics) == 0
    assert isinstance(diagnostics, dict)


@mock.patch.object(climt.Diagnostic, 'ensure_state_is_valid_input')
@mock.patch.object(climt.Diagnostic, '__call__')
def test_diagnostic_collection_calls_one_diagnostic(mock_call, mock_ensure):
    mock_call.return_value = ({'foo': 50.})
    mock_ensure.return_value = None
    diagnostic_collection = DiagnosticCollection([Diagnostic()])
    state = {'air_temperature': 273.15}
    diagnostic_collection.ensure_state_is_valid_input(state)
    assert mock_ensure.called
    diagnostics = diagnostic_collection(state)
    assert mock_call.called
    assert diagnostics == {'foo': 50.}


def test_empty_monitor_collection():
    # mainly we're testing that nothing errors
    monitor_collection = MonitorCollection([])
    state = {'air_temperature': 273.15}
    monitor_collection.ensure_state_is_valid_input(state)
    monitor_collection.store(state)


@mock.patch.object(climt.Monitor, 'ensure_state_is_valid_input')
@mock.patch.object(climt.Monitor, 'store')
def test_monitor_collection_calls_one_monitor(mock_store, mock_ensure):
    mock_store.return_value = None
    mock_ensure.return_value = None
    monitor_collection = MonitorCollection([Monitor()])
    state = {'air_temperature': 273.15}
    monitor_collection.ensure_state_is_valid_input(state)
    assert mock_ensure.called
    monitor_collection.store(state)
    assert mock_store.called


@mock.patch.object(climt.Monitor, 'ensure_state_is_valid_input')
@mock.patch.object(climt.Monitor, 'store')
def test_monitor_collection_calls_two_monitors(mock_store, mock_ensure):
    mock_store.return_value = None
    mock_ensure.return_value = None
    monitor_collection = MonitorCollection([Monitor(), Monitor()])
    state = {'air_temperature': 273.15}
    monitor_collection.ensure_state_is_valid_input(state)
    assert mock_ensure.called
    assert mock_ensure.call_count == 2
    monitor_collection.store(state)
    assert mock_store.called
    assert mock_store.call_count == 2


if __name__ == '__main__':
    pytest.main([__file__])
