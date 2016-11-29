import pytest
import mock
from climt import (
    Prognostic, Diagnostic, Monitor, PrognosticCollection, DiagnosticCollection,
    MonitorCollection, SharedKeyException,
)


class MockPrognostic(Prognostic):

    def ensure_state_is_valid_input(self, state):
        return

    def __call__(self, state):
        return {}, {}


class MockDiagnostic(Diagnostic):

    def ensure_state_is_valid_input(self, state):
        return

    def __call__(self, state):
        return {}


class MockMonitor(Monitor):
    def ensure_state_is_valid_input(self, state):
        return

    def store(self, state):
        return


def test_empty_prognostic_collection():
    prognostic_collection = PrognosticCollection([])
    state = {'air_temperature': 273.15}
    prognostic_collection.ensure_state_is_valid_input(state)
    tendencies, diagnostics = prognostic_collection(state)
    assert len(tendencies) == 0
    assert len(diagnostics) == 0
    assert isinstance(tendencies, dict)
    assert isinstance(diagnostics, dict)


@mock.patch.object(MockPrognostic, 'ensure_state_is_valid_input')
@mock.patch.object(MockPrognostic, '__call__')
def test_prognostic_collection_calls_one_prognostic(mock_call, mock_ensure):
    mock_call.return_value = ({'air_temperature': 0.5}, {'foo': 50.})
    mock_ensure.return_value = None
    prognostic_collection = PrognosticCollection([MockPrognostic()])
    state = {'air_temperature': 273.15}
    prognostic_collection.ensure_state_is_valid_input(state)
    assert mock_ensure.called
    tendencies, diagnostics = prognostic_collection(state)
    assert mock_call.called
    assert tendencies == {'air_temperature': 0.5}
    assert diagnostics == {'foo': 50.}


@mock.patch.object(MockPrognostic, 'ensure_state_is_valid_input')
@mock.patch.object(MockPrognostic, '__call__')
def test_prognostic_collection_calls_two_prognostics(mock_call, mock_ensure):
    mock_call.return_value = ({'air_temperature': 0.5}, {})
    prognostic_collection = PrognosticCollection(
        [MockPrognostic(), MockPrognostic()])
    state = {'air_temperature': 273.15}
    prognostic_collection.ensure_state_is_valid_input(state)
    assert mock_ensure.called
    assert mock_ensure.call_count == 2
    tendencies, diagnostics = prognostic_collection(state)
    assert mock_call.called
    assert mock_call.call_count == 2
    assert tendencies == {'air_temperature': 1.}
    assert diagnostics == {}


@mock.patch.object(MockPrognostic, 'ensure_state_is_valid_input')
@mock.patch.object(MockPrognostic, '__call__')
def test_prognostic_errors_when_overlapping_diagnostics(mock_call, mock_ensure):
    """Test that when two Prognostic objects in a collection return the same
    diagnostic, a SharedKeyException is raised."""
    mock_call.return_value = ({'air_temperature': 0.5}, {'foo': 1.})
    prognostic_collection = PrognosticCollection(
        [MockPrognostic(), MockPrognostic()])
    state = {'air_temperature': 273.15}
    try:
        tendencies, diagnostics = prognostic_collection(state)
    except SharedKeyException:
        pass
    except Exception as err:
        raise err
    else:
        raise AssertionError


def test_empty_diagnostic_collection():
    diagnostic_collection = DiagnosticCollection([])
    state = {'air_temperature': 273.15}
    diagnostic_collection.ensure_state_is_valid_input(state)
    diagnostics = diagnostic_collection(state)
    assert len(diagnostics) == 0
    assert isinstance(diagnostics, dict)


@mock.patch.object(MockDiagnostic, 'ensure_state_is_valid_input')
@mock.patch.object(MockDiagnostic, '__call__')
def test_diagnostic_collection_calls_one_diagnostic(mock_call, mock_ensure):
    mock_call.return_value = ({'foo': 50.})
    mock_ensure.return_value = None
    diagnostic_collection = DiagnosticCollection([MockDiagnostic()])
    state = {'air_temperature': 273.15}
    diagnostic_collection.ensure_state_is_valid_input(state)
    assert mock_ensure.called
    diagnostics = diagnostic_collection(state)
    assert mock_call.called
    assert diagnostics == {'foo': 50.}


@mock.patch.object(MockPrognostic, 'ensure_state_is_valid_input')
@mock.patch.object(MockPrognostic, '__call__')
def test_diagnostic_errors_when_overlapping_diagnostics(mock_call, mock_ensure):
    """Test that when two Prognostic objects in a collection return the same
    diagnostic, a SharedKeyException is raised."""
    mock_call.return_value = ({'air_temperature': 0.5}, {'foo': 1.})
    diagnostic_collection = PrognosticCollection(
        [MockPrognostic(), MockPrognostic()])
    state = {'air_temperature': 273.15}
    try:
        diagnostics = diagnostic_collection(state)
    except SharedKeyException:
        pass
    except Exception as err:
        raise err
    else:
        raise AssertionError


def test_empty_monitor_collection():
    # mainly we're testing that nothing errors
    monitor_collection = MonitorCollection([])
    state = {'air_temperature': 273.15}
    monitor_collection.ensure_state_is_valid_input(state)
    monitor_collection.store(state)


@mock.patch.object(MockMonitor, 'ensure_state_is_valid_input')
@mock.patch.object(MockMonitor, 'store')
def test_monitor_collection_calls_one_monitor(mock_store, mock_ensure):
    mock_store.return_value = None
    mock_ensure.return_value = None
    monitor_collection = MonitorCollection([MockMonitor()])
    state = {'air_temperature': 273.15}
    monitor_collection.ensure_state_is_valid_input(state)
    assert mock_ensure.called
    monitor_collection.store(state)
    assert mock_store.called


@mock.patch.object(MockMonitor, 'ensure_state_is_valid_input')
@mock.patch.object(MockMonitor, 'store')
def test_monitor_collection_calls_two_monitors(mock_store, mock_ensure):
    mock_store.return_value = None
    mock_ensure.return_value = None
    monitor_collection = MonitorCollection([MockMonitor(), MockMonitor()])
    state = {'air_temperature': 273.15}
    monitor_collection.ensure_state_is_valid_input(state)
    assert mock_ensure.called
    assert mock_ensure.call_count == 2
    monitor_collection.store(state)
    assert mock_store.called
    assert mock_store.call_count == 2


def test_prognostic_collection_cannot_use_diagnostic():
    try:
        PrognosticCollection([MockDiagnostic()])
    except TypeError:
        pass
    except Exception as err:
        raise err
    else:
        raise AssertionError('TypeError should have been raised')


def test_diagnostic_collection_cannot_use_prognostic():
    try:
        DiagnosticCollection([MockPrognostic()])
    except TypeError:
        pass
    except Exception as err:
        raise err
    else:
        raise AssertionError('TypeError should have been raised')


@mock.patch.object(MockDiagnostic, '__call__')
def test_diagnostic_collection_call(mock_call):
    mock_call.return_value = {'foo': 5.}
    state = {'bar': 10.}
    diagnostics = DiagnosticCollection([MockDiagnostic()])
    new_state = diagnostics(state)
    assert list(state.keys()) == ['bar']
    assert state['bar'] == 10.
    assert list(new_state.keys()) == ['foo']
    assert new_state['foo'] == 5.


@mock.patch.object(MockDiagnostic, '__call__')
def test_diagnostic_collection_update(mock_call):
    mock_call.return_value = {'foo': 5.}
    state = {'bar': 10.}
    diagnostics = DiagnosticCollection([MockDiagnostic()])
    diagnostics.update_state(state)
    assert len(state.keys()) == 2
    assert 'foo' in state.keys()
    assert 'bar' in state.keys()
    assert state['foo'] == 5.
    assert state['bar'] == 10.


@mock.patch.object(MockDiagnostic, '__call__')
def test_diagnostic_collection_update_will_not_overwrite(mock_call):
    mock_call.return_value = {'foo': 5.}
    state = {'foo': 10.}
    diagnostics = DiagnosticCollection([MockDiagnostic()])
    try:
        diagnostics.update_state(state)
    except SharedKeyException:
        pass
    except Exception as err:
        raise err
    else:
        raise AssertionError('SharedKeyException should have been raised')


if __name__ == '__main__':
    pytest.main([__file__])
