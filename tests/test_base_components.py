import pytest
import mock
from climt import (
    Prognostic, Diagnostic, Monitor, PrognosticComposite, DiagnosticComposite,
    MonitorComposite, SharedKeyException,
)


def same_list(list1, list2):
    return (len(list1) == len(list2) and all(
        [item in list2 for item in list1] + [item in list1 for item in list2]))


class MockPrognostic(Prognostic):

    def __call__(self, state):
        return {}, {}


class MockDiagnostic(Diagnostic):

    def __call__(self, state):
        return {}


class MockMonitor(Monitor):

    def store(self, state):
        return


def test_empty_prognostic_composite():
    prognostic_composite = PrognosticComposite([])
    state = {'air_temperature': 273.15}
    tendencies, diagnostics = prognostic_composite(state)
    assert len(tendencies) == 0
    assert len(diagnostics) == 0
    assert isinstance(tendencies, dict)
    assert isinstance(diagnostics, dict)


@mock.patch.object(MockPrognostic, '__call__')
def test_prognostic_composite_calls_one_prognostic(mock_call):
    mock_call.return_value = ({'air_temperature': 0.5}, {'foo': 50.})
    prognostic_composite = PrognosticComposite([MockPrognostic()])
    state = {'air_temperature': 273.15}
    tendencies, diagnostics = prognostic_composite(state)
    assert mock_call.called
    assert tendencies == {'air_temperature': 0.5}
    assert diagnostics == {'foo': 50.}


@mock.patch.object(MockPrognostic, '__call__')
def test_prognostic_composite_calls_two_prognostics(mock_call):
    mock_call.return_value = ({'air_temperature': 0.5}, {})
    prognostic_composite = PrognosticComposite(
        [MockPrognostic(), MockPrognostic()])
    state = {'air_temperature': 273.15}
    tendencies, diagnostics = prognostic_composite(state)
    assert mock_call.called
    assert mock_call.call_count == 2
    assert tendencies == {'air_temperature': 1.}
    assert diagnostics == {}


@mock.patch.object(MockPrognostic, '__call__')
def test_prognostic_errors_when_overlapping_diagnostics(mock_call):
    """Test that when two Prognostic objects in a collection return the same
    diagnostic, a SharedKeyException is raised."""
    mock_call.return_value = ({'air_temperature': 0.5}, {'foo': 1.})
    prognostic_composite = PrognosticComposite(
        [MockPrognostic(), MockPrognostic()])
    state = {'air_temperature': 273.15}
    try:
        tendencies, diagnostics = prognostic_composite(state)
    except SharedKeyException:
        pass
    except Exception as err:
        raise err
    else:
        raise AssertionError


def test_empty_diagnostic_composite():
    diagnostic_composite = DiagnosticComposite([])
    state = {'air_temperature': 273.15}
    diagnostics = diagnostic_composite(state)
    assert len(diagnostics) == 0
    assert isinstance(diagnostics, dict)


@mock.patch.object(MockDiagnostic, '__call__')
def test_diagnostic_composite_calls_one_diagnostic(mock_call):
    mock_call.return_value = {'foo': 50.}
    diagnostic_composite = DiagnosticComposite([MockDiagnostic()])
    state = {'air_temperature': 273.15}
    diagnostics = diagnostic_composite(state)
    assert mock_call.called
    assert diagnostics == {'foo': 50.}


@mock.patch.object(MockPrognostic, '__call__')
def test_diagnostic_errors_when_overlapping_diagnostics(mock_call):
    """Test that when two Prognostic objects in a collection return the same
    diagnostic, a SharedKeyException is raised."""
    mock_call.return_value = ({'air_temperature': 0.5}, {'foo': 1.})
    diagnostic_composite = PrognosticComposite(
        [MockPrognostic(), MockPrognostic()])
    state = {'air_temperature': 273.15}
    try:
        diagnostics = diagnostic_composite(state)
    except SharedKeyException:
        pass
    except Exception as err:
        raise err
    else:
        raise AssertionError


def test_empty_monitor_collection():
    # mainly we're testing that nothing errors
    monitor_collection = MonitorComposite([])
    state = {'air_temperature': 273.15}
    monitor_collection.store(state)


@mock.patch.object(MockMonitor, 'store')
def test_monitor_collection_calls_one_monitor(mock_store):
    mock_store.return_value = None
    monitor_collection = MonitorComposite([MockMonitor()])
    state = {'air_temperature': 273.15}
    monitor_collection.store(state)
    assert mock_store.called


@mock.patch.object(MockMonitor, 'store')
def test_monitor_collection_calls_two_monitors(mock_store):
    mock_store.return_value = None
    monitor_collection = MonitorComposite([MockMonitor(), MockMonitor()])
    state = {'air_temperature': 273.15}
    monitor_collection.store(state)
    assert mock_store.called
    assert mock_store.call_count == 2


def test_prognostic_composite_cannot_use_diagnostic():
    try:
        PrognosticComposite([MockDiagnostic()])
    except TypeError:
        pass
    except Exception as err:
        raise err
    else:
        raise AssertionError('TypeError should have been raised')


def test_diagnostic_composite_cannot_use_prognostic():
    try:
        DiagnosticComposite([MockPrognostic()])
    except TypeError:
        pass
    except Exception as err:
        raise err
    else:
        raise AssertionError('TypeError should have been raised')


@mock.patch.object(MockDiagnostic, '__call__')
def test_diagnostic_composite_call(mock_call):
    mock_call.return_value = {'foo': 5.}
    state = {'bar': 10.}
    diagnostics = DiagnosticComposite([MockDiagnostic()])
    new_state = diagnostics(state)
    assert list(state.keys()) == ['bar']
    assert state['bar'] == 10.
    assert list(new_state.keys()) == ['foo']
    assert new_state['foo'] == 5.


@mock.patch.object(MockDiagnostic, '__call__')
def test_diagnostic_composite_update(mock_call):
    mock_call.return_value = {'foo': 5.}
    state = {'bar': 10.}
    diagnostics = DiagnosticComposite([MockDiagnostic()])
    diagnostics.update_state(state)
    assert len(state.keys()) == 2
    assert 'foo' in state.keys()
    assert 'bar' in state.keys()
    assert state['foo'] == 5.
    assert state['bar'] == 10.


@mock.patch.object(MockDiagnostic, '__call__')
def test_diagnostic_composite_update_will_not_overwrite(mock_call):
    mock_call.return_value = {'foo': 5.}
    state = {'foo': 10.}
    diagnostics = DiagnosticComposite([MockDiagnostic()])
    try:
        diagnostics.update_state(state)
    except SharedKeyException:
        pass
    except Exception as err:
        raise err
    else:
        raise AssertionError('SharedKeyException should have been raised')


def test_prognostic_composite_includes_attributes():
    prognostic = MockPrognostic()
    prognostic.inputs = ('input1',)
    prognostic.diagnostics = ('diagnostic1',)
    prognostic.tendencies = ('tendency1',)
    composite = PrognosticComposite([prognostic])
    assert composite.inputs == ('input1',)
    assert composite.diagnostics == ('diagnostic1',)
    assert composite.tendencies == ('tendency1',)


def test_prognostic_composite_includes_attributes_from_two():
    prognostic1 = MockPrognostic()
    prognostic1.inputs = ('input1',)
    prognostic1.diagnostics = ('diagnostic1',)
    prognostic1.tendencies = ('tendency1',)
    prognostic2 = MockPrognostic()
    prognostic2.inputs = ('input2',)
    prognostic2.diagnostics = ('diagnostic2',)
    prognostic2.tendencies = ('tendency2',)
    composite = PrognosticComposite([prognostic1, prognostic2])
    assert same_list(composite.inputs, ('input1', 'input2'))
    assert same_list(composite.diagnostics, ('diagnostic1', 'diagnostic2'))
    assert same_list(composite.tendencies, ('tendency1', 'tendency2'))


def test_prognostic_merges_attributes():
    prognostic1 = MockPrognostic()
    prognostic1.inputs = ('input1',)
    prognostic1.diagnostics = ('diagnostic1', 'diagnostic2')
    prognostic1.tendencies = ('tendency1', 'tendency2')
    prognostic2 = MockPrognostic()
    prognostic2.inputs = ('input1', 'input2')
    prognostic2.diagnostics = ('diagnostic1', 'diagnostic2')
    prognostic2.tendencies = ('tendency2',)
    composite = PrognosticComposite([prognostic1, prognostic2])
    assert same_list(composite.inputs, ('input1', 'input2'))
    assert same_list(composite.diagnostics, ('diagnostic1', 'diagnostic2'))
    assert same_list(composite.tendencies, ('tendency1', 'tendency2'))


def test_diagnostic_composite_includes_attributes():
    diagnostic = MockDiagnostic()
    diagnostic.inputs = ('input1',)
    diagnostic.diagnostics = ('diagnostic1',)
    composite = DiagnosticComposite([diagnostic])
    assert composite.inputs == ('input1',)
    assert composite.diagnostics == ('diagnostic1',)


def test_diagnostic_composite_includes_attributes_from_two():
    diagnostic1 = MockDiagnostic()
    diagnostic1.inputs = ('input1',)
    diagnostic1.diagnostics = ('diagnostic1',)
    diagnostic2 = MockDiagnostic()
    diagnostic2.inputs = ('input2',)
    diagnostic2.diagnostics = ('diagnostic2',)
    composite = DiagnosticComposite([diagnostic1, diagnostic2])
    assert same_list(composite.inputs, ('input1', 'input2'))
    assert same_list(composite.diagnostics, ('diagnostic1', 'diagnostic2'))


def test_diagnostic_merges_attributes():
    diagnostic1 = MockDiagnostic()
    diagnostic1.inputs = ('input1',)
    diagnostic1.diagnostics = ('diagnostic1', 'diagnostic2')
    diagnostic2 = MockDiagnostic()
    diagnostic2.inputs = ('input1', 'input2')
    diagnostic2.diagnostics = ('diagnostic1', 'diagnostic2')
    composite = DiagnosticComposite([diagnostic1, diagnostic2])
    assert same_list(composite.inputs, ('input1', 'input2'))
    assert same_list(composite.diagnostics, ('diagnostic1', 'diagnostic2'))

if __name__ == '__main__':
    pytest.main([__file__])
