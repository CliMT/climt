import pytest
from climt import set_prognostic_update_frequency, Prognostic
from datetime import datetime, timedelta


class MockPrognostic(Prognostic):

    def __init__(self):
        self._num_updates = 0

    def __call__(self, state):
        self._num_updates += 1
        return {}, {'num_updates': self._num_updates}


def test_set_prognostic_update_frequency_calls_initially():
    set_prognostic_update_frequency(MockPrognostic, timedelta(hours=1))
    prognostic = MockPrognostic()
    state = {'time': timedelta(hours=0)}
    tendencies, diagnostics = prognostic(state)
    assert len(tendencies) == 0
    assert len(diagnostics) == 1
    assert diagnostics['num_updates'] == 1


def test_set_prognostic_update_frequency_caches_result():
    set_prognostic_update_frequency(MockPrognostic, timedelta(hours=1))
    prognostic = MockPrognostic()
    state = {'time': timedelta(hours=0)}
    tendencies, diagnostics = prognostic(state)
    tendencies, diagnostics = prognostic(state)
    assert len(tendencies) == 0
    assert len(diagnostics) == 1
    assert diagnostics['num_updates'] == 1


def test_set_prognostic_update_frequency_caches_result_with_datetime():
    set_prognostic_update_frequency(MockPrognostic, timedelta(hours=1))
    prognostic = MockPrognostic()
    state = {'time': datetime(2000, 1, 1)}
    tendencies, diagnostics = prognostic(state)
    tendencies, diagnostics = prognostic(state)
    assert len(tendencies) == 0
    assert len(diagnostics) == 1
    assert diagnostics['num_updates'] == 1


def test_set_prognostic_update_frequency_updates_result_when_equal():
    set_prognostic_update_frequency(MockPrognostic, timedelta(hours=1))
    prognostic = MockPrognostic()
    state = {'time': timedelta(hours=0)}
    tendencies, diagnostics = prognostic({'time': timedelta(hours=0)})
    tendencies, diagnostics = prognostic({'time': timedelta(hours=1)})
    assert len(tendencies) == 0
    assert len(diagnostics) == 1
    assert diagnostics['num_updates'] == 2


def test_set_prognostic_update_frequency_updates_result_when_greater():
    set_prognostic_update_frequency(MockPrognostic, timedelta(hours=1))
    prognostic = MockPrognostic()
    state = {'time': timedelta(hours=0)}
    tendencies, diagnostics = prognostic({'time': timedelta(hours=0)})
    tendencies, diagnostics = prognostic({'time': timedelta(hours=2)})
    assert len(tendencies) == 0
    assert len(diagnostics) == 1
    assert diagnostics['num_updates'] == 2


if __name__ == '__main__':
    pytest.main([__file__])
