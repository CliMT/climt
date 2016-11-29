import pytest
import mock
from climt import Prognostic, Leapfrog
from datetime import timedelta
import numpy as np


class MockPrognostic(Prognostic):

    def ensure_state_is_valid_input(self, state):
        return

    def __call__(self):
        return {}, {}


@mock.patch.object(MockPrognostic, '__call__')
def test_leapfrog_float_no_change_one_step(mock_prognostic_call):
    mock_prognostic_call.return_value = ({'air_temperature': 0.}, {})
    state = {'air_temperature': 273.}
    timestep = timedelta(seconds=1.)
    time_stepper = Leapfrog([MockPrognostic()])
    new_state = time_stepper.step(state, timestep)
    assert state == {'air_temperature': 273.}
    assert new_state == {'air_temperature': 273.}


@mock.patch.object(MockPrognostic, '__call__')
def test_leapfrog_float_no_change_three_steps(mock_prognostic_call):
    mock_prognostic_call.return_value = ({'air_temperature': 0.}, {})
    state = {'air_temperature': 273.}
    timestep = timedelta(seconds=1.)
    time_stepper = Leapfrog([MockPrognostic()])
    new_state = time_stepper.step(state, timestep)
    assert state == {'air_temperature': 273.}
    assert new_state == {'air_temperature': 273.}
    state = new_state
    new_state = time_stepper.step(state, timestep)
    assert state == {'air_temperature': 273.}
    assert new_state == {'air_temperature': 273.}
    state = new_state
    new_state = time_stepper.step(state, timestep)
    assert state == {'air_temperature': 273.}
    assert new_state == {'air_temperature': 273.}


@mock.patch.object(MockPrognostic, '__call__')
def test_leapfrog_float_one_step(mock_prognostic_call):
    mock_prognostic_call.return_value = ({'air_temperature': 1.}, {})
    state = {'air_temperature': 273.}
    timestep = timedelta(seconds=1.)
    time_stepper = Leapfrog([MockPrognostic()])
    new_state = time_stepper.step(state, timestep)
    assert state == {'air_temperature': 273.}
    assert new_state == {'air_temperature': 274.}


@mock.patch.object(MockPrognostic, '__call__')
def test_leapfrog_float_three_steps(mock_prognostic_call):
    mock_prognostic_call.return_value = ({'air_temperature': 1.}, {})
    state = {'air_temperature': 273.}
    timestep = timedelta(seconds=1.)
    time_stepper = Leapfrog([MockPrognostic()])
    new_state = time_stepper.step(state, timestep)
    assert state == {'air_temperature': 273.}
    assert new_state == {'air_temperature': 274.}
    state = new_state
    new_state = time_stepper.step(new_state, timestep)
    assert state == {'air_temperature': 274.}
    assert new_state == {'air_temperature': 275.}
    state = new_state
    new_state = time_stepper.step(state, timestep)
    assert state == {'air_temperature': 275.}
    assert new_state == {'air_temperature': 276.}


@mock.patch.object(MockPrognostic, '__call__')
def test_leapfrog_float_two_steps_filtered(mock_prognostic_call):
    """Test that the Asselin filter is being correctly applied"""
    mock_prognostic_call.return_value = ({'air_temperature': 0.}, {})
    state = {'air_temperature': 273.}
    timestep = timedelta(seconds=1.)
    time_stepper = Leapfrog([MockPrognostic()], asselin_strength=0.5)
    new_state = time_stepper.step(state, timestep)
    assert state == {'air_temperature': 273.}
    assert new_state == {'air_temperature': 273.}
    state = new_state
    mock_prognostic_call.return_value = ({'air_temperature': 2.}, {})
    new_state = time_stepper.step(state, timestep)
    # Asselin filter modifies the current state
    assert state == {'air_temperature': 274.}
    assert new_state == {'air_temperature': 277.}


@mock.patch.object(MockPrognostic, '__call__')
def test_leapfrog_requires_same_timestep(mock_prognostic_call):
    """Test that the Asselin filter is being correctly applied"""
    mock_prognostic_call.return_value = ({'air_temperature': 0.}, {})
    state = {'air_temperature': 273.}
    timestep = timedelta(seconds=1.)
    time_stepper = Leapfrog([MockPrognostic()], asselin_strength=0.5)
    state = time_stepper.step(state, timedelta(seconds=1.))
    try:
        time_stepper.step(state, timedelta(seconds=2.))
    except ValueError:
        pass
    except Exception as err:
        raise err
    else:
        raise AssertionError('Leapfrog must require timestep to be constant')


@mock.patch.object(MockPrognostic, '__call__')
def test_leapfrog_array_no_change_one_step(mock_prognostic_call):
    mock_prognostic_call.return_value = (
        {'air_temperature': np.zeros((3, 3))}, {})
    state = {'air_temperature': np.ones((3, 3))*273.}
    timestep = timedelta(seconds=1.)
    time_stepper = Leapfrog([MockPrognostic()])
    new_state = time_stepper.step(state, timestep)
    assert list(state.keys()) == ['air_temperature']
    assert (state['air_temperature'] == np.ones((3, 3))*273.).all()
    assert list(new_state.keys()) == ['air_temperature']
    assert (new_state['air_temperature'] == np.ones((3, 3))*273.).all()


@mock.patch.object(MockPrognostic, '__call__')
def test_leapfrog_array_no_change_three_steps(mock_prognostic_call):
    mock_prognostic_call.return_value = (
        {'air_temperature': np.ones((3, 3))*0.}, {})
    state = {'air_temperature': np.ones((3, 3))*273.}
    timestep = timedelta(seconds=1.)
    time_stepper = Leapfrog([MockPrognostic()])
    new_state = time_stepper.step(state, timestep)
    assert list(state.keys()) == ['air_temperature']
    assert (state['air_temperature'] == np.ones((3, 3))*273.).all()
    assert list(new_state.keys()) == ['air_temperature']
    assert (new_state['air_temperature'] == np.ones((3, 3))*273.).all()
    state = new_state
    new_state = time_stepper.step(state, timestep)
    assert list(state.keys()) == ['air_temperature']
    assert (state['air_temperature'] == np.ones((3, 3))*273.).all()
    assert list(new_state.keys()) == ['air_temperature']
    assert (new_state['air_temperature'] == np.ones((3, 3))*273.).all()
    state = new_state
    new_state = time_stepper.step(state, timestep)
    assert list(state.keys()) == ['air_temperature']
    assert (state['air_temperature'] == np.ones((3, 3))*273.).all()
    assert list(new_state.keys()) == ['air_temperature']
    assert (new_state['air_temperature'] == np.ones((3, 3))*273.).all()


@mock.patch.object(MockPrognostic, '__call__')
def test_leapfrog_array_one_step(mock_prognostic_call):
    mock_prognostic_call.return_value = (
        {'air_temperature': np.ones((3, 3))*1.}, {})
    state = {'air_temperature': np.ones((3, 3))*273.}
    timestep = timedelta(seconds=1.)
    time_stepper = Leapfrog([MockPrognostic()])
    new_state = time_stepper.step(state, timestep)
    assert list(state.keys()) == ['air_temperature']
    assert (state['air_temperature'] == np.ones((3, 3))*273.).all()
    assert list(new_state.keys()) == ['air_temperature']
    assert (new_state['air_temperature'] == np.ones((3, 3))*274.).all()


@mock.patch.object(MockPrognostic, '__call__')
def test_leapfrog_array_three_steps(mock_prognostic_call):
    mock_prognostic_call.return_value = (
        {'air_temperature': np.ones((3, 3))*1.}, {})
    state = {'air_temperature': np.ones((3, 3))*273.}
    timestep = timedelta(seconds=1.)
    time_stepper = Leapfrog([MockPrognostic()])
    new_state = time_stepper.step(state, timestep)
    assert list(state.keys()) == ['air_temperature']
    assert (state['air_temperature'] == np.ones((3, 3))*273.).all()
    assert list(new_state.keys()) == ['air_temperature']
    assert (new_state['air_temperature'] == np.ones((3, 3))*274.).all()
    state = new_state
    new_state = time_stepper.step(new_state, timestep)
    assert list(state.keys()) == ['air_temperature']
    assert (state['air_temperature'] == np.ones((3, 3))*274.).all()
    assert list(new_state.keys()) == ['air_temperature']
    assert (new_state['air_temperature'] == np.ones((3, 3))*275.).all()
    state = new_state
    new_state = time_stepper.step(state, timestep)
    assert list(state.keys()) == ['air_temperature']
    assert (state['air_temperature'] == np.ones((3, 3))*275.).all()
    assert list(new_state.keys()) == ['air_temperature']
    assert (new_state['air_temperature'] == np.ones((3, 3))*276.).all()


@mock.patch.object(MockPrognostic, '__call__')
def test_leapfrog_array_two_steps_filtered(mock_prognostic_call):
    """Test that the Asselin filter is being correctly applied"""
    mock_prognostic_call.return_value = (
        {'air_temperature': np.ones((3, 3))*0.}, {})
    state = {'air_temperature': np.ones((3, 3))*273.}
    timestep = timedelta(seconds=1.)
    time_stepper = Leapfrog([MockPrognostic()], asselin_strength=0.5)
    new_state = time_stepper.step(state, timestep)
    assert list(state.keys()) == ['air_temperature']
    assert (state['air_temperature'] == np.ones((3, 3))*273.).all()
    assert list(new_state.keys()) == ['air_temperature']
    assert (new_state['air_temperature'] == np.ones((3, 3))*273.).all()
    state = new_state
    mock_prognostic_call.return_value = (
        {'air_temperature': np.ones((3, 3))*2.}, {})
    new_state = time_stepper.step(state, timestep)
    # Asselin filter modifies the current state
    assert list(state.keys()) == ['air_temperature']
    assert (state['air_temperature'] == np.ones((3, 3))*274.).all()
    assert list(new_state.keys()) == ['air_temperature']
    assert (new_state['air_temperature'] == np.ones((3, 3))*277.).all()


@mock.patch.object(MockPrognostic, '__call__')
def test_leapfrog_array_two_steps_filtered_williams(mock_prognostic_call):
    """Test that the Asselin filter is being correctly applied with a
    Williams factor of alpha=0.5"""
    mock_prognostic_call.return_value = (
        {'air_temperature': np.ones((3, 3))*0.}, {})
    state = {'air_temperature': np.ones((3, 3))*273.}
    timestep = timedelta(seconds=1.)
    time_stepper = Leapfrog([MockPrognostic()], asselin_strength=0.5, alpha=0.5)
    new_state = time_stepper.step(state, timestep)
    assert list(state.keys()) == ['air_temperature']
    assert (state['air_temperature'] == np.ones((3, 3))*273.).all()
    assert list(new_state.keys()) == ['air_temperature']
    assert (new_state['air_temperature'] == np.ones((3, 3))*273.).all()
    state = new_state
    mock_prognostic_call.return_value = (
        {'air_temperature': np.ones((3, 3))*2.}, {})
    new_state = time_stepper.step(state, timestep)
    # Asselin filter modifies the current state
    assert list(state.keys()) == ['air_temperature']
    assert (state['air_temperature'] == np.ones((3, 3))*273.5).all()
    assert list(new_state.keys()) == ['air_temperature']
    assert (new_state['air_temperature'] == np.ones((3, 3))*276.5).all()

if __name__ == '__main__':
    pytest.main([__file__])
