import pytest
import numpy as np
from climt import DataArray


def test_array_addition():
    a = DataArray(np.array([1., 2., 3.]))
    b = DataArray(np.array([2., 1., 3.]))
    result = a + b
    assert (result.values == np.array([3., 3., 6.])).all()
    assert len(result.attrs) == 0


def test_array_subtraction():
    a = DataArray(np.array([1., 2., 3.]))
    b = DataArray(np.array([2., 1., 3.]))
    result = a - b
    assert (result.values == np.array([-1., 1., 0.])).all()
    assert len(result.attrs) == 0


def test_array_addition_keeps_left_attr():
    a = DataArray(np.array([1., 2., 3.]), attrs={'units': 'K'})
    b = DataArray(np.array([2., 1., 3.]), attrs={'units': 'm/s', 'foo': 'bar'})
    result = a + b
    assert (result.values == np.array([3., 3., 6.])).all()
    assert len(result.attrs) == 1
    assert result.attrs['units'] == 'K'


def test_array_subtraction_keeps_left_attrs():
    a = DataArray(np.array([1., 2., 3.]), attrs={'units': 'm/s', 'foo': 'bar'})
    b = DataArray(np.array([2., 1., 3.]), attrs={'units': 'K'})
    result = a - b
    assert (result.values == np.array([-1., 1., 0.])).all()
    assert len(result.attrs) == 2
    assert result.attrs['units'] == 'm/s'
    assert result.attrs['foo'] == 'bar'

if __name__ == '__main__':
    pytest.main([__file__])
