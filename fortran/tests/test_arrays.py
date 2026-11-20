import cffi
import numpy as np
import pytest

from climt_fortran._arrays import (
    as_double_scalar_pointer, as_int_pointer, check_array, contiguous_arguments, copy_back_outputs)

ffi = cffi.FFI()


def test_scalar_pointers():
    assert ffi.typeof(as_int_pointer(ffi, 7)) is ffi.typeof("int32_t *")
    assert as_int_pointer(ffi, 7)[0] == 7
    assert as_double_scalar_pointer(ffi, 2.5)[0] == 2.5
    with pytest.raises(TypeError):
        as_int_pointer(ffi, 1.5)


def test_check_array_accepts_any_layout():
    float64 = np.dtype(np.float64)
    for array in (np.zeros((3, 2)), np.zeros((2, 3)).T, np.zeros((3, 4))[:, ::2]):
        check_array("x", array, float64, (3, 2), output=True)
    check_array("x", np.zeros(5), float64)  # no shape given: any shape


@pytest.mark.parametrize("bad, error, message", [
    (np.zeros((3, 2), np.float32), TypeError, "x: expected a float64 array, got float32"),
    (np.zeros((3, 2), ">f8"), TypeError, "x: expected a float64 array, got >f8"),
    ([0.0, 1.0], TypeError, "x: expected a float64 array, got list"),
    (None, TypeError, "x: expected a float64 array, got NoneType"),
    (np.zeros((2, 3)), ValueError, r"x: expected shape \(3, 2\), got \(2, 3\)"),
])
def test_check_array_rejects(bad, error, message):
    with pytest.raises(error, match=message):
        check_array("x", bad, np.dtype(np.float64), (3, 2))


def test_check_array_output_must_be_writeable():
    array = np.zeros((3, 2))
    array.flags.writeable = False
    check_array("x", array, np.dtype(np.float64), (3, 2))  # fine as an input
    with pytest.raises(ValueError, match="x: array is read-only"):
        check_array("x", array, np.dtype(np.float64), (3, 2), output=True)


def _read_only(array):
    array.flags.writeable = False
    return array


def test_contiguous_arguments_copies_only_when_needed():
    values = np.arange(6.0).reshape(3, 2)
    contiguous = values.copy()
    transposed = np.asfortranarray(values)
    strided = np.zeros((3, 4))[:, ::2]
    strided[...] = values
    read_only = _read_only(values.copy())
    copy_back = []
    arrays = [contiguous, transposed, read_only, strided]
    specs = [("a", (3, 2), False), ("b", (3, 2), False), ("c", (3, 2), False), ("d", (3, 2), True)]
    result = contiguous_arguments(arrays, specs, copy_back)
    assert result[0] is contiguous
    for array in result:
        assert array.flags.c_contiguous and array.flags.writeable
        np.testing.assert_array_equal(array, values)
    assert len(copy_back) == 1 and copy_back[0][0] is strided  # outputs only
    result[3][...] = -1.0
    copy_back_outputs(copy_back)
    assert np.all(strided == -1.0)


def test_contiguous_arguments_checks_every_array():
    good = ("a", (3, 2), False)
    with pytest.raises(TypeError, match="b: expected a float64"):
        contiguous_arguments([np.zeros((3, 2)), np.zeros((3, 2), np.float32)], [good, ("b", (3, 2), False)], [])
    with pytest.raises(ValueError, match="b: array is read-only"):
        contiguous_arguments([np.zeros((3, 2)), _read_only(np.zeros((3, 2)))], [good, ("b", (3, 2), True)], [])
