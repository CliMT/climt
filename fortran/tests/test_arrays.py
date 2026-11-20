import cffi
import numpy as np
import pytest

from climt_fortran._arrays import (
    as_double_pointer, as_double_scalar_pointer, as_int_pointer, check_array, copy_back_outputs, input_pointer, output_pointer)

ffi = cffi.FFI()


def test_double_pointer_reads_the_array():
    array = np.arange(6.0).reshape(2, 3)
    pointer = as_double_pointer(ffi, array)
    assert [pointer[i] for i in range(6)] == list(range(6))
    pointer[4] = -1.0
    assert array[1, 1] == -1.0


@pytest.mark.parametrize("bad", [np.zeros(3, np.float32), np.zeros(3, ">f8"), [0.0, 1.0]])
def test_double_pointer_rejects_other_types(bad):
    with pytest.raises(TypeError, match="float64"):
        as_double_pointer(ffi, bad)


@pytest.mark.parametrize("bad", [np.zeros((3, 2)).T, np.zeros(6)[::2]])
def test_double_pointer_rejects_non_contiguous(bad):
    with pytest.raises(ValueError, match="C-contiguous"):
        as_double_pointer(ffi, bad)


def test_double_pointer_rejects_read_only():
    array = np.zeros(3)
    array.flags.writeable = False
    with pytest.raises(ValueError, match="read-only"):
        as_double_pointer(ffi, array)


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


def test_input_pointer_copies_only_non_contiguous_arrays():
    array = np.arange(6.0).reshape(3, 2)
    assert input_pointer(ffi, "x", array, (3, 2))[3] == 3.0
    pointer = input_pointer(ffi, "x", array, (3, 2))
    pointer[0] = -1.0
    assert array[0, 0] == -1.0  # contiguous: the array itself
    transposed = np.arange(6.0).reshape(2, 3).T
    pointer = input_pointer(ffi, "x", transposed, (3, 2))
    assert [pointer[i] for i in range(6)] == list(transposed.ravel())
    with pytest.raises(TypeError, match="x: expected a float64"):
        input_pointer(ffi, "x", array.astype(np.float32), (3, 2))
    with pytest.raises(ValueError, match="x: expected shape"):
        input_pointer(ffi, "x", array, (2, 3))


def test_output_pointer_writes_back_non_contiguous_arrays():
    contiguous, strided = np.zeros((3, 2)), np.zeros((3, 4))[:, ::2]
    copy_back = []
    pointers = [output_pointer(ffi, name, array, (3, 2), copy_back)
                for name, array in (("a", contiguous), ("b", strided))]
    assert len(copy_back) == 1 and copy_back[0][0] is strided
    for pointer in pointers:
        for i in range(6):
            pointer[i] = i
    assert not strided.any()
    copy_back_outputs(copy_back)
    expected = np.arange(6.0).reshape(3, 2)
    np.testing.assert_array_equal(contiguous, expected)
    np.testing.assert_array_equal(strided, expected)


def test_output_pointer_rejects_read_only():
    for array in (np.zeros((3, 2)), np.zeros((3, 4))[:, ::2]):
        array.flags.writeable = False
        with pytest.raises(ValueError, match="read-only"):
            output_pointer(ffi, "x", array, (3, 2), [])
