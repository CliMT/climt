import cffi
import numpy as np
import pytest

from climt_fortran._arrays import as_double_pointer, as_double_scalar_pointer, as_int_pointer

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
