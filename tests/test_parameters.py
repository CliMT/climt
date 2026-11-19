from collections import namedtuple
from typing import NamedTuple

import numpy as np
import pytest

from climt._core.parameters import per_column


def test_per_column_broadcasts_scalars_and_keeps_arrays():
    P = namedtuple("P", "a b")
    out = per_column(P(a=2.0, b=np.arange(3.0)), 3)
    assert type(out) is P and out.a.shape == (3,) and out.a.dtype == np.float64
    np.testing.assert_array_equal(out.a, [2.0, 2.0, 2.0])
    np.testing.assert_array_equal(out.b, [0.0, 1.0, 2.0])


def test_per_column_gives_int64_for_int_annotated_fields():
    class Q(NamedTuple):
        n: int
        x: float

    out = per_column(Q(n=3, x=np.array([1.5])), 4)
    assert type(out) is Q
    assert out.n.dtype == np.int64 and out.n.shape == (4,)
    assert out.x.dtype == np.float64 and out.x.shape == (4,)
    np.testing.assert_array_equal(out.n, [3, 3, 3, 3])
    np.testing.assert_array_equal(out.x, [1.5, 1.5, 1.5, 1.5])


def test_per_column_broadcast_fields_are_read_only():
    P = namedtuple("P", "a")
    out = per_column(P(a=1.0), 5)
    assert not out.a.flags.writeable


def test_per_column_rejects_wrong_length():
    with pytest.raises(ValueError, match="b"):
        per_column(namedtuple("P", "a b")(a=1.0, b=np.ones(2)), 3)
