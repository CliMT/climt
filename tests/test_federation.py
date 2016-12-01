import pytest
from climt import Federation


def test_initialize_empty_federation():
    Federation(initial_state={})


def test_initialize_empty_federation_with_lists():
    Federation(initial_state={}, prognostic=[], diagnostic=[], implicit=None,
               monitor=[])

if __name__ == '__main__':
    pytest.main([__file__])
