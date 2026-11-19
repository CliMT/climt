import pytest
from climt import FortranUnavailableError


@pytest.fixture
def component():
    raise FortranUnavailableError("fixture needs the extensions")


def test_raises_in_body():
    raise FortranUnavailableError("body needs the extensions")


def test_raises_in_fixture(component):
    pass


def test_other_import_errors_still_fail():
    raise ImportError("an unrelated import problem")
